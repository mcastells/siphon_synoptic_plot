from datetime import datetime
import time

import cartopy as cart
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cmocean
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage

from cftime import num2pydate
from metpy.units import units
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS

# Set the date/time to plot. CFSR/CFSv2 is hourly data, so make sure you set this at a whole hour, i.e. 0 minutes and 0 seconds
dt = datetime(2011, 5, 28, 0)
# Set the region to plot. [west, east, south, north]
plot_extent = [260, 330, -70, -40]
# Set the region of the subset of data to retrieve. I'm plotting with orthographic projection at high latitude, so I extend this region in west/east direction to avoid gaps at the poleward corners of the plot. [west, east, south, north]
data_extent = [235, 355, -70, -40]

# If the date/time chosen is between 1 Jan 1979 and 31 Mar 2011, CFSR reanalysis data is retrieved. If it is between 1 Apr 2011 and present, CFSv2 operational analysis data is retrieved.
cfsr_dates = [datetime(1979, 1, 1, 0), datetime(2011, 3, 31, 23, 59)]
cfsrv2_dates = [datetime(2011, 4, 1, 0), datetime.utcnow()]

# The API calls to retrieve data often result in server errors. These are errors occurring on NOAA's servers and are out of our control. Our code reattempts the API calls until they return successfully. The following variables determine the number of seconds to wait between attempts and how many attempts to try before giving up and stopping the program. The delay between attempts is to avoid spamming the servers.
retry_delay = 10
retry_limit = 100

# Helper function for finding proper time variable
def find_time_var(var, time_basename='time'):
    for coord_name in var.coordinates.split():
        if coord_name.startswith(time_basename):
            return coord_name
    raise ValueError('No time variable found for ' + var.name)

def retrieve_ncss(cat_url, ds_name, dt):
    ncss = None
    count = 0
    while ncss is None and count < retry_limit:
        count += 1
        try:
            print(f'Retrieving catalog from {cat_url}...')
            cat = TDSCatalog(cat_url)
            #print(f'{cat.datasets = }\n')
            ds = cat.datasets[ds_name]
            print(f'Retrieving NCSS dataset for {ds_name}...')
            ncss = ds.subset()
        except Exception as e:
            print(f'Attempt #{count}: {e}')
            time.sleep(retry_delay)
    if ncss is not None:
        return ncss
    else:
        raise Exception(f'Retry limit ({retry_limit}) reached.')

def retrieve_data(ncss, query):
    data = None
    count = 0
    while data is None and count < retry_limit:
        count += 1
        try:
            print(f'Querying dataset...')
            # Obtain the data we've queried for
            data = ncss.get_data(query)
        except Exception as e:
            print(f'Attempt #{count}: {e}')
            time.sleep(retry_delay)
    if data is not None:
        return data
    else:
        raise Exception(f'Retry limit ({retry_limit}) reached.')

def set_query_params(query, dt, extent):
    query.time(dt)
    query.lonlat_box(west=extent[0], east=extent[1], south=extent[2], north=extent[3])
    query.accept('netcdf')
    return query

if dt > cfsr_dates[0] and dt < cfsr_dates[1]:
    base_url = 'https://www.ncei.noaa.gov/thredds/catalog/model-cfs_reanl_ts/'
    cat_url = f'{base_url}{dt:%Y}/{dt:%Y%m}/catalog.xml'
    ds_title = 'CFS Reanalysis'
elif dt > cfsrv2_dates[0] and dt < cfsrv2_dates[1]:
    base_url = 'https://www.ncei.noaa.gov/thredds/catalog/model-cfs_v2_anl_ts/'
    cat_url = f'{base_url}{dt:%Y}/{dt:%Y%m}/catalog.xml'
    ds_title = 'CFSv2 Operational Analysis'
#base_url = 'https://www.ncei.noaa.gov/thredds/model-gfs-g4-anl-files-old/'
#cat = TDSCatalog(f'{base_url}{dt:%Y%m}/{dt:%Y%m%d}/catalog.xml')
#ncss = cat.datasets[f'gfsanl_4_{dt:%Y%m%d}_{dt:%H}00_000.grb2'].subset()

ds_name = f'prmsl.l.gdas.{dt:%Y%m}.grib2'
ncss = retrieve_ncss(cat_url, ds_name, dt)
query = ncss.query()
query = set_query_params(query, dt, data_extent)
query.variables('Pressure_reduced_to_MSL_msl')
mslp_data = retrieve_data(ncss, query)

ds_name = f'tmp2m.l.gdas.{dt:%Y%m}.grib2'
ncss = retrieve_ncss(cat_url, ds_name, dt)
query = ncss.query()
query = set_query_params(query, dt, data_extent)
query.variables('Temperature_height_above_ground')
temp_data = retrieve_data(ncss, query)

ds_name = f'wnd10m.l.gdas.{dt:%Y%m}.grib2'
ncss = retrieve_ncss(cat_url, ds_name, dt)
query = ncss.query()
query = set_query_params(query, dt, data_extent)
query.variables('u-component_of_wind_height_above_ground',
                'v-component_of_wind_height_above_ground')
wind_data = retrieve_data(ncss, query)

ds_name = f'ocnsst.l.gdas.{dt:%Y%m}.grib2'
ncss = retrieve_ncss(cat_url, ds_name, dt)
query = ncss.query()
query = set_query_params(query, dt, data_extent)
query.variables('Potential_temperature_depth_below_sea_1_Hour_Average')
sst_data = retrieve_data(ncss, query)

def get_data_at_time(data, var_names, dt):
    time_var = data.variables[find_time_var(data.variables[var_names[0]])]
    # Convert number of hours since the reference time into an actual date
    data_time = num2pydate(time_var[:].squeeze(), time_var.units)
    # Get the indices where the specified date/time occurs. Each model run contains a set of model output date/times, so there can be multiple occurences of model outputs with the date/time we are looking for.
    data_time_indices = np.where(data_time == dt)
    # We get the *last* occurance of the date/time we are looking for (from the latest model run) because this would be presumably be the most accurate data available.
    last_time_index = (data_time_indices[0][-1], data_time_indices[1][-1])
    data_at_time = []
    for var_name in var_names:
        data_at_time.append(data.variables[var_name][:].squeeze()[last_time_index])
    lat = data.variables['lat'][:].squeeze()
    lon = data.variables['lon'][:].squeeze()
    return data_at_time, lat, lon

data_at_time, mslp_lat, mslp_lon = get_data_at_time(mslp_data, ['Pressure_reduced_to_MSL_msl'], dt)
mslp = units.Pa * data_at_time[0]

data_at_time, temp_lat, temp_lon = get_data_at_time(temp_data, ['Temperature_height_above_ground'], dt)
temp = units.K * data_at_time[0]

data_at_time, wind_lat, wind_lon = get_data_at_time(wind_data, ['u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground'], dt)
u_wind_10m = units('m/s') * data_at_time[0]
v_wind_10m = units('m/s') * data_at_time[1]

data_at_time, sst_lat, sst_lon = get_data_at_time(sst_data, ['Potential_temperature_depth_below_sea_1_Hour_Average'], dt)
sst = units.K * data_at_time[0]

# Convert winds to knots
u_wind_10m.ito('kt')
v_wind_10m.ito('kt')

# Combine 1D latitude and longitudes into a 2D grid of locations
mslp_lon_2d, mslp_lat_2d = np.meshgrid(mslp_lon, mslp_lat)
temp_lon_2d, temp_lat_2d = np.meshgrid(temp_lon, temp_lat)
wind_lon_2d, wind_lat_2d = np.meshgrid(wind_lon, wind_lat)
sst_lon_2d, sst_lat_2d = np.meshgrid(sst_lon, sst_lat)

# Smooth MSLP a little
# Be sure to only put in a 2D lat/lon or Y/X array for smoothing
smooth_mslp = ndimage.gaussian_filter(mslp, sigma=3, order=0) * units.Pa
smooth_mslp.ito('hPa')

# Smooth temp a little
smooth_temp = ndimage.gaussian_filter(temp, sigma=1, order=0) * units.K

# Set Projection of Data
datacrs = ccrs.PlateCarree()

# Set Projection of Plot
center_lon = (plot_extent[0] + plot_extent[1])/2.0
center_lat = (plot_extent[2] + plot_extent[3])/2.0
projection = ccrs.Orthographic(center_lon, center_lat)

# Create new figure
plt.style.use('dark_background')
fig = plt.figure(figsize=(25.60, 10.80))

# Add the map and set the extent
ax = plt.subplot(111, projection=projection)
plt.title(f'{ds_title} - MSLP (hPa), 2m Temperature (C), Wind Barbs (kt)'
          f' {dt:%d %B %Y %H:%MZ}', fontsize=16)
ax.set_extent(plot_extent)
ax.coastlines()
ax.add_feature(cart.feature.LAND, edgecolor='k', facecolor='#aaa')
gl = ax.gridlines(crs=datacrs, draw_labels=True, alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

# Add state boundaries to plot
states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                name='admin_1_states_provinces_lakes',
                                                scale='50m', facecolor='none')
ax.add_feature(states_provinces, edgecolor='black', linewidth=1)

# Add country borders to plot
country_borders = cfeature.NaturalEarthFeature(category='cultural',
                                               name='admin_0_countries',
                                               scale='50m', facecolor='none')
ax.add_feature(country_borders, edgecolor='black', linewidth=1)

# Plot SST Contours
cs3 = ax.contourf(sst_lon_2d, sst_lat_2d, sst.to(units('degC')), levels=100, cmap=cmocean.cm.thermal, transform=datacrs)
cbar = fig.colorbar(cs3).set_label(f'sea surface temperature (potential temperature at 5m depth) (C)', rotation=-90, va='bottom')

# Plot MSLP Contours
clev_mslp = np.arange(0, 1200, 4)
cs = ax.contour(mslp_lon_2d, mslp_lat_2d, smooth_mslp, clev_mslp, linewidths=1.5, linestyles='solid', transform=datacrs)
plt.clabel(cs, fontsize=10, inline=1, inline_spacing=10, fmt='%i', rightside_up=True, use_clabeltext=True)

# Plot 2m Temperature Contours
clevtemp = np.arange(-60, 101, 4)
cs2 = ax.contour(temp_lon_2d, temp_lat_2d, smooth_temp.to(units('degC')), clevtemp, colors='tab:red', linewidths=1.75, linestyles='dotted', transform=datacrs)
plt.clabel(cs2, fontsize=10, inline=1, inline_spacing=10, fmt='%i', rightside_up=True, use_clabeltext=True)

# Plot 10m Wind Barbs
barbs = ax.barbs(wind_lon_2d, wind_lat_2d, u_wind_10m.magnitude, v_wind_10m.magnitude, barbcolor='white', flagcolor='white', length=6, regrid_shape=20, pivot='middle', transform=datacrs)

plt.show()