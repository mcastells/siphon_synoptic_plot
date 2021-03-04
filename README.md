# README #

This code is based on an [example from Unidata](https://unidata.github.io/python-gallery/examples/MSLP_temp_winds.html) and is modified to use CFSR and CFSv2 data instead of GFS data and include SST in the plot.

### What is this repository for? ###

* Learning to use [Siphon](https://unidata.github.io/siphon/latest/index.html) to access subsets of data hosted on a THREDDS Data Server (TDS).

### How do I get set up? ###

* Install dependencies:
    conda install -c conda-forge matplotlib scipy cartopy metpy siphon
* Set the variables near the top of the file as desired for date/time and lat/lon extent
* Run the code:
    python3 main.py

### Notes ###

* The CFSR/CFSv2 data has a separate file for each variable, so a netcdf subset (NCSS) request needs to be made for each variable. There seems to be inconsistencies in the lat/lon gridding between the different variables, so the sets of lat/lon are kept separate for each variable.
* Although a specific date/time is specified in the NCSS request, the entire month of data is returned for some reason. The data returned has 4 dimensions: date/time of model run, date/time of model output, lat, lon. There can be multiple model runs that contain a model output for the date/time we are looking for. The get_data_at_time function gets the *last* occurance of this model output (from the latest model run) because this would presumably be the most accurate data available.
* At time of writing, siphon calls to retrieve data from the NOAA NCEI THREDDS server often results in server errors. The code is designed to repeat these attempts until it gets a successful result.
* To modify this code to plot different variables, you will need to determine the filenames for the variables you want and the variable names that are contained in those files. The files can be viewed by going to the URL for any specific month of data (e.g., https://www.ncei.noaa.gov/thredds/catalog/model-cfs_v2_anl_ts/2021/202101/catalog.html). The filenames are abbreviations for the variables that they contain. A translation of abbreviations to actual variable names can be found here in the "CFSR Hourly Timeseries" PDF at the bottom of this page: https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/climate-forecast-system-version2-cfsv2. The best way I have found to get the actual variable names used in the files is to enter the code that retrieves the GRIB file object, set a breakpoint after that line, run the code in debugging mode and then inspect that object. It will have an attribute called 'variables' that contains the list of the variable names in that file. You can then use those names as needed to finish writing your code.

### Who do I talk to? ###

* Marc Castells
* mcastells@coaps.fsu.edu