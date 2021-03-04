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

### Who do I talk to? ###

* Marc Castells
* mcastells@coaps.fsu.edu