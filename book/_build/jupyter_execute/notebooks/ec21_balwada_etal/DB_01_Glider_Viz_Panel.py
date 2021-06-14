#!/usr/bin/env python
# coding: utf-8

# # Glider Viz Panel

# ## Authors
# **Dhruv Balwada, Scott Henderson, Alison R Gray**
# - Author1 = {"name": "Dhruv Balwada", "affiliation": "School of Oceanography, University of Washington", "email": "dbalwada@uw.edu", "orcid": "0000-0001-6632-0187"}
# - Author2 = {"name": "Scott Henderson", "affiliation": "Earth and Space Sciences, University of Washington", "email": "scottyh@uw.edu", "orcid": "0000-0003-0624-4965"}
# - Author3 = {"name": "Alison R Gray", "affiliation": "School of Oceanography, University of Washington", "email": "argray@uw.edu", "orcid": "0000-0002-1644-7654"}  

# ## Purpose
# Many oceanographic observational platforms, such as gliders, Argo floats, or ships, collect measurements on complex spatio-temporal paths. These sampling patterns are sometimes determined by logistical choices, e.g. due to weather, or sometimes dictated by the underlying oceanic flow causing the platforms to drift around. Deriving insights from these observations is challenging, and researchers almost always need to manually inspect the data. One such manual task might be to distinguish between signals that result due to sampling patterns vs signals that are an actual representation of the environment being studied.
# For example, a feature drifting by a  stationary platform, which samples the same location repeatedly, might look similar to a slowly evolving feature that a platform traverses through. These distinctions can be made by carefully inspecting the path of the instrument in reference to background measurements, like those from satellites, and this task can be made easier by using interactive visualizations. 
# 
# This notebook presents a workflow for setting up a visualization dashboard to interact with and analyze ocean glider data in combination with other available background measurements (eg. sea surface height, climatology, etc) that might be available in the geographical region sampled by the glider. This interactive visualization dashboard allows for easier data exploration and accelerates the rate at which insights can be derived. The primary purpose of this notebook is  not to simply build a standalone visualization dashboard for our specific data set, but rather provide an example that can be easily adapted to the specific use cases that individual researchers might have.
# 
# We visualize data for a particular glider campaign that was conducted in the Southern Ocean from May to August 2019, where two SeaGliders sampled the top 1000m of the ocean water column (a profile every 4 hours) in close proximity to each other. Sea Surface Height (SSH) and Finite Size Lypunov Exponents (FSLE, a measure of strain in the flow) are available from satellite based measurements and help provide background context in the sampling region. 
# 
# ## Technical contributions
# The main technical contributions of this notebook are as follows:
# - demonstrates how [GliderTools](https://glidertools.readthedocs.io/en/latest/) in combination with Scipy interpolation routines can be used to grid data collected by gliders, which is stored as point measurements.
# - demonstrates how [Xarray](http://xarray.pydata.org/) and the [Holoviz](https://holoviz.org/) visualization ecosystem, such as panel, hvplot, geoviews etc, can be used in combination to easily setup an interactive visualization tool to analyze glider data. 
# 
# ## Methodology
# The glider data is a collection of point measurements in space and time, where each measured variable has an associated position (longitude, latitude, depth) and time stamp. The glider data can be loaded and quality controlled using [GliderTools](https://glidertools.readthedocs.io/en/latest/). For making visualization easy and more interpretable this data is first gridded onto two regular grids; one is a time vs depth grid (T-grid) and the other other is an along track distance vs depth grid (D-grid). The T-grid helps visualize the duration of time over which different signals are observed, while the D-grid helps visualize the spatial size of the objects in the platform following coordinate. The gridding also helps colocate variables that might be slightly offset due to the different sampling frequencies of different sensors.  This gridded data set is then visualized along with the SSH and FSLE data sets, which are already available on a uniform grid.
# 
# The visualization dashboard has three main parts. The first part is a set of widgets that allow the user to make choices, such as selecting which variable to plot, the range of time or distance, choice of colormaps etc. The second part is a spatial plot of the surface properties (SSH or FSLE) with the glider tracks overlaid on it, where the selected glider and time slice to plot are based on the widgets. The third part is a plot of the glider section, where the exact range of time or distances plotted and the variables plotted are based on the variables. 
# 
# The visualization is made possible using the [Holoviz](https://holoviz.org/) libraries. Each widget selection sets the value of a certain parameter using the [param](https://param.holoviz.org/index.html) library, where parameters could be the variable to plot, the choice of colormap, etc. The plots are made possible using the libraries like [HoloViews](https://holoviews.org/), [hvPlot](https://hvplot.holoviz.org/), or [GeoViews](http://geoviews.org/). The parameter values are relayed to the plots using [Panel](https://panel.holoviz.org/), which also controls of the layout of the different components and manages the updates when interactions take place.
# 
# 
# ## Results
# The main purpose of this notebook is to demonstrate how a dashboard to explore glider data can be built using open source packages. The code is demonstrated in detail below and will hopefully serve as a useful example for others to explore their own datasets. The final cell of the notebook should result in dashboard that can be displayed in a notebook cell or in a standalone browser window. Some of the features of this dashboard are showcased in the gif below. The dashboard can also be accessed directly, without running all the cells in this notebook, by using the link to the 'Binder panel' in the repo readme; this makes the dashboard particularly appealing for sharing with collaborators that do not want to run the code.
# 
# ![Animation showing how to use some of the different features of the glider dashboard.](glider_dashboard.gif)
# 
# ## Funding
# 
# - Award1 = {"agency": "US National Science Foundation", "award_code": "OCE-1756882", "award_URL": "https://www.nsf.gov/awardsearch/showAward?AWD_ID=1756882"}
# 
# ## Keywords
# 
# keywords=["gliders", "interactive visualization", "ocean tracers", "submesoscale variability"]
# 
# ## Citation
# Balwada, Henderson, & Gray 2021. Interactive visualization tools for ocean glider data. Accessed at https://github.com/dhruvbalwada/ec2021_balwada_etal
# 
# 
# ## Suggested next steps
# 
# - We recommend that readers try to adapt this dashboard for their own datasets! 
# - To generalize the current dashboard code to use with additional datasets, it would be useful to add an option to read remote web-optimized formats instead of storing preprocessed netcdf datasets with the dashboard code. 
# - A few possible (non-essential) user experience improvements were identified (see wishlist on https://github.com/dhruvbalwada/glider-panel-demo), which can incorporated in future releases.
# 
# ## Acknowledgements 
# This project was supported by University of Washington's [eScience Institute](https://escience.washington.edu/) as part of the the 2021 [Winter Incubator](https://escience.washington.edu/winter-2021-incubator-projects/), and benefited from conversations with [Rob Fatland](https://escience.washington.edu/people/rob-fatland/) and [Don Setiawan](https://escience.washington.edu/people/landung-don-setiawan/) during the incubator. 
# We would also like to thank [Lily Dove](https://www.gps.caltech.edu/people/lilian-a-lily-dove?back_url=%2Fpeople%3Fcategory%3D11) and [Andrew Thompson](http://web.gps.caltech.edu/~andrewt/) from Caltech who processed a lot of the glider data that is used here, and a version of this data in a kriging mapped format can be found at https://www.ncei.noaa.gov/archive/accession/0228185 and https://www.ncei.noaa.gov/archive/accession/0228187. 

# # Setup
# 
# ## Library import
# Import all the required Python libraries.

# In[ ]:


# For data manipulation
import numpy as np
import xarray as xr
import pandas as pd
import glidertools as gt
import os
from scipy.interpolate import griddata

# For Visualization
import panel as pn
import holoviews as hv
from holoviews import opts
import geoviews as gv
import param
import matplotlib.pyplot as plt

## Import hvplot apis for xarray and pandas
import hvplot.xarray
import hvplot.pandas  


# # Data import
# 
# ## Import glider data
# 
# Two gliders were deployed in the Southern Ocean from May-August 2019 as part of an experiment called SOGOS (Southern Ocean Glider Observations of the Submesoscales). Here we use the data from these gliders, which is provided in the `data` folder. This data is stored in single netcdf files for each sensor, where each file contains the data for all the glider dives. This data can be opened directly using xarray.
# 
# Alternatively, often glider data is provided in single netcdf file for each dive. This dive data could be loaded using the instructions at https://glidertools.readthedocs.io/en/latest/loading.html. 

# In[ ]:


# locate the data folder
data_folder = './data'


# In[ ]:


# open glider files

ds_CTD_659 = xr.load_dataset(os.path.join(data_folder , 'sg659', 'CTD_659.nc'))
ds_CTD_660 = xr.load_dataset(os.path.join(data_folder , 'sg660', 'CTD_660.nc'))

ds_O2_659 = xr.load_dataset(os.path.join(data_folder , 'sg659', 'O2_659.nc'))
ds_O2_660 = xr.load_dataset(os.path.join(data_folder , 'sg660', 'O2_660.nc'))

ds_Chl_659 = xr.load_dataset(os.path.join(data_folder , 'sg659', 'Chl_659.nc'))
ds_Chl_660 = xr.load_dataset(os.path.join(data_folder , 'sg660', 'Chl_660.nc'))


# These data files are stored as 1D arrays of measurements at each observation point, where the location (longitude, latitude, and depth) and time of measurement for the observation point are also part of the data set.

# In[ ]:


# print to see what the data format is
ds_CTD_659


# ## Import surface data 
# 
# The sea surface height (SSH) and finite scale Lyapunov exponent (FSLE) datasets were obtained from the Copernicus and Aviso websites respectively, and manually cut for the region and time period of the glider deployments. Other data sets can also be similarly accessed; some community datasets could also be directly accessed using different webservices.  

# In[ ]:


# open SSH and FSLE files
ds_ssh = xr.open_dataset(os.path.join(data_folder, 'SSH_sogos.nc'))
ds_fsle = xr.open_dataset(os.path.join(data_folder, 'FSLE_sogos.nc'))


# In[ ]:


# print to see what the data format is
ds_ssh


# # Data processing and analysis
# 
# We first show how the data sets are processed to a form that is easily digestable by the visualization libraries. Then we show how the visualization libraries can be used to easily setup an interactive dashboard.

# ## Data processing

# ### Surface Data Processing

# In[ ]:


def datetime2ytd(time):
    """" Return time in YTD format from datetime format."""
    return  (time - np.datetime64('2019-01-01'))/np.timedelta64(1, 'D')
    


# In[ ]:


# Create variables for quiver plot
# Quiver plot requires the vectors to be defined in a very specific format.
ds_ssh['mag'] = np.sqrt(ds_ssh.ugos**2 + ds_ssh.vgos**2)
ds_ssh['angle'] = (np.pi/2.) - np.arctan2(ds_ssh.ugos/ds_ssh['mag'], 
                                          ds_ssh.vgos/ds_ssh['mag'])

# Create a new coordinate with time in year day units, as it is easier to work with.
# In future version handling of regular datetime formats could be introduced.
ds_ssh = ds_ssh.assign_coords(days = datetime2ytd(ds_ssh.time))
ds_fsle = ds_fsle.assign_coords(days = datetime2ytd(ds_fsle.time))

del ds_ssh.attrs['_NCProperties']
# need to delete this attribute because of the issue:
# https://github.com/pydata/xarray/issues/2822


# ### Glider Data Processing

# In[ ]:


# Convert glider time axis also to year day units, 
# so it matches that units for surface properties.

ds_CTD_659['days'] = datetime2ytd(ds_CTD_659.time)
ds_O2_659['days']  = datetime2ytd(ds_O2_659.time)
ds_Chl_659['days'] = datetime2ytd(ds_Chl_659.time)

ds_CTD_660['days'] = datetime2ytd(ds_CTD_660.time)
ds_O2_660['days']  = datetime2ytd(ds_O2_660.time)
ds_Chl_660['days'] = datetime2ytd(ds_Chl_660.time)


# In[ ]:


# Calculate along track distance
dXdist = gt.utils.distance(ds_CTD_659.longitude, ds_CTD_659.latitude)/1e3 # Convert to km
ds_CTD_659['distance'] = xr.DataArray(np.nancumsum(dXdist), 
                                       dims=ds_CTD_659.dims,
                                       coords=ds_CTD_659.coords)

dXdist = gt.utils.distance(ds_CTD_660.longitude, ds_CTD_660.latitude)/1e3
ds_CTD_660['distance'] = xr.DataArray(np.nancumsum(dXdist), 
                                       dims=ds_CTD_660.dims,
                                       coords=ds_CTD_660.coords)


# In[ ]:


# Group and average locations by dives 
# This makes plotting of locations on a map much faster, as
# there are less points to plot.
# These are used is only for plotting on a 2D map,
# where the depth coordinate is compressed.
ds_659_locs = xr.Dataset()
ds_660_locs = xr.Dataset()

ds_659_diveav = ds_CTD_659.groupby('dives').mean()
ds_660_diveav = ds_CTD_660.groupby('dives').mean()

ds_659_locs['longitude'] = ds_659_diveav.longitude
ds_659_locs['latitude']  = ds_659_diveav.latitude
ds_659_locs['days']      = ds_659_diveav.days
ds_659_locs['distance']  = ds_659_diveav.distance

ds_660_locs['longitude'] = ds_660_diveav.longitude
ds_660_locs['latitude']  = ds_660_diveav.latitude
ds_660_locs['days']      = ds_660_diveav.days
ds_660_locs['distance']  = ds_660_diveav.distance

# convert to pandas dataframe as it is much easier to handle in holoviz for traj data.
ds_659_locs = ds_659_locs.to_dataframe()
ds_660_locs = ds_660_locs.to_dataframe()


# In[ ]:


# Estimate additional derived variables
# Here we estimate density from the CTD measurements 
ds_CTD_659['potdens'] = gt.physics.potential_density(ds_CTD_659.salinity, 
                                                     ds_CTD_659.temperature, 
                                                     ds_CTD_659.pressure, 
                                                     ds_CTD_659.latitude, 
                                                     ds_CTD_659.longitude)

ds_CTD_660['potdens'] = gt.physics.potential_density(ds_CTD_660.salinity, 
                                                     ds_CTD_660.temperature, 
                                                     ds_CTD_660.pressure,
                                                     ds_CTD_660.latitude, 
                                                     ds_CTD_660.longitude)

# we can add mixed layer depth, N2 etc in the future versions


# Now that we have calculated some extra variables and new coordinates, we will go to the step of gridding the data onto a regular grid. 
# 
# In addition to the above tasks additional quality control procedures can be introduced at this stage (prior to gridding), using the QC procedures that are part of GliderTools (https://glidertools.readthedocs.io/en/latest/quality_control.html). 

# In[ ]:


# Make functions that put the point measurements from the glider 
# onto a regular grid. 
# There are many ways this can be done. Here we choose a simple linear interpolation 
# in time and pressure/depth or along-track distance and pressure/depth.

# Note this is different from how Glidertools goes gridding at the moment.
# These functions might be absorbed into Glidertools in future releases.

def interp_pres_time(ds_glid, var): 
    """ Return data variable interpolated onto a pressure-time grid.
    
    Keyword argument.
    ds_glid -- dataset of the glider data
    var     -- variable that needs to be interpolated
    """
    pres_ug = ds_glid.pressure
    time_ug = ds_glid.days
    
    # convert to points values
    points = np.stack([time_ug.values, pres_ug.values], axis=1)
    values = ds_glid[var].values
    
    # remove nans
    non_nan = np.logical_and(np.logical_and(
                                    ~np.isnan(points[:,0]), 
                                    ~np.isnan(points[:,1])),
                              ~np.isnan(values))
    
    points =points[non_nan,:]
    values =values[non_nan]
    
    # define grid 
    # In the future this can be made into an input from the users
    pres_grid = np.linspace(0,1000,251) 
    time_grid = np.arange(119, 207, 2/24)
    grid_p, grid_t = np.meshgrid(pres_grid, time_grid)
    
    temp_grided = griddata(points, 
                           values, 
                           (grid_t, grid_p), 
                           method='linear', 
                           rescale=True)
    
    return xr.DataArray(temp_grided.T, 
                        dims=["pressure", "time"],
                        coords={"pressure":pres_grid, "time":time_grid}
                       ).rename(var)


def interp_pres_dist(ds_glid, var): 
    """ Return data variable interpolated onto a pressure-along track distance grid.
    
    Keyword argument.
    ds_glid -- dataset of the glider data
    var     -- variable that needs to be interpolated
    """
    pres_ug = ds_glid.pressure
    dist_ug = ds_glid.distance
    
    # convert to points values
    points = np.stack([dist_ug.values, pres_ug.values],axis=1)
    values = ds_glid[var].values
    
    # remove nans
    non_nan = np.logical_and(np.logical_and(
                                    ~np.isnan(points[:,0]), 
                                    ~np.isnan(points[:,1])),
                             ~np.isnan(values))
    
    points =points[non_nan,:]
    values =values[non_nan]
    
    # define grid
    # In the future this can be made into an input from the users    
    pres_grid = np.linspace(0,1000,251)
    dist_grid = np.arange(0, dist_ug.max().values, 3)
    grid_p, grid_d = np.meshgrid(pres_grid, dist_grid)
    
    temp_grided = griddata(points, 
                           values, 
                           (grid_d, grid_p), 
                           method='linear', 
                           rescale=True)
    
    return xr.DataArray(temp_grided.T, 
                        dims=["pressure", "distance"],
                        coords={"pressure":pres_grid, "distance":dist_grid}
                       ).rename(var)


def convert_glider_time_pres(ds_glid, vars_convert= ['temperature','salinity','potdens','spice']):
    """ Return a dataset gridded onto pressure vs time grid.
    
    This is a helper function to apply gridding to multiple glider variables.
    """    
    ds_grid = xr.Dataset()
    
    for v in vars_convert:
        ds_grid[v] = interp_pres_time(ds_glid, v)
        print('Gridded ' + v)
        
    return ds_grid

def convert_glider_dist_pres(ds_glid, vars_convert= ['temperature','salinity','potdens','spice']):
    """ Return a dataset gridded onto pressure vs along track distance grid.
    
    This is a helper function to apply gridding to multiple glider variables.
    """ 
    ds_grid = xr.Dataset()
    
    for v in vars_convert:
        ds_grid[v] = interp_pres_dist(ds_glid, v)
        print('Gridded ' + v)
    
    return ds_grid


# In[ ]:


# Some variables measured by different sensors might be at different points. 
# Here we use a simple interpolation, from numpy, to collocate the Oxygen and Chlorophyll
# measurements to the CTD data point.

ds_CTD_659['oxygen'] = xr.DataArray(np.interp(ds_CTD_659.days, ds_O2_659.days, ds_O2_659.oxygen),
                                    dims = ds_CTD_659.dims, 
                                    coords = ds_CTD_659.coords
                                   ).rename('oxygen')
ds_CTD_659['Chl'] = xr.DataArray(np.interp(ds_CTD_659.days, ds_Chl_659.days, ds_Chl_659.Chl),
                                 dims = ds_CTD_659.dims, 
                                 coords = ds_CTD_659.coords
                                ).rename('Chl')
ds_CTD_660['oxygen'] = xr.DataArray(np.interp(ds_CTD_660.days, ds_O2_660.days, ds_O2_660.oxygen),
                                    dims = ds_CTD_660.dims, 
                                    coords = ds_CTD_660.coords
                                   ).rename('oxygen')
ds_CTD_660['Chl'] = xr.DataArray(np.interp(ds_CTD_660.days, ds_Chl_660.days, ds_Chl_660.Chl),
                                 dims = ds_CTD_660.dims, 
                                 coords = ds_CTD_660.coords
                                ).rename('Chl')


# **Note:** The next 4 cells can take a lot of time to run as a lot of compute heavy interpolations need to be done. 
# At this point the `load_flag` is set to 1, which will bypass the code in these cells without doing anything. The data that is generated here has been saved and provided with the repo, which will be loaded in the 5th cell below. 
# However, incase you want to run these set the `load_flag` to 0. 
# These cells will run with some waiting on local machines; while I have managed to get these to run on Binder, it often crashes due to the 2gb limit on Binder memory. 

# In[ ]:


get_ipython().run_cell_magic('time', '', "# convert from point data to gridded data\n# This cell will take the most time to run \n# (~10mins on laptop, ~20 mins on Binder). \n\nload_flag = 1 # Set this to 1 for loading the data instead of running the cells below, or 0 if you want to see what these cells do\nif load_flag == 0:\n    ds_659_Tgrid = convert_glider_time_pres(ds_CTD_659, vars_convert= ['temperature','salinity','potdens', 'oxygen', 'Chl'])\n    ds_660_Tgrid = convert_glider_time_pres(ds_CTD_660, vars_convert= ['temperature','salinity','potdens', 'oxygen', 'Chl'])\n\n    ds_659_Dgrid = convert_glider_dist_pres(ds_CTD_659, vars_convert= ['temperature','salinity','potdens', 'oxygen', 'Chl'])\n    ds_660_Dgrid = convert_glider_dist_pres(ds_CTD_660, vars_convert= ['temperature','salinity','potdens', 'oxygen', 'Chl'])\n\n\n# Alternatively users can separate the data processing and data \n# visualization sections into separate notebooks, and the gridded output can be generated once \n# and saved as netcdf files (using xarray's .to_netcdf() option). These netcdf files can then \n# directly be read into a visualization only notebook. \n# An example of how to do this is available at: https://github.com/dhruvbalwada/glider-panel-demo\n# This is how the data sets were saved, which are being loaded if load_flag is set to 1. ")


# In[ ]:


# Estimate an anomaly field based on time mean. 
# This is just an additional variable that we were interested in looking at.
# Could be defined in more complex ways too, like choose climatology as mean.
if load_flag == 0:
    ds_659_Tgrid_anomaly = ds_659_Tgrid - ds_659_Tgrid.mean('time')
    ds_660_Tgrid_anomaly = ds_660_Tgrid - ds_660_Tgrid.mean('time')

    ds_659_Dgrid_anomaly = ds_659_Dgrid - ds_659_Dgrid.mean('distance')
    ds_660_Dgrid_anomaly = ds_660_Dgrid - ds_660_Dgrid.mean('distance')


# In[ ]:


# Estimate the distance axis that goes with the time axis
# The gridding to a time axis was done for a uniform time grid,
# so the associated distance axis will likely be non-uniform.
if load_flag == 0:
    ds_659_Tgrid_loc = convert_glider_time_pres(ds_CTD_659, vars_convert= ['latitude','longitude'])
    ds_660_Tgrid_loc = convert_glider_time_pres(ds_CTD_660, vars_convert= ['latitude','longitude'])

    dXdist = gt.utils.distance(ds_659_Tgrid_loc.longitude.mean('pressure'), 
                               ds_659_Tgrid_loc.latitude.mean('pressure'))/1e3
    ds_659_Tgrid['distance'] = np.nancumsum(dXdist)
    ds_659_Tgrid_anomaly['distance'] = np.nancumsum(dXdist)

    dXdist = gt.utils.distance(ds_660_Tgrid_loc.longitude.mean('pressure'), 
                               ds_660_Tgrid_loc.latitude.mean('pressure'))/1e3

    ds_660_Tgrid['distance'] = np.nancumsum(dXdist)
    ds_660_Tgrid_anomaly['distance'] = np.nancumsum(dXdist)


# In[ ]:


# Estimate the time axis that goes with the along track distance 
# Similar to above cell, but now the non-uniform time axis that goes 
# with the uniform distance axis. 
if load_flag == 0:
    temp = convert_glider_dist_pres(ds_CTD_659, vars_convert=['days'])
    ds_659_Dgrid['time'] = temp.days.mean('pressure').values

    temp = convert_glider_dist_pres(ds_CTD_660, vars_convert=['days'])
    ds_660_Dgrid['time'] = temp.days.mean('pressure').values


# In[ ]:


# Incase you don't have the time to run the above 4 cells
# Or if binder keeps crashing trying to run the above 4 cells.

if load_flag == 1:
    ds_659_Tgrid = xr.open_dataset(os.path.join(data_folder, '659_Tgrid.nc'))
    ds_660_Tgrid = xr.open_dataset(os.path.join(data_folder, '660_Tgrid.nc'))

    ds_659_Tgrid_anomaly = xr.open_dataset(os.path.join(data_folder, '659_Tgrid_anomaly.nc'))
    ds_660_Tgrid_anomaly = xr.open_dataset(os.path.join(data_folder, '660_Tgrid_anomaly.nc'))

    ds_659_Dgrid = xr.open_dataset(os.path.join(data_folder, '659_Dgrid.nc'))
    ds_660_Dgrid = xr.open_dataset(os.path.join(data_folder, '660_Dgrid.nc'))

    ds_659_Dgrid_anomaly = xr.open_dataset(os.path.join(data_folder, '659_Dgrid_anomaly.nc'))
    ds_660_Dgrid_anomaly = xr.open_dataset(os.path.join(data_folder, '660_Dgrid_anomaly.nc'))


# At this point in the notebook all the data processing steps have been executed, and the data sets are in a format that is ready to be visualized. 

# ## Setting up the interactive dashboard
# 
# This is the main contribution of this submission. Here we show how Holoviz libraries can be used to create an interactive visualization of glider data. The visualization choices here are based on our particular use case, and we accordingly chose the variables that will be plotted and the widgets that are available. However, this notebook should be viewed more as an example that can be modified to the particular visualization use case that others might have. 

# In[ ]:


# Create variable maps 
# These are dictionaries linking variable names to particular data sets

# The different gliders
glider_nums = ['sg659', 'sg660']

# The different surface variables
surface_var_map = {
              'SSH' : ds_ssh['adt'],
              'SSHA': ds_ssh['sla'],
              'FSLE': ds_fsle['fsle_max']
              }

# The different variables in a particular glider data set
glider_vars = list(ds_659_Tgrid.keys()) # just need to do once bevause all glider data sets here have same variables

# The different colormaps available
cmap_options = plt.colormaps()

# Dictionary linking variables to default properties
# Here we only define a default colormap, but other defaults 
# can be added.
var_select_map = {
            'oxygen': {'cmap_sel': 'YlOrBr' },
            'Chl': {'cmap_sel': 'Greens' },
            'salinity': {'cmap_sel': 'YlGnBu'},
            'temperature': {'cmap_sel': 'RdBu_r'},
            'potdens': {'cmap_sel': 'Purples' }
             }
# For future versions would be nice if some of these things could come from the attributes. 

# Different data sets for each glider
glider_map = {
            'sg659': {'Time grid': ds_659_Tgrid, 'Distance grid': ds_659_Dgrid, 'loc': ds_659_locs},
            'sg660': {'Time grid': ds_660_Tgrid, 'Distance grid': ds_660_Dgrid, 'loc': ds_660_locs},
             }

glider_map_anom = {
            'sg659': {'Time grid': ds_659_Tgrid_anomaly, 'Distance grid': ds_659_Dgrid_anomaly},
            'sg660': {'Time grid': ds_660_Tgrid_anomaly, 'Distance grid': ds_660_Dgrid_anomaly},
             }


# The dashboard will be made using three classes that will correspond to the:
# - Widgets: `GliderParams`: containing definitions of all the widgets.
# - Trajectory plot: `GliderTrajectoryPlot`: sets up the plot with the glider trajectories overlaid on the surface variable plots and bahtymetry.
# - Glider section plot: `GliderVerticalSectionPlot`: sets up the plot for the glider section.
# 
# The plotting classes will inherit the class with the widgets, and then a combined class (`GliderCombinedPlot`) will inherit the plotting classes. The dashboard will be an object of this combined class, which can be passed to panel for being layed out.

# In[ ]:


class GliderParams(param.Parameterized):
    """ Class containting all the parameters for the widgets, and some default methods. """
    
    surface_var       = param.Selector(surface_var_map.keys(), default='SSH',
                                       label='Surface Field', precedence=0)
    glider_num        = param.Selector(glider_map.keys(), default='sg659',
                                       label='Glider Num', precedence=0)
    time_slider       = param.Range(label='Days in 2019', 
                                    bounds=(119, 205), 
                                    default=(119, 135), precedence=3)
    alpha_slider      = param.Magnitude(label='Transparency', precedence=4)
    glider_grid       = param.Selector(['Time grid', 'Distance grid'], default='Time grid', 
                                       label='Grid Type', precedence=0)
    glider_var        = param.Selector(glider_vars, default='temperature', 
                                       label='Glider Variable', precedence=1)
    var_colormap      = param.Selector(default='RdBu_r', objects=cmap_options, 
                                       label='Glider Section Colormap', precedence=2)
    distance_slider   = param.Range(label='Along Track Distance',
                                    bounds=(0, 2.2e3), default=(0, 400), 
                                    precedence=-1) # start with a negative precedence, in accordance with default being Tgrid
    anomaly_boolean   = param.Boolean(default=False, label='Anomaly', precedence=3)
    density_boolean   = param.Boolean(default=True, label='Show Density Contours', precedence=4)
    density_range     = param.Range(label='Density range', bounds=(1026.8, 1027.9), default=(1026.8, 1027.9),precedence=10)
    density_gradation = param.Integer(label='Density levels', default=11, bounds=(2, 21),precedence=10)
    
    def _set_tools(self, plot, element):
        """ Method to not revome default active toolbars. """
        plot.state.toolbar.active_drag = None
        plot.state.toolbar.active_inspect = None
    
    @param.depends('glider_var', watch=True)
    def _update_colormap(self):
        """ Update default colormap choices with changing variables. """
        self.var_colormap = var_select_map[self.glider_var]['cmap_sel']
    
    # The next couple of methods toggle the widgets visible or not. 
    @param.depends('density_boolean', watch=True)
    def _update_density_widgets(self):
        """ Remove density widgets when not being used. """ 
        if self.density_boolean:
            self.param.density_range.precedence=10
            self.param.density_gradation.precedence=10
        else:
            self.param.density_range.precedence=-1
            self.param.density_gradation.precedence=-1
            
    @param.depends('glider_grid', watch=True)
    def _update_grid_widgets(self):
        """ Show time or distance widgets based on grid choice. """
        if self.glider_grid == 'Time grid':
            self.param.time_slider.precedence=3
            self.param.distance_slider.precedence=-1
        else: 
            self.param.time_slider.precedence=-1
            self.param.distance_slider.precedence=3


# In[ ]:


class GliderTrajectoryPlot(GliderParams):
    """ Class containing the setup for the trajectory plot. """
    
    @param.depends('glider_num', 'time_slider', 'distance_slider')
    def plot_traj(self):
        """ Plot glider trajectories. """
        time_rng = self.time_slider
        dist_rng = self.distance_slider
        
        ###
        # For the selected glider do the proper time vs distance conversion
        # but for the unselected glider and surface plots always stick to the
        # corresponding time 
        ds = glider_map[self.glider_num]['loc']
        if self.glider_grid=='Time grid':
            ds_tsel = ds.loc[(ds.days>=time_rng[0]) & (ds.days<=time_rng[1])]
            dsel = (ds_tsel.iloc[0].distance, ds_tsel.iloc[-1].distance)
            self.distance_slider = dsel
        else:
            ds_tsel = ds.loc[(ds.distance>=dist_rng[0]) & (ds.distance<=dist_rng[1])]
            if int(ds_tsel.iloc[-1].days)<=205: # since the netcdf files for surface fields don't have a 206
                tsel = (int(ds_tsel.iloc[0].days), int(ds_tsel.iloc[-1].days))
            else:
                tsel = (int(ds_tsel.iloc[0].days), int(205))
            self.time_slider = tsel
        ###
        
        time_rng = self.time_slider # make sure time_rng has the most up to date values
        
        traj = {}
        for glid in glider_nums:
            ds = glider_map[glid]['loc']
            ds_tsel = ds.loc[(ds.days>time_rng[0]) & (ds.days<time_rng[1])]
        
            traj[glid] = ds_tsel.hvplot.points(geo=True,  x='longitude', y='latitude', 
                                               hover=True, hover_cols=['days'], 
                                               size=1)
        traj[self.glider_num].opts(size=2.5)
        
        return traj['sg659']*traj['sg660']
    
    def surf_tiles(self):
        """ Plot bathymetry tile. """
        gebco_tiles = 'https://tiles.arcgis.com/tiles/C8EMgrsFcRFL6LrL/arcgis/rest/services/GEBCO_basemap_NCEI/MapServer/tile/{Z}/{Y}/{X}'
        return gv.WMTS( gebco_tiles )
    
    @param.depends('time_slider')
    def surf_vec(self):
        """ Plot velocity vectors as quivers. """
        time_sel = self.time_slider[1] # show map for last day on time slider
        
        return ds_ssh.where(ds_ssh.days==time_sel, drop=True).squeeze('time'
                    ).hvplot.vectorfield(x='longitude', y='latitude', 
                                         angle='angle', mag='mag',
                                         geo=True, hover=False).opts(magnitude='mag')
    
    @param.depends('surface_var', 'time_slider', 'alpha_slider')
    def plot_surface(self):
        """ Plot surface variable of choice. """
        time_sel = self.time_slider[1] # show map for last day on time slider
        
        ds_all = surface_var_map[self.surface_var]
        ds = ds_all.where(ds_all.days==time_sel, drop=True).squeeze('time')
        if self.surface_var == 'FSLE':    
            surf_plot = ds.hvplot.image(geo=True)
            surf_plot.opts(clim=(-0.6,0), cmap='Blues_r', clabel='FSLE')
        elif self.surface_var == 'SSH':
            surf_plot = ds.hvplot.image(geo=True)
            surf_plot.opts(clim=(-1,0), cmap='cividis', clabel='SSH')
        else: 
            surf_plot = ds.hvplot.image(geo=True)
            surf_plot.opts(clim=(-0.3,0.3), cmap='RdBu_r', clabel='SSHA')
        
        surf_plot.opts(frame_width=450, alpha=self.alpha_slider, tools=['hover'], hooks=[self._set_tools])

        return surf_plot
        
   
    def view(self):
        """ Make lat-lon plot with surface vars, bathymetry, and glider tracks. """
        return (hv.DynamicMap(self.plot_surface)
                * hv.DynamicMap(self.surf_tiles)
                * hv.DynamicMap(self.surf_vec)
                * hv.DynamicMap(self.plot_traj))


# In[ ]:


class GliderVerticalSectionPlot(GliderParams):
    """ Class containing the setup for the glider section plot """
    
    @param.depends('density_range', 'density_gradation', 'glider_grid','glider_num')
    def density_contours(self):
        """ Plot the density contours. """
        #print('in contour')
        contour = glider_map[self.glider_num][self.glider_grid]['potdens'].hvplot.contour(
                        flip_yaxis=True, levels=np.linspace(self.density_range[0],
                        self.density_range[1],self.density_gradation)
                        ).opts(tools=[])
        return contour
    
    @param.depends('anomaly_boolean', 'glider_grid', 'glider_num', 'glider_var')
    def glider_image(self):
        """ Plot the image for the glider section. """
        # Change the data set if wanting to plot anomaly
        if self.anomaly_boolean:
            glid_ds = glider_map_anom
        else:
            glid_ds = glider_map

        # plot the image in Distance or Time
        if self.glider_grid=='Distance grid':
            image = hv.Image( 
                    (glid_ds[self.glider_num][self.glider_grid].distance, 
                     glid_ds[self.glider_num][self.glider_grid].pressure,
                     glid_ds[self.glider_num][self.glider_grid][self.glider_var]), 
                    ['Distance [km]', 'Pressure [dBar]'], self.glider_var)

        else:
            image = hv.Image( 
                    (glid_ds[self.glider_num][self.glider_grid].time, 
                     glid_ds[self.glider_num][self.glider_grid].pressure,
                     glid_ds[self.glider_num][self.glider_grid][self.glider_var]), 
                    ['Time [days]', 'Pressure [dBar]'], self.glider_var)
        
        # estimate the color range so that outliers don't create problems
        bin_range = np.nanpercentile(glid_ds[self.glider_num][self.glider_grid][self.glider_var], [.5,99.5])
        
        # set properties for image like colorbar etcs.
        image = image.opts(opts.Image(
                        colorbar=True,
                        cmap=self.var_colormap,
                        invert_yaxis=True,
                        clim=(bin_range[0], bin_range[1]),
                        width=800,
                        tools=['hover'], hooks=[self._set_tools]
                        ))
        return image

    def viewable(self):
        """ Make combined plot of the glider section. """
        image = hv.DynamicMap(self.glider_image).hist()
                
        title = (str(self.time_slider[0]) +'-'+ str(self.time_slider[1])
                 + ' Days & ' + str(int(self.distance_slider[0]))+'-'
                 + str(int(self.distance_slider[1])) + ' km')
        
        if self.glider_grid=='Distance grid':
            image.opts(opts.Image(xlim = self.distance_slider, title=title))
        else:
            image.opts(opts.Image(xlim = self.time_slider, title=title))

        # add the density contours or not. 
        if self.density_boolean:
            return image*hv.DynamicMap(self.density_contours)
        else:
            return image
            


# In[ ]:


class GliderCombinedPlot(GliderTrajectoryPlot, GliderVerticalSectionPlot):
    """ Acombined class that inherits the above classes and will be used to define the main object. """
    pass


# In[ ]:


dashboard = GliderCombinedPlot()


# In[ ]:


text_title = '## Southern Ocean Glider Observations of the Submesoscales\n Interactive dashboard to explore glider data collected in the Southern Ocean (zoom out on top panel to see exact location) during May-August 2019'
text_tip = '*Tip: The box select tool (from tools on top right of this plot) can be used to select a range on the histogram plot to the right, which adjusts the color limits for the section plot*'
dashboard_panel = pn.Column(pn.Row(
                                pn.Param(dashboard.param, name=''),
                                pn.Column(pn.panel(text_title),dashboard.view())), 
                            pn.Column(
                                dashboard.viewable, 
                                pn.panel(text_tip))) 


# You will need to comment out and uncomment the different lines below based on how you want to view the dashboard. By default we have it setup to be able to view the dashboard directly as a standalone app using the Binder panel button.

# In[ ]:


# To render in this notebook and make available as a panel app
dashboard_panel.servable()

# To run locally in a standalone browser window (this may not work on Binder)
#dashboard_panel.show()

# To display in notebook cell below
#dashboard_panel


# In[ ]:




