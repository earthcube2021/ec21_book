#!/usr/bin/env python
# coding: utf-8

# # Leveraging satellite data and pangeo to investigate the role of Gulf Stream frontal dynamics in ocean productivity

# ## Authors
# 
# - Author1 = {"name": "Patrick Clifton Gray", "affiliation": "Duke University Marine Lab", "email": "patrick.c.gray@duke.edu", "orcid": "0000-0002-8997-5255"}
# - Author2 = {"name": "David W. Johnston", "affiliation": "Duke University Marine Lab", "email": "david.johnston@duke.edu", "orcid": "0000-0003-2424-036X"}
#     

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# 
# <div class="toc">
#     <ul class="toc-item">
#         <li>
#             <span>
#                 <a href="#Leveraging-satellite-data-and-pangeo-to-investigate-the-role-of-Gulf-Stream-frontal-dynamics-in-ocean-productivity" data-toc-modified-id="Leveraging satellite data and pangeo to investigate the role of Gulf Stream frontal dynamics in ocean productivity"><span class="toc-item-num">1&nbsp;&nbsp;</span>Leveraging satellite data and pangeo to investigate the role of Gulf Stream frontal dynamics in ocean productivity</a>
#             </span>
#             <ul class="toc-item">
#                 <li>
#                     <span>
#                         <a href="#Authors" data-toc-modified-id="Author(s)-1.1">
#                             <span class="toc-item-num">1.1&nbsp;&nbsp;</span>
#                             Authors
#                         </a>
#                     </span>
#                 </li>
#                 <li><span><a href="#Purpose" data-toc-modified-id="Purpose-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Purpose</a></span></li>
#                 <li><span><a href="#Technical-contributions" data-toc-modified-id="Technical-contributions-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>Technical contributions</a></span></li>
#                 <li><span><a href="#Methodology" data-toc-modified-id="Methodology-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>Methodology</a></span></li>
#                 <li><span><a href="#Results" data-toc-modified-id="Results-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>Results</a></span></li>
#                 <li><span><a href="#Funding" data-toc-modified-id="Funding-1.6"><span class="toc-item-num">1.6&nbsp;&nbsp;</span>Funding</a></span></li>
#                 <li><span><a href="#Keywords" data-toc-modified-id="Keywords-1.7"><span class="toc-item-num">1.7&nbsp;&nbsp;</span>Keywords</a></span></li>
#                 <li><span><a href="#Citation" data-toc-modified-id="Citation-1.8"><span class="toc-item-num">1.8&nbsp;&nbsp;</span>Citation</a></span></li>
#                 <li><span><a href="#Acknowledgements" data-toc-modified-id="Acknowledgements-1.9"><span class="toc-item-num">1.9&nbsp;&nbsp;</span>Acknowledgements</a></span></li>
#             </ul>
#         </li>
#         <li>
#             <span>
#                 <a href="#Setup" data-toc-modified-id="Setup-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Setup</a>
#             </span>
#             <ul class="toc-item">
#                 <li><span><a href="#Library-import" data-toc-modified-id="Library-import-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>Library import</a></span></li>
#                 <li><span><a href="#Local-library-import" data-toc-modified-id="Local-library-import-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>Local library import</a></span></li>
#             </ul>
#         </li>
#         <li>
#             <span>
#                 <a href="#Data-preparation" data-toc-modified-id="Data-preparation-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Data preparation</a>
#             </span>
#             <ul class="toc-item">
#                 <li><span><a href="#Data-download" data-toc-modified-id="Data-download-3.1"><span class="toc-item-num">3.1&nbsp;&nbsp;</span>Data download</a></span></li>
#                 <li><span><a href="#Data-import" data-toc-modified-id="Data-import-3.2"><span class="toc-item-num">3.2&nbsp;&nbsp;</span>Data import</a></span></li>
#                 <li><span><a href="#Gulf-Stream-front-location" data-toc-modified-id="Gulf-Stream-front-location-3.3"><span class="toc-item-num">3.3&nbsp;&nbsp;</span>Gulf Stream front location</a></span></li>
#                 <li><span><a href="#Selecting-data-along-the-front" data-toc-modified-id="Selecting-data-along-the-front-3.4"><span class="toc-item-num">3.4&nbsp;&nbsp;</span>Selecting data along the front</a></span></li>
#             </ul>
#         </li>
#         <li>
#             <span>
#                 <a href="#Data-processing-and-analysis" data-toc-modified-id="Data-processing-and-analysis-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Data processing and analysis</a>
#             </span>
#             <ul class="toc-item">
#                 <li><span><a href="#Visualizing-all-datasets" data-toc-modified-id="Visualizing-all-datasets-4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>Visualizing all datasets</a></span></li>
#                 <li><span><a href="#Seasonality-along-the-Gulf-Stream" data-toc-modified-id="Seasonality-along-the-Gulf-Stream-4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>Seasonality along the Gulf Stream</a></span></li>
#                 <li><span><a href="#The-full-time-series" data-toc-modified-id="The-full-time-series-4.3"><span class="toc-item-num">4.3&nbsp;&nbsp;</span>The full time series</a></span></li>
#             </ul>
#         </li>
#         <li>
#             <span>
#                 <a href="#Suggested-next-steps" data-toc-modified-id="Suggested-next-steps-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>Suggested next steps</a>
#             </span>
#         </li>
#         <li>
#             <span>
#                 <a href="#References" data-toc-modified-id="References-6"><span class="toc-item-num">6&nbsp;&nbsp;</span>References</a>
#             </span>
#         </li>
#     </ul>   
# </div>

# **In a vessel floating on the Gulf Stream one sees nothing of the current and knows nothing but what experience tells him; but to be anchored in its depths far out of the sight of land, and to see the mighty torrent rushing past at a speed of miles per hour, day after day and day after day, one begins to think that all the wonders of the earth combined can not equal this one river in the ocean.” - J. E. Pillsbury, The Gulf Stream (1891)**
# 
# 
# ## Purpose
# 
# The Gulf Stream is a dominant physical feature in the North Atlantic, influencing the marine ecosystem across trophic levels, yet the physical-biological interactions along this western boundary current are poorly understood. In this work we analyze Gulf Stream physical and biological patterns over multiple years with particular emphasis on how the current's front impacts biology across space and time. We ask specifically how phytoplankton patterns along the front correspond to seasonality and wind with speculation on possible nutrient limitations. This work primarily uses satellite-based sea surface temperature (SST), chlorophyll-a (chla), and sea surface height (SSH) data. Using `xarray` and running in `jupyter` in the cloud allows efficient analysis across many years of satellite data and `holoviz` enables interactive visualization and data exploration. Basing this work on a [pangeo Docker image](https://github.com/pangeo-data/pangeo-docker-images) allows easy reproducibility and a robust set of packages for scientific computing. This work aims to provide insight into the overall role of Gulf Stream frontal features and dynamics in ocean productivity and increase the ease of access, manipulation, and visualization of these highly dimensional and complex datasets.
# 
# ## Technical contributions
# - collation and visualization of a variety of oceanographic datasets permitting investigation of the Gulf Stream front
# - derivation of the Gulf Stream front location every month for a decade from Sea Surface Height data
# - demonstration of an array of scientific computing tools for oceanography analysis
# - investigation of the physical and biological properties of the Gulf Stream front over a decade
# 
# ## Methodology
# The goal of this analysis is to investigate the physical and biological trends and relationships along the Gulf Stream front using SST, wind, and chla data. For comparison we also investigate that same data on the coastal and open ocean side of the front. The first half of the notebook pulls in and visualizes all the different datasets and then derives the front location based on the 0.25m isoline in sea surface height data. The final step in the first half of the notebook is to select data from the main datasets (SST, wind, chla) using these front lines at every time step.
# 
# The second half of the notebook then analyzes this along-front data using a variety of visualizations and interprets the results.
# 
# ## Results
# 
# Results are discussed in detail throughout the analysis section. Some main takeaways are that for nearly all of the Gulf Stream, from Florida out into the North Atlantic, there is a fall bloom that is comparable or greater than the spring bloom and this is true on both sides of the front. 
# 
# The northern region (further downstream) has a major increase in chla in the spring - the well studied spring bloom - and the region further south (upstream) around the South Atlantic Bight (SAB) has a relatively large increase throughout the whole winter, starting mid fall and staying high until the spring. The region around Cape Hatteras falls in the middle of these two patterns, but with generally higher chla, particularly on the coastal side, likely due to proximity to the coast and generally large outflows from Pamlico Sound and Chesapeake Bay injecting nutrients into the water. The sustained spikes in chla off Cape Hatteras near the front on the coastal side are likely drivers of the increased biodiversity in that region. It is also worth noting that the northerly section of the front (downstream) has a more notable mid-winter decrease, possibly when light begins to be limiting, and this is not as pronounced in the southerly section of the front (upstream). Though the light difference is not extreme. There are 9.5 hours of sunlight on the winter solstice at 37° N and 10.6 hours at 27° N.
# 
# The satellite imagery makes it clear that the Gulf Stream is a major boundary between the coastal and the open ocean ecosystems. A large increase in chla beyond the Gulf Stream front, in the Sargasso Sea, is evident from late fall to early spring. The spring bloom in the North Atlantic is also clear where it peaks in April and May. 
# 
# The SST data matches what one might expect, decreasing throughout the winter and increasing in summer with the Gulf Stream itself staying as the warmest water throughout all seasons. We can also see a few eddy footprints in SST and SSH even though this data is resampled to monthly time steps.
# 
# Overall the analysis generally agrees with the understanding of phytoplankton dynamics in the Atlantic ([Siegel et al 2002](https://doi.org/10.1126/science.1069174), [Fischer et al 2015](https://doi.org/10.5670/oceanog.2014.26)) but it doesn't follow exactly along with the commonly studied bloom in the North Atlantic. Fall production is more intense along the Gulf Stream compared to the North Atlantic, potentially supplied with nutrients by mesoscale upwelling ([Pascual et al 2014](https://doi.org/10.1002/2014GL062569)) and ageostrophic circuation at the front ([Levy et al 2018](https://doi.org/10.1038/s41467-018-07059-3)) along the front and supplied with seed populations from further south. This is more in agreement with ideas of winter biomass accumulation despite low light conditions (e.g. [Boss and Behrenfeld 2010](https://doi.org/10.1029/2010GL044174)). Wind patterns are noisy and don't lead to clear agreement or disagreement with the general assumptions that winter winds drive deeper mixing which inject nutrients supplying the bloom and stronger winds in the spring deepen mixed layers which delay the start of the subpolar spring ([Ueyama and Monger 2005](https://doi.org/10.4319/lo.2005.50.6.1820), [Henson et al 2009](https://doi.org/10.1029/2008JC005139)), but this is an avenue for further investigation with additional wind datasets.
# 
# ## Funding
# 
# - Award1 = {"agency": "US National Aeronautics and Space Administration", "award_code": "80NSSC19K1366", "award_URL": "https://nspires.nasaprs.com/external/viewrepositorydocument/cmdocumentid=701992/solicitationId=%7B913A7DEE-2747-6539-130C-0AB1E2322F42%7D/viewSolicitationDocument=1/Updated%20ESD%20FINESST19%20SELECTIONS%207.24.19.pdf"}
# 
# ## Keywords
# keywords=["gulf stream", "satellite oceanography", "xarray", "physical-biological", "ocean productivity"]
# 
# ## Citation
# Gray, P.C., & Johnston D.W., 2021. Leveraging satellite data and pangeo to investigate the role of Gulf Stream frontal dynamics in ocean productivity. Accessed 5/15/2021 at https://github.com/patrickcgray/gulf_stream_productivity_dynamics
# 
# ## Acknowledgements 
# 
# We acknowledge helpful conversations with many collaborators and mentors as well as the substantial effort of the open source community developing xarray, zarr, holoviz, geopandas, jupyter and the many other libraries we relied on for this analysis.
# 
# The template is licensed under a <a href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License.</a>

# # Setup
# 
# ## Library import

# In[1]:


# Data manipulation
import pandas as pd
import geopandas as gpd
import numpy as np
import xarray as xr

# time
import datetime

# downloading data from GDrive
from google_drive_downloader import GoogleDriveDownloader as gdd

# Options for pandas
pd.options.display.max_columns = 50
pd.options.display.max_rows = 30

# Visualizations
import hvplot.xarray
import hvplot.pandas
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# geographic and geometry tools
import cartopy.crs as crs
from shapely.ops import unary_union
from shapely.ops import transform
from shapely.geometry import Polygon
from affine import Affine

# there are many divide by zero warnings that we're going to suppress
# but advised to comment this line out when actively coding
import warnings
warnings.filterwarnings('ignore')


# ## Local library import

# In[2]:


# this is a package from https://github.com/GeoscienceAustralia/dea-notebooks that allows contours to be pulled from xarray datasets via geopanda points
from dea_spatial import subpixel_contours


# # Data preparation

# Data is all stored as zarr and downloaded from Google Drive in the next code block. This is done to ensure speed of operation for this notebook. Current data providers are too slow for use in an ephemeral Binder environment and some require manual downloading. All data can be downloaded from the original repositories via the information below if needed. All data use the same bounds, the latitude range is from 44° North to 26° North and the longitude range is from 82° West to 66° West.
# 
# #### Altimetry data 
# is from AVISO distributed via Copernicus and uses multiple altimeters to generate a 0.25 degree daily SSH and (sea level anomaly) SLA product. The data can be downloaded [here](https://resources.marine.copernicus.eu/?option=com_csw&task=results). After navigating to that site the user must manually search for SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046 and specify the bounds and time period (January 1st, 2010 through December 31st 2019).
# 
# #### Ocean Color data 
# is from the Ocean Colour Climate Change Initiative's multi-sensor global satellite chlorophyll-a product we use the CCI_ALL-v5.0-8DAY product which can viewed [here](https://www.oceancolour.org/about) and our exact dataset can be downloaded with [this link](https://www.oceancolour.org/thredds/ncss/CCI_ALL-v5.0-5DAY?var=chlor_a&north=42&west=-80&east=-70&south=30&disableProjSubset=on&horizStride=1&time_start=2014-01-01T00%3A00%3A00Z&time_end=2020-12-31T00%3A00%3A00Z&timeStride=1).
# 
# #### Sea Surface Temperature 
# is the GHRSST Level 4 OSPO Global Nighttime Foundation Sea Surface Temperature Analysis product. More can be found on this product [here](https://podaac.jpl.nasa.gov/dataset/Geo_Polar_Blended_Night-OSPO-L4-GLOB-v1.0) and all data can be downloaded [here](https://thredds.jpl.nasa.gov/thredds/catalog_ghrsst_gds2.html?dataset=Geo_Polar_Blended_Night-OSPO-L4-GLOB-v1.0).
# 
# #### Wind speed data 
# is Metop-A ASCAT, 0.25°, Global, Near Real Time, 2009-present (8 Day) available for download [here](https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdQAwind8day.graph). And for this exact dataset [this link](https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdQAwind8day.nc?x_wind%5B(2010-04-01T00:00:00Z):1:(2020-01-01T00:00:00Z)%5D%5B(10.0):1:(10.0)%5D%5B(26):1:(44)%5D%5B(278):1:(294)%5D,y_wind%5B(2010-04-01T00:00:00Z):1:(2020-01-01T00:00:00Z)%5D%5B(10.0):1:(10.0)%5D%5B(26):1:(44)%5D%5B(278):1:(294)%5D) can be used.

# ## Data download 
# Download all the necessary files from Google Drive. This should take approximately 1-2 minutes to download the ~1GB of data,

# In[3]:


file_ids = [
    '1KVDlXdaqPCi_gomHn6QK76NYmF7RONT5',
    '1MRUr0wAhVFIdrFjUg0YcDjaKaABKs09o',
    '1gw7BG7wlZfU_BGNDf3QH1UARWUjwsIjz',
    '1Cwl0gLUHFTwkvdSEm1KZiAsUU7jZyrAN'
]
paths = [
    'aviso.zip',
    'chla.zip',
    'sst.zip',
    'winds.zip'
]

for f, p in zip(file_ids, paths):
    gdd.download_file_from_google_drive(file_id=f,
                                    dest_path='./data/'+ p,
                                    showsize=True,
                                    unzip=True)


# ## Data import
# 
# Here we pull all data into the notebook from local storage and do some preliminary visualization and data prep.
# 
# The chlorophyll-a data, altimetry, and wind data cover 2010 through the end of 2019. The sea surface temp dataset covers 2017 through 2019.

# ### Chlorophyll-a
# from the Ocean Colour Climate Change Initiative project. This is 5 day averaged product generated by merging data from many sources including MODIS Aqua, Sentinel-3, MERIS and more. 

# In[4]:


chla_ds = xr.open_zarr('data/chla_long.zarr')
chla_ds


# Check out how the data is chunked and its dimensions.

# In[5]:


chla_ds.chlor_a


# Visualize the first time step of the data

# In[6]:


fig,ax = plt.subplots(figsize=(12,8))
chla_ds = chla_ds.sel(lat=slice(44,26),lon=slice(-82,-66))
chla_ds.chlor_a[0].plot(ax=ax, cmap='jet', norm=LogNorm(vmin=0.01, vmax=10))


# And now use hvplot to make the plotting interactive and add the ability to move through each time step with the scroll bar on the right.

# In[7]:


proj = crs.Orthographic(-90, 30)

chla_ds.chlor_a.hvplot.quadmesh(
    'lon', 'lat', projection=proj, project=True,
    cmap='jet', dynamic=True, coastline='10m', 
    frame_width=300, logz=True, clim=(0.01,20), rasterize=True)


# ## Sea Surface Height 
# via AVISO which uses multiple altimeters to generate a .25 degree daily SSH and SLA product.
# 
# Note that the data variable we're using in this `xarray` dataset is called adt which stands for absolute dynamic topography and is our SSH variable.

# In[8]:


ssh_ds = xr.open_zarr('data/aviso.zarr')
ssh_ds


# Again inspect chunking and data details. 

# In[9]:


ssh_ds.adt


# Visualize one time step of the absolute dynamic topography (ADT) variable in this dataset.

# In[10]:


fig, ax = plt.subplots(figsize=(12,9), subplot_kw=dict(projection=crs.PlateCarree()))
ssh_ds.adt[0].plot(ax=ax, cmap='bwr')
ax.coastlines(resolution='10m')
ax.set_ylim(26,44)
ax.set_xlim(-82,-66)


# Now use hvplot to allow interactive visualization of this dataset

# In[11]:


proj = crs.Orthographic(-90, 30)
ssh_ds.adt.hvplot.quadmesh(
    'longitude', 'latitude', projection=proj, project=True,
    cmap='bwr', dynamic=True, coastline='10m', 
    frame_width=300, clim=(-1,1), rasterize=True)


# We'll be using contours below (aka isolines in the sea surface height) to determine the location of the Gulf Stream front so let's see what that looks like now

# In[12]:


fig, ax = plt.subplots(figsize=(12,9), subplot_kw=dict(projection=crs.PlateCarree()))
ssh_ds.adt[0].plot.contourf(ax=ax, vmin=-1,vmax=1,center=0.25)
ax.coastlines(resolution='10m')
ax.set_ylim(26,44)
ax.set_xlim(-82,-66)


# ## Sea Surface Temperature
# via NOAA and the Group for High Resolution Sea Surface Temperature (GHRSST).

# In[13]:


sst_ds = xr.open_zarr('data/sst.zarr')
sst_ds


# In[14]:


sst_ds.analysed_sst


# Let's take a quick look at a static and interactive plot of SST

# In[15]:


fig,ax = plt.subplots(figsize=(12,8))
sst_ds.analysed_sst[0].plot(ax=ax, vmin=278, vmax=298, cmap='inferno')


# In[16]:


proj = crs.Orthographic(-90, 30)

sst_ds.analysed_sst.hvplot.quadmesh(
    'lon', 'lat', projection=proj, project=True,
    cmap='inferno', dynamic=True, coastline='10m', 
    frame_width=500, rasterize=True)


# Now let's take the mean of a year of daily SST data just to clearly see the thermal footprint of the Gulf Stream

# In[17]:


fig,ax = plt.subplots(figsize=(12,8))
sst_ds.analysed_sst[0:365].mean(dim='time', skipna=True).plot(ax=ax, vmin=284, vmax=300, cmap='inferno')
ax.set_title('Mean SST of the Gulf Stream')


# ## Sea Wind Speed
# As our last dataset we'll pull in data from EUMETSAT's Metop-A satellite. ASCAT is a microwave scatterometer designed to measure surface winds over the global ocean and these are 8 day composites of wind speed.

# In[18]:


wind_ds = xr.open_zarr('data/winds.zarr')
wind_ds


# Note above that the longitude goes from 0 to 360 as opposed to -180 through 180 as was the case with the other datasets. We'll process it here to match the others.

# In[19]:


wind_ds['longitude'] = wind_ds['longitude'] -360


# This dataset comes with x and y wind speeds, let's combine them by adding the absolute value of each just to get an idea of total wind speed in a given location.

# In[20]:


wind_ds['wind_total'] = abs(wind_ds.y_wind) + abs(wind_ds.y_wind)


# Let's take a look at the average wind speed across this region.

# In[21]:


fig, ax = plt.subplots(figsize=(12,9), subplot_kw=dict(projection=crs.PlateCarree()))
wind_ds.wind_total[0].plot(ax=ax, cmap='coolwarm')
ax.coastlines(resolution='10m')
ax.set_ylim(26,44)
ax.set_xlim(-82,-66)


# This data is relatively noisy and satellite tracks are somewhat visible, but we'll resample to monthly below so the noise isn't as much of an issue. As a test let's look at the average wind in each month for the entire dataset.

# In[22]:


wind_ds.groupby('time.month').median(dim='time').wind_total.hvplot.quadmesh(
    'longitude', 'latitude', projection=proj, project=True,
    cmap='bwr', dynamic=True, coastline='10m', 
    frame_width=300, clim=(0,12), rasterize=True)


# ## Gulf Stream front location
# 
# Now as a final step in our data prep before we get into the bulk of the analysis we'll find the location of the front at each time step by finding the 0.25 contour line in sea surface height as done in [Andres 2016](http://dx.doi.org/10.1002/2016GL069966).

# As a quick overview we're going to:
# 1. Get the monthly Gulf Stream front location via the SSH data based on the 25cm contour.
# 2. Turn those contours into a geodataframe of monthly front lines
# 3. Use these front lines to grab data from all our xarray DataSets to create monthly along-front DataSets

# First we'll resample our SSH data from daily to monthly time steps.

# In[23]:


ssh_ds_monthly = ssh_ds.resample(time="1M").mean()


# Then we'll run the subpixel_contours() function which takes in an xarray DataArray and related geospatial metadata and returns a geodataframe of the contours at each time step.
# 
# To do this we need to define an affine transformation for the SSH data since it doesn't exist in the dataset currently. Learn more [here](https://www.perrygeo.com/python-affine-transforms.html). The basic format of an affine transformation is (a, b, c, d, e, f):
# 
#     a = width of a pixel
#     b = row rotation (typically zero)
#     c = x-coordinate of the upper-left corner of the upper-left pixel
#     d = column rotation (typically zero)
#     e = height of a pixel (typically negative)
#     f = y-coordinate of the of the upper-left corner of the upper-left pixel

# We define the CRS as 4326 which is WGS84 so it just uses lat and lon. Thus the width of a pixel is .25 degrees.

# In[24]:


ssh_affine = Affine(0.25, 0.0, -81.875-.25/2, 0.0, 0.25, 26.125-.25/2)
ssh_affine


# In[25]:


# this function needs all the data in memory so load it in
ssh_ds_monthly.adt.load()

gdf = subpixel_contours(ssh_ds_monthly.adt, [0.25], crs='EPSG:4326', min_vertices=50, affine=ssh_affine, verbose=True)


# In[26]:


gdf.plot()
plt.title('All front lines from subpixel_contours()')
plt.xlabel('longitude')
plt.ylabel('latitude')


# As can be seen above, this approach to finding the front does pick up a few anomalous frontal eddies in the South Atlantic Bight so we're going to require all lines pass south of 28 degrees latitude and east of -70 degrees longitude.

# In[27]:


print(len(gdf))
gdf = gdf.cx[:,26:28].cx[-70:-66,:]
print(len(gdf))
gdf.plot()
plt.title('Filtered front lines')
plt.xlabel('longitude')
plt.ylabel('latitude')


# Note that we lose 13 months out of the 10 years of data because of this constraint.
# 
# And finally because the Gulf Stream is so close to the Florida coast around 26 degrees north, we'll clip it to 28 N and above and to west of -66.5 longitude so it isn't on the edge of our data.

# In[28]:


# Create a custom polygon to clip with
polygon = Polygon([(-82, 28.2), (-82, 44), (-66.5, 44), (-66.5, 28.2), (-82, 28.2)])
poly_gdf = gpd.GeoDataFrame([1], geometry=[polygon], crs=gdf.crs)
poly_gdf.plot()
plt.title('clip to this area')
plt.xlabel('longitude')
plt.ylabel('latitude')


# In[29]:


gdf = gpd.clip(gdf, polygon)
gdf.plot()
plt.title('Filtered and clipped front lines')
plt.xlabel('longitude')
plt.ylabel('latitude')


# That seems to nicely solve the issue

# In[30]:


gdf.hvplot.paths('Longitude', 'Latitude', geo=True, color='red', alpha=0.8,
                       tiles='ESRI', frame_width=500, title='filtered and clipped front lines')


# Let's quickly check to see if our contour matches exactly with the contours from the xarray contour function.
# 

# In[31]:


fig,ax= plt.subplots(1,2, figsize=(16,8))
ssh_ds_monthly.adt[0].plot.contourf(ax=ax[0], vmin=-1,vmax=1,center=0.25,levels=5)
ssh_ds_monthly.adt[0].plot.contourf(ax=ax[1], vmin=-1,vmax=1,center=0.25,levels=5)
gdf.head(1).plot(ax=ax[1], color='black')
ax[0].set_title('contours from xarray')
ax[1].set_title('adding in our newly calculated contour')


# ## Selecting data along the front 

# Now let's turn our contour lines into points so that we can grab the xarray values at each point.

# In[32]:


lon,lat = gdf.iloc[0].geometry.coords.xy
df = pd.DataFrame(list(zip(lat,lon)), columns=['lat', 'lon'])

points_gdf = gpd.GeoDataFrame(
    df, geometry=gpd.points_from_xy(df.lon, df.lat))


# In[33]:


points_gdf.hvplot.points(geo=True, color='red', alpha=0.8,
                       tiles='ESRI', frame_width=500, title='points from a contour line')


# These points are irregularly spaced as a product of the contour algorithm and we want to have them all regularly spaced along the line. So we'll make 90 evenly spaced points on each front line.
# 
# We also want a coastal line and an open ocean line that run parallel to the front line but allow us to compare differences along the front and nearby waters. We'll offset lines ~20 km from our front line on either side.

# In[34]:


points_list = []
points_list_coastal = []
points_list_sargasso = []

# number of points along the line
n = 90
for i, line in gdf.iterrows():
    # create equally spaced points along the line
    distances = np.linspace(0, line.geometry.length, n)
    # turn them into shapely points
    points = [line.geometry.interpolate(distance) for distance in distances]
    # join them into multipoint objects
    multipoint = unary_union(points)
    # move north and west
    coastal_points = []
    for point in points:
        coastal_points.append(transform(lambda x, y: (x-0.2, y+0.2), point))
    multipoint_coastal = unary_union(coastal_points)
    # move south and east
    sargasso_points = []
    for point in points:
        sargasso_points.append(transform(lambda x, y: (x+0.2, y-0.2), point))
    multipoint_sargasso = unary_union(sargasso_points)
    
    points_list.append(multipoint)
    points_list_coastal.append(multipoint_coastal)
    points_list_sargasso.append(multipoint_sargasso)


# Turn those points into a geodataframe

# In[35]:


points_gdf = gpd.GeoDataFrame(gdf['time'], geometry=points_list)
points_gdf.head()


# In[36]:


coastal_points_gdf = gpd.GeoDataFrame(gdf['time'], geometry=points_list_coastal)
coastal_points_gdf.head()


# In[37]:


sargasso_points_gdf = gpd.GeoDataFrame(gdf['time'], geometry=points_list_sargasso)
sargasso_points_gdf.head()


# Display 10 frontal lines against the backdrop of SSH

# In[38]:


ssh_ds_monthly.adt[0].hvplot.quadmesh(
    'longitude', 'latitude',
    cmap='bwr', dynamic=True, coastline='10m', 
    frame_width=300, clim=(-1,1), datashade=False, title='evenly spaced points from front contours') * \
points_gdf.head(10).hvplot.points(geo=True, color='black', alpha=0.7)


# Now let's check out the coastal, frontal, and sargasso lines together

# In[39]:


ssh_ds_monthly.adt[0].hvplot.quadmesh(
    'longitude', 'latitude',
    cmap='bwr', dynamic=True, coastline='10m', 
    frame_width=300, clim=(-1,1), datashade=False, title='points from coastal, front, and sargasso lines') * \
coastal_points_gdf.head(1).hvplot.points(geo=True, color='green', alpha=0.7, label='coastal') * \
sargasso_points_gdf.head(1).hvplot.points(geo=True, color='blue', alpha=0.7, label='sargasso') * \
points_gdf.head(1).hvplot.points(geo=True, color='black', alpha=0.7, label='front')


# Now let's run this same process but for each time step.
# 
# Create an xarray DataArray each day and each frontal line

# In[40]:


# there could be a better way to do this than iterating through the gdf
x_list = []
y_list = []
for i in range(len(points_gdf)):
    points = np.array([[p.x, p.y] for p in points_gdf.iloc[i].geometry])
    x_list.append(points[:,0])
    y_list.append(points[:,1])
x_list = np.array(x_list)
y_list = np.array(y_list)
x_list.shape, y_list.shape


# Using those points and time steps let's create a few DataArrays that we'll used to select from the primary DataSets

# In[41]:


target_lon = xr.DataArray(x_list, dims=["time", "points"])
target_lat = xr.DataArray(y_list, dims=["time", "points"])
target_time = xr.DataArray(pd.to_datetime(points_gdf['time']).values, dims="time")
target_lat


# Now let's resample all the datasets to monthly averages

# In[42]:


chla_ds_monthly = chla_ds.resample(time="1M").mean()
sst_ds_monthly = sst_ds.resample(time="1M").mean()
wind_ds_monthly = wind_ds.resample(time="1M").mean()


# And now run select the nearest neighbor in time and space to each point

# In[43]:


ssh_ds_frontal = ssh_ds_monthly.adt.sel(time=target_time, latitude=target_lat, longitude=target_lon, method="nearest")
chla_ds_frontal = chla_ds_monthly.chlor_a.sel(time=target_time, lat=target_lat, lon=target_lon, method="nearest")
sst_ds_frontal = sst_ds_monthly.analysed_sst.sel(time=target_time, lat=target_lat, lon=target_lon, method="nearest")
wind_ds_frontal = wind_ds_monthly.wind_total.sel(time=target_time, latitude=target_lat, longitude=target_lon, method="nearest")


# As a sanity check, let's take a look at the SSH along the points in this front line which was derived from the 0.25 contour line just to ensure that is is approximately 0.25 still.

# In[44]:


ssh_ds_frontal[0].plot()
plt.title('SSH at points processed from 0.25m contour line')


# Note that this does introduce a little error since we've interpolated the contour and it isn't exactly on the .25m line

# In[45]:


chla_ds_frontal


# In[46]:


chla_ds_frontal[0].plot()
plt.title('chla along the front')
plt.xlabel('points from up to downstream')


# In[47]:


sst_ds_frontal[0].plot()
plt.title('SST along the front')
plt.xlabel('points from up to downstream')


# In[48]:


wind_ds_frontal[0].plot()
plt.title('wind along the front')
plt.xlabel('points from up to downstream')


# And finally we want to do the exact same thing for the coastal and sargasso lines. We'll write a quick function for this and apply it to the coastal and sargasso lines.

# In[49]:


def frontal_dataset_generation(line_gdf, ssh_monthly, chla_monthly, sst_monthly, wind_monthly):
    x_list = []
    y_list = []
    for i in range(len(line_gdf)):
        points = np.array([[p.x, p.y] for p in line_gdf.iloc[i].geometry])
        x_list.append(points[:,0])
        y_list.append(points[:,1])
    x_list = np.array(x_list)
    y_list = np.array(y_list)

    target_lon = xr.DataArray(x_list, dims=["time", "points"])
    target_lat = xr.DataArray(y_list, dims=["time", "points"])
    target_time = xr.DataArray(pd.to_datetime(points_gdf['time']).values, dims="time")

    ssh_line = ssh_monthly.adt.sel(time=target_time, latitude=target_lat, longitude=target_lon, method="nearest")
    chla_line = chla_monthly.chlor_a.sel(time=target_time, lat=target_lat, lon=target_lon, method="nearest")
    sst_line = sst_monthly.analysed_sst.sel(time=target_time, lat=target_lat, lon=target_lon, method="nearest")
    wind_line = wind_monthly.wind_total.sel(time=target_time, latitude=target_lat, longitude=target_lon, method="nearest")
    
    return(ssh_line,chla_line,sst_line,wind_line)


# In[50]:


ssh_ds_sargasso,chla_ds_sargasso,sst_ds_sargasso,wind_ds_sargasso = frontal_dataset_generation(sargasso_points_gdf, ssh_ds_monthly, chla_ds_monthly, sst_ds_monthly, wind_ds_monthly)

ssh_ds_coastal,chla_ds_coastal,sst_ds_coastal,wind_ds_coastal = frontal_dataset_generation(coastal_points_gdf, ssh_ds_monthly, chla_ds_monthly, sst_ds_monthly, wind_ds_monthly)


# **Now our data is ready let's analyze!**
# 
# Let's just quickly look at a single time step from each dataset. We'll represent the coastal line with green, frontal line with black, and open ocean (sargasso sea) line with blue.

# In[51]:


fig, ax = plt.subplots()

sst_ds_coastal[0].plot(ax=ax, color='green', label='coastal')
sst_ds_frontal[0].plot(ax=ax, color='black', label='frontal')
sst_ds_sargasso[0].plot(ax=ax, color='blue', label='sargasso')

ax.legend(loc='best')

plt.title('SST along the front')
plt.xlabel('points from up to downstream')


# In[52]:


fig, ax = plt.subplots()

chla_ds_coastal[0].plot(ax=ax, color='green', label='coastal')
chla_ds_frontal[0].plot(ax=ax, color='black', label='frontal')
chla_ds_sargasso[0].plot(ax=ax, color='blue', label='sargasso')

ax.legend(loc='best')

plt.title('Chla along the front')
plt.xlabel('points from up to downstream')


# In[53]:


fig, ax = plt.subplots()

wind_ds_coastal[0].plot(ax=ax, color='green', label='coastal')
wind_ds_frontal[0].plot(ax=ax, color='black', label='frontal')
wind_ds_sargasso[0].plot(ax=ax, color='blue', label='sargasso')

ax.legend(loc='best')

plt.title('Wind along the front')
plt.xlabel('points from up to downstream')


# # Data processing and analysis
# 
# **Now we've made it through all the data processing steps and can begin analysis of this frontal dataset we've created.**

# ## Visualizing all datasets

# First let's just examine the data quickly to ensure there isn't anything obviously erroneous and to see if we find any clear trends. These charts show chla, SST, and SSH with the front plotted over top of them for all of 2018.

# In[54]:


# add a datetime object to the front gdf so that we can select by time
points_gdf['dt_time'] = pd.to_datetime(points_gdf['time'])


# In[55]:


for i in range(12,24):
    fig,ax = plt.subplots(1,4,figsize=(30,8))
    # there must be a better way to convert from cftime to datetime in xarray
    current_dt = datetime.datetime.strptime(str(sst_ds_monthly.time[i].dt.strftime("%Y, %m, %d").values),"%Y, %m, %d")
    
    sst_ds_monthly.analysed_sst[i].plot(ax=ax[1], cmap='inferno', vmin=290, vmax=304)
    # have to select from these based on SST because they are much longer time series
    chla_ds_monthly.chlor_a.sel(time=current_dt, method="nearest").plot(ax=ax[0], cmap='jet', norm=LogNorm(vmin=0.01, vmax=10))
    ssh_ds_monthly.adt.sel(time=current_dt, method="nearest").plot(ax=ax[2], cmap='bwr', vmin=-1, vmax=1)
    wind_ds_monthly.wind_total.sel(time=current_dt, method="nearest").plot(ax=ax[3], cmap='coolwarm', vmin=0, vmax=9)
    month_gdf = points_gdf[points_gdf['dt_time'] == current_dt]
    month_gdf.plot(ax=ax[0], color='black', alpha=0.9, markersize=10)
    month_gdf.plot(ax=ax[1], color='black', alpha=0.9, markersize=10)
    month_gdf.plot(ax=ax[2], color='black', alpha=0.9, markersize=10)
    month_gdf.plot(ax=ax[3], color='black', alpha=0.9, markersize=10)
    plt.show()


# **These figures really make it clear how much the GS is major boundary between the coastal and the open ocean ecosystems.**
# 
# A large increase in chlorophyll-a beyond the Gulf Stream front, aka in the Sargasso Sea, is evident from late fall to early spring. The spring bloom in the North Atlantic is also clear where it peaks in April and May. Note that the scale is logarithmic so the red colors are two orders of magnitude highers than teal. 
# 
# The SST data matches what one might expect, decreasing throughout the winter and increasing in summer with the Gulf Stream itself staying as the warmest water throughout all seasons. We can also see a few eddies in SST and SSH even though this data is resampled to monthly time steps.

# ## Seasonality along the Gulf Stream

# Now let's look at the chla on the front itself over time. Make sure you use the scroll bar to explore the data over time.

# In[56]:


chla_ds_frontal.hvplot.line(x='points', value_label='time', ylim=(0,1), color='black', label='frontal') * chla_ds_coastal.hvplot.line(x='points', value_label='time', ylim=(0,1), color='green', label='coastal') * chla_ds_sargasso.hvplot.line(x='points', value_label='time', ylim=(0,1), color='blue', label='sargasso', xlabel='points up to downstream')


# **Note sustained spike in chla off Cape Hatteras near the front on the coastal side.**
# 
# The spring bloom is also noticable in all three lines from points 60 to 90.
# 
# Moving on here to SST:

# In[57]:


sst_ds_frontal.hvplot.line(x='points', value_label='time', ylim=(290,304), color='black', label='frontal') * sst_ds_coastal.hvplot.line(x='points', value_label='time', ylim=(290,304), color='green', label='coastal') *sst_ds_sargasso.hvplot.line(x='points', value_label='time', ylim=(290,304), color='blue', label='sargasso', xlabel='points up to downstream')


# SST follows a relatively predictable pattern with most variation seemingly tied to seasonality.
# 
# Moving to wind now:

# In[58]:


wind_ds_coastal.hvplot.line(x='points', value_label='time', ylim=(0,15), xlim=(1,90), color='green', label='coastal') * wind_ds_frontal.hvplot.line(x='points', value_label='time', ylim=(0,15), xlim=(1,90),  color='black', label='frontal') * wind_ds_sargasso.hvplot.line(x='points', value_label='time', ylim=(0,15), xlim=(1,90), color='blue', label='sargasso', xlabel='points up to downstream')


# Winds are highly variable, but show some seasonality with an increase in the winter. This is a part of most descriptions of the spring bloom, winter storms drive mixing which adds sufficient nutrients into the euphotic zone for the major spring increase once sunlight/stratification returns.

# You can see some hints at the dynamics at play in these plots, but let's look at it all at once via a Hovmöller diagram with the spatial points on the x axis, time on the y, and the color showing the value of the variable.
# 
# In these plots Point 0 is the furthest point upstream around Florida and Point 90 is downstream in the North Atlantic.

# In[59]:


fig, ax = plt.subplots(1,3,figsize=(20,8))

# cut the last two dates off so it doesn't strech it to 2020
chla_ds_coastal[:-2].plot(ax=ax[0],norm=LogNorm(vmin=0.05, vmax=5), cmap='jet', xlim=(2,90))
chla_ds_frontal[:-2].plot(ax=ax[1],norm=LogNorm(vmin=0.05, vmax=5), cmap='jet', xlim=(2,90))
chla_ds_sargasso[:-2].plot(ax=ax[2],norm=LogNorm(vmin=0.05, vmax=5), cmap='jet', xlim=(2,90))

ax[0].set_title('coastal chla')
ax[0].set_xlabel('points up to downstream')
ax[1].set_title('frontal chla')
ax[1].set_xlabel('points up to downstream')
ax[2].set_title('sargasso chla')
ax[2].set_xlabel('points up to downstream')

plt.show()


# Let's look one more time at chla over just a few years

# In[60]:


fig, ax = plt.subplots(1,3,figsize=(20,8))

# cut the last two dates off so it doesn't strech it to 2020
chla_ds_coastal[:36].plot(ax=ax[0],norm=LogNorm(vmin=0.05, vmax=5), cmap='jet', xlim=(2,90))
chla_ds_frontal[:36].plot(ax=ax[1],norm=LogNorm(vmin=0.05, vmax=5), cmap='jet', xlim=(2,90))
chla_ds_sargasso[:36].plot(ax=ax[2],norm=LogNorm(vmin=0.05, vmax=5), cmap='jet', xlim=(2,90))

ax[0].set_title('coastal chla')
ax[0].set_xlabel('points up to downstream')
ax[1].set_title('frontal chla')
ax[1].set_xlabel('points up to downstream')
ax[2].set_title('sargasso chla')
ax[2].set_xlabel('points up to downstream')

plt.show()


# In[61]:


fig, ax = plt.subplots(1,3,figsize=(20,5))

sst_ds_coastal.plot(ax=ax[0], vmin=293, vmax=303, cmap='inferno', xlim=(2,90))
sst_ds_frontal.plot(ax=ax[1],vmin=293, vmax=303, cmap='inferno', xlim=(2,90))
sst_ds_sargasso.plot(ax=ax[2],vmin=293, vmax=303, cmap='inferno', xlim=(2,90))

ax[0].set_title('coastal SST')
ax[0].set_xlabel('points up to downstream')
ax[1].set_title('frontal SST')
ax[1].set_xlabel('points up to downstream')
ax[2].set_title('sargasso SST')
ax[2].set_xlabel('points up to downstream')

plt.show()


# In[62]:


fig, ax = plt.subplots(1,3,figsize=(20,8))

wind_ds_coastal.plot(ax=ax[0], vmin=3, vmax=12, cmap='coolwarm', xlim=(2,90))
wind_ds_frontal.plot(ax=ax[1],vmin=3, vmax=12, cmap='coolwarm', xlim=(2,90))
wind_ds_sargasso.plot(ax=ax[2],vmin=3, vmax=12, cmap='coolwarm', xlim=(2,90))

ax[0].set_title('coastal wind')
ax[0].set_xlabel('points up to downstream')
ax[1].set_title('frontal wind')
ax[1].set_xlabel('points up to downstream')
ax[2].set_title('sargasso wind')
ax[2].set_xlabel('points up to downstream')

plt.show()


# The intense seasonality of SST is clearly apparent and a dominant factor but chla is much more variable, influenced by more than just the temperature and light exposure.
# 
# We can show more quantatively that chla is more variable via the coefficient of variation here. 
# 
# **For reference point 20 is approx Charleston SC, point 40 is approx Cape Lookout NC, point 50 is Cape Hatteras NC, and point 60 and beyond are offshore.**

# In[63]:


(chla_ds_coastal.std(dim='time')/chla_ds_frontal.mean(dim='time')).plot(xlim=(2,90), color='green',label='coastal',)
(chla_ds_frontal.std(dim='time')/chla_ds_frontal.mean(dim='time')).plot(xlim=(2,90), color='black',label='frontal',)
(chla_ds_sargasso.std(dim='time')/chla_ds_frontal.mean(dim='time')).plot(xlim=(2,90), color='blue', label='sargasso')

plt.title('coefficient of variation of chla at each point')
plt.legend()
plt.xlabel('points up to downstream')
plt.ylim(0,3)


# You can see that right around point 45 which is nearshore and near Cape Hatteras is where the physical structure of GS connects with an area that has major bursts of productivity.

# In[64]:


(sst_ds_coastal.std(dim='time')/sst_ds_frontal.mean(dim='time')).plot(xlim=(2,90), color='green',label='coastal')
(sst_ds_frontal.std(dim='time')/sst_ds_frontal.mean(dim='time')).plot(xlim=(2,90), color='black',label='frontal')
(sst_ds_sargasso.std(dim='time')/sst_ds_frontal.mean(dim='time')).plot(xlim=(2,90), color='blue',label='sargasso')

plt.title('coefficient of variation of SST at each point')
plt.legend()
plt.xlabel('points up to downstream')
#plt.ylim(0,3)


# As you can see chla is 50-100 times more variable than SST and the temperature of the front is much more variable downstream which makes sense given the more intense seasonality there.
# 
# ## The full time series
# 
# Now let's look at the full time series of chla in the various segments. First we'll look at the offshore segment in the North Atlantic.

# In[65]:


fig,ax = plt.subplots(figsize=(12,7))
chla_ds_frontal[:,60:].median('points').plot(ax=ax,color='black', label='frontal')
chla_ds_coastal[:,60:].median('points').plot(ax=ax,color='green', label='coastal')
chla_ds_sargasso[:,60:].median('points').plot(ax=ax,color='blue', label='sargasso')
ax.set_ylim(0,1.2)

ax.set_title('avg chla of offshore section over time')
ax.legend()


# And then the section around Cape Lookout and Cape Hatteras

# In[66]:


fig,ax = plt.subplots(figsize=(12,7))
chla_ds_frontal[:,40:50].median('points').plot(ax=ax,color='black', label='frontal')
chla_ds_coastal[:,40:50].median('points').plot(ax=ax,color='green', label='coastal')
chla_ds_sargasso[:,40:50].median('points').plot(ax=ax,color='blue', label='sargasso')
ax.set_ylim(0,1.2)

ax.set_title('avg chla of Cape Hatteras section over time')
ax.legend()


# In[67]:


fig,ax = plt.subplots(figsize=(12,7))
chla_ds_frontal[:,20:40].median('points').plot(ax=ax,color='black', label='frontal')
chla_ds_coastal[:,20:40].median('points').plot(ax=ax,color='green', label='coastal')
chla_ds_sargasso[:,20:40].median('points').plot(ax=ax,color='blue', label='sargasso')
ax.set_ylim(0,1.2)

ax.set_title('avg chla of South Atlantic Bight section over time')
ax.legend()


# **A couple major takeaways from this analysis thus far.**
# 
# The offshore region further downstream seems to have a major increase in the spring - the well studied spring bloom - and the region further south (upstream) around the South Atlantic Bight (SAB) has a relatively large increase throughout the whole winter starting mid fall and staying high until the spring. The region around Cape Hatteras falls in the middle of these two patterns but with generally higher chla, particularly on the coastal side, which makes sense given its proximity to the coast in this region and generally large outflows from Pamlico Sound and Chesapeake Bay injecting nutrients into the water. It is also worth noting that the more northerly section of the front does seem to have a more notable dip mid-winter, possibly when light begins to be limiting, and this is not as pronounced in the southerly section of the front. Though the light difference is not extreme, it is 9.5 hours on the winter solstice at 37 degrees N and 10.6 at 27 degrees N.

# In[68]:


fig,ax = plt.subplots(figsize=(12,7))
chla_ds_frontal.median('points').plot(ax=ax,color='black', label='frontal')
chla_ds_coastal.median('points').plot(ax=ax,color='green', label='coastal')
chla_ds_sargasso.median('points').plot(ax=ax,color='blue', label='sargasso')
ax.set_ylim(0,0.6)

ax.set_title('avg chla of entire line over time')
ax.legend()


# It is interesting to see that for nearly all of the Gulf Stream from from Florida out into the North Atlantic there is a fall bloom that is comparable or greater than the spring bloom and this is true of both sides of the front. It decreases around mid-winter, likely due to light limitation (depending on the region), and then increases again as light begins to increase. At least that is our general idea!

# **Inspecting the SST data here in the same fashion to see if there are any surprising trends.**

# In[69]:


fig,ax = plt.subplots(figsize=(12,7))
sst_ds_frontal[:,60:].mean('points').plot(ax=ax,color='black', label='frontal')
sst_ds_coastal[:,60:].mean('points').plot(ax=ax,color='green', label='coastal')
sst_ds_sargasso[:,60:].mean('points').plot(ax=ax,color='blue', label='sargasso')
ax.set_ylim(290,303)

ax.set_title('avg SST of offshore section over time')
ax.legend()


# In[70]:


fig,ax = plt.subplots(figsize=(12,7))
sst_ds_frontal[:,40:50].mean('points').plot(ax=ax,color='black', label='frontal')
sst_ds_coastal[:,40:50].mean('points').plot(ax=ax,color='green', label='coastal')
sst_ds_sargasso[:,40:50].mean('points').plot(ax=ax,color='blue', label='sargasso')
ax.set_ylim(290,303)

ax.set_title('avg SST of Cape Hatteras section over time')
ax.legend()


# In[71]:


fig,ax = plt.subplots(figsize=(12,7))
sst_ds_frontal[:,20:40].mean('points').plot(ax=ax,color='black', label='frontal')
sst_ds_coastal[:,20:40].mean('points').plot(ax=ax,color='green', label='coastal')
sst_ds_sargasso[:,20:40].mean('points').plot(ax=ax,color='blue', label='sargasso')
ax.set_ylim(290,303)

ax.set_title('avg SST of South Atlantic Bight section over time')
ax.legend()


# **Finally let's see how the wind behaves in relation to chla along the front**

# In[72]:


fig,ax = plt.subplots(figsize=(12,7))
wind_ds_frontal[:,:,60:].median('points').plot(ax=ax,color='black', label='frontal')
wind_ds_coastal[:,:,60:].median('points').plot(ax=ax,color='green', label='coastal')
wind_ds_sargasso[:,:,60:].median('points').plot(ax=ax,color='blue', label='sargasso')
ax.set_ylim(0,10)

ax.set_title('avg wind of offshore section over time')
ax.legend()


# In[73]:


fig,ax = plt.subplots(figsize=(12,7))
wind_ds_frontal[:,:,60:].median('points').groupby('time.month').median(dim='time').plot(ax=ax,color='black', label='frontal')
wind_ds_coastal[:,:,60:].median('points').groupby('time.month').median(dim='time').plot(ax=ax,color='green', label='coastal')
wind_ds_sargasso[:,:,60:].median('points').groupby('time.month').median(dim='time').plot(ax=ax,color='blue', label='sargasso')
ax.set_ylim(0,10)

ax.set_title('avg monthly wind of offshore section')
ax.legend()


# In[74]:


fig,ax = plt.subplots(figsize=(12,7))
wind_ds_frontal[:,:,40:50].median('points').groupby('time.month').median(dim='time').plot(ax=ax,color='black', label='frontal')
wind_ds_coastal[:,:,40:50].median('points').groupby('time.month').median(dim='time').plot(ax=ax,color='green', label='coastal')
wind_ds_sargasso[:,:,40:50].median('points').groupby('time.month').median(dim='time').plot(ax=ax,color='blue', label='sargasso')
ax.set_ylim(0,10)

ax.set_title('avg monthly wind of Cape Hatteras section')
ax.legend()


# In[75]:


fig,ax = plt.subplots(figsize=(12,7))
wind_ds_frontal[:,:,20:40].median('points').groupby('time.month').median(dim='time').plot(ax=ax,color='black', label='frontal')
wind_ds_coastal[:,:,20:40].median('points').groupby('time.month').median(dim='time').plot(ax=ax,color='green', label='coastal')
wind_ds_sargasso[:,:,20:40].median('points').groupby('time.month').median(dim='time').plot(ax=ax,color='blue', label='sargasso')
ax.set_ylim(0,10)

ax.set_title('avg monthly wind of South Atlantic Bight section')
ax.legend()


# All three sections generally show a peak in early summer and in the fall but the patterns aren't 100% clear and the data itself may be a bit suspect. 
# 
# # Suggested next steps
# 
#  * This work both agrees with many previous observations and modeling studies on the phenology of the Atlantic and also seems to differ in the exact timing and intensity of chla. It would be worthwhile to compare these findings more directly to data further from the front in the Sargasso Sea and the North Atlantic to try to differentiate patterns.
#  * Longer time series of chla and SSH would allow for an investigation into impacts of decadal cycles such as the North Atlantic Oscillation which may play a role of similar magnitude to the front dynamics over large scales.
#  * The intense spikes of productivity off Cape Hatter and Cape Lookout were intriguing and warrant further analysis. Pairing this productivity data with observations of higher trophic levels and top predators may help illuminate why this area is such a hot spot for biodiversity.
#  * An improved wind product with higher spatial resolution could improve future work and allow investigation to look more directly at the impact of wind on chla. This would be useful to do on a time lagged basis.

# # References
# 
#  * Andres M. (2016) On the recent destabilization of the Gulf Stream path downstream of Cape Hatteras. Geo Res Letters. https://doi.org/10.1002/2016GL069966
#  * Boss E., Behrenfeld M. (2010) In situ evaluation of the initiation of the North Atlantic phytoplankton bloom. Geo Res Letters. https://doi.org/10.1029/2010GL044174.
#  * Fischer A.D., Moberg E.A., Alexander H., Brownlee E.F., Hunter-Cevera K.R., Pitz K.J., Rosengard S.Z., Sosik H.M. (2014) Sixty Years of Sverdrup: A Retrospective of Progress in the Study of Phytoplankton Blooms. Oceanography. https://doi.org/10.5670/oceanog.2014.26
#  * Henson S.A., Dunne J.P., Sarmiento J.L. (2009) Decadal variability in North Atlantic phytoplankton blooms. JGR: Oceans. https://doi.org/10.1029/2008JC005139
#  * Lévy, M., Franks, P.J.S. & Smith, K.S. (2018) The role of submesoscale currents in structuring marine ecosystems. Nat Commun. https://doi.org/10.1038/s41467-018-07059-3
#  * Pascual A., Ruiz S., Nardelli B.B., Guinehut S., Iudicone S., Tintoré J. (2014) Net primary production in the Gulf Stream sustained by quasi-geostrophic vertical exchanges. Geo Res Letters. https://doi.org/10.1002/2014GL062569
#  * Siegel D. A., Doney S. C., Yoder J. A. (2002) The North Atlantic Spring Phytoplankton Bloom and Sverdrup's Critical Depth Hypothesis. Science. https://doi.org/10.1126/science.1069174
#  * Ueyama R., Monger B.C. (2005) Wind-induced modulation of seasonal phytoplankton blooms in the North Atlantic derived from satellite observations. Limn and Oceano. https://doi.org/10.4319/lo.2005.50.6.1820

# In[ ]:




