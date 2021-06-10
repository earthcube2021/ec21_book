#!/usr/bin/env python
# coding: utf-8

# # Frequency-Domain Analysis of Large Datasets

# ## Author(s)
# **Paige Martin and Ryan Abernathey**
# 
# - Author1 = {"name": "Paige Martin", "affiliation": "Australian National University/Lamont-Doherty Earth Observatory", "email": "paigemar@umich.edu", "orcid": "0000-0003-3538-633X"}
# - Author2 = {"name": "Ryan Abernathey", "affiliation": "Lamont-Doherty Earth Observatory", "email": "rpa@ldeo.columbia.edu", "orcid": "0000-0001-5999-4917"}

# ## Purpose
# 
# Climate model datasets are typically stored as global snapshots, i.e. chunked in time rather than space. For many workflows, this chunking works well (e.g. computations across spatial domains at every point in time). However, this storage format can create serious challenges for processing long time series at each point in space, as is the case for frequency-domain analysis. For large datasets with frequent (e.g. daily) output, it is not feasible to process each spatial point as a single time series, even with the help of distributed computing such as Dask.
# 
# This notebook provides an example scientific workflow for performing frequency-domain analysis on large datasets. Specifically, this notebook presents a workflow for computing the power spectrum of sea surface temperature in the [Community Earth System Model](https://www.cesm.ucar.edu) (CESM). While we carry out computations on a specific model, the main goal of this notebook is not to interpret scientific results from our computations, but rather to provide a working example of frequency-domain analysis on large datasets that others could follow.
# 
# Because the goal of this notebook is to provide an example of a workflow that works on large datasets, we have chosen to use a large dataset (CESM) that is available on Pangeo Cloud. This notebook is therefore developed for a Jupyter Hub environment that can access data stored on the Pangeo Cloud (e.g. [PangeoBinder](https://binder.pangeo.io)).
# 
# 
# ## Technical contributions
# 
# - demonstrates how to quickly rechunk data, e.g. from chunks in time to chunks in space, using the package [Rechunker](https://rechunker.readthedocs.io/en/latest/)
# - demonstrates how to easily perform Fourier analysis using the package [xrft](https://xrft.readthedocs.io/en/latest/)
# - shows that all of these steps, with the use of [Xarray](http://xarray.pydata.org/en/stable/) and [Dask](https://dask.org), can be taken with large datasets
# 
# 
# ## Methodology
# 
# The notebook follows three main steps:
# 
# 1) Rechunk the data. We begin by rechunking the CESM sea surface temperature (SST) output from global, daily snapshots to chunks in space and 5-year chunks in time. This step is accomplished using the package `rechunker`. 
# 
# 2) Fourier analysis. Next we compute the power spectrum of SST using the package `xrft`, which nicely integrates with `Xarray` and `Dask`. Within `xrft`, we are also able to easily apply detrending and windowing functions to our data, and also account for the fact that our data are real (with no imaginary components).
# 
# 3) Visualize the data. Last, we average over various frequency bands to show the spatial distribution of the SST power spectrum as global maps.
# 
# Between each of these steps, we would typically write out the intermediate data. Specifically, we would write out the rechunked data, as well as the processed power spectra. Being able to write out data at intermediate steps is crucial to this workflow. However, due to the inability to write out data from Binder, we instead import previously rechunked data and use a spatial subset to compute and plot the power spectrum in this notebook. We still include all code necessary to run every step if a user wishes to run this notebook elsewhere that allows for data to be written out.
# 
# 
# ## Results
# This notebook presents a feasible example for performing frequency-domain analysis on large datasets. Specifically, this notebook demonstrates how to quickly rechunk ~500GB of data from chunking in time to chunking in space. It also demonstrates how to pair the `xrft` library with the `rechunker` library to perform frequency-domain spectral analysis (here power spectra) and obtain interpretable results. We finish with a few sample plots to round out the workflow. This notebook is meant to serve as an example for others who wish to perform similar types of analysis.
# 
# 
# ## Funding
# 
# - Award1 = {"agency": "Gordon and Betty Moore Foundation", "award_code": "", "award_URL": "https://www.moore.org"}
# 
# ## Keywords
# 
# keywords=["frequency-domain", "Pangeo", "spectral analysis", "rechunking", "cloud computing"]
# 
# ## Citation
# 
# Martin and Abernathey 2021. Frequency-Domain Analysis of Large Datasets. Accessed at https://github.com/paigem/EC2021_Martin_and_Abernathey.
# 
# ## Acknowledgements 
# 
# We thank the Pangeo community for developing and maintaining most of the packages used in this notebook. We also acknowledge Pangeo Cloud, which provides the computing power for this analysis. 

# # Setup
# 
# ## Library import

# In[ ]:


# Reading in data
import intake
import gcsfs
import zarr
import os

# Data manipulation
import xarray as xr
import dask.array as dsa
from rechunker import rechunk
import xrft

# Distributed computing
from dask_gateway import Gateway
from dask.distributed import Client

# Visualization
import matplotlib.pyplot as plt


# # Data import, processing, and analysis

# ## Step 1: Open and rechunk the original data
# 
# This dataset is stored on the [Pangeo Cloud](https://catalog.pangeo.io/) and contains daily output for 41 model years.

# ### Access the data

# In[ ]:


# Access CESM POP2 control output from Pangeo Cloud data catalog
cat = intake.open_catalog("https://raw.githubusercontent.com/pangeo-data/pangeo-datastore/master/intake-catalogs/ocean/CESM_POP.yaml")
item = cat['CESM_POP_hires_control']  # Specify CESM high resolution control run
ds_orig = item.to_dask().reset_coords(drop=True)  # drop unneeded coordinates for efficiency
ds_orig


# We are interested in working with SST.
# Let's make a quick plot to see what the data look like.
# Pre-coarsening makes the plot faster to load.

# In[ ]:


ds_orig.SST[0].coarsen(nlat=5, nlon=5).mean().plot(figsize=(18, 10))
plt.title('SST snapshot ($^{\circ}$C)',fontsize=18)


# Note the chunk structure of the original data: contiguous in the spatial dimension and chunked in the time dimension.
# This is not optimal for frequency-domain spectral analysis.

# In[ ]:


ds_orig.SST


# ## Rechunk the data
# 
# Rechunker allows us to transform the chunk structure of the dataset using Dask distributed computing in order to facilitate time-domain analysis of the dataset.
# Since we are only interested in specific variables, we rechunk a single variable (in this case SST) at a time.
# 
# 
# <div class="alert alert-info">
# 
# **Note:** 
# Rechunking requires the ability to write hundreds of GB to cloud storage.
# If you are running this notebook in a binder, your environment has read access to Google Cloud Storage, but not write access.
# (Providing anonymous write access would be a major security threat.)
# Therefore, below we display the code needed for rechunking in a markdown cell, rather than an executable cell.
# Instead, in Step 2 below we point to a rechunked dataset that has previously been written and then continue with the Fourier transform steps.
# </div>
# 
# Users running on Pangeo Cloud should be able to run this code with no changes, since credentials are automatically populated.
# Users who have their own Google Cloud account could modify this code by providing valid credentials and a path to their own storage buckets.

# ```python
# # Open data as Zarr group
# gcs = gcsfs.GCSFileSystem(requester_pays=True) # Connect to Google Cloud Storage (GCS)
# mapper = gcs.get_mapper(item.urlpath) # gcs.get_mapper() is like getting the pathname of the data on GCS
# zgroup = zarr.open_consolidated(mapper) # Open the dataset in Zarr format
# 
# # Select SST variable
# varname = 'SST' # name of variable you wish to rechunk (SST in this case)
# array = zgroup[varname] # select only the SST variable from the Zarr group
# 
# # Write intermediary data to a temporary path. Items stored here get deleted every 7 days.
# scratch_path = os.environ['PANGEO_SCRATCH']
# 
# # Define options needed by the `rechunk()` function of `rechunker`.
# max_mem = '1GB' # request memory for dask workers
# target_chunks = (365*5, 90, 180) # chunk in time (5 years) and space
# 
# # Set paths to temporary and target storage locations and set names of output data
# tmp_path = f'{scratch_path}/CESM_POP_hires_control/{varname}_tmp.zarr'
# target_path = f'paigem-pangeo/CESM_POP_hires_control/{varname}_target.zarr'
# ```

# ```python
# # Delete arrays of the same name that already exist at those paths
# def clear_targets():
#     for path in tmp_path,target_path:
#         try:
#             gcs.rm(path + '/.zarray')
#         except FileNotFoundError:
#             pass
# clear_targets()
# ```

# ```python
# # Create mappings to temporary and target storage buckets
# store_tmp = gcs.get_mapper(tmp_path)
# store_target = gcs.get_mapper(target_path)
# ```

# ```python
# # Call the `rechunk()` function. Below it shows the type of the source, intermediate, and target data.
# r = rechunk(array, target_chunks, max_mem,
#             store_target, temp_store=store_tmp, executor='dask')
# ```

# ### (Set up Dask cluster)
# 
# Note that to successfully run `rechunk()` it is necessary to spin up a Dask cluster. To see how we do so on Pangeo Cloud using Dask Gateway, see the following section on computing the power spectrum.

# #### Execute the `rechunk()` function
# 
# Executing this function will perform the rechunking and save the rechunked SST as intermediate data in the storage bucket designated by `target_path` defined above.

# ```python
# r.execute(retries=10) # `retries=10` sets the number of times a Dask worker will retry if computing fails (default is 0)
# ```

# ***

# ## Step 2: Compute the power spectrum of SST
# 
# The power spectrum is defined as follows, where a hat denotes a Fourier transform and the star denotes a complex conjugate.
# 
# $$
# \widehat{SST}^* \widehat{SST}
# $$

# Here we read in a previously rechunked SST field, which was obtained by following the same steps as shown above.

# In[ ]:


# Read previously rechunked data from my GCS storage bucket
gcs = gcsfs.GCSFileSystem(requester_pays=True)
varname = 'SST' 
target_path = f'pangeo-paigem/CESM_POP_hires_control/{varname}_target.zarr'
T_rechunked = dsa.from_zarr(gcs.get_mapper(target_path)) # load as dask file
T_rechunked


# ### Work with a 10-year subset of the rechunked data
# 
# Though we only show the workflow for one 10-year subset below, when running this workflow for scientific output we compute the power spectrum across several (in this case 10-year) windows. By averaging over all windows in our final plots, we are able to smooth out the naturally noisy Fourier Transform while also getting better statistics. 

# In[ ]:


# Create 10-year subset over first 10 years of output
yr1 = 0
yr2 = 10

T_rechunked = T_rechunked[365*yr1:365*yr2,:,:]

# Convert from dask to xarray DataArray
Txr_rechunked = xr.DataArray(T_rechunked,dims=['time','nlat','nlon'])
Txr_rechunked


# ### Define function to take power spectrum

# In[ ]:


def take_power_spectrum(var,real_arg):
    
    var = var.chunk({'time':None}) # there can be no chunking in the time dimension
    var_filled = var.fillna(0) # fill NaNs with zeros
    
    # Take power spectrum in the time domain, setting time to be a real dimension, with both a linear detrend and a windowing function
    var_hat = xrft.power_spectrum(var_filled,dim='time',real=real_arg,detrend='linear',window=True)
    
    return var_hat


# ### Compute the power spectrum of SST

# In[ ]:


# Call the
T_calc_power_spectrum = take_power_spectrum(Txr_rechunked,'time')

# Take only the real output, and immediately coarsen to 0.5 degree grid to reduce the memory usage
T_power_spectrum = T_power_spectrum.real.coarsen(nlat=5, nlon=5).mean()
T_power_spectrum


# ### Load Some Results
# 
# So far, all of our calculations have been "lazy".
# No computation has actually happened yet.
# In order to load some data, we will subset the data to a size that can fit into our notebook's memory.
# We select a region in the North Atlantic.

# In[ ]:


region = dict(nlat=slice(290, 380), nlon=slice(50, 220))
T_power_spectrum_NAtl = T_power_spectrum.isel(**region)
T_power_spectrum_NAtl


# #### Start Dask cluster
# 
# Now that we are ready to compute, we will start a Dask cluster.
# If running this notebook on Pangeo binder, this Dask cluster will spin up on Pangeo's Cloud allocation.
# We use Dask Gateway here (as opposed to, e.g., Dask Kubernetes), as Dask Gateway provides flexibility with cluster scaling, and is the Dask library installed in the Pangeo Cloud infrastructure.
# This notebook can also be run with other Dask configurations.
# This Dask cluster sets each worker to use 8 GB of memory, and allows for up to 20 workers to spin up (using adaptive scaling). 
# Note that this may take several minutes to spin up - once you see nonzero numbers below, you are ready to continue.

# In[ ]:


from dask_gateway import Gateway
gateway = Gateway()
options = gateway.cluster_options()
options.worker_memory = 8 # assign each worker 8 GB of memory
cluster = gateway.new_cluster(options)
cluster.adapt(minimum=1, maximum=20) # use adaptive scaling to allow up to 20 workers
cluster # print out cluster information - you can use this to view when the cluster spins up


# In[ ]:


# Start the client
from dask.distributed import Client
client = Client(cluster)
client # the link below can be used to access the Dask dashboard


# Calling `load` triggers the cluster to scale up elastically and finish the computation.

# In[ ]:


T_power_spectrum_NAtl.load()


# Re-mask the data:

# In[ ]:


mask = (ds_orig.SST[0].notnull().astype(int).coarsen(nlon=5, nlat=5).mean() > 0.2).isel(**region).load()
T_power_spectrum_NAtl = T_power_spectrum_NAtl.where(mask)


# ### Write to Zarr
# 
# At this point in our workflow, to do the full globe, we would write out the power spectrum computation back to disk. However, we again include this cell as Markdown only since Binder does not have write access.

# ```python
# # Set path to save power spectrum
# url = f'paigem-pangeo/CESM_POP_hires_control/SST_power_spectrum_yr{yr1}_{yr2}_earthcube_test.zarr'
# 
# # Save to Zarr
# T_power_spectrum_reset.to_dataset(name='SST_power_spectrum').to_zarr(gcs.get_mapper(url)) # need to convert xarray DataArray to Dataset first
# ```

# ***

# ## Step 3: Plot Results
# 
# ### Area Average Power Spectrum

# In[ ]:


T_power_spectrum_NAtl.mean(dim=('nlon', 'nlat')).plot(yscale='log', xscale='log')
plt.title('Power spectrum of SST in North Atlantic region')


# ### Maps of different frequency band averages

# In[ ]:


bands = {
    "> 1 year": slice(1/(3650), 1/400),
    "~ 1 year": slice(1/370, 1/360),
    "100 days - 1 year": slice(1/360, 1/100),
    "10 days - 100 days": slice(1/100, 1/10),
    "< 10 days": slice(1/10, None)
}


# In[ ]:


for name, band_slice in bands.items():
    plt.figure()
    T_power_spectrum_NAtl.sel(freq_time=band_slice).mean('freq_time').plot(robust=True)
    plt.title(name)


# In[ ]:




