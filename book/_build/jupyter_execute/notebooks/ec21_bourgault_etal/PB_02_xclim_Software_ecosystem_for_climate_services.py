#!/usr/bin/env python
# coding: utf-8

# # xclim: Software ecosystem for climate services

# ## Author(s)
# List authors, their current affiliations,  up-to-date contact information, and ORCID if available. Add as many author lines as you need.
# 
# - Author1 = {"name": "Pascal Bourgault", "affiliation": "Ouranos Inc", "email": "bourgault.pascal@ouranos.ca", "orcid": "0000-0003-1192-0403"}
# - Author2 = {"name": "Travis Logan", "affiliation": "Ouranos Inc", "email": "logan.travis@ouranos.ca", "orcid": "0000-0002-2212-9580"}
# - Author3 = {"name": "David Huard", "affiliation": "Ouranos Inc", "email": "huard.david@ouranos.ca", "orcid": "0000-0003-0311-5498"}
# - Author4 = {"name": "Trevor J. Smith", "affiliation": "Ouranos Inc", "email": "smith.trevorj@ouranos.ca", "orcid": "0000-0001-5393-8359"}

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#xclim:-Software-ecosystem-for-climate-services" data-toc-modified-id="xclim:-Software-ecosystem-for-climate-services-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>xclim: Software ecosystem for climate services</a></span><ul class="toc-item"><li><span><a href="#Author(s)" data-toc-modified-id="Author(s)-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>Author(s)</a></span></li><li><span><a href="#Purpose" data-toc-modified-id="Purpose-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Purpose</a></span></li><li><span><a href="#Technical-contributions" data-toc-modified-id="Technical-contributions-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>Technical contributions</a></span></li><li><span><a href="#Methodology" data-toc-modified-id="Methodology-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>Methodology</a></span></li><li><span><a href="#Results" data-toc-modified-id="Results-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>Results</a></span></li><li><span><a href="#Funding" data-toc-modified-id="Funding-1.6"><span class="toc-item-num">1.6&nbsp;&nbsp;</span>Funding</a></span></li><li><span><a href="#Keywords" data-toc-modified-id="Keywords-1.7"><span class="toc-item-num">1.7&nbsp;&nbsp;</span>Keywords</a></span></li><li><span><a href="#Citation" data-toc-modified-id="Citation-1.8"><span class="toc-item-num">1.8&nbsp;&nbsp;</span>Citation</a></span></li><li><span><a href="#Suggested-next-steps" data-toc-modified-id="Suggested-next-steps-1.9"><span class="toc-item-num">1.9&nbsp;&nbsp;</span>Suggested next steps</a></span></li><li><span><a href="#Acknowledgements" data-toc-modified-id="Acknowledgements-1.10"><span class="toc-item-num">1.10&nbsp;&nbsp;</span>Acknowledgements</a></span></li></ul></li><li><span><a href="#Setup" data-toc-modified-id="Setup-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Setup</a></span><ul class="toc-item"><li><span><a href="#Library-import" data-toc-modified-id="Library-import-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>Library import</a></span></li></ul></li><li><span><a href="#Data-import" data-toc-modified-id="Data-import-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Data import</a></span></li><li><span><a href="#Data-processing-and-analysis" data-toc-modified-id="Data-processing-and-analysis-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Data processing and analysis</a></span><ul class="toc-item"><li><span><a href="#Climate-indices-computation-and-ensemble-statistics" data-toc-modified-id="Climate-indices-computation-and-ensemble-statistics-4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>Climate indices computation and ensemble statistics</a></span><ul class="toc-item"><li><span><a href="#Creating-the-ensemble-dataset" data-toc-modified-id="Creating-the-ensemble-dataset-4.1.1"><span class="toc-item-num">4.1.1&nbsp;&nbsp;</span>Creating the ensemble dataset</a></span></li><li><span><a href="#Compute-climate-indicators" data-toc-modified-id="Compute-climate-indicators-4.1.2"><span class="toc-item-num">4.1.2&nbsp;&nbsp;</span>Compute climate indicators</a></span></li><li><span><a href="#Ensemble-statistics" data-toc-modified-id="Ensemble-statistics-4.1.3"><span class="toc-item-num">4.1.3&nbsp;&nbsp;</span>Ensemble statistics</a></span></li></ul></li><li><span><a href="#Bias-adjustment" data-toc-modified-id="Bias-adjustment-4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>Bias-adjustment</a></span><ul class="toc-item"><li><span><a href="#Simple-quantile-mapping" data-toc-modified-id="Simple-quantile-mapping-4.2.1"><span class="toc-item-num">4.2.1&nbsp;&nbsp;</span>Simple quantile mapping</a></span></li></ul></li><li><span><a href="#Finch-:-xclim-as-a-service" data-toc-modified-id="Finch-:-xclim-as-a-service-4.3"><span class="toc-item-num">4.3&nbsp;&nbsp;</span>Finch : xclim as a service</a></span><ul class="toc-item"><li><ul class="toc-item"><li><span><a href="#What-is-a-WPS?" data-toc-modified-id="What-is-a-WPS?-4.3.0.1"><span class="toc-item-num">4.3.0.1&nbsp;&nbsp;</span>What is a WPS?</a></span></li></ul></li><li><span><a href="#Chain-computations" data-toc-modified-id="Chain-computations-4.3.1"><span class="toc-item-num">4.3.1&nbsp;&nbsp;</span>Chain computations</a></span></li><li><span><a href="#Ensemble-statistics" data-toc-modified-id="Ensemble-statistics-4.3.2"><span class="toc-item-num">4.3.2&nbsp;&nbsp;</span>Ensemble statistics</a></span></li></ul></li></ul></li><li><span><a href="#References" data-toc-modified-id="References-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>References</a></span><ul class="toc-item"><li><ul class="toc-item"><li><span><a href="#Bias-adjustment" data-toc-modified-id="Bias-adjustment-5.0.1"><span class="toc-item-num">5.0.1&nbsp;&nbsp;</span>Bias-adjustment</a></span></li></ul></li></ul></li></ul></div>

# ## Purpose
# This notebook presents [`xclim`](https://xclim.readthedocs.io/en/stable/), a python package allowing easy manipulation and analysis of N-D climate dataset. It goes in details on what tools xclim provides, and how it simplifies and standardizes climate analyses workflows. The notebook also demontrates how xclim functionality can be exposed as a Web Processing Service ([WPS](http://opengeospatial.github.io/e-learning/wps/text/basic-main.html)) through [`finch`](https://pavics-sdi.readthedocs.io/projects/finch/en/latest/). This notebook was developped for Python 3.7+, but expects a special python environment and a running instance of finch. The configuration is optimized for running on `binder`.
# 
# 
# ## Technical contributions
# - Development of a climate analysis library regrouping many different tools, with special attention to community standards (especially [CF](https://cfconventions.org/)) and computationally efficient implementations.
#   * Numerically efficient climate indices algorithms yielding standard-compliant outputs;
#   * Model ensemble statistics and robustness metrics
#   * Bias-adjustment algorithms
#   * Other tools (not shown in this notebook) (subsetting, spatial analogs, calendar conversions)
# - Development of a web service exposing climate analysis processes through OGC's WPS protocol.
# 
# 
# ## Methodology
# This notebook shows two examples of climate analysis workflows:
# 
# 1. Calculating climate indices on an ensemble of simulation data, then computing ensemble statistics;
# 2. Bias-adjustment of model data.
# 3. Calculating climate indices through a web service (`finch`).
# 
# The examples interleave code and markdown cells to describe what is being done. They also include some graphics created with  holoviews, in order to give a clearer idea of what was achieved.
# 
# ## Results
# This notebook demontrates how the use of xclim can help simplify and standardize climate analysis workflows.
# 
# ## Funding
# 
# This project benefitted from direct and indirect funding from multiple sources, which are cited in the "Acknowledgements" section below.
# 
# ## Keywords
# 
# keywords=["climate", "xarray", "bias-adjustment", "wps"]
# 
# ## Citation
# 
# Bourgault, Pascal et al. Ouranos, 2021. xclim: Software ecosystem for climate services. Accessed 2021/05/15 at https://github.com/Ouranosinc/xclim/blob/earthcube-nb/docs/notebooks/PB_01_xclim_Software_ecosystem_for_climate_services.ipynb
# 
# 
# ## Suggested next steps
# 
# After reading through this notebook, we recommend looking at xclim's documentation, especially the [examples](https://xclim.readthedocs.io/en/stable/notebooks/index.html) page, where other notebooks extend the examples shown here. Xclim offers other features that are not shown here, but might be of interest to readers : [a spatial analogs submodule](https://xclim.readthedocs.io/en/stable/api.html#module-xclim.analog) and [internationalization utilities](https://xclim.readthedocs.io/en/stable/internationalization.html), among others. The [list of indicators](https://xclim.readthedocs.io/en/stable/indicators.html) is also a good stop, and is regularly updated with new algorithms. It covers most indicators from [ECAD](https://xclim.readthedocs.io/en/stable/indicators.html#icclim-indices) to more complex and specialized one, like the [Fire Weather Indexes](https://xclim.readthedocs.io/en/stable/indices.html#fire-weather-indices-submodule).
# 
# The web service demontrated here is part of the [`bird-house`](http://bird-house.github.io/) ecosystem and project. It is used by different platforms, such as [PAVICS](https://pavics.ouranos.ca/) and [climatedata.ca](https://climatedata.ca).
# 
# 
# ## Acknowledgements
# 
# Much of the development work of xclim was done to support climate services delivery at [Ouranos](https://ouranos.ca) and the [Canadian Centre for Climate Services](https://www.canada.ca/en/environment-climate-change/services/climate-change/canadian-centre-climate-services.html), and was funded by Environment and Climate Change Canada (ECCC). It has become a core component of the Power Analytics and Visualization for Climate Science (PAVICS) platform (Ouranos and CRIM 2021). PAVICS has been funded by [CANARIE](https://canarie.ca) and the Québec Fonds Vert, and through other projects also benefits from the support of the Canadian Foundation for Innovation (CFI) and the Fonds de Recherche du Québec (FRQ). PAVICS is a joint effort of Ouranos and the Centre de Recherche en Informatique de Montréal (CRIM). xclim itself
# 

# # Setup
# 
# ## Library import

# In[1]:


import xclim
# Also import some submodules for direct access
import xclim.ensembles  # Ensemble creation and statistics
import xclim.sdba  # Bias-adjustment

# Data manipulation
import xarray as xr
import xclim.testing

# For interacting with finch WPS
import birdy

# Visualizations and display
from pprint import pprint
from IPython.display import display
import matplotlib as mpl
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar

# For handling file paths
from pathlib import Path

# Some little figure tweaking
mpl.style.use('seaborn-deep')
mpl.style.use('seaborn-darkgrid')
get_ipython().run_line_magic('matplotlib', 'inline')
mpl.rcParams['figure.figsize'] = (12, 5)


# # Data import
# 
# In the first part, we use the "Ouranos standard ensemble of bias-adjusted climate scenarios version 1.0", bias-adjusted and downscaled data from Ouranos. It includes 11 simulations from different climate models, all produced in the RCP8.5 experiment of the CMIP5 project. More information can be found on [this page](https://pavics.ouranos.ca/datasets.html).

# In[2]:


get_ipython().system('du -hs data/*.nc')


# The second example, dealing with bias-adjustment, will use two sources. The model data comes from the [CanESM2 model](https://www.canada.ca/en/environment-climate-change/services/climate-change/science-research-data/modeling-projections-analysis/centre-modelling-analysis/models/second-generation-earth-system-model.html) of CCCma. The reference time series are extracted from the Adjusted and Homogenized Canadian Climate Data ([AHCCD](https://open.canada.ca/data/en/dataset/9c4ebc00-3ea4-4fe0-8bf2-66cfe1cddd1d)), from Environment and Climate Change Canada. Both datasets are distributed under the Open Government Licence (Canada). 
# We chose three weather stations representative of a reasonable range of climatic conditions and extracted time series over those locations from the model data.
# 
# Finally, the third section will use two datasets. A first one extracted from ECMWF's [ERA5 reanalysis](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5) includes time series extracted on a few locations for multiple variables. The second dataset is a subset of Natural Resources Canada (NRCAN)'s  daily dataset derived from observational data. See https://cfs.nrcan.gc.ca/projects/3/4.

# # Data processing and analysis
# 
# ## Climate indices computation and ensemble statistics
# 
# A large part of climate data analytics is based on the computation and comparison of climate indices. Xclim's main goal is to make it easier for researchers to process climate data and calculate those indices. It does so by providing a large library of climate indices, grouped in categories (CMIP's "realms") for convenience, but also by providing tools for all the small tweaks and pre- and post-processing steps that such workflows require.
# 
# xclim is built on [`xarray`](xarray.pydata.org/), which itself makes it easy to scale computations using [`dask`](https://dask.org). While many indices are conceptually  simple, we aim to provide implementations that balance performance and code complexity. This allows for ease in understanding and extending the code base.
# 
# In the following steps we will use xclim to:
# 
# - create an ensemble dataset, 
# - compute a few climate indices, 
# - calculate ensemble statistics on the results.
# 
# ### Creating the ensemble dataset
# 
# In xclim, an "ensemble dataset" is simply an xarray `Dataset` where multiple realizations of the same climate have been concatenated along a "realization" dimension. In our case here, as presented above, the realizations are in fact simulations from different models, but all have been bias-adjusted with the same reference. All ensemble-related functions are in the [`ensembles`](https://xclim.readthedocs.io/en/stable/api.html#module-xclim.ensembles) module.
# 
# Given a list of files, creating the ensemble dataset is straightforward:

# In[3]:


# Get a list of files in the `earthcube_data`:
# The sort is not important here, it only ensures a consistent notebook output for automated checks
files = sorted(list(Path('data').glob('*1950-2100.nc')))

# Open files into the ensemble
ens = xclim.ensembles.create_ensemble(files)
ens


# (Note: similarly to `xarray.open_mfdataset`, used internally, the dataset attributes are actually from the first file only.)
# 
# But how is this different from `xarray.concat(files, concat_dim='realization)`? While only slighlty more shorter to write, the main difference is the handling of incompatible calendars. Different models often use different [calendars](http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#calendar) and this causes issues when merging/concatening their data in xarray. Xclim provides [`xclim.core.calendar.convert_calendar`](https://xclim.readthedocs.io/en/stable/api.html#xclim.core.calendar.convert_calendar) to handle this. For example, in this ensemble, we have a simulation of the `HadGEM2` model which uses the "360_day" calendar:

# In[4]:


ds360 = xr.open_dataset('data/HadGEM2-CC_rcp85_singlept_1950-2100.nc')
ds360.time[59].data


# In[5]:


# Convert to the `noleap` calendar, see doc for details on arguments.
ds365 = xclim.core.calendar.convert_calendar(ds360, 'noleap', align_on='date')
ds365.time[59].data


# This conversion dropped the invalid days (for example: the 29th and 30th of February 1950) and convert the dtype of the time axis. This does alter the data, but it makes it possible to concatenate different model simulations into a single `DataArray`.
# 
# In `create_ensemble`, all members are converted to a common calendar. By default, it uses the normal numpy/pandas data type (a [range-limited](https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#timestamp-limitations) calendar equivalent to `proleptic_gregorian`), but for this model ensemble, it's better to use `noleap` (all years have 365 days), as it is the calendar used in most models included here.

# In[6]:


ens = xclim.ensembles.create_ensemble(files, calendar='noleap')
ens


# ### Compute climate indicators
# 
# Climate indicators are stored in `xclim.indicators.[realm]` and act as functions, taking `xarray.DataArray` as input variables and returning one (or more) `DataArray`. Most indicators also take keyword arguments to control various parameters. Under the hood, indicators are python objects, storing information on what the computation does, its inputs and its outputs. For example, let's look at the `tropical_nights` indicator:

# In[7]:


TN = xclim.indicators.atmos.tropical_nights

# Print general information:
print('Title: ', TN.title)
print('General description: ', TN.abstract, '\n')

# Info on the inputs
print('Inputs description:')
pprint(TN.parameters)
print()

# Info on the output
print('Output description:')
pprint(TN.cf_attrs)


# We use the indicator by calling it like a function with the variable and parameters. We will look at the _seasonal_ number of "tropical nights" with a threshold of 5°C. For frequency arguments, we use the same syntax as [pandas](https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#dateoffset-objects). Here `QS-DEC` will resample to seasons of three (3) months, starting in December (so DJF, MMA, JJA, SON).

# In[8]:


out = xclim.indicators.atmos.tropical_nights(tasmin=ens.tasmin, thresh='5 degC', freq='QS-DEC')
display(out.data)
out.isel(time=slice(0, 2)).load()


# Let's talk about the many things that happened here:
# 
# 1. **Input checks**
# 
#     Before the real computation, xclim performs some checks on the inputs. These can differ between indicators, but the most common are CF convention and input frequency checks. The former looks at the `standard_name` and `cell_methods` attributes of the input variable and raises warnings when they do not fit with what it expects. Here, the `standard_name` ("air_temperature") is correct, but the variable is missing the `cell_methods` attribute that would normally indicate that this is indeed the *minimum* temperature within each days.
# 
#     The second type of checks looks at the data itself. Most indicators expect data at a **daily** frequency, and if the input has a different frequency, the indicator checks will raise an error.
# 
# 
# 2. **Unit handling**
# 
#     Another check on the inputs is performed to ensure the correct units were passed. Above, in the `parameters` dictionary, we saw that `tropical_nights` expects `tasmin` and `thresh` to  both be "[temperature]" values. In the function, conversion are made to ensure compatible units. Here, `tasmin` is in Kelvin and `thresh` was given in Celsius. Under the hood, xclim uses [`pint`](pint.readthedocs.io/) to handle this, but with small modifications made so that CF-compliant unit strings are understood.
# 
# 
# 3. **Missing values handling**
# 
#     After the computation, for indicators involving a resampling operation, periods where at least one value was missing (i.e. was a `np.NaN`) are masked out. Here the first period is masked for all realization since there was no data for the month of december 1949. This is the default behavior and it can be changed, as we will show further down.
# 
# 
# 4. **Metadata formatting**
# 
#    The output's attributes were formatted to include some information from the passed parameters. Notice how the description includes the given threshold and frequency (*Seasonal number of tropical nights [...]*).
# 
# 
# 5. **Lazy computing**
# 
#    This is not a feature inherent to xclim, but, by default, `create_ensemble` will open the ensemble dataset using the `dask` backend, one "chunk" for each file. This means that all further calculations are _lazy_, ie: python remembers the operations but doesn't perform them as long as the result is not requested. This is why we needed to call `.load()` to see the first values of the output. This behaviour from `xarray` is extremely useful for computationally-intensive workflows as it allows `dask` to efficiently manage the memory and CPU ressources. Xclim aims to implement its indicators in a lazy and efficient way.
# 
# 
# Let's compute a few indicators and merge them into another ensemble dataset. This time, we will change the missing values handling to a more relaxed method, where periods are masked when at least 5% of the elements are missing. Moreover, we will set an option so that CF-checks warnings are ignored.
# 
# If you want to know more about the indicators we are computing, you can use `help(indicator)` or `indicator?` (in ipython/jupyter) to print their docstrings, which include most information shown above.

# In[9]:


# Set options using a "context". Options set here are reset to their defaults as soon as we exit the "with" block
with xclim.set_options(cf_compliance='log', check_missing='pct', missing_options={'pct': {'tolerance': 0.05}}):

    # Intermediate computation : Approximate mean daily temp by taking the mean of tasmin and tasmax
    tas = xclim.indicators.atmos.tg(tasmin=ens.tasmin, tasmax=ens.tasmax)

    # First day below, first day of the year where tasmin goes below "thresh" for at least "window" days
    # Note that we resample yearly, but starting in July, so that the full winter is included
    out_fda = xclim.indicators.atmos.first_day_below(tasmin=ens.tasmin, thresh='-5 degC', window=3, freq='AS-JUL')

    # Number frost days : total number of days where tasmin < 0 degC
    out_nfd = xclim.indicators.atmos.frost_days(tasmin=ens.tasmin, freq='AS-JUL')

    # Solid precipitation accumulation : total thickness of solid precipitation (estimated by pr when tas < 0°C)
    out_spa = xclim.indicators.atmos.solid_precip_accumulation(pr=ens.pr, tas=tas, freq='AS-JUL')

    # Cold spell freq : Number of spells where tas is under "thresh" for at least "window" days.
    out_csf = xclim.indicators.atmos.cold_spell_frequency(tas=tas, thresh='-10 degC', window=3, freq='AS-JUL')

out = xr.merge([out_fda, out_nfd, out_spa, out_csf])
out


# ### Ensemble statistics
# 
# All this is quite cool, but we still have a large amount of data to manage. In order to get an idea of the climatic evolution of our indicators and of the model uncertainty associated with it, we can compute ensemble percentiles. Here we extract the 10th, 50th (median) and 90th percentile of the distribution made up of the 11 members of the ensemble.

# In[10]:


out_perc = xclim.ensembles.ensemble_percentiles(out, values=[10, 50, 90], split=False)
out_perc


# However, the `first_day_below` indicator is problematic. The data is in a "day of year" (doy) format, an integer starting at `1` on January 1st and going to `365` on December 31st (in our *no leap* calendar). There might be years where the first day below 0°C happened _after_ December 31st, meaning that taking the percentiles directly will give incorrect results. The trick is simply to convert the "day of year" data to a format where the numerical order reflects the temporal order within each period.
# 
# The next graph show the _incorrect_ version. Notice how for the years where some members have values in January ($0 < doy \le 31$), those "small" numerical values are classified into the **10th** percentile, while they are in fact _later dates_, so should be considered so and classified as the **90th** percentile.

# In[11]:


out.first_day_below.plot(hue='realization', color='grey', alpha=0.5)
out_perc.first_day_below.plot(hue='percentiles')
plt.title('Incorrect percentiles computation');


# In[12]:


# Convert "doys" to the number of days elapsed since the time coordinate (which is 1st of July)
fda_days_since = xclim.core.calendar.doy_to_days_since(out.first_day_below)

# Take percentiles now that the numerical ordering is correct
fda_days_since_perc = xclim.ensembles.ensemble_percentiles(fda_days_since, values=[10, 50, 90], split=False)

# Convert back to doys for meaningful values
fda_doy_perc = xclim.core.calendar.days_since_to_doy(fda_days_since_perc)

# Replace the bad data with the good one
out_perc['first_day_below'] = fda_doy_perc


# Voilà! Lets enjoy our work!

# In[13]:


with ProgressBar():
    out_perc.load()


# In[14]:


for name, var in out_perc.data_vars.items():
    plt.figure()
    var.plot(hue='percentiles')
    plt.title(name)


# Notice that the doy conversion trick did work for `first_day_below`, but plotting this kind of variable is not straightforward, and not the point of this notebook.
# 
# ## Bias-adjustment
# 
# Besides climate indicators, xclim comes with the `sdba` (statistical downscaling and bias-adjustment) module. The module tries to follow a "modular" approach instead of implementing published methods as black box functions. It comes with some pre- and post-processing function in `sdba.processing` and a few adjustment objects in `sdba.adjustment`. Adjustment algorithms all conform to the train - adjust scheme, formalized within `Adjustment` classes. Given a reference time series (`ref`), historical simulations (`hist`) and simulations to be adjusted (`sim`), any bias-adjustment method would be applied by first estimating the adjustment factors between the historical simulation and the observations series, and then applying these factors to `sim`, which could be a future simulation.
# 
# Most algorithms implemented also perform the adjustment separately on temporal groups. For example, for a *Quantile Mapping* method, one would compute quantiles and adjustment factors for each day of the year individually or within a moving window, so that seasonal variations do not interfere with the result.
# 
# ### Simple quantile mapping
# 
# Let's start here with a simple *Empirical Quantile Mapping* adjustment for `pr`. Instead of adjusting a future period, we adjust the historical timeseries, so we can have a clear appreciation of the adjustment results. This adjustment method is based on [1].
# 
# **Multiplicative** and **Additive** modes : In most of the literature, bias-adjustment involving `pr` is done in a *multiplicative* mode, in opposition to the *additive* mode used for temperatures. In xclim.sdba, this is controlled with the `kind` keyword, which defaults to `'+'` (additive). This means that `sim` values will be multiplied to the adjustment factors in the case `kind='*'`, and added (or subtracted) with `'+'`.

# In[15]:


# We are opening the datasets with xclim.testing.open_dataset
# a convenience function for accessing xclim's test and example datasets, stored in a github repository.
dsref = xclim.testing.open_dataset('sdba/ahccd_1950-2013.nc')
dssim = xclim.testing.open_dataset('sdba/CanESM2_1950-2100.nc')

# Extract 30 years of the daily data
# The data is made of 3 time series, we take only one "location"
# `ref` is given in Celsius, we must convert to Kelvin to fit with hist and sim
ref = dsref.pr.sel(time=slice('1980', '2010'), location='Vancouver')
ref = xclim.core.units.convert_units_to(ref, 'mm/d')
hist = dssim.pr.sel(time=slice('1980', '2010'), location='Vancouver')
hist = xclim.core.units.convert_units_to(hist, 'mm/d')

# We want to group data on a day-of-year basis with a 31-day moving window
group_doy_31 = xclim.sdba.Grouper('time.dayofyear', window=31)

# Create the adjustment object
EQM = xclim.sdba.EmpiricalQuantileMapping(nquantiles=15, group=group_doy_31, kind='*')

# Train it
EQM.train(ref, hist)

# Adjust hist
# And we estimate correct adjustment factors by interpolating linearly between quantiles
scen = EQM.adjust(hist, interp='linear')


# In[16]:


# Plot the season cycle
fig, ax = plt.subplots()
ref.groupby('time.dayofyear').mean().plot(ax=ax, label='Reference')
hist.groupby('time.dayofyear').mean().plot(ax=ax, label='Model - raw')
scen.groupby('time.dayofyear').mean().plot(ax=ax, label='Model - adjusted')
ax.legend();


# To inspect the trained adjustment data, we can access it with `EQM.ds`:

# In[17]:


EQM.ds


# More complex adjustment workflows can be created using the multiple tools exposed by xclim in `sdba.processing` and `sdba.detrending`. A more complete description of those is a bit too complex for the purposes of this notebook, but we invite interested users to look at [the documentation](https://xclim.readthedocs.io/en/stable/notebooks/sdba.html).
# 
# Let's make another simple example, using the `DetrendedQuantileMapping` method. Compared to the `EmpiricalQuantileMapping` method, this one will normalize and detrend the inputs before computing the adjustment factors. This ensures a better quantile-to-quantile comparison when the climate change signal is strong. In reality, the same steps could all be performed explicitly, making use of xclim's modularity, as is shown in the doc examples. It is based on [2].

# In[18]:


# Extract the data.
ref = xclim.core.units.convert_units_to(
    dsref.tasmax.sel(time=slice('1980', '2010'), location='Vancouver'),
    'K'
)
hist = xclim.core.units.convert_units_to(
    dssim.tasmax.sel(time=slice('1980', '2010'), location='Vancouver'),
    'K'
)
sim = xclim.core.units.convert_units_to(
    dssim.tasmax.sel(location='Vancouver'),
    'K'
)


# Adjustment
DQM = xclim.sdba.DetrendedQuantileMapping(nquantiles=15, group=group_doy_31, kind='+')
DQM.train(ref, hist)

scen = DQM.adjust(sim, interp='linear', detrend=1)


# In[19]:


fig, ax = plt.subplots()
ref.groupby('time.dayofyear').mean().plot(ax=ax, label='Reference'),
hist.groupby('time.dayofyear').mean().plot(ax=ax, label='Simulation - raw')
scen.sel(time=slice('1981', '2010')).groupby('time.dayofyear').mean().plot(ax=ax, label='Simulation - adjusted')
ax.legend();


# Just for fun, let's compute the annual mean of the daily maximum temperature and plot the three time series:

# In[20]:


with xclim.set_options(cf_compliance='log'):
    ref_tg = xclim.indicators.atmos.tx_mean(tasmax=ref)
    sim_tg = xclim.indicators.atmos.tx_mean(tasmax=sim)
    scen_tg = xclim.indicators.atmos.tx_mean(tasmax=scen)

fig, ax = plt.subplots()
ref_tg.plot(ax=ax, label='Reference')
sim_tg.plot(ax=ax, label='Model - Raw')
scen_tg.plot(ax=ax, label='Model - Adjusted')
ax.legend();


# ## Finch : xclim as a service
# 
# This final part of the notebook demonstrates the use of `finch` to compute indicators over the web using the Web Processing Service (WPS) standard. When running this notebook on binder, a local server instance of finch is already running on port 5000. We will interact with it using [`birdy`](https://birdy.readthedocs.io/en/latest/) a lightweight WPS client built upon [`OWSLib`](https://geopython.github.io/OWSLib/).
# 
# #### What is a WPS?
# 
# Standardized by the OGC, WPS is a standard describing web requests and responses to and from geospatial processing services. The standard defines how process execution can be requested, and how the output from that process is handled. `finch` provides such a service, bundling xclim's indicator and some other tools into individual processes.
# 
# It's possible with WPS to either embed input data in the request, or pass a link to the data. In the latter case, the server will download data locally before performing its calculations. It's also possible to use DAP links so the web server only streams the data it actually needs. In the following example, we use a subset of the ERA5 reanalysis data, available on xclim's testdata GitHub repository.

# In[21]:


from birdy import WPSClient

wps = WPSClient('http://localhost:5000')
dataurl = 'https://github.com/Ouranosinc/xclim-testdata/raw/main/ERA5/daily_surface_cancities_1990-1993.nc'


# In[22]:


# As in xclim, processes and their inputs are documented 
# by parsing the indicator's attributes:
help(wps.growing_degree_days)


# In[23]:


# Through birdy, the indicator call looks a lot like the one from xclim
# but passing an url instead of an xarray object.
# Also, the `variable` argument tells which variable from the dataset to use.
response = wps.growing_degree_days(
    dataurl,
    thresh='5 degC',
    freq='MS',
    variable='tas'
)

response.get()


# The response we received is a list of URLs pointing to output files. Birdy makes it easy to open links to datasets directly in xarray:

# In[24]:


out = response.get(asobj=True).output_netcdf
out.growing_degree_days.plot(hue='location');


# ### Chain computations
# 
# As the goal of `finch` is to free the user of most, if not all, computation needs, it includes many other processes for all common climate analysis tasks. For example, here we will take maps of `tasmax` and `tasmin` over southern Québec, subset them with a square bounding box of lon/lat, compute `tas` as the mean of the other two and then finally, compute the growing season length with custom parameters.
# 
# All processes will be fed with the output of the previous one, so that no intermediate result is ever downloaded to the client. Only the last output is directly downloaded to file using python's `urllib`.

# In[25]:


tasmax_url = "https://github.com/Ouranosinc/xclim-testdata/raw/main/NRCANdaily/nrcan_canada_daily_tasmax_1990.nc"
tasmin_url = "https://github.com/Ouranosinc/xclim-testdata/raw/main/NRCANdaily/nrcan_canada_daily_tasmin_1990.nc"
pr_url = "https://github.com/Ouranosinc/xclim-testdata/raw/main/NRCANdaily/nrcan_canada_daily_pr_1990.nc"

# Subset the maps
resp_tx = wps.subset_bbox(tasmax_url, lon0=-72, lon1=-68, lat0=46, lat1=50)
resp_tn = wps.subset_bbox(tasmin_url, lon0=-72, lon1=-68, lat0=46, lat1=50)
resp_pr = wps.subset_bbox(pr_url, lon0=-72, lon1=-68, lat0=46, lat1=50)

# Compute tas
# We don't need to provide variable names because variables have the same names as the argument names
resp_tg = wps.tg(tasmin=resp_tn.get().output, tasmax=resp_tx.get().output)

# Compute growing season length, as defined by the period between
# Start/End when tas is over/under 10°C for 6 consecutive days
resp_lpr = wps.liquid_precip_ratio(
    pr=resp_pr.get().output,
    tas=resp_tg.get().output_netcdf,
    freq='YS'
)


# In[26]:


# Plot the liquid precipitation ratio of 1990
out_lpr = resp_lpr.get(asobj=True).output_netcdf
out_lpr.liquid_precip_ratio.plot()


# In[27]:


# or we can plot it as a map!
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Set a CRS transformation:
crs = ccrs.PlateCarree()

# Add some map elements
ax = plt.subplot(projection=crs)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE)
out_lpr.liquid_precip_ratio.plot(ax=ax, transform=crs, cmap="rainbow")

plt.show()


# Output attributes are formatted the same as when xclim is used locally:

# In[28]:


out_lpr.liquid_precip_ratio


# ### Ensemble statistics
# 
# In addition to the indicators and the subsetting processes, finch wraps xclim's ensemble statistics functionality in all-in-one processes. All indicator are provided in several ensemble-subset processes, where the inputs are subsetted before the indicator is computed and ensemble statistics taken. For example, `ensemble_gridpoint_growing_season_length`, which computes `growing_season_length` over a single gridpoint for all members of the specified ensemble. However, at least for the current version, ensemble datasets are hardcoded in the configuration of the server instance and one can't send their own list of members.
# 
# In order not to overload Ouranos' servers where these current hardcoded datasets are stored, no example of this is run by default here. But interested users are invited to uncomment the following lines to see for themselves.
# 
# The ensemble data used by this example is a bias-adjusted and downscaled product from the Pacific Climate Impacts Consortium (PCIC). It is described [here](https://www.pacificclimate.org/data/statistically-downscaled-climate-scenarios).

# In[29]:


help(wps.ensemble_grid_point_max_n_day_precipitation_amount)


# In[30]:


## Get maximal 5 days precipitation amount for Niagara Falls between 2040 and 2070
## for 12 rcp45 datasets bias-adjusted and downscaled with BCCAQv2 by PCIC
## Careful : this call will take a long time to succeed!

#resp_ens = wps.ensemble_grid_point_max_n_day_precipitation_amount(
#    lat=43.08,
#    lon=-79.07,
#    start_date='2040-01-01',
#    end_date='2070-12-31',
#    ensemble_percentiles=[10, 50, 90],
#    rcp='rcp45',
#    models='PCIC12',
#    window=5,
#    freq='YS',
#    data_validation='warn',
#)


# In[31]:


#resp_ens.get(asobj=True).output


# This concludes our xclim and finch tour! We hope these tools can be useful for your workflows. The development of both is always on-going, visit us at our GitHub pages for any questions or if you want to contribute! 
# 
# [`xclim`](https://github.com/Ouranosinc/xclim/) 
# 
# [`finch`](https://github.com/bird-house/finch/)

# # References
# 
# ### Bias-adjustment
# [1] Dequé, M. (2007). Frequency of precipitation and temperature extremes over France in an anthropogenic scenario: Model results and statistical correction according to observed values. Global and Planetary Change, 57(1–2), 16–26. https://doi.org/10.1016/j.gloplacha.2006.11.030
# 
# [2] Cannon, A. J., Sobie, S. R., & Murdock, T. Q. (2015). Bias correction of GCM precipitation by quantile mapping: How well do methods preserve changes in quantiles and extremes? Journal of Climate, 28(17), 6938–6959. https://doi.org/10.1175/JCLI-D-14-00754.1
# 
# ### Software used
# 
# **xarray 0.18.2** : Stephan Hoyer, Joe Hamman, Maximilian Roos, keewis, Deepak Cherian, Clark Fitzgerald, … Benoit Bovy. (2021, May 19). pydata/xarray: v0.18.2 (Version v0.18.2). Zenodo. http://doi.org/10.5281/zenodo.4774304
# 
# **numpy** : Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 0.1038/s41586-020-2649-2
# 
# **dask** : Dask Development Team (2016). Dask: Library for dynamic task scheduling
# URL https://dask.org
# 
# **cftime 1.4** : Jeff Whitaker, Constantine Khrulev, Spencer Clark, Joe Hamman, Bill Little, Ryan May, … Stefan. (2021, January 31). Unidata/cftime: version 1.4.0 release (Version v1.4.0rel). Zenodo. http://doi.org/10.5281/zenodo.4484141
# 

# In[ ]:




