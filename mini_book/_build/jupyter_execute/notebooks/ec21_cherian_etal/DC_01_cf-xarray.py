#!/usr/bin/env python
# coding: utf-8

# # cf-xarray: Scale your analysis across datasets with less data wrangling and more metadata handling
# 
# <img src="https://github.com/xarray-contrib/cf-xarray/blob/main/doc/_static/full-logo.png?raw=true" width="40%" align="center">
# 

# ## Author(s)
# 
# - Author1 = {"name": "Mattia Almansi", "affiliation": "National Oceanography
#   Centre", "email": "mattia.almansi@noc.ac.uk", "orcid": "0000-0001-6849-3647"}
# - Author2 = {"name": "Deepak Cherian", "affiliation": "National Center for
#   Atmospheric Research", "email": "deepak@cherian.net", "orcid":
#   "0000-0002-6861-8734"}
# - Author3 = {"name": "Pascal Bourgault", "affiliation": "Ouranos Inc.", "email":
#   "bourgault.pascal@ouranos.ca", "orcid": "0000-0003-1192-0403"}
# 

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#cf-xarray:-Scale-your-analysis-across-datasets-with-less-data-wrangling-and-more-metadata-handling" data-toc-modified-id="cf-xarray:-Scale-your-analysis-across-datasets-with-less-data-wrangling-and-more-metadata-handling-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>cf-xarray: Scale your analysis across datasets with less data wrangling and more metadata handling</a></span><ul class="toc-item"><li><span><a href="#Author(s)" data-toc-modified-id="Author(s)-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>Author(s)</a></span></li><li><span><a href="#Purpose" data-toc-modified-id="Purpose-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Purpose</a></span></li><li><span><a href="#Technical-contributions" data-toc-modified-id="Technical-contributions-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>Technical contributions</a></span></li><li><span><a href="#Methodology" data-toc-modified-id="Methodology-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>Methodology</a></span></li><li><span><a href="#Results" data-toc-modified-id="Results-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>Results</a></span></li><li><span><a href="#Funding" data-toc-modified-id="Funding-1.6"><span class="toc-item-num">1.6&nbsp;&nbsp;</span>Funding</a></span></li><li><span><a href="#Keywords" data-toc-modified-id="Keywords-1.7"><span class="toc-item-num">1.7&nbsp;&nbsp;</span>Keywords</a></span></li><li><span><a href="#Citation" data-toc-modified-id="Citation-1.8"><span class="toc-item-num">1.8&nbsp;&nbsp;</span>Citation</a></span></li><li><span><a href="#Work-In-Progress---improvements" data-toc-modified-id="Work-In-Progress---improvements-1.9"><span class="toc-item-num">1.9&nbsp;&nbsp;</span>Work In Progress - improvements</a></span></li><li><span><a href="#Suggested-next-steps" data-toc-modified-id="Suggested-next-steps-1.10"><span class="toc-item-num">1.10&nbsp;&nbsp;</span>Suggested next steps</a></span></li><li><span><a href="#Acknowledgements" data-toc-modified-id="Acknowledgements-1.11"><span class="toc-item-num">1.11&nbsp;&nbsp;</span>Acknowledgements</a></span></li></ul></li><li><span><a href="#Setup" data-toc-modified-id="Setup-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Setup</a></span><ul class="toc-item"><li><span><a href="#Library-import" data-toc-modified-id="Library-import-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>Library import</a></span></li></ul></li><li><span><a href="#Parameter-definitions" data-toc-modified-id="Parameter-definitions-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Parameter definitions</a></span></li><li><span><a href="#Data-import" data-toc-modified-id="Data-import-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Data import</a></span><ul class="toc-item"><li><span><a href="#MOM6" data-toc-modified-id="MOM6-4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>MOM6</a></span></li><li><span><a href="#CMIP6-Ocean" data-toc-modified-id="CMIP6-Ocean-4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>CMIP6 Ocean</a></span></li><li><span><a href="#CMIP6-Ice" data-toc-modified-id="CMIP6-Ice-4.3"><span class="toc-item-num">4.3&nbsp;&nbsp;</span>CMIP6 Ice</a></span></li><li><span><a href="#NCEP" data-toc-modified-id="NCEP-4.4"><span class="toc-item-num">4.4&nbsp;&nbsp;</span>NCEP</a></span></li></ul></li><li><span><a href="#Data-processing-and-analysis" data-toc-modified-id="Data-processing-and-analysis-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>Data processing and analysis</a></span><ul class="toc-item"><li><span><a href="#Overview-of-cf-xarray" data-toc-modified-id="Overview-of-cf-xarray-5.1"><span class="toc-item-num">5.1&nbsp;&nbsp;</span>Overview of cf-xarray</a></span><ul class="toc-item"><li><span><a href="#.cf-is-an-entrypoint-for-cf-xarray-functionality" data-toc-modified-id=".cf-is-an-entrypoint-for-cf-xarray-functionality-5.1.1"><span class="toc-item-num">5.1.1&nbsp;&nbsp;</span><code>.cf</code> is an entrypoint for cf-xarray functionality</a></span></li><li><span><a href="#Use-CF-metadata-in-standard-xarray-methods" data-toc-modified-id="Use-CF-metadata-in-standard-xarray-methods-5.1.2"><span class="toc-item-num">5.1.2&nbsp;&nbsp;</span>Use CF metadata in standard xarray methods</a></span></li><li><span><a href="#Dictionaries-mapping-CF-keys-to-variable-names" data-toc-modified-id="Dictionaries-mapping-CF-keys-to-variable-names-5.1.3"><span class="toc-item-num">5.1.3&nbsp;&nbsp;</span>Dictionaries mapping CF keys to variable names</a></span></li><li><span><a href="#Dealing-with-incomplete-metadata" data-toc-modified-id="Dealing-with-incomplete-metadata-5.1.4"><span class="toc-item-num">5.1.4&nbsp;&nbsp;</span>Dealing with incomplete metadata</a></span></li><li><span><a href="#Indexing-using-CF-keys" data-toc-modified-id="Indexing-using-CF-keys-5.1.5"><span class="toc-item-num">5.1.5&nbsp;&nbsp;</span>Indexing using CF keys</a></span></li><li><span><a href="#Automagic-plotting" data-toc-modified-id="Automagic-plotting-5.1.6"><span class="toc-item-num">5.1.6&nbsp;&nbsp;</span>Automagic plotting</a></span></li><li><span><a href="#CF-keys-expansion" data-toc-modified-id="CF-keys-expansion-5.1.7"><span class="toc-item-num">5.1.7&nbsp;&nbsp;</span>CF keys expansion</a></span></li></ul></li><li><span><a href="#Use-case-examples" data-toc-modified-id="Use-case-examples-5.2"><span class="toc-item-num">5.2&nbsp;&nbsp;</span>Use case examples</a></span><ul class="toc-item"><li><span><a href="#Seamlessly-extract-statistics-from-CF-compliant-datasets" data-toc-modified-id="Seamlessly-extract-statistics-from-CF-compliant-datasets-5.2.1"><span class="toc-item-num">5.2.1&nbsp;&nbsp;</span>Seamlessly extract statistics from CF-compliant datasets</a></span></li><li><span><a href="#Standardize-datasets" data-toc-modified-id="Standardize-datasets-5.2.2"><span class="toc-item-num">5.2.2&nbsp;&nbsp;</span>Standardize datasets</a></span></li><li><span><a href="#Regridding-with-xESMF" data-toc-modified-id="Regridding-with-xESMF-5.2.3"><span class="toc-item-num">5.2.3&nbsp;&nbsp;</span>Regridding with xESMF</a></span></li></ul></li></ul></li></ul></div>
# 

# ## Purpose
# 
# This notebook demonstrates how the **cf-xarray** Python package (Cherian et
# al, 2021) helps climate data scientists to process several CF-compliant datasets
# from a variety of sources. Under the hood, CF-xarray decodes and makes use of
# the widely adopted [Climate and Forecast (CF) conventions](cfconventions.org/).
# Therefore, workflows integrating cf-xarray do not need knowledge of arbitrary
# dataset-specific metadata such as variable names.
# 
# Xarray (Hoyer et al, 2021) is a Python package that enables easy and convenient
# labelled data analytics by allowing users to leverage metadata such as dimension
# names and coordinate labels. Xarray provides two core data structures:
# 
# - DataArray: a container wrapping a single multidimensional array with metadata.
# - Dataset: a dict-like container of multiple DataArrays.
# 
# cf-xarray uses Xarray's plugin interface, or "accessor", to provide extensive
# functionality on both Datasets and DataArrays under the `.cf` namespace.
# 
# For example, the zonal average of an Xarray Dataset `ds` is seamlessly
# calculated as `ds.cf.mean("longitude")` on any CF-compliant dataset, regardless
# of the actual name of the "longitude" variable (e.g., `"lon"`, `"lon_rho"`,
# `"long"`, ...).
# 
# ## Technical contributions
# 
# Development of cf-xarray, an extension library that
# 
# 1. adds awareness of the Climate and Forecast (CF) conventions to core Xarray
#    functionality.
# 1. Provides utility functions with minimal dependencies, allowing easy
#    integration of CF-aware functionality in other packages such as xesmf (for
#    regridding).
# 
# ## Methodology
# 
# This notebook can be executed in a pre-configured interactive environment at:
# https://binder.pangeo.io/v2/gh/malmans2/cf-xarray-earthcube/main?filepath=DC_01_cf-xarray.ipynb
# 
# ## Results
# 
# This notebook contains use case examples demonstrating the following
# functionalities of cf-xarray:
# 
# 1. Seamless analysis of various CF-compliant datasets.
# 2. Standardization of datasets to comply with CF conventions.
# 3. Inference of grid-cell coordinates and bounds using CF conventions.
# 4. Integration with other libraries (xESMF).
# 
# ## Funding
# 
# N/A
# 
# ## Keywords
# 
# Include up to 5 keywords, using the template below.
# 
# keywords=["cf-conventions", "xarray", "netcdf"]
# 
# ## Citation
# 
# Almansi, M. and Cherian, D. and Bourgault, P. (2021). cf-xarray: Scale your
# analysis across datasets with less data wrangling and more metadata handling.
# _2021 EarthCube Annual Meeting_. Accessed 14/5/2021 at
# https://github.com/malmans2/cf-xarray-earthcube
# 
# ## Work In Progress - improvements
# 
# N/A
# 
# ## Suggested next steps
# 
# It is common for datasets to not be perfectly CF-compliant. Here we work around
# these deficiencies using `assign_coordinates_and_cell_measures`. `cf_xarray` is
# considering adding heuristics to guess such metadata attributes, possibly using
# other metadata conventions such as SGRID.
# 
# ## Acknowledgements
# 
# We acknowledge contributions from all cf-xarray contributors:
# https://github.com/xarray-contrib/cf-xarray/graphs/contributors. We also
# acknowledge MetPy for providing inspiration and various criteria for identifying
# CF variables. Discussion with Jon Thielen was instrumental in development of
# cf-xarray. We also acknowledge the Pangeo project for collating and providing an
# immense amount of datasets on the cloud that motivated this work, as well as
# provoking, enabling, and fostering discussions that led to the development of
# this project.
# 

# # Setup
# 
# ## Library import
# 

# In[ ]:


# The package demonstrated here
import cf_xarray as cfxr

# For parallelization
import dask

# For loading shared data
import intake

# Visualizations
import matplotlib.pyplot as plt

# For basic data manipulation
import numpy as np
import xarray as xr

# For regridding
import xesmf as xe

# silence a minor warning
dask.config.set(**{"array.slicing.split_large_chunks": False})


# # Parameter definitions
# 

# In[ ]:


# Paths and urls pointing to data
MOM6_GRID_PATH = "./data/ocean_grid_sym_OM4_05.nc"
MOM6_DATA_URL = "http://35.188.34.63:8080/thredds/dodsC/OM4p5/ocean_monthly_z.200301-200712.nc4"
CMIP6_OCE_CATALOG = "https://storage.googleapis.com/cmip6/pangeo-cmip6.json"
CMIP6_OCE_EXPERIMENT = dict(
    table_id="Omon",
    grid_label="gn",
    source_id="ACCESS-CM2",
    # source_id="GFDL-CM4",
    member_id="r1i1p1f1",
    experiment_id=["historical"],
    variable_id=["thetao", "volcello"],
)
CMIP6_ICE_PATH = "data/sic_SImon_CCCma-CanESM5_ssp245_r13i1p2f1_2020.nc"
NCEP_PATH = "data/air_temperature.nc"


# # Data import
# 
# This notebook uses a variety of publicly available data:
# 
# - `ds_mom6`: Data from the Modular Ocean Model - v6 (Adcroft et al, 2019)
# - `ds_cmip6_oce`: Ocean Data from the Climate Model Intercomparison Project -
#   Phase 6
# - `ds_cmip6_ice_ice`: Ice Data from the Climate Model Intercomparison Project -
#   Phase 6
# - `ds_ncep`: Data from the National Centers for Atmospheric Prediction
#   Reanalysis (Kalnay et al, 1996)
# 

# In[ ]:


def assign_coordinates_and_cell_measures(ds):

    """
    Functions to add missing CF metadata (coordinates and cell measures).
    Fully CF-compliant datasets do not need this pre-processing.
    Functions to automatically assign missing coordinates
    and measures metadata will be implemented in cf_xarray.
    See https://github.com/xarray-contrib/cf-xarray/issues/201

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset to modify
    """

    for varname, variable in ds.data_vars.items():

        # Add coordinates attribute when the dimensions
        # of a coordinate variable are a subset of those of
        # a variable
        coordinates = []
        for coord in sum(ds.cf.coordinates.values(), []):
            if set(ds[coord].dims) <= set(variable.dims):
                coordinates.append(coord)
        # sets an attribute like "geolon geolat"
        if coordinates:
            variable.attrs["coordinates"] = " ".join(coordinates)
        else:
            variable.attrs.pop("coordinates", None)

        # Add cell_measures attribute when appropriate measures are available
        cell_measures = {}
        possible_measures = {
            "cell_thickness",
            "cell_area",
            "ocean_volume",
        } & set(ds.cf.standard_names)
        for stdname in possible_measures:
            key = stdname.split("_")[-1]
            value = ds.cf.standard_names[stdname]
            for measure in value:
                if (
                    set(ds[measure].dims) <= set(variable.dims)
                    and measure != varname
                ):
                    cell_measures[key] = measure

        if cell_measures:
            # sets an attribute like "area: areacello volume: volcello"
            variable.attrs["cell_measures"] = " ".join(
                [f"{k}: {v}" for k, v in cell_measures.items()]
            )
        else:
            variable.attrs.pop("cell_measures", None)


# ## MOM6
# 
# Read grid and data variables for a MOM6 ocean model simulation
# 

# In[ ]:


# Open grid and data variables, then merge them together to one dataset
grid = xr.open_dataset(MOM6_GRID_PATH, chunks={})
ds = xr.open_dataset(
    MOM6_DATA_URL,
    chunks={"time": 1},
)
ds_mom6 = xr.merge([grid, ds], compat="override")

# Illustrate the equivalent of a curvilinear grid case,
# where axes and coordinates are different
axes = ["xh", "xq", "yh", "yq"]
ds_mom6 = ds_mom6.drop_vars(axes)
ds_mom6 = ds_mom6.assign_coords({axis: ds_mom6[axis] for axis in axes})
ds_mom6 = ds_mom6.set_coords(
    [
        var
        for var in ds_mom6.variables
        for prefix in ["geo"]
        if var.startswith(prefix)
    ]
)
assign_coordinates_and_cell_measures(ds_mom6)


# ## CMIP6 Ocean
# 
# Read a historical CMIP6 simulation from the ACCESS-OM2 climate modelling system.
# 

# In[ ]:


# Use intake-esm to access data on Pangeo Cloud
col = intake.open_esm_datastore(CMIP6_OCE_CATALOG)
cat = col().search(**CMIP6_OCE_EXPERIMENT)

ddict = cat.to_dataset_dict(
    zarr_kwargs={
        "consolidated": True,
        "decode_times": True,
        "use_cftime": True,
    }
)
_, ds_cmip6_oce = ddict.popitem()
assign_coordinates_and_cell_measures(ds_cmip6_oce)


# ## CMIP6 Ice
# 

# In[ ]:


ds_cmip6_ice = xr.open_dataset(CMIP6_ICE_PATH)


# ## NCEP
# 

# In[ ]:


ds_ncep = xr.tutorial.open_dataset(NCEP_PATH)


# # Data processing and analysis
# 

# ## Overview of cf-xarray
# 
# cf-xarray uses Xarray's plugin interface, or "accessors", to provide extensive
# functionality on both Datasets and DataArrays under the `.cf` namespace.
# 

# ### `.cf` is an entrypoint for cf-xarray functionality
# 
# When `cf_xarray` is imported, the `cf` accessor is automatically added to any
# `xarray` object. cf-xarray is able to wrap most of Xarray's functions. The repr
# for `Dataset.cf` prints out a list of detected "CF names" and the corresponding
# dataset-specific variable names. cf-xarray parses the attributes associated with
# each variable in the Dataset to build these mappings
# 

# In[ ]:


# After `import cf_xarray`, the cf_xarray accessor has been added to the xarray object
ds_mom6.cf


# ### Use CF metadata in standard xarray methods
# 
# With standard `xarray` syntax, one would specify the variable names on the
# right-hand side of the mappings printed above. With cf_xarray, one can instead
# use the standardized "CF names" on the left-hand side.
# 
# For example, the next two cells show two ways of calculating an average along
# the vertical dimension
# 

# In[ ]:


# xarray way:
ds_mom6.mean(["z_i", "z_l"])


# In[ ]:


# cf_xarray way:
# By calling .cf.mean we can provide "Z" which is then rewritten to ["z_i", "z_l"]
# This statement is entirely equivalent to ds_mom6.mean(["z_i", "z_l"])
ds_mom6.cf.mean("Z")


# `cf_xarray` knows that `z_i` and `z_l` correspond `Z` axes because the two
# variables have a CF-compliant attribute `cartesian_axis: Z`. A full list of
# criteria used by cf-xarray is documented
# [here](https://cf-xarray.readthedocs.io/en/latest/criteria.html).
# 

# In[ ]:


for var_name in ["z_i", "z_l"]:
    print(f"{var_name}: {ds_mom6[var_name].attrs['cartesian_axis']}")


# ### Dictionaries mapping CF keys to variable names
# 
# The example object contains variables lying on staggered grids. Therefore, a CF
# key can be associated with multiple variables. `cf_xarray` provides several
# properties that return dictionaries mapping CF keys to lists of variable names,
# such as:
# 
# - `.cf.axes`
# - `.cf.coordinates`
# - `.cf.cell_measures`
# - `.cf.standard_names`
# - `.cf.bounds`
# 

# In[ ]:


# maps "axes" to variable names
ds_mom6.cf.axes


# In[ ]:


# maps CF standard name to variable name
ds.cf.standard_names


# ### Dealing with incomplete metadata
# 
# The usefulness of cf-xarray fully depends on the amount of CF-compliant metadata
# present in a dataset. Although many datasets have incomplete metadata, in most
# cases one can guess appropriate metadata by looking at variable names.
# 
# We can use `Dataset.cf.guess_coord_axis` to identify, guess, and add missing CF
# metadata for "axes" ("X", "Y", "Z", "T") and "coordinates" ("latitude",
# "longitude", "time"). It does so by using regular expressions to parse variable
# names and make reasonable guesses.
# 
# We also demonstrate the `verbose` mode so that the user can double check
# cf-xarray's inferences.
# 
# First, note that no `X` or `Y` axes variables have been detected in the current
# dataset:
# 

# In[ ]:


# The "axes" property maps axes names 'X', 'Y', 'Z', 'T' to variable names in the dataset
# here the metadata only identify 'Z' and 'T'
ds_mom6.cf.axes


# Now we will ask `cf_xarray` to autoguess more axes:
# 

# In[ ]:


ds_mom6 = ds_mom6.cf.guess_coord_axis(verbose=True)


# In[ ]:


# The `X` and `Y` axes variables that have been detected are sensible!
ds_mom6.cf.axes


# ### Indexing using CF keys
# 
# CF metadata precisely describes the physical quantities being represented by all
# variables. More importantly, CF conventions also describe links between
# different variables in a dataset.
# 
# Here we examine the "sea floor depth" variable. First, we pick it out using
# standard xarray syntax:
# 

# In[ ]:


xr_da = ds_mom6["deptho"]
xr_da


# Note that, in the output above, the `cell_methods` attribute indicates that the
# `areacello` variable contains the appopriate "cell area" for this variable (see
# under "Attributes"). However, this variable is not associated with the DataArray
# (see under "Coordinates").
# 
# Now we pick out the variable using the `cf` accessor and the appropriate
# standard name:
# 

# In[ ]:


cf_da = ds_mom6.cf["sea_floor_depth_below_geoid"]
cf_da


# Now notice that `areacello` is present under "Coordinates". This is because
# `cf_xarray` decodes CF metadata linking variables with each other (e.g.,
# `coordinates`, `cell_measures`, `ancillary_variables`).
# 
# As opposed to `xr_da`, `cf_da` extracted in the previous cell contains all
# `cell_measures` associated with the variable extracted.
# 

# In[ ]:


additional_coords = set(cf_da.coords) - set(xr_da.coords)
print("Cell measure extracted by cf_xarray:", additional_coords)


# ### Automagic plotting
# 
# `cf_xarray` automagically sets some optional keyword arguments for plotting
# functions. As opposed to `xarray`, in the example below `cf_xarray` assigns the
# appropriate coordinates to the plot axes (i.e., longitude and latitude).
# `cf_xarray` does so by parsing the `"coordinates"` attribute to identify the
# appropriate latitude and longitude variables:
# 

# In[ ]:


fig, (xr_ax, cf_ax) = plt.subplots(1, 2, figsize=(12, 4))

# left: xarray plot
xr_da.plot(ax=xr_ax)
xr_ax.set_title("xarray")

# right: cf_xarray plot
cf_da.cf.plot(ax=cf_ax)
cf_ax.set_title("cf_xarray")

fig.tight_layout()


# ### CF keys expansion
# 
# As mentioned above, the example dataset is characterized by multiple dimensions
# associated with the same spatial axes. Such information is decoded by
# `cf_xarray` and is used under the hood of wrapped functions. In the example
# below, the CF Axes keys (i.e., "X", "Y", and "Z") are expanded and multiple
# dimensions are sliced at once:
# 

# In[ ]:


ds_mom6_sliced = ds_mom6.cf.isel(
    X=slice(10), Y=slice(10), Z=slice(10), T=slice(10)
)
print("Original dataset sizes:", dict(ds_mom6.sizes))
print("  Sliced dataset sizes:", dict(ds_mom6_sliced.sizes))


# ## Use case examples
# 

# ### Seamlessly extract statistics from CF-compliant datasets
# 
# `cf_xarray` allows one to use the same code on a wide variety of CF-compliant
# objects that each has their own nomenclature.
# 
# There are two approaches to leveraging cf-xarray in applications.
# 
# #### Access variables through the `.cf` interface
# 
# In the example below, we define a function that uses many `cf_xarray` features,
# then we apply to objects with different dimension and coordinate names. All
# cf_xarray functionality is accessed using the `.cf` accessor
# 

# In[ ]:


def compute_top_10m_temp_anomaly(ds):
    """
    Compute the volume weighted temperature anomaly from the climatology.

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset to analyze

    Returns
    -------
    DataArray
    """

    # Compute and plot line
    with xr.set_options(keep_attrs=True):
        # Extract ocean potential temperature
        da = ds.cf["sea_water_potential_temperature"]
        # Fill missing cell volumes with zeros
        da = da.cf.assign_coords(volume=da.cf.coords["volume"].fillna(0))
        # Select temperature in the top 10m in 2003
        da = da.cf.sel(T="2003", Z=slice(0, 10))
        # Compute volume-weighted mean temperature
        da = da.cf.weighted("volume").mean(["X", "Y", "Z"])
        # Calculate an anomaly relative to the time mean
        da = da - da.cf.mean("T")

    # Update metadata
    da.attrs["standard_name"] += "_anomaly"
    da.attrs["long_name"] += " Anomaly"

    return da.squeeze(drop=True)


# Run the function on two different datasets and compare the results
compute_top_10m_temp_anomaly(ds_mom6).cf.plot(label="ds_mom6")
compute_top_10m_temp_anomaly(ds_cmip6_oce).cf.plot(label="ds_cmip6_oce")
_ = plt.legend()


# #### Standardize datasets using cf_xarray
# 
# Alternatively, `cf_xarray` provides utility functions to rename variables and
# dimensions in one object to match another object. Matching variables/dimensions
# are determined using CF metadata. One might choose to use this approach of
# standardizing datasets prior to passing them through a data processing pipeline.
# 
# Here we illustrate the `rename_like`
# [feature](https://cf-xarray.readthedocs.io/en/latest/generated/xarray.DataArray.cf.rename_like.html#xarray.DataArray.cf.rename_like).
# cf_xarray also supports renaming datasets through `.cf.rename`
# 

# In[ ]:


mom6_da = ds_mom6.cf["sea_water_potential_temperature"]
cmip6_da = ds_cmip6_oce.cf["sea_water_potential_temperature"]
renamed_mom6_da = mom6_da.cf.rename_like(cmip6_da)
print("        MOM6 dimensions:", mom6_da.dims)
print("       CMIP6 dimensions:", cmip6_da.dims)
print("renamed MOM6 dimensions:", renamed_mom6_da.dims)


# ### Regridding with xESMF
# 
# `cf-xarray` is used by [xESMF](https://pangeo-xesmf.readthedocs.io), a
# regridding package wrapping the powerful Fortran libray "ESMF". This example
# will regrid sea ice concentration data extracted from the CCCma-CanESM5 model, a
# participant in the CMIP6 experiment.
# 
# Our original sea ice data is on a tripolar grid. The target grid is a regular
# grid used by the NCEP reanalysis dataset that ships with Xarray. This regridding
# problem requires providing grid cell corners to xESMF in a specific
# CF-compatible format. Here we illustrate utility functions provided by cf_xarray
# to make this task easy and convenient.
# 

# In[ ]:


# Let's look at the grid shape itself and the data for one time step
fig, axs = plt.subplots(ncols=2, figsize=(12, 4))

# Notice how with .cf we will use the same keyword arguments
# Although here we explicitely pass the coordinate standard names for the plot axes,
# cf_xarray default scatter plot would produce the same results.
scatter_kwargs = dict(x="longitude", y="latitude", s=0.1)

# CMIP6: Input grid
ds_cmip6_ice.cf.plot.scatter(**scatter_kwargs, ax=axs[0])
axs[0].set_title(
    "The input horizontal grid points as seen on a lat/lon map."
    "\nOnly the northern hemisphere is shown."
)
axs[0].set_ylim(0, 90)

# NCEP: Target grid
ds_ncep.cf.plot.scatter(**scatter_kwargs, ax=axs[1])
axs[1].set_title("The target horizontal grid points")
axs[1].set_ylim(0, 90)

fig.tight_layout()


# We will regrid the sea ice data using the "conservative" method. Please refer to
# the xESMF documentation for details about the
# [different algorithms](https://xesmf.readthedocs.io/en/latest/notebooks/Compare_algorithms.html).
# The important information here is that the "conservative" regridding algorithms
# requires the grid points coordinates, but also the grid _corners_ coordinates.
# While the target grid doesn't provide them, they are easily computable with the
# help of `cf_xarray`:
# 

# In[ ]:


# Make a reasonable guess of the bounds of the spatial coordinates:
ds_ncep = ds_ncep.cf.add_bounds(["latitude", "longitude"])
ds_ncep


# This was easy since the grid is regular (i.e., latitude and longitude are 1D).
# Inferring bounds of 2D grids is
# [not yet supported](https://github.com/xarray-contrib/cf-xarray/issues/163) by
# cf-xarray.
# 
# Luckily, our sea ice data includes the corner coordinates in the
# `vertices_latitude` and `vertices_longitude` variables. However, xESMF expects a
# format that is different from the CF convention followed here. No worries,
# `cf_xarray` has an helper method just for this:
# 

# In[ ]:


# Get the bounds variable and convert them to "vertices" format
# Order=None, means that we do not know if the bounds are
# listed clockwise or counterclockwise, so we ask cf_xarray to try both.
lat_corners = cfxr.bounds_to_vertices(
    ds_cmip6_ice.vertices_latitude, "vertices", order=None
)
lon_corners = cfxr.bounds_to_vertices(
    ds_cmip6_ice.vertices_longitude, "vertices", order=None
)

# We are using special variable names "lon_b" and "lat_b" for easier detection by xESMF
ds_in = ds_cmip6_ice.assign(lon_b=lon_corners, lat_b=lat_corners)
ds_in


# Finally, the regridding is performed with xESMF. Under the hood, it uses
# cf-xarray to get the coordinates and their bounds, so we do not need to worry
# about renaming. The only exception is the input grid's corners, where we used
# hardcoded variable names because of the difference between xESMF's and CF's
# syntaxes for 2D grid bounds.
# 

# In[ ]:


# Regrid
regridder = xe.Regridder(ds_in, ds_ncep, "conservative")
sic_reg = regridder(ds_in.cf["sea_ice_area_fraction"])

# Plot the results
sic_reg.isel(time=0).plot()
_ = plt.title("Regridded sic data (Jan 2020)")


# # Conclusion
# 
# This notebook was a quick walkthrough of a few core cf-xarray features that
# enable scaling analysis pipelines _across_ datasets. It is fairly common for
# datasets to not have consistent terminology, and be imperfectly tagged with CF
# attributes. cf_xaray both allows you to leverage the presence of attributes, and
# provides utility functions to quickly fix imperfect tagging. For more see the
# [documentation](https://cf-xarray.readthedocs.io/en/latest/), specifically the
# [introductory notebook](https://cf-xarray.readthedocs.io/en/latest/examples/introduction.html).
# 

# # References
# 

# Adcroft, Alistair, Whit G Anderson, V Balaji, Chris Blanton, Mitchell Bushuk,
# Carolina O Dufour, John P Dunne, Stephen M Griffies, Robert Hallberg, Matthew J
# Harrison, Isaac M Held, Malte Jansen, Jasmin G John, John P Krasting, Amy R
# Langenhorst, Sonya Legg, Zhi Liang, Colleen McHugh, Aparna Radhakrishnan,
# Brandon G Reichl, Anthony Rosati, Bonita L Samuels, Andrew Shao, Ronald J
# Stouffer, Michael Winton, Andrew T Wittenberg, Baoqiang Xiang, Niki Zadeh, and
# Rong Zhang, October 2019: The GFDL Global Ocean and Sea Ice Model OM4.0: Model
# Description and Simulation Features. Journal of Advances in Modeling Earth
# Systems, 11(10), DOI:10.1029/2019MS001726.
# 
# Cherian, Deepak, Almansi, Mattia, Bourgault, Pascal, keewis, Kent, Julia, Thyng,
# Kristen, … Chauhan, Subhendra Singh. (2021, May 11). xarray-contrib/cf-xarray:
# (Version v0.5.2). Zenodo. http://doi.org/10.5281/zenodo.4749736
# 
# Hoyer, Stephan, Hamman, Joe, Roos, Maximilian, keewis, Cherian, Deepak,
# Fitzgerald, Clark … Bovy, Benoit. (2021, May 19). pydata/xarray: v0.18.2
# (Version v0.18.2). Zenodo. http://doi.org/10.5281/zenodo.4774304
# 
# Kalnay, E., Kanamitsu, M., Kistler, R., Collins, W., Deaven, D., Gandin, L.,
# Iredell, M., Saha, S., White, G., Woollen, J., Zhu, Y., Chelliah, M., Ebisuzaki,
# W., Higgins, W., Janowiak, J., Mo, K. C., Ropelewski, C., Wang, J., Leetmaa, A.,
# Reynolds, R., Jenne, R., & Joseph, D. (1996). The NCEP/NCAR 40-Year Reanalysis
# Project, Bulletin of the American Meteorological Society, 77(3), 437-472.
# 
