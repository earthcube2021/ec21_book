#!/usr/bin/env python
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc">
#   <ul class="toc-item">
#     <li>
#       <span><a href="#WH_01_Spatio-Temporal_Deep_Analogs_2021" data-toc-modified-id="Template-Notebook-for-EarthCube---Long-Version-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Weiming Hu's Notebook for EarthCube 2021: Machine Learning Guided Weather Analogs</a></span>
#       <ul class="toc-item">
#         <li><span><a href="#Authors" data-toc-modified-id="Author(s)-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>Authors</a></span></li>
#         <li><span><a href="#Purpose" data-toc-modified-id="Purpose-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Purpose</a></span></li>
#         <li><span><a href="#Technical-Contributions" data-toc-modified-id="Technical-contributions-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>Technical Contributions</a></span></li>
#         <li><span><a href="#Methodology" data-toc-modified-id="Methodology-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>Methodology</a></span></li>
#         <li><span><a href="#Data" data-toc-modified-id="Data-1.45"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>Data</a></span></li>
#         <li><span><a href="#Results" data-toc-modified-id="Results-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>Results</a></span></li>
#         <li><span><a href="#Funding" data-toc-modified-id="Funding-1.6"><span class="toc-item-num">1.6&nbsp;&nbsp;</span>Funding</a></span></li>
#         <li><span><a href="#Keywords" data-toc-modified-id="Keywords-1.7"><span class="toc-item-num">1.7&nbsp;&nbsp;</span>Keywords</a></span></li>
#         <li><span><a href="#Citation" data-toc-modified-id="Citation-1.8"><span class="toc-item-num">1.8&nbsp;&nbsp;</span>Citation</a></span></li>
#         <li><span><a href="#Acknowledgements" data-toc-modified-id="Acknowledgements-1.11"><span class="toc-item-num">1.11&nbsp;&nbsp;</span>Acknowledgements</a></span></li>
#       </ul>
#     </li>
#     <li>
#       <span><a href="#Setup" data-toc-modified-id="Setup-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Setup</a></span>
#     </li>
#     <li><span><a href="#Parameter-Definitions" data-toc-modified-id="Parameter-definitions-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Parameter Definitions</a></span></li>
#     <li><span><a href="#Data-Import" data-toc-modified-id="Data-import-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Data Import</a></span></li>
#     <li>
#       <span><a href="#Data-Analysis" data-toc-modified-id="Data-processing-and-analysis-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>Data Analysis</a></span>
#       <ul class="toc-item">
#         <li><span><a href="#Ensemble-Visualization" data-toc-modified-id="The-10-rules-5.1"><span class="toc-item-num">5.1&nbsp;&nbsp;</span>Ensemble Visualization</a></span></li>
#         <li><span><a href="#Latent-Feature-Visualization" data-toc-modified-id="The-10-rules-5.2"><span class="toc-item-num">5.1&nbsp;&nbsp;</span>Latent Feature Visualization</a></span></li>
#       </ul>
#     </li>
#     <li><span><a href="#Summary-and-Future-Work" data-toc-modified-id="Summary-6"><span class="toc-item-num">6&nbsp;&nbsp;</span>Summary and Future Work</a></span></li>
#     <li><span><a href="#References" data-toc-modified-id="References-7"><span class="toc-item-num">7&nbsp;&nbsp;</span>References</a></span></li>
#   </ul>
# </div>
# 
# 

# # WH_01_Spatio-Temporal_Deep_Analogs_2021

# ## Authors
# 
# - Author1 = {"name": "Weiming Hu", "affiliation": "Department of Geography, The Pennsylvania State University", "email": "weiming@psu.edu", "orcid": "0000-0003-4501-1435"}
# - Author2 = {"name": "Guido Cervone", "affiliation": "Department of Geography, The Pennsylvania State University", "email": "cervone@psu.edu"}
# - Author3 = {"name": "George S. Young", "affiliation": "Department of Meteorology and Atmospheric Science, The Pennsylvania State University", "email": "g3y@psu.edu"}
# - Author4 = {"name": "Luca Delle Monache", "affiliation": "Scripps Institution of Oceanography, University of California, San Diego", "email": "ldellemonache@ucsd.edu"}

# ## Purpose
# 
# This notebook presents the **Spatio-Temporal Deep Analog**, an [Analog Ensemble](https://weiming-hu.github.io/AnalogsEnsemble/) (AnEn) technique with Machine Learning powered weather similarity metric.
# 
# It, for the first time, integrates Machine Learning techniques with the Analog Ensemble technique. While the Analog Ensemble uses a multivariate distance function to define weather similarity, it has several limitation regarding weight optimization and model updates in the search history. This notebook showcases how a neural network can be trained and used as a powerful weather similarity metric that leads to an increased flexibility and an improved accuracy compared to the original weather similarity.
# 
# This work has been inspired by the recent progress in Computer Vision, especially in face recognition and identification like the [FaceNet](https://www.cv-foundation.org/openaccess/content_cvpr_2015/html/Schroff_FaceNet_A_Unified_2015_CVPR_paper.html) from Google.

# ## Technical Contributions
# 
# The development of Deep Analog (DA) spans across multiple years of research and several projects. A summary of the technical contributions is provided below:
# 
# 1. [Parallel Analog Ensemble](https://weiming-hu.github.io/AnalogsEnsemble): This open-source project aims to make the generation of weather analogs faster and easier. It is implemented in C++ and it provides documented interface for development.
# 2. [RAnEn](https://weiming-hu.github.io/AnalogsEnsemble/R/): The R package is implemented on top of the PAnEn. This package is mainly developed for researchers familiar with R and who want to quickly prototype a research idea.
# 3. [PyAnEn](https://github.com/Weiming-Hu/PyAnEn): The Python package is developed mainly for post processing predictions generated from the PAnEn. It aims at providing the functionality to analyze and verify large weather prediction data on supercomputers utilizing distributed computing and parallel processing.
# 4. [Deep Analog](https://github.com/Weiming-Hu/DeepAnalogs): The Python package is developed for building and training a neural network to be used as the weather similarity metric with the PAnEn.

# ## Methodology
# 
# ### Analog Ensemble
# 
# Analog Ensemble is a technique to generate ensemble forecasts from a single run of deterministic weather model and the corresponding historical observations. It has several major advantages:
# 
# 1. **Computation**: Only a single run of a deterministic model is needed to generate ensembles.
# 1. **Uncertainty Quantification**: The spread of the generated ensemble is well correlated with the prediction accuracy, providing a good estimate of forecast uncertainty.
# 1. **Decoupled Prediction**: The AnEn can be extended to predict variables that are not directly modeled by weather models, as long as an observation history is obtained, e.g., the power generation from a solar farm.
# 
# It is assumed that, given a static weather model, similar forecasts are associated with similar errors. Therefore, by using the associated historical observations for the most similar forecasts, the error of the target forecast can be effectively corrected. The key issue is how we can identify good-quality weather analogs.
# 
# The core to Analog Ensemble is a weather similarity metric:
# 
# $\left \| F_{t},A_{{t}'} \right \| = \sum_{i=1}^{N_{v}} {\frac{\omega_{i}}{\sigma_{f_{i}}} \sqrt{\sum_{j = -\tilde{t}}^{\tilde{t}} {\left ( F_{i,t+j}-A_{i,{t}'+j} \right )^{2}}}}$,
# 
# where:
# 
# - $F_{t}$ is the model prediction valid at the weather model initialization time stamp $t$ at a specific location and lead time.
# - $A_{{t}'}$ is the historical repository of the weather model from the search space at the same location and lead time, but with a different model initialization time ${t}'$.
# - $N_{v}$ is the number of physical variables used during forecast similarity calculation.
# - $\omega_{i}$ is the weight for each physical variable which suggests the relevant importance of the physical variable with respect to the others.
# - $\sigma_{f_{i}}$ is the standard deviation for the physical variable $i$ calculated from the historical forecasts at the same location and Forecast Lead Time (FLT).
# - $\tilde{t}$ equals to half of the time window size of the FLTs to be compared so that weather analogs are identified within a very small time window.
# - $F_{i,t+j}$ is the value of the current forecast for the physical variable $i$ at the valid time $t+j$.
# - $A_{i,{t}'+j}$ is the value of the historical forecast for the physical variable $i$ at the valid time ${t}'+j$.
# 
# <img src="figures/AnEn-schema.png" alt="AnEn Schema" style="width: 600px;"/>
# 
# The above sketch provides a pictorial representation of generating a four-member forecast ensemble.
# 
# - The top arrow represents an operational run of a deterministic weather model with the grey shaded area representing the history.
# - The bottom arrow represents the historical observations of the variable of prediction, associated with the model. For example, the bottom arrow would represents surface temperature records if the predictand is temperature. Since it is the historical archive of observations, it only overlaps with the grey shaded area from the first arrow.
# 
# Generating a four-member ensemble with the AnEn consists of four steps:
# 
# 
# 1. A multi-variate target forecast is retrieved from the deterministic model.
# 1. The similarity measure is calculated between the target forecast and each of the historical forecasts based on the weather similarity metric. Four candidate forecasts with the highest similarity (the lowest distances) are identified as weather analogs to the target forecasts.
# 1. Observations associated with the four candidate forecasts are retrieved.
# 1. The historical observations become forecast ensemble members in the final analog ensemble.
# 
# However, as many scientific toolbox, Analog Ensemble has its own challenges and limitation:
# 
# 1. **Grid-to-Grid Comparison**: Weather similarity is defined on a single grid basis, meaning that forecasts are compared independently at each grid. It limits the weather analogs from detecting spatial patterns.
# 1. **Model Updates**: An important assumption of the Analog Ensemble is that similar weather forecasts have similar forecast errors. Thus, using the observations can effectively correct for this error. However, operational models are constantly changing and upgrading. Model behaviors are subject to changes across time which could negatively affect the correction.
# 1. **Weight Optimization**: Predictor weights, $\omega_{i}$, needs to be determined as priori. This is usually done via extensive grid search which is computational expensive and error-prone. It also limits the number of predictors used from the weather model to only a few while, in practice, many more variables are modeled by weather models.

# ### Deep Analog
# 
# To address the above challenges, we propose a renovated weather similarity metric by using Machine Learning, specifically a Convolutional Long Short-Term Memory neural network as the weather similarity metric.
# 
# <img src="figures/DA_Structure.png" alt="ConvLSTM" style="width: 1200px;"/>

# The above figure shows the structure of a particular embedding network designed for [WRF NAM NMM](https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/north-american-mesoscale-forecast-system-nam).
# 
# - The input is a four-dimensional data structure with [height x width x variables x Time window]. The color coded image is an example of the solar irradiance field from the weather model.
# - The embedding model has two convolutional layers and one convolutional LSTM layer, each one followed by a max-pooling layer.
# - The output of the embedding network is a one-dimensional vector with a length of 120 values, referred to as the latent features or the hidden features.
# 
# The embedding network is applied to the entire forecast dataset to transform all forecast variables into a latent space, specifically a 120-dimensional latent space. Finally, forecast clusters, calculated based on the Euclidean distance, identified in the latent space is treated as weather analogs.

# ### Reverse Analog as Training
# 
# <img src="figures/ReverseAnalogs.png" alt="Reverse Analog" style="width: 700px;"/>
# 
# How do we train an effective Machine Learning similarity metric?
# 
# We propose the Reverse Analog procedure in order to train an effective embedding network.
# 
# We **ARE NOT**:
# 
# 1. training a neural network that directly predict a variable of interest;
# 1. training a compression neural network that only seeks to compress the original forecasts with fewer variables.
# 
# Instead, we **ARE**:
# 
# 1. training a neural network that transform the forecasts into a latent space so that analogs can be more effectively found with better accuracy;
# 1. training a neural network that learns a relationship between a spatial forecasts and the corresponding observations;
# 1. training a neural network that, by looking at the forecasts, it learns whether the observations are actually similar.
# 
# The training process is illustrated below:
# 
# 1. Given an observation, e.g., the observed solar irradiance, from a target date, find the historical dates that have the most similar observations.
# 1. The associated forecasts on those most similar dates are queried.
# 1. The network is trained to produce a latent space that places similar forecasts closer to each other, dissimilar forecasts further away from each other. Specifically, the blue pairs are more similar than the grey pairs.
# 
# The core idea is that, the neural network is trained with **both forecasts and observations** in order to generate effective embeddings.

# ## Data
# 
# A summmar of the research data is provided below:
# 
# - **Weather Model**: [GFS](https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs) with 3 different resolutions, 1°, 0.50°, and 0.25°, and [WRF NAM NMM](https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/north-american-mesoscale-forecast-system-nam) with 12 km resolution. Results from WRF NAM NMM is the focus of this notebook.
# - **Observation**: Solar irradiance reaching the surface. Observations are collectecd from the [Surface Radiation Budget Network](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00540).
# - **Geographic Location**: Center Pensylvania
# - **Prediction Period**: 2018/01/01 to 2019/10/31
# - **Search Period**: 2015/01/01 to 2017/12/31

# ## Results
# 
# <img src="figures/SingleStation_MAE_Bias_CRPS.png" alt="ConvLSTM" style="width: 700px;"/>

# The above figure shows the prediction error of solar irradiance from various configurations of Analog Ensemble and weather models. Specifically,
# 
# - `NWP` stands for Numerical Weather Prediction, namely the uncalibrated prediction from the numerical weather prediction model.
# - `AnEn Spatial` stands for the Analog Ensemble with nearby forecasts as additional predictors.
# - `AnEn` stands for the conventional Analog Ensemble with forecasts on a single grid.
# - `DA IS` stands for the Deep Analog with no convolutional layers (spatial information) and trained on forecasts from only one grid.
# - `DA SSE` stands for the Deep Analog with no convolutional layers (spatial information) and trained on forecasts from nearby grids.
# - `DA Spatial` stands for the Deep Analog with convolutional layers (spatial information).
# 
# The verification metrics include Mean Absolute Error (MAE), Bias, and Continuous Ranked Probability Score (CRPS). CRPS can be understood as the ensemble version of Root-Mean-Square Error (RMSE). It is calculated using the package [properscoring](https://pypi.org/project/properscoring/).
# 
# Several key observations:
# 
# 1. AnEn and its variants all outperforms the baseline forecasts demonstrating its capability in correcting weather models.
# 2. Comparing `AnEn Spatial` and `AnEn`, it suggests that simply using forecasts from nearby locations as additional predictors is not an effective way of exploiting spatial information.
# 3. `DA Spatial` outperforms `AnEn` demonstrating the added capability to correct weather forecasts with a Machine Learning simialrity metric.
# 4. `DA Spatial` outperforms other variants of DA, suggesting the convolutional layers are essential to extract spatial information.

# <img src="figures/MultiStation_NAM_NMM.png" alt="ConvLSTM" style="width: 900px;"/>

# The above figures visualize the spatial error when searching nearby forecasts. The question they try to address is, can you find better analogs with nearby points. By applying a spatial mask. `DA Spatial shows` a much smoother error surface when nearby forecasts are compared. In the other cases, forecasts are compared on a single grid, therefore the error surface is more disjointed. The red `x` indicates the best search grid, meaning that if you compare the analog selected across all the grids in the domain, they will be the best grids to search.   
# 
# Specifically,
# 
# - `AnEn` and `DA IS` both do not exploit spatial information. Prediction error drastically increases when forecasts from distant regions are used.
# - `DA SSE` tries to train an embedding model using data from nearby forecasts, but it is shown that simply adding to the training data without a change in the model architecture is not enough to achieve effective improvement in prediction accuracy.
# - `DA Spatial`, by far, shows the best prediction accuracy. The smooth error surface centered at the close vincinity of the location of interest suggest its capability of capturing the spatial features.

# <img src="figures/MAE_vs_Distance.png" alt="ConvLSTM" style="width: 800px;"/>

# The above figure assesses the effectiveness of different similarity metrics when a distant forecasts is used. It is desired to observe a trend that, while the error should increase as the compared forecasts comes from a more distant location, the speed of the increase is low so that the AnEn can benefit from searching nearby locations.
# 
# In all four panels of the above figure, `DA Spatial` (shown in red) always has the lowest prediction error across all bins whereas the alternatives typically show a faster increase of error when nearby forecasts are compared. This is due to the grid-by-grid comparison without any spatial information. As for `DA Spatial`, the nearby forecasts over a spatial mask still appears largely similar to the target forecast when there is only a small geographic offset. This is similar to the sliding windows at two consecutive steps.

# ## Funding
# 
# Include references to awards that supported this research. Add as many award references as you need.
# 
# - Award1 = {"agency": "US National Science Foundation", "award_code": "1639707", "award_URL": "https://www.nsf.gov/awardsearch/showAward?AWD_ID=1639707"}
# - Award2 = {"agency": "Department of Geography, The Pennsylvania State University"}
# 
# ## Keywords
# 
# keywords=["Ensemble Forecast", "Analog Ensemble", "Machine Learning", "Renewable Energy", "Solar Irradiance"]
# 
# ## Citation
# 
# Hu, W., Cervone, G., Young, G., Delle Monache, L., 2021. Machine Learning Guided Weather Analogs. Jupyter Notebook. Accessed on 2021/6/11 from https://github.com/Weiming-Hu/EarthCube2021
# 
# ## Acknowledgements
# 
# We thank the support from the [RADICAL group](http://radical.rutgers.edu/) at Rutgers Univeristy through their [Ensemble Toolkit](https://radical-cybertools.github.io/entk/index.html), which makes generating analog ensembles on supercomputers possible.
# 
# The notebook is licensed under a <a href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License.</a>

# # Setup
# 
# **Only run the following cell if you are on Binder!**
# 
# Some packages are not properly set up on Binder. Therefore, I need to manually fix them. The following cell only needs to be run **once**.

# In[ ]:


get_ipython().run_cell_magic('bash', '', 'pip uninstall pillow -y\npip install pillow')


# The following libraries are needed to run the cells. Please see comments for their purposes.

# In[ ]:


# Numpy is a very common package for handling high-dimensional data structure.
import numpy as np

# Pandas is a very common package for dealing with tables.
import pandas as pd

# PyAnEn is the package to help reading analog results.
import PyAnEn.IO as AnEnIO

# The script provides mostly visualization functions.
from functions import *

# Autoreload extension. This force the session to reload the package in case of code changes.
if 'autoreload' not in get_ipython().extension_manager.loaded:
    get_ipython().run_line_magic('load_ext', 'autoreload')
    
get_ipython().run_line_magic('autoreload', '2')


# # Parameter Definitions

# In[ ]:


# This is the NetCDF file that has DA predictions.
da_file = 'data/DA_Spatial.nc'

# This is the NetCDF file that has AnEn predictions.
anen_file = 'data/AnEn.nc'

# This is the NetCDF file that has WRF NAM NMM predictions.
#
# There were original 48 parameters. Due to file size limit, only 2 parameters have been kept in the sample data file.
#
# They are
#
# 1. downward shortwave radiation flux [260087_0_surface_dswrf]
# 2. total cloud cover [228164_0_atmosphere_tcc].
#
fcst_file = 'data/NAM_NMM.nc'

# This is the observations for solar irradiance.
obs_file = 'data/SURFRAD_PennState.nc'

# This is the latent features for WRF NAM NMM forecasts.
embedding_file = 'data/NAM_NMM_Spatial.nc'


# # Data Import
# 
# The rest of the notebook visualizes already generated analog ensembles. Training the deep network and generating analog ensembles are skipped because of the computational cost. Typically, I used [Deep Analog](https://github.com/Weiming-Hu/DeepAnalogs) to train the network and used [Parallel Analog Ensemble](https://weiming-hu.github.io/AnalogsEnsemble/) to generate analog ensembles. Please feel free to contact Weiming if you are interested in integrating the methodology to your scientific workflow.
# 
# Please run the following cell to download a sample of the data.

# In[ ]:


get_ipython().run_cell_magic('bash', '', 'if test -f "EarthCube2021.tar.gz"; then\n    echo "Tarball exists. Skip downloading."\nelse\n    echo "Downloading data tarball from https://prosecco.geog.psu.edu/EarthCube2021/EarthCube2021.tar.gz ..."\n    wget https://prosecco.geog.psu.edu/EarthCube2021/EarthCube2021.tar.gz 1> log 2> log\nfi\n\ntar xf EarthCube2021.tar.gz\necho "Done!"')


# In the following cell, data are loaded into memory. Please read the comments for detailed explanations.

# In[ ]:


da_ds = AnEnIO.open_dataset(da_file, decode=True)


# In[ ]:


da_ds


# Internally, `PyAnEn.IO` uses [xarray](http://xarray.pydata.org/en/stable/) which is a labeled data structure for NetCDF data. It is very intuitive to use. The above cell prints a description of the loaded data. The most important information shown above is:
# 
# - `analogs`: A four-dimensional array for analog predictions generated by Deep Analog. There are 21 ensemble member and 647 test times. There are three lead times, but all our subsequent analysis is only for the second lead time, as an example. There is only one station, which is at the center Pennsylvania.
# 
# Please note the naming convention: the first part, `da` stands for the prediction technqiue, and the second part, `ds` standards for `xarray` dataset.
# 
# The following code snippet plots the location of the station as an example to show how one would typically access the data in an `xarray` dataset.

# In[ ]:


fig, ax = plt.subplots(1, 1, figsize=(10, 5), subplot_kw={'projection': cartopy.crs.PlateCarree()})

ax.add_feature(cartopy.feature.OCEAN, linewidth=1)
ax.add_feature(cartopy.feature.LAKES, linewidth=1)
ax.add_feature(cartopy.feature.STATES, linewidth=1)

ax.scatter(da_ds['Xs'], da_ds['Ys'], marker='*', c='red', s=100, label='SURFRAD')
ax.legend(loc='upper right')
fig.show()


# Similarly, we load the other datasets.

# In[ ]:


# This dataset stores analog predictions from AnEn
anen_ds = AnEnIO.open_dataset(anen_file, decode=True)

# This dataset stores WRF NAM NMM forecasts
fcst_ds = AnEnIO.open_dataset(fcst_file, decode=True)

# This dataset stores SURFRAD observations
obs_ds = AnEnIO.open_dataset(obs_file, decode=True)

# This dataset stores the transformed WRF NAM NMM forecasts, namely the 120 latent features generated from DA Spatial
embeddings_ds = AnEnIO.open_dataset(embedding_file, decode=True).squeeze('num_stations').sel(num_flts=3600 * 18)


# In[ ]:


print('There are {} features in the transformed latest space for WRF NAM NMM.'.format(embeddings_ds.dims['num_parameters']))


# These latent features are generated by running each individual forecasts into a trained Deep Analog embedding network.
# 
# Since the training of such a network and the generation of the embeddings are very time consuming, the results are simply provided here and we will only be looking at forecasts at 1 PM which typically correspond to the highest solar irradiance throughout the day.

# # Data Analysis
# 
# ## Ensemble Visualization
# 
# In this section, we focus on visualizing ensembles generated from `AnEn` and `DA Spatial`.
# 
# Several key differences between `AnEn` and `DA Spatial` to keep in mind while reading the figures:
# 
# 1. `AnEn` determines similar weather by looking at WRF NAM NMM forecasts at a single grid location.
# 1. `DA Spatial` determines similar weather by looking at WRF NAM NMM forecasts with a spatial domain. The network is also trained with observations so that `DA Spatial` is designed to learn what the observations would be given a spatial forecasts and trying to cluster forecasts based on the associated observations.

# In[ ]:


# Select a date for which the forecast ensemble will be plotted.
#
# The date must be between 2018/01/01 and 2019/10/31.
#
target_date = pd.to_datetime('2019/05/05')


# In[ ]:


# If you want to extract the dates of similar weather from AnEn, use this code
# anen_times = get_analog_forecast_times(anen_ds, time=target_date)[:21]

# If you want to extract the dates of similar weather from DA Spatial, use this code
anen_times = get_analog_forecast_times(da_ds, time=target_date)[:21]

################################################################################
# 21 is the number of analog members, namely the number of most similar dates. #
# It is better not to change this number as it might affect visualization.     #
################################################################################

print('{} is most similar to the following dates:'.format(target_date.strftime('%Y/%m/%d')))
print(np.sort([dt.strftime('%Y/%m/%d') for dt in anen_times]))


# In[ ]:


# Combine the target date and all the similar dates. The first date is plotted in a slightly bigger panel.
plot_times = np.concatenate([np.array([target_date]), anen_times])

plot_analog_forecasts(obs_ds, fcst_ds, plot_times, obs_ds)


# The above figure shows the following information:
# 
# - Each panel is a WRF NAM NMM forecasts from a particular date labeled at the uppper left of each panel.
# - The observed solar irradiance is shown following the date.
# - Points are color coded with WRF NAM NMM predictions for solar irradiance.
# - `AnEn` determines weather similarity based on a single grid point from WRF NAM NMM (shown in star). However, `DA Spatial` determines weather similarity based on a spatial domain and the associated observations.

# ## Latent Feature Visualization
# 
# In this section, we visualize the latent features in its relationship with the original forecast variable. We try to see correlation between the learned features and the original weather variables.

# In[ ]:


# Define the time range to plot.
# We are plotting the entire test period.
#
datetimes = pd.date_range('2018/01/01', '2019/10/31')

# Some dates are missing from the WRF NAM NMM, so we need to make sure only plot dates that exist.
datetimes = np.array([dt for dt in datetimes if dt in embeddings_ds['Times']])

# Query the embeddings
embeddings = embeddings_ds['Data'].sel(num_times=datetimes).data.transpose()

print('You have extracted embeddings on {} dates.'.format(embeddings.shape[1]))
print('Each embedding has {} latent features.'.format(embeddings.shape[0]))


# In[ ]:


# Select which latent feature to use when creating the scatter plot
latent_feature = 115

# The point index of cluster centers to create in the scatter plot.
# The visualization function supports three cluster centers.
#
center1 = 170
center2 = 190
center3 = 10

# The number of nearest neighbors to visualize
n_nbs = 10


# In[ ]:


fig = plot_feature(latent_feature, center1, center2, center3, n_nbs, embeddings, datetimes, fcst_ds, obs_ds)


# The above figure shows the following information:
# 
# - The left panel shows the scatter plot of the selected latent feature with respect to the observed solar irradiance. It is a good visual diagnosis of the correlation between a single latent feature and the predictand.
# - Each point in the scatter plot represent a embedding WRF NAM NMM forecast on a particular date.
# - Three clusters are created based on the scatter plot, controlled by `center1`, `center2`, and `center3`.
# - Panels on the right are WRF NAM NMM associated with the clustered points. Please note the corresponding labels in the legend.

# # Summary and Future Work
# 
# This notebook has proposed the Deep Analog, a weather similarity metric driven by a Conv-LSTM neural network. Results show that Deep Analog outperforms Analog Ensemble when predicting solar irradiance. The integration of spatial information allows Deep Analog to be associated with a smoother error surface. Weight optimization is also significantly faster while, with Analog Ensemble, weight optimization would be almost impossible for more than a dozen predictors.
# 
# We suggest future research should further investigate the efficacy of such a Machine Learning guided technique on other meteorological variables, for example wind. Analog forecasting is also helpful in downscaling climate and reginal models. It would be interesting to study how this architecture can be used for model downscaling.
# 
# # References
# 
# 1. Delle Monache, Luca, et al. "Probabilistic weather prediction with an analog ensemble." Monthly Weather Review 141.10 (2013): 3498-3516.
# 1. Hu, Weiming, and Guido Cervone. "Dynamically Optimized Unstructured Grid (DOUG) for Analog Ensemble of numerical weather predictions using evolutionary algorithms." Computers & Geosciences 133 (2019): 104299.
# 1. Hu, Weiming, et al. "PARALLEL ANALOG ENSEMBLE–THE POWER OF WEATHER ANALOGS." NCAR Technical Notes NCAR/TN-564+ PROC: 1.
# 1. Schroff, Florian, Dmitry Kalenichenko, and James Philbin. "Facenet: A unified embedding for face recognition and clustering." Proceedings of the IEEE conference on computer vision and pattern recognition. 2015.
# 1. Shi, Xingjian, et al. "Convolutional LSTM network: A machine learning approach for precipitation nowcasting." arXiv preprint arXiv:1506.04214 (2015).
# 1. Hu, Weiming, et al. "Weather Analogs with a Machine Learning Similarity Metric for Renewable Resource Forecasting." arXiv preprint arXiv:2103.04530 (2021).
