#!/usr/bin/env python
# coding: utf-8

# # Automated Machine Learning for Earth Science via AutoGluon

# ## Authors

# - Author1 = {"name": "Xingjian Shi", "affiliation": "Amazon Web Services", "email": "xjshi@amazon.com", "orcid": ""}
# - Author2 = {"name": "Wen-ming Ye", "affiliation": "Amazon Web Services", "email": "wye@amazon.com", "orcid": ""}
# - Author3 = {"name": "Nick Erickson", "affiliation": "Amazon Web Services", "email": "neerick@amazon.com", "orcid": ""}
# - Author4 = {"name": "Jonas Mueller", "affiliation": "Amazon Web Services", "email": "jonasmue@amazon.com", "orcid": ""}
# - Author5 = {"name": "Alexander Shirkov", "affiliation": "Amazon Web Services", "email": "ashyrkou@amazon.com", "orcid": ""}
# - Author6 = {"name": "Zhi Zhang", "affiliation": "Amazon Web Services", "email": "zhiz@amazon.com", "orcid": ""}
# - Author7 = {"name": "Mu Li", "affiliation": "Amazon Web Services", "email": "mli@amazon.com", "orcid": ""}
# - Author8 = {"name": "Alexander Smola", "affiliation": "Amazon Web Services", "email": "alex@smola.org", "orcid": ""}

# ## Table of Contents
# * [Purpose](#purpose)
# * [Setup](#setup)
# * [Forest Cover Type Classification](#forest-cover-type-classification)
#     * [Train Model with One Line](#train-model-with-one-line)
#     * [Evaluation and Prediction](#evaluation-and-prediction)
#     * [Load the Predictor](#load-the-predictor)
#     * [Feature Importance](#feature-importance)
#     * [Achieve Better Performance](#achieve-better-performance)
# * [Solar Radiation Prediction](#solar-radiation-predictions)
# * [More Information](#more-information)

# ## Purpose

# In this notebook, we introduce [AutoGluon](https://github.com/awslabs/autogluon) to the Earth science community. AutoGluon is an automated machine learning toolkit that enables users to solve machine learning problems with a single line of code. Many earth science problems involve tabular-like datasets. With AutoGluon, you can feed in the **raw** data table and specify the `label` column. AutoGluon will deliver a model that has reasonable performance in a short period of time. In addition, with AutoGluon, you can also analyze the importance of each feature column with a single line of code. In the following, we illustrate how to use AutoGluon to build machine learning models for two Earth Science problems.

# ## Setup
# 
# We have pre-installed [AutoGluon](https://github.com/awslabs/autogluon) via pip. Here, we will fix the random seed.

# In[1]:


# Uncomment below to install autogluon
# !python3 -m pip install autogluon
import random
import numpy as np
random.seed(123)
np.random.seed(123)


# ## Forest Cover Type Classification

# In the first example, we will predict the forest cover type (the predominant kind of tree cover) from strictly cartographic variables. The dataset is downloaded from [Kaggle Forest Cover Type Prediction](https://www.kaggle.com/c/forest-cover-type-prediction). Study area of the dataset includes four wilderness areas located in the Roosevelt National Forest of northern Colorado. The actual forest cover type for a given 30 x 30 meter cell was determined from US Forest Service (USFS) Region 2 Resource Information System data. Independent variables were then derived from data obtained from the US Geological Survey and USFS. The data is in raw form and contains binary columns of data for qualitative independent variables such as wilderness areas and soil type. Let's first download the dataset.

# In[2]:


get_ipython().system('wget https://deep-earth.s3.amazonaws.com/datasets/earthcube2021_demo/forest-cover-type-prediction.zip -O forest-cover-type-prediction.zip')
get_ipython().system('unzip -o forest-cover-type-prediction.zip -d forest-cover-type-prediction')


# Here, we load and visualize the dataset. We will split the dataset to 80% training and 20% development for the purpose of reporting the score on the development data. Also, for the purpose of demonstration, we will subsample the dataset to 5000 samples.

# In[3]:


import pandas as pd
from sklearn.model_selection import train_test_split
df = pd.read_csv('forest-cover-type-prediction/train.csv.zip')
df = df.drop('Id', 1)
df = df.sample(5000, random_state=100)
train_df, dev_df = train_test_split(df, random_state=100)


# By visualizing the dataset, we can see that there are 54 feature columns and 1 label column called `"Cover_Type"`.

# In[4]:


train_df.head(5)


# ### Train Model with One Line

# Next, we train a model in AutoGluon with a single line of code. We will just need to specify the label column before calling `.fit()`. Here, the label column is `Cover_Type`. AutoGluno will inference the problem type automatically. In our example, it can correctly figure out that it is a "multiclass" classification problem and output the model with the best accuracy. Internally, it will also figure out the feature type automatically.

# In[5]:


import autogluon
from autogluon.tabular import TabularPredictor
predictor = TabularPredictor(label='Cover_Type', path='ag_ec2021_demo').fit(train_df)


# We can visualize the performance of each model with `predictor.leaderboard()`. Internally, AutoGluon trains a diverse set of different tabular models and computes a weighted ensemble to combine these models.

# In[6]:


predictor.leaderboard()


# ### Evaluation and Prediction

# We can also evaluate the model performance on the heldout predictor dataset by calling `.evaluate()`.

# In[7]:


predictor.evaluate(dev_df)


# To get the prediction, you may just use  `predictor.predict()`.

# In[8]:


predictions = predictor.predict(dev_df)
predictions


# For classification problems, we can also use `.predict_proba` to get the probability.

# In[9]:


probs = predictor.predict_proba(dev_df)
probs.head(5)


# ### Load the Predictor

# Loading a AutoGluon model is straight-forward. We can directly call `.load()`

# In[10]:


predictor_loaded = TabularPredictor.load('ag_ec2021_demo')
predictor_loaded.evaluate(dev_df)


# ### Feature Importance

# AutoGluon offers a built-in method for calculating the relative importance of each feature based on [permutation-shuffling](https://scikit-learn.org/stable/modules/permutation_importance.html). In the following, we calculate the feature importance and print the top-10 important features. Here, `importance` means the importance score and the other values give you an understanding of the statistical significance of the calculated score.

# In[11]:


importance = predictor.feature_importance(dev_df, subsample_size=500)
importance.head(10)


# From the results, we can see that `Elevation` is the most important feature. `Horizontal_Distance_To_Roadways` is the 2nd most important feature.

# ### Achieve Better Performance

# The default behavior of AutoGluon is to compute a weighted ensemble of a diverse set of models. Usually, you can achieve better performance via stack ensembling. To achieve better performance based on automated stack ensembling, you can specify `presets="best_quality"` when calling `.fit()` in AutoGluon. For more details, you can also checkout our provided script. The detailed architecture is described in [1] and we also provide the following figure so you can know the general architecture.
# 
# <img src="https://deep-earth.s3.amazonaws.com/datasets/earthcube2021_demo/stacking.png" alt="screenshot" style="width: 500px;"/>
# 
# With `.fit(train_df, presets="best_quality")`, we are able to achieve 82/1692 in the competition. To reproduce our number, you may try the command mentioned in [link](https://github.com/sxjscience/EC2021_autogluon_notebook).
# 
# <img src="https://deep-earth.s3.amazonaws.com/datasets/earthcube2021_demo/forest_cover_type.png" alt="screenshot" style="width: 500px;"/>

# ## Solar Radiation Prediction

# In the second example, we will train model to predict the solar radiation. The orignal dataset is available in [Kaggle Solar Radiation Prediction](https://www.kaggle.com/dronio/SolarEnergy). The dataset contains such columns as: "wind direction", "wind speed", "humidity" and "temperature". The response parameter that is to be predicted is: "Solar_radiation". It contains measurements for the past 4 months and you have to predict the level of solar radiation. Let's download and load the dataset.

# In[12]:


get_ipython().system('wget https://deep-earth.s3.amazonaws.com/datasets/earthcube2021_demo/SolarPrediction.csv.zip -O SolarPrediction.csv.zip')


# In[13]:


import pandas as pd
df = pd.read_csv('SolarPrediction.csv.zip')
train_df, dev_df = train_test_split(df, random_state=100)


# In[14]:


train_df.head(10)


# Like in our previos example, we can directly train a predictor with a single `.fit()` call. The difference is that AutoGluon can automatically determine that it is a regression problem.

# In[15]:


predictor = TabularPredictor(label='Radiation', eval_metric='r2', path='ag_ec2021_demo2').fit(train_df)


# We can evaluate on the development set by calling `.evaluate()`. Here, we have specified the model to use [R2 score](https://en.wikipedia.org/wiki/Coefficient_of_determination) so it will report the R2.

# In[16]:


predictor.evaluate(dev_df)


# Similarly, we can also measure the feature importance.

# In[17]:


importance = predictor.feature_importance(dev_df)
importance


# ## More Information

# You may check our website for more information and tutorials: https://auto.gluon.ai/. We also support automatically train models with text, image, and multimodal tabular data.

# ## References

# 1. Erickson, Nick and Mueller, Jonas and Shirkov, Alexander and Zhang, Hang and Larroy, Pedro and Li, Mu and Smola, Alexander, AutoGluon-Tabular: Robust and Accurate AutoML for Structured Data, 2020, https://arxiv.org/pdf/2003.06505.pdf
