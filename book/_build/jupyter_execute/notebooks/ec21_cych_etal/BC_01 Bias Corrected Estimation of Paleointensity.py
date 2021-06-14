#!/usr/bin/env python
# coding: utf-8

# # BC_01 Bias Corrected Estimation of Paleointensity

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#BC_01-Bias-Corrected-Estimation-of-Paleointensity" data-toc-modified-id="BC_01-Bias-Corrected-Estimation-of-Paleointensity-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>BC_01 Bias Corrected Estimation of Paleointensity</a></span><ul class="toc-item"><li><span><a href="#Author(s)" data-toc-modified-id="Author(s)-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>Author(s)</a></span></li><li><span><a href="#Purpose" data-toc-modified-id="Purpose-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Purpose</a></span></li><li><span><a href="#Technical-contributions" data-toc-modified-id="Technical-contributions-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>Technical contributions</a></span></li><li><span><a href="#Methodology" data-toc-modified-id="Methodology-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>Methodology</a></span></li><li><span><a href="#Results" data-toc-modified-id="Results-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>Results</a></span></li><li><span><a href="#Funding" data-toc-modified-id="Funding-1.6"><span class="toc-item-num">1.6&nbsp;&nbsp;</span>Funding</a></span></li><li><span><a href="#Keywords" data-toc-modified-id="Keywords-1.7"><span class="toc-item-num">1.7&nbsp;&nbsp;</span>Keywords</a></span></li><li><span><a href="#Citation" data-toc-modified-id="Citation-1.8"><span class="toc-item-num">1.8&nbsp;&nbsp;</span>Citation</a></span></li><li><span><a href="#Acknowledgements" data-toc-modified-id="Acknowledgements-1.9"><span class="toc-item-num">1.9&nbsp;&nbsp;</span>Acknowledgements</a></span></li></ul></li><li><span><a href="#Setup" data-toc-modified-id="Setup-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Setup</a></span><ul class="toc-item"><li><span><a href="#Library-import" data-toc-modified-id="Library-import-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>Library import</a></span></li><li><span><a href="#Local-library-import" data-toc-modified-id="Local-library-import-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>Local library import</a></span></li></ul></li><li><span><a href="#Parameter-definitions" data-toc-modified-id="Parameter-definitions-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Parameter definitions</a></span></li><li><span><a href="#Data-import" data-toc-modified-id="Data-import-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Data import</a></span></li><li><span><a href="#Data-processing-and-analysis" data-toc-modified-id="Data-processing-and-analysis-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>Data processing and analysis</a></span><ul class="toc-item"><li><ul class="toc-item"><li><span><a href="#Looking-at-site-and-specimen-data" data-toc-modified-id="Looking-at-site-and-specimen-data-5.0.1"><span class="toc-item-num">5.0.1&nbsp;&nbsp;</span>Looking at site and specimen data</a></span></li><li><span><a href="#The-BiCEP-method--accounting-for-bias" data-toc-modified-id="The-BiCEP-method--accounting-for-bias-5.0.2"><span class="toc-item-num">5.0.2&nbsp;&nbsp;</span>The BiCEP method- accounting for bias</a></span></li><li><span><a href="#Data-manipulation-and-choosing-interpretations." data-toc-modified-id="Data-manipulation-and-choosing-interpretations.-5.0.3"><span class="toc-item-num">5.0.3&nbsp;&nbsp;</span>Data manipulation and choosing interpretations.</a></span></li><li><span><a href="#Saving-data-to-MagIC-format." data-toc-modified-id="Saving-data-to-MagIC-format.-5.0.4"><span class="toc-item-num">5.0.4&nbsp;&nbsp;</span>Saving data to MagIC format.</a></span></li><li><span><a href="#GUI" data-toc-modified-id="GUI-5.0.5"><span class="toc-item-num">5.0.5&nbsp;&nbsp;</span>GUI</a></span></li><li><span><a href="#GUI-Documentation" data-toc-modified-id="GUI-Documentation-5.0.6"><span class="toc-item-num">5.0.6&nbsp;&nbsp;</span>GUI Documentation</a></span></li><li><span><a href="#Data-References-and-Expected-Field-Values" data-toc-modified-id="Data-References-and-Expected-Field-Values-5.0.7"><span class="toc-item-num">5.0.7&nbsp;&nbsp;</span>Data References and Expected Field Values</a></span></li></ul></li></ul></li><li><span><a href="#References" data-toc-modified-id="References-6"><span class="toc-item-num">6&nbsp;&nbsp;</span>References</a></span></li></ul></div>

# ## Author(s)

# In[1]:


Author1={"name":"Brendan Cych","affiliation":"Scripps Institution of Oceanography, UCSD","email":"bcych@ucsd.edu","orcid":"0000-0003-2387-3544"}
Author2={"name":"Matthias Morzfeld","affiliation":"Scripps Institution of Oceanography, UCSD","email":"matti@ucsd.edu","orcid":"0000-0003-2257-8930"}
Author3={"name":"Lisa Tauxe","affiliation":"Scripps Institution of Oceanography, UCSD","email":"ltauxe@ucsd.edu","orcid":"0000-0002-4837-8200"}


# ## Purpose
# Paleointensities are estimates of the Earth's ancient magnetic field strength, obtained from materials including rocks and pottery sherds which cooled in that field. Paleointensity can be obtained by comparing the magnetization lost on cooling a specimen in zero field, to the magnetization gained on cooling in a known lab field. The theoretical relationship between these two magnetizations should be proportional at all temperatures for an ideal specimen, with the ratio of magnetization lost/magnetization gained being the same as the ratio of ancient field/lab field.
# 
# In natural samples, the assumptions of the paleointensity experiment are frequently violated, leading to non linear behaviour on the Arai plot (Nagata et al, 1963). Traditionally, paleomagnetists use "selection criteria" to exclude specimens from the analysis if they believe they violate the experimental assumptions. These criteria are generally based on an empirical description of non-ideal behaviour, with a specimen being excluded if it fails to meet a threshold value. There is little agreement between labs on which set of criteria to use.
# 
# Bias Corrected Estimation of Paleointensity (BiCEP) is a new method for analyzing Thellier-type paleointensity data. Instead of excluding specimens from the analysis, BiCEP assumes that the paleointensity estimate for each specimen has some amount of bias based on the values of the curvature criterion of Paterson (2011). In this way we can estimate the unbiased value of the paleointensity. 
# 
# BiCEP differs from other methods of paleointensity analysis, because it reduces the number of experiments that need to be excluded from the analysis. The credible intervals produced by BiCEP are more reliable than the standard deviation traditional selection criteria, with the precision on the estimate starting low, and increasing as more specimens are measured. It also allows for weighting of results from different specimens, as better results contribute more to the estimate than worse results. BiCEP can be used to analyze new paleointensity data, and to reanalyze old paleointensity data using a consistent methodology which significantly reduces the number of subjective choices needed to make a paleointensity estimate.
# 
# 
# ## Technical contributions
# This notebook
# - Uses the BiCEP_functions library to read in, plot and analyze paleointensity data.
# 
# - Uses a Markov Chain Monte Carlo (MCMC) sampler using the `pystan` package (http://mc-stan.org) to determine the relationship between the curvature ($\vec{k}$) parameter (Paterson, 2011) for specimens which were magnetized in the same field.
# 
# - Allows you to save these data into the MagIC database format (https://www2.earthref.org/MagIC)
# 
# - Demonstrates a GUI which allows you to perform the previous steps more easily. 
# 
# 
# ## Methodology
# - Thellier Type Paleointensity data must be imported using the MagIC format. A measurements.txt, specimens.txt, samples.txt and sites.txt are required. These are then saved as an internal csv format for use with BiCEP.
# 
# - Paleointensity data can be plotted and some limited analysis using different selection criteria can be performed.
# 
# - The BiCEP method can be fit to a collection of specimens, either sample or site, and the results can be plotted. Reanalysis is possible. For more information on the paleointensity method- see the citation for more detailed information on how the method operates.
# 
# - The results of the BiCEP method can be saved to MagIC format. 
# 
# ## Results
# The BiCEP method yields a similar number of accurate paleointensity estimates to the looser sets of traditional selection criteria analyzed, but with higher precision and accuracy more akin to tighter sets of criteria which yield fewer estimates. See Figure 7 in Cych et al (in press) for more information.
# 
# ## Funding

# In[2]:


Award1 = {"agency": "National Science Foundation", "award_code": "EAR1547263", "award_URL": "award_URL"}
Award2 = {"agency": "National Science Foundation", "award_code": "EAR1827263", "award_URL": "award_URL"}


# ## Keywords
# Include up to 5 keywords, using the template below.
# 
# keywords=["Data Analysis", "Paleomagnetism", "Paleointensity", "Markov Chain Monte Carlo", "Bayesian Statistics"]
# 
# ## Citation
# Cych et al (In Press) doi: 10.1002/essoar.10506403.1 https://www.essoar.org/doi/10.1002/essoar.10506403.1
# 
# ## Acknowledgements 
# We are deeply grateful to Lennart de Groot and Greig Paterson for their very helpful reviews and for the advice and guidance given by Andrew Roberts, David Heslop and Joseph Wilson.

# # Setup
# 
# ## Library import

# In[3]:


import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
from IPython.display import Image


# ## Local library import

# In[4]:


import sys
sys.path.append('modules_and_scripts/') 
#The BiCEP_functions library imports all other needed libraries
import BiCEP_functions as BiCEP


# # Parameter definitions

# In[5]:


B_exp=36.4 #The expected field for the site we will be using BiCEP with today


# 
# # Data import
# 
# We can produce a data file that BiCEP can read internally from the magic format using the `generate_arai_plot_table` function in BiCEP_GUI. We will specify an output name of `example_data.csv`. For your convenience, we have included a set of MagIC formatted data to convert to the BiCEP format. This may take a few minutes

# In[6]:


BiCEP.generate_arai_plot_table('example_data')


# When imported, BiCEP data is stored internally as a set of nested objects, one for the dataset which contains a object for each site (a set of samples expected to have been magnetized in the same field), or each sample (an object collected in the field) which are described using a SpecimenCollection object. Each of these SpecimenCollection objects contains an object for each specimen, which is a subsample that a set of independent measurements are made on. To make the data class, we can initialize a `ThellierData` object using our `example_data.csv` file. This may take a few seconds to create all of the objects. Printing our ThellierData object allows us to see the number of sites and specimens.

# In[9]:


data=BiCEP.ThellierData('example_data.csv')
print(data)


# # Data processing and analysis

# ### Looking at site and specimen data
# We can index for our sites and specimens using a key for the site name. Let's look at site "hw126" and print out the object to see the specimens. This site is a lava flow erupted in Hawaii in 1935 and originally analyzed by Cromwell et al (2015). Because this eruption occurred in the age of instrumentation, we know that the original field this lava flow was erupted in was approximately 36.4 μT.

# In[10]:


site=data['hw126']
print(site)


# Let's take a look at some of these data from one of the specimens from this site. We can use the `plot_arai()` method on our specimen to plot the Thellier data on an Arai plot (Nagata & Arai, 1963), which shows the magnitude of the magnetization lost by cooling in zero field against the magnetization gained when cooling in a known field. We can also plot a Zijderveld plot, which plots the magnetization vectors of the zero field steps. In the zijderveld plot, black circles represent the component of the magnetization in the x,y plane, whereas red squares represent the component in the x,z plane.

# In[19]:


specimen=site['hw126a1']
fig,ax=plt.subplots(1,2,figsize=(12,4))
specimen.plot_arai(ax[0])
specimen.plot_zijd(ax[1])
ax[0].set_title('Arai Plot')
ax[1].set_title('Zijderveld Plot');


# Both our Arai plot and Zijderveld plot are a straight line. Because we're plotting magnetization lost in zero field, against magnetization gained in a known field, the ratio of our two lines is the ratio of the Ancient field to the Lab field which is a little under 2. We can check the lab field using the `B_lab` attribute of the specimen.

# In[20]:


specimen.B_lab


# The lab field was 20 μT. If we expect a value of 36.4 μT, then a slope of slightly under 2 is what we might expect. Let's look at another specimen.

# In[29]:


specimen=site['hw126a2']
fig,ax=plt.subplots(1,2,figsize=(12,4))
specimen.plot_arai(ax[0])
specimen.plot_zijd(ax[1])
ax[0].set_title('Arai Plot')
ax[1].set_title('Zijderveld Plot');


# The data for this specimen do not follow a straight line. Although the line slope appears to start out at a reasonable value, it quickly becomes shallower at higher temperatures. If we were to fit a line to these data we might end up with an underestimate of the paleointensity, and there is no reason to suggest which temperature range to pick for our intensity data. Some strict selection criteria might exclude this specimen, others might include it, but could include an incorrect temperature range, which could give a biased paleointensity estimate.

# ### The BiCEP method- accounting for bias
# 
# BiCEP works differently. It fits a circle to a scaled set of paleointensity data and uses the tangent joining the circle center and the origin. BiCEP generates many circle fits from a probability distribution. Let's look at the circle fits that BiCEP samples from the posterior distribution, and the median tangent which we use to fit the line. BiCEP fits are performed at a site level for a reason which will become apparent later. We can use the method `site.BiCEP_fit()` to fit the BiCEP method (this may take several minutes depending on the speed of your processor and memory available) and plot the circle and tangent fits to the data.

# In[30]:


#There is also a model_circle_slow which fits these circles using a slower model,
#use this if you are having problems with sampling.
site.BiCEP_fit(model=BiCEP.model_circle_fast)
#This cell may generate warnings from the internal pystan code. 
#The main error to watch for is one about R_hat, which indicates the sampler hasn't converged, a "bad" run


# We can use the `specimen.plot_circ()` method to draw a circle under the Arai plot data

# In[75]:


fig,ax=plt.subplots(1,2,figsize=(12,4))
specimen.plot_circ(ax[0],tangent=True)
specimen.plot_arai(ax[0])
specimen.plot_zijd(ax[1])
ax[0].set_title('Arai Plot')
ax[1].set_title('Zijderveld Plot');


# We notice that this tangent is slightly less than the initial slope, how does this help us? 
# 
# BiCEP assumes that the curvature, $\vec{k}$ (1/radius with sign depending on direction of curvature, Paterson, 2011) of the circle fit is proportional to the calculated ancient field ($B_{anc}$). For $\vec{k}=0$, our "circle" has infinite radius, and so becomes a perfect line. BiCEP tries to fit a linear relationship between $B_{anc}$ and $\vec{k}$ for all specimens in a site. The $B_{anc}$ value can be extrapolated (or interpolated) to $\vec{k}$=0 to give us our unbiased paleointensity estimate. In the case of the plotted specimen hw126a2, BiCEP predicts that this specimen **should** underestimate the site paleointensity slightly, due to its curvature.
# 
# We can plot the relationship between $B_{anc}$ and $\vec{k}$ using the method `site.regplot()`. We plot the expected field value as a red line.

# In[76]:


fig,ax=plt.subplots()
site.regplot(ax)

#x limits, y limits and y label are not set by this method as it is used slightly differently in our GUI.

ax.set_xlim(-0.2,1.2)
ax.set_ylim(15,50);
ax.set_ylabel('$B_{anc}$ ($\mu$T)')
ax.axhline(36.4,color='r',zorder=0);


# The green dots with black errorbars are our individual specimen estimates for $\vec{k}$ and $B_{anc}$ and their 95% confidence intervals. The blue lines are our BiCEP fits to these data. Notice that these lines cross the $\vec{k}=0$ line close to the expected field value, yielding an accurate and precise estimate of the paleointensity, without excluding any specimens! We can get a sense of the distribution of estimates for $B_{anc}$ using the `site.histplot()` method.

# In[19]:


fig,ax=plt.subplots()
site.histplot(ax)

ax.set_xlim(30,40)
ax.axvline(36.4,color='r');


# We can see that our estimate for $B_{anc}$ is close to the median value of this distribution, and well within the 95% confidence interval denoted by the black line. Our estimate is accurate and precise, with a 95% confidence interval of ~4 μT

# ### Data manipulation and choosing interpretations.
# 
# You might notice that some of our specimens have larger error bars on their $\vec{k}$ estimates than others in our BiCEP fitting figure. Let's investigate what's going on with these specimens.

# In[77]:


fig,ax=plt.subplots(1,2,figsize=(12,4));
hw126a7=site['hw126a7']
hw126a8=site['hw126a8']
hw126a7.plot_circ(ax[0],tangent=True)
hw126a7.plot_arai(ax[0])
ax[0].set_title('hw126a7')

hw126a8.plot_circ(ax[1],tangent=True)
hw126a8.plot_arai(ax[1])
ax[1].set_title('hw126a8');


# For these specimens, our analysis might be slightly different. Our specimen "hw126a7" has a pTRM check (triangle) which is very displaced from the line at the 300 °C temperature step. This is a repeat of the 300 °C in-field measurement after heating to 350 °C in zero field. It indicates a chemical alteration of the specimen between 300 °C and 350 °C. We can use this by calculating the DRAT statistic (Selkin & Tauxe, 2000), using the attribute `specimen.drat`, and noticing that it is high (criteria vary, but 15 is considered a high value by any metric). Although we have already demonstrated the effectiveness of BiCEP without using this statistic, we have evidence that we should only be using the temperatures up to 300 °C for this specimen. We can do this using the method `specimen.change_temps()`, and we save this change using the method `specimen.save_changes()`. Note that our attributes, such as `specimen.drat` are affected by the change in temperatures.
# 
# Our specimen "hw126a8" has a very early failed pTRM check, and fitting a circle seems like it would probably be inappropriate. We can elect to exclude this specimen as no part of the line is likely to give a good paleointensity estimate. This can be done by setting the attribute `specimen.active=False` Ultimately, due to its large uncertainty in the circle fit, it has very little effect on the estimate of $B_{anc}$ in the first place.
# 
# The MAD statistic of Kirshvink (1980) and the DANG Statistic (Tanaka & Kobayashi, 2003, Tauxe & Staudigel, 2004) are also available using the methods `specimen.mad` and `specimen.dang`.

# In[79]:


#Printing the DRAT parameter of the specimen to 
print('hw126a7 DRAT: %2.1f'%hw126a7.drat)
print('hw126a8 DRAT: %2.1f'%hw126a8.drat)

lowerTemp=0
upperTemp=300
lowerTempK,upperTempK=lowerTemp+273,upperTemp+273

hw126a7.change_temps(lowerTempK,upperTempK)
hw126a7.save_changes()
hw126a8.active=False


# We now run the BiCEP method again with these changes in place, and plot up our site results with the results for these two specimens.

# In[80]:


site.BiCEP_fit(model=BiCEP.model_circle_fast)


# In[83]:


#set up subplots
fig,ax=plt.subplots(2,2,figsize=(12,8))

#Plot our Arai plots with circle fits
hw126a7.plot_circ(ax[0,0])
hw126a7.plot_arai(ax[0,0])

hw126a8.plot_circ(ax[0,1])
hw126a8.plot_arai(ax[0,1])


#Plot our plot of the BiCEP fit.
site.regplot(ax[1,0])
ax[1,0].set_xlim(-0.2,1.2)
ax[1,0].set_ylim(15,50);
ax[1,0].set_ylabel('$B_{anc}$ ($\mu$T)')
ax[1,0].axhline(36.4,color='r',zorder=0);


#Plot our histogram
site.histplot(ax[1,1])
ax[1,1].set_xlim(30,40)
ax[1,1].axvline(36.4,color='r');


# We notice from our histogram that doing this data manipulation has very slightly improved our accuracy and precision- but the overall distribution has not changed vary much. We notice from our Arai plot for specimen hw126a7 that the uncertainty in our circle fits has become much larger by excluding these temperature steps. This is because we assume the uncertainty in the circle fit for the whole Arai plot, rather than just the segment we're looking at. In this way, we are penalized for excluding data from the Arai plot, without having to use some kind of "length of the line" criterion. 

# ### Saving data to MagIC format.
# 
# Once we have run our site- we can save our results directly into the MagIC formatted data file. This will save results of the BiCEP method for our site and all specimens in that site. We can do this using the method
# `site.save_magic_tables()`

# In[84]:


site.save_magic_tables()


# ### GUI
# As we think that using these functions for every site will be cumbersome, we present a GUI which performs all of these steps. The GUI can save figures produced, and also produces diagnostics about the sampler, as well as important diagnostic information about the sampler. We can use the `run_gui()` function to initialize the GUI. Note that the `%matplotlib widget` magic command is required to have the BiCEP GUI display plots correctly. Above figures can be set to the normal plotting style when rerunning cells using the `%matplotlib inline` magic command. 

# In[7]:


Image('readme_image/GUI_layout.png')


# ### GUI Documentation

# On launch you should have a GUI with the above layout:
# 
# 1. File selection button. Press select, choose your file, and press select again. Then press "Run" to import the data to the GUI. You cannot then select a new file.
# 
# 2. Site and specimen dropdowns. These dropdown menus allow you choose a particular paleointensity experiment.
# 
# 3. Minimum and maximum temperature steps (in Celcius) to use for the paleointensity experiment. We recommend using the Zijderveld plot (7.) and pTRM checks to choose which set of temperatures to use. By default, we use all temperature steps to make a paleointensity estimate. Currently it is required to make an estimate for all specimens.
# 
# 4. Statistics about the direction and alteration of the ChRM used for paleointensity estimation. These may help with choosing which set of temperature steps to use. See the standard paleointensity definitons (Paterson et al, 2014, https://earthref.org/PmagPy/SPD/DL/SPD_v1.1.pdf) for more information on these statistics. In addition to these statistics, we present the worst R_hat diagnostic for a specimen. If R_hat>1.1 or R_hat<0.9, it may indicate an issue with the sampler (see 13.). In this case, this box will show up as red, and the specimen may be excluded using the checkbox (8.))
# 
# 4. Arai plot with zero field first steps plotted as red circles, in field first steps plotted as blue circles, pTRM checks plotted as triangles, and pTRM tail checks plotted as squares. Additivity checks are not currently plotted. Circle fits from the BiCEP method will be plotted as green lines under the Arai plot after the site fit (9) has been performed. All plots can be rescaled using the "move" button (3rd symbol from the bottom on left side of plot) and right clicking and dragging, or the "zoom" button (2nd symbol from the bottom) and left clicking and dragging to zoom in, or right clicking and dragging to zoom out. The "home" button (second symbol from the top) resets the plot axis, as does changing the temperatures.
# 
# 5. Zijderveld plot of the data, with x,y plotted as black circles and x,z plotted as red squares.
# 
# 7. Saves the temperatures used for that specimen to RAM. This must be done for each specimen individually to change temperatures before running the sampler (9.). By default, all temperature steps are used for every specimen.
# 
# 8. Checkbox for excluding a specimen. This should only be done if there is no reasonable interpretation of the specimen (e.g. alteration at low temperature, not demagnetizing to the origin). 
# 
# 9. The "Process Site Data" button performs the BiCEP method on that site and calculates the site level paleointensity. Depending on the speed of your computer and the sampler parameters used (10), this may take a while to run, especially for sites with many specimens. Please be patient.
# 
# 10. Parameters for the MCMC sampler for the BiCEP method. The "n samples" slider increases the number of samples used in the MCMC sampler. Smaller numbers will take less time to run but result in lower accuracy in the resulting probability distribution. The "Sampler" button changes the parameterization of the MCMC sampler slightly (mathematically, the model is the same, but the parameters being sampled from are specified slightly differently). The "Slow, more accurate" sampler is much slower than the "Fast, less accurate" sampler, but generally (though not always) results in better sampler diagnostics than the "Fast, less accurate" sampler, particularly for sites with small numbers of specimens.
# 
# 11. Plot of the estimated paleointensity for each specimen against Arai plot curvature. The currently displayed specimen in the Arai and Zijderveld plots has a red circle around it in this plot. The blue lines are samples from the posterior distribution for the relationship between specimen level paleointensity and curvature. The y intercept is the estimated site level paleointensity.
# 
# 12. Histogram of the site level paleointensities sampled from the posterior distribution. This corresponds to the distribution of intercepts of the blue lines in (12.).
# 
# 13. Diagnostics for the MCMC sampler (see Cych et al, in prep. or the Stan Documentation at https://mc-stan.org/docs/2_26/reference-manual/notation-for-samples-chains-and-draws.html, https://mc-stan.org/docs/2_26/reference-manual/effective-sample-size-section.html). 0.9<R_hat<1.1 and n_eff>1000 is desired, with R_hat=1.00 and n_eff>10000 being ideal. Tweak the sampler parameters (10.) or measure more specimens if these parameters give poor results (indicated by an amber color for n_eff<1000 or a red color for bad R_hat). Also displayed here is the 95% credible interval for the site and the Category (see Cych et al once again for an explanation). The color of the category box indicates how to proceed. Green (Category A or B): Accept site, Amber (Category C or D): Measure more specimens, Red (Category D): Ignore site.
# 
# 14. Saves figures to file. Currently the Zijderveld plot and Arai plot have to be saved together (as do both site plots).
# 
# 15. Saves the results from the BiCEP method to the MagIC tables (site and specimen tables are altered).

# In[85]:


get_ipython().run_line_magic('matplotlib', 'widget')
BiCEP.run_gui();


# ### Data References and Expected Field Values
# We follow with a list of DOIs for sites included in our data file as well as the expected field values at those sites. This can be useful for comparison when using BiCEP_GUI.

# In[23]:


site_info=pd.read_csv('site_info.txt')
site_info


# # References
# List relevant references:
# 
# Cromwell, G., Tauxe, L., Staudigel, H., & Ron, H.(2015). Paleointensity estimates from historic and modern hawaiian lava flows using glassy basalt as a primary source material. Phys. Earth Planet. Int., 241, 44–56.  doi:89710.1016/j.pepi.2014.12.007
# 
# Kirschvink, J. L. (1980), The least-squares line and plane and the analysis of palaeomagnetic data,Geophys. J. R. Astr. Soc.,62, 699–718, doi:10.1111/j.1365-246X.1980.tb02601.x
# 
# Nagata, T., Arai, Y., & Momose, K.  (1963).  Secular variation of the geomagnetic total force during the last 5000 years. J. Geophys. Res.,68(18), 5277-5281.   doi:0.1029/j.2156-2202.1963.tb00005.x
# 
# Paterson, G. A. (2011), A simple test for the presence of multidomain behaviour during paleointensity experiments,J. Geophys. Res.,116, B10,104, doi:10.1029/2011JB008369.
# 
# Paterson,  G.  A.,  L.  Tauxe,  A.  J.  Biggin,  R.  Shaar,  and  L.  C.  Jonestrask  (2014),  On  im-proving  the  selection  of  thellier-type  paleointensity  data,Geochem.  Geophys.  Geosyst.,  doi:10.1002/2013GC005135.
# 
# Selkin, P. A., and L. Tauxe (2000), Long-term variations in palaeointensity, Phil.  Trans.  R.  Soc.London,358, 1065–1088, doi:10.1098/rsta.2000.0574.
# 
# Tanaka,  H.,  and  T.  Kobayashi  (2003),  Paleomagnetism  of  the  late  Quaternary  Ontake  Volcano,Japan:  directions, intensities, and excursions,Earth Planets Space,55, 189–202.
# 
# Tauxe, L., and H. Staudigel (2004), Strength of the geomagnetic field in the Cretaceous Normal Superchron:  New data from submarine basaltic glass of the Troodos Ophiolite,Geochem. Geophys.Geosyst.,5, Q02H06, doi:10.1029/2003GC000635
# 
# Zijderveld, J. D. A.(1967).A. c. demagnetization of rocks:  Analysis of results. In D. W. Collinson, K. M. Creer, & S. K. Runcorn (Eds.), Developments in solid earth geophysics (Vol. 3, pp. 254–286). Elsevier. doi:102110.1016/B978-1-4832-2894-5.50049-5
