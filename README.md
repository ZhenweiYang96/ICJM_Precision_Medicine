# Title: Personalized Biopsy Schedules Using an Interval-censored Cause-specific Joint Model

**Author**: Zhenwei Yang 

**Affiliation**: Department of Biostatistics, Erasmus Medical Center

**DOI**: 10.48550/ARXIV.2209.00105

****

## Content

* ICJM Precision Medicine.Rproj
* [Rscript](#Rscript)
* [Output](#Output)
* [img](#Visualized_example)

****

### Rscript

* Data analysis

|File name| Description|
|:----------:|:--------------|
|ICJM1_inits| initials of the coefficients in the longitudinal part of ICJM1 (only include PSA)|
|ICJM1.R | fit the ICJM1 on the Canary PASS data|
|ICJM2_inits | initials of the coefficients in the longitudinal part of ICJM2 (PSA + core ratio)|
|ICJM2_inits_prepare.R | fit the longitudinal part of ICJM2 on the PASS data and extract the coefficients as the inits|
|ICJM2.R | fit the ICJM2 on the Canary PASS data|

* Data preprocessing

|File name| Description|
|:----------:|:--------------|
|combine data.R | pool all the raw information into one dataset|
|dataframe for core ratio.R | clean up the info of the core ratio|

* Functions

|File name| Description|
|:----------:|:--------------|
|[notrt]dynpred_Risk prediction function.R| the function used in scheduling to predict patient-specific risk of cancer progression based on the dynamic PSA, refer to Equation (6)
|dynpred_Risk prediction function.R | the function used in evaluation to predict patient-specific risk of cancer progression based on the dynamic PSA right after biopsies being done, refer to Equation (8)|
|Final personalized schedule evaluation function.R | a function to evaluate the personalized scheduling in the simulation study|
|function.R | help and auxiliary functions|
|ICCSJoint model function JAGS.R | a function to fit the ICJM|
|screening schedule planning function.R | a function to plan personalized schedules in the future at a stage|
|test screening schedule.R | a function to test schedules with a certain assigned risk threshold|
|true risk derivation.R | a function to calculate the true risk of cancer progression in simulation with true coefficients |

* Simulation

|File name| Description|
|:----------:|:--------------|
|seed | including the seeds to generate the data & fit the model, the seeds for individual's personalized scheduling, and model test|
|Dataset generation.R | training and test sets generation based on ICJM1|
|Estimation model.R | model fit based on 200 simualted training sets|
|Model test.R | evaluate model performance of 200 simulated datasets|
|Schedules (non progressing).R | generate personalized scheduling for patients who were not detected with cancer progression in the 200 testsets|
|Schedules (progressing).R | generate personalized scheduling for patients who were detected with cancer progression in the 200 testsets|

* Tables&figures

|File name| Description|
|:----------:|:--------------|
|Fig1.ppt | for Figure 1 in the manuscript|
|Main text_figures.R | code to generate figures in the manuscript|
|Main text_tables.R | code for the tables in the manuscript|
|Supplementary materials_figures.R | code to generate figures in the supplementary materials|
|Supplementary materials_tables.R | code for the tables in the supplementary materials|


### Output

* Data analysis

|File name| Description|
|:----------:|:--------------|
|ICJM1.R | model results (including the MCMC samples) of the ICJM1|
|ICJM2.R | model results (including the MCMC samples) of the ICJM2|

* Tables&figures


|File name| Description|
|:----------:|:--------------|
|Main text | all figures used in the manuscripts|
|Supplementary| all figures used in the supplementary materials|

### Package dependencies

- rajgs
- mcmcplots
- GLMMadaptive
- ggplot2
- tidyverse
- splines
- future
- mcmcse
- mvtnorm
- JMbayes2
- JMbayes
- MASS
- doParallel
- truncnorm
- Matrix
- latex2exp
- cowplot

> [!Note]
> - Please make sure to install above-mentioned packages by `install.packages()` before running the R code

### Visualized example

- Previous biopsy: year 0.5; current time: year 1.5:

![](https://github.com/ZhenweiYang96/ICJM_Precision_Medicine/blob/main/img/Schedule_1.5.png)

- Previous biopsy: year 0.5; current time: year 2:

![](https://github.com/ZhenweiYang96/ICJM_Precision_Medicine/blob/main/img/Schedule_2.png)

- Previous biopsy: year 0.5; current time: year 2.5:

![](https://github.com/ZhenweiYang96/ICJM_Precision_Medicine/blob/main/img/Schedule_2.5.png)

- Previous biopsy: year 2.5; current time: year 3:

![](https://github.com/ZhenweiYang96/ICJM_Precision_Medicine/blob/main/img/Schedule_3.png)