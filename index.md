**UNDER CONSTRUCTION**

# Singapore’s SARS-CoV-2 Infection Trend In 2020: A Lower Bound Estimate

By Julian Chia, 1<sup>st</sup> March 2022

_In remembrance of the 2<sup>nd</sup> Anniversary of Singapore's COVID-19 epidemic._

## Abstract
Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2) infection causes Coronavirus Disease 2019 (COVID-19). Knowledge on the SARS-CoV-2 infection trend is, however, lacking. This viral infection is invisible to the naked eye and is challenging to profile in real-time. Its closest indication is its documented COVID-19 epidemic trend. Fortunately, that is published daily and globally due to the pandemic situation of COVID-19. Leveraging the empirical nature and availability of these COVID-19 epidemic trends, this paper posits that these trends are, in fact, lower-bound estimates of the respective localities COVID-19 epidemic situation, and each is interrelated to the lower-bound SARS-CoV-2 infection trend in their locality. A model based on this posit is developed and applied to Singapore. The predicted Local SARS-CoV-2 infection trends provided a novel reference to understanding the Singapore COVID-19 epidemic that was previously not possible. They evidenced the window of opportunity where Singapore could have mitigated its COVID-19 epidemic via its Circuit Breaker (CB) that it had missed. They evidenced Singapore's additional tighter CB measures and extended CB dateline, both implemented during the CB on 21<sup>st</sup> April 2020, were timely and effective. They derived the population of imminent COVID-19 individuals that the empirical Local COVID-19 epidemic trend had undocumented daily. These undocumented populations are sizable and a possible factor for the COVID-19 epidemic and its protracted recovery. Finally, these SARS-CoV-2 trends provided circumstantial evidence that Singapore's COVID-19 epidemic originated from COVID-19 cases imported into Singapore.

## Keywords
SARS-CoV-2 infection trend· COVID-19 epidemic trend· COVID-19 Confirmation Period · Normal Distribution · Cartesian Product · Backcasting · Forecasting · Lower-bound · Modelling · Python Programming  · 

## 1 Introduction

Singapore reported its first two Coronavirus Disease 2019 (COVID-19) [1] cases on 23<sup>rd</sup> January 2020 [2]. This disease onset evolved from several independent sporadic outbreaks into an epidemic within that year. Publicly accessible daily COVID-19 Case Reports and Situation Reports issued by the Ministry of Health of Singapore [3] describe the extent of this epidemic. The _number of daily confirmed COVID-19 cases_ is their unit of measurement of daily COVID-19 prevalence; they are the daily cumulative positive real-time reverse transcription-polymerase chain reaction (rt-PCR) tests and serologic assays [4]. Its daily compilation constructs the Singapore COVID-19 epidemic curves [5] for Imported and Local cases.

![Figure 1](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_1_SG_COVID19_%20Epidemic_trends.png?raw=true)
**Figure 1:** Singapore’s COVID-19 Epidemic Curves from 23<sup>rd</sup> January to 18<sup>th</sup> August 2020.

Figure 1 illustrates these Imported and Local COVID-19 epidemic trends of Singapore (from 23<sup>rd</sup> January to 18<sup>th</sup> August 2020). It evidences the pervasiveness of the COVID-19 epidemic in Singapore that caused the implementation of a nationwide Circuit Breaker (CB), alternatively known as a Lockdown, that year. Although these trend curves are informative, they do not describe the events that caused COVID-19, i.e. they do not show when Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2) infection had occurred. The SARS-CoV-2 infection trend curve fulfils this function, yet it is unreported. A means to monitor SARS-CoV-2 transmission in real-time is unavailable. Also, COVID-19 patients could not see and did not know when they had contracted SARS-CoV-2.
  
There is unanimous consensus that the SARS-CoV-2 infection trend is not identical to the COVID-19 epidemic trend. One reason for this is that SARS-CoV-2 infection is a precursor of COVID-19. SARS-CoV-2 infection takes time to cause COVID-19 and takes time to be discovered and confirmed [6] by the patients, their doctors and contact tracers. Furthermore, COVID-19 confirmation can occur during the incubation, symptomatic and recovery periods of a SARS-CoV-2 infection as positive rt-PCR results are obtainable during these periods regardless of the virus transmissibility [7,8]. Presymptomatic SARS-CoV-2 transmissions, which occurs during the incubation period, have been evidenced in Singapore [9], China [10], and the USA [11]. In 2020, the incubation period of SARS-CoV-2 variants averaged at 5 to 6 days and ranged from 1 to 14 days [12,13]. Thereon, positive rt-PCR have evidenced too. A multi-centre USA study of 70,406 unique COVID-19 patients [14] found that a majority of COVID-19 patients yielded positive rt-PCR, i.e. shed SARS-CoV-2 ribonucleic acid (RNA), three weeks after their first positive SARS-CoV-2 rt-PCR test. Another study that reviewed 77 COVID-19 reports [15] found that the duration of SARS-CoV-2 RNA shedding can range from a minimum of 1 day to a maximum of 83 days. According to WHO [16], rt-PCR positivity generally appears up to 3 weeks or more for mild to moderate COVID-19 and a more extended period for severe COVID-19. These findings evidence the amount of lag that a COVID-19 epidemic trend could have against a SARS-CoV-2 infection trend.
	
The second manner where the SARS-CoV-2 infection trend is dissimilar to the COVID-19 epidemic curve is in their population. The  SARS-CoV-2 infectee population, in reality, is more than the documented COVID-19 epidemic population. The under-reporting of asymptomatic and presymptomatic SARS-CoV-2 transmissions cause this phenomenon [17]. In Wanzhou, China, a 14-weeks COVID-19 mass testing program [18] found that asymptomatic and presymptomatic SARS-CoV-2 transmission accounted for 75.9% of all SARS-CoV-2; the abundance of close contacts before symptom onset or diagnosis facilitated them. Using a model that incorporates daily testing information fit to the COVID-19 case and serology data from New York City, researchers estimated that presymptomatic and asymptomatic SARS-CoV-2 transmissions together comprised at least 50% of the infections at the outbreak peak [19]. A decision analytical model study estimated that asymptomatic and presymptomatic SARS-CoV-2 accounted for >50% of all transmissions [20]. These studies exemplify the extent of under-reporting possibility occurring between the documented COVID-19 epidemic population against the actual SARS-CoV-2 infectee population.
	
Succinctly put, the actual SARS-CoV-2 infection trend will always appear ahead and exhibit a larger integral than its empirical COVID-19 epidemic trend. Its modelling is non-trivial to perform; the number of transmission factors for SARS-CoV-2 is many, their quantification remains challenging to implement, and SARS-CoV-2 transmission is opportunistic. To simplify its modelling while retaining an empirical basis, I posit that the officially published Local COVID-19 epidemic trend is, in fact, a lower-bound estimate of the Local COVID-19 epidemic and that it is a result of the lower-bound Local SARS-CoV-2 infection trend. Such a postulation is the basis of this paper.
	
The following sections present a model that utilises the abovementioned postulation, statistics, Cartesian-product, computation algorithm and empirical Local COVID-19 epidemic data to estimate the lower-bound Local SARS-CoV-2 infection trend of Singapore from January to August 2020. Its results contain the findings on the estimation of Singapore's daily mean Local COVID-19 confirmation period trend, lower-bound daily Local COVID-19 epidemic trend, and lower-bound daily Local SARS-CoV-2 infection trend. Their discussions give new insights into Singapore's COVID-19 epidemic situation in 2020.
  
## 2 The Model

### 2.1 Assumptions

The primary assumptions of the model are:
1. Local SARS-CoV-2 infections that are presymptomatic, asymptomatic or symptomatic but are undocumented are negligible. Therefore, an empirical Local COVID-19 epidemic trend is a reasonable lower-bound estimate of the Local COVID-19 epidemic trend and is relatable to the lower-bound Local SARS-CoV-2 infection trend that precedes it.
2. The Local COVID-19 confirmation event always lags behind its Local SARS-Cov-2 infection event. The duration between these two events is called the _COVID-19 confirmation period (CCP)_. It follows the Normal/Gaussian distribution theory. Daily, the _CCP_ probability density function is:

   <img src="https://render.githubusercontent.com/render/math?math={\color{black} P(X) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp^{-\frac{(X-\mu)^2}{2\sigma^2}}}"> ..........(1)

   Here, _X є (-∞,∞)_ is the random variate of Eqn(1) and denotes the daily _CCP_, _μ_ denotes its daily mean, _σ_ denotes its daily standard deviation, and _π_ denotes the pi constant. _μ_ is not constant over the Local COVID-19 population history, i.e. _μ_ is a function of days.

### 2.2 Hypothesis

The _backcasting_ of a documented Local COVID-19 epidemic curve using _μ_ and _σ_ yields a lower-bound estimate of the SARS-CoV-2 infection trend. Vice versa, the _forecasting_ of the lower-bound Local SARS-CoV-2 infection trend with _μ_ and _σ_ yields an estimate of the lower-bound Local COVID-19 epidemic trend, which is the empirical Local COVID-19 epidemic trend.

### 2.3 Methodology

_Backcasting_ is a statistical computation algorithm performed from the Local COVID-19 confirmation event day. It involves: 
1. Estimating the _CCP_ of each Local COVID-19 patient with Normal distribution theory. 
2. Determining the SARS-CoV-2 infection event day of each Local COVID-19 patient by subtracting the estimated _CCP_ from the COVID-19 confirmation event.
3. Getting a histogram of all the estimated SARS-CoV-2 infection events. This histogram estimates the lower-bound Local SARS-CoV-2 infection trend of the Local COVID-19 population.

_Forecasting_ similarly is a statistical computation algorithm. But unlike backcasting, its execution is from the SARS-CoV-2 infection event day. It involves:
1. Estimating the _CCP_ of each Local SARS-CoV-2 infection event (this step is identical to backcasting).
2. Determining the Local COVID-19 confirmation event day of each Local SARS-CoV-2 infectee by adding the estimated _CCP_ to the Local SARS-CoV-2 event.
3. Getting a histogram of all the estimated Local COVID-19 confirmation events. This histogram yields the lower-bound Local COVID-19 epidemic trend that should resemble the empirical Local COVID-19 epidemic trend.

In _backcasting_ and _forecasting_, the profile of the estimated lower-bound Local SARS-CoV-2 infection and COVID-19 epidemic trends rely on the values of the _CCP_, which in turn are a function of _μ_ and _σ_. Empirical data on the _μ_ and _σ_ of Singapore’s Local COVID-19 epidemic trend are unpublished. As such, they have to be estimated.

The procedure developed to estimate _μ_ involves:
1. Iteratively perform the backcasting and forecasting algorithms for a range of constant _μ_ scenarios (termed as _μ<sub>c</sub>_).
2. Select _μ<sub>c</sub>_ as a probable value of _μ_ for a given day when the number of estimated Local COVID-19 confirmation events that day is identical to its empirical data.
3. Complete this selection procedure upon reaching 300 samples of the _μ_ estimates per day for the Local COVID-19 epidemic population.
The approximation of the range of _μ<sub>c</sub>_ in step 1 of this procedure is determined using empirical CCP data. Figure 2 illustrates the CCP from a sample of Local COVID-19 cases obtained from [3] and various national news sources. Its cumulative probability-density distribution is Gaussian-like with a mean and standard deviation of 17.531 days and 6.044 days, respectively. Assuming this empirical mean value indicates the maximum range of _μ<sub>c</sub>_ while letting a day be its minimum, then _μ<sub>c</sub>_=[1,18] since the ceiling value of 17.531 days is 18 days. Furthermore, _σ_ is assumed to be a quadratic function of _μ_:

   <img src="https://render.githubusercontent.com/render/math?math={\color{black} \sigma = f(\mu) = a\mu^2 %2B b\mu %2B c}"> ..........(2)

   where _a_ = -0.008665, _b_ = 0.483888, _c_ = 0.0 and _μ_=_μ<sub>c</sub>_. Figure 2 illustrates the Normal Cumulative Probability Distributions of _μ<sub>c</sub>_=[1,18].

The estimation of _μ_ via the abovementioned procedure will, unfortunately, continue indefinitely when one or more elements of _μ_ are unpredictable, i.e. when the element(s) of _μ_ have zero value. In such an eventuality, a _Resemblance Algorithm_ is to complete the estimation of _μ_. Firstly, the Cartesian Product of _μ<sub>c</sub>_=[1,18] replaces the element(s) of _μ_ that has zero value. Doing so allows the exploration of every possible sequencing of _μ<sub>c</sub>_ in these missing _μ_ elements. Next, the backcasting and forecasting of each possible sequencing of these _μ_ estimates the Local COVID-19 epidemic trend. Finally, the sequence of _μ_ that yields the Local COVID-19 epidemic trend that best resembles its empirical counterpart is its _μ_. The measurement of resemblance is via a _cumulative-absolute-difference_ (_CAD_) criterion:

<img src="https://render.githubusercontent.com/render/math?math={\color{black} CAD = \displaystyle\sum_{d=0}^{d_{max}} |{T_e-T_m}|}"> ..........(3)
   
and its _weighted_ counterpart (_WCAD_) criterion:

<img src="https://render.githubusercontent.com/render/math?math={\color{black} WCAD = \displaystyle\sum_{d=0}^{d_{max}}(\frac{T_e}{p} * |{T_e-T_m}|)}"> ..........(4)
   
Here, _T<sub>e</sub>_ and _T<sub>m</sub>_, respectively, denote the estimated and documented daily number of Local COVID-19 cases, _d_ denotes the day, _d<sub>max</sub>_ denotes its maximum, and _p_ denotes the population of the empirical Local COVID-19 epidemic. A complete resemblance occurs when the criterion → 0. The opposite is true when their value → ∞.

![Figure 2](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_2_SG_CCP_empirical_%26_model.png?raw=true)
**Figure 2:** Illustration of the _CCP_ of a sample of Singapore’s Local COVID-19 cases, their cumulative probability distribution, and the Gaussian cumulative probability distribution parameters used to estimate the _μ_ of Singapore’s Local COVID-19 epidemic trend.

### 2.4 Computation

The implementation of the methodology is by the Python3 scripting language [21] and optimized libraries such as NumPy [22] and SciPy [23]. The visualization of their results is through Matplotlib [24]. Their source codes are in [25]. All computations are by a workstation installed with an overclocked Intel® Core™ i9-7960X CPU comprising 32 logical cores and 94.0GB of DDR4 RAM.

A strategy to achieve high computation efficacy is executing an instance of the _concurrent.futures.ProcessPoolExecutor_ class of Python3 within a nested logical _while_-loop structure. This arrangement provides a continuous-concurrent stream of computation using every available logical core of the CPU.  Also, large three-dimensional instances of the NumPy _ndarray_ class facilitated data parallelism within each CPU logical core. The 1<sup>st</sup>, 2<sup>nd</sup> and 3<sup>rd</sup> dimensions of these ndarrays, respectively, represent the number of iterations performed in each CPU logical-core, the range of _μ<sub>c</sub>_, and the Local COVID-19 epidemic population. On this workstation, the optimum _ndarray_ size to operate 28 logical cores is 150x18x55136= 148,867,200 elements.

The generation of pseudo-randomness in the results are by the NumPy Permuted Congruential Generator 64-bit with DXSM (_PCG64DXSM_) class and the NumPy _random.SeedSequence_ class. Also, three unique random seeds had primed the _random.SeedSequence_ class, and thus the _PCG64DXSM_ bit generator, to understand possible variance in their computed results.
	
## 3 Results

### 3.1 The Statistical _μ_ Estimates

Figure 3 shows that the statistical estimation of _μ_ for each random seed took 14 to 15hrs or ~2.48x10<sup>6</sup> iterations to complete ~84.58% of the Local COVID-19 population history of 214 days. This result meant that 181 days of the 214 days achieved the sampling quota of 300 estimates per day. Of the remaining 15.42% (or 33 days of the 214 days) that are incomplete, 6.54-7.48% (or 14-16 days of the 214 days) achieved no estimate of _μ_. Beyond these thresholds, the statistical estimation of _μ_ became unproductive and discontinued. Figure 4 shows the mean of the estimated _μ_, i.e. _μ<sub>mean</sub>_, of each unique random seed. The low variance in these _μ<sub>mean</sub>_ trends indicates that the sampling quota of the _μ_ estimates is sufficiently large.

Figure 5 shows the _μ<sub>mean</sub>_ trend of all the _μ_ estimates combined. ~85.05% of the 214 days (or 182 days) achieved the sample quota of _μ_ while ~14.95% of 214 days (or 32 days) did not. Of which, ~6.07% of the 214 days (or 13 days) had no estimates of _μ_. The completion of the _μ<sub>mean</sub>_ trend via the _Resemblance Algorithm_ is needed and presented in the next section.

![Figure 3](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_3_%CE%BC_computation.png?raw=true)
**Figure 3:** Illustration of the duration and iterations to estimate _μ_. “Completed” denotes achieving a sample quota of 300 estimates per day. “Incomplete” denotes not achieving the sample quota. “No Estimate” is a subset of “Incomplete” without an estimation of _μ_.

![Figure 4](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_4_%CE%BC_estimates.png?raw=true)
**Figure 4:** The mean of the respective estimates of _μ_ for Random Seeds 1, 2 and 3. μmean=0 denotes no estimate of _μ_.

![Figure 5](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_5_%CE%BCmean_incomplete.png?raw=true)
**Figure 5:** The mean of the combined estimates of _μ_ for Random Seeds 1, 2 and 3. _μ<sub>mean</sub>_=0 denotes no estimate of _μ_.

### 3.2 The Estimated _μ<sub>mean</sub>_ Trend

The Cartesian product of _μ_=[1,18] for thirteen missing elements of the _μ<sub>mean</sub>_ trend in Figure 5 is needed. Accordingly, 18<sup>13</sup>=2.082296487×10¹⁶ possible sequencing of these missing μmean elements, as well as the _backcasting-forecasting_ and the _CAD-WCAD_ treatments of their constituted μmean trends, needs computing. The completion of this Big-O time complexity is too computationally intensive to achieve on a workstation. Addressing this issue requires a reduction of the problem size. To this end, the following _Modified Resemblance Algorithm_ is implemented:
1. The Cartesian product of _μ_=[1,18] shall not exceed five missing _μ<sub>mean</sub>_ elements for each computation. This decision discretizes the computation problem into three reasonably sized prediction steps (since 13days//5days=3 and 18<sup>5</sup>=1,889,568 iterations of the _backcasting-forecasting_ and the _CAD-WCAD_ treatments of their constituted _μ<sub>mean</sub>_ trends takes an hour or two to complete). The selection of these five missing _μ<sub>mean</sub>_ elements is according to whether they have the five highest empirical Local COVID-19 case counts and thus reordered. This reordering is performed on the three missing _μ<sub>mean</sub>_ elements of the last prediction step too.
2. In the 1st prediction step, the values of the eight missing _μ<sub>mean</sub>_ elements that are not selected are made constant for every possible scenario of _μ<sub>c</sub>_=[1,18]. In subsequent prediction steps, the _μ<sub>mean</sub>_ values from the iteration with the least _CAD_ and _WCAD_ scores in its previous prediction step replaces the unselected missing _μ_ elements.
3. In each prediction step, the _μ<sub>mean</sub>_ trend of the iteration with the least _CAD_ and _WCAD_ score, respectively, are carried over to the next prediction step. In the final prediction step, the _μ<sub>mean</sub>_ trend with the least _CAD_ and _WCAD_ score, respectively, are selected.
Consequently, Step0 (the 1<sup>st</sup> prediction step) performs 18<sup>5</sup>x18=1,889,568x18=34,012,224 iterations, Step1 performs 18<sup>5</sup>x2=1,889,568x2= 3,779,136 iterations and Step2 performs 18<sup>3</sup>x2=5,832x2=11,664 iterations, of _backcasting-forecasting_ and _CAD-WCAD_ treatments of their constituted _μ<sub>mean</sub>_ trends. These computations are performed for three unique random seeds to quantify the effects of statistical variances. Therefore in total, 37,803,024×3=113,409,072 iterations of _backcasting-forecasting_ and the _CAD-WCAD_ treatments of their constituted  _μ<sub>mean</sub>_ trends are solved. Figure 6 shows the value assigned to the unselected missing _μ<sub>mean</sub>_ elements via the _CAD_ and _WCAD_ criteria for different Random Seeds from Step0 can be similar and dissimilar. Table 1 evidence small improvement gains in _CAD_ and _WCAD_ step after step of the _Modified Resemblance Algorithm_, the values of _WCAD_ are one order smaller than _CAD_, and Seed2 yielded the least _CAD_ while Seed3 yielded the least _WCAD_.

Figure 7 illustrates the  _μ<sub>mean</sub>_ trends of Figure 5 after the _Modified Resemblance Algorithm_ treatment for three Random Seeds. Their 14 days windowed Simple Moving Averages (SMA) show:
1. Singapore started with the SMA daily _CCP_ of 9 days until ~1<sup>st</sup> March 2020. 
2. In the next five weeks, i.e. leading into the CB, this duration decreased to 4 days. 
3. Throughout the CB until 12<sup>th</sup> June 2020, this duration ranged from 4 to 8 days. This variation appears cyclical over a 4 to 5 weeks period.
4. Over the last 4 to 5 weeks, the SMA daily _CCP_ ranged between 5 to 9 days.
The 14 days window of the SMA reflects the self-isolation/quarantine period mandated by Singapore's Stay-Home-Notice (SHN) Order [26]. These SMA trends of  _μ<sub>mean</sub>_ evidenced the collective effort by Singapore to quickly confirm COVID-19 (given its SHN) started a month after the confirmation of its 1<sup>st</sup> Local COVID-19 case. This success then became periodic, fluctuating per month, during and post CB.

![Figure 6](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_6_MRA_Step0_leastCAD_leastWCAD_scores.png?raw=true)
**Figure 6:** The _least-CAD_ and _least-WCAD_ scores obtained in Step 0 of the _Modified Resemblance Algorithm_ for three unique random seeds.

![Table 1](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Table1_leastCAD_leastWCAD_of_Steps_0_1_2.png?raw=true)

**Table 1:** The _least-CAD_ and _least-WCAD_ scores of Steps 0, 1 and 2 for three unique random seeds.

![Figure 7](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_7_%CE%BCmean_completed.png?raw=true)
**Figure 7:** The mean of the estimates of _μ_, i.e. _μ<sub>mean</sub>_, for Random Seeds 1, 2 and 3 combined (see Figure 5) completed with the _μ<sub>mean</sub>_ predicted by the _Modified Resemblance Algorithm_. Included is also their 14 days windowed Simple Moving Averages.

### 3.3 The Estimated Lower-bound Local COVID-19 Epidemic Trends

The estimated lower-bound Local COVID-19 epidemic trends resembled their empirical counterpart (see Figures 8 and 9). According to Table 2, the _CAD_ scores from these estimated trends are ~21% of the Local COVID-19 population. This result means that these estimated trends predicted or resembled ~79% of the empirical data.

An advantage of the _least-WCAD_ criterion appears to be its ability to capture the peak of the Local COVID-19 epidemic trend (see Table 3). It predicted 1474 peak cases. This amount is closer to the actual 1426 cases than the 1226 cases estimated by the _least-CAD_ criterion. Moreover, it correctly predicted the day of the peak event. The _least-CAD_ criterion predictions were a day or two later.

With the strong resemblances to the empirical Local COVID-19 epidemic trend achieved by the estimated lower-bound Local COVID-19 epidemic trends, credence in the estimated lower-bound Local SARS-CoV-2 trends by the _least-CAD_ and _least-WCAD_ criteria is reasonable and presented next.

![Figure 8](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_8_lowerbound_Local_COVID19_epidemic_trends_leastCAD.png?raw=true)
**Figure 8:** Comparison of the lower-bound Local COVID-19 epidemic trends estimated by the _least-CAD_ criterion against its empirical counterpart.

![Figure 9](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_9_lowerbound_Local_COVID19_epidemic_trends_leastWCAD.png?raw=true)
**Figure 9:** Comparison of the lower-bound Local COVID-19 epidemic trends estimated by the _least-WCAD_ criterion against its empirical counterpart.

![Table 2](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Table2_leastCAD_leastWCAD_of_COVID19_epidemic_trends.png?raw=true)
**Table 2:** The _CAD_ and _WCAD_ scores of the lower-bound Local COVID-19 epidemic trends estimated by the _least-CAD_ and _least-WCAD_ criteria.

![Table 3](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Table3_peak_COVID19_cases_estimated_vs_empirical.png?raw=true)
**Table 3:** Data on the estimated and empirical peak daily number of confirmed COVID-19 cases.

### 3.4 The Estimated Lower-bound Local SARS-CoV-2 Infection Trends
## 4 Discussions
## 5 Conculsions
## 6 References



## Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

## Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/JulianChia/lowerboundSARSCOV2/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

## Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.

