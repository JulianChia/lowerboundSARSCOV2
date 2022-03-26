---
layout: default
---

**By Julian Chia, 1<sup>st</sup> March 2022**

_In remembrance of the 2<sup>nd</sup> Anniversary of Singapore's COVID-19 epidemic._

## Summary
1. A model to "extract" lower-bound estimates of Singapore's daily SARS-CoV-2 infection trend from its daily COVID-19 epidemic trend and a few statistical parameters from the confirmation period of a sample of its COVID-19 cases is successfully demonstrated.

2. Key findings:

   1. Singapore had an early window of opportunity to mitigate its COVID-19 epidemic with its Circuit Breaker, but it was unseized. An extended Circuit Breaker with tighter mitigation measures to quell Singapore's COVID-19 epidemic then needed implementation even having reduced its period to confirm COVID-19 cases.

   2. Many people with SARS-CoV-2 who ultimately had COVID-19 remained unidentified until a later point in time. I believe this factor is one of the reasons for the epidemic nature of COVID-19 in Singapore and its protracted recovery.

   3. The influx of Imported COVID-19 individuals in March 2020 caused Singapore's COVID-19 epidemic

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

   <img src="https://render.githubusercontent.com/render/math?math={\color{black} P(X) = \displaystyle\frac{1}{\sqrt{2\pi\sigma^2}} \exp^{-\frac{(X-\mu)^2}{2\sigma^2}}}"> ..........(1)

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
1. Iteratively perform the _backcasting_ and _forecasting_ algorithms for a range of constant _μ_ scenarios (termed as _μ<sub>c</sub>_).
2. Select _μ<sub>c</sub>_ as a probable value of _μ_ for a given day when the number of estimated Local COVID-19 confirmation events that day is identical to its empirical data.
3. Complete this selection procedure upon reaching 300 samples of the _μ_ estimates per day for the Local COVID-19 epidemic population.
The approximation of the range of _μ<sub>c</sub>_ in step 1 of this procedure is determined using empirical _CCP_ data. Figure 2 illustrates the _CCP_ from a sample of Local COVID-19 cases obtained from [3] and various national news sources. Its cumulative probability-density distribution is Gaussian-like with a mean and standard deviation of 17.531 days and 6.044 days, respectively. Assuming this empirical mean value indicates the maximum range of _μ<sub>c</sub>_ while letting a day be its minimum, then _μ<sub>c</sub>_=[1,18] since the ceiling value of 17.531 days is 18 days. Furthermore, _σ_ is assumed to be a quadratic function of _μ_:

   <img src="https://render.githubusercontent.com/render/math?math={\color{black} \sigma = f(\mu) = a\mu^2 %2B b\mu %2B c}"> ..........(2)

   where _a_ = -0.008665, _b_ = 0.483888, _c_ = 0.0 and _μ_=_μ<sub>c</sub>_. Figure 2 illustrates the Normal Cumulative Probability Distributions of _μ<sub>c</sub>_=[1,18].

The estimation of _μ_ via the abovementioned procedure will, unfortunately, continue indefinitely when one or more elements of _μ_ are unpredictable, i.e. when the element(s) of _μ_ have zero value. In such an eventuality, a _Resemblance Algorithm_ is to complete the estimation of _μ_. Firstly, the Cartesian Product of _μ<sub>c</sub>_=[1,18] replaces the element(s) of _μ_ that has zero value. Doing so allows the exploration of every possible sequencing of _μ<sub>c</sub>_ in these missing _μ_ elements. Next, the backcasting and forecasting of each possible sequencing of these _μ_ estimates the Local COVID-19 epidemic trend. Finally, the sequence of _μ_ that yields the Local COVID-19 epidemic trend that best resembles its empirical counterpart is its _μ_. The measurement of resemblance is via a _cumulative-absolute-difference_ (_CAD_) criterion:

<img src="https://render.githubusercontent.com/render/math?math={\color{black} CAD = \displaystyle\sum_{d=0}^{d_{max}} \left|{T_e-T_m}\right|}"> ..........(3)

and its _weighted_ counterpart (_WCAD_) criterion:

<img src="https://render.githubusercontent.com/render/math?math={\color{black} WCAD = \displaystyle\sum_{d=0}^{d_{max}} \left( \frac{T_e}{p} * \left|{T_e-T_m} \right| \right)}"> ..........(4)

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

Figures 10 and 11 plot the estimated lower-bound daily Local SARS-CoV-2 infections alongside the empirical Local COVID-19 epidemic trend. Table 4 presents data on their peak SARS-CoV-2 event. Together, they show that :

1. Local SARS-CoV-2 infection trends precede the Local COVID-19 epidemic trend.
2. The profile of the daily Local SARS-CoV-2 infection trends foreshadowed the growth and reduction of the daily Local COVID-19 cases.
3. Local SARS-CoV-2 infections had peaked on the same day as when the Local COVID-19 epidemic peaked, i.e. 20<sup>th</sup> April 2020.  The _least-WCAD_ criterion conservatively predicted that 2044 (±124 to ±183) daily SARS-CoV-2 infections had occurred that day. This amount is 618 (±124 to ±183) cases greater than the 1426 COVID-19 cases documented that day.
4. During the CB, a secondary peak SARS-CoV-2 infection event had occurred on 1<sup>st</sup> May 2020. A minimum of 1205 (±86 to ±129) to 1225 (±88 to ±132) Local SARS-CoV-2 infections occurred that day.

![Figure 10](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_10_lowerbound_Local_SARSCOV2_infection_trends_leastCAD.png?raw=true)
**Figure 10:** The lower-bound Local SARS-CoV-2 infection trends estimated by the _least-CAD_ criterion contrasted against the empirical Local COVID-19 epidemic trend.

![Figure 11](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_11_lowerbound_Local_SARSCOV2_infection_trends_leastWCAD.png?raw=true)
**Figure 11:** The lower-bound Local SARS-CoV-2 infection trends estimated by the _least-WCAD_ criterion contrasted against the empirical Local COVID-19 epidemic trend.

![Table 4](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Table4_peak_SARSCOV2_cases_estimated.png?raw=true)
**Table 4:** Data on the predicted peak Local SARS-CoV-2 Infection event.

## 4 Discussions

The lower-bound Local SARS-CoV-2 infection trends in Section 3.4 provide a novel reference to understand the Local COVID-19 epidemic of Singapore. For example, Singapore could have mitigated the COVID-19 epidemic in 2020 if it had implemented its CB a week or two sooner. Figures 10 and 11 supports such a view as it shows that the origin of the exponential transmission of SARS-CoV-2 in Singapore predated the CB announcement and implementation by approximately two and three weeks, respectively. Also, given that Singapore had announced on 21<sup>st</sup> April 2020 further tighter CB measures until 4<sup>th</sup> May 2020 and extended the CB dateline to 1<sup>st</sup> June 2020 [27], in hindsight, these mitigation measures do seem appropriate in scale and timing. This conclusion is arrived at because this announcement predates 1<sup>st</sup> May 2020, the day sustain reduction in the transmission rate of SARS-CoV-2 in Singapore began to show signs.

Figure 12 plots the cumulative number of lower-bound Local SARS-CoV-2 infections leading to four events in Singapore's COVID-19 epidemic history. Namely, (a) the announcement of the CB, (b) the start of the CB, (c) when the daily Local COVID-19 epidemic peaked (this day also marks when daily Local SARS-CoV-2 infections peaked) and (d) when the 2nd highest peak in daily Local SARS-CoV-2 infection occurred. The figure also plots the cumulative number of empirical Local COVID-19 epidemic cases. Included are also the cumulative case counts on the day of these events. Together, these pieces of information illustrate the significant amount of imminent COVID-19 individuals undocumented by the Local COVID-19 epidemic trend. For example, at the CB announcement, the cumulative number of Local COVID-19 cases was only 570. However, the reality then was that 1800(±15) to 1807(±39) SARS-CoV-2 infectees, who shall develop COVID-19 within a week or later and who shall require medical resources, had already existed. Thus, a minimum of 1230(±15) to 1237(±39) COVID-19 individuals were still undocumented that day. These undocumented imminent COVID-19 individuals nearly doubled when the CB started. There were at least 2070(±15) to 2070(±39) such undocumented individuals that day. By the time the COVID-19 epidemic trend had peaked, the population of these undocumented individuals more than doubled to reach 4816(±15) to 4837(±39) such individuals. Despite decreasing to 3208(±15) to 3425(±39) ~1.5 weeks later, this amount of undocumented imminent COVID-19 individuals is still significant. Figure 13 also evidence that these undocumented individuals existed throughout the Local COVID-19 epidemic trend. The inadequacy of the Local COVID-19 trend to document all carriers of SARS-CoV-2 that are guaranteed to develop COVID-19 in real-time is evident. It is a possible factor for the epidemic and its protracted recovery. The lower-bound Local SARS-CoV-2 trend avoids this inadequacy and is advantageous in this regard.

![Figure 12](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_12_SARSCOV2_vs_COVID19_population.png?raw=true)
**Figure 12:** The evolving population of the estimated lower-bound Local SARS-CoV-2 infection trends and the empirical Local COVID-19 epidemic trend leading to (a) the announcement of the CB, (b) the start of the CB, (c) when the daily Local COVID-19 epidemic trend peaked (it is also when the daily Local SARS-CoV-2 infection trend peaked) and (d) the 2<sup>nd</sup> highest peak of the daily Local SARS-CoV-2 infection trend.

![Figure 13](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_13_undocumented_imminent_COVID19_population_trend.png?raw=true)
**Figure 13:** Differences in the population of the estimated lower-bound Local SARS-CoV-2 infection trends and the empirical Local COVID-19 epidemic trend, i.e. the undocumented imminent Local COVID-19 population. Note, its demise towards 18<sup>th</sup> August, though technically correct, is an artefact of the termination effect of the data used. If post 18<sup>th</sup> August 2020 COVID-19 epidemic data are included in the modelling, the sizable population of undocumented imminent Local COVID-19 individuals will have propagated.

Another use of the lower-bound Local SARS-CoV-2 trends is in the forensic analysis of the Local COVID-19 epidemic origin. Figure 14 compares them against the empirical Imported COVID-19 epidemic trend. The steep rise in Local SARS-CoV-2 infections occurred when Imported COVID-19 cases rose in Singapore. These two events occurred concurrently; circumstantially, they appear related. The likelihood that the influx of Imported COVID-19 individuals into Singapore in March 2020 caused the exponential transmission of Local SARS-CoV-2 that led to the COVID-19 epidemic seems plausible. Comparisons of the arrival dates by the Imported COVID-19 individuals against the lower-bound Local SARS-CoV-2 infection trends can evaluate this plausibility. If the Imported COVID-19 individuals did not cause the Local COVID-19 epidemic, their arrival in Singapore would never precede the start of the exponential transmissions of Local SARS-CoV-2. The opposite rationale applies if the arrival of the Imported COVID-19 individuals to Singapore precedes the exponential Local SARS-CoV-2 infections.

According to published records [3], there were 569 Imported COVID-19 cases by 19th April 2020. Of which, 110 had their arrival date published while the remaining 459 did not and need to be estimated. Figure 15 shows that the cumulative probabilistic distribution of the confirmation period of the Imported cases can be modelled by Gaussian Kernel Density Estimation (_KDE<sub>gauss</sub>_) theory better than Normal Distribution theory. Thus, subtracting the Press-Release Date of each of the 459 Imported COVID-19 cases by the COVID-19 confirmation period predicted by _KDE<sub>gauss</sub>_ could reasonably approximate their arrival dates.

The daily arrival of the Imported COVID-19 cases (both empirical and estimated) are in Figure 16. A sizable sample of these Imported COVID-19 cases precedes by a week or two the start of the exponential rise in daily Local SARS-CoV-2 infections in March 2020. This lead time certainly gave opportunities for SARS-CoV-2 transmission to Singapore residents. The ease to flout the SHN by overseas arrivals at that time [28] further increases the dangers posed by these lead times to transmit SARS-CoV-2.

A prelude to the above events is that 11<sup>th</sup> March 2020 marked the pronouncement by the Director-General of The World Health Organisation that COVID-19 has evolved into a pandemic. Globally, the number of COVID-19 cases outside China had increased 13-fold, and the number of countries with COVID-19 cases had increased 3-fold [29]. In Singapore, this phenomenon resulted in an exodus of Singaporeans returning from many COVID-19 plagued countries that month. Minister Lawrence Wong [30], one of the co-chairs of Singapore’s COVID-19 Multi-Ministry Task Force, made known in a Parliament meeting on 25th March 2020 that the daily number of Imported COVID-19 cases would grow significantly. He described that there were about 200,000 Singaporeans abroad and that Singapore should be prepared to receive these returnees when they choose to weather out the COVID-19 pandemic in Singapore. Singapore had been witnessing ~1200 returnees from the UK and US daily. He feared the end of this high trend of returnees did not seem to be in sight. According to Figure 12, Singapore’s Local SARS-CoV-2 population that day had reached at least a mean of 589 with a standard deviation of 9.849, while the officially reported Local COVID-19 population was only 264. That is, the dangers posed by at least 325 imminent COVID19 individuals were still unrecognised that day.

Viewing all this information together with the circumstantial evidence presented in Figures 14 and 16 and their reasoning, it is not plausible to not attribute the surge of SARS-CoV-2 infections in Singapore, hence Singapore’s COVID-19 epidemic, to the influx of Imported COVID-19 individuals in March 2020.

![Figure 14](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_14_Imported_COVID19_vs_SARSCOV2_CAD%26WCAD.png?raw=true)
**Figure 14:** Singapore’s lower-bound daily Local SARS-CoV-2 infection trends vs its Imported COVID-19 epidemic trend.

![Figure 15](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_15_Imported_covid19_confirmation_period_new2.png?raw=true)
**Figure 15:** The cumulative probabilistic distribution of the COVID-19 confirmation period of the Imported COVID-19 cases.

![Figure 16](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_16_Imported_COVID19_arrival_vs_SARSCOV2_CAD%26WCAD.png?raw=true)
**Figure 16:** Singapore’s lower-bound daily Local SARS-CoV-2 infection trends vs the arrival trend of its Imported COVID-19 cases.

## 5 Conculsions

A novel approach to estimate the lower-bound Local SARS-CoV-2 infection trend of Singapore between January to August 2020 is presented. It involved first modelling Singapore's daily mean COVID-19 confirmation periods, i.e. μmean, and the lower-bound Local COVID-19 epidemic trend. The former showed Singapore's ability to reduce its daily COVID-19 confirmation periods. Yet, this measure alone was not sufficient to avert a COVID-19 epidemic. The latter achieved resemblance to the empirical lower-bound Local COVID-19 epidemic trend. These resemblances gave credence to the estimated lower-bound Local SARS-CoV-2 infection trends. These trends showed Singapore missing an early window of opportunity to mitigate its COVID-19 epidemic with its Circuit Breaker. As such, an extended Circuit Breaker with tighter mitigation measures to quell Singapore's COVID-19 epidemic needed implementation. They also provided the ability to quantify the daily imminent Local COVID-19 population within a week or more before their confirmation by the Local COVID-19 epidemic trend. This undocumented population is sizable and evidenced the COVID-19 epidemic trend weakness to report the "real" COVID-19 population (i.e. those infected with SARS-CoV-2 that are, and also have yet to be, confirmed with COVID-19) promptly; a possible factor for the COVID-19 epidemic and its protracted recovery. These trends also elucidated Singapore's COVID-19 epidemic origin. They provided circumstantial evidence that the exponential surge of SARS-CoV-2 in Singapore, hence Singapore's COVID-19 epidemic, was caused by the influx of Imported COVID-19 individuals in March 2020.

## 6 References

1. Gorbalenya AE, Baker SC, Baric RS, de Groot RJ, Drosten C, Gulyaeva AA, et al. The species Severe acute respiratory syndromerelated coronavirus: classifying 2019-nCoV and naming it SARS-CoV-2.  Nature Microbiology, 2020-03;5(4): 536–544. [https://doi.org/10.1038/s41564-020-0695-z](https://doi.org/10.1038/s41564-020-0695-z)
2. Confirmed imported case of novel coronavirus infection in Singapore; multi-ministry taskforce ramps up precautionary measures, Ministry of Health of Singapore, Updates on COVID-19 (Coronavirus Disease 2019) Local Situation, 23rd Jan 2020. [https://www.moh.gov.sg/news-highlights/details/confirmed-imported-case-of-novel-coronavirus-infection-in-singapore-multi-ministry-taskforce-ramps-up-precautionary-measures](https://www.moh.gov.sg/news-highlights/details/confirmed-imported-case-of-novel-coronavirus-infection-in-singapore-multi-ministry-taskforce-ramps-up-precautionary-measures)
3. Updates on COVID-19 (Coronavirus Disease 2019) Local Situation, Ministry of Health of Singapore. [https://www.moh.gov.sg/covid-19](https://www.moh.gov.sg/covid-19)
4. Ng Y, Li Z, Chua YX, et al. Evaluation of the effectiveness of surveillance and containment measures for the first 100 patients with COVID-19 in Singapore—January 2–February 29, 2020. MMWR Morb Mortal Wkly Rep 2020;69:307–11. [https://doi.org/10.15585/mmwr.mm6911e1](https://doi.org/10.15585/mmwr.mm6911e1)
5. Singapore COVID-19 Epidemic Trends:   Click on the “Number of Cases” tab in  [https://covidsitrep.moh.gov.sg/](https://covidsitrep.moh.gov.sg/)
6. Wu Z, McGoogan JM. Characteristics of and Important Lessons From the Coronavirus Disease 2019 (COVID-19) Outbreak in China: Summary of a Report of 72 314 Cases From the Chinese Center for Disease Control and Prevention. JAMA. 2020;323(13):1239–1242. [https://doi.org/10.1001/jama.2020.2648](https://doi.org/10.1001/jama.2020.2648)
7. Atkinson, B. and Petersen, E., SARS-CoV-2 shedding and infectivity. The Lancet, 2020-04-15;395(10233):1339-1340. [https://doi.org/10.1016/S0140-6736(20)30868-0](https://doi.org/10.1016/S0140-6736(20)30868-0)
8. Bullard J, Dust K, Funk D, Strong JE, Alexander D, Garnett L, et al. Predicting infectious SARS-CoV-2 from diagnostic samples. Clin Infect Dis. 2020-11-15;71(10):2663–2666, [https://doi.org/10.1093/cid/ciaa638](https://doi.org/10.1093/cid/ciaa638)
9. Wei WE, Li Z, Chiew CJ, Yong SE, Toh MP, Lee VJ. Presymptomatic Transmission of SARS-CoV-2 — Singapore, January 23–March 16, 2020. MMWR Morb Mortal Wkly Rep 2020;69:411–415. DOI: [https://dx.doi.org/10.15585/mmwr.mm6914e1](https://dx.doi.org/10.15585/mmwr.mm6914e1)
10. Qian G, Yang N, Ma AHY, et al, A COVID-19 Transmission within a family cluster by presymptomatic infectors in China. Clin Infect Dis, 2020-07-28;71(15):861-862. [https://doi.org/10.1093/cid/ciaa316](https://doi.org/10.1093/cid/ciaa316).
11. Kimball A, Hatfield KM, Arons M, et al. Asymptomatic and presymptomatic SARS-CoV-2 infections in residents of a long-term care skilled nursing facility—King County, Washington, March 2020. MMWR Morb Mortal Wkly Rep 2020. Epub March 27, 2020. [https://www.cdc.gov/mmwr/volumes/69/wr/mm6913e1.htm](https://www.cdc.gov/mmwr/volumes/69/wr/mm6913e1.htm)
12. How long does it take to develop symptoms?, World Health Organisation, Coronavirus disease (COVID-19) Q&A, 12 Oct 2020, [https://www.who.int/emergencies/diseases/novel-coronavirus-2019/question-and-answers-hub/q-a-detail/coronavirus-disease-covid-19](https://www.who.int/emergencies/diseases/novel-coronavirus-2019/question-and-answers-hub/q-a-detail/coronavirus-disease-covid-19)
13. Lauer SA, Grantz KH, Bi Q, et al, The Incubation Period of Coronavirus Disease 2019 (COVID-19) From Publicly Reported Confirmed Cases: Estimation and Application, Annals of Internal Medicine, 2020-05-5. [https://doi.org/10.7326/M20-0504](https://doi.org/10.7326/M20-0504)
14. Ridgway, J., Shah, N., & Robicsek, A., Prolonged shedding of severe acute respiratory coronavirus virus 2 (SARS-CoV-2) RNA among patients with coronavirus disease 2019 (COVID-19). Infect Control Hosp Epidemiol, 2020-06-24;41(10):1235-1236. [https://doi.org/10.1017/ice.2020.307](https://doi.org/10.1017/ice.2020.307)
15. Fontana, L., Villamagna, A., Sikka, M., & McGregor, J., Understanding viral shedding of severe acute respiratory coronavirus virus 2 (SARS-CoV-2): Review of current literature. Infect Control Hosp Epidemiol, 2020-10-20;FirstView:1-10. [https://doi.org/10.1017/ice.2020.1273](https://doi.org/10.1017/ice.2020.1273)
16. Transmission of SARS-CoV-2: implications for infection prevention precautions, World Health Organization, Scientific Brief - 9 July 2020. [https://www.who.int/news-room/commentaries/detail/transmission-of-sars-cov-2-implications-for-infection-prevention-precautions](https://www.who.int/news-room/commentaries/detail/transmission-of-sars-cov-2-implications-for-infection-prevention-precautions)
17. Qiu J, Covert coronavirus infections could be seeding new outbreaks. Nat News 2020-03-20. [https://doi.org/10.1038/d41586-020-00822-x](https://doi.org/10.1038/d41586-020-00822-x)
18. Shi, Q., Hu, Y., Peng, B. et al. Effective control of SARS-CoV-2 transmission in Wanzhou, China. Nat Med 2021;27:86-93. [https://doi.org/10.1038/s41591-020-01178-5](https://doi.org/10.1038/s41591-020-01178-5)
19. Subramanian R, He Q, Pascual M, Quantifying asymptomatic infection and transmission of COVID-19 in New York City using observed cases, serology, and testing capacity. PNAS 2021-03-02;118(9):e2019716118. [https://doi.org/10.1073/pnas.2019716118](https://doi.org/10.1073/pnas.2019716118)
20. Johansson MA, Quandelacy TM, Kada S, et al. SARS-CoV-2 Transmission From People Without COVID-19 Symptoms. JAMA Netw Open. 2021;4(1):e2035057. [https://doi.org/10.1001/jamanetworkopen.2020.35057](https://doi.org/10.1001/jamanetworkopen.2020.35057)
21. Python 3, [https://www.python.org/](https://www.python.org/)
22. NumPy, [https://numpy.org/](https://numpy.org/)
23. SciPy, [https://scipy.org/](https://scipy.org/)
24. Matplotlib, [https://matplotlib.org/stable/index.html](https://matplotlib.org/stable/index.html)
25. Model’s source code in Github repository [https://github.com/JulianChia/lowerboundSARSCOV2](https://github.com/JulianChia/lowerboundSARSCOV2)
26. Implementation of new Stay-Home-Notice,  Ministry of Health of Singapore. [https://www.moh.gov.sg/news-highlights/details/implementation-of-new-stay-home-notice](https://www.moh.gov.sg/news-highlights/details/implementation-of-new-stay-home-notice)
27. PM Lee announced tighter measures as well as an extension to the Circuit Breaker until 1 June 2020. Published on 21 Apr 2020. [https://www.gov.sg/article/pm-lees-address-on-the-covid-19-situation-in-singapore-21-april-2020](https://www.gov.sg/article/pm-lees-address-on-the-covid-19-situation-in-singapore-21-april-2020)
28. Shaffiq A. Man who breached Covid-19 stay-home notice to eat bak kut teh convicted; DPP argues for 10- to 12-week jail term, 16 Apr 2020, [https://www.straitstimes.com/singapore/courts-crime/man-who-breached-stay-home-notice-to-eat-bak-kut-teh-convicted-dpp-argues-for?cx_testId=20&cx_testVariant=cx_1&cx_artPos=4#cxrecs_s~](https://www.straitstimes.com/singapore/courts-crime/man-who-breached-stay-home-notice-to-eat-bak-kut-teh-convicted-dpp-argues-for?cx_testId=20&cx_testVariant=cx_1&cx_artPos=4#cxrecs_s~)
29. Cucinotta D, Vanelli M. WHO Declares COVID-19 a Pandemic. Acta Biomed. 2020 Mar 19;91(1):157-160. doi: 10.23750/abm.v91i1.9397. PMID: 32191675; PMCID: PMC7569573. [https://pubmed.ncbi.nlm.nih.gov/32191675/](https://pubmed.ncbi.nlm.nih.gov/32191675/)
30. In Parliament: Lawrence Wong on ‘critical phase’ of COVID-19 fight as Singaporeans return home. 25th Mar 2020. [https://www.youtube.com/watch?v=\_-3u63Ve-C4](https://www.youtube.com/watch?v=_-3u63Ve-C4)

