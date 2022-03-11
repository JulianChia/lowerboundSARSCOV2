# Singapore’s SARS-CoV-2 Infection Trend In 2020: A Lower Bound Estimate

By Julian Chia, 1<sup>st</sup> March 2022

_In remembrance of the 2nd Anniversary of Singapore's COVID-19 epidemic._

## Abstract
Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2) infection causes Coronavirus Disease 2019 (COVID-19). Knowledge on the SARS-CoV-2 infection trend is, however, lacking. This viral infection is invisible to the naked eye and is challenging to profile in real-time. Its closest indication is its documented COVID-19 epidemic trend. Fortunately, that is published daily and globally due to the pandemic situation of COVID-19. Leveraging the empirical nature and availability of these COVID-19 epidemic trends, this paper posits that these trends are, in fact, lower-bound estimates of the respective localities COVID-19 epidemic situation, and each is interrelated to the lower-bound SARS-CoV-2 infection trend in their locality. A model based on this posit is developed and applied to Singapore. The predicted Local SARS-CoV-2 infection trends provided a novel reference to understanding the Singapore COVID-19 epidemic that was previously not possible. They evidenced the window of opportunity where Singapore could have mitigated its COVID-19 epidemic via its Circuit Breaker (CB) that it had missed. They evidenced Singapore's additional tighter CB measures and extended CB dateline, both implemented during the CB on 21<sup>st</sup> April 2020, were timely and effective. They derived the population of imminent COVID-19 individuals that the empirical Local COVID-19 epidemic trend had undocumented daily. These undocumented populations are sizable and a possible factor for the COVID-19 epidemic and its protracted recovery. Finally, these SARS-CoV-2 trends provided circumstantial evidence that Singapore's COVID-19 epidemic originated from COVID-19 cases imported into Singapore.

## Keywords
SARS-CoV-2 infection trend· COVID-19 epidemic trend· COVID-19 Confirmation Period · Normal Distribution · Cartesian Product · Backcasting · Forecasting · Lower-bound · Modelling · Python Programming  · 

## 1 Introduction

Singapore reported its first two Coronavirus Disease 2019 (COVID-19) [1] cases on 23<sup>rd</sup> January 2020 [2]. This disease onset evolved from several independent sporadic outbreaks into an epidemic within that year. Publicly accessible daily COVID-19 Case Reports and Situation Reports issued by the Ministry of Health of Singapore [3] describe the extent of this epidemic. The _number of daily confirmed COVID-19 cases_ is their unit of measurement of daily COVID-19 prevalence; they are the daily cumulative positive real-time reverse transcription-polymerase chain reaction (rt-PCR) tests and serologic assays [4]. Its daily compilation constructs the Singapore COVID-19 epidemic curves [5] for Imported and Local cases.

[Figure 1](https://github.com/JulianChia/lowerboundSARSCOV2/blob/main/1_Figures/Figure_1_SG_COVID19_%20Epidemic_trends.png?raw=true)
  
Figure 1 illustrates these Imported and Local COVID-19 epidemic trends of Singapore (from 23<sup>rd</sup> January to 18<sup>th</sup> August 2020). It evidences the pervasiveness of the COVID-19 epidemic in Singapore that caused the implementation of a nationwide Circuit Breaker (CB), alternatively known as a Lockdown, that year. Although these trend curves are informative, they do not describe the events that caused COVID-19, i.e. they do not show when Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2) infection had occurred. The SARS-CoV-2 infection trend curve fulfils this function, yet it is unreported. A means to monitor SARS-CoV-2 transmission in real-time is unavailable. Also, COVID-19 patients could not see and did not know when they had contracted SARS-CoV-2.
  
There is unanimous consensus that the SARS-CoV-2 infection trend is not identical to the COVID-19 epidemic trend. One reason for this is that SARS-CoV-2 infection is a precursor of COVID-19. SARS-CoV-2 infection takes time to cause COVID-19 and takes time to be discovered and confirmed [6] by the patients, their doctors and contact tracers. Furthermore, COVID-19 confirmation can occur during the incubation, symptomatic and recovery periods of a SARS-CoV-2 infection as positive rt-PCR results are obtainable during these periods regardless of the virus transmissibility [7,8]. Presymptomatic SARS-CoV-2 transmissions, which occurs during the incubation period, have been evidenced in Singapore [9], China [10], and the USA [11]. In 2020, the incubation period of SARS-CoV-2 variants averaged at 5 to 6 days and ranged from 1 to 14 days [12,13]. Thereon, positive rt-PCR have evidenced too. A multi-centre USA study of 70,406 unique COVID-19 patients [14] found that a majority of COVID-19 patients yielded positive rt-PCR, i.e. shed SARS-CoV-2 ribonucleic acid (RNA), three weeks after their first positive SARS-CoV-2 rt-PCR test. Another study that reviewed 77 COVID-19 reports [15] found that the duration of SARS-CoV-2 RNA shedding can range from a minimum of 1 day to a maximum of 83 days. According to WHO [16], rt-PCR positivity generally appears up to 3 weeks or more for mild to moderate COVID-19 and a more extended period for severe COVID-19. These findings evidence the amount of lag that a COVID-19 epidemic trend could have against a SARS-CoV-2 infection trend.
	
The second manner where the SARS-CoV-2 infection trend is dissimilar to the COVID-19 epidemic curve is in their population. The  SARS-CoV-2 infectee population, in reality, is more than the documented COVID-19 epidemic population. The under-reporting of asymptomatic and presymptomatic SARS-CoV-2 transmissions cause this phenomenon [17]. In Wanzhou, China, a 14-weeks COVID-19 mass testing program [18] found that asymptomatic and presymptomatic SARS-CoV-2 transmission accounted for 75.9% of all SARS-CoV-2; the abundance of close contacts before symptom onset or diagnosis facilitated them. Using a model that incorporates daily testing information fit to the COVID-19 case and serology data from New York City, researchers estimated that presymptomatic and asymptomatic SARS-CoV-2 transmissions together comprised at least 50% of the infections at the outbreak peak [19]. A decision analytical model study estimated that asymptomatic and presymptomatic SARS-CoV-2 accounted for >50% of all transmissions [20]. These studies exemplify the extent of under-reporting possibility occurring between the documented COVID-19 epidemic population against the actual SARS-CoV-2 infectee population.
	
Succinctly put, the actual SARS-CoV-2 infection trend will always appear ahead and exhibit a larger integral than its empirical COVID-19 epidemic trend. Its modelling is non-trivial to perform; the number of transmission factors for SARS-CoV-2 is many, their quantification remains challenging to implement, and SARS-CoV-2 transmission is opportunistic. To simplify its modelling while retaining an empirical basis, I posit that the officially published Local COVID-19 epidemic trend is, in fact, a lower-bound estimate of the Local COVID-19 epidemic and that it is a result of the lower-bound Local SARS-CoV-2 infection trend. Such a postulation is the basis of this paper.
	
The following sections present a model that utilises the abovementioned postulation, statistics, Cartesian-product, computation algorithm and empirical Local COVID-19 epidemic data to estimate the lower-bound Local SARS-CoV-2 infection trend of Singapore from January to August 2020. Its results contain the findings on the estimation of Singapore's daily mean Local COVID-19 confirmation period trend, lower-bound daily Local COVID-19 epidemic trend, and lower-bound daily Local SARS-CoV-2 infection trend. Their discussions give new insights into Singapore's COVID-19 epidemic situation in 2020.
  
## 2 The Model

The primary assumptions of the model are:

1. Local SARS-CoV-2 infections that are presymptomatic, asymptomatic or symptomatic but are undocumented are negligible. Therefore, an empirical Local COVID-19 epidemic trend is a reasonable lower-bound estimate of the Local COVID-19 epidemic trend and is relatable to the lower-bound Local SARS-CoV-2 infection trend that precedes it.

2. The Local COVID-19 confirmation event always lags behind its Local SARS-Cov-2 infection event. The duration between these two events is called the _COVID-19 confirmation period (CCP)_. It follows the Normal/Gaussian distribution theory. Daily, the _CCP_ probability density function is:
                             (1)
			     
 Here, _X є_ (-_∞_,_∞_) is the random variate of Eqn(1) and denotes the daily CCP, _μ_ denotes its daily mean, _σ_ denotes its daily standard deviation, and _π_ denotes the pi constant. _μ_ is not constant over the Local COVID-19 population history, i.e. _μ_ is a function of days.

### 2.1 Assumptions
### 2.2 Hypothesis
### 2.3 Methodology
### 2.4 Computation
## 3 Results
### 3.1 The Statistical μ Estimates
### 3.2 The Estimated μ<sub>mean</sub> Trend
### 3.3 The Estimated Lower-bound Local COVID-19 Epidemic Trends
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

