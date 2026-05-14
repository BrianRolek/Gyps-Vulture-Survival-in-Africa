This repository contains supplementary files (data and scripts) for implementation of survival analyses of *Gyps* vultures in eastern Africa. Analyses included multi-event capture-recapture survival models implemented using R, NIMBLE, and NIMBLEECOLOGY. The accompanying manuscript is:

B. W. Rolek, L. Dunn, C. J. Kendall, C. J. W. McClure, S. Thomsett, D. Muteti, M. Z. Virani, E. R. Buechley, R. Buij. 2026. Low survival of *Gyps* vultures in eastern Africa coinciding with interventions to address poisoning.

A full workflow is provided in docs/index.html or can be viewed at <https://brianrolek.github.io/Gyps-Vulture-Survival-in-Africa/>.

Data are included in the "data" folder as .rdata files. R and NIMBLE code are included in the "R" folder. Use load("data\\data.RData") to access the data in R and the data are saved within the list objects "datl" and "constl". Additional details are included within the manuscript and appendix.

datl : a list object

- *y*: the observation state of each *Gyps* vulture (rows) during each time step (columns). Observation states can range between one and five and are fully described in the methods. A matrix of dimensions *nind* x *ntime*.

- *first_age*: the estimated age of each vulture during first capture. NAs indicate that the vulture was a subadult that could not be accurately aged. The NAs are imputed by using a submodel for age class. A vector of length *nind*.

constl: a list object

- *nind*: the total number of individual *Gyps* vultures included in analysis.

- *ntime*: the total number of discrete time steps where each step is one month.

- *nyears*: the total number of years.

- *f*: The monthly time step when a vulture is captured, tagged, and released (sometimes release occurred after rehabilitation). A vector of length *nind*.

- *last*: The last monthly time step when a vulture was observed as a fatality or the time step when a vulture disappeared from tracking. A vector of length *nind*.

- *period.cat*: an explanatory variable describing which study period a vulture was observed within (0 = 2009-2011, 1= 2017-2024). A vector of length *ntime*.

- *year.cont*: year as a continuous explanatory variable, centered so the median is zero and scaled between -1 and 1. A vector of length *ntime*.

- *rehabbed:* a binary explanatory variable where 0 = not rehabilitated and 1 = rehabilitated. A vector of length *nind*.

- *tag.age.sc*: Age of transmitter for each vulture during each time step. A matrix of dimensions *nind* x *ntime*. Scaled and centered.

- *sp*: a binary explanatory variable of species where 0 = African white-backed vulture and 1 = Ruppell's vulture. A vector of length *nind*.

- *known*: a binary index where 1 = known age and 2 = unknown age. A vector of length *nind*.

- *y.first*: the observation state during the first time interval for each vulture. A vector of length *nind*.

  Note that the binary explanatory variables *period.cat*, *rehabbed*, and *sp* are converted within the model so that values of 0 = -1 and 1=1, thereby centering the covariate on zero and allowing the intercept to be interpreted as the overall mean.
