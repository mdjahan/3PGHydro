# 3PG-Hydro V1.1
3PG-Hydro is available as an R-package, Coding by Anja Nölte & Marc Djahangard.

3PG-Hydro is an update of the original 3PG Forest Growth model by Landsberg and Waring (1997) (https://doi.org/10.1016/S0378-1127(97)00026-1). 3PG-Hydro calculates on a daily timestep, includes a soil-water-model, seperated soil evaporation calculation as well as a snow routine. Further information in the publication on 3PG-Hydro by Yousefpour and Djahangard (2021) (https://doi.org/10.3390/f12121729).

V1.1

We integrated a dynamic leaf fall and leaf growth model for decidious trees as well as dynamic CO2-concentration function that prdouces yearly values for historical period as well as climate change scenarios RCP2.6, RCP4.5 & RCP8.5. Moreover, 3PG-Hydro integrates the stand volume & height calculation method as described in Forrester et al. (2021) (https://doi.org/10.1007/s10342-021-01370-3). Therefore, the parameters from Forrester et al. (2021) can be used.
Moreover, V1.1 includes calculation of soil respiration to provide the output of net ecosystem exchange (NEE) adpated from Meyer et al. (2018) (https://doi.org/10.1016/j.foreco.2018.01.034) and a module for simulating a bark beetle attack adapted from Meyer et al. (2017) (https://doi.org/10.1016/j.foreco.2017.03.019).

V1.0

Schema 3PG-Hydro (as in Yousefpour & Djahangard (2021)):
![schema1](https://user-images.githubusercontent.com/122866605/213461877-833bb89d-31b6-4ca1-99d3-ef30e7fcff4f.png)
