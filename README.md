# 3-PG-Hydro V2
New code will be released soon.
Major update of version two is the integration of 3-PGmix (Forrester, 2014 (https://doi.org/10.1016/j.ecolmodel.2013.12.021); Forrester et al., 2014 (https://doi.org/10.1186/s40663-014-0017-0); Forrester & Tang, 2016 (https://doi.org/10.1016/j.ecolmodel.2015.07.010); Trotsiuk et al., 2020 (https://doi.org/10.1111/2041-210X.13474)).

## V1.1
3-PG-Hydro is available as an R-package, Coding by Anja Nölte & Marc Djahangard.

3-PG-Hydro is an update of the original 3-PG Forest Growth model by Landsberg and Waring (1997) (https://doi.org/10.1016/S0378-1127(97)00026-1). 3-PG-Hydro calculates on a daily timestep, includes a soil-water-model, seperated soil evaporation calculation as well as a snow routine. Further information in the publication on 3-PG-Hydro by Yousefpour and Djahangard (2021) (https://doi.org/10.3390/f12121729).

## V1.1

We integrated a dynamic leaf fall and leaf growth model for decidious trees based on the Growing-Degreee-Days method adapted from Fu et al. (2014) (https://doi.org/10.1371/journal.pone.0109544). Moroever, we developed a dynamic CO2-concentration function that prdouces yearly values for the historical period as well as climate change scenarios RCP2.6, RCP4.5 & RCP8.5. Moreover, 3-PG-Hydro integrates the stand volume & height calculation method as described in Forrester et al. (2021) (https://doi.org/10.1007/s10342-021-01370-3). Therefore, the parameters from Forrester et al. (2021) can be used.
Moreover, V1.1 includes calculation of soil respiration to provide the output of net ecosystem exchange (NEE) adpated from Meyer et al. (2018) (https://doi.org/10.1016/j.foreco.2018.01.034) and a module for simulating a bark beetle attack adapted from Meyer et al. (2017) (https://doi.org/10.1016/j.foreco.2017.03.019).

## V1.0

Schema 3-PG-Hydro (as in Yousefpour & Djahangard (2021)):
![schema1](https://user-images.githubusercontent.com/122866605/213461877-833bb89d-31b6-4ca1-99d3-ef30e7fcff4f.png)
