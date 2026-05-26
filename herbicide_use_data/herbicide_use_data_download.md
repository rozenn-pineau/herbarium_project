### County-level herbicide data use, 1992 - 2019

Dowloading the herbicide use information from 1992 to 2019 from the [Pesticide National Synthesis Project](https://water.usgs.gov/nawqa/pnsp/usage/maps/show_map.php?year=1992&map=GLYPHOSATE&hilo=L).

```
for year in {1992..2019}
do
    wget "https://water.usgs.gov/nawqa/pnsp/usage/maps/county-level/PesticideUseEstimates/EPest.county.estimates.${year}.txt"
done
```



### State-level herbicide use data 1992-2017

From : https://www.sciencebase.gov/catalog/item/5e95c13a82ce172707f25252 

This data release provides state-level estimates of annual agricultural use of pesticide compounds by major crop or crop group for states 
in the conterminous United States, for the time period 1992-2017, compiled from data used to make county-level estimates by means of methods 
described in Thelin and Stone (2013) and Baker and Stone (2015). The source of these data is the same as the published county-level pesticide use 
estimates for 1992-2009 (Stone, 2013), estimates for 2008-2012 (Baker and Stone, 2015), and estimates for 2013-17 (Wieben, 2019). 
County-level by-crop estimates are not published because of the increased uncertainty in estimating the geographic distribution of 
compounds applied to specific crops. High-acreage crops (corn, soybeans, wheat, cotton, rice, and alfalfa) are individually aggregated to state 
level while low-acreage crops are combined into groups (vegetables and fruit, orchards and grapes, pasture and hay, and other crops) prior to
aggregating to the state level.

This data release contains two tables of state-level annual agricultural pesticide use estimates by crop or crop group 
(one for low estimates and one for high estimates) and associated metadata. These data were used to produce annual time-series charts 
for individual pesticide by crop or crop group for 1992-2017 available on the Pesticide National Synthesis Project (PNSP) webpage: 
https://doi.org/doi:10.5066/F7NP22KM.

