# BBS_Seasonal
Playing with a BBS model that also accounts for the seasonal time of survey

The BBS field methods do a reasonably good job of standardizing the seasonal day of a survey. Observers try hard to survey their routes at approximately the same time each year. 
In addition, looking at the existing data shows there's been little to know change in survey day over time. 

However, it is possible that the changes in abundance over time may be confounded with changes in survey date, particularly in some relatively data-sparse regions. 
It is also possible that the bird's phenology has changed over time. 

This project represents a start to explicitly model the observation-component of a day-of-year-effects (season, within-year, not sure which term makes the most sense).
 It also represents an example of a more general approach to adding flexible effort covariates to the BBS model.
