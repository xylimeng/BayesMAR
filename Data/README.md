# Read-world data
Data used in real data application, collected on Aug 9th 2018 

## Producer Price Index 
- tp.csv: Raw data shown in Figure 2.
- diffp,csv: Lagged data of order one shown in Figure 3, left column.
- original source: https://fred.stlouisfed.org/series/PPIACO

## Treasury Bill Rate
- tr.csv: Raw data shown in Figure 2.
- diffr.csv: Lagged data of order one shown in Figure 3, left column.
- original source: https://fred.stlouisfed.org/series/TB3MS


## Unemployment Rate
- tu.csv: Raw data shown in Figure 2.
- diffu.csv: Lagged data of order one shown in Figure 3, left column.
- original source: https://fred.stlouisfed.org/series/LRUN64TTUSQ156N

## TBR_BayesMAR_forecast.RData
- preliminary results for demonstration

## CRPS_tem
- samples to calculate CRPS for BayesMAR-MAP at step #6
- CRPS for BayesMAR-BMA could be calculated in similar way, while requiring samples at each time point for all order
