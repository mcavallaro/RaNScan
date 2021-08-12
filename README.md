
# Cluster detection with Random Neighbourhood Covering

This repository contains code used in reference [1]. Please cite [1] if you find it useful.

As a demonstration, we consider the synthetic data in `case.df` which have been simulated from baseline `b.matrix` (see `Simulation_experiment4baseline.Rmd`).

```
> case.df = read.table(file='Data/Simulation_experiment_case_df.csv',  sep=',', header=T)> 
> b.matrix = 
> df.cases2 = 
```

```
> head(case.df)
  row SAMPLE_DT_numeric postcode latitude longitude population        y           x
1  483                 0  BB120EZ 53.79724 -2.264622        157 5982.253  -67.346355
2  488                 0  BB2 1HN 53.74896 -2.496710        113 5976.884  -74.560339
3  555                 0  BB8 7AR 53.86487 -2.163336        109 5989.774  -63.955308
4  973                 0  BN274EW 50.88265  0.263503        234 5658.151    9.803753
5 1079                 0  BS106DD 51.50526 -2.603483         37 5727.385  -92.748101
6 1142                 0  BS247EQ 51.34984 -2.917601        164 5710.103 -105.092027
```

```
> b.matrix[1:6,1:6]
                 0           1           2           3
AL1 1TA 0.008866186 0.009163252 0.009446268 0.009711113
AL1 1UB 0.005344003 0.005523056 0.005693641 0.005853273
AL1 2JT 0.002793456 0.002887052 0.002976221 0.003059666
```

```
> head(df.cases2)
        Postcode latitude longitude n.cases
AL1 1TA  AL1 1TA 51.74220 -0.320578       4
AL1 1UB  AL1 1UB 51.73746 -0.316577       1
AL1 2JT  AL1 2JT 51.73965 -0.340915       2
AL1 4XG  AL1 4XG 51.76646 -0.317093       6
AL1 5DF  AL1 5DF 51.75159 -0.306273       4
AL1 5JQ  AL1 5JQ 51.74941 -0.310961       4
```


Create  100,000 cyclinders to cover the detected cases:
```
> cylinders = CreateCylinders(observation.matrix = sim2, baseline.matrix = b.matrix,
emmtype = 'sim', week.range = c(0,99), n.cylinders = 100000, coord.df=df.cases2)
```

Compute the warning score for each case:
```
> source('R/Init2.R')
> source('R/surveillance_utilis.R)
> case.df[,'warning.score'] = apply(case.df, 1, FUN=warning.score, cylinders)
```

[1] M. Cavallaro, J. Coelho, D. Ready, V. Decraene, T. Lamagni, N. D. McCarthy, D. Todkill, M. J. Keeling,
Cluster detection with random neighbourhood covering: application to invasive Group A Streptococcal disease, 2021


