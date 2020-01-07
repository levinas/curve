

This program fits dose response curves for one or multiple drug screening experiments using a Hill Slope model with 3 bounded parameters [1].
The program depends on the [Uno data loader](https://github.com/ECP-CANDLE/Benchmarks/blob/master/Pilot1/Uno/uno_data.py) in the [CANDLE](https://github.com/ECP-CANDLE) project.

The following command computes the results on an example melanoma cell line treated with paclitaxel.

```
python curve.py --cell LOXIMVI --drug paclitaxel
```

![LOXIMVI treated with paclitaxel](/figs/LOXIMVI-paclitaxel.png?raw=true)

The parameters and metrics of the model are as follows:
* Einf - fraction of cells not susceptible to drug
* EC50 - drug dose to have half target receptors bound
* HS - hill slope binding cooperativity
* AUC - area under growth curve for a fixed dose range: [4, 10] (-log10(M))
* IC50 - drug dose to have 50% growth
* EC50se - standard error of the estimated EC50
* R2fit - R2 score between the unclipped, real growths and fitted growth values
* AUC1 - area under growth curve for the measured dose range in a study
* AAC1 - area above growth curve for the measured dose range in a study
* DDS1 - a more robust variant of AUC1 normalized against dose range [2]

![LOXIMVI treated with paclitaxel](/figs/LOXIMVI-paclitaxel-table.png?raw=true)

For all cell line datasets, drug sensitivity is quantified by dose response values that measure the ratio of surviving treated to untreated cells after exposure to a given drug at a given concentration.
To facilitate comparison across cell lines, all response values were linearly rescaled to share a common range (from 0 to 100).
PharmacoDB incorporated multiple dose-independent metrics to allow integration of heterogeneous drug response data from the different studies.
We applied the same method to dose response data from five datasets: NCI60, CCLE, CTRP, GDSC, and gCSI.

We look at two confidence metrics to evaluate the quality of the fit: EC50se and R2fit.
In the examples in Table Y, NCI60 rep 2, CTRP rep1 and CTRP rep2 are pathological cases of curve fitting as they have high standard deviation of estimated EC50 or low R2 score of fit or both.
IC50 is not always applicable as certain drugs are not effective against cell lines regardless of concentration.

To process the integrated multi-study dose response data frame, use the default command:
```
python curve.py
```

### Options
```
usage: curve.py [-h] [-c CELL] [-d DRUG] [-s STUDY [STUDY ...]] [--reps REPS]
                [--start START] [--count COUNT] [--out OUT]

optional arguments:
  -h, --help            show this help message and exit
  -c CELL, --cell CELL  cell line name
  -d DRUG, --drug DRUG  drug name
  -s STUDY [STUDY ...], --study STUDY [STUDY ...]
                        list of data sources/studies
  --reps REPS           maximum of replicates to include in plots
  --start START         start index (0-based)
  --count COUNT         number of experiments
  --out OUT             prefix of output file
```

For more details, see the examples in the [Jupyter notebook](https://github.com/levinas/curve/blob/master/response-curve-fit.ipynb). 

### References
1. Petr Smirnov, Victor Kofia, Alexander Maru, Mark Freeman, Chantal Ho, Nehme El-Hachem, George-Alexandru Adam, Wail Ba-alawi, Zhaleh Safikhani, and Benjamin Haibe-Kains. Pharmacodb: an integrative database for mining in vitro anticancer drug screening studies. Nucleic acids research, 46(D1):D994–D1002, 2017.
2. Bhagwan Yadav, Tea Pemovska, Agnieszka Szwajda, Evgeny Kulesskiy, Mika Kontro, Riikka Karjalainen, Muntasir Mamun Majumder, Disha Malani, Astrid Murumägi, Jonathan Knowles, et al. Quantitative scoring of differential drug sensitivity for individually optimized anticancer therapies. Scientific reports, 4:5193, 2014.
