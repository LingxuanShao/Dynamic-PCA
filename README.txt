This directory contains the online supplements  (R codes).

## Datasets
The datasets used in this paper can be loaded from the following links from World Bank Website
Military expenditures: {https://data.worldbank.org/indicator/MS.MIL.XPND.GD.ZS}
Gross capital formation: {https://data.worldbank.org/indicator/NE.GDI.TOTL.ZS}
Imports of goods and services: {https://data.worldbank.org/indicator/NE.IMP.GNFS.ZS}
Exports of goods and services: {https://data.worldbank.org/indicator/NE.EXP.GNFS.ZS}
The four csv files **API_MS.MIL.XPND.GD.ZS_DS2_en_csv_v2_2055780.csv**, **API_NE.EXP.GNFS.ZS_DS2_en_csv_v2_2058633.csv**, **API_NE.GDI.TOTL.ZS_DS2_en_csv_v2_2252056.csv** and **API_NE.IMP.GNFS.ZS_DS2_en_csv_v2_2059780.csv** include the data from 1960 to 2019



## DCA code
The folder **DCA_code** provides necessary R code to implement the numerical experiments in this paper
The R script **simulation_onestep_unrolling.R** contains the assessment of our proposed method in Section 5.1
The R script **simulation_unrolling_by_iteration.R** contains the assessment of unrolling by iteration method in Section 5.1
The R script **simulation_gradient_descent.R** contains the assessment of gradient descent method in Section 5.1
The R script **simulation_switch_and_non-switch_case.R** contains the comparison simulation in Section 5.2
The csv file **Indicator.csv** marks the target countries and locations in the data application in Section 6
The R script **data_application.R** contains the data application in Section 6
