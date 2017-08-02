# CMV Serology project

This repository contains data, scripts/programs, and output of a study into the population dynamics of antibodies against cytomegalovirus 

> Infectious reactivation of cytomegalovirus explaining age- and sex-specific patterns of seroprevalence' by Michiel van Boven, Jan van de Kassteele, Marjolein korndewal, Christiaan van Dorp, Mirjam Kretzschmar, Fiona van der Klis, Hester de Melker, Ann Vossen, and Debbie van Baarle 

See [bioRxiv](http://biorxiv.org/content/early/2017/01/24/102491) for a preprint.

Specifically, the directory contains:

1. The serological data (`data/cmvdata.csv`). Each line is a participant, and columns include age in months (`lftinmnd2`), sex, ethnicity (`nl`), and antibody measurements (last column). *TODO: currently missing... what about cmvdatanew.csv*

2. The contact matrix and demographic composition of the population (`data/contact_intensities_aggregated.csv`). See for notation 

> [van de Kasteele et al.](http://projecteuclid.org/euclid.aoas/1491616883), Efficient estimation of age-specific social contact rates between men and women, Annals of Applied Statistics, Volume 11, Number 1 (2017), 320-339.

3. An R script ('multinomiallogit.R') for the multinomial logistic model. This script still needs to be streamlined (08/12/2016).

4. A stripped-down Mathematica 10.0 program ('estimation.nb') for the analyses with the transmission model.

5. Output of the multinomial logistic regression analyses ('results/p_neg_M.csv' for the posterior of the seroprevalence of the negative component in males, 'results/p_boo.csv' for the the posterior of the seroprevalence of the component of increased antibodies, etcetera). Figures 1-3 in the manuscript are made using this output. See the R script for details.

6. Output of the MCMC output of the transmission model analyses ('results/scenario1.csv' and 'results/scenario2.csv'). See the Mathematica program for details.

Updated 15/06/2017: 

7. An R/Stan model for alternative analyses of the data using the transmission model (`cmvmodel.stan` and `cmvmodel.R`). The code is drafted by Christiaan van Dorp, Sophia de Jong, and Michiel van Boven. This method is much faster than the Mathematica program, and yields identical results.

Updated 20/062017: 

8. csv files (in the `results` folder) containing posterior parameter estimates (5000 samples) for each of models A-G described in the revised manuscript (Table 1). Figures 5-6 and Table 2 are based on these outputs.
