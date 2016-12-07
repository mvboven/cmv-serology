This directory contains data, scripts/programs, and output of a study into the population dynamics of antibodies against cytomegalovirus ('Infectious reactivation of cytomegalovirus explaining age- and sex-specific patterns of seroprevalence' by Michiel van Boven, Jan van de Kassteele, Marjolein korndewal, Christiaan van Dorp, Mirjam Kretzschmar, Fiona van der Klis, Hester de Melker, Ann Vossen, and Debbie van Baarle).

Specifically, the directory contains:

1-The serological data ('data.csv'). Each line is a participant, and columns include age in months ('lftinmnd2'), sex, ethnicity ('nl'), and antibody measurements (last column).

2-The contact matrix and demographic composition of the population ('contact intensities aggregated.csv'). See van de Kasteele et al (2016) Efficient estimation of age-specific social contact rates between men and women, accepted for publication in the Annals of Applied Statistics.

3-An R script ('xxx.r') for the multinomial logistic model. mention version R and script.

4-A stripped-down Mathematica 10.0 program ('estimation (07122017).nb') for the analyses with the transmission model.  

5-Output of the multinomial logistic regression analyses ('p_neg_M.csv' for the posterior of the seroprevalence of the negative component in males, 'p_boo.csv' for the the posterior of the seroprevalence of the component of increased antibodies, etcetera). Figures 1-3 in the manuscript are made using this output. See the R script for details.

6-Output of the MCMC output of the transmission model analyses ('scenario1.csv' and 'scenario2.csv'). See the Mathemtica program details.

7-A directory with a Python/STAN model for alternative analyses of the data using the transmission model. The code is musch faster than the Mathematica code and yields identical results, but is formally an approximation to the exact analyses in the Mathematica program. [TO BE COMPLETED]
