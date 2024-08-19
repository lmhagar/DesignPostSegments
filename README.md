# DesignPostSegments
Supplementary code for Design of Posterior Analyses with Sampling Distribution Segments manuscript

The code files for this manuscript are divided into three groups. The code files for each group should be run in order.

Group 1: code for the numerical studies with the gamma example in Sections 5.1 and 6
- 01-sec5p1-setup.R: obtains the data generation processes and informative priors via simulated data
- 02-sec5p1-functions.R: defines functions used to conduct the numerical studies with this example
- 03-sec5p1-studies.R: conducts the numerical studies and generates the figures and numbers in the text for these sections
- JAGS_gamma.txt: model to approximate the posteriors for the gamma model using MCMC

Group 2: code for the numerical studies with the regression example in Sections 5.2 and Appendix D
- 04-sec5p2-functions.R: defines functions used to conduct the numerical studies with this example
- 05-sec5p3-studies.R: conducts the numerical studies and generates the figures and numbers in the text for these sections

Group 3: code for the numerical studies with the Weibull example in Appendix E
- 06-appe-setup.R: obtains the data generation processes and informative priors via simulated data
- 07-appe-functions.R: defines functions used to conduct the numerical studies with this example
- 08-appe-studies.R: conducts the numerical studies and generates the figures and numbers in the text for this appendix
- JAGS_weibull.txt: model to approximate the posteriors for the Weibull model using MCMC
