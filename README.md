### Correlated-Synthetic-Controls-

## Abstract
Synthetic Control methods have recently gained considerable attention in applications with only one treated unit.  Their popularity is partly based on the key insight that we can predict good synthetic counterfactuals for our treated unit. However, this insight of predicting counterfactuals is generalisable to microeconometric settings where we often observe many treated units. We propose the Correlated Synthetic Controls (CSC) estimator for such situations: intuitively, it creates synthetic controls that are correlated across individuals with similar observables. When treatment assignment is correlated with unobservables, we show that the CSC estimator has more desirable theoretical properties than the difference-in-difference estimator. We also utilise CSC in practice to obtain heterogeneous treatment effects in the well-known Mariel Boatlift example, leveraging additional information from the PSID.

## Contents of the folder
This is the GitHub repository for my MPhil thesis at the University of Oxford. It contains the coding used for the simulation in the paper (Section 5) and the Empirical Application (Section 6).

The main branch contains three folders:
* Simulation - this folder contains all of the coding that was conducted for the Simulation in Section 5 of the paper. 
* PSID - here we can find the coding used for automatically extracting waves and variables from PSID. Shout out to @floswald (Florian Oswald) for package psidR (Check https://github.com/floswald)
* Mariel Boatlift - this is the coding for the results in the Empirical Application of Section 6 of the paper. It uses the data extracted from PSID.

## Future work
In the future, we hope to turn the Correlated Synthetic Controls into an R package

## Disclaimer
Huge thanks to @fditraglia and @maxkasy for excellent supervision! All errors are my own and please contact me when you spot them.
