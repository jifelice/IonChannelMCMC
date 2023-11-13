# IonChannelMCMC

This repository contains the R code to explore the parameter space of an ion channel kinetic model.

## Contents
The following source codes are necessary for performing the exploration and must be placed in a source directory that is defined in the second line of the code: AJUSTEcall.R and AJUSTEcallTT.R with the variable: "directorioFuentes". To perform stage 1/stage 2 the codes AJUSTEcall.R/AJUSTEcallTT.R must be executed. They call other subprograms that are also listed below.

_**AJUSTEcall.R**_

Runs stage 1. Arguments are read from the file _**Parametros.R**_ in /additional\_files directory. 

_**AJUSTEcallTT.R**_

Runs stage 2. Arguments are read from the file _**ParametrosTT.R**_ in /additional\_files directory.

_**AJUSTEcoreDIAG.R**_

Calculates the Po vs t curves.

_**AJUSTEcoreMCMC.R**_

Performs _N<sub>c</sub>_ cycles of Monte Carlo steps at stage 1.

_**AJUSTEcoreMCMCTT.R**_

Performs _N<sub>c</sub>_ cycles of Monte Carlo steps at stage 2.

_**AJUSTEdensfunction.R**_

Calculates $d^2$

_**AJUSTEmatrizHVCN1.R**_

Calculates de W matrix for the model considered as an example.

_**AJUSTEinforme.R**_ and _**AJUSTEinformeTT.R**_ save to .txt files the initial values used in the simulations for stage1 and stage2, respectively. 
