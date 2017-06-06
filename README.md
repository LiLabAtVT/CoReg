# CoReg
## Installation
### 1. Install CoReg package through github
To install CoReg through github, first you will need to install 'devtools' package for R:
```R
install.packages("devtools")
```
Load the devtools package:
```R
library(devtools)
```
Install CoReg package through github:
```R
install_github("LiLabAtVT/CoReg")
```
## Quick start
This tutorial gives an example of identifying co-regulatory modules and performing rewiring simulation in Arabidopsis transcription network using CoReg package
### 1. Load package and network
Load CoReg package in R:
```R
library(CoReg)
```
Load the arabidopsis example network(use `loadNetwork()` function if you want to load you own network for analysis):
```R
athNet<-data(athNet)
```
Below is an example of using `loadNetwork()` to load other network. Note that network file should be formatted as two-column edge list. First column represents the ID of transcription factor and second column represents the ID of target gene.
```R
athNet<-loadNetwork("araNet.csv",",")
```
### 2. Identify co-regulatory modules
Identify co-regulatory modules in arabidopsis network:
```R
athRes<-CoReg(athNet)
```
Show module finding result
```R
athRes$module
```
Show similarity matrix
```R
athRes$similarity_matrix
```
### 3. Perform simulation analysis
#### a. Rewiring recall score
`networkSim()` function performs rewiring simulation on Arabidopsis network and get rewiring recall score. CoReg will be compared to other three module-finding methods: label propagation, edge betweeness and walk trap. 50 nodes will be duplicated and rewiring probabilities will be 0.3 and 0.5. 
```R
simRes<-networkSim(athNet,nDup = 50, dDup = 10, c(0.3,0.5),c("lp","wt","eb"),2)
```
Show simulation result
```R
print(simRes)
simRes$evalResult
```
Plot the rewiring recall score curves
```R
plot(simRes)
```
#### b. auROC analysis
`computeAuROC()` function computes simulation-based AUC value and ROC curves. The network simulation process is the same as `networkSim()`. The AUC and ROC for CoReg + each of three similarity indices are computed:
```R
auROCres <- computeAuROC(athNet,nDup=50,dDup=10,rewProb=0.5,simMethods=c("jaccard","geometric","invlogweighted","wt"))
```
Show AUC values
```R
re$AUC
```
Plot ROC curves
```R
plot(auROCres)
```
