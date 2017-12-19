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
### 3. Perform rewiring simulation analysis on Arabidopsis network
#### a. Rewiring recall score
`rewSim()` function performs rewiring simulation on Arabidopsis network and get rewiring recall score. CoReg will be compared to other three module-finding methods: label propagation, edge betweeness and walk trap. 50 nodes will be duplicated and rewiring probabilities will be 0.3 and 0.5. 
```R
simRes<-rewSim(athNet,nDup = 50, dDup = 10, c(0.3,0.5),c("lp","wt","eb"),2)
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
### 4. Perform evaluation of different module-finding methods on simulated network
#### a. Compute Normalized Mutual Information (NMI) score
`netSimAndEval` calls function `generateSimNet()` to generate simulated network(s) with pre-specified modular structure, and then different module-finding algorithms are run to identify modules. The correlation between pre-specified modules and algorithm identified modules is calculated using NMI score. The following example generates a simulated network with 5 pre-specified modules. Each module has 10 regulators and each regulator has 20 targets. There are 100 other nodes in the network which do not have outgoing edges (auxiliary nodes). The co-regulation probability for the simulated network is 0.5. After network is constructed, four clustering methods (as specified by `testMethods`) will be run to identify modules in the simulated network. 
```R
re<-netSimAndEval(10,5,20,100,0.5,testMethods=c("coregJac","lp","wt","eb"))
```
See the summary of simulated network(s)
```R
summary(re)
```
Plot the evaluation result
```R
plot(re)
```
#### b. Generate a simulated network
CoReg also provides the functionality of generating the simulated network for other use. The simulated network is returned as an edge list, represented by a two column matrix in R. In the example below, we use the same paramters as we have in the previous example to generate a simulated network
```R
re<-generateSimNet(10,5,20,100,0.5)

# This is the edge list
re$el

# This is the ground-truth module partition
re$modulePartition
```
## Citation
Please cite the following paper if CoReg is used in a publication:\
Song Q, Grene R, Heath LS, Li S: Identification of regulatory modules in genome scale transcription regulatory networks. BMC Syst Biol 2017, 11:140.\
https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-017-0493-2
