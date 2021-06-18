### jDMR: a heuristic DMR caller for population-level WGBS data

jDMR is a heuristic and fast DMR caller for population level epigenomic studies and also control/treatment experiments.The method uses two main approaches to find DMRs. The first method assumes that differentially methylated regions (DMRs) can only occur as clusters of cytosines. The algorithm starts by identifying cytosine clusters from the raw sequence and uses them as units of analysis to assess whether these regions are differentially methylated between individuals. The second method uses a grid approach by splitting the genome into sliding windows/bins of user defined sizes. Both methods, produces methylation state calls for each cluster or sliding window based on a Hidden Markov Model (HMM). Following which all samples are merged togther in the form of a matrix and non-polymorphic patterns are filtered to obtain DMRs.


##### Installing from the Github

To install from GitHub (development version), follow the steps given below. 

##### Step 1 — Install a last version of R (>=3.6)

##### Step 2 — In R, please install all dependencies and execute the following commands:
 - 2.1.  install.packages("devtools")
 - 2.2.  library(devtools)
 - 2.3.  devtools::install_github("jlab-code/jDMR")

------------------------------------------------------------------------

### How to use jDMR

Please open the [vignette](https://github.com/jlab-code/jDMR/blob/master/vignettes/jDMR-tutorial.pdf) file.


### Contributors:

- Rashmi Hazarika - rashmi.hazarika@tum.de
- Yadi Shahryary - y.shahryary@tum.de
- Frank Johannes - frank@johanneslab.org