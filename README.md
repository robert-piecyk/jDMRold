### jDMR: a heuristic DMR caller for population-level WGBS data

We have developed jDMR, a heuristic and fast DMR caller for population level epigenomic studies and also control/treatment experiments.The method uses two main approaches to find DMRs. The first method assumes that differentially methylated regions (DMRs) can only occur as clusters of cytosines. The algorithm starts by identifying cytosine clusters from the raw sequence and uses them as units of analysis to assess whether these regions are differentially methylated between individuals. The second method uses a grid approach by splitting the genome into sliding windows/bins of user defined sizes. Both methods, produces methylation state calls for each cluster or sliding window based on a Hidden Markov Model (HMM). Following which all samples are merged togther in the form of a matrix and non-polymorphic patterns are filtered to obtain DMRs.

### Documentation

1. [jDMR tutorial](docs/jDMR-tutorial.pdf)

### Contributors:

- Rashmi Hazarika - rashmi.hazarika@tum.de
- Yadi Shahryary - y.shahryary@tum.de
- Frank Johannes - frank@johanneslab.org