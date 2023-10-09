# Bayesian belief network for calculating imperilment of fishes based on three R's framework
Repository contains supporting files for manuscript, Dunn et al. 2023. (in press at Ecosphere), "Using resiliency, redundancy, and representation in a Bayesian belief network to assess imperilment of riverine fishes." This manuscript presents a Bayesian belief network based on the "three Rs" framework for calculating imperilment using concepts of resiliency and redundancy within representative settings for Piebald Madtom. The U.S. Fish and Wildlife Service recently adopted the three Rs framework to implement Species Status Assessments for specis petitioned for listing under the U.S. Endangered Species Act. All data and analyses are preliminary and meant to serve as a proof of concept to demonstrate how the model works. 

Description of files:
- piebaldmadtom_inputs_cr.neta: Netica file containing Bayesian belief network that replicates results in Dunn et al. in press at Ecosphere. Netica is a proprietary software that may not be available to manuscript readers, so we developed supplemental R code that replicates Netica's Bayesian belief network:
- bbn_inputs_pmt_cr.xlsx: Excel sheet containing model inputs imported into R by script. Individual sheets contain default prior probabilities, observed values for Piebald Madtom that inform prior probabilities, and conditional probability tables within interior child nodes of Bayesian belief network.
- bbn_piebaldmadtom_cr.R: R script that replicates Bayesian Belief network performed in piebaldmadtom_inputs_cr.neta. Code users need to change the local directory to import data inputs.
- piebaldmadtom_results.xlsx: an output Excel file generated by bbn_piebaldmadtom_cr.R with results matching those presented in Dunn et al. in press at Ecosphere
- dunn_etal-ecosphere_metadata.docx: metadata for data contained in bbn_inputs_pmt_cr.xlsx.
