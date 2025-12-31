# Benthic-assemblages modelling

Joint species distribution modelling of Antarctic benthic assemblages.

This repository contains the model implementation used to analyse temporal and depth-related patterns in benthic community assembly in Potter Cove (Antarctica). The analysis is based on a Joint Species Distribution Model (JSDM) implemented under the Hierarchical Modelling of Species Communities (HMSC) framework.

The model integrates species occurrence and abundance (percent cover) data with environmental variables, species functional traits, and phylogenetic relationships to infer the ecological filtering processes shaping benthic assemblages through time.


# File description

The repository is organized in a modular way to facilitate reproducibility of the analyses and results presented in the manuscript. Below is a brief description of the main files included in the Model/ folder:

 - stan_model.stan: contains the Stan code specifying the JSDM used in our study, including the hierarchical structure, environmental effects, species functional traits, and phylogenetic relationships.

- run_model.R: R script used to prepare the data, compile the Stan model, and run the Bayesian inference within the HMSC framework. This script controls the model settings, sampling parameters, and storage of model outputs.

- model_prediction.R: contains all scripts required to generate model predictions, post-processing analyses, and the figures presented in the manuscript. *Running this script allows full replication of the results and visualizations based on the model outputs*.

Each script starts with a description of the required input files and includes detailed comments guiding the user through the workflow. For a complete explanation of the modelling framework, assumptions, and prediction procedures, we strongly recommend consulting the Extended Methods section of the Supplementary Information.

Benthic abundance data belong to Laboratorio de Ecosistemas Marinos y Polares (ECOMARES; Instituto de Diversidad y Ecología Animal, Universidad Nacional de Córdoba - Consejo Nacional de Ciencia y Tecnología).
*Requests to use the data should be addressed to Dr. Ricardo Sahade (rsahade@unc.edu.ar)*. Trait data and taxonomic matrix are available in Supplementary Information A and D, and should be cited accordingly.
