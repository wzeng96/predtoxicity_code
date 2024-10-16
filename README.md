# Predicting Chemical Toxicity by Applying a Hierarchical Bayesian Approach with Priors to the Tox21 Assay Data

Repo for the chemical toxicity prediction using Tox21 assay data

## Description

This project uses Tox21 assay to predict chemical toxicity. We are interedted to see if assay-level metadata, such as sex or organisms, has an affect in chemical toxicity, to do this ridge regression, naive bayesian and hierarchical bayesian models are applied for prediction. Feature selection was done before applying the models, coefficients from the models are used for feature importance.
* "feature_selection" folder contains code and results for reducing the features from 208 to 40;
* "predictive_models" folder contains code for training and testing the three machine learning models;
* "feature_importance" folder contains code and results for using the coefficients to generate feature importance plots.

## Getting Started

### Dependencies and Installing

* See "conda_env" for the dependencies and python packages used for this project.

### Executing program

* See the "stepwise" document for the details on running this project.

## Authors

Contributors names and contact info

Wenyu Zeng  
[@wz34](linkedin.com/in/wz96)

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments

Inspiration of readme structure:
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
