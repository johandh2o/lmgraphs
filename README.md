# Simulation Study: Causal Estimation under Missing Data Mechanisms

This repository contains the R code and supporting materials for a simulation study comparing multiple causal effect estimators —FATE, NATE, Multiple Imputation, and the Missing Indicator Method (MIM)— under various missing data scenarios.

## 📄 Description

The simulation mimics a causal system with:
- Binary exposure `A`
- Outcome `Y`
- Confounder `X` subject to missingness
- Four missingness mechanisms, generating different missing data rates

## 🔍 Estimators Evaluated

- **FATE (Full Average Treatment Effect)** using AIPW and pseudo-outcomes  
- **NATE (Natural Average Treatment Effect)** integrating across missingness strata  
- **Multiple Imputation (MI)** using `mice::pmm`  
- **Missing Indicator Method (MIM)** with `X=0` substitution  

## 📁 Files

- `UAI2.R`: Main script to reproduce the simulations and generate the final plot
- `/results`: For saving generated plots
- `/data`: For storing simulated data (not included in current script)
