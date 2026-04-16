# CCEBA: Enhanced Bat Optimizer for Multi-threshold LN Image Segmentation

This repository provides the MATLAB implementation of CCEBA, the enhanced bat optimizer proposed in the manuscript:

**An Enhanced Bat Optimizer with Elite Selection and Crisscross Strategies for Multi-threshold Image Segmentation of Lupus Nephritis**

CCEBA integrates elite selection and crisscross strategies to improve the exploration-exploitation balance of the original bat optimizer. The code includes the optimizer implementation and scripts used for benchmark optimization and lupus nephritis image segmentation experiments.

## Environment

The experiments in the manuscript were conducted using:

- MATLAB R2021a
- Windows Server 2016
- Intel Xeon Silver 4210R CPU at 2.40 GHz

## Notes

- All stochastic optimizers were executed for 20 independent runs.
- The lupus nephritis images used in the manuscript are not included because of data privacy restrictions.

## Usage

Run the corresponding MATLAB experiment scripts for benchmark optimization or image segmentation after adding the required folders to the MATLAB path.

```matlab
addpath(genpath(pwd));
