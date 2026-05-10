# Modified-Tamed-Euler-Experiments

One may refer to this project for the codes related to the paper "A modified tamed scheme for stochastic differential equations with superlinear drifts", which can be accessed by https://arxiv.org/abs/2507.09475.

The experiments are done with MATLAB R2025b with seeds fixed.

The MTE_xd.m files give the code of strong and weak convergence experiments, in which tamed Euler, truncated Euler and modified tamed Euler (with random batch method or not) are considered. For 1-d case, we consider the multiplicative noise, while for the 2-d and 20-d, we consider the additive noise. More details can be found in the paper.

The SGLDE_xd.m file gives the code of KL distance experiment for SGLD. Here the KL distance is calculated by definition (i.e., sample and then use kernel density to estimate the distribution).

Remark: 1. the paper related (to the different schemes) are given;

        2. the random batch is to split the coefficients into linear term and superlinear term, and do some unbiased estimation, like E(x-x^3) = 1/2 E(2x)+1/2 E(-2x^3).

## Repository layout

- `programs/`: MATLAB experiment scripts
  - `MTE_1d.m`
  - `MTE_2d.m`
  - `MTE_20d.m`
  - `SGLDE_2d.m`
- `run_program.m`: Root launcher that adds `programs/` to MATLAB path and runs a selected script

## How to run

Option 1 (recommended, from repository root):

```matlab
run_program('MTE_1d')
run_program('MTE_2d')
run_program('MTE_20d')
run_program('SGLDE_2d')
```

Option 2 (run directly from `programs/`):

```matlab
cd programs
MTE_1d
```
