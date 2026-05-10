# Modified-Tamed-Euler-Experiments

One may refer to this project for the codes related to the paper "A modified tamed scheme for stochastic differential equations with superlinear drifts", which may be accessed by https://arxiv.org/abs/2507.09475.

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
