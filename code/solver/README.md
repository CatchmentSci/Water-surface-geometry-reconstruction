# Solver

| File | Purpose |
|---|---|
| `KLT_synthetic_truth_test_v3.m` | Constructs synthetic water-surface cases from saved inputs, runs or resumes the solver, and records recovery metrics and checkpoints. |
| `KLT_wse_solver_paths_Apoint_block_jacobi_v5.m` | Performs the block-Jacobi water-surface-elevation inversion used by the synthetic-truth workflow. |

`KLT_synthetic_truth_test_v3` expects a case-definition CSV and a solver-input MAT file. Each MAT file must contain `app_in`, `camA_fullmodel`, `camA_first_fullmodel`, `xyzA_wse`, `xyzB_wse`, `aa`, and `wse_map`; `globalPolarity` is optional.

The workflow was tested using MATLAB R2024a with Parallel Computing Toolbox. It requires all five files in `../dependencies`; see that directory's README for their roles.
