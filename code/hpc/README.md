# HPC workflows

These Bash scripts submit and manage the synthetic KLT experiments under Slurm using MATLAB R2024a.

| File | Purpose |
|---|---|
| `submit_klt_syn.sh` | Submits the configured set of synthetic cases as individual Slurm jobs. |
| `resubmit_missing_klt_jobs.sh` | Checks whether named cases are active and resubmits missing jobs. |
| `run_one_klt_case.sh` | Runs one case selected by the exported `KLT_S_ID` value. |
| `run_one_klt_case_ramlog.sh` | Runs one case while recording MATLAB and descendant-process memory use. |

The `real` subdirectory contains the operational real-video workflow:

| File | Purpose |
|---|---|
| `real/submit_klt_analysis_batch.sh` | Resolves the configured video cases and submits one Slurm job per case. |
| `real/resubmit_missing_klt_analysis_jobs.sh` | Checks and resubmits missing real-video jobs. |
| `real/run_one_klt_analysis_case.sh` | Loads a video solver input/checkpoint and runs the block-Jacobi solver. |
| `real/klt_analysis_case_lookup.tsv` | Maps video identifiers to sweep/discharge values for submission. |

The scripts reflect the Newcastle University Comet configuration and contain user- and system-specific absolute paths. Before running them elsewhere, update the Slurm account and partition, working directories, log locations, project root, and MATLAB module name.

The case CSV and solver-input MAT files are part of the associated Zenodo dataset rather than this Git repository. The four top-level scripts are byte-for-byte copies of the operational files under `/mnt/nfs/home/nmp65/slurm_commands/synthetic` on Comet; the three scripts and lookup under `real` match `/mnt/nfs/home/nmp65/slurm_commands/real`.

`run_one_klt_case_ramlog.sh`, selected by the operational `submit_klt_syn.sh`, uses Excel-visible row numbers from `synthetic_cases.csv` (row 1 is the header) and maps them as follows:

| Excel rows | Solver input |
|---|---|
| 2-10 | `solver_inputs_0pt88.mat` |
| 11-19 | `solver_inputs_1pt50.mat` |
| 20-31 | `solver_inputs_1pt70.mat` |
| 32-40 | `solver_inputs_1pt10.mat` |

The active `S_IDS` value in `submit_klt_syn.sh` records the selected operational run and does not define the complete dataset. `run_one_klt_case.sh` and `resubmit_missing_klt_jobs.sh` are also retained exactly from Comet; the standard worker contains an earlier row-range configuration and should not be substituted for the RAM-logging worker when reproducing the final 39-case mapping.

The real workflow additionally requires `sweep_limits.csv`, the selected `<identifier>_solver_inputs.mat` files, and any checkpoints used for resumption. These data files are in the associated Zenodo archive. The lookup table is included here because the unmodified Comet submission scripts expect it beside the scripts; the Zenodo copy is identical.

Comet's `synthetic/testing` directory contains six MATLAB/Slurm environment diagnostics. They are intentionally excluded because neither output-generation workflow calls them.
