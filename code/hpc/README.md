# HPC workflows

These Bash scripts submit and manage the synthetic KLT experiments under Slurm using MATLAB R2024a.

| File | Purpose |
|---|---|
| `submit_klt_syn.sh` | Submits the configured set of synthetic cases as individual Slurm jobs. |
| `resubmit_missing_klt_jobs.sh` | Checks whether named cases are active and resubmits missing jobs. |
| `run_one_klt_case.sh` | Runs one case selected by the exported `KLT_S_ID` value. |
| `run_one_klt_case_ramlog.sh` | Runs one case while recording MATLAB and descendant-process memory use. |

The scripts reflect the Newcastle University Comet configuration and contain user- and system-specific absolute paths. Before running them elsewhere, update the Slurm account and partition, working directories, log locations, project root, and MATLAB module name.

The case CSV and solver-input MAT files are part of the associated Zenodo dataset rather than this Git repository.

Both worker scripts use Excel-visible row numbers from `synthetic_cases.csv` (row 1 is the header) and map them to the four current solver inputs as follows:

| Excel rows | Solver input |
|---|---|
| 2-10 | `solver_inputs_0pt88.mat` |
| 11-19 | `solver_inputs_1pt50.mat` |
| 20-31 | `solver_inputs_1pt70.mat` |
| 32-40 | `solver_inputs_1pt10.mat` |

The `S_IDS` arrays in the submission scripts are deliberately editable job selections; they do not define the complete 39-case dataset.
