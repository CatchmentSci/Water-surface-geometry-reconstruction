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
