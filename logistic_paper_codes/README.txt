Logistic SGLASSO: Prediction with Correlated Predictor Groups
Study-specific run code, separate from paper_codes/, which accompanies:
Yüzbaşı, B. and Cao, J. (2026). Collinear Groupwise Selection via
Scaled Group Lasso. The American Statistician, 1–23.
Advance online publication, 17 September 2026.
https://doi.org/10.1080/00031305.2026.2709494
Repository location: https://github.com/byuzbasi/sglasso/tree/main/logistic_paper_codes

CONTENTS AND SCOPE
Code only: launchers, necessary unchanged R/Rcpp numerical sources, textual
configurations (including original seeds and folds), manifests and instructions.
No measurements, fitted models, predictions, result archives, manuscript PDFs,
figures, compiled binaries, local validation receipts or deployment paths.
Research code is retained to reproduce the published calculations; substituting
the general package API would not reproduce every study-specific audit.
SOURCE_MAP.csv identifies all unchanged scientific source files.
GPL-3 or later applies to the software; DATA_NOTICE.txt retains original data
and annotation attribution, without distributing their measurements here.

INPUTS
Original measurements: GEO GSE84727 (development) and GSE80417 (external test).
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84727
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80417
Original annotation: https://github.com/waterlandlab/CoRSIV-Methylation-based-SZ-Risk-Score
Data paper: https://doi.org/10.1038/s41398-021-01496-3
Fitting requires the separately distributed, processed CoRSIVSZ_v1.rds file,
not the raw GEO matrix. It is checked against the package's size/SHA-256 metadata
in configuration/CoRSIVSZ_v1.dcf. No hosting endpoint is invented here.
Obtain the exact asset from its separate distribution notice before external
validation or fitting. See sglasso::load_CoRSIVSZ and download_CoRSIVSZ.
No automatic download, imputation, sample filtering or grouping is performed.
Simulation fitting needs no data file. AUC/MCC inference consumes the verified
external run that YOU generate; existing predictions are not bundled.
External refits and resulting intervals may differ across numerical environments;
the original manuscript results are not silently replaced by new runs.

REQUIREMENTS
R, C++17 and an R compilation toolchain; installed digest, jsonlite, knitr,
Rcpp, RcppArmadillo, adelie, grpreg, logistf, mltools and pROC.
No dependencies are installed automatically. Historical versions are recorded
in configuration/software_versions.R. New runs record their own versions.
Numerical execution uses Unix forks (macOS/Linux); Windows is not supported.
Local validation is required in the execution environment before production.

COMMANDS (run from this directory; replace /absolute paths before use)
cd /absolute/path/to/sglasso/logistic_paper_codes
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
Rscript --vanilla reproduce.R --mode=verify
Rscript --vanilla reproduce.R --mode=validate --data=/absolute/CoRSIVSZ_v1.rds --output=/absolute/new_checks

Full numerical workload: run only explicitly in your own terminal after validation.
Rscript --vanilla reproduce.R --mode=production --study=simulation --cores=56 --output=/absolute/new_runs --validation=/absolute/new_checks/VALIDATED.rds
Rscript --vanilla reproduce.R --mode=production --study=external --data=/absolute/CoRSIVSZ_v1.rds --output=/absolute/new_runs --validation=/absolute/new_checks/VALIDATED.rds
Rscript --vanilla reproduce.R --mode=production --study=auc --external=/absolute/new_runs/portable_external_production_v2 --output=/absolute/new_runs --validation=/absolute/new_checks/VALIDATED.rds
Rscript --vanilla reproduce.R --mode=production --study=mcc --external=/absolute/new_runs/portable_external_production_v2 --output=/absolute/new_runs --validation=/absolute/new_checks/VALIDATED.rds

Resume: use the identical command and inputs plus --resume. --max-units=1 allows
an intentional pause; it does not change the scientific grid. Only signature-
matching, numerically validated checkpoints count as completed.
Verify a run using the same --study, --output and input flags, replace
--mode=production with --mode=verify-run --stage=production.
Follow stdout and each run's progress.json/progress.tsv. ETA excludes final
aggregation/verification; SLURM remaining time is reported separately if present.
Output must be outside this code directory; existing outputs are not overwritten.
Do not remove stale locks automatically; inspect them before recovery.

WORKLOAD
Simulation: 400 tasks (eight scenarios, 50 repetitions each), all six audit variants
per task. Approximately three hours was observed with 56 workers, not guaranteed.
External: five development CV folds and refit; approximately 164 minutes observed
sequentially. --cores applies only to simulation tasks, not external CV or inference.
AUC/MCC: 10000 paired, class-stratified draws, fixed trained models. Threshold for
MCC is 0.5. Six intervals and five paired contrasts retain the original adjustment
rules. Small smoke tests use reduced original designs/draw counts, not study evidence.
For the 16-person external smoke fixture only, the MCC count validator uses
its eight cases/eight controls instead of the production 353/322 margins.
The original MCC checker and all numerical criteria are unchanged in production.

The separate code-and-data journal supplement is NOT copied into this folder.
No Git commits, pushes or data releases are performed by these run scripts.
