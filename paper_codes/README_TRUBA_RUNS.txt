TRUBA ALL RUNS
==============

Bu dosya, hakem revizyonu icin gercek veri, ana simulasyon,
tuning sensitivity ve scalability kosularini TRUBA'da ayni proje
klasoru altinda calistirmak icin hazirlandi.


1) Beklenen server klasor yapisi
--------------------------------

TRUBA'da proje kok dizini:

  /arf/scratch/byuzbasi/sglasso

Beklenen yapi:

  /arf/scratch/byuzbasi/sglasso/
  |-- DESCRIPTION
  |-- NAMESPACE
  |-- README.md
  |-- R/
  |-- src/
  |-- man/
  |-- data/
  |-- paper_codes/
  |-- logs/
  |-- results/
  `-- tmp/

paper_codes klasoru mutlaka kok dizinin icinde bulunmalidir:

  /arf/scratch/byuzbasi/sglasso/paper_codes

paper_codes icindeki beklenen yapi:

  paper_codes/
  |-- 00_setup/
  |-- 01_real_data/
  |-- 02_main_simulation/
  |-- 03_tuning_sensitivity/
  |-- 04_scalability/
  |-- 99_optional_plotting/
  |-- README_TRUBA_RUNS.txt
  |-- run_sim_slasso_all.sh
  |-- run_sim_slasso_debug.sh
  `-- submit_all_runs.sh


2) GitHub'dan temiz kurulum
---------------------------

En temiz yol, eski server klasorunu yedekleyip repo'yu yeniden klonlamaktir:

  cd /arf/scratch/byuzbasi
  mv sglasso sglasso_old_$(date +%Y%m%d_%H%M)
  git clone https://github.com/byuzbasi/sglasso.git sglasso
  cd /arf/scratch/byuzbasi/sglasso

Cikti klasorlerini olusturun:

  mkdir -p logs results/real_data results/main_simulation results/tuning_sensitivity results/scalability tmp


3) Ilk kontrol komutlari
------------------------

Kok dizinde oldugunuzu kontrol edin:

  pwd

Beklenen cevap:

  /arf/scratch/byuzbasi/sglasso

Klasor yapisini kontrol edin:

  ls
  ls paper_codes

Kok dizinde DESCRIPTION dosyasi yoksa veya paper_codes yoksa, server
klasoru eski/yanlis yapidadir.

Son GitHub commit'ini kontrol edin:

  git log --oneline -3

Beklenen son commitlerden biri:

  2833658 Fix paper code run paths


4) R module ve kutuphane kontrolu
---------------------------------

TRUBA'da R 4.3 module yuklenir:

  module purge
  module load apps/R/4.3.0-gcc-11.3.1

Tum .sh dosyalari varsayilan olarak su R kutuphanesini kullanir:

  /arf/home/byuzbasi/R/x86_64-pc-linux-gnu-library/4.3

Manuel kontrol:

  export R_LIBS_USER=/arf/home/byuzbasi/R/x86_64-pc-linux-gnu-library/4.3
  export R_LIBS=/arf/home/byuzbasi/R/x86_64-pc-linux-gnu-library/4.3
  Rscript -e 'cat(.libPaths(), sep="\n")'

Gerekli paketlerin bulundugunu kontrol edin:

  Rscript -e 'pkgs <- c("sglasso","simstudy","gglasso","grpreg","grpnet","adelie","data.table","dplyr","readr","foreach","doParallel","doRNG","ggplot2","caret","Metrics","matrixcalc","tibble"); print(sapply(pkgs, find.package))'

sglasso GitHub kokunden kurulacak sekilde hazirlandi:

  remotes::install_github("byuzbasi/sglasso")

Paket versiyon kontrolu:

  Rscript -e 'library(sglasso); packageVersion("sglasso"); find.package("sglasso")'

Beklenen versiyon:

  1.1.12


5) Run script kontrolu
----------------------

Scriptlerin parse/syntax kontrolu:

  bash -n paper_codes/01_real_data/run_real_slasso_all.sh
  bash -n paper_codes/02_main_simulation/run_sim_slasso_all.sh
  bash -n paper_codes/03_tuning_sensitivity/run_sglasso_tuning_sensitivity.sh
  bash -n paper_codes/04_scalability/run_scalability_analysis.sh
  bash -n paper_codes/submit_all_runs.sh

R scriptlerinin dosya olarak bulundugunu kontrol edin:

  test -f paper_codes/01_real_data/run_real_slasso_all.R
  test -f paper_codes/02_main_simulation/run_sim_slasso.all.R
  test -f paper_codes/03_tuning_sensitivity/run_sglasso_tuning_sensitivity.R
  test -f paper_codes/04_scalability/run_scalability_analysis.R


6) Kosulari gonderme
--------------------

Kok dizinden calisin:

  cd /arf/scratch/byuzbasi/sglasso

Tek tek gondermek icin:

  sbatch paper_codes/01_real_data/run_real_slasso_all.sh
  sbatch paper_codes/02_main_simulation/run_sim_slasso_all.sh
  sbatch paper_codes/03_tuning_sensitivity/run_sglasso_tuning_sensitivity.sh
  sbatch paper_codes/04_scalability/run_scalability_analysis.sh

Hepsini arka arkaya gondermek icin:

  bash paper_codes/submit_all_runs.sh

Onemli:
  .sh dosyalari paper_codes altinda dursa da varsayilan calisma dizini
  repo kokudur. Bu nedenle loglar ve sonuclar kokteki logs/ ve results/
  altina yazilir.


7) Opsiyonel ayarlar
--------------------

Ana simulasyon tekrar sayisi:

  sbatch --export=ALL,REPEATNUM=50 paper_codes/02_main_simulation/run_sim_slasso_all.sh

Tuning sensitivity tekrar sayisi:

  sbatch --export=ALL,NREP=10 paper_codes/03_tuning_sensitivity/run_sglasso_tuning_sensitivity.sh

Scalability d degeri:

  sbatch --export=ALL,NREP=10,D_VALUE=0.5 paper_codes/04_scalability/run_scalability_analysis.sh

R module override:

  sbatch --export=ALL,R_MODULE=apps/R/4.4.0-gcc-12.2.0 paper_codes/04_scalability/run_scalability_analysis.sh

R kutuphane yolu override:

  sbatch --export=ALL,R_LIBS_USER=/arf/home/byuzbasi/R/x86_64-pc-linux-gnu-library/4.3 paper_codes/02_main_simulation/run_sim_slasso_all.sh


8) Cikti yapisi
---------------

Tum ciktilar kokteki results klasoru altinda ayrilir:

  results/real_data/
  results/main_simulation/
  results/tuning_sensitivity/
  results/scalability/

Log dosyalari:

  logs/

tmp dosyalari:

  tmp/

Beklenen onemli cikti dosyalari:

  results/real_data/ALL_summary.csv
  results/main_simulation/tables/paper_summary_overall.csv
  results/main_simulation/tables/all_sim_slasso_results.csv
  results/tuning_sensitivity/sglasso_tuning_summary_by_scenario.csv
  results/tuning_sensitivity/summary_tune_scenario_table.csv
  results/scalability/scalability_summary.csv


9) Job izleme
-------------

Aktif joblari gormek:

  squeue -u byuzbasi

Job detayini gormek:

  scontrol show job JOBID

Loglari canli izlemek:

  tail -f logs/*JOBID*.out
  tail -f logs/*JOBID*.err

Son 160 satiri okumak:

  tail -n 160 logs/*JOBID*.out
  tail -n 160 logs/*JOBID*.err

Ornek:

  tail -n 160 logs/*5908409*.out
  tail -n 160 logs/*5908409*.err

Tamamlanan job durumunu gormek:

  sacct -j JOBID --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS

Butun kendi joblarini gormek:

  sacct -u byuzbasi --starttime today --format=JobID,JobName,State,ExitCode,Elapsed

Job iptal etmek:

  scancel JOBID


10) Hata ayiklama
-----------------

Hata:

  Cannot find run_scalability_analysis.R.

Neden:
  Eski klasor yapisindan calisiliyor olabilir. Yeni yapida script
  paper_codes/04_scalability/ altindadir.

Kontrol:

  pwd
  ls
  ls paper_codes/04_scalability

Dogru submit:

  sbatch paper_codes/04_scalability/run_scalability_analysis.sh


Hata:

  Cannot find sglasso_sim_function_full_tuning.R.

Kontrol:

  ls paper_codes/02_main_simulation/sglasso_sim_function_full_tuning.R
  ls paper_codes/03_tuning_sensitivity/sglasso_sim_function_full_tuning.R
  ls paper_codes/04_scalability/sglasso_sim_function_full_tuning.R


Hata:

  Package not installed: PACKAGE_NAME

Kontrol:

  Rscript -e 'find.package("PACKAGE_NAME")'
  Rscript -e 'cat(.libPaths(), sep="\n")'

R_LIBS_USER dogru degilse:

  export R_LIBS_USER=/arf/home/byuzbasi/R/x86_64-pc-linux-gnu-library/4.3
  export R_LIBS=/arf/home/byuzbasi/R/x86_64-pc-linux-gnu-library/4.3


Hata:

  there is no package called 'sglasso'

Kontrol:

  Rscript -e 'find.package("sglasso")'

Gerekirse yeniden kurulum:

  Rscript -e 'remotes::install_github("byuzbasi/sglasso", upgrade="never", force=TRUE)'


11) Kisa test kosulari
----------------------

Scalability kisa test:

  sbatch --export=ALL,P_GRID=1000,J_GRID=200,NREP=1,NLAMBDA=3,OUTDIR=results/scalability_test paper_codes/04_scalability/run_scalability_analysis.sh

Tuning sensitivity kisa test:

  sbatch --export=ALL,NREP=1,CORES=1,OUTDIR=results/tuning_sensitivity_test paper_codes/03_tuning_sensitivity/run_sglasso_tuning_sensitivity.sh

Ana simulasyon kisa test:

  sbatch --export=ALL,REPEATNUM=1,OUTDIR=results/main_simulation_test paper_codes/02_main_simulation/run_sim_slasso_all.sh

Gercek veri kisa test:

  sbatch --export=ALL,NREP=1,OUTDIR=results/real_data_test paper_codes/01_real_data/run_real_slasso_all.sh


12) Notlar
----------

- Account, mail ve kullaniciya ozel scratch yolu hard-code edilmedi.
- BLAS/OpenMP thread sayilari .sh dosyalarinda 1'e sabitlendi.
- Gercek veri ve tuning sensitivity kosularinda varsayilan cekirdek
  sayisi SLURM_CPUS_PER_TASK degerinden turetilir.
- OUTDIR verilmezse tum sonuclar results/ altindaki ilgili klasore yazilir.
- Server klasoru eski duzendeyse en temiz cozum repo'yu yeniden klonlamaktir.
