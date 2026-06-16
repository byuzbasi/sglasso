TRUBA ALL RUNS
==============

Bu klasor, hakem revizyonu icin kullanilacak gercek veri, ana simulasyon,
tuning sensitivity ve scalability kosularini tek yerde toplar.

Calistirma yeri
---------------

TRUBA'da bu klasorun kok dizinine gelin:

  cd /arf/scratch/byuzbasi/sglasso

SLURM log ve cikti klasorlerini olusturun:

  mkdir -p logs results/real_data results/main_simulation results/tuning_sensitivity results/scalability tmp

Ardindan kosulari sirayla gonderin:

  sbatch paper_codes/01_real_data/run_real_slasso_all.sh
  sbatch paper_codes/02_main_simulation/run_sim_slasso_all.sh
  sbatch paper_codes/03_tuning_sensitivity/run_sglasso_tuning_sensitivity.sh
  sbatch paper_codes/04_scalability/run_scalability_analysis.sh

Hepsini arka arkaya gondermek icin:

  bash paper_codes/submit_all_runs.sh

Cikti yapisi
------------

Tum ciktilar tek results klasoru altinda ayrilir:

  results/real_data/
  results/main_simulation/
  results/tuning_sensitivity/
  results/scalability/

Log dosyalari:

  logs/

Kosular
-------

01_real_data:
  Gercek veri uygulamasi. Bardet, Birthwt ve GenAtHum veri setlerini calistirir.
  Varsayilan cikti: results/real_data/

02_main_simulation:
  Ana simulasyon kosusu.
  Icerir:
    sglasso_sim_function_full_tuning.R
    run_sim_slasso.all.R
    run_sim_slasso_all.sh
  Varsayilan cikti: results/main_simulation/

03_tuning_sensitivity:
  Hakem istegi icin SGLASSO tuning kriterleri duyarlilik analizi.
  Varsayilan cikti: results/tuning_sensitivity/

04_scalability:
  Hakem istegi icin sure/olceklenebilirlik analizi.
  Varsayilan cikti: results/scalability/

Opsiyonel ayarlar
-----------------

Parametreler sbatch ile degistirilebilir:

  sbatch --export=ALL,REPEATNUM=50 paper_codes/02_main_simulation/run_sim_slasso_all.sh
  sbatch --export=ALL,NREP=10 paper_codes/03_tuning_sensitivity/run_sglasso_tuning_sensitivity.sh
  sbatch --export=ALL,NREP=10,D_VALUE=0.5 paper_codes/04_scalability/run_scalability_analysis.sh

R modul adi farkliysa:

  sbatch --export=ALL,R_MODULE=apps/R/4.4.0-gcc-12.2.0 paper_codes/04_scalability/run_scalability_analysis.sh

Kutuphane yollari
-----------------

Tum .sh dosyalari varsayilan olarak R 4.3 kutuphanesini kullanir:

  /arf/home/byuzbasi/R/x86_64-pc-linux-gnu-library/4.3

sglasso, simstudy, gglasso, grpreg, grpnet, adelie, data.table, dplyr,
readr, foreach, doParallel, doRNG, ggplot2, caret, Metrics, matrixcalc
ve tibble bu kutuphanede gorundugu icin yeni bir kutuphane dizinine
kurulum yapmadan calisabilir. MASS R modul kutuphanesinden yuklenir.

Gerekirse yol override edilebilir:

  sbatch --export=ALL,R_LIBS_USER=/path/to/R/4.3 paper_codes/02_main_simulation/run_sim_slasso_all.sh

Paket kontrolu
--------------

sglasso paketi GitHub repo kokunden kurulacak sekilde hazirlandi:

  remotes::install_github("byuzbasi/sglasso")

Notlar
-----

- .sh dosyalari paper_codes altinda dursa da varsayilan calisma dizini repo kokudur.
  Bu nedenle loglar ve sonuclar /arf/scratch/byuzbasi/sglasso/logs ve
  /arf/scratch/byuzbasi/sglasso/results altina yazilir.
- OUTDIR verilmezse tum sonuclar yukaridaki results alt klasorlerine yazilir.
- Account, mail ve kullaniciya ozel scratch yolu hard-code edilmedi.



tail -n 160 logs/*5908409*.out
tail -n 160 logs/*5908409*.err
