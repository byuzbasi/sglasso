# Small end-to-end execution tests only; never a production numerical study.
portable_tests <- function(bundle,out) {
  portable_assert(!dir.exists(out),"Validation output exists; select a new versioned directory")
  dir.create(out,recursive=TRUE)
  before <- portable_inventory(bundle)
  relocated <- file.path(out,"relocated supplement with spaces")
  dir.create(relocated)
  for(f in c(before$file,"MANIFEST.csv","MANIFEST.sha256")) {
    dest <- file.path(relocated,f); dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    portable_assert(file.copy(file.path(bundle,f),dest,overwrite=FALSE),"Relocation copy failed")
  }
  checks <- list(); timings <- list()
  record <- function(name,ok) {
    checks[[name]] <<- data.frame(check=name,passed=isTRUE(ok))
    portable_assert(isTRUE(ok),paste("Test failed:",name))
  }
  invoke <- function(args,label) {
    start <- Sys.time(); path <- file.path(out,paste0(label,".log"))
    args <- c(args,paste0("--data=",paper_inputs$data))
    if(any(grepl("^--study=(auc|mcc)$",args))) args <- c(args,
      paste0("--external=",file.path(out,"runs/portable_external_smoke_v2")))
    # Use the public executable, not an in-session substitute for the launcher.
    result <- suppressWarnings(system2(Sys.which("Rscript"),
      shQuote(c("--vanilla",file.path(relocated,"reproduce.R"),args)),stdout=TRUE,stderr=TRUE))
    writeLines(result,path)
    status <- attr(result,"status")
    error <- if(is.null(status) || status==0L) NULL else paste(tail(result,2L),collapse="; ")
    elapsed <- as.numeric(difftime(Sys.time(),start,units="secs"))
    timings[[label]] <<- data.frame(test=label,seconds=elapsed,error=if(is.null(error)) "" else error)
    cat("SMALL EXECUTION TEST:",label,sprintf("%.2fs",elapsed),if(is.null(error)) "passed" else error,"\n"); flush.console()
    list(result=result,error=error)
  }
  for(f in before$file[endsWith(before$file,".R")]) parse(file.path(relocated,f))
  record("all_R_sources_parse",TRUE)
  record("relocated_manifest",is.null(invoke("--mode=verify","verify")$error))
  runs <- file.path(out,"runs")
  # Configuration generation only: no production fit or resampling call.
  simulation <- paper_configuration(relocated,"simulation")
  e <- portable_simulation(relocated,runs,"production")
  cfg <- e$lsg_configuration_v7("production")
  design <- e$lsg_design_for_stage_v7(file.path(relocated,"research"),"production")
  tasks <- e$lsg_task_grid_v7(design,cfg)
  record("original_400_tasks_seeds_and_settings",identical(cfg,simulation$configuration) &&
    identical(design,simulation$design) && identical(tasks,simulation$tasks))
  external <- portable_external(relocated,runs,"production","portable_external_production_v2")
  original <- paper_configuration(relocated,"external")
  record("original_external_settings_and_folds",identical(external$spec$configuration,original$configuration) &&
    identical(external$spec$fold,original$fold))
  for(study in c("simulation","external","auc","mcc")) {
    prefix <- c("--mode=smoke",paste0("--study=",study),paste0("--output=",runs),"--cores=1","--interval=1")
    record(paste0(study,"_one_unit_pause"),is.null(invoke(c(prefix,"--max-units=1"),paste0(study,"_pause"))$error))
    run <- file.path(runs,paste0("portable_",study,"_smoke_v2"))
    progress <- jsonlite::read_json(file.path(run,"progress.json"),simplifyVector=TRUE)
    record(paste0(study,"_pause_progress"),progress$completed==1L && progress$pending>0L &&
      progress$running==0L && progress$failed==0L && progress$phase=="paused")
    record(paste0(study,"_no_overwrite"),!is.null(invoke(prefix,paste0(study,"_no_overwrite"))$error))
    record(paste0(study,"_resume_complete"),is.null(invoke(c(prefix,"--resume"),paste0(study,"_resume"))$error))
    progress <- jsonlite::read_json(file.path(run,"progress.json"),simplifyVector=TRUE)
    record(paste0(study,"_complete_progress"),progress$completed==progress$total && progress$failed==0L &&
      progress$running==0L && progress$pending==0L && progress$phase=="completed" &&
      file.exists(file.path(run,"progress.tsv")) && !is.null(progress$heartbeat_utc))
    hashes <- portable_inventory(run)
    record(paste0(study,"_verify_run"),is.null(invoke(c("--mode=verify-run",paste0("--study=",study),
      "--stage=smoke",paste0("--output=",runs)),paste0(study,"_verify"))$error))
    record(paste0(study,"_complete_resume_no_refit"),is.null(invoke(c(prefix,"--resume"),paste0(study,"_complete_resume"))$error))
    # Data, shards and final manifests must not change on completed resume.
    stable <- hashes[!grepl("^progress[.]",hashes$file),]
    record(paste0(study,"_complete_outputs_immutable"),all(vapply(seq_len(nrow(stable)),function(i)
      portable_hash(file.path(run,stable$file[i]))==stable$sha256[i],logical(1))))
  }
  # Corrupt a COPY of a completed checkpoint: rejection must not overwrite it.
  invalid <- file.path(out,"invalid")
  src <- file.path(runs,"portable_mcc_smoke_v2")
  dst <- file.path(invalid,"portable_mcc_smoke_v2"); dir.create(dst,recursive=TRUE)
  for(f in portable_inventory(src)$file) {
    p <- file.path(dst,f); dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE); file.copy(file.path(src,f),p)
  }
  shard <- file.path(dst,"shards/chunk_0001.rds")
  bad <- readRDS(shard); bad$signature <- "deliberately_invalid_test_signature"; saveRDS(bad,shard)
  bad_hash <- portable_hash(shard)
  record("invalid_shard_rejected",!is.null(invoke(c("--mode=smoke","--study=mcc","--resume",paste0("--output=",invalid)),"invalid_shard")$error))
  record("invalid_shard_preserved",portable_hash(shard)==bad_hash)
  altered <- file.path(invalid,"portable_mcc_smoke_v2/specification.rds")
  bad_spec <- readRDS(altered); bad_spec$seed_base <- bad_spec$seed_base+1L; saveRDS(bad_spec,altered)
  record("changed_specification_rejected",!is.null(invoke(c("--mode=smoke","--study=mcc","--resume",paste0("--output=",invalid)),"invalid_specification")$error))
  fake <- file.path(out,"wrong_runtime_test_receipt.rds")
  saveRDS(list(accepted=TRUE,checks=data.frame(passed=TRUE),
    distribution_manifest_sha256=portable_hash(file.path(relocated,"MANIFEST.csv")),
    R="deliberately_invalid_R",platform=R.version$platform,
    packages=c(digest=as.character(packageVersion("digest")))),fake)
  guarded <- invoke(c("--mode=production","--study=mcc",paste0("--output=",out,"/not_executed"),
    "--max-units=1",paste0("--validation=",fake)),"production_runtime_guard")
  record("production_wrong_runtime_rejected_before_compute",!is.null(guarded$error) &&
    !dir.exists(file.path(out,"not_executed")))
  record("inputs_unchanged",identical(before,portable_inventory(bundle)) && identical(before,portable_inventory(relocated)))
  tab <- do.call(rbind,checks); timing <- do.call(rbind,timings)
  write.csv(tab,file.path(out,"checks.csv"),row.names=FALSE)
  write.csv(timing,file.path(out,"timings.csv"),row.names=FALSE)
  receipt <- list(schema="portable_supplement_validation_v2",accepted=all(tab$passed),checks=tab,
    distribution_manifest_sha256=portable_hash(file.path(bundle,"MANIFEST.csv")),
    production_executed=FALSE,production_model_fits=0L,relocation_with_spaces=TRUE,
    platform=R.version$platform,R=R.version.string,
    packages=vapply(c("Rcpp","RcppArmadillo","digest","adelie","grpreg","logistf","mltools","jsonlite","knitr","pROC"),
      function(p) as.character(packageVersion(p)),""))
  saveRDS(receipt,file.path(out,"VALIDATED.rds"))
  write.csv(portable_inventory(out),file.path(out,"MANIFEST.csv"),row.names=FALSE)
  cat("LOCAL VALIDATED: small tests only; production NOT executed. Receipt:",file.path(out,"VALIDATED.rds"),"\n")
  invisible(receipt)
}
