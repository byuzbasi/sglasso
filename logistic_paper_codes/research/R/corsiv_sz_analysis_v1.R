# Analysis-only adapter: frozen fitting/selection functions are inherited unchanged.
csa_context <- function(root) {
  b <- new.env(parent = globalenv())
  files <- c("R/genathum_binary_workflow_v1.R", "R/genathum_binary_io_v1.R",
    "R/genathum_binary_firth_v3.R", paste0("R/all_binary_fair_tuning_v", 2:4, ".R"),
    unlist(lapply(c(1,2,3,4,6,7), function(v) paste0("R/all_binary_", c("workflow","io"), "_v", v, ".R"))),
    "R/all_binary_workflow_v8.R", "R/gse25066_binary_workflow_v1.R")
  for (f in files) sys.source(file.path(root, f), b)
  for (name in c("source_inventory", "atomic_rds", "valid_shard", "finalize", "verify"))
    assign(paste0("base_", name), get(paste0("gse_", name), b), b)
  bindings <- c(gse_load_data="csa_data", gse_source_inventory="csa_inventory",
    gse_specification="csa_spec", gse_atomic_rds="csa_atomic", gse_valid_shard="csa_valid",
    gse_progress="csa_progress", gse_finalize="csa_finalize", gse_verify="csa_verify")
  for (name in names(bindings)) {
    f <- get(bindings[[name]], mode="function"); environment(f) <- b
    assign(name, f, b)
  }
  b
}

csa_data <- function(root, stage) {
  out <- file.path(root, "data_processed/corsiv_sz_complete_groups_v1")
  m <- readRDS(file.path(out,"MANIFEST.rds"))
  for (i in seq_len(nrow(m))) {
    f <- file.path(out,m$file[i])
    gse_assert(file.exists(f) && file.info(f)$size == m$bytes[i] &&
      digest::digest(file=f,algo="sha256") == m$sha256[i], "Prepared input checksum failure")
  }
  path <- file.path(out,"data.rds"); x <- readRDS(path)
  gse_assert(identical(x$schema_version,"corsiv_sz_complete_groups_v1") &&
    identical(dim(x$train$X),c(847L,1107L)) && identical(dim(x$test$X),c(675L,1107L)) &&
    length(x$group_name)==409L && sum(x$train$y)==414L && sum(x$test$y)==353L &&
    identical(colnames(x$train$X), x$probe_id) && identical(colnames(x$test$X),x$probe_id) &&
    !anyDuplicated(c(x$train$sample_id,x$test$sample_id)), "CoRSIV identity changed")
  if (stage == "smoke") {
    # Small execution test, not evidence about external-cohort performance.
    sizes <- tabulate(x$group)
    chosen <- vapply(c(2L,3L,4L), function(s) which(sizes==s)[1L], integer(1))
    chosen <- sort(chosen); cols <- which(x$group %in% chosen)
    for (part in c("train","test")) {
      z <- x[[part]]; n <- if (part=="train") 24L else 8L
      rows <- unlist(lapply(0:1,function(y) head(which(z$y==y),n)))
      for (key in c("y","sample_id","geo_sample_id","source")) z[[key]] <- z[[key]][rows]
      z$metadata <- z$metadata[rows,,drop=FALSE]; z$X <- z$X[rows,cols,drop=FALSE]
      x[[part]] <- z
    }
    x$group <- match(x$group[cols], chosen); x$group_name <- x$group_name[chosen]
    x$probe_id <- x$probe_id[cols]
  }
  x$processed_sha256 <- digest::digest(file=path,algo="sha256")
  x
}

csa_inventory <- function(root) {
  inherited <- base_source_inventory(root)
  own <- c("R/corsiv_sz_analysis_v1.R", "scripts/112_run_corsiv_sz_v1.R",
    "tests/test_corsiv_sz_analysis_v1.R", "CORSIV_SZ_ANALYSIS_PROTOCOL_V1.md")
  paths <- file.path(root,own)
  rbind(inherited,data.frame(file=own,bytes=unname(file.info(paths)$size),
    sha256=unname(vapply(paths,digest::digest,"",file=TRUE,algo="sha256"))))
}

csa_spec <- function(root,e,stage,version) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  data <- gse_load_data(root,stage); cfg <- allb_configuration_v8(e,stage)
  cfg$schema_version <- "corsiv_sz_analysis_configuration_v1"
  cfg$endpoint_id <- "schizophrenia_vs_control"; cfg$event_label <- "schizophrenia"
  cfg$reference_label <- "control"; cfg$outer_train_fraction <- NULL
  cfg$repeats <- 1L; cfg$inner_folds <- if(stage=="smoke") 2L else 5L
  cfg$seed_base <- if(stage=="smoke") 1010500L else 1010501L
  cfg$response_threshold <- .5
  identity <- list(schema_version="corsiv_sz_analysis_v1",stage=stage,version=version,
    configuration=cfg,fold=gse_fold_assignment(data$train$y,cfg$inner_folds,cfg$seed_base),
    data_identity=list(processed_sha256=data$processed_sha256,
      train_sample_id=data$train$sample_id,test_sample_id=data$test$sample_id,
      train_response_sha256=digest::digest(data$train$y,algo="sha256"),
      test_response_sha256=digest::digest(data$test$y,algo="sha256"),
      group_sha256=digest::digest(data$group,algo="sha256")),
    source_inventory=gse_source_inventory(root),runtime=list(R=R.version.string,
      platform=R.version$platform, rng_kind=RNGkind(), packages=vapply(c("Rcpp","RcppArmadillo","digest",
        "adelie","grpreg","logistf","mltools","jsonlite"),function(p) as.character(packageVersion(p)),"")))
  identity$scientific_signature <- digest::digest(identity,algo="sha256")
  identity
}

csa_atomic <- function(object,path) {
  base_atomic_rds(object,path)
  if (basename(dirname(path))=="shards") base_atomic_rds(list(bytes=file.info(path)$size,
    sha256=digest::digest(file=path,algo="sha256")),paste0(path,".receipt.rds"))
  invisible(path)
}

csa_valid <- function(path,spec,unit) {
  if(!file.exists(path) || !file.exists(paste0(path,".receipt.rds"))) return(FALSE)
  tryCatch({
    r <- readRDS(paste0(path,".receipt.rds"))
    identical(file.info(path)$size,r$bytes) &&
      identical(digest::digest(file=path,algo="sha256"),r$sha256) && base_valid_shard(path,spec,unit)
  },error=function(e) FALSE)
}

csa_progress <- function(out,spec,status,phase,started,current=NULL) {
  units <- names(status)
  valid <- vapply(units,function(u) gse_valid_shard(gse_shard(out,u),spec,u),logical(1))
  fresh <- units[valid & status=="passed" & startsWith(units,"fold_")]
  duration <- vapply(fresh,function(u) readRDS(gse_shard(out,u))$runtime_seconds,0.0)
  pending_folds <- sum(!valid & startsWith(units,"fold_"))
  eta <- if(length(duration)>=2L && pending_folds>0L) mean(duration)*pending_folds else NULL
  elapsed <- as.numeric(difftime(Sys.time(),started,units="secs")); now<-Sys.time()
  running <- as.integer(!is.null(current) && !valid[[current]])
  end <- suppressWarnings(as.numeric(Sys.getenv("SLURM_JOB_END_TIME","")))
  p <- list(run_id=spec$version,slurm_job_id=Sys.getenv("SLURM_JOB_ID","local"),phase=phase,
    work_unit="validated_cv_fold_or_full_development_refit",total=length(units),completed=sum(valid),
    running=running,failed=sum(status=="failed"),pending=length(units)-sum(valid)-running-sum(status=="failed"),
    percent_completed=100*mean(valid),elapsed_seconds=elapsed,
    throughput_units_per_second=if(elapsed>0) sum(valid & status=="passed")/elapsed else 0,
    eta_seconds=eta,estimated_completion_utc=if(is.null(eta)) NULL else format(now+eta,tz="UTC",usetz=TRUE),
    eta_status=if(is.null(eta)) "estimating_or_unknown" else "fresh_invocation_fold_durations",
    eta_scope="remaining_cv_folds_including_current_excludes_refit_finalization_compilation",
    slurm_remaining_seconds=if(is.finite(end)) max(0,end-as.numeric(now)) else NULL,
    heartbeat_utc=format(now,tz="UTC",usetz=TRUE),current_work_unit=current,
    most_recently_completed_work_unit=if(any(valid)) tail(units[valid],1L) else NULL)
  tmp<-tempfile("progress_",out); jsonlite::write_json(p,tmp,auto_unbox=TRUE,null="null",pretty=TRUE)
  gse_assert(file.rename(tmp,file.path(out,"progress.json")),"Progress rename failed")
  h<-file.path(out,"progress.tsv")
  row<-as.data.frame(p[c("heartbeat_utc","phase","total","completed","running","failed","pending",
    "percent_completed","elapsed_seconds","throughput_units_per_second")])
  row$eta_seconds<-if(is.null(eta)) NA_real_ else eta
  write.table(row,h,sep="\t",quote=FALSE,row.names=FALSE,col.names=!file.exists(h),append=file.exists(h))
  cat(p$heartbeat_utc,phase,sum(valid),"/",length(units),"validated; current:",current,
    "; ETA:",if(is.null(eta)) "unknown" else round(eta),"\n"); flush.console()
  invisible(p)
}

csa_finalize <- function(e,data,spec,out,selected) {
  if(!file.exists(file.path(out,"COMPLETE.rds"))) base_finalize(e,data,spec,out,selected)
  rows <- readRDS(file.path(out,"final/selected_refit_rows.rds"))
  grid <- spec$configuration$fit_configuration$fair_common_lambda_relative_grid
  range <- data.frame(method=rows$method,point_type=rows$point_type,
    lambda_relative=rows$lambda_relative_to_reference,
    common_lower_selected=rows$point_type=="finite" & abs(rows$lambda_relative_to_reference-tail(grid,1))<1e-10,
    policy="Frozen finite-grid endpoint allowed; beyond-grid performance unassessed")
  gse_atomic_rds_or_identical(range,file.path(out,"final/lambda_range_report.rds"))
  response <- readRDS(file.path(out,"final/results.rds"))
  response$accuracy <- 1-response$classification_error
  gse_atomic_rds_or_identical(response,file.path(out,"final/response_summary.rds"))
  paths<-c(file.path(out,"specification.rds"),list.files(file.path(out,"shards"),full.names=TRUE),
    list.files(file.path(out,"final"),full.names=TRUE),file.path(out,"COMPLETE.rds"))
  m<-data.frame(file=substring(paths,nchar(out)+2L),bytes=unname(file.info(paths)$size),
    sha256=unname(vapply(paths,digest::digest,"",file=TRUE,algo="sha256")))
  gse_atomic_rds_or_identical(m,file.path(out,"CHECKPOINT_MANIFEST.rds"))
}

csa_verify <- function(root,e,spec) {
  base_verify(root,e,spec)
  out<-gse_output(root,spec$version); m<-readRDS(file.path(out,"CHECKPOINT_MANIFEST.rds"))
  for(i in seq_len(nrow(m))) {
    p<-file.path(out,m$file[i]); gse_assert(file.exists(p) && file.info(p)$size==m$bytes[i] &&
      digest::digest(file=p,algo="sha256")==m$sha256[i],"Checkpoint manifest failed")
  }
  expected <- readRDS(file.path(out,"final/results.rds"))
  expected$accuracy <- 1-expected$classification_error
  gse_assert(identical(expected,readRDS(file.path(out,"final/response_summary.rds"))),
    "Accuracy summary reconstruction failed")
  invisible(TRUE)
}
