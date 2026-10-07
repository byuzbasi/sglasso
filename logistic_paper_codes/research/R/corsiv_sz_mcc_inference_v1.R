# Post-estimation MCC adapter; imports immutable cai_* integrity utilities.
cmi_score <- function(y, prediction) {
  cai_assert(length(y)==length(prediction) && all(y %in% 0:1) &&
    all(prediction %in% 0:1), "Invalid binary score inputs")
  counts <- c(TP=sum(y==1 & prediction==1), FP=sum(y==0 & prediction==1),
    TN=sum(y==0 & prediction==0), FN=sum(y==1 & prediction==0))
  margins <- c(counts[1]+counts[2],counts[1]+counts[4],
    counts[3]+counts[2],counts[3]+counts[4])
  list(counts=counts, mcc=as.numeric(do.call(mltools::mcc,as.list(counts))),
    degenerate=any(margins==0))
}
cmi_inputs <- function(root) {
  x <- cai_inputs(root)
  x$prediction <- vapply(x$rocs,function(r) as.integer(r$predictor>=.5),integer(length(x$y)))
  original <- readRDS(file.path(root,"outputs/study/corsiv_sz_external_v1/final/response_summary.rds"))
  x$observed <- vapply(seq_along(x$methods),function(j) cmi_score(x$y,x$prediction[,j])$mcc,0.0)
  cai_assert(max(abs(x$observed-original$mcc))<1e-12,"Frozen MCC reconstruction failed")
  x
}
cmi_spec <- function(root,x,stage) {
  files <- c("R/corsiv_sz_auc_inference_v1.R","R/corsiv_sz_mcc_inference_v1.R",
    "scripts/114_bootstrap_corsiv_mcc_v1.R","tests/test_corsiv_sz_mcc_inference_v1.R",
    "CORSIV_SZ_MCC_INFERENCE_PROTOCOL_V1.md")
  s <- list(schema="corsiv_mcc_inference_v1",stage=stage,
    run_id=if(stage=="smoke") "corsiv_sz_mcc_smoke_v1" else "corsiv_sz_mcc_b10000_v1",
    boot_n=if(stage=="smoke") 50L else 10000L,
    chunk_size=if(stage=="smoke") 25L else 100L,
    seed_base=1010610L,rng_kind=c("Mersenne-Twister","Inversion","Rejection"),
    threshold=.5,positive_rule="probability >= threshold",stratified=TRUE,paired=TRUE,
    interval="percentile",quantile_type=7L,confidence=.95,
    difference_confidence_bonferroni=.99,comparison_family_size=5L,
    zero_denominator="mltools_zero_retained_with_flag",model_fits=0L,
    methods=x$methods,input_identity=x$identity,sources=cai_manifest(root,files),
    runtime=list(R=R.version.string,platform=R.version$platform,
      packages=vapply(c("mltools","pROC","digest","jsonlite","knitr"),
        function(p) as.character(utils::packageVersion(p)),"")))
  s$signature <- digest::digest(s,algo="sha256"); s
}
cmi_units <- function(s) sprintf("chunk_%04d",seq_len(ceiling(s$boot_n/s$chunk_size)))
cmi_ids <- function(u,s) {
  k <- as.integer(sub("chunk_","",u))
  seq.int((k-1L)*s$chunk_size+1L,min(k*s$chunk_size,s$boot_n))
}
cmi_draw <- function(x,s,id) {
  do.call(RNGkind,as.list(s$rng_kind)); set.seed(s$seed_base+id)
  # The same indices feed ALL methods: do not sample methods independently.
  unlist(lapply(0:1,function(a) {i <- which(x$y==a); sample(i,length(i),replace=TRUE)}),use.names=FALSE)
}
cmi_compute <- function(u,x,s) {
  ids <- cmi_ids(u,s); n <- length(ids)
  counts <- array(0,c(n,6L,4L),dimnames=list(NULL,x$methods,c("TP","FP","TN","FN")))
  mcc <- matrix(NA_real_,n,6L,dimnames=list(NULL,x$methods))
  degenerate <- matrix(FALSE,n,6L,dimnames=dimnames(mcc)); hashes <- character(n)
  for(i in seq_along(ids)) {
    rows <- cmi_draw(x,s,ids[i]); hashes[i] <- digest::digest(rows,algo="sha256")
    for(j in 1:6) {
      z <- cmi_score(x$y[rows],x$prediction[rows,j])
      counts[i,j,] <- z$counts; mcc[i,j] <- z$mcc; degenerate[i,j] <- z$degenerate
    }
  }
  list(ids=ids,counts=counts,mcc=mcc,degenerate=degenerate,index_sha256=hashes)
}
cmi_check <- function(v,u,x,s) {
  n <- length(cmi_ids(u,s)); c <- v$counts
  cai_assert(identical(v$ids,cmi_ids(u,s)) && identical(dim(c),c(n,6L,4L)) &&
    identical(dim(v$mcc),c(n,6L)) && identical(dim(v$degenerate),c(n,6L)) &&
    identical(colnames(v$mcc),x$methods) && is.logical(v$degenerate) &&
    all(is.finite(c)) && all(c>=0 & c==floor(c)) &&
    all(c[,,1]+c[,,4]==353) && all(c[,,2]+c[,,3]==322),"Invalid bootstrap counts")
  denominator <- sqrt((c[,,1]+c[,,2])*(c[,,1]+c[,,4])*(c[,,3]+c[,,2])*(c[,,3]+c[,,4]))
  zero <- denominator==0; denominator[zero] <- 1
  expected <- (c[,,1]*c[,,3]-c[,,2]*c[,,4])/denominator
  cai_assert(all(is.finite(v$mcc)) && all(abs(v$mcc)<=1+1e-14) &&
    max(abs(expected-v$mcc))<1e-14 && all(zero==v$degenerate) &&
    length(v$index_sha256)==n && all(grepl("^[a-f0-9]{64}$",v$index_sha256)),"MCC/degeneracy validation failed")
  invisible(TRUE)
}
cmi_valid <- function(out,u,x,s) {
  p <- cai_path(out,u)
  if(!file.exists(p) || !file.exists(paste0(p,".receipt.rds"))) return(FALSE)
  tryCatch({r <- readRDS(paste0(p,".receipt.rds")); z <- readRDS(p)
    identical(r$sha256,cai_hash(p)) && r$bytes==file.info(p)$size &&
      identical(z$signature,s$signature) && identical(z$unit,u) &&
      isTRUE(cmi_check(z$value,u,x,s))},error=function(e) FALSE)
}
cmi_one <- function(out,u,x,s,work=function() cmi_compute(u,x,s)) {
  if(cmi_valid(out,u,x,s)) return(invisible("resumed"))
  p <- cai_path(out,u)
  cai_assert(!file.exists(p) && !file.exists(paste0(p,".receipt.rds")),"Invalid existing shard; preserve and inspect")
  t <- Sys.time()
  v <- withCallingHandlers(work(),warning=function(w) stop(conditionMessage(w),call.=FALSE))
  cmi_check(v,u,x,s)
  cai_save(list(signature=s$signature,unit=u,value=v,
    seconds=as.numeric(difftime(Sys.time(),t,units="secs"))),p)
  cai_save(list(bytes=file.info(p)$size,sha256=cai_hash(p)),paste0(p,".receipt.rds"))
  invisible("passed")
}
cmi_tables <- function(out,x,s) {
  zs <- lapply(cmi_units(s),function(u) readRDS(cai_path(out,u))$value)
  ids <- unlist(lapply(zs,`[[`,"ids")); cai_assert(identical(ids,seq_len(s$boot_n)),"Missing/duplicated draws")
  b <- do.call(rbind,lapply(zs,`[[`,"mcc")); g <- do.call(rbind,lapply(zs,`[[`,"degenerate"))
  q <- function(v,p) as.numeric(quantile(v,p,type=s$quantile_type,names=FALSE,na.rm=FALSE))
  intervals <- do.call(rbind,lapply(1:6,function(j) {
    ci <- q(b[,j],c(.025,.975)); data.frame(method=x$methods[j],mcc=x$observed[j],
      lower95=ci[1],upper95=ci[2],degenerate_draws=sum(g[,j]),boot_n=s$boot_n)
  }))
  contrasts <- do.call(rbind,lapply(2:6,function(j) {
    delta <- b[,1]-b[,j]; ci <- q(delta,c(.025,.975,.005,.995))
    data.frame(reference=x$methods[1],comparator=x$methods[j],
      mcc_difference=x$observed[1]-x$observed[j],lower95=ci[1],upper95=ci[2],
      lower99_bonferroni=ci[3],upper99_bonferroni=ci[4],
      either_degenerate_draws=sum(g[,1] | g[,j]),boot_n=s$boot_n)
  }))
  list(mcc_intervals=intervals,paired_mcc_differences=contrasts)
}
cmi_finish <- function(out,x,s) {
  tables <- cmi_tables(out,x,s)
  for(n in names(tables)) {
    cai_same_or_save(tables[[n]],file.path(out,"final",paste0(n,".rds")))
    p <- file.path(out,"final",paste0(n,".tex"))
    lines <- strsplit(as.character(knitr::kable(tables[[n]],format="latex",booktabs=TRUE,
      digits=6,row.names=FALSE)),"\n",fixed=TRUE)[[1]]
    if(file.exists(p)) cai_assert(identical(readLines(p),lines),"Changed table") else writeLines(lines,p)
  }
  files <- c("specification.rds",paste0("shards/",list.files(file.path(out,"shards"))),
    paste0("final/",list.files(file.path(out,"final"))))
  cai_same_or_save(cai_manifest(out,files),file.path(out,"MANIFEST.rds"))
  cai_same_or_save(list(signature=s$signature,boot_n=s$boot_n,units=length(cmi_units(s)),
    model_fits=0L,manifest_sha256=cai_hash(file.path(out,"MANIFEST.rds"))),file.path(out,"COMPLETE.rds"))
}
cmi_verify <- function(out,x,s) {
  cai_assert(identical(readRDS(file.path(out,"specification.rds")),s),"Changed specification")
  cai_assert(all(vapply(cmi_units(s),function(u) cmi_valid(out,u,x,s),logical(1))),"Incomplete/invalid chunks")
  cai_check_manifest(out,readRDS(file.path(out,"MANIFEST.rds")))
  done <- readRDS(file.path(out,"COMPLETE.rds"))
  cai_assert(identical(done$signature,s$signature) && done$boot_n==s$boot_n &&
    done$units==length(cmi_units(s)) && done$model_fits==0L &&
    identical(done$manifest_sha256,cai_hash(file.path(out,"MANIFEST.rds"))),"Invalid completion")
  tables <- cmi_tables(out,x,s)
  for(n in names(tables)) cai_assert(identical(tables[[n]],readRDS(file.path(out,"final",paste0(n,".rds")))),"Changed summary")
  invisible(TRUE)
}
cmi_progress <- function(out,s,status,start,phase,durations=numeric(),current=NULL) {
  now <- Sys.time(); elapsed <- as.numeric(difftime(now,start,units="secs"))
  completed <- sum(status %in% c("passed","resumed")); failed <- sum(status=="failed")
  running <- as.integer(!is.null(current)); remaining <- sum(status=="pending")
  eta <- if(!remaining) 0 else if(length(durations)) mean(durations)*remaining else NULL
  deadline <- suppressWarnings(as.numeric(Sys.getenv("SLURM_JOB_END_TIME","")))
  p <- list(run_id=s$run_id,slurm_job_id=Sys.getenv("SLURM_JOB_ID","local"),
    phase=phase,work_unit=paste(s$chunk_size,"paired bootstrap draws, all six methods"),
    total=length(status),completed=completed,running=running,failed=failed,pending=remaining-running,
    percent_completed=100*completed/length(status),elapsed_seconds=elapsed,
    throughput_units_per_second=if(elapsed>0) sum(status=="passed")/elapsed else 0,
    eta_seconds=eta,estimated_completion_utc=if(is.null(eta)) NULL else format(now+eta,tz="UTC",usetz=TRUE),
    eta_status=if(is.null(eta)) "estimating" else "measured_current_invocation",
    eta_scope="remaining_resampling_only_excludes_final_aggregation_verification",
    slurm_remaining_seconds=if(is.finite(deadline)) max(0,deadline-as.numeric(now)) else NULL,
    heartbeat_utc=format(now,tz="UTC",usetz=TRUE),current_work_unit=current,
    last_completed_work_unit=if(completed) tail(names(status)[status %in% c("passed","resumed")],1L) else NULL)
  tmp <- tempfile("cmi_progress_",out); jsonlite::write_json(p,tmp,auto_unbox=TRUE,null="null",pretty=TRUE)
  cai_assert(file.rename(tmp,file.path(out,"progress.json")),"Progress rename failed")
  row <- as.data.frame(p[c("heartbeat_utc","phase","total","completed","running","failed","pending",
    "percent_completed","elapsed_seconds","throughput_units_per_second")])
  row$eta_seconds <- if(is.null(eta)) NA_real_ else eta
  hist <- file.path(out,"progress.tsv")
  write.table(row,hist,sep="\t",quote=FALSE,row.names=FALSE,col.names=!file.exists(hist),append=file.exists(hist))
  cat(p$heartbeat_utc,phase,completed,"/",length(status),"validated; ETA:",
    if(is.null(eta)) "estimating" else round(eta),"seconds\n"); flush.console()
}
