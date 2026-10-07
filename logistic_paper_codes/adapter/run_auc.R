run <- function() {
  cai_assert(!dir.exists(out) || resume,"Output exists; use --resume or --mode=verify")
  dir.create(file.path(out,"shards"),recursive=TRUE,showWarnings=FALSE)
  lock <- file.path(out,"RUNNING.lock")
  cai_assert(dir.create(lock,showWarnings=FALSE),"Active/stale lock; inspect, do not delete blindly")
  on.exit(unlink(lock,recursive=TRUE))
  cai_same_or_save(spec,file.path(out,"specification.rds"))
  if(file.exists(file.path(out,"COMPLETE.rds"))) {
    cai_verify(out,input,spec); cat("Already completed; verified without resampling.\n"); return(invisible(NULL))
  }
  status <- setNames(rep("pending",11L),cai_units()); durations <- numeric()
  for(u in names(status)) if(cai_valid(cai_path(out,u),u,input,spec)) status[u] <- "resumed"
  started <- Sys.time(); fresh <- 0L
  cai_progress(out,spec,status,started,"preflight_complete")
  for(u in names(status)[status=="pending"]) {
    cai_progress(out,spec,status,started,"resampling",u,durations)
    child <- parallel::mcparallel(tryCatch({cai_one(u,out,input,spec); list(ok=TRUE)},
      error=function(e) list(ok=FALSE,message=conditionMessage(e))),silent=TRUE)
    repeat {
      got <- parallel::mccollect(child,wait=FALSE)
      if(!is.null(got)) break
      Sys.sleep(interval)
      cai_progress(out,spec,status,started,"resampling",u,durations)
    }
    result <- got[[1L]]
    if(!is.list(result) || !isTRUE(result$ok) || !cai_valid(cai_path(out,u),u,input,spec)) {
      status[u] <- "failed"
      cai_save(list(unit=u,result=result),file.path(out,"failures",paste0(u,"_",Sys.getpid(),"_",format(Sys.time(),"%Y%m%dT%H%M%S"),".rds")))
      cai_progress(out,spec,status,started,"failed",durations=durations)
      stop("Inference unit failed; records preserved: ",u,"; ",if(is.list(result)) result$message else "worker failure")
    }
    status[u] <- "passed"; fresh <- fresh+1L
    durations[u] <- readRDS(cai_path(out,u))$runtime_seconds
    cai_progress(out,spec,status,started,"unit_validated",durations=durations)
    if(fresh>=limit && any(status=="pending")) {
      cai_progress(out,spec,status,started,"paused",durations=durations)
      cat("Planned pause. Resume with --resume; max-units is an execution-only cap.\n")
      return(invisible(NULL))
    }
  }
  cai_progress(out,spec,status,started,"aggregating_and_verifying",durations=durations)
  cai_finalize(out,input,spec); cai_verify(out,input,spec)
  cai_progress(out,spec,status,started,"completed",durations=durations)
  cat("Complete: TRUE\nOutput:",out,"\n")
}
run()
