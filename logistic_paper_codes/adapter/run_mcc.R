run <- function() {
  cai_assert(!dir.exists(out) || resume,"Output exists; use --resume or verify")
  dir.create(file.path(out,"shards"),recursive=TRUE,showWarnings=FALSE)
  lock <- file.path(out,"RUNNING.lock")
  cai_assert(dir.create(lock,showWarnings=FALSE),"Active/stale lock: inspect before resuming")
  on.exit(unlink(lock,recursive=TRUE))
  cai_same_or_save(s,file.path(out,"specification.rds"))
  if(file.exists(file.path(out,"COMPLETE.rds"))) {
    cmi_verify(out,x,s); cat("Already complete; verified without resampling.\n"); return(invisible(NULL))
  }
  status <- setNames(rep("pending",length(cmi_units(s))),cmi_units(s)); durations <- numeric()
  for(u in names(status)) if(cmi_valid(out,u,x,s)) status[u] <- "resumed"
  start <- Sys.time(); fresh <- 0L
  cmi_progress(out,s,status,start,"preflight_complete")
  for(u in names(status)[status=="pending"]) {
    cmi_progress(out,s,status,start,"resampling",durations,u)
    tryCatch(cmi_one(out,u,x,s),error=function(e) {
      status[u] <<- "failed"
      cai_save(list(unit=u,message=conditionMessage(e)),file.path(out,"failures",
        paste0(u,"_",Sys.getpid(),"_",format(Sys.time(),"%Y%m%dT%H%M%S"),".rds")))
      cmi_progress(out,s,status,start,"failed",durations); stop(e)
    })
    status[u] <- "passed"; durations[u] <- readRDS(cai_path(out,u))$seconds; fresh <- fresh+1L
    cmi_progress(out,s,status,start,"unit_validated",durations)
    if(fresh>=limit && any(status=="pending")) {
      cmi_progress(out,s,status,start,"paused",durations)
      cat("Planned pause; resume with --resume (omit --max-chunks to finish).\n"); return(invisible(NULL))
    }
  }
  cmi_progress(out,s,status,start,"aggregating_and_verifying",durations)
  cmi_finish(out,x,s); cmi_verify(out,x,s)
  cmi_progress(out,s,status,start,"completed",durations)
  cat("Complete: TRUE\nOutput:",out,"\n")
}
run()
