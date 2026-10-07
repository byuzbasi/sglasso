run<-function() {
  if(dir.exists(out) && !resume) stop("Output exists; use --resume or --mode=verify")
  dir.create(file.path(out,"shards"),recursive=TRUE,showWarnings=FALSE)
  lock<-file.path(out,"RUNNING.lock")
  if(!dir.create(lock,showWarnings=FALSE)) stop("Active/stale lock; inspect before recovery")
  on.exit(unlink(lock,recursive=TRUE))
  b$gse_atomic_rds_or_identical(spec,file.path(out,"specification.rds"))
  if(file.exists(file.path(out,"CHECKPOINT_MANIFEST.rds"))) {
    b$gse_verify(root,e,spec); cat("Already complete; verified without refitting.\n"); return(invisible(NULL))
  }
  units<-c(paste0("fold_",seq_len(spec$configuration$inner_folds)),"refit")
  status<-setNames(rep("pending",length(units)),units)
  for(u in units) if(b$gse_valid_shard(b$gse_shard(out,u),spec,u)) status[u]<-"resumed"
  start<-Sys.time(); done<-0L
  b$gse_progress(out,spec,status,"preflight_complete",start)
  for(k in seq_len(spec$configuration$inner_folds)) {
    u<-paste0("fold_",k)
    if(status[u]=="resumed") next
    b$gse_run_unit(u,out,spec,function() b$gse_fit_fold(e,data,spec,k),status,start,interval)
    status[u]<-"passed"; done<-done+1L
    if(done>=max_units) {b$gse_progress(out,spec,status,"paused",start); cat("Planned pause. Use identical flags plus --resume.\n"); return(invisible(NULL))}
  }
  selected<-b$gse_select_from_folds(e,data,spec,out)
  if(!all(vapply(selected$selections,function(x) nrow(x)==1L && x$eligible[1L] %in% TRUE,logical(1)))) stop("Invalid CV selection")
  if(status["refit"]!="resumed") {
    b$gse_run_unit("refit",out,spec,function() b$gse_fit_refit(e,data,spec,selected),status,start,interval)
    status["refit"]<-"passed"
  }
  b$gse_progress(out,spec,status,"aggregating_and_verifying",start)
  b$gse_finalize(e,data,spec,out,selected); b$gse_verify(root,e,spec)
  b$gse_progress(out,spec,status,"completed",start)
  cat("Complete: TRUE\nOutput:",out,"\n")
}
run()
