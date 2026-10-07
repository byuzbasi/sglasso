# Run-code distribution; data are supplied explicitly, never bundled.
here <- dirname(sys.frame(1L)$ofile)
sys.source(file.path(here,"base.R"),environment())
sys.source(file.path(here,"corsivsz_loader.R"),environment())
paper_inputs <- new.env(parent=emptyenv())
paper_configuration <- function(bundle,name) {
  z <- parse(file.path(bundle,"configuration",paste0(name,".R")))
  portable_assert(length(z)==1L,"Invalid configuration file")
  eval(z[[1L]],envir=baseenv())
}
portable_verify <- function(bundle) {
  portable_require("digest")
  m <- read.csv(file.path(bundle,"MANIFEST.csv"),stringsAsFactors=FALSE)
  portable_assert(isTRUE(all.equal(m,portable_inventory(bundle),tolerance=0,check.attributes=FALSE)),
                  "Code manifest/size/checksum mismatch")
  portable_assert(identical(strsplit(readLines(file.path(bundle,"MANIFEST.sha256"))," +")[[1]][1],
                            portable_hash(file.path(bundle,"MANIFEST.csv"))),"Manifest checksum mismatch")
  portable_assert(!any(nzchar(Sys.readlink(file.path(bundle,m$file)))),"Symlink files refused")
  map <- read.csv(file.path(bundle,"SOURCE_MAP.csv"))
  portable_assert(all(vapply(seq_len(nrow(map)),function(i)
    portable_hash(file.path(bundle,map$file[i]))==map$original_sha256[i],logical(1))),
    "Original numerical source changed")
  sim <- paper_configuration(bundle,"simulation")
  portable_assert(nrow(sim$tasks)==400L && nrow(sim$design)==8L &&
    identical(sim$configuration$alpha_grid,seq(0,1,.1)) &&
    identical(sim$configuration$d_grid,seq(0,1,.1)),"Simulation settings invalid")
  cat("VERIFIED: code-only manifest and original research sources/settings. Fits=0; draws=0.\n")
  invisible(TRUE)
}
portable_data <- function(bundle,stage) {
  path <- paper_inputs$data
  portable_assert(length(path)==1L && file.exists(path),
                  "Supply --data=/absolute/path/CoRSIVSZ_v1.rds; no automatic download")
  d <- as.list(read.dcf(file.path(bundle,"configuration/CoRSIVSZ_v1.dcf"))[1,])
  original_data <- .corsivsz_read_verified(path,d)
  x <- list(train=original_data$development,test=original_data$external,
            group=original_data$group,group_name=original_data$group_name,probe_id=original_data$probe_id)
  original <- paper_configuration(bundle,"external")$data_identity
  portable_assert(identical(x$train$sample_id,original$train_sample_id) &&
    identical(x$test$sample_id,original$test_sample_id) &&
    identical(digest::digest(x$train$y,algo="sha256"),original$train_response_sha256) &&
    identical(digest::digest(x$test$y,algo="sha256"),original$test_response_sha256) &&
    identical(digest::digest(x$group,algo="sha256"),original$group_sha256),"Cohort identities changed")
  if(stage=="smoke") {
    chosen <- sort(vapply(c(2L,3L,4L),function(s) which(tabulate(x$group)==s)[1L],integer(1)))
    cols <- which(x$group %in% chosen)
    for(part in c("train","test")) {
      z <- x[[part]]; n <- if(part=="train") 24L else 8L
      rows <- unlist(lapply(0:1,function(y) head(which(z$y==y),n)))
      for(key in c("y","sample_id","geo_sample_id")) z[[key]] <- z[[key]][rows]
      z$X <- z$X[rows,cols,drop=FALSE]; x[[part]] <- z
    }
    x$group <- match(x$group[cols],chosen); x$group_name <- x$group_name[chosen]; x$probe_id <- x$probe_id[cols]
  }
  x$processed_sha256 <- portable_hash(path); x
}
portable_inference_input <- function(bundle,stage) {
  out <- paper_inputs$external
  portable_assert(length(out)==1L && dir.exists(out),"Supply --external=/absolute/path/to/a/verified/external/run")
  source(file.path(bundle,"research/R/corsiv_sz_auc_inference_v1.R"),local=environment())
  m <- readRDS(file.path(out,"CHECKPOINT_MANIFEST.rds")); cai_check_manifest(out,m)
  spec <- readRDS(file.path(out,"specification.rds")); complete <- readRDS(file.path(out,"COMPLETE.rds"))
  portable_assert(spec$stage==stage && complete$validated_units==spec$configuration$inner_folds+1L &&
    complete$method_rows==6L && identical(complete$scientific_signature,spec$scientific_signature),
    "External run is incomplete, unverified or the wrong stage")
  z <- readRDS(file.path(out,"final/predictions.rds"))
  r <- readRDS(file.path(out,"final/response_summary.rds")); methods <- r$method
  ids <- spec$data_identity$test_sample_id
  portable_assert(length(methods)==6L && !anyDuplicated(methods) && !anyDuplicated(ids) &&
    nrow(z)==6L*length(ids) && setequal(unique(z$method),methods),"Invalid fixed predictions")
  parts <- lapply(methods,function(method) {
    p <- z[z$method==method,,drop=FALSE]
    portable_assert(nrow(p)==length(ids) && !anyDuplicated(p$sample_id) && setequal(p$sample_id,ids),"Unpaired IDs")
    p[match(ids,p$sample_id),,drop=FALSE]
  })
  y <- parts[[1L]]$y
  portable_assert(all(y %in% 0:1) && all(table(factor(y,levels=0:1))>0) &&
    all(vapply(parts,function(p) identical(p$y,y) && all(is.finite(p$probability)) &&
      all(p$probability>=0 & p$probability<=1),logical(1))),"Invalid responses or probabilities")
  if(stage=="production") {
    original <- paper_configuration(bundle,"external")
    portable_assert(identical(ids,original$data_identity$test_sample_id) &&
      length(y)==675L && sum(y==1)==353L && sum(y==0)==322L &&
      identical(spec$configuration,original$configuration) && identical(spec$fold,original$fold),
      "Production inference requires the original external cohort/settings/folds")
  }
  rocs <- setNames(lapply(parts,function(p) pROC::roc(y,p$probability,levels=c(0,1),direction="<",quiet=TRUE,na.rm=FALSE)),methods)
  auc <- vapply(rocs,function(r) as.numeric(pROC::auc(r)),0.0)
  prediction <- vapply(parts,function(p) as.integer(p$probability>=.5),integer(length(y)))
  observed <- vapply(seq_along(methods),function(j) mltools::mcc(preds=prediction[,j],actuals=y),0.0)
  portable_assert(max(abs(auc-r$auc))<1e-12 && max(abs(observed-r$mcc))<1e-12,"AUC/MCC reconstruction failed")
  list(rocs=rocs,methods=methods,auc=auc,ids=ids,y=y,prediction=prediction,observed=observed,
       identity=list(parent_signature=spec$scientific_signature,
         manifest_sha256=portable_hash(file.path(out,"CHECKPOINT_MANIFEST.rds")),
         predictions_sha256=portable_hash(file.path(out,"final/predictions.rds"))))
}
# Only the two hardcoded cohort margins in the original MCC checker are adapted
# for the 16-person smoke fixture. Production uses the original checker verbatim.
paper_smoke_mcc_check <- function(f) {
  replace <- function(z) {
    if(is.numeric(z) && length(z)==1L && z==353) return(quote(sum(x$y==1)))
    if(is.numeric(z) && length(z)==1L && z==322) return(quote(sum(x$y==0)))
    if(is.call(z)) for(i in seq_along(z)[-1L]) z[i] <- list(replace(z[[i]]))
    z
  }
  body(f) <- replace(body(f)); f
}
portable_main <- function(bundle,args) {
  portable_require(c("digest","knitr","pROC","mltools"))
  portable_assert(all(grepl("^--(mode|study|stage|output|cores|max-units|interval|validation|data|external)=|^--resume$",args)),"Unknown option")
  arg <- function(n,d) {
    x <- args[startsWith(args,paste0("--",n,"="))]; portable_assert(length(x)<=1L,"Duplicate option")
    if(length(x)) sub("^[^=]+=","",x) else d
  }
  paper_inputs$data <- arg("data",""); paper_inputs$external <- arg("external","")
  mode <- arg("mode","verify"); study <- arg("study","simulation")
  portable_assert(mode %in% c("verify","preflight","validate","smoke","production","verify-run"),"Invalid mode")
  portable_assert(study %in% c("simulation","external","auc","mcc"),"Invalid study")
  portable_verify(bundle)
  if(mode=="verify") return(invisible(TRUE))
  out <- portable_output(bundle,arg("output",""))
  if(mode=="preflight") {
    if(study=="external") portable_data(bundle,"production")
    if(study %in% c("auc","mcc")) portable_inference_input(bundle,"production")
    cat("PREFLIGHT: code/settings and required inputs checked; no fits or resampling.\n")
    return(invisible(TRUE))
  }
  portable_require(c("jsonlite","Rcpp","RcppArmadillo","adelie","grpreg","logistf"))
  portable_assert(.Platform$OS.type=="unix","Numerical execution requires Unix forks")
  for(n in c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"))
    portable_assert(Sys.getenv(n)=="1",paste("Set",n,"=1 before starting R"))
  if(mode=="validate") {sys.source(file.path(bundle,"adapter/tests.R"),environment()); return(portable_tests(bundle,out))}
  stage <- if(mode=="smoke") "smoke" else if(mode=="production") "production" else arg("stage","production")
  portable_assert(stage %in% c("smoke","production"),"Invalid stage")
  cores <- as.integer(arg("cores","1")); interval <- as.numeric(arg("interval","30"))
  limit <- as.numeric(arg("max-units","Inf")); resume <- "--resume" %in% args
  portable_assert(!is.na(cores) && cores>=1 && cores<=56 && is.finite(interval) && interval>=1 && interval<=60 &&
    !is.na(limit) && limit>=1 && (is.infinite(limit) || limit==floor(limit)),"Invalid execution options")
  if(stage=="smoke") portable_assert(cores<=2,"Smoke tests allow at most two workers")
  allocated <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  portable_assert(is.na(allocated) || cores<=allocated,"Workers exceed allocated CPUs")
  if(mode=="production") {
    receipt <- arg("validation",""); portable_assert(file.exists(receipt),"Local validation receipt required")
    r <- readRDS(receipt)
    portable_assert(isTRUE(r$accepted) && all(r$checks$passed) &&
      identical(r$distribution_manifest_sha256,portable_hash(file.path(bundle,"MANIFEST.csv"))),"Code validation mismatch")
    portable_assert(identical(r$R,R.version.string) && identical(r$platform,R.version$platform) &&
      identical(r$packages,vapply(names(r$packages),function(p) as.character(packageVersion(p)),"")),"Validation runtime differs")
  }
  version <- paste0("portable_",study,"_",stage,"_v2"); run_out <- file.path(out,version)
  portable_assert(mode=="verify-run" || !dir.exists(run_out) || resume,"Output exists; use identical flags plus --resume")
  root <- file.path(bundle,"research")
  if(study=="simulation") {
    e <- portable_simulation(bundle,out,stage)
    if(mode=="verify-run") return(e$lsg_verify_v7(root,stage,version))
    answer <- e$lsg_run_v19(root,version,cores=cores,max_seconds=if(stage=="smoke") 120 else 18000,
      task_limit=limit,interval=interval,stage=stage)
    portable_assert(!any(vapply(answer$outcomes,function(x) x$status %in% c("gate_failed","error"),logical(1))),"Numerical gate failed; evidence preserved")
    return(invisible(answer))
  }
  scope <- new.env(parent=environment())
  if(study=="external") {
    x <- portable_external(bundle,out,stage,version)
    for(n in names(x)) scope[[n]] <- x[[n]]
    scope$root <- root; scope$out <- run_out; scope$resume <- resume; scope$interval <- interval; scope$max_units <- limit
    if(mode=="verify-run") return(x$b$gse_verify(root,x$e,x$spec))
  } else {
    for(f in c("corsiv_sz_auc_inference_v1.R","corsiv_sz_mcc_inference_v1.R")) sys.source(file.path(root,"R",f),scope)
    if(stage=="smoke" && study=="mcc") scope$cmi_check <- paper_smoke_mcc_check(scope$cmi_check)
    input <- portable_inference_input(bundle,stage); spec <- paper_configuration(bundle,study)
    spec$stage <- stage; spec$run_id <- version; spec$input_identity <- input$identity; spec$sources <- portable_inventory(bundle)
    spec$runtime <- list(R=R.version.string,platform=R.version$platform,packages=vapply(
      if(study=="auc") c("pROC","digest","jsonlite","knitr") else c("mltools","pROC","digest","jsonlite","knitr"),
      function(p) as.character(packageVersion(p)),""))
    if(stage=="smoke") {spec$boot_n <- if(study=="auc") 25L else 50L; if(study=="mcc") spec$chunk_size <- 25L}
    spec$signature <- digest::digest(spec,algo="sha256")
    scope$input <- scope$x <- input; scope$spec <- scope$s <- spec
    scope$out <- run_out; scope$resume <- resume; scope$interval <- interval; scope$limit <- limit
    if(mode=="verify-run") return(if(study=="auc") scope$cai_verify(run_out,input,spec) else scope$cmi_verify(run_out,input,spec))
  }
  sys.source(file.path(bundle,"adapter",paste0("run_",study,".R")),scope)
}
