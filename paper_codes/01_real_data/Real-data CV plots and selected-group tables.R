# ============================================================
# Real-data CV plots and selected-group tables
# Datasets: Birthwt, bardet, GenAtHum
# Methods: GLASSO, ADELIE, GENET, SGLASSO, GSCAD, GMCP
#
# This script is for full-data CV visualization and selected-group tables.
# It does not perform train/test splitting and is not used for the repeated
# real-data prediction-performance tables.
# ============================================================

library(grpreg)
library(sglasso)
library(kableExtra)

# Do NOT call library(grpnet) or library(adelie)
# Use grpnet:: and adelie:: to avoid S3 method conflicts.

# ------------------------------------------------------------
# 1. General settings
# ------------------------------------------------------------

#set.seed(2025)

nfolds <- 5
maxit <- 1e8
eps <- 1e-4
lambda_min_ratio <- 0.01
nlambda <- 100
alpha_seq <- seq(0.1, 0.9, 0.1)
nd <- 11
standardize <- TRUE
sglasso_screen <- "SSR_fast"
sglasso_transform <- "lazy"

base_output_dir <- file.path("results", "real_data")
cv_plot_dir <- file.path(base_output_dir, "cv_plots")
selected_group_dir <- file.path(base_output_dir, "selected_groups")

dir.create(cv_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(selected_group_dir, recursive = TRUE, showWarnings = FALSE)

check_symbol <- "$\\checkmark$"
dash_symbol  <- "-"

# ------------------------------------------------------------
# 2. Helper functions for CV plots
# ------------------------------------------------------------

plot_cv_generic <- function(
    log_lambda,
    cv_error,
    cv_se,
    lambda_min,
    lambda_1se = NULL,
    group_counts = NULL,
    title_text = "",
    xlab = expression(log(lambda)),
    ylab = "Cross-validation error",
    col = "red",
    reverse_x = TRUE,
    selected = TRUE,
    vertical.line = TRUE,
    cex.lab = 1.4,
    cex.axis = 1.2,
    cex.main = 1.4
) {

  ind <- is.finite(log_lambda) &
    is.finite(cv_error) &
    is.finite(cv_se)

  x <- log_lambda[ind]
  y <- cv_error[ind]
  se <- cv_se[ind]

  if (!is.null(group_counts)) {
    group_counts <- group_counts[ind]
  }

  ylim <- range(c(y - se, y + se), na.rm = TRUE)

  xlim <- range(x, na.rm = TRUE)
  if (reverse_x) {
    xlim <- rev(xlim)
  }

  plot(
    x, y,
    type = "n",
    xlab = xlab,
    ylab = ylab,
    xlim = xlim,
    ylim = ylim,
    las = 1,
    bty = "n",
    cex.lab = cex.lab,
    cex.axis = cex.axis
  )

  arrows(
    x0 = x,
    x1 = x,
    y0 = y - se,
    y1 = y + se,
    code = 3,
    angle = 90,
    col = "gray80",
    length = 0.05
  )

  points(x, y, col = col, pch = 19, cex = 0.55)
  lines(x, y, col = "black", lwd = 0.6)

  if (vertical.line) {
    abline(v = log(lambda_min), lty = 2, lwd = 0.9)

    if (!is.null(lambda_1se)) {
      abline(v = log(lambda_1se), lty = 3, lwd = 0.9)
    }
  }

  if (selected && !is.null(group_counts)) {
    axis(
      side = 3,
      at = x,
      labels = group_counts,
      tick = FALSE,
      line = -0.5,
      cex.axis = 1.0
    )

    mtext(
      "Groups selected",
      side = 3,
      line = 1.5,
      cex = 0.8
    )
  }

  title(
    title_text,
    adj = 0,
    line = 2.5,
    cex.main = cex.main
  )
}

get_grpreg_group_counts <- function(cvmod, group_membership, tol = 1e-8) {

  beta_path <- coef(cvmod$fit)

  if (is.list(beta_path)) {
    beta_path <- do.call(cbind, beta_path)
  }

  beta_path <- as.matrix(beta_path)
  beta_path <- beta_path[-1, , drop = FALSE]

  sapply(seq_len(ncol(beta_path)), function(j) {
    beta_j <- beta_path[, j]
    active_groups <- tapply(abs(beta_j) > tol, group_membership, any)
    sum(active_groups)
  })
}

get_adelie_group_counts <- function(cv_adelie, group_membership, tol = 1e-8) {

  beta_path <- cv_adelie$grpnet.fit$state$betas
  beta_path <- as.matrix(beta_path)

  sapply(seq_len(nrow(beta_path)), function(j) {
    beta_j <- beta_path[j, ]
    active_groups <- tapply(abs(beta_j) > tol, group_membership, any)
    sum(active_groups)
  })
}

get_genet_group_counts <- function(cv_genet, group_membership, tol = 1e-8) {

  beta_path <- cv_genet$grpnet.fit$beta
  beta_path <- as.matrix(beta_path)

  if (nrow(beta_path) != length(group_membership) &&
      ncol(beta_path) == length(group_membership)) {
    beta_path <- t(beta_path)
  }

  sapply(seq_len(ncol(beta_path)), function(j) {
    beta_j <- beta_path[, j]
    active_groups <- tapply(abs(beta_j) > tol, group_membership, any)
    sum(active_groups)
  })
}

get_sglasso_group_counts <- function(cv_sglasso, X) {

  pred_ns <- predict(
    cv_sglasso$fit,
    newx = X,
    lambda = cv_sglasso$lambda,
    d = cv_sglasso$d.min,
    type = "groups"
  )

  if (is.list(pred_ns)) {
    return(sapply(pred_ns, length))
  }

  if (is.matrix(pred_ns) || is.data.frame(pred_ns)) {
    return(apply(pred_ns, 2, function(z) length(unique(z))))
  }

  rep(as.numeric(pred_ns), length(cv_sglasso$lambda))
}

# ------------------------------------------------------------
# 3. Helper functions for selected-group tables
# ------------------------------------------------------------

get_selected_groups_from_beta <- function(beta, group, tol = 1e-6) {

  beta <- as.numeric(beta)

  if (length(beta) == length(group) + 1) {
    beta <- beta[-1]
  }

  if (length(beta) > length(group)) {
    beta <- tail(beta, length(group))
  }

  if (length(beta) != length(group)) {
    stop("beta and group lengths do not match.")
  }

  active <- abs(beta) > tol

  as.numeric(names(which(tapply(active, group, any))))
}

make_selected_group_table <- function(selection_list, group_membership,
                                      group_names = NULL) {

  all_groups <- sort(unique(group_membership))

  if (is.null(group_names)) {
    group_names <- paste0("Group ", all_groups)
  }

  tab <- data.frame(
    Group = group_names,
    stringsAsFactors = FALSE
  )

  for (m in names(selection_list)) {
    tab[[m]] <- ifelse(
      all_groups %in% selection_list[[m]],
      check_symbol,
      dash_symbol
    )
  }

  tab
}

export_selected_group_table <- function(tab, dataset_name, label_name) {

  latex_tab <- tab
  colnames(latex_tab)[1] <- "Group"

  latex_out <- latex_tab |>
    kbl(
      format = "latex",
      booktabs = TRUE,
      longtable = dataset_name == "GenAtHum",
      align = paste0(
        "l",
        paste(rep("c", ncol(latex_tab) - 1), collapse = "")
      ),
      caption = paste0(
        "Selected groups by each method for the ",
        dataset_name,
        " data set."
      ),
      label = label_name,
      escape = FALSE,
      linesep = ""
    ) |>
    kable_styling(
      latex_options = if (dataset_name == "GenAtHum") {
        c("repeat_header")
      } else {
        c("hold_position", "scale_down")
      },
      font_size = if (dataset_name == "GenAtHum") 6 else 8,
      full_width = FALSE
    ) |>
    row_spec(0, bold = TRUE) |>
    footnote(
      general = paste(
        "A check mark indicates that the corresponding group was selected by the method.",
        "A dash indicates that the group was not selected."
      ),
      general_title = "\\\\textit{Note: } ",
      footnote_as_chunk = FALSE,
      threeparttable = TRUE,
      escape = FALSE
    )

  cat(
    latex_out,
    file = file.path(selected_group_dir, paste0("Selected_Groups_", dataset_name, ".tex"))
  )
}

# ------------------------------------------------------------
# 4. Main function for one dataset
# ------------------------------------------------------------

run_realdata_cv_and_selection <- function(
    X,
    y,
    group,
    dataset_name,
    group_names = NULL
) {

  set.seed(26042026)

  cat("\n============================================================\n")
  cat("Dataset:", dataset_name, "\n")
  cat("============================================================\n")

  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)

  group_membership <- as.numeric(group)
  adelie_groups <- which(!duplicated(group_membership))

  foldid <- sample(rep(seq_len(nfolds), length.out = n))

  # ----------------------------------------------------------
  # CV fitting
  # ----------------------------------------------------------

  cv_glasso <- grpreg::cv.grpreg(
    X = X, y = y, group = group_membership, fold = foldid,
    alpha = 1, penalty = "grLasso",
    nlambda = nlambda, lambda.min = lambda_min_ratio,
    eps = eps, max.iter = maxit,
    dfmax = ncol(X), gmax = length(unique(group_membership))
  )

  cv_gscad <- grpreg::cv.grpreg(
    X = X, y = y, group = group_membership, fold = foldid,
    alpha = 1, penalty = "grSCAD",
    nlambda = nlambda, lambda.min = lambda_min_ratio,
    eps = eps, max.iter = maxit,
    dfmax = ncol(X), gmax = length(unique(group_membership))
  )

  cv_gmcp <- grpreg::cv.grpreg(
    X = X, y = y, group = group_membership, fold = foldid,
    alpha = 1, penalty = "grMCP",
    nlambda = nlambda, lambda.min = lambda_min_ratio,
    eps = eps, max.iter = maxit,
    dfmax = ncol(X), gmax = length(unique(group_membership))
  )

  adelie_fits <- lapply(alpha_seq, function(a) {
    adelie::cv.grpnet(
      X = X,
      glm = adelie::glm.gaussian(y = y),
      groups = adelie_groups,
      foldid = foldid,
      alpha = a,
      min_ratio = lambda_min_ratio,
      lmda_path_size = nlambda,
      n_threads = 1,
      tol = eps,
      max_iters = maxit,
      irls_tol = eps,
      irls_max_iters = maxit,
      newton_tol = eps,
      newton_max_iters = maxit,
      early_exit = FALSE
    )
  })

  adelie_cvm <- sapply(adelie_fits, function(fit) {
    min(fit$cvm, na.rm = TRUE)
  })

  adelie_best_ind <- which.min(adelie_cvm)
  adelie_best_alpha <- alpha_seq[adelie_best_ind]
  cv_adelie <- adelie_fits[[adelie_best_ind]]

  genet_fits <- lapply(alpha_seq, function(a) {
    grpnet::cv.grpnet(
      x = X,
      y = y,
      group = group_membership,
      foldid = foldid,
      type.measure = "mse",
      alpha = a,
      thresh = eps,
      maxit = maxit,
      lambda.min.ratio = lambda_min_ratio,
      verbose = FALSE
    )
  })

  genet_cvm <- sapply(genet_fits, function(fit) {
    min(fit$cvm, na.rm = TRUE)
  })

  genet_best_ind <- which.min(genet_cvm)
  genet_best_alpha <- alpha_seq[genet_best_ind]
  cv_genet <- genet_fits[[genet_best_ind]]

  sglasso_fits <- lapply(alpha_seq, function(a) {
    sglasso::cv.sglasso(
      X = X,
      Y = y,
      group = group_membership,
      fold = foldid,
      nlambda = nlambda,
      nd = nd,
      alpha = a,
      lambda.min.ratio = lambda_min_ratio,
      eps = eps,
      max_iter = maxit,
      standardize = standardize,
      screen = sglasso_screen,
      transform = sglasso_transform
    )
  })

  sglasso_cve <- sapply(sglasso_fits, function(fit) {
    min(fit$cve, na.rm = TRUE)
  })

  sglasso_best_ind <- which.min(sglasso_cve)
  sglasso_best_alpha <- alpha_seq[sglasso_best_ind]
  cv_sglasso <- sglasso_fits[[sglasso_best_ind]]

  # ----------------------------------------------------------
  # CV plot
  # ----------------------------------------------------------

  setEPS()
  postscript(
    file = file.path(cv_plot_dir, paste0("CV_Plots_", dataset_name, ".eps")),
    horizontal = FALSE,
    onefile = FALSE,
    paper = "special",
    height = 14,
    width = 14
  )

  layout(
    matrix(
      c(1, 2,
        3, 4,
        5, 6),
      nrow = 3,
      byrow = TRUE
    )
  )

  par(mar = c(4.5, 4.5, 3.5, 1))

  plot_cv_generic(
    log_lambda = log(cv_glasso$lambda),
    cv_error = cv_glasso$cve,
    cv_se = cv_glasso$cvse,
    lambda_min = cv_glasso$lambda.min,
    group_counts = get_grpreg_group_counts(cv_glasso, group_membership),
    title_text = expression(bold("GLASSO"))
  )

  adelie_group_counts <- get_adelie_group_counts(
    cv_adelie,
    group_membership
  )

  plot_cv_generic(
    log_lambda = log(cv_adelie$lambda),
    cv_error = cv_adelie$cvm,
    cv_se = cv_adelie$cvsd,
    lambda_min = cv_adelie$lambda.min,
    group_counts = adelie_group_counts,
    title_text = bquote(
      bold("ADELIE ") *
        (alpha == .(adelie_best_alpha))
    )
  )

  plot_cv_generic(
    log_lambda = log(cv_gscad$lambda),
    cv_error = cv_gscad$cve,
    cv_se = cv_gscad$cvse,
    lambda_min = cv_gscad$lambda.min,
    group_counts = get_grpreg_group_counts(cv_gscad, group_membership),
    title_text = expression(bold("GSCAD"))
  )

  genet_group_counts <- get_genet_group_counts(
    cv_genet,
    group_membership
  )

  plot_cv_generic(
    log_lambda = log(cv_genet$lambda),
    cv_error = cv_genet$cvm,
    cv_se = cv_genet$cvsd,
    lambda_min = cv_genet$lambda.min,
    group_counts = genet_group_counts,
    title_text = bquote(
      bold("GENET ") *
        (alpha == .(genet_best_alpha))
    )
  )

  plot_cv_generic(
    log_lambda = log(cv_gmcp$lambda),
    cv_error = cv_gmcp$cve,
    cv_se = cv_gmcp$cvse,
    lambda_min = cv_gmcp$lambda.min,
    group_counts = get_grpreg_group_counts(cv_gmcp, group_membership),
    title_text = expression(bold("GMCP"))
  )

  d_ind <- cv_sglasso$min_ind["col"]

  sglasso_group_counts <- get_sglasso_group_counts(
    cv_sglasso,
    X
  )

  plot_cv_generic(
    log_lambda = log(cv_sglasso$lambda),
    cv_error = cv_sglasso$cve[, d_ind],
    cv_se = cv_sglasso$cvse[, d_ind],
    lambda_min = cv_sglasso$lambda.min,
    group_counts = sglasso_group_counts,
    title_text = bquote(
      bold("SGLASSO ") *
        (alpha == .(sglasso_best_alpha)) * ", " *
        (d == .(round(cv_sglasso$d.min, 2)))
    )
  )

  dev.off()

  # ----------------------------------------------------------
  # Selected groups from full-data optimal fits
  # ----------------------------------------------------------

  fit_glasso_full <- grpreg::grpreg(
    X = X, y = y, group = group_membership,
    alpha = 1, penalty = "grLasso",
    lambda = cv_glasso$lambda.min,
    eps = eps, max.iter = maxit
  )

  fit_gscad_full <- grpreg::grpreg(
    X = X, y = y, group = group_membership,
    alpha = 1, penalty = "grSCAD",
    lambda = cv_gscad$lambda.min,
    eps = eps, max.iter = maxit
  )

  fit_gmcp_full <- grpreg::grpreg(
    X = X, y = y, group = group_membership,
    alpha = 1, penalty = "grMCP",
    lambda = cv_gmcp$lambda.min,
    eps = eps, max.iter = maxit
  )

  fit_genet_full <- grpnet::grpnet(
    x = X,
    y = y,
    group = group_membership,
    alpha = genet_best_alpha,
    lambda = cv_genet$lambda.min,
    thresh = eps,
    maxit = maxit
  )

  fit_sglasso_full <- sglasso::sglasso(
    X = X,
    Y = y,
    group = group_membership,
    d = cv_sglasso$d.min,
    alpha = sglasso_best_alpha,
    lambda = cv_sglasso$lambda.min,
    eps = eps,
    max_iter = maxit,
    standardize = standardize,
    screen = sglasso_screen,
    transform = sglasso_transform
  )

  adelie_beta_path <- as.matrix(cv_adelie$grpnet.fit$state$betas)
  adelie_lambda_index <- which.min(
    abs(cv_adelie$lambda - cv_adelie$lambda.min)
  )
  beta_adelie <- adelie_beta_path[adelie_lambda_index, ]

  selected_groups_list <- list(
    SGLASSO = get_selected_groups_from_beta(
      coef(fit_sglasso_full),
      group_membership
    ),
    ADELIE = get_selected_groups_from_beta(
      beta_adelie,
      group_membership
    ),
    GLASSO = get_selected_groups_from_beta(
      coef(fit_glasso_full),
      group_membership
    ),
    GENET = get_selected_groups_from_beta(
      fit_genet_full$beta,
      group_membership
    ),
    GSCAD = get_selected_groups_from_beta(
      coef(fit_gscad_full),
      group_membership
    ),
    GMCP = get_selected_groups_from_beta(
      coef(fit_gmcp_full),
      group_membership
    )
  )

  Selected_Groups_Table <- make_selected_group_table(
    selected_groups_list,
    group_membership = group_membership,
    group_names = group_names
  )

  write.csv(
    Selected_Groups_Table,
    file.path(selected_group_dir, paste0("Selected_Groups_", dataset_name, ".csv")),
    row.names = FALSE
  )

  export_selected_group_table(
    Selected_Groups_Table,
    dataset_name = dataset_name,
    label_name = paste0("selected_groups_", dataset_name)
  )

  invisible(
    list(
      cv = list(
        GLASSO = cv_glasso,
        ADELIE = cv_adelie,
        GENET = cv_genet,
        SGLASSO = cv_sglasso,
        GSCAD = cv_gscad,
        GMCP = cv_gmcp
      ),
      selected_groups = Selected_Groups_Table
    )
  )
}

# ============================================================
# 5. Load real datasets and run analyses
# ============================================================

# ------------------------------------------------------------
# Birthwt data from grpreg
# ------------------------------------------------------------

data(Birthwt, package = "grpreg")

X_birthwt <- as.matrix(Birthwt$X)
y_birthwt <- Birthwt$bwt
group_birthwt <- Birthwt$group

result_birthwt <- run_realdata_cv_and_selection(
  X = X_birthwt,
  y = y_birthwt,
  group = group_birthwt,
  dataset_name = "Birthwt"
)

# ------------------------------------------------------------
# Bardet data from gglasso
# ------------------------------------------------------------

if (!requireNamespace("gglasso", quietly = TRUE)) {
  stop("Package 'gglasso' is required for the Bardet data set.", call. = FALSE)
}

data(bardet, package = "gglasso")

X_bardet <- as.matrix(bardet$x)
y_bardet <- bardet$y
group_bardet <- rep(1:20,each=5)

result_bardet <- run_realdata_cv_and_selection(
  X = X_bardet,
  y = y_bardet,
  group = group_bardet,
  dataset_name = "bardet"
)

# ------------------------------------------------------------
# Gene Atlas Human data from sglasso
# ------------------------------------------------------------

data(GenAtHum, package = "sglasso")

X_genathum <- as.matrix(GenAtHum$X)
y_genathum <- GenAtHum$y
group_genathum <- GenAtHum$group

genathum_group_names <- if (!is.null(GenAtHum$groups_name)) {
  GenAtHum$groups_name
} else {
  paste0("Group ", sort(unique(group_genathum)))
}

sanitize_latex <- function(x) {

  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x)
  x <- gsub("~", "\\\\textasciitilde{}", x)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x <- gsub("−", "-", x)

  x
}

genathum_group_names <- sanitize_latex(genathum_group_names)


result_genathum <- run_realdata_cv_and_selection(
  X = X_genathum,
  y = y_genathum,
  group = group_genathum,
  dataset_name = "GenAtHum",
  group_names = genathum_group_names
)

# ------------------------------------------------------------
# Save all results
# ------------------------------------------------------------

all_cv_selection_results <- list(
  Birthwt = result_birthwt,
  bardet = result_bardet,
  GenAtHum = result_genathum
)

saveRDS(
  all_cv_selection_results,
  file = file.path(selected_group_dir, "RealData_CV_SelectedGroups_Results.rds")
)
