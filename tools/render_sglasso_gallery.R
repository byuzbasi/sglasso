#!/usr/bin/env Rscript
# Documentation only: one tiny Gaussian fit and verified, already fitted ROC curves.
# Usage: Rscript --vanilla tools/render_sglasso_gallery.R \
#   --external=/path/to/completed/external/run --output=/path/to/new/gallery
# Optional: --library=/path/to/an/already/installed/R/library
args <- commandArgs(trailingOnly = TRUE)
option <- function(key) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit) != 1L) stop("Provide exactly one --", key, "= argument.")
  sub(paste0("^--", key, "="), "", hit)
}
if (any(grepl("^--library=", args))) {
  .libPaths(c(normalizePath(option("library"), mustWork = TRUE), .libPaths()))
}
for (p in c("sglasso", "digest", "pROC", "svglite")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing installed package: ", p)
}
stopifnot(utils::packageVersion("sglasso") >= "1.2.0")
external <- normalizePath(option("external"), mustWork = TRUE)
out <- option("output")
if (file.exists(out)) stop("Output already exists; choose a new versioned directory.")
hash <- function(f) digest::digest(file = f, algo = "sha256", serialize = FALSE)
manifest_file <- file.path(external, "CHECKPOINT_MANIFEST.rds")
manifest <- readRDS(manifest_file)
stopifnot(all(c("file", "bytes", "sha256") %in% names(manifest)))
stopifnot(!anyDuplicated(manifest$file), !any(grepl("(^/|(^|/)\\.\\.(/|$))", manifest$file)))
paths <- file.path(external, manifest$file)
stopifnot(all(file.exists(paths)), all(file.info(paths)$size == manifest$bytes))
stopifnot(identical(unname(vapply(paths, hash, character(1))), manifest$sha256))
complete <- readRDS(file.path(external, "COMPLETE.rds"))
spec <- readRDS(file.path(external, "specification.rds"))
stopifnot(complete$validated_units == 6L, complete$method_rows == 6L,
          identical(complete$scientific_signature, spec$scientific_signature),
          identical(complete$manifest_sha256, hash(file.path(external, "final/manifest.rds"))))
pred <- readRDS(file.path(external, "final/predictions.rds"))
summary <- readRDS(file.path(external, "final/response_summary.rds"))
stopifnot(nrow(pred) == 4050L, nrow(summary) == 6L,
          length(unique(pred$method)) == 6L,
          all(is.finite(pred$probability)), all(pred$probability >= 0 & pred$probability <= 1))
original <- c("Logistic SGLASSO", "Logistic Group Elastic Net (adelie)",
              "Logistic Group Lasso (grpreg)", "Logistic Group MCP (grpreg)",
              "Logistic Group SCAD (grpreg)")
labels <- c("SGLASSO", "Group ENET", "Group Lasso", "Group MCP", "Group SCAD")
palette <- c("#087F8C", "#7963B2", "#C06E32", "#395B85", "#A4446B")
by_method <- split(pred, pred$method)
reference <- by_method[[original[1]]]
stopifnot(nrow(reference) == 675L, !anyDuplicated(reference$sample_id),
          sum(reference$y == 1) == 353L, sum(reference$y == 0) == 322L)
for (z in by_method) {
  stopifnot(identical(z$sample_id, reference$sample_id), identical(z$y, reference$y))
}
rocs <- lapply(original, function(m) {
  z <- by_method[[m]]
  r <- pROC::roc(z$y, z$probability, levels = c(0, 1), direction = "<", quiet = TRUE)
  a <- summary$auc[match(m, summary$method)]
  stopifnot(is.finite(a), abs(as.numeric(pROC::auc(r)) - a) < 1e-12)
  r
})
dir.create(out, recursive = TRUE)
draw <- function(name, fun, width = 9, height = 5.5) {
  svglite::svglite(file.path(out, name), width = width, height = height, bg = "white")
  tryCatch(fun(), finally = grDevices::dev.off())
}
set.seed(19)
x <- matrix(rnorm(80 * 6), 80, 6)
group <- rep(1:3, each = 2)
y <- x[, 1] - x[, 3] + rnorm(80)
lambda <- c(0.3, 0.1, 0.03)
fit <- sglasso::sglasso(x, y, group, family = "gaussian", lambda = lambda,
                       d = c(0, 0.5), alpha = 0.5, screen = "none")
coefs <- vapply(lambda, function(l) stats::coef(fit, lambda = l, d = 0.5)[-1], numeric(6))
stopifnot(all(is.finite(coefs)))
draw("sglasso-quickstart-path.svg", function() {
  par(mar = c(4.3, 4.4, 3.5, 1), fg = "#203448", col.axis = "#203448", las = 1)
  matplot(lambda, t(coefs), type = "n", log = "x", xlab = "Penalty strength (lambda)",
          ylab = "Coefficient", main = "Grouped coefficient paths", bty = "l")
  abline(h = 0, col = "#D9E5EC", lty = 2)
  for (j in 1:6) lines(lambda, coefs[j, ], col = palette[group[j]], lwd = 2,
                        lty = if (j %% 2) 1 else 2, type = "b", pch = if (j %% 2) 16 else 1)
  legend("bottomleft", legend = paste("Group", 1:3), col = palette[1:3], lwd = 2,
         bty = "n", horiz = TRUE, cex = 0.85)
  mtext("Illustrative Gaussian example | n = 80, p = 6 | alpha = 0.5, d = 0.5",
        side = 3, line = 0.2, cex = 0.8, col = "#53697C")
})
draw("sglasso-external-roc.svg", function() {
  par(mar = c(4.3, 4.4, 3.5, 1), fg = "#203448", col.axis = "#203448", las = 1)
  plot(c(0, 1), c(0, 1), type = "n", xlab = "False-positive rate (1 - specificity)",
       ylab = "True-positive rate (sensitivity)", main = "Independent external-test ROC", bty = "l")
  abline(0, 1, col = "#B7C7D1", lty = 2)
  for (i in seq_along(rocs)) lines(1 - rocs[[i]]$specificities, rocs[[i]]$sensitivities,
                                  col = palette[i], lwd = 2, lty = i)
  legend("bottomright", legend = sprintf("%s  |  AUC %.3f", labels,
         vapply(rocs, function(r) as.numeric(pROC::auc(r)), numeric(1))),
         col = palette, lwd = 2, lty = 1:5, bty = "n", cex = 0.85)
  mtext("CoRSIVSZ | 675 independent-test participants | fixed fitted predictions",
        side = 3, line = 0.2, cex = 0.8, col = "#53697C")
})
saveRDS(list(schema = "sglasso_gallery_v1", seed = 19L, tiny_model_fits = 1L,
             external_model_fits = 0L, bootstrap_draws = 0L,
             external_signature = complete$scientific_signature,
             checks = c(manifest = TRUE, alignment = TRUE, auc_reconstruction = TRUE,
                        equal_roc_line_width = TRUE, finite_coefficients = TRUE),
             package_versions = vapply(c("sglasso", "digest", "pROC", "svglite"),
               function(p) as.character(utils::packageVersion(p)), character(1)),
             r_version = R.version.string), file.path(out, "VALIDATED.rds"))
files <- list.files(out, full.names = TRUE)
utils::write.csv(data.frame(file = basename(files), bytes = file.info(files)$size,
  sha256 = vapply(files, hash, character(1))), file.path(out, "MANIFEST.csv"), row.names = FALSE)
cat("VERIFIED: completed external manifest, 675 aligned participants, five AUC reconstructions.\n")
cat("Tiny illustrative model fits: 1; external refits and bootstrap draws: 0.\nOutput:", out, "\n")
