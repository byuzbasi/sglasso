list(schema = "corsiv_auc_inference_v1", stage = "production", 
    run_id = "corsiv_sz_auc_inference_b10000_v1", boot_n = 10000L, 
    conf_level = 0x1.e666666666666p-1, seed_ci = 1010601L, seed_comparison = 1010602L, 
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"
    ), stratified = TRUE, paired = TRUE, alternative = "two.sided", 
    direction = "<", levels = c(0x0p+0, 0x1p+0), percent = FALSE, 
    adjustment = "holm", multiplicity_family = "five_sglasso_contrasts_separately_per_test_method", 
    interpretation = "exploratory_conditional_on_trained_models_and_class_counts", 
    model_fits = 0L, methods = c("Logistic SGLASSO", "Logistic SGLASSO (d=0 boundary)", 
    "Logistic Group Elastic Net (adelie)", "Logistic Group Lasso (grpreg)", 
    "Logistic Group MCP (grpreg)", "Logistic Group SCAD (grpreg)"
    ), runtime = list(R = "R version 4.6.0 (2026-04-24)", platform = "aarch64-apple-darwin23", 
        packages = c(pROC = "1.19.0.1", digest = "0.6.39", jsonlite = "2.0.0", 
        knitr = "1.51")))
