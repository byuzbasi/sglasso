list(schema = "corsiv_mcc_inference_v1", stage = "production", 
    run_id = "corsiv_sz_mcc_b10000_v1", boot_n = 10000L, chunk_size = 100L, 
    seed_base = 1010610L, rng_kind = c("Mersenne-Twister", "Inversion", 
    "Rejection"), threshold = 0x1p-1, positive_rule = "probability >= threshold", 
    stratified = TRUE, paired = TRUE, interval = "percentile", 
    quantile_type = 7L, confidence = 0x1.e666666666666p-1, difference_confidence_bonferroni = 0x1.fae147ae147aep-1, 
    comparison_family_size = 5L, zero_denominator = "mltools_zero_retained_with_flag", 
    model_fits = 0L, methods = c("Logistic SGLASSO", "Logistic SGLASSO (d=0 boundary)", 
    "Logistic Group Elastic Net (adelie)", "Logistic Group Lasso (grpreg)", 
    "Logistic Group MCP (grpreg)", "Logistic Group SCAD (grpreg)"
    ), runtime = list(R = "R version 4.6.0 (2026-04-24)", platform = "aarch64-apple-darwin23", 
        packages = c(mltools = "0.3.5", pROC = "1.19.0.1", digest = "0.6.39", 
        jsonlite = "2.0.0", knitr = "1.51")))
