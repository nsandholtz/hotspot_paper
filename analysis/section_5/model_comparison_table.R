library(xtable)

acquisition_fits = readRDS("./analysis/section_4/model_output/acquisition_fits.rds")
acquisition_fits_threshold = readRDS("./analysis/section_5/model_output/acquisition_fits_threshold.rds")
acquisition_fits_prospect = readRDS("./analysis/section_5/model_output/acquisition_fits_prospect.rds")

model_comparison_table = data.frame(naive_out_ll = unlist(lapply(acquisition_fits, function(x) x$out_log_lik)),
                                    threshold_out_ll = unlist(lapply(acquisition_fits_threshold, function(x) x$out_log_lik)),
                                    prospect_out_ll = unlist(lapply(acquisition_fits_prospect, function(x) x$out_log_lik)))

model_comparison_table_full = rbind(model_comparison_table, apply(model_comparison_table, 2, sum))
xtable::xtable(model_comparison_table_full)

                                    