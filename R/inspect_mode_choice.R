library(apollo)
data(apollo_modeChoiceData)
d <- apollo_modeChoiceData
cat("Total rows:", nrow(d), "  IDs:", length(unique(d$ID)), "\n")

cat("\n=== RP vs SP split ===\n")
cat("RP rows:", sum(d$RP), "  SP rows:", sum(d$SP), "\n")

d_sp <- d[d$SP == 1, ]
cat("\n=== SP only ===\n")
cat("Rows:", nrow(d_sp), "  IDs:", length(unique(d_sp$ID)),
    "  Tasks/ID range:", paste(range(table(d_sp$ID)), collapse="-"), "\n")
cat("Availability in SP:\n")
for (a in c("av_car","av_bus","av_air","av_rail")) {
  cat(sprintf("  %s = 1 in %d (%.1f%%)\n", a, sum(d_sp[[a]]), 100*mean(d_sp[[a]])))
}
cat("All 4 available:", sum(d_sp$av_car & d_sp$av_bus & d_sp$av_air & d_sp$av_rail),
    " (", round(100*mean(d_sp$av_car & d_sp$av_bus & d_sp$av_air & d_sp$av_rail), 1), "%)\n")

d_full <- d_sp[d_sp$av_car == 1 & d_sp$av_bus == 1 & d_sp$av_air == 1 & d_sp$av_rail == 1, ]
cat("\n=== SP, all 4 available ===\n")
cat("Rows:", nrow(d_full), "  IDs:", length(unique(d_full$ID)), "\n")
tasks_per_id <- table(d_full$ID)
cat("Tasks/ID summary:\n"); print(summary(as.integer(tasks_per_id)))
cat("Choice distribution:\n"); print(table(d_full$choice))
