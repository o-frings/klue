# dev/make_figures.R
# Regenerates the three Results figures (fig_reliability, fig_ic, fig_time) from
# the output data files and writes them next to the manuscript, which includes
# them from the repository root. Run from the repo root:  Rscript dev/make_figures.R
# Deterministic: no randomness, no seed needed. It also prints the plotted
# percentages and per-rung mean seconds.
#
# Style: the shared thesis figure style, dev/theme_thesis.R (a verbatim copy of
#  Thesis/figure-style/theme_thesis.R; its STYLE.md has the rules). The text
#  width of the manuscript and of the thesis is 155 mm. Every figure is saved
#  at the width it is included at (fig_ic, fig_reliability at \linewidth =
#  FIG_W_FULL; fig_time at 0.75\linewidth = FIG_W_NARROW), so LaTeX applies no
#  scaling and the type prints at the theme's 8/7 pt. Panel letters are bold
#  lower case (a, b) above the panels.
#
# Palette: ONE shared vector of seven slots, each a (colour, line type, marker)
#  triple, so no series is told apart by colour alone. Slot 1 is black (the
#  clustering starts of this paper), slot 4 grey; slots 2, 3, 5, 6 and 7 are
#  the thesis categorical colours in their fixed order (blue, red, purple,
#  green, orange), which pass the colour-vision checks for all pairs. Colour
#  follows the entity: clustering+ML is slot 1, clustering+EM slot 2, Apollo's
#  default search slot 3 and gmnl slot 4; the three other Apollo arms take
#  slots 5-7.
#
# Reproducibility:
#  - fig_ic            : computed from output/main_blocked.csv
#  - fig_reliability   : (a) output/estimator_blocked.rds (direct ML, per-start
#                        rate by K*), (b) output/stress_replicate_blocked.rds
#                        (share of the 10 seeds per rung within 0.5 LL of the best),
#                        seven arms: klue_ml, klue_em, gmnl_lc, apollo_searchStart
#                        (Apollo's default search), apollo_ss_published (Apollo's
#                        reduced search), apollo_lcEM (Apollo's EM routine),
#                        apollo_clust (Apollo's estimator from the clustering starts)
#  - fig_time          : wall-clock per fit, means over the 10 seeds of each rung
#                        from output/stress_replicate_blocked.rds ($secs), for the
#                        same arms. Apollo's reduced search is the exception: the
#                        ladder's stored times for apollo_ss_published are stale
#                        (dev/stress_replicate_blocked.R), so its times are the
#                        seed-1 runs in output/apollo_smartstart_diag.csv
#                        (column bounded_secs).
#  Points are dodged horizontally (offsets below) so coinciding values stay visible.

if (!file.exists("dev/make_figures.R")) stop("run from the repository root")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("dev/theme_thesis.R")
inp  <- function(f) f           # inputs under output/
outp <- function(f) f           # PDFs in the repository root

## ===== shared design ========================================================
H_PANEL <- 72             # height of the plot row (mm)
H_LEG   <- 26             # height of a one-column legend of seven arms (mm)
LW_LINE <- 0.55           # line width of every series
PT_SIZE <- 1.7            # marker size of every series

col <- c("black", PAL_CAT[["blue"]], PAL_CAT[["red"]], "grey55",
         PAL_CAT[["purple"]], PAL_CAT[["green"]], PAL_CAT[["orange"]])
lty <- c("solid", "dashed", "dotdash", "dotted", "longdash", "twodash", "22")
pch <- c(19, 17, 15, 4, 6, 0, 18)

# Line-and-marker series: one row per point (x, y, series); slots give the
# (colour, line type, marker) of each series label, in legend order
series_plot <- function(d, labels, slots, legend_cols = 1) {
  d$series <- factor(d$series, levels = labels)
  sc <- function(v) setNames(v[slots], labels)
  key <- guide_legend(ncol = legend_cols, byrow = TRUE)
  ggplot(d, aes(x, y, colour = series, linetype = series, shape = series, group = series)) +
    geom_line(linewidth = LW_LINE) +
    geom_point(size = PT_SIZE) +
    scale_colour_manual(values = sc(col), name = NULL) +
    scale_linetype_manual(values = sc(lty), name = NULL) +
    scale_shape_manual(values = sc(pch), name = NULL) +
    guides(colour = key, linetype = key, shape = key) +
    theme_thesis() +
    theme(legend.key.width = unit(24, "pt"), legend.justification = "left")
}
long <- function(m, x, labels) data.frame(x = rep(x, ncol(m)), y = as.vector(m),
                                          series = rep(labels, each = nrow(m)))
tagged <- function(p) p + plot_annotation(tag_levels = "a")

## the seven head-to-head arms, in table order, with their palette slots
ARMS <- data.frame(
  key  = c("klue_ml", "klue_em", "gmnl_lc", "apollo_searchStart",
           "apollo_ss_published", "apollo_lcEM", "apollo_clust"),
  lab  = c("clustering+ML", "clustering+EM", "gmnl", "Apollo's default search",
           "Apollo's reduced search", "Apollo's EM routine",
           "Apollo's estimator from the clustering starts"),
  slot = c(1, 2, 4, 3, 5, 6, 7), stringsAsFactors = FALSE)
DODGE <- seq(-0.15, 0.15, length.out = nrow(ARMS))   # horizontal offsets per arm
RUNGS <- c("moderate", "hard", "very hard")
arms_long <- function(m) {             # rungs x arms matrix, dodged per arm
  d <- long(m, 1:3, ARMS$lab)
  d$x <- d$x + rep(DODGE, each = 3)
  d
}
rung_axis <- scale_x_continuous(breaks = 1:3, labels = RUNGS, limits = c(0.75, 3.25))

## ===== fig_ic : information-criterion accuracy by K* and by separation =====
d <- read.csv(inp("output/main_blocked.csv"))
accBy <- function(var, vals, crit) sapply(vals, function(v)
  100*mean(d[[paste0(crit,"_correct")]][d[[var]]==v]))
K <- 1:5; kap <- c(0.50,0.75,1.00,1.25,1.50)
icK <- cbind(BIC=accBy("true_K",K,"bic"), AIC=accBy("true_K",K,"aic"), ICL=accBy("true_K",K,"icl"))
d2 <- d[d$true_K>=2,]; accBy2 <- function(crit) sapply(kap, function(v) 100*mean(d2[[paste0(crit,"_correct")]][d2$kappa==v]))
icK2 <- cbind(BIC=accBy2("bic"), AIC=accBy2("aic"), ICL=accBy2("icl"))
IC_SLOTS <- c(1, 3, 2)                 # BIC black, AIC red, ICL blue
ic_y <- scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))
p_ic_a <- series_plot(long(icK, K, colnames(icK)), colnames(icK), IC_SLOTS, legend_cols = 3) +
  scale_x_continuous(breaks = K) + ic_y +
  labs(title = "By number of classes", x = expression("Number of classes "*K^"*"),
       y = "Correct enumeration (%)")
p_ic_b <- series_plot(long(icK2, kap, colnames(icK2)), colnames(icK2), IC_SLOTS, legend_cols = 3) +
  scale_x_continuous(breaks = kap) + ic_y +
  labs(title = expression("By separation ("*K^"*" >= 2*")"),
       x = expression("Inter-class separation "*kappa), y = "Correct enumeration (%)")
save_thesis_fig(tagged(p_ic_a + p_ic_b + plot_layout(guides = "collect")) &
                  theme(legend.position = "bottom", legend.justification = "center"),
                outp("fig_ic.pdf"), width_mm = FIG_W_FULL, height_mm = H_PANEL + 10)

## ===== fig_reliability : (a) per-start by K*  (b) reliability by difficulty =====
e <- readRDS(inp("output/estimator_blocked.rds"))
aggE <- aggregate(at_global_rate ~ strategy + K, e[e$estimator=="ml",], mean)
getE <- function(st) sapply(c(3,4,5), function(k) 100*aggE$at_global_rate[aggE$strategy==st & aggE$K==k])
ya <- cbind(getE("clustering"), getE("random_partition"), getE("perturbation"), getE("random"))
s <- readRDS(inp("output/stress_replicate_blocked.rds")); rungs <- c("moderate","hard","very_hard")
stopifnot(all(vapply(s, function(z) all(ARMS$key %in% names(z$gap)), TRUE)))
relB <- sapply(ARMS$key, function(m) sapply(rungs, function(rg){
  g <- sapply(s[sapply(s,function(c)c$rung==rg)], function(c) c$gap[[m]]); g<-g[is.finite(g)]
  stopifnot(length(g) == 10); 100*mean(g >= -0.5) }))
cat("panel (b), % of seeds within 0.5 LL of the best:\n"); print(round(relB))
STARTS <- c("clustering", "random-partition", "pooled-MNL perturbation", "diffuse random")
p_rel_a <- series_plot(long(ya, 3:5, STARTS), STARTS, 1:4) +
  scale_x_continuous(breaks = 3:5) + ic_y +
  labs(title = "Per-start, direct ML", x = expression("Number of classes "*K^"*"),
       y = "Per-start reach-global (%)")
p_rel_b <- series_plot(arms_long(relB), ARMS$lab, ARMS$slot) +
  rung_axis + ic_y +
  labs(title = "Package head-to-head", x = "Data difficulty",
       y = "Seeds within 0.5 LL of best (%)")
save_thesis_fig(tagged(p_rel_a + p_rel_b) & theme(legend.position = "bottom"),
                outp("fig_reliability.pdf"), width_mm = FIG_W_FULL, height_mm = H_PANEL + H_LEG)

## ===== fig_time : estimation time per fit (log) by difficulty =====
ts <- sapply(ARMS$key, function(m) sapply(rungs, function(rg){
  v <- sapply(s[sapply(s,function(c)c$rung==rg)], function(c) c$secs[[m]]); v<-v[is.finite(v)]; mean(v) }))
ss <- read.csv(inp("output/apollo_smartstart_diag.csv"))
stopifnot(identical(ss$rung, rungs))
ts[, "apollo_ss_published"] <- ss$bounded_secs          # seed-1 runs (see header)
cat("fig_time, seconds per fit:\n"); print(round(ts, 1))
stopifnot(max(ts) < 1e4, min(ts) > 1)
YT <- c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000)
p_time <- series_plot(arms_long(ts), ARMS$lab, ARMS$slot) +
  rung_axis +
  scale_y_log10(limits = c(1, 1e4), breaks = YT,
                labels = paste0(format(YT, scientific = FALSE, trim = TRUE), "s")) +
  labs(x = "Data difficulty", y = "Estimation time per fit (s, log scale)")
save_thesis_fig(p_time, outp("fig_time.pdf"), width_mm = FIG_W_NARROW, height_mm = 80 + H_LEG)
cat("regenerated fig_ic.pdf, fig_reliability.pdf, fig_time.pdf at print size in the repository root\n")
