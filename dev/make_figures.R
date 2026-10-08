# dev/make_figures.R
# Regenerates the three Results figures (fig_reliability, fig_ic, fig_time) from
# the output data files and writes them next to the manuscript, which includes
# them from the repository root. Run from the repo root:  Rscript dev/make_figures.R
# Deterministic: no randomness, no seed needed. It also prints the plotted
# percentages and per-rung mean seconds.
#
# Sizing convention (constant effective font size across figures):
#  The manuscript is a 10pt article on letter paper with 20mm margins, so
#  \linewidth = 8.5in - 2*(20/25.4)in = 6.925in. Every figure is saved at the
#  exact physical width it is included at (fig_ic, fig_reliability at
#  \linewidth; fig_time at 0.62\linewidth), so LaTeX applies NO scaling and the
#  shared pointsize = 9 prints at a true 9pt in every figure. All axis/label/
#  title text uses cex = 1 (9pt; titles bold); legends use the single shared
#  LEG_CEX. Shared palette, line widths, grid style, and margins throughout.
#  Heights are free: fig_reliability and fig_time carry their legends in a strip
#  below the plots, because seven series leave no empty zone inside the panels.
#
# Palette: ONE shared vector of seven slots, each a (colour, line type, marker)
#  triple, so no series is told apart by colour alone. Slots 1-4 are the
#  original palette (black, #2c7fb8, #d95f0e, grey55); slots 5-7 extend it with
#  #4a3aa7, #008300, #c51b7d. The chromatic slots in the order 2, 3, 5, 6, 7
#  pass the dataviz validator on adjacent pairs (light, white surface; worst
#  CVD dE 11.3, worst normal-vision dE 29.2; all-pairs normal vision >= 17.8);
#  pairs it cannot separate under colour-vision deficiency are separated by
#  line type and marker. Colour follows the entity: clustering+ML is slot 1,
#  clustering+EM slot 2, Apollo's default search slot 3 and gmnl slot 4, as in
#  the previous figure; the three new Apollo arms take slots 5-7.
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
inp  <- function(f) f           # inputs under output/
outp <- function(f) f           # PDFs in the repository root

## ===== shared design system ================================================
LW      <- 6.925          # \linewidth in inches (letter, 20mm margins)
H_PANEL <- 3.30           # height of the plot row of the two-panel figures
H_LEG   <- 1.05           # height of the legend strip below fig_reliability / fig_time
PS      <- 9              # pointsize: true printed size, no LaTeX scaling
LEG_CEX <- 0.85           # single legend size shared by all figures
GRID    <- function(at) abline(h = at, col = "grey92", lwd = 0.6)

col <- c("black", "#2c7fb8", "#d95f0e", "grey55", "#4a3aa7", "#008300", "#c51b7d")
lty <- c("solid", "dashed", "dotdash", "dotted", "longdash", "twodash", "22")  # = 1, 2, 4, 3, 5, 6 and short dashes
pch <- c(19, 17, 15, 4, 6, 0, 18)

new_fig <- function(file, width, height) {
  pdf(outp(file), width = width, height = height, pointsize = PS, useDingbats = FALSE)
  par(mgp = c(2.4, 0.7, 0), las = 1, cex.main = 1)
}
legend_strip <- function(labels, slots, ncol = 1) {
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("top", labels, col = col[slots], lty = lty[slots], pch = pch[slots], lwd = 2,
         bty = "n", cex = LEG_CEX, ncol = ncol, seg.len = 3)
}

## the seven head-to-head arms, in table order, with their palette slots
ARMS <- data.frame(
  key  = c("klue_ml", "klue_em", "gmnl_lc", "apollo_searchStart",
           "apollo_ss_published", "apollo_lcEM", "apollo_clust"),
  lab  = c("clustering+ML", "clustering+EM", "gmnl", "Apollo's default search",
           "Apollo's reduced search", "Apollo's EM routine",
           "Apollo's estimator from the clustering starts"),
  slot = c(1, 2, 4, 3, 5, 6, 7), stringsAsFactors = FALSE)
DODGE <- seq(-0.15, 0.15, length.out = nrow(ARMS))   # horizontal offsets per arm

## ===== fig_ic : information-criterion accuracy by K* and by separation =====
d <- read.csv(inp("output/main_blocked.csv"))
accBy <- function(var, vals, crit) sapply(vals, function(v)
  100*mean(d[[paste0(crit,"_correct")]][d[[var]]==v]))
K <- 1:5; kap <- c(0.50,0.75,1.00,1.25,1.50)
icK <- cbind(BIC=accBy("true_K",K,"bic"), AIC=accBy("true_K",K,"aic"), ICL=accBy("true_K",K,"icl"))
d2 <- d[d$true_K>=2,]; accBy2 <- function(crit) sapply(kap, function(v) 100*mean(d2[[paste0(crit,"_correct")]][d2$kappa==v]))
icK2 <- cbind(BIC=accBy2("bic"), AIC=accBy2("aic"), ICL=accBy2("icl"))
cic <- c("black","#d95f0e","#2c7fb8"); lic <- c(1,4,2); pic <- c(19,15,17)
new_fig("fig_ic.pdf", width = LW, height = H_PANEL)
par(mfrow = c(1,2), mar = c(4.0, 4.0, 2.2, 0.8))
matplot(K, icK, type="n", xaxt="n", ylim=c(0,100),
        xlab=expression("number of classes "*K^"*"), ylab="correct enumeration (%)", main="(a) by number of classes")
GRID(seq(0, 100, 20))
matlines(K, icK, type="b", lwd=2, col=cic, lty=lic, pch=pic)
axis(1,at=K); legend("right", c("BIC","AIC","ICL"), col=cic, lty=lic, pch=pic, lwd=2, bty="n", cex=LEG_CEX)
matplot(kap, icK2, type="n", xaxt="n", ylim=c(0,100),
        xlab=expression("inter-class separation "*kappa), ylab="correct enumeration (%)",
        main=expression("(b) by separation ("*K^"*">=2*")"))
GRID(seq(0, 100, 20))
matlines(kap, icK2, type="b", lwd=2, col=cic, lty=lic, pch=pic)
axis(1,at=kap); legend("right", c("BIC","AIC","ICL"), col=cic, lty=lic, pch=pic, lwd=2, bty="n", cex=LEG_CEX)
dev.off()

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
new_fig("fig_reliability.pdf", width = LW, height = H_PANEL + H_LEG)
layout(rbind(c(1, 2), c(3, 4)), heights = c(H_PANEL, H_LEG))
par(mar = c(4.0, 4.0, 2.2, 0.8))
matplot(3:5, ya, type="n", xaxt="n", ylim=c(0,100),
        xlab=expression("number of classes "*K^"*"), ylab="per-start reach-global (%)", main="(a) Per-start, direct ML")
GRID(seq(0, 100, 20))
matlines(3:5, ya, type="b", lwd=2, col=col[1:4], lty=lty[1:4], pch=pch[1:4])
axis(1,at=3:5)
par(mar = c(4.0, 4.0, 2.2, 0.8))
plot(NA, xlim = c(0.75, 3.25), ylim = c(0, 100), xaxt = "n",
     xlab = "data difficulty", ylab = "seeds within 0.5 LL of best (%)", main = "(b) Package head-to-head")
GRID(seq(0, 100, 20))
for (i in seq_len(nrow(ARMS))) {
  k <- ARMS$slot[i]
  lines(1:3 + DODGE[i], relB[, i], type = "b", lwd = 2, col = col[k], lty = lty[k], pch = pch[k])
}
axis(1, at = 1:3, labels = gsub("_", " ", rungs))
legend_strip(c("clustering", "random-partition", "pooled-MNL perturbation", "diffuse random"), 1:4)
legend_strip(ARMS$lab, ARMS$slot)
dev.off()

## ===== fig_time : estimation time per fit (log) by difficulty =====
ts <- sapply(ARMS$key, function(m) sapply(rungs, function(rg){
  v <- sapply(s[sapply(s,function(c)c$rung==rg)], function(c) c$secs[[m]]); v<-v[is.finite(v)]; mean(v) }))
ss <- read.csv(inp("output/apollo_smartstart_diag.csv"))
stopifnot(identical(ss$rung, rungs))
ts[, "apollo_ss_published"] <- ss$bounded_secs          # seed-1 runs (see header)
cat("fig_time, seconds per fit:\n"); print(round(ts, 1))
stopifnot(max(ts) < 1e4, min(ts) > 1)
YT <- c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000)
new_fig("fig_time.pdf", width = 0.62 * LW, height = 3.20 + H_LEG)
layout(matrix(1:2, ncol = 1), heights = c(3.20, H_LEG))
par(mar = c(4.0, 5.3, 1.2, 0.8))
plot(NA, log = "y", xaxt = "n", yaxt = "n", xlab = "data difficulty", ylab = "",
     ylim = c(1, 1e4), xlim = c(0.75, 3.25))
GRID(YT)
for (i in seq_len(nrow(ARMS))) {
  k <- ARMS$slot[i]
  lines(1:3 + DODGE[i], ts[, i], type = "b", lwd = 2, col = col[k], lty = lty[k], pch = pch[k])
}
axis(1, at = 1:3, labels = gsub("_", " ", rungs))
axis(2, at = YT, labels = paste0(format(YT, scientific = FALSE, trim = TRUE), "s"))
title(ylab = "estimation time per fit (s, log scale)", line = 4.1)  # clear of the wide "10000s" tick label
legend_strip(ARMS$lab, ARMS$slot)
dev.off()
cat("regenerated fig_ic.pdf, fig_reliability.pdf, fig_time.pdf at print size (9pt effective) in the repository root\n")
