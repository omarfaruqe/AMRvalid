# Draw the simulation study figures from saved results.
#
# Usage (from the package root):
#   Rscript simulation/make_figures.R
# Reads simulation/results/replicates.csv and summary.csv; writes
# fig_optimism.png and fig_coverage.png. run_study.R calls this at the end.
#
# Colours (colourblind-safe, validated for protan, deutan and tritan vision):
#   orange = the naive choice (random split, isolate-level bootstrap)
#   blue   = the population-aware choice (lineage-held-out, lineage bootstrap)

out_dir <- file.path("simulation", "results")
reps <- utils::read.csv(file.path(out_dir, "replicates.csv"), stringsAsFactors = FALSE)
summ <- utils::read.csv(file.path(out_dir, "summary.csv"), stringsAsFactors = FALSE)

col_naive <- "#eb6834"
col_aware <- "#2a78d6"
ink <- "#2b2b2b"
ink_muted <- "#6b6b6b"
grid_col <- "#e6e6e3"

setup_par <- function() {
  graphics::par(mar = c(6.5, 5, 3, 1), col.axis = ink_muted, col.lab = ink, col.main = ink, fg = ink_muted,
                family = "sans", las = 1, mgp = c(3.2, 0.7, 0), tcl = -0.3)
}
short <- function(s) sub(" lineage effect", "", s)
scen_order <- unique(summ$scenario[order(summ$lineage_sd)])

# Figure 1: optimism -------------------------------------------------------
grDevices::png(file.path(out_dir, "fig_optimism.png"), width = 2200, height = 1200, res = 220)
setup_par()
cells <- unique(summ[order(summ$lineage_sd, summ$predictor, summ$design),
                     c("scenario", "predictor", "design")])
pos <- integer(nrow(cells))
gap_after <- 0
for (i in seq_len(nrow(cells))) {
  if (i > 1 && cells$scenario[i] != cells$scenario[i - 1]) gap_after <- gap_after + 1
  if (i > 1 && cells$predictor[i] != cells$predictor[i - 1] && cells$scenario[i] == cells$scenario[i - 1]) {
    gap_after <- gap_after + 0.4
  }
  pos[i] <- i + gap_after
}
groups <- lapply(seq_len(nrow(cells)), function(i) {
  reps$optimism[reps$scenario == cells$scenario[i] & reps$predictor == cells$predictor[i] &
                  reps$design == cells$design[i]]
})
ylim <- range(unlist(groups)) + c(-0.03, 0.08)
graphics::plot.new()
graphics::plot.window(xlim = c(0.4, max(pos) + 0.6), ylim = ylim)
graphics::abline(h = pretty(ylim), col = grid_col, lwd = 1)
graphics::abline(h = 0, col = ink_muted, lty = 2, lwd = 1.2)
graphics::boxplot(groups, at = pos, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.7,
                  col = ifelse(cells$design == "random", col_naive, col_aware),
                  border = ink, medcol = "white", medlwd = 2, whisklty = 1, whiskcol = ink_muted,
                  staplecol = ink_muted)
graphics::axis(2, at = pretty(ylim), lwd = 0, lwd.ticks = 0)
graphics::mtext("Internal AUC minus AUC on new lineages", side = 2, line = 3.2, col = ink, las = 0)
pred_mid <- tapply(pos, paste(cells$scenario, cells$predictor), mean)
pred_lab <- sub(".* ", "", names(pred_mid))
graphics::axis(1, at = pred_mid, labels = pred_lab, lwd = 0, line = -0.4, col.axis = ink)
scen_mid <- tapply(pos, cells$scenario, mean)[scen_order]
graphics::mtext(paste(short(scen_order), "lineage effect"), side = 1, at = scen_mid, line = 2.6, col = ink, font = 2)
graphics::legend("topleft", fill = c(col_naive, col_aware), border = NA, bty = "n", cex = 0.85,
                 legend = c("random split", "lineage-held-out split"), text.col = ink)
graphics::title(main = "Random splits overstate the lineage-memorising predictor", adj = 0, cex.main = 1,
                font.main = 1)
grDevices::dev.off()

# Figure 2: coverage -------------------------------------------------------
grDevices::png(file.path(out_dir, "fig_coverage.png"), width = 2200, height = 1200, res = 220)
setup_par()
lin <- summ[summ$design == "lineage", ]
lin <- lin[order(lin$lineage_sd, lin$predictor), ]
cov <- rbind(lin$coverage_isolate_bootstrap, lin$coverage_lineage_bootstrap)
space <- c(0.08, 0.6)
mids <- graphics::barplot(cov, beside = TRUE, ylim = c(0, 1.12), axes = FALSE, names.arg = rep("", ncol(cov)),
                          col = c(col_naive, col_aware), border = NA, space = space)
graphics::abline(h = seq(0, 1, 0.2), col = grid_col, lwd = 1)
graphics::barplot(cov, beside = TRUE, add = TRUE, axes = FALSE, names.arg = rep("", ncol(cov)),
                  col = c(col_naive, col_aware), border = NA, space = space)
graphics::abline(h = 0.95, col = ink_muted, lty = 2, lwd = 1.2)
graphics::text(max(mids) + 0.9, 0.95, "95% target", pos = 3, cex = 0.75, col = ink_muted, xpd = TRUE)
graphics::text(mids, cov - 0.01, sprintf("%.0f%%", 100 * cov), pos = 1, cex = 0.7, col = "white", font = 2)
graphics::axis(2, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20), "%"), lwd = 0)
graphics::mtext("Intervals containing the true AUC", side = 2, line = 3.2, col = ink, las = 0)
graphics::axis(1, at = colMeans(mids), labels = lin$predictor, lwd = 0, line = -0.6, col.axis = ink)
scen_mid <- tapply(colMeans(mids), lin$scenario, mean)[scen_order]
graphics::mtext(paste(short(scen_order), "lineage effect"), side = 1, at = scen_mid, line = 2.2, col = ink, font = 2)
graphics::legend("topleft", fill = c(col_naive, col_aware), border = NA, bty = "n", cex = 0.85, horiz = TRUE,
                 legend = c("isolate-level bootstrap", "lineage-level bootstrap"), text.col = ink)
graphics::title(main = "Lineage-level intervals keep close to their nominal coverage (lineage-held-out design)",
                adj = 0, cex.main = 1, font.main = 1)
grDevices::dev.off()

cat("Figures written to", out_dir, "\n")
