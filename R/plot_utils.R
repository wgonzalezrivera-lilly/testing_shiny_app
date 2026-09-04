plot_finemap_pip <- function(
  data,
  position_range,
  pip_threshold,
  credible_set_threshold
) {
  visible <- filter_finemap_data(data, position_range)
  if (nrow(visible) == 0) {
    plot.new()
    text(0.5, 0.5, "No variants in the selected position range.", cex = 1.1)
    return(invisible(NULL))
  }

  ordered <- data[order(-data$PIP, data$POS), , drop = FALSE]
  ordered$cumulative_pip <- cumsum(ordered$PIP) / sum(ordered$PIP)
  credible_set_end <- which(ordered$cumulative_pip >= credible_set_threshold)[1]
  if (is.na(credible_set_end)) credible_set_end <- nrow(ordered)
  ordered$in_credible_set <- seq_len(nrow(ordered)) <= credible_set_end
  visible$in_credible_set <- visible$SNP %in%
    ordered$SNP[ordered$in_credible_set]

  point_color <- ifelse(
    visible$in_credible_set,
    "#C2410C",
    ifelse(visible$PIP >= pip_threshold, "#1D4ED8", "#94A3B8")
  )
  plot(
    visible$POS,
    visible$PIP,
    pch = 21,
    bg = point_color,
    col = "white",
    cex = 1.05,
    xlab = "Position on chromosome 22 (bp)",
    ylab = "Posterior inclusion probability (PIP)",
    main = "Fine-mapping posterior inclusion probabilities",
    yaxs = "i",
    xaxs = "i",
    ylim = c(0, 1)
  )
  abline(h = pip_threshold, col = "#1D4ED8", lty = 2)
  legend(
    "topright",
    legend = c("Credible set", "Above PIP threshold", "Other variants"),
    pt.bg = c("#C2410C", "#1D4ED8", "#94A3B8"),
    pch = 21,
    pt.cex = 1.1,
    bty = "n"
  )
}
