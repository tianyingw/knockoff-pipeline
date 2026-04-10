plot_manhattan <- function(result_all, out_path, title = 'Manhattan plot (W statistic)') {
  manhattan_data <- data.frame(
    CHR = result_all[,"chr"],
    BP = (result_all[,"start"] + result_all[,"end"]) %/% 2,
    P = result_all[,"W_KS"]
  )
  manhattan_data$SNP <- paste0(manhattan_data$CHR, ":", manhattan_data$BP)
  # manhattan_data$P[!is.finite(manhattan_data$P)] <- 0
  finite_vals <- manhattan_data$P[is.finite(manhattan_data$P)]
  max_finite <- max(finite_vals, na.rm = TRUE)
  manhattan_data$P[!is.finite(manhattan_data$P)] <- max_finite
  
  threshold <- result_all[,"W_Threshold"][1]
  png(file.path(out_path, title), width = 1200, height = 600, res = 150) 
  manhattan(
      manhattan_data,
      chr = 'CHR',
      bp = 'BP',
      p = 'P',
      snp = 'SNP',
      logp = FALSE,
      ylab = 'W statistic',
      col = c("blue4", "orange3"),
      suggestiveline = FALSE,
      genomewideline = FALSE,
      main = "Manhattan Plot",
      cex = 0.7,
      cex.axis = 1.2,
      cex.lab = 1.4,
      cex.main = 1.6,
      las = 1
  )
  if (!is.null(threshold) && is.finite(threshold)) {
    abline(h = threshold, col = "#e41a1c", lwd = 2, lty = 2)
  }
  dev.off()
}