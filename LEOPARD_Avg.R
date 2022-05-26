rm(list = ls())

suppressMessages(require(data.table))
suppressMessages(require(BEDMatrix))
suppressMessages(require(optparse))

options(stringsAsFactors = F)
option_list = list(
  make_option("--weights_prefix", action = "store", default = NA, type = "character"),
  make_option("--rep", action = "store", default = NA, type = "numeric"),
  make_option("--out", action = "store", default = NA, type = "character")
)
opt = parse_args(OptionParser(option_list=option_list))


cat("********************************************************************* \n")
cat("* Cross-population Weighting (X-Wing) \n")
cat("* Version 1.0.0 \n")
cat("* Step3: Linear combination for polygenic score using gwas summary statistics data (LEOPARD) \n")
cat("* (C) Jiacheng Miao and Hanmin Guo \n")
cat("* University of Wisconsin–Madison and Tsinghua University \n")
cat("* https://github.com/qlu-lab/X-Wing \n")
cat("* GNU General Public License v3 \n")
cat("********************************************************************* \n \n")

cat("Options in effect: \n")
cat("Rscript LEOPARD_Avg.R \\ \n")
cat(paste0("--weights_prefix ", opt$weights_prefix, " \\ \n"))
cat(paste0("--rep ", opt$rep, " \\ \n"))
cat(paste0("--out ", opt$out, " \n \n"))

# Input
weights_prefix <- opt$weights_prefix
rep <-  opt$rep
out <- opt$out


cat("#### Read the linear combination weights for each replications ####\n")
cat("--- Calculating the averaged linear combination weights \n")


cal_avg_rel_weights <- function(path){
  weights_file <- fread(path)
  weights_non_0 <-  ifelse(weights_file$Weights < 0, 0, weights_file$Weights)
  rel_weights <- weights_non_0/(sum(weights_non_0))
  return(rel_weights)
}

rel_weights_list <- lapply(c(1:rep), function(i) cal_avg_rel_weights(paste0(weights_prefix, i, ".txt")))
avg_weights <- colMeans(as.data.frame(do.call(rbind, rel_weights_list)))

df_out <- data.frame(Path = fread(paste0(weights_prefix, "1.txt"))$Path, Weights = avg_weights)
rownames(df_out) <- NULL
fwrite(df_out, out,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t")
cat("--- The LEOPARD averaged linear combination weights are:  \n")
print(df_out)
cat("\n")
cat("#### Finsh LEOPARD calculation for the averaged linear combination weights! ####\n")
cat("\n")