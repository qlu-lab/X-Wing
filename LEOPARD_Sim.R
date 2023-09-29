rm(list = ls())

suppressMessages(require(data.table))
suppressMessages(require(BEDMatrix))
suppressMessages(require(optparse))

options(stringsAsFactors = F)
option_list = list(
  make_option("--sumstats", action = "store", default = NA, type = "character"),
  make_option("--n_gwas", action = "store", default = NA, type = "numeric"),
  make_option("--train_prop", action = "store", default = NA, type = "numeric"),
  make_option("--seed", action = "store", default = NA, type = "numeric"),
  make_option("--ref_prefix", action = "store", default = NA, type = "character"),
  make_option("--rep", action = "store", default = NA, type = "numeric"),
  make_option("--out_prefix", action = "store", default = NA, type = "character")
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
# Print the input

cat("Options in effect: \n")
cat("Rscript LEOPARD_Sim.R \\ \n")
cat(paste0("--sumstats ", opt$sumstats, " \\ \n"))
cat(paste0("--n_gwas ", opt$n_gwas, " \\ \n"))
cat(paste0("--train_prop ", opt$train_prop, " \\ \n"))
cat(paste0("--ref_prefix ", opt$ref_prefix, " \\ \n"))
if (!is.na(opt$seed)){
  cat(paste0("--seed ", opt$seed, " \\ \n"))
}
cat(paste0("--rep ", opt$rep, " \\ \n"))
cat(paste0("--out_prefix ", opt$out_prefix, " \n \n"))


# Input
sumstats_path <- opt$sumstats
N <- as.numeric(opt$n_gwas)
train_prop <- as.numeric(opt$train_prop)
seed <- as.numeric(opt$seed)
ref <-  opt$ref_prefix
rep <-  opt$rep
out <- opt$out_prefix

if(is.na(seed)){
  seed <- sample(1:(42*10^5), 1)
}

# Data preparation
cat("#### Begin Sampling GWAS summary statistics for traning and validing set ####\n")
cat(paste0("--- Propotion of traning set is ", train_prop," \n"))
cat("--- Reading the data \n")
ref_bed <- suppressMessages(BEDMatrix(paste0(ref, ".bed")))
ref_bim <- fread(paste0(ref, ".bim"))
sumstats <- fread(sumstats_path)
N_train <- floor(N*train_prop)
N_valid <- N - N_train

# Take the overlap between reference panel and GWAS sumstats
snps_ovp <- intersect(sumstats$SNP, ref_bim$V2)
ref_bim_ovp <- ref_bim[match(snps_ovp, ref_bim$V2), ]
sumstats_ovp <- sumstats[match(snps_ovp, sumstats$SNP), ]
# standardized the genotype
ref_bed_ovp <- scale(ref_bed[, match(snps_ovp, ref_bim$V2)])

# Align the A1/A2 for the sumstats
## Remove the SNPs with umatched A1/A2
# Here I should also add A1/A2
ind_right <- which(sumstats_ovp$A1 == ref_bim_ovp$V5 | sumstats_ovp$A1 == ref_bim_ovp$V6) # 0, need to change it in order to make it works better
ref_bim_ovp <- ref_bim_ovp[ind_right, ]
sumstats_ovp <- sumstats_ovp[ind_right, ]
ref_bed_ovp <- ref_bed_ovp[, ind_right]

## If A1 for sumstats = A2 in reference panel, filp the sign of beta in sumstats
ind_mis <- which(sumstats_ovp$A1 == ref_bim_ovp$V6) 
sumstats_ovp$BETA[ind_mis] <- -sumstats_ovp$BETA[ind_mis]
sumstats_ovp$A1[ind_mis] <- ref_bim_ovp$V5[ind_mis]
sumstats_ovp$A2[ind_mis] <- ref_bim_ovp$V6[ind_mis]

# xTy in the full data
Z <- qnorm(sumstats_ovp$P/2, lower.tail=FALSE )
sumstats_ovp$xTy <- Z * sqrt(N) * sign(sumstats_ovp$BETA)

# Sample the training data
cat("--- Initializing sampling \n")
for (i in 1:rep){
  set.seed(seed * i)
  cat(paste0("--- Sampling for replication ", i," \n"))
  cat(paste0("--- Current seed is ", seed * i," \n"))
  mu <- (sumstats_ovp$xTy/N)
  Sigma <- sqrt((N_valid/N_train/N/nrow(ref_bed_ovp))) * t(as.matrix(ref_bed_ovp)) %*% rnorm(nrow(ref_bed_ovp)) # Here is where the randomness comes from

  # training data sumstats
  rho_train <- mu + Sigma
  Z_train <- sqrt(N_train * rho_train ^2 / (1 - rho_train^2)) * sign(rho_train)
  P_train <- 2*pnorm(abs(Z_train), 0, 1, lower.tail = F)
  out_train <- sumstats_ovp[, 1:7]
  out_train$P <- P_train
  out_train$BETA <- sign(rho_train) * abs(qnorm(out_train$P/2)) / sqrt(N_train)

  # validation data sumstats
  rho_valid <- (sumstats_ovp$xTy - rho_train * N_train) / N_valid
  Z_valid <- sqrt(N_valid * rho_valid ^2 / (1 - rho_valid^2)) * sign(rho_valid)
  P_valid <- 2*pnorm(abs(Z_valid), 0, 1, lower.tail = F)
  out_valid <- sumstats_ovp[, 1:7]
  out_valid$P <- P_valid
  out_valid$BETA <- sign(rho_valid) * abs(qnorm(out_valid$P/2)) / sqrt(N_valid)

  # Remove NA from data
  ind_noNA <- intersect(which(!(is.na(out_train$BETA) )), which(!(is.na(out_valid$BETA) )))

  fwrite(out_train[ind_noNA, ], paste0(out, "_rep", i, "_train.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
  fwrite(out_valid[ind_noNA, ], paste0(out, "_rep", i, "_valid.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
  fwrite(data.frame(N_train, N_valid), paste0(out, "_rep", i,"_train_valid_N.txt"), sep = "\t", quote = F, row.names = F, col.names = T) 
}
cat("#### Finsh LEOPARD GWAS summary statistics sampling! ####\n")
cat("\n")
