rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(BEDMatrix))
options(stringsAsFactors = F)
option_list = list(
  make_option("--beta_file ", action = "store", default = NA, type = "character"),
  make_option("--valid_file", action = "store", default = NA, type = "character"),
  make_option("--n_valid", action = "store", default = NA, type = "numeric"),
  make_option("--ref_prefix", action = "store", default = NA, type = "character"),
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
# Print the input

cat("Options in effect: \n")
cat("Rscript LEOPARD_Weights.R \\ \n")
cat(paste0("--beta_file ", opt$beta_file, " \\ \n"))
cat(paste0("--valid_file ", opt$valid_file, " \\ \n"))
cat(paste0("--n_valid ", opt$n_valid, " \\ \n"))
cat(paste0("--ref_prefix ", opt$ref_prefix, " \\ \n"))
cat(paste0("--out ", opt$out, " \n \n"))

cat("### Begin estimating linear combination weights! ###\n")

beta_file <- unlist(strsplit(opt$beta_file, ","))
valid_file <- opt$valid_file
n_valid <- opt$n_valid
ref_prefix <- opt$ref_prefix
out <- opt$out

# Read the data
cat("--- Reading the data \n")
beta <- c()
snps_union <- c()
for (i in 1:length(beta_file)){
  beta[[i]] <- fread(beta_file[i])
  snps_union <- c(snps_union, beta[[i]]$SNP)
}
snps_union <- unique(snps_union)

rho.valid <- fread(valid_file)

ref_frq <- fread(paste0(ref_prefix, ".frq"))
ref_bed <- suppressMessages(BEDMatrix(paste0(ref_prefix, ".bed")))
ref_bim <- fread(paste0(ref_prefix, ".bim"))
colnames(ref_bim) <- c("CHR", "SNP", "POS", "BP", "A1", "A2")
ref_fam <- fread(paste0(ref_prefix, ".fam"))
n_ref <- nrow(ref_bed)

# QC for input data
cat("--- QC for input data \n")
## Find overlapped SNPs acorss all data
snps_ovp <- intersect(ref_bim$SNP, snps_union)
snps_ovp <- intersect(rho.valid$SNP, snps_ovp)
plink_ind <- match(snps_ovp, ref_bim$SNP)
ref_bim <- ref_bim[plink_ind, ]
ref_frq <- ref_frq[plink_ind, ]
ref_bed <- ref_bed[, plink_ind]
rho.valid <- rho.valid[match(snps_ovp, rho.valid$SNP), ]
for (i in 1:length(beta_file)){
  beta[[i]] <- beta[[i]][match(snps_ovp, beta[[i]]$SNP), ]
}
cat("--- ", length(snps_ovp)," SNPs are matched in all input files \n",sep = "")

## Remove ambiguous SNPs
# replace T with A, replace G with C; A=1, C=2
ref_bim$A1[ref_bim$A1 == "T"] <- "A"
ref_bim$A1[ref_bim$A1 == "G"] <- "C"
ref_bim$A2[ref_bim$A2 == "T"] <- "A"
ref_bim$A2[ref_bim$A2 == "G"] <- "C"
ref_A1 <- ifelse(ref_bim$A1=="A",1,2)
ref_A2 <- ifelse(ref_bim$A2=="A",1,2)

rho.valid$A1[rho.valid$A1 == "T"] <- "A"
rho.valid$A1[rho.valid$A1 == "G"] <- "C"
rho.valid$A2[rho.valid$A2 == "T"] <- "A"
rho.valid$A2[rho.valid$A2 == "G"] <- "C"
rho.valid_A1 <- ifelse(rho.valid$A1=="A",1,2)
rho.valid_A2 <- ifelse(rho.valid$A2=="A",1,2)

snps_keep <- c()
for (i in 1:length(beta_file)){
  beta[[i]]$A1[beta[[i]]$A1 == "T"] <- "A"
  beta[[i]]$A1[beta[[i]]$A1 == "G"] <- "C"
  beta[[i]]$A2[beta[[i]]$A2 == "T"] <- "A"
  beta[[i]]$A2[beta[[i]]$A2 == "G"] <- "C"
  beta_A1 <- ifelse(beta[[i]]$A1=="A",1,2)
  beta_A2 <- ifelse(beta[[i]]$A2=="A",1,2)
  snps_rm <- (ref_A1+ref_A2)!=(rho.valid_A1+rho.valid_A2) | (ref_A1+ref_A2)!=(beta_A1+beta_A2) | (rho.valid_A1+rho.valid_A2)!=(beta_A1+beta_A2)
  snps_keep <- c(snps_keep, snps_ovp[!(snps_rm)])
}
snps_keep <- unique(snps_keep)
snps_keep <- snps_keep[!is.na(snps_keep)]
cat("--- ", length(snps_keep)," SNPs left after removing ", length(snps_ovp) - length(snps_keep)," ambiguous SNPs \n",sep = "")

# 
plink_ind2 <- match(snps_keep, ref_bim$SNP)
ref_bim <- ref_bim[plink_ind2, ]
ref_frq <- ref_frq[plink_ind2, ]
ref_bed <- ref_bed[, plink_ind2]
rho.valid <- rho.valid[match(snps_keep, rho.valid$SNP), ]
for (i in 1:length(beta_file)){
  beta[[i]] <- beta[[i]][match(snps_keep, beta[[i]]$SNP), ]
}


# 
rho.valid$Z <- abs(qnorm(rho.valid$P/2)) * sign(rho.valid$BETA)
rho.valid$rho <- sqrt((rho.valid$Z)^2 * n_valid) * sign(rho.valid$BETA)

# Calculate the SNP effects w.r.t standardized allele scale
BETA_STD <- data.frame(beta[[i]]$SNP)
for (i in 1:length(beta_file)){
  beta[[i]]$FREQ <- ifelse(beta[[i]]$A1 == ref_frq$A1, ref_frq$MAF, 1- ref_frq$MAF)
  beta[[i]]$BETA_STD <-  beta[[i]]$BETA * sqrt(2 * (1 - beta[[i]]$FREQ) * beta[[i]]$FREQ)
  BETA_STD <- cbind(BETA_STD, beta[[i]]$BETA_STD)
}

# Calcualte the PRS in reference panel using the input PRS weights
prs <- ref_fam[,2]
for (i in 1:length(beta_file)){
  bed_matrix <- data.matrix(ref_bed)
  bed_matrix[is.na(bed_matrix)] <- 0
  # Match the A1/A2
  beta[[i]]$BETA_prs <- ifelse(beta[[i]]$A1 == ref_bim$A1, beta[[i]]$BETA, -beta[[i]]$BETA)
  beta[[i]]$BETA_prs[is.na(beta[[i]]$BETA_prs)] <-  0

  effect.martrix <- data.matrix(beta[[i]]$BETA_prs)
  prs_value <-  bed_matrix %*% effect.martrix
  prs_value <- prs_value - mean(prs_value)
  prs <- cbind(prs, data.frame(prs_value))
}

# Calculate the linear combination weights
cat("--- Calculating the linear combination weights \n")
xbeta <- as.matrix(data.frame(prs[, -1]))

Q <- t(xbeta) %*% xbeta
Q.inv <- solve(Q)

# Make sure the sign of xTy is correct
xTy <- as.matrix(rho.valid$rho)
xTy[is.na(xTy)] <- 0
beta.tilde <- as.matrix(data.frame(BETA_STD[, -1]))
beta.tilde[is.na(beta.tilde)] <- 0
second <- t(beta.tilde) %*% xTy

weights <- Q.inv %*% second

weights_save <- data.frame(beta_file = beta_file, weights = weights)
colnames(weights_save)<- NULL

# print(weights_save)
df_out <- data.frame(weights_save)
colnames(df_out) <- c("Path", "Weights")
rownames(df_out) <- NULL
fwrite(df_out, out,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t")
cat("--- The LEOPARD linear combination weights in order are:  \n")
print(df_out)
cat("\n")
cat("### Finsh estimating linear combination weights! ###\n")
cat("\n")