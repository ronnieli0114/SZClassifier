# This script runs a t-test on the 
# logit-transformed beta values between
# cases (schizophrenia) and controls

library("car")
library("dplyr")
library("arrow")
library("parallel")

# load data
num.cores <- 12
data.dir <- "/projects/compbio/users/yli459/sz_classifier/data"
cpg.meth <- arrow::read_feather(file.path(data.dir, "GSE147221_methylation.ftr")) %>% as.data.frame()
subj.info <- read.table(file.path(data.dir, "GSE147221_subject_info.txt"), sep='\t', header=TRUE)
message('Loaded data')

# define function to test on one CpG
cpg.ttest <- function(probe) {
  
  result.list <- list()
  
  cpg.data <- subset(cpg.meth, cpg==probe) %>% select(-cpg) %>% t() %>% as.data.frame()
  colnames(cpg.data) <- c("beta")
  
  Xy <- merge(cpg.data, subj.info, by.x='row.names', by.y='sample_id')
  Xy$beta <- logit(Xy$beta)
  
  beta.case <- subset(Xy, target=='Case')$beta
  beta.control <- subset(Xy, target=='Control')$beta
  
  mean.diff <- mean(beta.case) - mean(beta.control)
  
  # conduct t-test
  ttest.res <- t.test(beta.case, beta.control)
  p.val <- ttest.res$p.value
  t <- ttest.res$statistic
  # print(sprintf("%s    pval=%s", probe, p.val))
  
  result.list[['cpg']] <- probe
  result.list[['mean_diff']] <- mean.diff
  result.list[["t"]] <- t
  result.list[['p_nominal']] <- p.val
  
  return(data.frame(result.list))
  
}

results <- mclapply(cpg.meth$cpg, cpg.ttest, mc.cores=num.cores)
results.df <- do.call(rbind, results)
results.df$p_adj <- p.adjust(results.df$p_nominal)

write.table(results.df, file.path(data.dir, 'SZ_vs_control_ttest.txt'), sep='\t', quote=F, row.names=F, col.names=T)
message('=== done ===')