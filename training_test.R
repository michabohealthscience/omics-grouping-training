#==============================================================================
# Training day
#==============================================================================
root <- 'https://raw.githubusercontent.com/michabohealthscience/training-fsa/main'

sample_metadata_all <- read.csv(file.path(root, 'data/HILIC_POS_male/0_sample_metadata.csv'))

head(sample_metadata_all)


hilic_pos_all <- read.csv(file.path(root, 'data/HILIC_POS_male/1_unfiltered.csv'))

# feature count
feature_c_all = nrow(hilic_pos_all)

# sample count
samp_c_all = ncol(hilic_pos_all)-1

feature_c_all
samp_c_all


hilic_pos_f <- read.csv(file.path(root, 'data/HILIC_POS_male/2_filtered.csv'))

# feature count
feature_c_f = nrow(hilic_pos_f)
# sample count
samp_c_f = ncol(hilic_pos_f)-1


c('All samples:', samp_c_all)
c('Samples removed:', samp_c_all-samp_c_f)
c('Remaining samples:', samp_c_f)

c('All features:', feature_c_all)
c('Features removed:', feature_c_all-feature_c_f)
c('remaining features:', feature_c_f)

# Read in the sample metadata and processed data matrix
sample_metadata <- read.csv(file.path(root, 'data/HILIC_POS_male/0_sample_metadata_filtered.csv'))
glog <- read.csv(file.path(root, 'data/HILIC_POS_male/5_glog.csv'), row.names = 1)

# transpose 
glog_t <- t(glog)

qc_names <- sample_metadata[,1][sample_metadata$Class=='QC']
sample_metadata_no_qcs <- sample_metadata[sample_metadata$Class!='QC',]

glog_t_no_qcs <- glog_t[!rownames(glog_t) %in% qc_names,]


pca_no_qcs <- prcomp(glog_t_no_qcs, center = TRUE, scale. = TRUE)

library(ggfortify)


autoplot(pca_no_qcs, x=1, y=2, data=sample_metadata_no_qcs, colour="test_substance",  shape="dose_group", frame=TRUE, frame.colour = 'test_substance')+
  scale_colour_manual(values=c("khaki","black","blue", 'red', "grey",'orange', 'magenta', 'yellow', 'green', 'brown', 'purple'))+
  scale_fill_manual(values=c("khaki","black","blue", 'red', "grey",'orange', 'magenta', 'yellow', 'green', 'brown', 'purple'))+
  scale_shape_manual(values=c(4,8,2))+
  theme_bw()




pqn <- read.csv(file.path(root, 'data/HILIC_POS_male/3_pqn.csv'))

# Get all the control samples (i.e. all those with dose 0)
control_samples <- sample_metadata[,1][sample_metadata$dose_group==0]

# Get the high dose group (2) for test substance 1
ts1_dose_2_samples <- sample_metadata[,1][sample_metadata$dose_group==2 & sample_metadata$test_substance==1]

# Select the 2000th feature as an example
f2000 <- pqn[2000,]


controldf <- data.frame(intensity = unlist(f2000[,control_samples]), control_treated = 'control')

treateddf <- data.frame(intensity = unlist(f2000[,ts1_dose_2_samples]), control_treated = 'treated')

# Perform a t-test for this feature subset based on the chosen samples
ttest_out <- t.test(treateddf$intensity, controldf$intensity)

ttest_out

ctdf <- rbind(controldf, treateddf)

boxplot(intensity ~ control_treated, data = ctdf, frame = FALSE,  main=f2000['feature_name'])


ttest_ts2 <- read.csv(file.path(root, 'data/ttest_ts2.csv'))

sum(ttest_ts2$pvalue<0.05)

ttest_ts2$qvalue <- p.adjust(ttest_ts2$pvalue, method = 'fdr')
sum(ttest_ts2$qvalue<0.05)

ttest_ts2$diff <- 'NONE'
ttest_ts2$diff[ttest_ts2$tstat<0 & ttest_ts2$qvalue<0.05] <- 'DOWN'
ttest_ts2$diff[ttest_ts2$tstat>0 & ttest_ts2$qvalue<0.05] <- 'UP'


ggplot(data=ttest_ts2, aes(x=log2(fc), y=-log10(qvalue), col=diff)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")




tstats <- read.csv(file.path(root, 'data/combined_male_tstats.csv'), row.names = 1)

tstats_scaled <- scale(tstats, center=FALSE, scale=TRUE)

library(pvclust)
pvclust_res <- pvclust(tstats_scaled, method.dist="euclidean", method.hclust = "ward.D2", nboot=100)


plot(pvclust_res, c("si", "au", "bp"), hang=-1)



library(pvclust)
unblinded_names <- read.csv(file.path(root, 'data/unblinded_ordered_names.csv'), header = TRUE)
pvclust_res_unblind <- pvclust_res
unblinded_names_sorted <- unblinded_names[match(pvclust_res$hclust$labels, unblinded_names$test_substance_dose), ]  
pvclust_res_unblind$hclust$labels <- unblinded_names_sorted$compound_name_dose_moa
plot(pvclust_res_unblind, c("si", "au", "bp"), hang=-1)


plot(pvclust_res_unblind, c("si", "au", "bp"), hang=-1)
