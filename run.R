package_list <- c("dplyr","tidyr","purrr","stringr","ggplot2","data.table","h2o","splitTools")
installed_packages <- installed.packages()
packages_to_install <- setdiff(package_list, installed_packages[,1])
if(length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}
# Load the packages
sapply(package_list, require, character.only = TRUE)

# initiate h2o instance
h2o.init(port = 18181 + sample(1:1000,1), nthreads = 16, max_mem_size = "256G")
# load the model fitting function
source("./model_fitting.R")
source("./assemble_internal_validation_results.R")

### fitting model requires 1.feature_data, 2. feature_type (a string separated by comma), 
### 3. sample_info, 5. output_path

### prepare inputs for full cfDNA model with external validation

train_frame_cfDNA <- train_frame[,!grepl("Train_",colnames(train_frame))]

feature_data_full_cfDNA <- rbind(train_frame_cfDNA, valid_frame)


fit_model(feature_data = feature_data_full_cfDNA, feature_type = "FSD,CNV, FSC, NP",
          sample_info = sample_info_full,
          output_path = "./cfDNA_output/")
h2o.removeAll()
# h2o.shutdown(prompt = F)

### assemble results for training and validation cohorts
train_score_final <- read.csv(paste0("./train_score.csv")) %>%
  mutate(Final_score = (l1_1 + l1_2 + l1_3+ l1_4+ l1_5+ l1_6+ l1_7+ l1_8+ l1_9+ l1_10)/10)
train_roc <- roc(train_score_final$Train_Group, train_score_final$Final_score, ci = T, levels = c("Healthy","Cancer"))

valid_score_final <- read.csv(paste0("./valid_score.csv")) %>%
  mutate(Final_score =  (l1_1 + l1_2 + l1_3+ l1_4+ l1_5+ l1_6+ l1_7+ l1_8+ l1_9+ l1_10)/10)
valid_roc <- roc(valid_score_final$Train_Group, valid_score_final$Final_score, ci = T, levels = c("Healthy","Cancer"))

print(train_roc$auc)
print(valid_roc$auc)