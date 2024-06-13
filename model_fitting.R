# the module takes in data and parameters to fit models

### load packages
library(dplyr)
library(tidyr)
library(stringr)
library(pROC)
library(h2o)

fit_model <- function(feature_data = feature_data,
                      feature_type = "FSD,CNV, FSC, NP", # feature_type is a list of feature types separated by commas
                      sample_info = sample_info, # sample_contains train/valid group information
                      output_path = "./output/"){
  
  ### initiate h2o
  # h2o.init(port = 18181 + sample(1:1000,1), nthreads = 16, max_mem_size = "256G")
  
  ### create output directory
  dir.create(output_path, recursive = T)
  
  ### label feature data; extract train/valid frames based on sample_info
  feature_data_labeled <- feature_data %>%
    inner_join(select(sample_info, ID, train, fold_column), by = c("ID")) %>%
    mutate(Train_Group = factor(Train_Group, levels = c("Healthy","Cancer")))
    
  train_frame <- feature_data_labeled %>%
    filter(train == "training")
  
  valid_frame <- feature_data_labeled %>%
    filter(train == "validation")
  
  ### create score frames for export
  train_score <- select(train_frame, ID, Train_Group) 
  valid_score <- select(valid_frame, ID, Train_Group)
  
  ### start fitting base models by specified feature type
  model_list <- list()
  predictor_list <- list()
  leaderboard_top <- list()
  cols <- colnames(feature_data)
  leaderboard_top_output <- data.frame()
  
  for(feature_type_tmp in unlist(str_split(feature_type,","))){
    
    predictor_list[[feature_type_tmp]] <- cols[grepl(paste0(feature_type_tmp,"_"), cols)]
  
    model_list[[feature_type_tmp]] <- h2o.automl(x = predictor_list[[feature_type_tmp]],
                                             y = "Train_Group",
                                             
                                             project_name = feature_type_tmp,
                                             
                                             fold_column = "fold_column", # pre-specified fold column to allow for reproducibility
                                             balance_classes = T,
                                             
                                             keep_cross_validation_predictions = T,
                                             keep_cross_validation_fold_assignment = T,
                                             
                                             stopping_metric = "AUC", 
                                             max_runtime_secs_per_model=1200, 
                                             
                                             include_algos = c("XgBoost","GLM","DRF","DeepLearning"),

                                             training_frame = as.h2o(train_frame), 

                                             max_models = 200, 
                                             seed = 99)
    
    ##### get leaderboard
    leaderboard_top[[feature_type_tmp]] <- model_list[[feature_type_tmp]]@leaderboard %>% as.data.frame() %>%
      mutate(Train_AUC = NA, Valid_AUC = NA, 
             Train_cutoff95 = NA, Train_sens95 = NA, Train_spec95 = NA, Valid_sens95 = NA, Valid_spec95 = NA,
             Train_cutoff98 = NA, Train_sens98 = NA, Train_spec98 = NA, Valid_sens98 = NA, Valid_spec98 = NA) %>%
      head(5)
    
    ##### save top5 base models by AUC in training frame (5XCV)
    for(i in 1:5){
      
      model_i <- h2o.getModel(leaderboard_top[[feature_type_tmp]]$model_id[i])
      varimp_i <- h2o.varimp(model_i) %>% as.data.frame()
      h2o.saveModel(model_i, path = paste0(output_path,"/Base_model_",feature_type_tmp,"_",i), 
                    export_cross_validation_predictions = T, force = T)
      write.csv(varimp_i, paste0(output_path,"/Base_model_",feature_type_tmp,"_",i,
                                 "/varimp_",feature_type_tmp,"_",i,".csv"), row.names = F)
      
      ####### save performance metrics
      leaderboard_top[[feature_type_tmp]]$Train_AUC[i] = h2o.auc(h2o.performance(model_i, xval = T))

      if(nrow(valid_frame) > 0){
        leaderboard_top[[feature_type_tmp]]$Valid_AUC[i] = h2o.auc(h2o.performance(model_i, newdata = as.h2o(valid_frame)))
      }
      
      ####### get predict score
      model_name <- paste0(feature_type_tmp,"_",i)
      train_score[, model_name] = as.numeric(as.data.frame(h2o.getFrame(model_i@model$cross_validation_holdout_predictions_frame_id$name))[,"Cancer"])
      train_score$fold_column = train_frame$fold_column
      
      ####### get roc metrics of train frame
      roc_tmp <- roc(train_score$Train_Group,train_score[,model_name],plot = F,levels = c("Healthy", "Cancer"))
      coords95_tmp <- coords(roc = roc_tmp, x = 'all', transpose = F, as.matrix = T) %>%
        as.data.frame() %>%
        filter(specificity >= 0.95) %>%
        filter(specificity <= min(specificity)) %>%
        arrange(desc(sensitivity))
      leaderboard_top[[feature_type_tmp]]$Train_cutoff95[i] <- coords95_tmp$threshold[1]
      leaderboard_top[[feature_type_tmp]]$Train_sens95[i] <- coords95_tmp$sensitivity[1]
      leaderboard_top[[feature_type_tmp]]$Train_spec95[i] <- coords95_tmp$specificity[1]
      
      coords98_tmp <- coords(roc = roc_tmp, x = 'all', transpose = F, as.matrix = T) %>%
        as.data.frame() %>%
        filter(specificity >= 0.98) %>%
        filter(specificity <= min(specificity)) %>%
        arrange(desc(sensitivity))
      leaderboard_top[[feature_type_tmp]]$Train_cutoff98[i] <- coords98_tmp$threshold[1]
      leaderboard_top[[feature_type_tmp]]$Train_sens98[i] <- coords98_tmp$sensitivity[1]
      leaderboard_top[[feature_type_tmp]]$Train_spec98[i] <- coords98_tmp$specificity[1]
      
      ####### get roc metrics of valid frame
      if(nrow(valid_frame) > 0){
        valid_score[, model_name] = as.numeric(as.data.frame(h2o.predict(model_i, newdata=as.h2o(valid_frame)))[,"Cancer"])
        roc_valid_tmp <- roc(valid_score$Train_Group,valid_score[,model_name],plot = F,levels = c("Healthy", "Cancer"))
        
        coords95_valid_tmp <- coords(roc = roc_valid_tmp, x = coords95_tmp$threshold[1], input = "threshold", ret = "all")
        leaderboard_top[[feature_type_tmp]]$Valid_sens95[i] <- coords95_valid_tmp$sensitivity[1]
        leaderboard_top[[feature_type_tmp]]$Valid_spec95[i] <- coords95_valid_tmp$specificity[1]
        
        coords98_valid_tmp <- coords(roc = roc_valid_tmp, x = coords98_tmp$threshold[1], input = "threshold", ret = "all")
        leaderboard_top[[feature_type_tmp]]$Valid_sens98[i] <- coords98_valid_tmp$sensitivity[1]
        leaderboard_top[[feature_type_tmp]]$Valid_spec98[i] <- coords98_valid_tmp$specificity[1]
      }
    } # end of top5 base model iteration
    
    tmp_output <- leaderboard_top[[feature_type_tmp]] %>% mutate(feature_type = feature_type_tmp)
    if(nrow(leaderboard_top_output)){
      leaderboard_top_output <- rbind(leaderboard_top_output, tmp_output)
    }else{
      leaderboard_top_output <- tmp_output
    }
  
  } # end of feature type iteration
  
  ##### export layer 1 results
  write.csv(leaderboard_top_output, paste0(output_path,"Layer1_metrics.csv"),row.names = F)
  write.csv(train_score, paste0(output_path,"Layer1_train_score.csv"),row.names = F)
  if(nrow(valid_frame) > 0){
    write.csv(valid_score, paste0(output_path,"Layer1_valid_score.csv"),row.names = F)
  }
  
  h2o.removeAll()

} # end of model fitting function

