library(ggplot2)
library(patchwork)



################################################################################
###############################Prepare data#####################################
################################################################################



#Set datasets, optimization methods and target accuracies to be considered
#in the plot. The method "Avg_GRM_MinMax is not included here because it
#works in a different way to the rest and it will be added manually later
datasets <- c("Maize", "Rice", "Rice_high_PS", "Sorghum", "Spruce", "Switchgrass")
targeted_methods <- c("CDtarg", "Rscoretarg")
untargeted_methods <- c("CD", "Rscore", "Avg_GRM_self")
methods <- c(untargeted_methods, targeted_methods)
tested_sizes <- c("0.95", "0.99")



#Initialize list for the optimal size selected by each method
selected_sizes <- vector("list", length(methods))
names(selected_sizes) <- methods

for (method in methods) {
  selected_sizes[[method]] <- vector("list", length(tested_sizes))
  names(selected_sizes[[method]]) <- tested_sizes
}


#Initialize list for the actual accuraciy obtained with the optimized training set sizes
achieved_accuracies <- vector("list", length(methods))
names(achieved_accuracies) <- methods

for (method in methods) {
  achieved_accuracies[[method]] <- vector("list", length(tested_sizes))
  names(achieved_accuracies[[method]]) <- tested_sizes
}


#initialize vectors for the optimal size and corresponding accuracy for 
#targeted Avg_GRM_MinMax
selected_sizes_Avg_GRM_MinMaxtarg <- c()
achieved_accuracies_Avg_GRM_MinMaxtarg <- c()

#for all datasets considered
for(dataset in datasets) {
  #load data:
  load(paste("../Data/", dataset, "_predictions.RData", sep = ""))
  #The list "predictions" contains for each optimization method the optimal
  #training set size (x) and the corresponding value of the evaluation metric (y)
  #Each optimization method has several optimal values for different target 
  #accuracies (0.90; 0.95; 0.99) except Avg_GRM_MinMaxtarg, that has a real optimal value
  load(paste("../Data/", dataset, "_interpolations_summary.RData", sep = ""))
  #in the interpolations summary there are several dataframes with the mean,
  #median and sd accross iterations of the interpolated accuracy for the 
  #optimized training set sizes
  load(paste("../Data/", dataset, "_training_optimization_results_merged.RData", sep = ""))
  #"Optimization resluts" contains a lot of results but it will only be used
  #here to get the size of the candidate set

  #To avoid overwriting, assign the loaded data to new variables with the
  #dataset name in it. This could also be done with lists.
  assign(paste("Entire_CS_size_", dataset, sep = ""),
         length(Optimization_results[[paste(dataset, "Entire_CS", sep = "")]][[1]][[1]]))
  
  assign(paste("predictions", dataset, sep = ""),
         predictions)
  
  assign(paste("table_mean_df", dataset, sep = ""),
         table_mean_df)
  
  assign(paste("table_median_df", dataset, sep = ""),
         table_median_df)
  
  assign(paste("table_sd_df", dataset, sep = ""),
         table_sd_df)
  
  
  #store the optimal size found by Avg_GRM_MinMaxtarg relative to the 
  #entire candidate set
  selected_sizes_Avg_GRM_MinMaxtarg <- c(selected_sizes_Avg_GRM_MinMaxtarg,
                                         predictions[["Avg_GRM_MinMaxtarg"]][["x"]][["opt"]]/
                                           get(paste("Entire_CS_size_", dataset, sep = "")))
  
  #store the interpolated accuracy for the optimal size found by Avg_GRM_MinMaxtarg
  #for all traits except the simulated one
  achieved_accuracies_Avg_GRM_MinMaxtarg <- c(achieved_accuracies_Avg_GRM_MinMaxtarg,
                                              unlist(get(paste("table_mean_df", dataset, sep = ""))["Avg_GRM_MinMaxtarg", -ncol(get(paste("table_mean_df", dataset, sep = "")))]))
  

  
  #store the optimal size and corresponding accuracies for the other methods
  for (method in methods) {
    
    for (size in tested_sizes) {
      
      selected_sizes[[method]][[size]] <- c(selected_sizes[[method]][[size]],
                                            predictions[[method]][["x"]][[size]]/get(paste("Entire_CS_size_", dataset, sep = "")))
      
      #use the mean
      achieved_accuracies[[method]][[size]] <- c(achieved_accuracies[[method]][[size]],
                                                 unlist(get(paste("table_mean_df", dataset, sep = ""))[paste(method, "_", size, "_relative", sep = ""), -ncol(get(paste("table_mean_df", dataset, sep = "")))]))
      

      
      
    }
    
  }
  

}





#convert list "achieved_accuracies" to dataframes for untargeted methods
for (size in tested_sizes) {

  assign(paste("df_", size, "_acc", sep = ""), c())

  for (method in untargeted_methods) {
    placeholder <- data.frame(Accuracy = achieved_accuracies[[method]][[size]]*100,
                              Method = method,
                              Target_accuracy = paste("TA: ", as.numeric(size)*100, " %", sep = ""),
                              Target = as.numeric(size)*100)


    assign(paste("df_", size, "_acc", sep = ""),
           rbind(get(paste("df_", size, "_acc", sep = "")), placeholder))
  }



}

#merge dataframes for the considered target accuracies
bigdf_untarg_acc <- rbind(df_0.95_acc, df_0.99_acc)

colnames(bigdf_untarg_acc)[1] <- "Percentage of maximum accuracy"




#convert list "achieved_accuracies" to dataframes for targeted methods
for (size in tested_sizes) {
  
  assign(paste("df_", size, "_acc", sep = ""), c())
  
  for (method in targeted_methods) {
    placeholder <- data.frame(Accuracy = achieved_accuracies[[method]][[size]]*100,
                              Method = method,
                              Target_accuracy = paste("TA: ", as.numeric(size)*100, " %", sep = ""),
                              Target = as.numeric(size)*100
                              )
    
    
    assign(paste("df_", size, "_acc", sep = ""),
           rbind(get(paste("df_", size, "_acc", sep = "")), placeholder))
  }
  
  
  
}


#merge dataframes for the considered target accuracies
bigdf_targ_acc <- rbind(df_0.95_acc, df_0.99_acc)

colnames(bigdf_targ_acc)[1] <- "Percentage of maximum accuracy"




#convert list "selected_sizes" to dataframes for untargeted methods
for (size in tested_sizes) {
  
  assign(paste("df_", size, "_size", sep = ""), c())
  
  for (method in untargeted_methods) {
    placeholder <- data.frame(TRS_size = selected_sizes[[method]][[size]]*100, 
                              Method = method,
                              Target_accuracy = paste("TA: ", as.numeric(size)*100, " %", sep = ""))
    
    
    assign(paste("df_", size, "_size", sep = ""), 
           rbind(get(paste("df_", size, "_size", sep = "")), placeholder))
  }
  
  
  
}


#merge dataframes for the considered target accuracies
bigdf_untarg <- rbind(df_0.95_size, df_0.99_size)

colnames(bigdf_untarg)[1] <- "TRS size (percentage of CS)"





#convert list "selected_sizes" to dataframes for targeted methods
for (size in tested_sizes) {
  
  assign(paste("df_", size, "_size", sep = ""), c())
  
  for (method in targeted_methods) {
    placeholder <- data.frame(TRS_size = selected_sizes[[method]][[size]]*100, 
                              Method = method,
                              Target_accuracy = paste("TA: ", as.numeric(size)*100, " %", sep = ""))
    
    
    assign(paste("df_", size, "_size", sep = ""), 
           rbind(get(paste("df_", size, "_size", sep = "")), placeholder))
  }
  
  
  
}




#merge dataframes for the considered target accuracies
bigdf_targ <- rbind(df_0.95_size, df_0.99_size)

colnames(bigdf_targ)[1] <- "TRS size (percentage of CS)"



#convert vectors with data for Avg_GRM_MinMax to dataframes
MinMaxtarg_df_acc <- data.frame(Accuracy = achieved_accuracies_Avg_GRM_MinMaxtarg*100,
                                Method = "Avg_GRM_MinMaxtarg",
                                Target_accuracy = "Opt")

colnames(MinMaxtarg_df_acc)[1] <- "Percentage of maximum accuracy"




MinMaxtarg_df_size <- data.frame(TRS_size = selected_sizes_Avg_GRM_MinMaxtarg*100, 
                                                    Method = "Avg_GRM_MinMaxtarg",
                                                    Target_accuracy = "Opt")

colnames(MinMaxtarg_df_size)[1] <- "TRS size (percentage of CS)"





#combine data for accuracy:
#For each method and target accuracy there will be one data point for each 
#dataset-trait combination (19 in total). Each data point is an interpolation
#done using the average of 40 repetitions in the cross-validation

#To combine dataframes, we need a new column stating if the method corresponds to  
#"Targeted" or "Untargeted" scenario
bigdf_untarg_acc <- cbind(bigdf_untarg_acc, "Untargeted")
colnames(bigdf_untarg_acc)[length(colnames(bigdf_untarg_acc))] <- c("Scenario")
bigdf_untarg_acc$Method <- gsub(pattern = "CD", replacement = "CDmean", x = bigdf_untarg_acc$Method)
#rename CD to CDmean

bigdf_targ_acc <- cbind(bigdf_targ_acc, "Targeted")
colnames(bigdf_targ_acc)[length(colnames(bigdf_targ_acc))] <- c("Scenario")
bigdf_targ_acc$Method <- gsub(pattern = "CDtarg", replacement = "CDmean", x = bigdf_targ_acc$Method)
bigdf_targ_acc$Method <- gsub(pattern = "Rscoretarg", replacement = "Rscore", x = bigdf_targ_acc$Method)
#Remove "targ" after targeted methods, it is no longer needed as we have a 
#column stating if method is targeted or untargeted

#Add targeted #Avg_GRM_MinMax
MinMaxtarg_df_acc <- cbind(MinMaxtarg_df_acc, "Targeted")
colnames(MinMaxtarg_df_acc)[length(colnames(MinMaxtarg_df_acc))] <- c("Scenario")
MinMaxtarg_df_acc$Target <- NA #Avg_GRM_MinMax doesn't need a target accuracy
MinMaxtarg_df_acc$Method <- gsub(pattern = "Avg_GRM_MinMaxtarg", replacement = "Avg_GRM_MinMax", x = MinMaxtarg_df_acc$Method)

#merge all dataframes
combine_df_acc <- rbind(bigdf_untarg_acc, bigdf_targ_acc, MinMaxtarg_df_acc)
  


#combine data for optimized training set size:
#For each method and target accuracy there will be one data point for each 
#dataset (6 in total)

#To combine dataframes, we need a new column stating if the method corresponds to  
#"Targeted" or "Untargeted" scenario
bigdf_untarg <- cbind(bigdf_untarg, "Untargeted")
colnames(bigdf_untarg)[length(colnames(bigdf_untarg))] <- c("Scenario")
bigdf_untarg$Method <- gsub(pattern = "CD", replacement = "CDmean", x = bigdf_untarg$Method)
#rename CD to CDmean

bigdf_targ <- cbind(bigdf_targ, "Targeted")
colnames(bigdf_targ)[length(colnames(bigdf_targ))] <- c("Scenario")
bigdf_targ$Method <- gsub(pattern = "CDtarg", replacement = "CDmean", x = bigdf_targ$Method)
bigdf_targ$Method <- gsub(pattern = "Rscoretarg", replacement = "Rscore", x = bigdf_targ$Method)
#Remove "targ" after targeted methods, it is no longer needed as we have a 
#column stating if method is targeted or untargeted


#Add targeted #Avg_GRM_MinMax
MinMaxtarg_df_size <- cbind(MinMaxtarg_df_size, "Targeted")
colnames(MinMaxtarg_df_size)[length(colnames(MinMaxtarg_df_size))] <- c("Scenario")
MinMaxtarg_df_size$Method <- gsub(pattern = "Avg_GRM_MinMaxtarg", replacement = "Avg_GRM_MinMax", x = MinMaxtarg_df_size$Method)

#merge all dataframes
combine_df <- rbind(bigdf_untarg, bigdf_targ, MinMaxtarg_df_size)













################################################################################
#################################Make plot######################################
################################################################################


#Make boxplots for the accuracy obtained with optimized training set sizes
 (combined_acc <- ggplot(combine_df_acc, aes(y=Method, x = `Percentage of maximum accuracy`)) + 
     geom_boxplot(aes(color = Method, fill= Method), alpha = 0.3)+
     #Make 6 different plots:
     #Two columns for targeted vs untargeted
     #3 rows for the two target accuracies considered plus the true optimal size
     #obtained by Avg_GRM_MinMax
     facet_grid(vars(`Target_accuracy`), vars(Scenario)
                , scales = "free_y"#scales = "free_y" is needed because otherwise 
                #ggplot tries to include all for methods in each subplot
     ) +
     #change theme
     theme_bw() +
     theme(axis.title.y=element_blank(),
           axis.text.y=element_blank(),
           axis.ticks.y=element_blank()) +
     geom_vline(data = combine_df_acc,  #add vertical line indicating target accuracy
                aes(xintercept = Target),
                linetype = "dashed", color = "blue") + 
     ggtitle("B) Concordance between real and target accuracy")) #add title
#There will be a warning about missing values because we don't have all methods
#in each subplot. It is normal




#Make boxplots for the obtained optimized training set sizes
 (combined_size <- ggplot(combine_df, aes(y=Method, x = `TRS size (percentage of CS)`)) + 
     geom_boxplot(aes(color = Method, fill= Method), alpha = 0.3) +
     facet_grid(vars(`Target_accuracy`), vars(Scenario)
                #Make 6 different plots:
                #Two columns for targeted vs untargeted
                #3 rows for the two target accuracies considered plus the true optimal size
                #obtained by Avg_GRM_MinMax
                , scales = "free_y"
     ) +
     #change theme
     theme_bw() +
     theme(axis.title.y=element_blank(),
           axis.text.y=element_blank(),
           axis.ticks.y=element_blank()) + 
     ggtitle("A) Selected TRS size")) #add title
 
 
 #Use patchwork package to combine the two boxplots into a single one
 Combined_all <- (combined_size) /
                 (combined_acc)
 
 
 #save the plot
 ggsave(
   filename = "All_boxplots_good.png",
   plot = Combined_all,
   device = png(),
   #using width and height you can control the proportions
   width = 8, 
   height = 6
 )
 
