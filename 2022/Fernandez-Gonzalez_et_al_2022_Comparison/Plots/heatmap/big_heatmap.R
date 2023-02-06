library(ComplexHeatmap)
library(circlize)



################################################################################
###############################Prepare data#####################################
################################################################################




#set datasets considered
dataset_names <- c("Maize", "Rice", "Soy", "Spruce", "Sorghum",
                   "Switchgrass", "Rice_high_PS")

#we will do the average of the accuracy of the tree models used
models <- c("GBLUP", "BayesB", "RKHS")


#traits <- c("HT", "YLD", "FT",  "R8", "DBH", "DE", "MO", "AN", "ST", "FP", "PC", "simulated1")

#initialize the matrix we will represent as heatmap
heatmap_matrix <- matrix(data = NA, nrow = 15, ncol = 0)

#the rownames are the methods
rownames(heatmap_matrix) <- c("StratSamp", "PAM", "Rscore", "CDmean", "OvClustCDmean", "WIClustCDmean", "Avg_GRM",       
                              "Avg_GRM_MinMax", "Avg_GRM_self", "Avg_GRMtarg", "Avg_GRM_MinMaxtarg", "Rscoretarg", 
                              "CDmeantarg", "OvClustCDmeantarg", "WIClustCDmeantarg" )



#for each dataset
for (dataset in dataset_names) {
  
    #load performance data for each model
    for (model in models) {
      
      
      if (file.exists(paste("../Data/AUC_", dataset, "_", model, ".RData", sep = ""))) {
        
        #the data is stored as a dataframe (summary_new) containing the performance
        #(measured as area under the curve, AUC) of each method for each trait
        #relative to random sampling
        load(paste("../Data/AUC_", dataset, "_", model, ".RData", sep = ""))
        
        #assign summary_new to a new object with the dataset in its name so that
        #it is not overwritten. We remove the columns with the method names.
        #That would usually be a problem, but the method names in "summary_new"
        #are in the same order as the ones we have set in rownames(heatmap_matrix)
        assign(paste("placeholder_", model, sep = ""), summary_new[,3:ncol(summary_new)])
        
      }
      
    }
    
    

    #elementwise mean of the matrices for the 3 different models
    placeholder_list <- list(placeholder_GBLUP,placeholder_BayesB,placeholder_RKHS)
    placeholder <- Reduce("+",placeholder_list)/length(placeholder_list)

    #add data to heatmap_matrix. New dataset-trait combinations are included
    #as new columns
    heatmap_matrix <- cbind(heatmap_matrix, placeholder)
    

  
}



save(heatmap_matrix, file = "Heatmap_matrix.RData")








################################################################################
#################################Make plot######################################
################################################################################

#Use ComplexHeatmap
#https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html

#In "heatmap_matrix" each dataset is present 4 times (once per trait) except 
#Rice_high_PS, which has 5 traits and therefore appears 5 times
#Store the dataset names in the correct order and the correct number of times
#in "dataset_groups"
dataset_groups <- c()

for (name in dataset_names) {
  dataset_groups <- c(dataset_groups, rep(c(name), 4))
}

#add Rice_high_PS the fifth time
dataset_groups <- c(dataset_groups, c("Rice_high_PS"))

#rename Rice_high_PS to RicePopStr because we decided to change its name
#in a later stage.
dataset_groups <- gsub("Rice_high_PS", "RicePopStr", dataset_groups)


#assign a number to each dataset and do the same thing we did before
dataset_groups_numeric <- rep(1:7, each = 4)
dataset_groups_numeric <- c(dataset_groups_numeric, 7)


#the first 9 methods are untargeted and the next 6 ones are targeted. 
#store that in "targ_untarg"
targ_untarg <- c(rep("Untargeted", 9), rep("Targeted", 6))

#same as before but numeric
targ_untarg_numeric <- c(rep(1, 9), rep(2, 6))


#above the heatmap we want to add blocks corresponging to the dataset names we 
#have stored in "dataset_groups". We can do that with HeatmapAnnotation
ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  foo = anno_block(gp = gpar(fill = 2:8), labels = unique(dataset_groups))
)



#Make nonlinear color scale to deal with extreme values:
#circlice package needed!
col_fun<-colorRamp2(c(-50, -10, 0, 10, 50), c("#2A18B5", 
                                              "#7D74C5", 
                                              "#DCDCDC", 
                                              "#C57474",
                                              "#B51818"))

colnames(heatmap_matrix) <- gsub("simulated1", "Simulated", colnames(heatmap_matrix))
#rename trait "simulated1" to Simulated


Heatmap(as.matrix(heatmap_matrix),
        col = col_fun, #color scale
        column_split = dataset_groups_numeric, #group columns by dataset (in numeric form)
        row_split = targ_untarg_numeric,  #group rows by targeted-untargeted (in numeric form)
        column_order = 1:ncol(heatmap_matrix), #do not reorder columns
       cluster_rows = FALSE, #we don't want any clustering
       cluster_columns = FALSE, #we don't want any clustering
       row_labels =c("StratSamp", "PAM", "Rscore", "CDmean", "OvClustCDmean", "WIClustCDmean", "Avg_GRM",       
                     "Avg_GRM_MinMax", "Avg_GRM_self", "Avg_GRM", "Avg_GRM_MinMax", "Rscore", 
                     "CDmean", "OvClustCDmean", "WIClustCDmean" ), #rownames
       top_annotation = ha, #add in the top the names for the blocks corresponging
       #to the different datasets
       rect_gp = gpar(col = "black", lwd = 1), #black borders for the cells
       left_annotation =  rowAnnotation( #add to the left names for the blocks
         #corresponding to targeted-untargeted methods
         foo = anno_block(gp = gpar(fill = 2:3), labels = unique(targ_untarg))
       )
      )

      
#save the heatmap
dev.copy(png, "Big_heatmap_AUC.png",
         units="cm", width=(1+sqrt(5))/2*16, height=16, res=300)
dev.off() 
dev.off() 
  

#The heatmap already had a good legend, but a better one can be made,
#controling its dimensions and title
lgd = Legend(col_fun = col_fun, title = "AUC Gain (%)", at = c(-50,-25,0,25,50), 
             legend_height = unit(12, "cm"),
             legend_width = unit(4, "cm"))
draw(lgd)

#save the legend
dev.copy(png, "Heatmap_legend.png",
         units="cm", width=8, height=20, res=300)
dev.off() 
    
#We have generated a png file with the heatmap and another one with the legend
#If we want to include the better legend into the heatmap, we can do it with 
#any image editing tool.




