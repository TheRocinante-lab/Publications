manhattan_fun <- function(df,FDR,ylims,titleplot){
  #filename <- paste("Gwas_",titleplot,sep= "")  

results <- df
    names(results) = c("SNP","CHR","BP","P")
    #head(results)
    sorting <- sort(results$P,decreasing=T)
    ### ADD FDR and Bonferroni Correction
    cutoffbonf=round(-log10(0.05/length(results$P)),3) ## Bonferroni
    cutoff <- FDR 
    signi <- as.character(results$SNP[which(results$P>cutoff)])
    #as.character(results$SNP[which(sorting>cutoff)])
    hlight <- signi

    results <- results %>% 
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>% 
    left_join(results, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>% 
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P < cutoff, "yes", "no"))
    
    axisdf <- results %>% group_by(CHR) %>% summarise(center=( max(BPcum) + min(BPcum) ) / 2 )
    
    figura <- ggplot(results, aes(x=BPcum, y=P)) +
      # Show all points
      geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
      labs(title=titleplot,x = "Chromosome",y="-log10(pvalues")
      #Add highlighted points
      if (length(hlight)<=0){
        figura +geom_hline(yintercept = cutoffbonf,col="blue",linetype="dashed",alpha=0.5) +
          geom_hline(yintercept = FDR, linetype="dashed",col="red")+
          theme_classic(base_size = 10) +
          theme(
            legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.x = element_text(angle=45,hjust=1)) +
          annotate("text",x =2500, y = 7.5, col="blue",size=2,label = "----Bonferroni") +
          annotate("text",x =2500, y = 7, col="red",size=2,label = "----FDR") 
        #ggsave(figura, file = filename, 
        #     device = "png",
        #       width = 9,
        #        height = 8)
      } else{figura + geom_point(data=subset(results, is_highlight=="yes"), color="red", size=2) +
          geom_label_repel(data=results[results$is_highlight=="yes",], aes(label=as.factor(signi),alpha=0.7), 
                           size=2, force=1.3)}+
       
    
      #geom_point(data=subset(results, is_highlight=="yes"), color="red", size=2) +
      
      
      #geom_label_repel(data=results[results$is_highlight=="yes",], aes(label=as.factor(signi),alpha=0.7), 
                      # size=2, force=1.3) +
      
      # add genome-wide sig and sugg lines
      geom_hline(yintercept = cutoffbonf,col="blue",linetype="dashed",alpha=0.5) +
      geom_hline(yintercept = FDR, linetype="dashed",col="red")+
      theme_classic(base_size = 10) +
      theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1)) +
        annotate("text",x =2500, y = 7.5, col="blue",size=4,label = "----Bonferroni") +
        annotate("text",x =2500, y = 7, col="red",size=4,label = "----FDR") 
        #ggsave(figura, file = filename, 
         #     device = "png",
        #       width = 9,
         #        height = 8)
    
}

#annotate("text",x =2600, y = cutoffbonf+0.3, col="blue",label = "Bonferroni") +
#annotate("text",x =2700, y = FDR-0.45, col="red",label = "FDR") 

#manhattan_fun(gwasresults_DTB,4.98,4.29,c(0,8),"DTB")
#FDR=4.29
#titleplot="SUM_T2_HT2"
#df <- gwasresults_SUMnaive
#ylims=c(0,8)

