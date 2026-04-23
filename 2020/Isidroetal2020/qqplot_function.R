quantileplot <- function(df,trait){
cols <- "gray10"
#cols <- c("#ff5b1e","#ff9900","#c8e717","#25d181","#45818e")


n <- length(df)
datafr <- data.frame(observed =sort(df), expected=sort(-log10(c(1:n)*(1/n))))

  log10Pe <- expression(paste("Expected -log"[10], plain(Pvalues)))
  log10Po <- expression(paste("Observed -log"[10], plain(Pvalues)))
  ggplot(datafr) +
    geom_point(aes(expected, observed), shape = 19, size = 2,col=sample(cols,1)) +
    geom_abline(intercept = 0, slope = 0.95, alpha = 1) +
    #geom_line(aes(expected, cupper), linetype = 2) +
    #geom_line(aes(expected, clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po)+
    scale_y_continuous(expand = c(0, 0), limits = c(0,8))+
    theme_classic()+
    theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(hjust=2),
    text=element_text(size=12, family="Times New Roman"))+
    ggtitle(trait)
}

##8B4500, #82B446
#quantileplot (gwasresults_Lodging$Lodging)

#geom_abline(intercept = -0.04, slope = 0.93, alpha = 0.5) +
