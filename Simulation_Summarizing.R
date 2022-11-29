
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(patchwork)
##################################################
#
#         Set Size Plot
#
##################################################
Set_Size_fun<-function(all_Set_data_dataframe){
  p_size_box = ggplot(data=all_Set_data_dataframe,aes(x = reorder(Method,Size,median), y =Size,fill=Method))+geom_boxplot( aes(fill=Method))+scale_y_continuous(limits = c(0,100))+scale_fill_manual(values=c("MESuSiE"="#023e8a","SuSiE"="#2a9d8f","Paintor"="#f4a261"))+ylab("Set Size")+xlab("")+facet_grid(h2~causal_num,labeller=label_parsed)+theme_bw()
  p_size_box = p_size_box + theme(axis.text.x = element_text( size = 18),
                                  axis.text.y = element_text( size = 18),  
                                  axis.title.x = element_text( size = 18,face="bold"),
                                  axis.title.y = element_text( size = 18,face="bold"),
                                  strip.text.x = element_text(size = 18),
                                  strip.text.y= element_text(size = 18),
                                  legend.text=element_text(size=18),
                                  legend.title=element_text(size=22,face="bold"),
                                  plot.title = element_text(size=16,hjust = 0.5)) 
  return(p_size_box)
}

##################################################
#
#         Set Power Plot
#
##################################################
Set_Power_fun<-function(power_summary){
  p_power_bar = ggplot(data=power_summary,aes(x = Method, y =Power_name,fill=Method))+geom_bar(stat = "identity", position = "dodge", aes(fill=Method),alpha = 1.2)+scale_fill_manual(values=c("MESuSiE"="#023e8a","SuSiE"="#2a9d8f","Paintor"="#f4a261"))+ylab("Power")+xlab("")+facet_grid(h2~causal_num,labeller=label_parsed)+ylim(c(0,1))+theme_bw()
  p_power_bar = p_power_bar + theme(axis.text.x = element_text( size = 18),
                                    axis.text.y = element_text( size = 18),  
                                    axis.title.x = element_text( size = 18,face="bold"),
                                    axis.title.y = element_text( size = 18,face="bold"),
                                    strip.text.x = element_text(size = 18),
                                    strip.text.y= element_text(size = 18),
                                    legend.text=element_text(size=18),
                                    legend.title=element_text(size=22,face="bold"),
                                    plot.title = element_text(size=16,hjust = 0.5)) 
  return(p_power_bar)
}

##################################################
#
#         ROC for either/shared signal detection 
#
##################################################
ROC_shared_fun<-function(all_ROC_data_dataframe,ylab_title){
  p = ggplot(all_ROC_data_dataframe,aes(x = Power, y = 1 - FDR )) + geom_line(aes(x = Power, y = 1 - FDR,linetype = Method,color = Method),size=0.8)+scale_color_manual(values=c("MESuSiE"="#023e8a","SuSiE"="#2a9d8f","Paintor"="#f4a261"))+ scale_linetype_manual(values=c("MESuSiE" = "solid", "SuSiE"="dashed","Paintor" = "twodash")) + xlab("Recall")+ylab("Precision")+theme_bw()+facet_grid(h2~causal_num,labeller=label_parsed)
  p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
  p = p + geom_point(data = all_ROC_data_dataframe[all_ROC_data_dataframe$Cutoff==0.9|all_ROC_data_dataframe$Cutoff==0.5,],aes(x=Power, y=1-FDR, group=Method,color = Method),size =3,shape = 21)+scale_colour_manual(values=c("MESuSiE"="#023e8a","SuSiE"="#2a9d8f","Paintor"="#f4a261"))
  p = p + geom_text_repel(data = all_ROC_data_dataframe[all_ROC_data_dataframe$Cutoff==0.9|all_ROC_data_dataframe$Cutoff==0.5,],aes(x=Power, y=1-FDR,label=Cutoff,color=Method),size=5,show.legend = FALSE)+scale_colour_manual(values=c("MESuSiE"="#023e8a","SuSiE"="#2a9d8f","Paintor"="#f4a261"))
  p = p + theme(axis.text.x = element_text( size = 18),
                axis.text.y = element_text( size = 18),  
                axis.title.x = element_text( size = 18,face="bold"),
                axis.title.y = element_text( size = 18,face="bold"),
                strip.text.x = element_text(size = 18),
                strip.text.y= element_text(size = 18),
                legend.text=element_text(size=18),
                legend.title=element_text(size=22,face="bold"),
                plot.title = element_text(size=16,hjust = 0.5))
  p = p + theme( strip.background = element_rect(fill = "gray"))
  p<-grid.arrange(p,left = textGrob(ylab_title, rot = 90,gp = gpar(fontface = "bold", cex = 2)))
  return(p)
}
###################################################
#
#
#       ROC for ancestry-specific signal detection
#
#
###################################################
ROC_ancestry_fun<-function(all_ROC_data_dataframe,ylab_title){
  p = ggplot(all_ROC_data_dataframe,aes(x = Power, y = 1 - FDR )) + geom_line(aes(x = Power, y = 1 - FDR,linetype = Method,color = Method),size=0.8)+scale_color_manual(values=c("MESuSiE WB"="#023e8a","MESuSiE BB"="#023e8a","SuSiE WB"="#2a9d8f","SuSiE BB"="#2a9d8f","Paintor WB"="#f4a261","Paintor BB"="#f4a261"))+ scale_linetype_manual(values=c("MESuSiE WB" = "dashed", "MESuSiE BB"="solid","SuSiE WB" = "dashed", "SuSiE BB"="solid","Paintor WB" = "dashed", "Paintor BB"="solid")) + xlab("Recall")+ylab("Precision")+theme_bw()+facet_grid(h2~causal_num,labeller=label_parsed)
  p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
  p = p + geom_point(data = all_ROC_data_dataframe[all_ROC_data_dataframe$Cutoff==0.9|all_ROC_data_dataframe$Cutoff==0.5,],aes(x=Power, y=1-FDR, group=Method,color = Method),size =3,shape = 21)+scale_color_manual(values=c("MESuSiE WB"="#023e8a","MESuSiE BB"="#023e8a","SuSiE WB"="#2a9d8f","SuSiE BB"="#2a9d8f","Paintor WB"="#f4a261","Paintor BB"="#f4a261"))
  p = p + geom_text_repel(data = all_ROC_data_dataframe[all_ROC_data_dataframe$Cutoff==0.9|all_ROC_data_dataframe$Cutoff==0.5,],aes(x=Power, y=1-FDR,label=Cutoff,color=Method),size=5,show.legend = FALSE)+scale_color_manual(values=c("MESuSiE WB"="#023e8a","MESuSiE BB"="#023e8a","SuSiE WB"="#2a9d8f","SuSiE BB"="#2a9d8f","Paintor WB"="#f4a261","Paintor BB"="#f4a261"))
  p = p + theme(axis.text.x = element_text( size = 18),
                axis.text.y = element_text( size = 18),  
                axis.title.x = element_text( size = 18,face="bold"),
                axis.title.y = element_text( size = 18,face="bold"),
                strip.text.x = element_text(size = 18),
                strip.text.y= element_text(size = 18),
                legend.text=element_text(size=18),
                legend.title=element_text(size=22,face="bold"),
                plot.title = element_text(size=16,hjust = 0.5))
  p<-grid.arrange(p,left = textGrob(ylab_title, rot = 90,gp = gpar(fontface = "bold", cex = 2)))
  return(p)
}


##################################################
#
#         FDR Power for either/shared signal detection 
#
##################################################

FDR_Power_shared_fun<-function(FDR_Power){
  p<-ggplot(FDR_Power, aes(FDR, Power, fill = Method))+geom_bar(stat="identity", position = "dodge")+theme_bw()+facet_grid(vars(h2),vars(causal_num),labeller=label_parsed)+geom_text(aes(x=FDR,group=Method,y=Power,label=Power),position = position_dodge(width = 1),vjust=-0.5,size = 3,fontface="bold")+ylab("Power")+xlab("FDR")+scale_fill_manual(values=c("MESuSiE"="#023e8a","SuSiE"="#2a9d8f","Paintor"="#f4a261"))+ylim(0,0.7)
  p = p + theme(axis.text.x = element_text( size = 18),
                axis.text.y = element_text( size = 18),  
                axis.title.x = element_text( size = 18,face="bold"),
                axis.title.y = element_text( size = 18,face="bold"),
                strip.text.x = element_text(size = 18),
                strip.text.y= element_text(size = 18),
                legend.text=element_text(size=18),
                legend.title=element_text(size=22,face="bold"),
                plot.title = element_text(size=16,hjust = 0.5))
  return(p)
}
###################################################
#
#     FDR Power for ancestry-specific signal detection
#
#####################################################
FDR_Power_ancestry_fun<-function(FDR_Power){
  p = ggplot(FDR_Power, aes(FDR, Power)) + geom_col_pattern(position = "dodge", aes(pattern =Method, fill = Method), colour= 'black', pattern_density  = 0.35,pattern_key_scale_factor = 1.3) +
    theme_bw() +facet_grid(vars(h2),vars(causal_num),labeller=label_parsed)+ylab("Power")+xlab("FDR")+
    scale_fill_manual(name = "Method",values = c("MESuSiE WB"="#023e8a","MESuSiE BB"="#023e8a","SuSiE WB"="#2a9d8f","SuSiE BB"="#2a9d8f","Paintor WB"="#f4a261","Paintor BB"="#f4a261")) + scale_pattern_manual(values=c("MESuSiE WB"="none","MESuSiE BB"="stripe","SuSiE WB"="none","SuSiE BB"="stripe","Paintor WB"="none","Paintor BB"="stripe"),guide="none")+
    guides(fill = guide_legend(override.aes =  list(
      pattern = c("none","stripe","none", "stripe","none", "stripe"),
      pattern_spacing = .01,
      pattern_angle = c( 0, 45,0, 45,0, 45) )))
  p = p+ theme(axis.text.x = element_text( size = 18),
               axis.text.y = element_text( size = 18),  
               axis.title.x = element_text( size = 18,face="bold"),
               axis.title.y = element_text( size = 18,face="bold"),
               strip.text.x = element_text(size = 18),
               strip.text.y= element_text(size = 18),
               legend.text=element_text(size=18),
               legend.title=element_text(size=22,face="bold"),
               plot.title = element_text(size=16,hjust = 0.5))
  return(p)
}
#######################################################
#
#   PIP calibration for shared signal
#
#########################################################
PIP_calibration_shared_fun<-function(PIP_calibration){
  p = ggplot(data = PIP_calibration,aes(x = mean_pip,y = obs_frq)) + theme_bw()+ geom_pointrange(aes(ymin = obs_min,ymax = obs_max),col="grey") +geom_point(col="black")+geom_abline(intercept = 0, slope = 1,col="red")+ylab("Observed Frequency")+xlab("Mean PIP")+facet_grid(vars(Method),vars(causal_num),labeller=label_parsed)			 
  p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
  p = p +	theme(axis.text.x = element_text( size = 18),
                axis.text.y = element_text( size = 18),  
                axis.title.x = element_text( size = 18),
                axis.title.y = element_text( size = 18),
                strip.text.x = element_text(size = 18),
                strip.text.y= element_text(size = 18),
                legend.text=element_text(size=18),
                legend.title=element_text(size=22,face="bold"),
                plot.title = element_text(size=16,hjust = 0.5))			 
  return(p)
}
#######################################################
#
#   PIP calibration for ancestry specific signal
#
#########################################################
PIP_calibration_ancestry_fun<-function(PIP_calibration){
  p = ggplot(data = PIP_calibration,aes(x = mean_pip,y = obs_frq)) + theme_bw()+ geom_pointrange(aes(ymin = obs_min,ymax = obs_max),col="grey") +geom_point(col="black")+geom_abline(intercept = 0, slope = 1,col="red")+ylab("Observed Frequency")+xlab("Mean PIP")+facet_grid(Method~causal_num,labeller=label_parsed)			 
  p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
  p = p +	theme(axis.text.x = element_text( size = 18),
                axis.text.y = element_text( size = 18),  
                axis.title.x = element_text( size = 18,face="bold"),
                axis.title.y = element_text( size = 18,face="bold"),
                strip.text.x = element_text(size = 18),
                strip.text.y= element_text(size = 18),
                legend.text=element_text(size=18),
                legend.title=element_text(size=22,face="bold"),
                plot.title = element_text(size=16,hjust = 0.5))		
  return(p)
}
#########################################################################
#
#             PIP scatter plot
#
#########################################################################
PIP_scatter_fun<-function(pip_all_dataframe_long){
p = ggplot(data=pip_all_dataframe_long,aes(x=MESuSiE,y=PIP)) +geom_point(alpha=0.3)+ geom_point(data=pip_all_dataframe_long[pip_all_dataframe_long$Signal!="None",],aes(x=MESuSiE,y=PIP,color=Signal),size=3)+scale_colour_manual(values=c("#F2C1B6","#ABDB9F", "#6162B0"))+geom_abline(intercept=0,slope=1, col="red")+ylab("PIP")+theme_bw()
p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
p = p+facet_grid(vars(Method),vars(causal_num),labeller=label_parsed)			 
p = p + theme(axis.text.x = element_text( size = 18),
              axis.text.y = element_text( size = 18),  
              axis.title.x = element_text( size = 18),
              axis.title.y = element_text( size = 18),
              strip.text.x = element_text(size = 18),
              strip.text.y= element_text(size = 18),
              legend.text=element_text(size=18),
              legend.title=element_text(size=22,face="bold"),
              plot.title = element_text(size=16,hjust = 0.5))	
return(p)
}
#########################################################################
#
#             PIP scatter plot shared only
#
#########################################################################
PIP_scatter_shared_fun<-function(pip_all_dataframe_long){
  p = ggplot(data=pip_all_dataframe_long,aes(x=MESuSiE,y=PIP)) +geom_point(alpha=0.3)+ geom_point(data=pip_all_dataframe_long[pip_all_dataframe_long$Signal!="None",],aes(x=MESuSiE,y=PIP,color=Signal),size=3)+scale_colour_manual(values=c( "#6162B0"))+geom_abline(intercept=0,slope=1, col="red")+ylab("PIP")+theme_bw()
  p = p + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0", "0.25", "0.5", "0.75", "1"))
  p = p+facet_grid(vars(Method),vars(causal_num),labeller=label_parsed)			 
  p = p + theme(axis.text.x = element_text( size = 18),
                axis.text.y = element_text( size = 18),  
                axis.title.x = element_text( size = 18),
                axis.title.y = element_text( size = 18),
                strip.text.x = element_text(size = 18),
                strip.text.y= element_text(size = 18),
                legend.text=element_text(size=18),
                legend.title=element_text(size=22,face="bold"),
                plot.title = element_text(size=16,hjust = 0.5))	
  return(p)
}




setwd("/net/fantasia/home/borang/Susie_Mult/Paper_Figure/summarized_data_new")
fig_dir<-"/net/fantasia/home/borang/Susie_Mult/Paper_Figure/figure_new/"
########################################################################
#
#           50% causal variants are shared
#
########################################################################
##########################################################
#
# Either ancestry 
# Set Size&Power | ROC | FDR Power 
#
##########################################################

      ###################
      #
      #Set Size&Power
      #
      ###################
      load("shared_50_either_ancestry_all_Set_data_dataframe_either.RData")
      load("shared_50_either_ancestry_power_summar_eithery.RData")
      p_size_box<-Set_Size_fun(all_Set_data_dataframe)
      p_power_bar<-Set_Power_fun(power_summary)
      size_power<-p_size_box/p_power_bar
      ggsave(paste0(fig_dir,"shared_50_set_size_power.jpeg"),size_power,height =10,width =15,dpi=300)
      ###################
      #
      #ROC
      #
      ###################
      load("shared_50_either_ancestry_all_ROC_data_dataframe_either.RData")
      p<-ROC_shared_fun(all_ROC_data_dataframe,"50% Shared Causal Signals")
      ggsave(paste0(fig_dir,"shared_50_either_ROC.jpeg"),p,height=6,width=14,dpi=300)
      ###################
      #
      #FDR&Power
      #
      ###################
      load("shared_50_either_ancestry_FDR_Power_either.RData")
      p<-FDR_Power_shared_fun(FDR_Power)
      ggsave(paste0(fig_dir,"Shared_50_either_FDR_Power.jpeg"),p,height=6,width=14,dpi=300)
      
##########################################################
#
# Shared Signal 
#  ROC | FDR Power | PIP calibration | PIP scatter plot 
#
##########################################################
      ###################
      #
      #ROC
      #
      ###################
      load(paste0(new_res_dir,"shared_50_shared_all_ROC_data_dataframe_shared.RData"))
      p<-ROC_shared_fun(all_ROC_data_dataframe,"50% Shared Causal Signals")
      ggsave(paste0(fig_dir,"shared_50_shared_ROC.jpeg"),p,height=6,width=14,dpi=300)
      ###################
      #
      #FDR&Power
      #
      ###################
      load(paste0(new_res_dir,"shared_50_shared_FDR_Power_shared.RData"))
      p<-FDR_Power_shared_fun(FDR_Power)+ylim(0,0.8)
      ggsave(paste0(fig_dir,"Shared_50_shared_FDR_Power.jpeg"),p,height=6,width=14,dpi=300)
      ###################
      #
      #PIP calibration
      #
      ###################      
      load( paste0(new_res_dir,"shared_50_shared_pip_calibration.RData"))
      p<-PIP_calibration_shared_fun(PIP_calibration)
      ggsave(paste0(fig_dir,"Shared_50_shared_PIP_calibration.jpeg"),p,height=9,width=9,dpi=300)
      ###################
      #
      #PIP scatter plot
      #
      ################### 
      load(paste0(new_res_dir,"shared_50_shared_pip.RData"))
      p<-PIP_scatter_fun(pip_all_dataframe_long[pip_all_dataframe_long$h2==1,])
      ggsave(paste0(fig_dir,"Shared_50_shared_PIP_h2_1.jpeg"),p,height=7,width=14,dpi=300)
      
      p<-PIP_scatter_fun(pip_all_dataframe_long[pip_all_dataframe_long$h2==2,])
      ggsave(paste0(fig_dir,"Shared_50_shared_PIP_h2_2.jpeg"),p,height=7,width=14,dpi=300)
      
##########################################################
#
# Ancestry-specific Signal 
#  ROC | FDR Power | PIP calibration  
#
##########################################################
      ###################
      #
      #ROC
      #
      ###################
      load(paste0(new_res_dir,"shared_50_shared_all_ROC_data_dataframe_ancestry.RData"))
      p<-ROC_ancestry_fun(all_ROC_data_dataframe,"50% Shared Causal Signals")
      ggsave(paste0(fig_dir,"shared_50_ancestry_ROC.jpeg"),p,height=6,width=14,dpi=300)
      ###################
      #
      #FDR&Power
      #
      ###################
      load(paste0(new_res_dir,"shared_50_ancestry_FDR_Power_ancestry.RData"))
      p<-FDR_Power_ancestry_fun(FDR_Power)
      ggsave(paste0(fig_dir,"shared_50_ancestry_FDR_Power.jpeg"),p,height=6,width=14,dpi=300)
      ###################
      #
      #PIP calibration
      #
      ###################      
      load(paste0(new_res_dir,"shared_50_ancestry_pip_calibration.RData"))
      p<-PIP_calibration_ancestry_fun(PIP_calibration)
      ggsave(paste0(fig_dir,"Shared_50_ancestry_PIP_calibration.jpeg"),p,height=9,width=9,dpi=300)
########################################################################
#
#           100% causal variants are shared
#
########################################################################
      ##########################################################
      #
      # 
      # Set Size&Power | ROC | FDR Power 
      #
      ##########################################################
      
      ###################
      #
      #Set Size&Power
      #
      ###################
    
      load(paste0(new_res_dir,"shared_100_either_ancestry_all_Set_data_dataframe_either.RData"))
      load(paste0(new_res_dir,"shared_100_either_ancestry_power_summar_eithery.RData"))
      
      p_size_box<-Set_Size_fun(all_Set_data_dataframe)
      p_power_bar<-Set_Power_fun(power_summary)
      size_power<-p_size_box/p_power_bar
      ggsave(paste0(fig_dir,"shared_100_set_size_power.jpeg"),size_power,height =10,width =15,dpi=300)
      ###################
      #
      #ROC
      #
      ###################
      load(paste0(new_res_dir,"shared_100_shared_all_ROC_data_dataframe_shared.RData"))
      p<-ROC_shared_fun(all_ROC_data_dataframe,"All Shared Causal Signals")
      ggsave(paste0(fig_dir,"shared_100_shared_ROC.jpeg"),p,height=6,width=14,dpi=300)
      ###################
      #
      #FDR&Power
      #
      ###################
      load(paste0(new_res_dir,"shared_100_shared_FDR_Power_shared.RData"))
      p<-FDR_Power_shared_fun(FDR_Power)+ylim(0,1)
      ggsave(paste0(fig_dir,"Shared_100_shared_FDR_Power.jpeg"),p,height=6,width=14,dpi=300)
      ###################
      #
      #PIP calibration
      #
      ###################      
      load( paste0(new_res_dir,"shared_100_shared_pip_calibration.RData"))
      p<-PIP_calibration_shared_fun(PIP_calibration)
      ggsave(paste0(fig_dir,"Shared_100_shared_PIP_calibration.jpeg"),p,height=9,width=9,dpi=300)
      ###################
      #
      #PIP scatter plot
      #
      ################### 
      load(paste0(new_res_dir,"shared_100_shared_pip.RData"))
      p<-PIP_scatter_shared_fun(pip_all_dataframe_long[pip_all_dataframe_long$h2==1,])
      ggsave(paste0(fig_dir,"Shared_100_shared_PIP_h2_1.jpeg"),p,height=7,width=14,dpi=300)
      
      p<-PIP_scatter_shared_fun(pip_all_dataframe_long[pip_all_dataframe_long$h2==2,])
      ggsave(paste0(fig_dir,"Shared_100_shared_PIP_h2_2.jpeg"),p,height=7,width=14,dpi=300)
########################################################################
#
#           50% causal variants are shared use sub LD
#
########################################################################      
      ##########################################################
      #
      # Shared Signal 
      #  ROC | FDR Power | PIP calibration | PIP scatter plot 
      #
      ##########################################################
      ###################
      #
      #ROC
      #
      ###################
      load(paste0(new_res_dir,"shared_50_sub_shared_all_ROC_data_dataframe_shared.RData"))
      p<-ROC_shared_fun(all_ROC_data_dataframe,"External LD")
      ggsave(paste0(fig_dir,"shared_50_sub_shared_ROC.jpeg"),p,height=6,width=14,dpi=300)
      ###################
      #
      #FDR&Power
      #
      ###################
      load(paste0(new_res_dir,"shared_50_sub_shared_FDR_Power_shared.RData"))
      p<-FDR_Power_shared_fun(FDR_Power)+ylim(0,0.9)
      ggsave(paste0(fig_dir,"Shared_50_sub_shared_FDR_Power.jpeg"),p,height=6,width=14,dpi=300)
      ##########################################################
      #
      # Ancestry Signal 
      #  ROC | FDR Power 
      #
      ##########################################################
      ###################
      #
      #ROC
      #
      ###################
      load(paste0(new_res_dir,"shared_50_sub_shared_all_ROC_data_dataframe_ancestry.RData"))
      p<-ROC_ancestry_fun(all_ROC_data_dataframe,"External LD")
      ggsave(paste0(fig_dir,"shared_50_sub_ancestry_ROC.jpeg"),p,height=6,width=14,dpi=300)
      ###################
      #
      #FDR&Power
      #
      ###################
      load(paste0(new_res_dir,"shared_50_sub_ancestry_FDR_Power_ancestry.RData"))
      p<-FDR_Power_ancestry_fun(FDR_Power)
      ggsave(paste0(fig_dir,"Shared_50_sub_ancestry_FDR_Power.jpeg"),p,height=6,width=14,dpi=300)
########################################################################
#
#           50% causal variants are shared with different sample size 
#
########################################################################   
      ##########################################################
      #
      # Shared Signal 
      #  ROC | FDR Power | PIP calibration | PIP scatter plot 
      #
      ##########################################################
      ###################
      #
      #ROC
      #
      ###################
      load(paste0(new_res_dir,"shared_50_lowcor_ROC_data_dataframe_shared.RData"))
      p<-ROC_shared_fun(all_ROC_data_dataframe,"50% Correlation")
      ggsave(paste0(fig_dir,"shared_50_lowcor_shared_ROC.jpeg"),p,height=6,width=14,dpi=300)
         ##########################################################
      #
      # Ancestry Signal 
      #  ROC
      #
      ##########################################################
      ###################
      #
      #ROC
      #
      ###################
      load(paste0(new_res_dir,"shared_50_lowcor_ancestry_all_ROC_data_dataframe_ancestry.RData"))
      p<-ROC_ancestry_fun(all_ROC_data_dataframe,"50% Correlation")
      ggsave(paste0(fig_dir,"shared_50_lowcor_ancestry_ROC.jpeg"),p,height=6,width=14,dpi=300)
      