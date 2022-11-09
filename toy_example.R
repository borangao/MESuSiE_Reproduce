
library(data.table)
library(snpStats)
geno_dir<-paste0("/net/fantasia/home/borang/Genotype/Fig1/")

WB_plink<-read.plink(paste0(geno_dir,"WB.bed"))
BB_plink<-read.plink(paste0(geno_dir,"BB.bed"))
WB_plink_geno<-as(WB_plink$genotypes,"numeric")
WB_plink_geno<-apply(WB_plink_geno,2,function(x){
  x[is.na(x)]=mean(x,na.rm=T)
  x = scale(x)
  return(x)
})
WB_plink$map$snp.name
######## SNPs used for reproduction "rs74063479"  "rs113274765" "rs3122053"   "rs3120831" "rs3008244"   "rs3008245"  
BB_plink_geno<-as(BB_plink$genotypes,"numeric")
BB_plink_geno<-apply(BB_plink_geno,2,function(x){
  x[is.na(x)]=mean(x,na.rm=T)
  x = scale(x)
  return(x)
})

BB_plink_geno[,BB_plink$map$allele.1!=WB_plink$map$allele.1]<-1*BB_plink_geno[,BB_plink$map$allele.1!=WB_plink$map$allele.1]

WB_plink_geno<-WB_plink_geno[,c(3,4,1,2,5,6)]
BB_plink_geno<-BB_plink_geno[,c(3,4,1,2,5,6)]
WB_cov<-cov2cor(crossprod(WB_plink_geno))
BB_cov<-cov2cor(crossprod(BB_plink_geno))


colnames(WB_cov)<-paste0("SNP",seq(1,6,1))
colnames(BB_cov)<-paste0("SNP",seq(1,6,1))

colnames(WB_plink_geno)<-paste0("SNP",seq(1,6,1))
colnames(BB_plink_geno)<-paste0("SNP",seq(1,6,1))


causal_snp_list_1<-c("SNP2","SNP4")
causal_snp_list_2<-c("SNP4","SNP6")


num_causal_SNP=2
library(mvtnorm)
WB_plink_causal<-matrix(WB_plink_geno[colnames(WB_plink_geno)%in%causal_snp_list_1],ncol=num_causal_SNP)
BB_plink_causal<-matrix(BB_plink_geno[colnames(BB_plink_geno)%in%causal_snp_list_2],ncol=num_causal_SNP)

beta_BB<-c(0.02,0.02)
beta_BB_all<-rep(0,ncol(BB_cov))
beta_BB_all[which(colnames(BB_cov)%in%causal_snp_list_1)]<-beta_BB
beta_BB_marginal<-BB_cov%*%beta_BB_all


beta_WB<-c(0.02,0.02)
beta_WB_all<-rep(0,ncol(WB_cov))
beta_WB_all[which(colnames(WB_cov)%in%causal_snp_list_2)]<-beta_WB
beta_WB_marginal<-WB_cov%*%beta_WB_all

set.seed(881130)
y_null_WB<-rnorm(nrow(WB_plink_geno),0,sqrt(1-var(WB_cov%*%beta_WB_all)))
y_null_WB<-y_null_WB-mean(y_null_WB)
y_null_BB<-rnorm(nrow(BB_plink_geno),0,sqrt(1-var(BB_cov%*%beta_BB_all)))
y_null_BB<-y_null_BB-mean(y_null_BB)
err_beta_WB<-t(WB_plink_geno)%*%y_null_WB/nrow(WB_plink_geno)
err_beta_BB<-t(BB_plink_geno)%*%y_null_BB/nrow(BB_plink_geno)

target_WB_N<-300000
target_BB_N<-300000


err_beta_WB_scale<-sqrt(nrow(WB_plink_geno)/target_WB_N)*err_beta_WB
err_beta_BB_scale<-sqrt(nrow(BB_plink_geno)/target_BB_N)*err_beta_BB

z_WB<-(beta_WB_marginal+err_beta_WB_scale)*sqrt(target_WB_N)
z_BB<-(beta_BB_marginal+err_beta_BB_scale)*sqrt(target_BB_N)


R_mat_list=list("WB" = WB_cov,"BB" = BB_cov)

summary_stat_1 = data.frame("SNP" = colnames(WB_cov), "Beta"=beta_WB_marginal+err_beta_WB_scale,"Se"=1/sqrt(target_WB_N), "Z" =z_WB,  "N" =target_WB_N )
summary_stat_2 = data.frame("SNP" = colnames(WB_cov), "Beta"=beta_BB_marginal+err_beta_BB_scale,"Se"=1/sqrt(target_BB_N), "Z" =z_BB,  "N" =target_BB_N )
summary_stat_sd_list = list("WB" = summary_stat_1,"BB"=summary_stat_2 )  
save(R_mat_list,summary_stat_sd_list,file="/net/fantasia/home/borang/Genotype/Fig1/Fig1.RData")
source("/net/fantasia/home/borang/Susie_Mult/meSuSie_ancestry/meSuSie_R6.R")
source("/net/fantasia/home/borang/Susie_Mult/meSuSie_ancestry/meSuSie_utility.R")
library(Rcpp)
sourceCpp("/net/fantasia/home/borang/Susie_Mult/meSuSie_ancestry/meSuSie_R6.cpp")
bayes_fac = 3
ancestry_weight=c(bayes_fac/(bayes_fac*2+1),bayes_fac/(bayes_fac*2+1),1/(bayes_fac*2+1))

susie_multi_unique<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =T,max_iter =100)
library(susieR)
susie_WB<-susie_rss(z_WB,WB_cov)
susie_BB<-susie_rss(z_BB,BB_cov)


wrk_dir<-geno_dir
data_dir<-paste0(wrk_dir,"summary_data/")
result_dir<-paste0(wrk_dir,"result/")
system(paste0("mkdir -p ",result_dir))
system(paste0("mkdir -p ",data_dir))
out_dir<-paste0(wrk_dir,"out/")
paintor_dir<-paste0(data_dir,"Paintor/")

system(paste0("mkdir -p ",out_dir))
system(paste0("mkdir -p ",paintor_dir))
z_paintor_name<-paste0(paintor_dir,"fig1")

paintor_z<-data.frame("CHR" = WB_plink$map$chromosome[1],"POS" =WB_plink$map$position ,"RSID" =paste0("SNP",seq(1,6,1)) ,"zscore_1"=z_WB,"zscore_2"=z_BB)

write.table(paintor_z,z_paintor_name,col.names = T,row.names = F,quote=F,sep=" ")

annotation_file<-matrix(rep(1,nrow(WB_plink_geno)),ncol=1)
colnames(annotation_file)<-"coding"
annotation_paintor_name<-paste0(paintor_dir,"fig1.annotations")
write.table(annotation_file,annotation_paintor_name,col.names = T,row.names = F,quote=F)

input_file_name<-paste0(paintor_dir,"fig1.input")
write.table(paste0("fig1"),input_file_name,col.names = F,row.names = F,quote=F)

ld_WB_paintor_name<-paste0(paintor_dir,"fig1.LD1")
ld_BB_paintor_name<-paste0(paintor_dir,"fig1.LD2")
write.table(WB_cov,ld_WB_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")
write.table(BB_cov,ld_BB_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")

system(paste0("/net/fantasia/home/borang/software/fine_mapping/PAINTOR_V3.0/PAINTOR -input ",input_file_name," -Zhead  zscore_1,zscore_2 -LDname LD1,LD2 -in ",input_dir," -out ",result_dir," -mcmc  -annotations coding -Gname ",paintor_suffix," -Lname ",paste0(paintor_suffix,"_bayes_factor")," -RESname ",paste0("mcmc.paintor")))


#############################################################
#
#
#
#					meSuSie
#
#
#
###############################################################
library(ggplot2)
library(cowplot)
Signal_vec<-paste0(c("None","BB","None","Shared Signal","None","WB"))
susie_plot_data<-data.frame(SNP=paste0("SNP",seq(1,6,1)),PIP=susie_multi_unique$pip,Signal =Signal_vec )
susie_plot_data$Signal<-factor(susie_plot_data$Signal)
levels(susie_plot_data$Signal)<-c(levels(susie_plot_data$Signal),"Either")

locuszoom_susie<-ggplot(susie_plot_data,aes(x=SNP,y=PIP,color=Signal))+theme_bw()+geom_point(size =3)+scale_color_manual("Signal",values=c("WB"="#F2C1B6","BB"="#ABDB9F","Shared Signal" = "#6162B0","Either" = "#DB6C79","None"="gray"))+geom_text(aes(label=ifelse(Signal!="None",SNP,"")),hjust=1.5,vjust=0)+labs(x="MESuSiE",y="PIP")+theme(legend.title = element_text(size = 20,face = "bold"),
                                                                                                                                                                                                                                                                                                                                                legend.text = element_text(size = 10))+scale_y_continuous( limits=c(0, 1.1))+theme(legend.title = element_text(size=24))+  theme(legend.text = element_text(size=14))+geom_hline(yintercept=0.5, linetype="dashed")
locuszoom_susie<-locuszoom_susie+theme(axis.text.x = element_text( size = 12),
                                       axis.text.y = element_text( size = 12),  
                                       axis.title.x = element_text( size = 18,face="bold"),
                                       axis.title.y = element_text( size = 18,face="bold"))

plot_legend<-get_legend(locuszoom_susie)
locuszoom_susie<-locuszoom_susie+theme(legend.position="none")
#ggsave("/net/fantasia/home/borang/Susie_Mult/csg_present/fig/ethnic_unique_susie_locuszoom.jpeg",locuszoom_susie,dpi=300)

#############################################################
#
#
#
#					SuSie WB
#
#
#
###############################################################
library(ggplot2)
library(cowplot)
susie_WB$pip>0.5

result = rep(0,length(susie_WB$pip))
result[which(susie_WB$pip>0.5)]<-1
susie_signal<-rep(0,6)
susie_signal[which(result==1)]<-1
susie_signal[susie_signal==0]<-"Non"
susie_signal[susie_signal==1]<-"WB"
paintor_data<-data.frame("SNP"=paste0("SNP",seq(1,6,1)),"PIP"=susie_WB$pip,"Signal"=susie_signal,"RSID"=paste0("SNP",seq(1,6,1)))


locuszoom_susie_WB<-ggplot(paintor_data,aes(x=SNP,y=PIP,color=Signal))+theme_bw()+geom_point(size =3)+scale_y_continuous( limits=c(0, 1.1))+scale_color_manual("",values=c("Non"="gray","WB"="#F2C1B6"))+geom_text(aes(label=ifelse(Signal!="Non",RSID,"")),hjust=1.5,vjust=0) +labs(x="SuSiE WB",y="PIP")+theme(legend.title = element_text(size = 24,face = "bold"),
                                                                                                                                                                                                                                                                                                                 legend.text = element_text(size = 14))+ guides(fill=guide_legend(title=""))+geom_hline(yintercept=0.5, linetype="dashed")

locuszoom_susie_WB<-locuszoom_susie_WB+theme(legend.position="none")+theme(axis.text.x = element_text( size = 12),
                                                                           axis.text.y = element_text( size = 12),  
                                                                           axis.title.x = element_text( size = 18,face="bold"),
                                                                           axis.title.y = element_text( size = 18,face="bold"))							
#ggsave("/net/fantasia/home/borang/Susie_Mult/csg_present/fig/susie_WB_locuszoom.jpeg",locuszoom_susie_WB,dpi=300)
#############################################################
#
#
#
#					SuSie BB
#
#
#
###############################################################
library(ggplot2)
library(cowplot)
susie_BB$pip>0.5

result = rep(0,length(susie_BB$pip))
result[which(susie_BB$pip>0.5)]<-1
susie_signal<-rep(0,6)
susie_signal[which(result==1)]<-1
susie_signal[susie_signal==0]<-"Non"
susie_signal[susie_signal==1]<-"BB"
paintor_data<-data.frame("SNP"=paste0("SNP",seq(1,6,1)),"PIP"=susie_BB$pip,"Signal"=susie_signal,"RSID"=paste0("SNP",seq(1,6,1)))


locuszoom_susie_BB<-ggplot(paintor_data,aes(x=SNP,y=PIP,color=Signal))+theme_bw()+geom_point(size =3)+scale_y_continuous( limits=c(0, 1.1))+scale_color_manual("",values=c("Non"="gray","BB"="#ABDB9F"))+geom_text(aes(label=ifelse(Signal!="Non",RSID,"")),hjust=1.5,vjust=0) +labs(x="SuSiE BB",y="PIP")+theme(legend.title = element_text(size = 24,face = "bold"),
                                                                                                                                                                                                                                                                                                                  legend.text = element_text(size = 14))+ guides(fill=guide_legend(title=""))+geom_hline(yintercept=0.5, linetype="dashed")

locuszoom_susie_BB<-locuszoom_susie_BB+theme(legend.position="none")+theme(axis.text.x = element_text( size = 12),
                                                                             axis.text.y = element_text( size = 12),  
                                                                             axis.title.x = element_text( size = 18,face="bold"),
                                                                             axis.title.y = element_text( size = 18,face="bold"))							
#ggsave("/net/fantasia/home/borang/Susie_Mult/csg_present/fig/susie_BB_locuszoom.jpeg",locuszoom_susie_BB,dpi=300)


#############################################################
#
#
#
#					Paintor
#
#
#
###############################################################
paintor_post<-read.table(paste0("/net/fantasia/home/borang/Susie_Mult/csg_present/result/fig1.mcmc.paintor"),header=T)
result = rep(0,length(paintor_post$Posterior_Prob))
result[which(paintor_post$Posterior_Prob>0.5)]<-1
paintor_signal<-rep(0,6)
paintor_signal[which(result==1)]<-1
paintor_signal[paintor_signal==0]<-"Non"
paintor_signal[paintor_signal==1]<-"Paintor"

paintor_data<-data.frame("SNP"=paste0("SNP",seq(1,6,1)),"PIP"=paintor_post$Posterior_Prob,"Signal"=paintor_signal,"RSID"=paste0("SNP",seq(1,6,1)))


locuszoom_paintor<-ggplot(paintor_data,aes(x=SNP,y=PIP,color=Signal))+theme_bw()+geom_point(size =3)+scale_y_continuous( limits=c(0, 1.1))+scale_color_manual("",values=c("Non"="lightgrey","Paintor"="#DB6C79"))+geom_text(aes(label=ifelse(Signal!="Non",RSID,"")),hjust=1.5,vjust=0) +labs(x="Paintor",y="PIP")+theme(legend.title = element_text(size = 24,face = "bold"),
                                                                                                                                                                                                                                                                                                                         legend.text = element_text(size = 14))+ guides(fill=guide_legend(title=""))+geom_hline(yintercept=0.5, linetype="dashed")
paintor_legend<-get_legend(locuszoom_paintor)						
locuszoom_paintor<-locuszoom_paintor+theme(legend.position="none")+theme(axis.text.x = element_text( size = 12),
                                                                         axis.text.y = element_text( size = 12),  
                                                                         axis.title.x = element_text( size = 18,face="bold"),
                                                                         axis.title.y = element_text( size = 18,face="bold"))								
#ggsave("/net/fantasia/home/borang/Susie_Mult/csg_present/fig/ethnic_unique_paintor_locuszoom.jpeg",locuszoom_paintor,dpi=300)

library(cowplot)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(grid)
library(snpStats)
library(ggrepel)

pip_plot<-grid.arrange(grid.arrange(locuszoom_susie_WB,locuszoom_susie,locuszoom_susie_BB,locuszoom_paintor,nrow =2),plot_legend,ncol=2,widths=c(9,3))
#ggsave("/net/fantasia/home/borang/Susie_Mult/csg_present/fig/all_method_locuszoom.jpeg",pip_plot,dpi=300)

#############################################################
#
#
#
#					Barplot
#
#
#
###############################################################

# WB_bar_data <- data.frame(SNP=paste0("SNP",seq(1,6,1)),Effect = beta_WB_all)
#  
# # Grouped
# #"WB"="#fadf63","BB"="#99d5c9"
# bar_WB<-ggplot(WB_bar_data, aes( y=Effect, x=SNP)) + 
#     geom_bar(position="dodge", stat="identity",fill = "#fadf63")+theme_bw()
# 	
# 
# ggsave("/net/fantasia/home/borang/Susie_Mult/simulation/example_fig/fig/WB_bar.jpeg",bar_WB,dpi=300)
# BB_bar_data <- data.frame(SNP=paste0("SNP",seq(1,6,1)),Effect = beta_BB_all)
#  
# # Grouped
# bar_BB<-ggplot(BB_bar_data, aes( y=Effect, x=SNP)) + 
#     geom_bar(position="dodge", stat="identity",fill ="#99d5c9")+theme_bw()


#ggsave("/net/fantasia/home/borang/Susie_Mult/csg_present/fig/BB_bar.jpeg",bar_BB,dpi=300)


WBR_bar_data <- data.frame(SNP=paste0("SNP",seq(1,6,1)),Effect = beta_WB_all,Signal =c("None","BB","None","Shared","None","WBR"))
BB_bar_data <- data.frame(SNP=paste0("SNP",seq(1,6,1)),Effect = beta_BB_all,Signal =c("None","BB","None","Shared","None","WBR"))

# Grouped
bar_WB<-ggplot(WBR_bar_data, aes( y=Effect, x=SNP,fill = Signal)) + xlab("")+ylab("Effect")+
  geom_bar(position="dodge", stat="identity")+theme_bw()+scale_fill_manual("Signal",values=c("WBR"="#F2C1B6","BB"="#ABDB9F","Shared" = "#6162B0","None"="gray"))+theme(legend.position="none")
bar_WB = bar_WB + 	theme(axis.text.x =  element_text( size = 8),
                         axis.text.y = element_text( size = 8))	+theme(axis.text.x = element_text( size = 10),
                                                                       axis.text.y = element_text( size = 10),  
                                                                       axis.title.x = element_text( size = 18,face="bold"),
                                                                       axis.title.y = element_text( size = 18,face="bold"))
bar_BB<-ggplot(BB_bar_data, aes( y=Effect, x=SNP,fill = Signal)) + xlab("")+ylab("Effect")+
  geom_bar(position="dodge", stat="identity")+theme_bw()+scale_fill_manual("Signal",values=c("WBR"="#F2C1B6","BB"="#ABDB9F","Shared" = "#6162B0","None"="gray"))+theme(legend.position="none")
bar_BB = bar_BB + 	theme(axis.text.x =  element_text( size = 8),
                           axis.text.y = element_text( size = 8))	+theme(axis.text.x = element_text( size = 10),
                                                                         axis.text.y = element_text( size = 10),  
                                                                         axis.title.x = element_text( size = 18,face="bold"),
                                                                         axis.title.y = element_text( size = 18,face="bold"))


#############################################################
#
#
#
#					LD plot
#
#
#
###############################################################

melt_WB_cov<-matrix(ncol=3,nrow=0)
melt_WB_cov<-data.frame(melt_WB_cov)

melt_BB_cov<-matrix(ncol=3,nrow=0)
melt_BB_cov<-data.frame(melt_BB_cov)
for(SNP1 in 1:6){
  for(SNP2 in SNP1:6){
    melt_WB_cov<-rbind(melt_WB_cov,c(paste0("SNP",SNP1),paste0("SNP",SNP2),WB_cov[SNP1,SNP2]))
    melt_BB_cov<-rbind(melt_BB_cov,c(paste0("SNP",SNP1),paste0("SNP",SNP2),BB_cov[SNP1,SNP2]))
    
    
  }
  
  
  
}
colnames(melt_WB_cov)<-c("SNP1","SNP2","Value")
melt_WB_cov$Value<-as.numeric(melt_WB_cov$Value)
colnames(melt_BB_cov)<-c("SNP1","SNP2","Value")
melt_BB_cov$Value<-as.numeric(melt_BB_cov$Value)

#ggcorrplot(round(WB_cov,2),colors = c("lightblue", "white", "orange"),type ="upper", lab =TRUE,show.diag=T)
#ggcorrplot(round(BB_cov,2),colors = c("lightblue", "white", "orange"),type ="upper", lab =TRUE,show.diag=T)

WB_cor<-ggplot(melt_WB_cov, aes(x = SNP1, y = SNP2, fill = Value)) +
  geom_tile(aes( fill = Value))+scale_fill_stepsn(
    colors = c( "purple","blue","lightskyblue", "green", "orange", "red"),
    breaks = seq(0, 0.8, by = 0.2),
    limits = c(-0.2, 1),
    show.limits = TRUE,
    na.value = 'grey50',
    name = ""
  )+theme_bw() +
  theme(legend.position = "none")+xlab("WB")+ylab("LD")
WB_cor<-WB_cor+geom_text(aes(label=round(Value,2)))+theme(axis.text.x = element_text( size = 10),
                                                          axis.text.y = element_text( size = 10),  
                                                          axis.title.x = element_text( size = 18,face="bold"),
                                                          axis.title.y = element_text( size = 18,face="bold"))

BB_cor<-ggplot(melt_BB_cov, aes(x = SNP1, y = SNP2, fill = Value)) +
  geom_tile(aes( fill = Value))+scale_fill_stepsn(
    colors = c( "purple","blue","lightskyblue", "green", "orange", "red"),
    breaks = seq(0, 0.8, by = 0.2),
    limits = c(-0.2, 1),
    show.limits = TRUE,
    na.value = 'grey50',
    name = ""
  )+theme_bw() +
  theme(legend.position = "none")+xlab("BB")+ylab("")
BB_cor<-BB_cor+geom_text(aes(label=round(Value,2)))+theme(axis.text.x = element_text( size = 10),
                                                            axis.text.y = element_text( size = 10),  
                                                            axis.title.x = element_text( size = 18,face="bold"),
                                                            axis.title.y = element_text( size = 18,face="bold"))
#BB_cor = BB_cor +  theme(legend.position = c(0.9,0.35),legend.background = element_rect(fill = "transparent") )
#############################################################
#
#
#
#					GWAS plot
#
#
#
###############################################################

gwas_data<-data.frame("SNP"=paste0("SNP",seq(1,6,1)),"Pvalue"=z_WB,"RSID"=paste0("SNP",seq(1,6,1)))

locuszoom_gwas_WB<-ggplot(gwas_data,aes(x=SNP,y=-1*log10(2*pnorm(-1*abs(Pvalue)))))+theme_bw()+geom_point(size =3)+geom_hline(yintercept=8, linetype="dashed", 
                                                                                                                              color = "red", size=2)+labs(x="",y="-log10 Pvalue") +theme(legend.title = element_text(size = 24,face = "bold"),
                                                                                                                                                                                         legend.text = element_text(size = 20))+theme(axis.text.x = element_text( size = 10),
                                                                                                                                                                                                                                      axis.text.y = element_text( size = 10),  
                                                                                                                                                                                                                                      axis.title.x = element_text( size = 18,face="bold"),
                                                                                                                                                                                                                                      axis.title.y = element_text( size = 18,face="bold"))
gwas_data<-data.frame("SNP"=paste0("SNP",seq(1,6,1)),"Pvalue"=z_BB,"RSID"=paste0("SNP",seq(1,6,1)))

locuszoom_gwas_BB<-ggplot(gwas_data,aes(x=SNP,y=-1*log10(2*pnorm(-1*abs(Pvalue)))))+theme_bw()+geom_point(size =3)+geom_hline(yintercept=8, linetype="dashed", 
                                                                                                                               color = "red", size=2)+labs(x="",y="-log10 Pvalue") +theme(legend.title = element_text(size = 24,face = "bold"),
                                                                                                                                                                                          legend.text = element_text(size = 20))+theme(axis.text.x = element_text( size = 10),
                                                                                                                                                                                                                                       axis.text.y = element_text( size = 10),  
                                                                                                                                                                                                                                       axis.title.x = element_text( size = 18,face="bold"),
                                                                                                                                                                                                                                       axis.title.y = element_text( size = 18,face="bold"))


library(ggpubr)
beta_gwas_test<-grid.arrange(bar_WB,bar_BB,ncol=2,nrow=1)
locuszoom_test<-grid.arrange(locuszoom_gwas_WB,locuszoom_gwas_BB,ncol=2,nrow=1)
LD_test<-grid.arrange(WB_cor,BB_cor,ncol=2,nrow=1)
labeled_plot<-ggarrange(beta_gwas_test, locuszoom_test, LD_test, labels = c("A", "B", "C"),  
                        nrow=3, ncol=1)
#ggsave("/net/fantasia/home/borang/Susie_Mult/csg_present/fig/annotate_figure.jpeg",labeled_plot,dpi=300)

pip_plot<-grid.arrange(grid.arrange(locuszoom_susie_WB,locuszoom_susie_BB,nrow =2),grid.arrange(locuszoom_susie,locuszoom_paintor,nrow =2),plot_legend,ncol=3,widths=c(5,5,2.5))
pip_plot_test<-ggarrange(pip_plot, labels = c("D"))
Fig_1<-grid.arrange(labeled_plot,pip_plot_test,ncol=2,widths=c(4,6))
#ggsave("/net/fantasia/home/borang/Susie_Mult/csg_present/fig/Fig_1.jpeg",Fig_1,dpi=300,height = 6, width =14)
ggsave("/net/fantasia/home/borang/Susie_Mult/Paper_Figure/figure/Fig_1.jpeg",Fig_1,dpi=300,height = 9, width =18)






########################################
#
#         Save RData file for website
#
########################################


