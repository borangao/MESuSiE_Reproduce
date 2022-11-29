library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(VennDiagram)
library(gridExtra)
library(DescTools)
load("/net/fantasia/home/borang/Susie_Mult/website_material/real_data/res_all.RData")
#res_all_subset<-res_all%>%filter(SuSiE_Either>0.5|MESuSiE_PIP_Shared>0.5|Paintor_PIP>0.5|MESuSiE_PIP_WB>0.5|MESuSiE_PIP_BB>0.5)
#write.csv(res_all_subset,"/net/fantasia/home/borang/Susie_Mult/Paper_Figure/summarized_data_new/Table_S2.csv", row.names = FALSE,quote=F)
################################################
#
#		Set Size/Z-score/eQTL
#
#
###############################################
				################################################
				#
				#		Set SiZe Part
				#
				#
				###############################################		
				###Summary Stat of Set SiZe
				res_all%>%group_by(Region) %>% summarise(across(c("MESuSiE_cs", "SuSiE_cs","Paintor_cs"), ~ sum(.x, na.rm = TRUE)))%>% summarise(median(MESuSiE_cs),median(SuSiE_cs),median(Paintor_cs))
				###Median set size by Trait
				all_sets_info<-data.frame(res_all%>%group_by(Trait,Region) %>% summarise(across(c("MESuSiE_cs", "SuSiE_cs","Paintor_cs"), ~ sum(.x, na.rm = TRUE)))) ###Median Set Size across all locus

				all_sets_info_long<-all_sets_info%>%pivot_longer(!(Trait|Region), names_to = "Method", values_to = "Count")
				all_sets_info_long$Method<-factor(all_sets_info_long$Method,levels=c("MESuSiE_cs","SuSiE_cs","Paintor_cs"))
				levels(all_sets_info_long$Method)<-c("MESuSiE","SuSiE","Paintor")

				p_set = ggplot(data =all_sets_info_long,aes(x = Trait, y=Count,fill=Method))+geom_boxplot(aes(x = Trait,middle = mean(Count,na.rm=T),fill=Method))+scale_y_continuous(limits = quantile(all_sets_info_long$Count, c(0, 0.9),na.rm =T))+scale_fill_manual(values=c("MESuSiE"="#023e8a","SuSiE"="#2a9d8f","Paintor"="#f4a261"))
				p_set =p_set + theme_bw() + xlab("") +ylab("Set Size")
				p_set= p_set+ theme(axis.text.x = element_text( size = 16),
								axis.text.y = element_text( size = 16),  
								axis.title.x = element_text( size = 18),
								axis.title.y = element_text( size = 18,face="bold"),
								strip.text.x = element_text(size = 18),
								strip.text.y= element_text(size = 18),
								legend.text=element_text(size=16),
								legend.title=element_text(size=18,face="bold"),
								plot.title = element_text(size=16,hjust = 0.5))
				p_set = p_set + theme(legend.position = "none")				
				################################################
				#
				#		Z-score Part
				#
				#
				###############################################		
					###Summary Stat of Z-scores
					res_all%>%filter(MESuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))
					res_all%>%filter(SuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))
					res_all%>%filter(Paintor_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))
					
					MESuSiE_cs_Z<-res_all%>%group_by(Trait) %>%filter(MESuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))
					SuSiE_cs_Z<-res_all%>%group_by(Trait) %>%filter(SuSiE_cs==1)%>%summarise(zmax =median(pmax(abs(zscore_WB),abs(zscore_BB))))
					Paintor_cs_Z<-res_all%>%group_by(Trait) %>%filter(Paintor_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))
					set_size_z_info<-data.frame(cbind(MESuSiE_cs_Z,SuSiE_cs_Z$zmax,Paintor_cs_Z$zmax))
					colnames(set_size_z_info)<-c("Trait",c("MESuSiE_Z","SuSiE_Z","Paintor_Z"))

					set_size_z_info_long<-set_size_z_info %>%pivot_longer(!(Trait), names_to = "Method", values_to = "Z")
					set_size_z_info_long$Method<-factor(set_size_z_info_long$Method)
					levels(set_size_z_info_long$Method)<-c("MESuSiE","SuSiE","Paintor")

					p_z = ggplot(data = set_size_z_info_long,aes(x = Trait, y=Z,fill=Method))+geom_bar( stat = "identity",position="dodge")+scale_fill_manual(values=c("MESuSiE"="#023e8a","SuSiE"="#2a9d8f","Paintor"="#f4a261"))
					p_z = p_z + theme_bw() + xlab("") +ylab("Median |Z|")
					p_z = p_z + theme(axis.text.x = element_text( size = 16),
									axis.text.y = element_text( size = 16),  
									axis.title.x = element_text( size = 18),
									axis.title.y = element_text( size = 18,face="bold"),
									strip.text.x = element_text(size = 18),
									strip.text.y= element_text(size = 18),
									legend.text=element_text(size=18),
									legend.title=element_text(size=22,face="bold"),
									plot.title = element_text(size=16,hjust = 0.5))
					p_z = p_z + theme(legend.position="none")+ theme(strip.background = element_blank(),strip.text.y = element_blank())
				################################################
				#
				#		eQTL Part
				#
				#
				###############################################	
				res_all$utr_comb = res_all$utr_5+res_all$utr_3
				ann_col_name<-c("pLOF","missense", "synonymous", "utr_comb", "promotor", "CRE", "eQTL")
				MESuSiE_PIP_ann = res_all%>%group_by(Trait,MESuSiE_cs)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))%>%group_by(Trait)%>%summarise(across(ann_col_name,~.x[MESuSiE_cs==1]/.x[MESuSiE_cs==0]))
				SuSiE_PIP_ann = res_all%>%group_by(Trait,SuSiE_cs)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))%>%group_by(Trait)%>%summarise(across(ann_col_name,~.x[SuSiE_cs==1]/.x[SuSiE_cs==0]))
				Paintor_PIP_ann = res_all%>%group_by(Trait,Paintor_cs)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))%>%group_by(Trait)%>%summarise(across(ann_col_name,~.x[Paintor_cs==1]/.x[Paintor_cs==0]))

				Either_Ancestry_cs_Ann<-rbind(MESuSiE_PIP_ann,SuSiE_PIP_ann,Paintor_PIP_ann)			
				Either_Ancestry_cs_Ann<-data.frame(c(rep("MESuSiE",4),rep("SuSiE",4),rep("Paintor",4)),Either_Ancestry_cs_Ann)
				colnames(Either_Ancestry_cs_Ann)[1]<-"Method"
				Either_Ancestry_cs_Ann_long<-Either_Ancestry_cs_Ann %>%pivot_longer(!c(Method,Trait), names_to = "Cat", values_to = "Prop")
				Either_Ancestry_cs_Ann_long_subset<-Either_Ancestry_cs_Ann_long[ Either_Ancestry_cs_Ann_long$Cat=="eQTL",]
				Either_Ancestry_cs_Ann_long$Method<-factor(Either_Ancestry_cs_Ann_long$Method,levels=c("MESuSiE","SuSiE","Paintor"))
				
				Either_Ancestry_cs_Ann_long_subset_eQTL<- Either_Ancestry_cs_Ann_long[ Either_Ancestry_cs_Ann_long$Cat=="eQTL",]
				p_eQTL = ggplot(data = Either_Ancestry_cs_Ann_long_subset_eQTL,aes(x = Trait, y=Prop,fill=Method))+geom_bar( stat = "identity",position="dodge")+scale_fill_manual(values=c("MESuSiE"="#023e8a","SuSiE"="#2a9d8f","Paintor"="#f4a261"))
				p_eQTL = p_eQTL+ theme_bw() + xlab("") +ylab("eQTL Fold Enrichment")
				p_eQTL = p_eQTL + theme(axis.text.x = element_text( size = 16),
								axis.text.y = element_text( size = 16),  
								axis.title.x = element_text( size = 18),
								axis.title.y = element_text( size = 18,face="bold"),
								strip.text.x = element_text(size = 18),
								strip.text.y= element_text(size = 18),
								legend.text=element_text(size=14),
								legend.title=element_text(size=16,face="bold"),
								plot.title = element_text(size=16,hjust = 0.5))
				p_eQTL = p_eQTL + theme(strip.background = element_blank(),strip.text.y = element_blank())+theme(legend.position="bottom")

				p_out<-p_set/p_z/p_eQTL

				ggsave("/net/fantasia/home/borang/Susie_Mult/real_data/summary_data_lipid/070722/Figure/real_data_sets_info.jpeg",p_out,dpi=300,width=6,height =12)
				

###############################################################
#
#         Venn diagram for 95% credible set SNPs
#
###############################################################
				for(trait_name in c("HDL","LDL","TC","TG")){
				  MESuSiE_SNP = res_all%>%filter(Trait ==trait_name,MESuSiE_cs==1)%>%dplyr::select(SNP)%>%pull(SNP)
				  SuSiE_SNP = res_all%>%filter(Trait ==trait_name,SuSiE_cs==1)%>%dplyr::select(SNP)%>%pull(SNP)
				  Paintor_SNP = res_all%>%filter(Trait ==trait_name,Paintor_cs==1)%>%dplyr::select(SNP)%>%pull(SNP)
				  list1 <- list("MESuSiE"=MESuSiE_SNP, "SuSiE"=SuSiE_SNP,"Paintor" = Paintor_SNP )
				  venn1 <- venn.diagram(list1, filename = NULL,print.mode=c("raw","percent"),col=c("#023e8a", "#2a9d8f", "#f4a261"),
				                        fill = c(alpha("#023e8a",0.6), alpha("#2a9d8f",0.6), alpha("#f4a261",0.6)),main = trait_name,main.fontface="bold",main.cex=2)
				  assign(paste0("venn_",trait_name),venn1)
				}
				
				venn_out<-grid.arrange(gTree(children=venn_HDL),gTree(children=venn_LDL),gTree(children=venn_TC),gTree(children=venn_TG),ncol=2)
################################################################
#
#
#     Functional Annotation enrichment for 95% credible set SNPS
#
#
###############################################################
				
				ann_col_name<-c("pLOF","missense", "synonymous", "utr_comb", "promotor", "CRE", "eQTL")
				MESuSiE_PIP_ann = res_all%>%group_by(MESuSiE_cs)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))%>%summarise(across(ann_col_name,~.x[MESuSiE_cs==1]/.x[MESuSiE_cs==0]))
				SuSiE_PIP_ann = res_all%>%group_by(SuSiE_cs)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))%>%summarise(across(ann_col_name,~.x[SuSiE_cs==1]/.x[SuSiE_cs==0]))
				Paintor_PIP_ann = res_all%>%group_by(Paintor_cs)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))%>%summarise(across(ann_col_name,~.x[Paintor_cs==1]/.x[Paintor_cs==0]))
				
				
				Either_Ancestry_cs_Ann<-rbind(MESuSiE_PIP_ann,SuSiE_PIP_ann,Paintor_PIP_ann)			
				Either_Ancestry_cs_Ann<-data.frame(Method = c("MESuSiE","SuSiE","Paintor"),Either_Ancestry_cs_Ann)
				Either_Ancestry_cs_Ann_long<-Either_Ancestry_cs_Ann %>%pivot_longer(!c(Method), names_to = "Cat", values_to = "Prop")
				Either_Ancestry_cs_Ann_long$Method<-factor(Either_Ancestry_cs_Ann_long$Method,levels=c("MESuSiE","SuSiE","Paintor"))
				Either_Ancestry_cs_Ann_long_subset<-Either_Ancestry_cs_Ann_long%>%filter(!(Cat=="eQTL"))
				Either_Ancestry_cs_Ann_long_subset$Cat<-factor(Either_Ancestry_cs_Ann_long_subset$Cat,levels = c("pLOF","missense","synonymous","utr_comb","promotor","CRE"))
				levels(Either_Ancestry_cs_Ann_long_subset$Cat)<-c("pLOF","Missense","Synonymous","Utr","Promotor","CRE")
				Either_Ancestry_cs_Ann_long_subset$Prop<-round(Either_Ancestry_cs_Ann_long_subset$Prop,2)
				
				p = ggplot(data = Either_Ancestry_cs_Ann_long_subset,aes(x=Cat,y=Prop,fill=Method))+geom_col( position="dodge")
				p = p + xlab("")+ylab("Fold Enrichment") + theme_bw() + geom_hline(yintercept = 1, linetype = "dashed")+scale_fill_manual(values=c("MESuSiE"="#023e8a","SuSiE"="#2a9d8f","Paintor"="#f4a261"))
				p = p + geom_text(aes(x=Cat,group=Method,y=Prop,label=Prop),position = position_dodge(width = 1),vjust=-0.5,size = 4,fontface="bold")
				p = p + theme(axis.text.x = element_text( size = 16,angle=90),
				              axis.text.y = element_text( size = 16),  
				              axis.title.x = element_text( size = 18),
				              axis.title.y = element_text( size = 18,face="bold"),
				              strip.text.x = element_text(size = 18),
				              strip.text.y= element_text(size = 18),
				              legend.text=element_text(size=18),
				              legend.title=element_text(size=22,face="bold"),
				              plot.title = element_text(size=16,hjust = 0.5))
				
				ggsave("/net/fantasia/home/borang/Susie_Mult/real_data/summary_data_lipid/070722/Figure/real_data_annotation_info_combined.jpeg",p,dpi=300,width=12,height =6)

				
				
############################################################################
#
#
#                       Proportion of Signal Plot
#
#
############################################################################
			res_all_signal_number<-res_all%>%group_by(Trait)%>%summarise(Paintor_Shared_n = sum(Paintor_PIP>0.5),SuSiE_Shared_n = sum(SuSiE_Shared>0.5),SuSiE_WB_n = sum(SuSiE_WB>0.5&SuSiE_BB<0.5),SuSiE_BB_n = sum(SuSiE_WB<0.5&SuSiE_BB>0.5),MESuSiE_Shared_n = sum(MESuSiE_PIP_Shared>0.5),MESuSiE_WB_n = sum(MESuSiE_PIP_WB>0.5),MESuSiE_BB_n = sum(MESuSiE_PIP_BB>0.5))
				for(trait_name in c("HDL","LDL","TC","TG")){
				  res_all_signal_number_subset<-res_all_signal_number[res_all_signal_number$Trait==trait_name,]
				susie_signal_num_data<-data.frame(Cat = c("Shared","EUR","AFR"),
				                                  Num = c(sum(res_all_signal_number_subset$SuSiE_Shared_n),
				                                          sum(res_all_signal_number_subset$SuSiE_WB_n),
				                                          sum(res_all_signal_number_subset$SuSiE_BB_n)))
				
				# Compute the position of labels
				susie_signal_num_data <- susie_signal_num_data %>% 
				  arrange(desc(Cat)) %>%
				  mutate(prop = Num / sum(susie_signal_num_data$Num) *100) %>%
				  mutate(ypos = cumsum(prop)- 0.5*prop )
				susie_signal_num_data$label = paste0(susie_signal_num_data$Cat," ",susie_signal_num_data$Num)
				
				mesusie_signal_num_data<-data.frame(Cat = c("Shared","EUR","AFR"),
				                                    Num = c(sum(res_all_signal_number_subset$MESuSiE_Shared_n),
				                                            sum(res_all_signal_number_subset$MESuSiE_WB_n),
				                                            sum(res_all_signal_number_subset$MESuSiE_BB_n)))
				
				# Compute the position of labels
				mesusie_signal_num_data <- mesusie_signal_num_data %>% 
				  arrange(desc(Cat)) %>%
				  mutate(prop = Num / sum(mesusie_signal_num_data$Num) *100) %>%
				  mutate(ypos = cumsum(prop)- 0.5*prop )
				
				mesusie_signal_num_data$label = paste0(mesusie_signal_num_data$Cat," ",mesusie_signal_num_data$Num)
				
				
				paintor_signal_num_data<-data.frame(Cat = c("Either"),
				                                    Num = c(sum(res_all_signal_number_subset$Paintor_Shared_n)))
				# Compute the position of labels
				paintor_signal_num_data <- paintor_signal_num_data %>% 
				  arrange(desc(Cat)) %>%
				  mutate(prop = Num / sum(paintor_signal_num_data$Num) *100) %>%
				  mutate(ypos = cumsum(prop)- 0.5*prop )
				
				paintor_signal_num_data$label = paste0(paintor_signal_num_data$Cat," ",paintor_signal_num_data$Num)
				
				# Basic piechart
				susie_signal_pie_chart<-ggplot(susie_signal_num_data, aes(x="", y=prop, fill=Cat)) +
				  geom_bar(stat="identity", width=1, color="white") +
				  coord_polar("y", start=0) +
				  theme_void() + 
				  theme(legend.position="none") +
				  geom_text(aes(y = ypos, label = label),  size=4,color="white") +
				  scale_fill_manual(values=c("#ABDB9F","#F2C1B6","#6162B0"))+ggtitle(ifelse(trait_name=="HDL","SuSiE",""))+theme(plot.title = element_text(size = 20,face="bold",hjust = 0.5))
				
				mesusie_signal_pie_chart<-ggplot(mesusie_signal_num_data, aes(x="", y=prop, fill=Cat)) +
				  geom_bar(stat="identity", width=1, color="white") +
				  coord_polar("y", start=0) +
				  theme_void() + 
				  theme(legend.position="none") +
				  geom_text(aes(y = ypos, label = label), size=4,color="white") +
				  scale_fill_manual(values=c("#ABDB9F","#F2C1B6","#6162B0"))+ggtitle(ifelse(trait_name=="HDL","MESuSiE",""))+theme(plot.title = element_text(size = 20,face="bold",hjust = 0.5))
				
				
				paintor_signal_pie_chart<-ggplot(paintor_signal_num_data, aes(x="", y=prop, fill=Cat)) +
				  geom_bar(stat="identity", width=1, color="white") +
				  coord_polar("y", start=0) +
				  theme_void() + 
				  theme(legend.position="none") +
				  geom_text(aes(y = ypos, label = label),  size=4,color="white") +
				  scale_fill_manual(values=c("gray50"))+ggtitle(ifelse(trait_name=="HDL","Paintor",""))+theme(plot.title = element_text(size = 20,face="bold",hjust = 0.5))
				
				p_pie<-(mesusie_signal_pie_chart  + susie_signal_pie_chart+paintor_signal_pie_chart) 
				assign(paste0(trait_name,"_pie"),p_pie)
				}
		
				p_pie_out<-HDL_pie/LDL_pie/TC_pie/TG_pie
				
			#	ggsave(paste0("/net/fantasia/home/borang/Susie_Mult/real_data/summary_data_lipid/070722/Figure_New/signal_number.jpeg"),p_pie_out,width = 8,height =12,dpi=300)
				
###########################################################################################
#
#
#         Correlation of ancestry-specific signal by MESuSiE and SuSiE
#
#
##########################################################################################
				res_all_add<-res_all%>%filter(MESuSiE_PIP_BB>0.5|MESuSiE_PIP_WB>0.5|(SuSiE_WB>0.5&SuSiE_BB<0.5)|(SuSiE_BB>0.5&SuSiE_WB<0.5))
				res_all_add$MESuSiE_SuSiE_Cat=3
				res_all_add$MESuSiE_SuSiE_Cat[(res_all_add$MESuSiE_PIP_BB>0.5|res_all_add$MESuSiE_PIP_WB>0.5)&(res_all_add$SuSiE_WB<0.5&res_all_add$SuSiE_BB<0.5)]<-2
				res_all_add$MESuSiE_SuSiE_Cat[(res_all_add$MESuSiE_PIP_BB<0.5&res_all_add$MESuSiE_PIP_WB<0.5)]<-1
				res_all_add$MESuSiE_SuSiE_Cat<-factor(res_all_add$MESuSiE_SuSiE_Cat)
				levels(res_all_add$MESuSiE_SuSiE_Cat)<-c("SuSiE","MESuSiE","Shared")
				
					
				betaplot<-ggplot(res_all_add, aes(x = Beta_WB, y =  Beta_BB))+geom_smooth(aes(x = Beta_WB, y =  Beta_BB,color = MESuSiE_SuSiE_Cat,shape=MESuSiE_SuSiE_Cat),method="lm",se=FALSE,col="red",formula = y ~ x)+  stat_cor()+scale_shape_manual(values=c(16,17,18))+scale_size_manual(values=c(0.75,1.25,1.75))+scale_colour_manual(values=c("gray", "#2a9d8f","#f4a261"))+xlab("UKBB GWAS")+ylab("AA GWAS")+
				  geom_point(data=res_all_add,aes(x = Beta_WB, y =  Beta_BB, color = MESuSiE_SuSiE_Cat,shape=MESuSiE_SuSiE_Cat,size=MESuSiE_SuSiE_Cat))+theme_bw()+
				  facet_grid(vars(Trait),vars(MESuSiE_SuSiE_Cat))
				betaplot = betaplot + theme(axis.text.x = element_text( size = 16),
				                            axis.text.y = element_text( size = 16),  
				                            axis.title.x = element_text( size = 18),
				                            axis.title.y = element_text( size = 18),
				                            strip.text.x = element_text(size = 18),
				                            strip.text.y= element_text(size = 18),
				                            legend.text=element_text(size=18),
				                            legend.title=element_text(size=22,face="bold"),
				                            plot.title = element_text(size=16,hjust = 0.5))
				betaplot = betaplot + theme(legend.position="none")
				ggsave(paste0("/net/fantasia/home/borang/Susie_Mult/real_data/summary_data_lipid/070722/Figure_New/cor_beta_ancestry_specific.jpeg"),betaplot,width = 12,height =12,dpi=300)
				
			  res_all_add%>%group_by(MESuSiE_SuSiE_Cat,Trait)%>%summarise(cor(Beta_WB,Beta_BB))
				res_all_add%>%group_by(MESuSiE_SuSiE_Cat)%>%summarise(cor(Beta_WB,Beta_BB))
#################################################################
#
#
#               MAF and Conservative Score Part
#
#
##################################################################
				################################################
				#
				#		Phylop Score Part
				#
				#
				###############################################		
				
				conservative_score<-fread("/net/fantasia/home/borang/Susie_Mult/real_data/summary_data_lipid/070722/Annotation/conservative.txt")
				
				res_all$phylop<- conservative_score$phylop[match(res_all$SNP,conservative_score$Rsid)]
				res_all_Conservative<-res_all
				res_all_Conservative$MESuSiE_Conservative<-0
				res_all_Conservative$MESuSiE_Conservative[res_all_Conservative$MESuSiE_PIP_WB>0.5|res_all_Conservative$MESuSiE_PIP_BB>0.5]<-1
				res_all_Conservative$MESuSiE_Conservative[res_all_Conservative$MESuSiE_PIP_Shared>0.5]<-2
				res_all_Conservative$MESuSiE_Conservative<-factor(res_all_Conservative$MESuSiE_Conservative)
				levels(res_all_Conservative$MESuSiE_Conservative)<-c("None","Ancestry-specific","Shared")
				
				res_all_Conservative$SuSiE_Conservative<-0
				res_all_Conservative$SuSiE_Conservative[res_all_Conservative$SuSiE_WB>0.5|res_all_Conservative$SuSiE_BB>0.5]<-1
				res_all_Conservative$SuSiE_Conservative[res_all_Conservative$SuSiE_Shared>0.5]<-2
				res_all_Conservative$SuSiE_Conservative<-factor(res_all_Conservative$SuSiE_Conservative)
				levels(res_all_Conservative$SuSiE_Conservative)<-c("None","Ancestry-specific","Shared")
				
				mesusie_data<-data.frame(Method = "MESuSiE",Phylop = res_all_Conservative$phylop,Cat =res_all_Conservative$MESuSiE_Conservative )
				susie_data<-data.frame(Method = "SuSiE",Phylop = res_all_Conservative$phylop,Cat =res_all_Conservative$SuSiE_Conservative )
				plot_data<-rbind(mesusie_data,susie_data)
				plot_data%>%group_by(Method,Cat)%>%summarise(median_phylop = median(Phylop,na.rm=T))
				plot_data%>%group_by(Method,Cat)%>%summarise(mean_phylop = mean(Phylop,na.rm=T))
				plot_data$Signal<-factor(plot_data$Cat)
				
				JonckheereTerpstraTest(list(phylop_none$phylop,phylop_ancestry$phylop,phylop_shared$phylop),nperm=100)
				
				plot_data$Phylop_rescale<-unlist(lapply(plot_data$Phylop,function(x){
				  ifelse(x>=0,log(x+1),-log(-1*(x-1)))
				}))
				plot_data%>%group_by(Method,Cat)%>%summarise(mean(Phylop_rescale,na.rm=T),median(Phylop_rescale,na.rm=T))
				
			
			df <-plot_data %>%group_by(Method,Signal) %>%summarize(ymin = min(Phylop_rescale,na.rm = T), lower = quantile(Phylop_rescale, .25,na.rm = T),  middle = mean(Phylop_rescale,na.rm = T), upper = quantile(Phylop_rescale, .75,na.rm = T),ymax = max(Phylop_rescale,na.rm = T)) 
			
			p_set_phylop = df %>% ggplot(aes(x = factor(Method), fill=Signal)) + geom_boxplot(fatten = 6,aes(ymin = ymin, ymax = ymax, lower = lower, upper = upper, middle = middle), stat = 'identity')+scale_fill_manual(values=c("None"="#DAFFED","Ancestry-specific"="#fffbdb","Shared" = "#7776bc"))
			p_set_phylop =p_set_phylop + theme_bw() + xlab("") +ylab("Phylop Score")
			p_set_phylop= p_set_phylop+ theme(axis.text.x = element_text( size = 16),
			                                  axis.text.y = element_text( size = 16),  
			                                  axis.title.x = element_text( size = 18),
			                                  axis.title.y = element_text( size = 18,face="bold"),
			                                  strip.text.x = element_text(size = 18),
			                                  strip.text.y= element_text(size = 18),
			                                  legend.text=element_text(size=16),
			                                  legend.title=element_text(size=18,face="bold"),
			                                  plot.title = element_text(size=16,hjust = 0.5))
			p_set_phylop= p_set_phylop+ theme(legend.position="bottom")		
			
			ggsave("/net/fantasia/home/borang/Susie_Mult/real_data/summary_data_lipid/070722/Figure_New/phylop.jpeg",p_set_phylop,dpi=300,width=8,height =8)
			
			################################################
			#
			#		MAF Difference Part
			#
			#
			###############################################					
			
			
			
			maf_shared<-  data.frame(res_all)%>%filter(MESuSiE_PIP_Shared>0.5)%>%summarise(MAF_diff = MAF_WB-MAF_BB)
			maf_ancestry<-   data.frame(res_all)%>%filter(MESuSiE_PIP_WB>0.5|MESuSiE_PIP_BB>0.5)%>%summarise(MAF_diff =MAF_WB-MAF_BB)
			
		
			maf_data_plot<-data.frame("MAF" = c(maf_shared$MAF_diff,maf_ancestry$MAF_diff),"Signal"=c(rep("Shared",length(maf_shared$MAF_diff)),rep("Ancestry-specific",length(maf_ancestry$MAF_diff))))
			maf_data_plot$Signal<-factor(maf_data_plot$Signal)
			maf_data_plot$Signal<-factor(maf_data_plot$Signal,levels = c(levels(maf_data_plot$Signal),"None"))
			p_set = ggplot(data =maf_data_plot,aes(x = factor(Signal), y=MAF,fill=Signal))+geom_boxplot(aes(x = factor(Signal),fill=Signal)) +scale_fill_manual(values=c("None"="#DAFFED","Ancestry-specific"="#fffbdb","Shared" = "#7776bc"))
			p_set = p_set + stat_compare_means()
			p_set =p_set + theme_bw() + xlab("") +ylab("MAF difference")
			p_set= p_set+ theme(axis.text.x = element_text( size = 16),
			                    axis.text.y = element_text( size = 16),  
			                    axis.title.x = element_text( size = 18),
			                    axis.title.y = element_text( size = 18,face="bold"),
			                    strip.text.x = element_text(size = 18),
			                    strip.text.y= element_text(size = 18),
			                    legend.text=element_text(size=16),
			                    legend.title=element_text(size=18,face="bold"),
			                    plot.title = element_text(size=16,hjust = 0.5))
			p_set= p_set+ theme(legend.position="none")	
			p_out<-p_set+p_set_phylop+ plot_layout(guides = "collect") & theme(legend.position = "bottom")
		#	ggsave("/net/fantasia/home/borang/Susie_Mult/real_data/summary_data_lipid/070722/Figure_New/MAF_diff.jpeg",p_out,dpi=300,width=10,height =6)
			wilcox.test(maf_shared$MAF_diff, maf_ancestry$MAF_diff, alternative = "two.sided")

			
			
###########################################################
#
#
#                 Example Case HDL APOE Gene case
#
#
##########################################################
			###Function used for locuszoom plot
gwas_plot_fun<-function(res,pip,r2,lead_SNP,xlab_name,ylab_name,yintercept){
			  
			  
			  data_plot = data.frame(SNP = res$SNP,r2 = r2, POS = res$POS,pip = pip)
			  data_plot = data_plot[order(data_plot$r2),]
			  p_manhattan = ggplot() + geom_point(data =data_plot[-(data_plot$SNP==lead_SNP),], aes(x = POS, y = pip, color = r2),size = 1.5)
			  p_manhattan = p_manhattan + geom_point(data =data_plot[data_plot$SNP==lead_SNP,], aes(x = POS, y = pip),shape = 3,size = 2.5,color="red")
			  p_manhattan = p_manhattan + geom_text(data =data_plot[data_plot$SNP==lead_SNP,], mapping=aes(x=POS, y=pip, label=SNP),vjust=1.2, size=3,show.legend = FALSE)
			  p_manhattan = p_manhattan + theme_bw()+scale_color_stepsn(
			    colors = c("navy", "lightskyblue", "green", "orange", "red"),
			    breaks = seq(0.2, 0.8, by = 0.2),
			    limits = c(0, 1),
			    show.limits = TRUE,
			    na.value = 'grey50',
			    name = expression(R^2)
			  )
			  p_manhattan = p_manhattan + geom_hline(
			    yintercept =yintercept,
			    linetype = "dashed",
			    color = "grey50",
			    size = 0.5
			  ) + geom_vline(
			    xintercept = data_plot$POS[data_plot$SNP==lead_SNP],
			    linetype = "dashed",
			    color = "grey50",
			    size = 0.5
			  ) 
			  p_manhattan= p_manhattan+xlab(xlab_name)+ylab(ylab_name)
			  p_manhattan= p_manhattan+guides(fill=guide_legend(title=as.expression(bquote(R^2))))
			  return(p_manhattan)
			}
			###Function used for PIP plot	
finemap_plot_fun<-function(res,pip,r2,cat_label,lead_SNP,xlab_name,ylab_name,yintercept){
			  
			  
			  data_plot = data.frame(SNP = res$SNP,r2 = r2, POS = res$POS,pip = pip,cat =cat_label )
			  data_plot$cat<-factor(data_plot$cat)
			  levels(data_plot$cat) = list("Non" = "0","EU" = "1","AA" = "2","Shared" = "3","Paintor"="4")
			  data_plot = data_plot[order(data_plot$r2),]
			  p_manhattan = ggplot() + geom_point(data =data_plot[-(data_plot$SNP==lead_SNP),], aes(x = POS, y = pip, color = r2,shape = cat,size=cat))+scale_shape_manual(name="Category",drop=FALSE,values=c( 16, 17,15,3,18))+scale_size_manual(values=c( 1.5, 2.5,2.5,2.5,2.5))+guides( size = FALSE)
			  p_manhattan = p_manhattan + geom_point(data =data_plot[data_plot$SNP==lead_SNP,], aes(x = POS, y = pip),shape = 3,size = 2.5,color="red")
			  p_manhattan = p_manhattan + geom_text(data =data_plot[data_plot$SNP==lead_SNP,], mapping=aes(x=POS, y=pip, label=SNP),vjust=1.2, size=3,show.legend = FALSE)
			  p_manhattan = p_manhattan + theme_bw()+scale_color_stepsn(
			    colors = c("navy", "lightskyblue", "green", "orange", "red"),
			    breaks = seq(0.2, 0.8, by = 0.2),
			    limits = c(0, 1),
			    show.limits = TRUE,
			    na.value = 'grey50',
			    name = expression(R^2)
			  )
			  p_manhattan = p_manhattan + geom_hline(
			    yintercept =yintercept,
			    linetype = "dashed",
			    color = "grey50",
			    size = 0.5
			  ) + geom_vline(
			    xintercept = data_plot$POS[data_plot$SNP==lead_SNP],
			    linetype = "dashed",
			    color = "grey50",
			    size = 0.5
			  ) 
			  p_manhattan= p_manhattan+xlab(xlab_name)+ylab(ylab_name)
			  p_manhattan= p_manhattan+guides(fill=guide_legend(title=as.expression(bquote(R^2))))
			  return(p_manhattan)
}			

##############################################################
#HDL APOE case region 99
##########################################################
load("/net/fantasia/home/borang/Susie_Mult/website_material/real_data/HDL_APOE.RData")
lead_SNP<-which.max(abs(candidate_region$zscore_WB))
candidate_region$r2_EU = unname(unlist((WB_LD[,..lead_SNP])^2))
candidate_region$r2_AA = unname(unlist((BB_LD[,..lead_SNP])^2))

region_select<-candidate_region%>%filter(SuSiE_WB>0.5|SuSiE_BB>0.5|SuSiE_Shared>0.5|Paintor_PIP>0.5|MESuSiE_PIP_Shared>0.5|MESuSiE_PIP_WB>0.5|MESuSiE_PIP_BB>0.5)%>%summarise(MIN_POS = min(POS),MAX_POS = max(POS))
candidate_region<-candidate_region%>%filter(POS>=region_select$MIN_POS,POS<=region_select$MAX_POS)
lead_SNP<-which.max(abs(candidate_region$zscore_WB))
library(cowplot)
p_EU<-gwas_plot_fun(candidate_region,-log10(2*pnorm(-abs(candidate_region$zscore_WB))),candidate_region$r2_EU,candidate_region$SNP[lead_SNP],"UKBB GWAS","-log10(Pvalue)",-log10(5e-8))
p_AA<-gwas_plot_fun(candidate_region,-log10(2*pnorm(-abs(candidate_region$zscore_BB))),candidate_region$r2_AA,candidate_region$SNP[lead_SNP],"AA GWAS","-log10(Pvalue)",-log10(5e-8))


####Category Setting
SuSiE_cat = rep(0,nrow(candidate_region))
SuSiE_cat[candidate_region$SuSiE_WB>0.5&candidate_region$SuSiE_BB>0.5]<-3
SuSiE_cat[candidate_region$SuSiE_WB>0.5&candidate_region$SuSiE_BB<0.5]<-1
SuSiE_cat[candidate_region$SuSiE_WB<0.5&candidate_region$SuSiE_BB>0.5]<-2
SuSiE_cat<-factor(SuSiE_cat)

Paintor_cat = rep(0,nrow(candidate_region))
Paintor_cat[candidate_region$Paintor_PIP>0.5]<-4
Paintor_cat<-factor(Paintor_cat)

MESuSiE_cat = rep(0,nrow(candidate_region))
MESuSiE_cat[candidate_region$MESuSiE_PIP_WB>0.5]<-1
MESuSiE_cat[candidate_region$MESuSiE_PIP_BB>0.5]<-2
MESuSiE_cat[candidate_region$MESuSiE_PIP_Shared>0.5]<-3
MESuSiE_cat<-factor(MESuSiE_cat)


#################################


p_EU_SuSiE<-finemap_plot_fun(candidate_region,candidate_region$SuSiE_WB,candidate_region$r2_EU,SuSiE_cat,candidate_region$SNP[lead_SNP],"SuSiE UKBB","PIP",0.5)
p_AA_SuSiE<-finemap_plot_fun(candidate_region,candidate_region$SuSiE_BB,candidate_region$r2_AA,SuSiE_cat,candidate_region$SNP[lead_SNP],"SuSiE AA","PIP",0.5)

p_MESuSiE<-finemap_plot_fun(candidate_region,candidate_region$MESuSiE_PIP_Either,candidate_region$r2_EU,MESuSiE_cat,candidate_region$SNP[lead_SNP],"MESuSiE","PIP",0.5)
p_Paintor<-finemap_plot_fun(candidate_region,candidate_region$Paintor_PIP,candidate_region$r2_EU,Paintor_cat,candidate_region$SNP[lead_SNP],"Paintor","PIP",0.5)

library(data.table)
Gene_List<-fread("/net/fantasia/home/borang/Susie_Mult/simulation/simu_0120/data/Gencode_GRCh37_Genes_UniqueList2021.txt",header=T)
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))


plot.range <- c(min(candidate_region$POS),max(candidate_region$POS))

p2 <- ggplot(data = Gene_List_sub_coding) + 
  geom_linerange(aes(x = Gene, ymin = Start, ymax = End)) +
  coord_flip() + ylab("") + ylim(plot.range) + 
  geom_text(aes(x = Gene, y = Start, label = Gene), fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom", 
        panel.grid.major.y = element_blank()) + 
  expand_limits(y=c(-1, 1))


combined_plot<-((p_EU/p_EU_SuSiE/p_MESuSiE/p2)|(p_AA/p_AA_SuSiE/p_Paintor/p2))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
#ggsave("/net/fantasia/home/borang/Susie_Mult/real_data/summary_data_lipid/070722/Figure_New/HDL_APOE.jpeg", combined_plot ,width=12,height=8,dpi=300)


#########################################
#TC ARIC4 Gene rs6601924
#####################################

load("/net/fantasia/home/borang/Susie_Mult/website_material/real_data/TC_ARIC4.RData")
lead_SNP<-which.max(abs(candidate_region$zscore_WB))
candidate_region$r2_EU = unname(unlist((WB_LD[,..lead_SNP])^2))
candidate_region$r2_AA = unname(unlist((BB_LD[,..lead_SNP])^2))

library(cowplot)
p_EU<-gwas_plot_fun(candidate_region,-log10(2*pnorm(-abs(candidate_region$zscore_WB))),candidate_region$r2_EU,candidate_region$SNP[lead_SNP],"UKBB GWAS","-log10(Pvalue)",-log10(5e-8))
p_AA<-gwas_plot_fun(candidate_region,-log10(2*pnorm(-abs(candidate_region$zscore_BB))),candidate_region$r2_AA,candidate_region$SNP[lead_SNP],"AA GWAS","-log10(Pvalue)",-log10(5e-8))


####Category Setting
SuSiE_cat = rep(0,nrow(candidate_region))
SuSiE_cat[candidate_region$SuSiE_WB>0.5&candidate_region$SuSiE_BB>0.5]<-3
SuSiE_cat[candidate_region$SuSiE_WB>0.5&candidate_region$SuSiE_BB<0.5]<-1
SuSiE_cat[candidate_region$SuSiE_WB<0.5&candidate_region$SuSiE_BB>0.5]<-2
SuSiE_cat<-factor(SuSiE_cat)

Paintor_cat = rep(0,nrow(candidate_region))
Paintor_cat[candidate_region$Paintor_PIP>0.5]<-4
Paintor_cat<-factor(Paintor_cat)

MESuSiE_cat = rep(0,nrow(candidate_region))
MESuSiE_cat[candidate_region$MESuSiE_PIP_WB>0.5]<-1
MESuSiE_cat[candidate_region$MESuSiE_PIP_BB>0.5]<-2
MESuSiE_cat[candidate_region$MESuSiE_PIP_Shared>0.5]<-3
MESuSiE_cat<-factor(MESuSiE_cat)

##########################################
#
#		Summary for plot
#
#################################



p_EU_SuSiE<-finemap_plot_fun(candidate_region,candidate_region$SuSiE_WB,candidate_region$r2_EU,SuSiE_cat,candidate_region$SNP[lead_SNP],"SuSiE UKBB","PIP",0.5)
p_AA_SuSiE<-finemap_plot_fun(candidate_region,candidate_region$SuSiE_BB,candidate_region$r2_AA,SuSiE_cat,candidate_region$SNP[lead_SNP],"SuSiE AA","PIP",0.5)

p_MESuSiE<-finemap_plot_fun(candidate_region,candidate_region$MESuSiE_PIP_Either,candidate_region$r2_EU,MESuSiE_cat,candidate_region$SNP[lead_SNP],"MESuSiE","PIP",0.5)
p_Paintor<-finemap_plot_fun(candidate_region,candidate_region$Paintor_PIP,candidate_region$r2_EU,Paintor_cat,candidate_region$SNP[lead_SNP],"Paintor","PIP",0.5)

library(data.table)
Gene_List<-fread("/net/fantasia/home/borang/Susie_Mult/simulation/simu_0120/data/Gencode_GRCh37_Genes_UniqueList2021.txt",header=T)
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))


plot.range <- c(min(c(candidate_region$POS,Gene_List_sub_coding$Start)),max(c(candidate_region$POS,Gene_List_sub_coding$End)))

p2 <- ggplot(data = Gene_List_sub_coding) + 
  geom_linerange(aes(x = Gene, ymin = Start, ymax = End)) +
  coord_flip() + ylab("") + ylim(plot.range) + 
  geom_text(aes(x = Gene, y = Start, label = Gene), fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom", 
        panel.grid.major.y = element_blank()) + 
  expand_limits(y=c(-1, 1))


combined_plot<-((p_EU/p_EU_SuSiE/p_MESuSiE/p2)|(p_AA/p_AA_SuSiE/p_Paintor/p2))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
#ggsave(paste0("/net/fantasia/home/borang/Susie_Mult/real_data/summary_data_lipid/070722/Figure_New/TC_ARIC4.jpeg"), combined_plot ,width=12,height=8,dpi=300)




