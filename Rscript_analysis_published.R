#analysis for portugal D.mel populations
##written by Chaimae and edited by SKHsu
##20190801_V3

rm(list=ls())
library(edgeR)
library(topGO)
library(MESS)
library(WGCNA)
setwd("/Your/working/directory/")
####import function####
cont_table=function(query,background,classifyer){
  p1=length(Reduce(intersect,list(query,background,classifyer)))
  q1=length(intersect(query,background))-p1
  p0=length(setdiff(intersect(background,classifyer),intersect(query,classifyer)))
  q0=length(setdiff(background,query))-p0
  return(matrix(c(p1,p0,q1,q0),2,2))
}

####step1: input RNASeq count table####
counts = read.csv("count_table_published.csv", header = T, stringsAsFactors=F,row.names = 1)
counts_use=counts[apply(cpm(counts),1,function(x) !sum(x<1)>=1),]

evo=rep(c("C","H","B"),5)
#substr(colnames(counts_HF_PM),1,1)#translate labels to evolutionary states (B:ancestral; H: evolved)
y2=DGEList(counts=counts_use,group = evo)
y2=calcNormFactors(y2)

####Step2: Variance decoposition####
plotMDS(y2)
pca=prcomp(t(log(cpm(y2))))
ve=pca$sdev^2/sum(pca$sdev^2)
png("./Figure2.png",width = 8.7,height = 8.7,units = "cm",res=600,pointsize = 8)
plot(pca$x,col=c("forestgreen","royalblue","salmon")[as.factor(evo)],pch=19,asp=1,
     xlab=paste0("PC1 (",round(ve[1]*100,2),"%)"),ylab=paste0("PC2 (",round(ve[2]*100,2),"%)"))
legend("topleft",legend = c("hot-evolved","cold-evolved","ancestral"),col=c("salmon","royalblue","forestgreen"),pch=19)
dev.off()
####Step3: linear modeling and DE analysis####
ModelDesign=model.matrix(~0+evo)
DGE2=estimateDisp(y2,design = ModelDesign,robust = T)
GLM=glmFit(DGE2,design = ModelDesign)
mycontrast=makeContrasts("HB"=evoH-evoB,"CB"=evoC-evoB,"lab"=((evoH-evoB)+(evoC-evoB))/2,"interaction"=(evoH-evoB)-(evoC-evoB),
                         levels = ModelDesign)
LRT_res_HB=glmLRT(GLM,contrast = mycontrast[,"HB"])
res_table_HB=LRT_res_HB$table
res_table_HB$padj=p.adjust(res_table_HB$PValue,method = "BH")
LRT_res_CB=glmLRT(GLM,contrast = mycontrast[,"CB"])
res_table_CB=LRT_res_CB$table
res_table_CB$padj=p.adjust(res_table_CB$PValue,method = "BH")
LRT_res_LB=glmLRT(GLM,contrast = mycontrast[,"lab"])
res_table_LB=LRT_res_LB$table
res_table_LB$padj=p.adjust(res_table_LB$PValue,method = "BH")
LRT_res_interaction=glmLRT(GLM,contrast = mycontrast[,"interaction"])
res_table_interaction=LRT_res_interaction$table
res_table_interaction$padj=p.adjust(res_table_interaction$PValue,method = "BH")
write.csv(res_table_interaction,"./DE_test_TA.csv",quote = F)
write.csv(res_table_LB,"./DE_test_LB.csv",quote = F)

####step4: GO analysis####
background=rownames(y2)
query_ID5=list(lab_up=rownames(y2)[res_table_LB$padj<0.05&res_table_interaction$padj>0.05&res_table_LB$logFC>0],
               lab_down=rownames(y2)[res_table_LB$padj<0.05&res_table_interaction$padj>0.05&res_table_LB$logFC<0],
               TA_HC=rownames(y2)[res_table_LB$padj>0.05&res_table_interaction$padj<0.05&res_table_interaction$logFC>0],
               TA_CH=rownames(y2)[res_table_LB$padj>0.05&res_table_interaction$padj<0.05&res_table_interaction$logFC<0]
)

GO_res_table5=list()
for (i in 1:4){
  tmp5=factor(as.integer(rownames(y2)%in%query_ID5[[i]]))
  names(tmp5)=rownames(y2)#genelist
  tgd5=new( "topGOdata", ontology="BP", allGenes = tmp5, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")#data preparation#
  resTopGO.classic=runTest(tgd5, algorithm = "classic", statistic = "Fisher")#enrichment test#
  resTopGO.weight01=runTest(tgd5, algorithm = "weight01", statistic = "Fisher")
  tmp_res5=GenTable(tgd5,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)#analysis of results#
  GO_res_table5[[i]]=tmp_res5
}


####Step5: Comparison to clinal variation####
NA_Dmel_29=read.table("./mel-cline-gene-count.deseq2.SH-NH.exp.chr",
                      header = F,stringsAsFactors = F)
rownames(NA_Dmel_29)=NA_Dmel_29$V1
query_ID_NA_Dmel_29=list(HC=NA_Dmel_29$V1[NA_Dmel_29$V3>0&NA_Dmel_29$V7<0.05],
                         CH=NA_Dmel_29$V1[NA_Dmel_29$V3<0&NA_Dmel_29$V7<0.05])

EA_Dmel=read.csv("./hutter_2008_Dmel_22.csv",header = T,stringsAsFactors = F)
EA_Dmel=EA_Dmel[!duplicated(EA_Dmel$FBgn_ID),]
rownames(EA_Dmel)=EA_Dmel$FBgn_ID
query_ID_EA_Dmel=list(EU_Afr=read.table("./EU_Afr.txt",stringsAsFactors = F)[,1],
                      Afr_EU=read.table("./Afr_EU.txt",stringsAsFactors = F)[,1])

fisher.test(cont_table(query_ID5$TA_HC,background,query_ID_NA_Dmel_29$HC),alternative = "greater")
fisher.test(cont_table(query_ID5$TA_CH,background,query_ID_NA_Dmel_29$CH),alternative = "greater")
fisher.test(cont_table(query_ID5$TA_HC,background,query_ID_EA_Dmel$Afr_EU),alternative = "greater")
fisher.test(cont_table(query_ID5$TA_CH,background,query_ID_EA_Dmel$EU_Afr),alternative = "greater")

cor.test(res_table_interaction$logFC,NA_Dmel_29[background,]$V3,method = "spearman")
cor.test(res_table_interaction[EA_Dmel$FBgn_ID,]$logFC,EA_Dmel$log2FC,method = "spearman")
cor.test(NA_Dmel_29[EA_Dmel$FBgn_ID,3],EA_Dmel$log2FC,method = "spearman")

####Step6: Comparison to Manenti et al. 2018####
fluc_gene=read.csv("./fluctuating_genes.csv",header = T,stringsAsFactors = F)[,1]
fisher.test(cont_table(unlist(query_ID5[c(1,2)]),background,fluc_gene),alternative = "greater")

####Step7: WGCNA####
allowWGCNAThreads(24)
expr_dat=t(log(cpm(y2)))

powers = 1:20
sft_dat = pickSoftThreshold(expr_dat, powerVector = powers, verbose = 5)

plot(sft_dat$fitIndices[,1], -sign(sft_dat$fitIndices[,3])*sft_dat$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_dat$fitIndices[,1], -sign(sft_dat$fitIndices[,3])*sft_dat$fitIndices[,2],
     labels=powers,cex=1,col="red")
abline(h=0.90,col="red")

plot(sft_dat$fitIndices[,1], sft_dat$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_dat$fitIndices[,1], sft_dat$fitIndices[,5], labels=powers, cex=1,col="red")

net_dat = blockwiseModules(expr_dat, power = 6,
                           TOMType = "signed",networkType = "signed", minModuleSize = 100,
                           reassignThreshold = 1e-6, mergeCutHeight = 0.15,
                           numericLabels = TRUE, pamRespectsDendro = T,
                           saveTOMs = F,
                           verbose = 3)

table(net_dat$colors)
mergedColors = labels2colors(net_dat$colors)

plotDendroAndColors(net_dat$dendrograms[[2]], mergedColors[net_dat$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

mod_dat_GO_res=list()
for (i in sort(unique(net_dat$colors))){
  idx=net_dat$colors==i
  tmp=factor(as.integer(idx))
  names(tmp)=rownames(y2)#genelist#
  tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")#data preparation#
  resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")#enrichment test#
  resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)#analysis of results#
  mod_dat_GO_res[[i+1]]=tmp_res
}

names(mod_dat_GO_res)=paste0("Module",0:20)
lapply(mod_dat_GO_res,function(x) x[1:3,c(1:2,8)])

mod_dat_enr_res=c()
for (i in sort(unique(net_dat$colors))){
  idx=net_dat$colors==i
  tmp=sapply(query_ID5,function(x) fisher.test(cont_table(x,background,background[idx]),alternative = "greater")$p.value)
  mod_dat_enr_res=rbind(mod_dat_enr_res,tmp)
}
mod_dat_enr_res_adj=matrix(p.adjust(mod_dat_enr_res,method = "BH"),21,4,byrow = F)
rownames(mod_dat_enr_res_adj)=sort(unique(net_dat$colors))
colnames(mod_dat_enr_res_adj)=names(query_ID5)
mod_dat_enr_res_adj[apply(mod_dat_enr_res_adj,1,function(x) any(x<0.05)),]

mod_dat_enr_res_NA=c()
for (i in sort(unique(net_dat$colors))){
  idx=net_dat$colors==i
  tmp=sapply(query_ID_NA_Dmel_29,function(x) fisher.test(cont_table(x,background,background[idx]),alternative = "greater")$p.value)
  mod_dat_enr_res_NA=rbind(mod_dat_enr_res_NA,tmp)
}
mod_dat_enr_res_NA_adj=matrix(p.adjust(mod_dat_enr_res_NA,method = "BH"),21,2,byrow = F)
rownames(mod_dat_enr_res_NA_adj)=sort(unique(net_dat$colors))
colnames(mod_dat_enr_res_NA_adj)=names(query_ID_NA_Dmel_29)
mod_dat_enr_res_NA_adj[apply(mod_dat_enr_res_NA_adj,1,function(x) any(x<0.05)),]


lab_evo=c()
for (i in 1:21){
  lab_evo=c(lab_evo,paste(colnames(mod_dat_enr_res_adj)[apply(mod_dat_enr_res_adj,1,function(x) x<0.05)[,i]], collapse = "/"))
}
nat_evo=c()
for (i in 1:21){
  nat_evo=c(nat_evo,paste(colnames(mod_dat_enr_res_NA_adj)[apply(mod_dat_enr_res_NA_adj,1,function(x) x<0.05)[,i]], collapse = "/"))
}

out_tab=data.frame("name"=names(mod_dat_GO_res),"N.genes"=as.numeric(table(net_dat$colors)),
                   "top_term"=sapply(mod_dat_GO_res,function(x) x[1,2]),
                   "evolution_lab"=lab_evo,"evolution_nature"=nat_evo)

write.table(out_tab,"./WGCNA_res_table.txt",col.names = T,row.names = F,quote = F,sep="\t")

scl_expr_dat=apply(expr_dat,2,scale)
scl_avg_expr_dat=apply(scl_expr_dat,2,function(x) tapply(x,evo,mean))

png("./FigureS1.png",width = 20,height = 12,units = "cm",res=600,pointsize = 8)
par(mfrow=c(3,7))
for (i in 0:20){
  idx=net_dat$colors==i
  plot(NA,xlim=c(0.5,3.5),ylim=c(-2,2),xlab="pop",ylab="normalized mean expression",xaxt="n",
       main=paste0("Module",i),col.main=ifelse(apply(mod_dat_enr_res_adj,1,function(x) any(x<0.05))[i+1],"red","black"))
  axis(1,at=1:3,labels = c("B","C","H"))
  apply(scl_avg_expr_dat[,idx],2,function(x) lines(x,col=alpha("grey20",0.2)))
}
dev.off()


####Step8: Visualization####
png("./Figure3.png",width = 8.7,height = 8.7,units = "cm",res=600,pointsize = 6)
plot(res_table_HB$logFC,res_table_CB$logFC,asp=1,pch=19,col="grey",cex=1,cex.lab=1.4,xlim=c(-4,4),ylim=c(-4,4),
     xlab="Evolutionary response: Hot",ylab="Evolutionary response: Cold")
abline(v=0,h=0,col="grey70",lty=2)
abline(a=0,b=1,col="grey70",lty=2)
abline(a=0,b=-1,col="grey70",lty=2)
points(res_table_HB$logFC[res_table_LB$padj<0.05&res_table_interaction$padj>0.05],
       res_table_CB$logFC[res_table_LB$padj<0.05&res_table_interaction$padj>0.05],col="forestgreen",pch=19,cex=1)
points(res_table_HB$logFC[res_table_interaction$padj<0.05&res_table_LB$padj>0.05],
       res_table_CB$logFC[res_table_interaction$padj<0.05&res_table_LB$padj>0.05],col="gold",pch=19,cex=1)
legend("topleft",legend = c("concordant evolution: 541","divergent evolution: 203"), pch=19, col=c("forestgreen","gold"))
dev.off()

GOI1=intersect(genesInTerm(tgd5,whichGO = "GO:0032504")[[1]],query_ID5$lab_up)
GOI2=intersect(genesInTerm(tgd5,whichGO = "GO:0006120")[[1]],query_ID5$lab_down)
GOI3=intersect(genesInTerm(tgd5,whichGO = "GO:0034605")[[1]],query_ID5$TA_HC)
GOI4=intersect(genesInTerm(tgd5,whichGO = "GO:0019731")[[1]],query_ID5$TA_CH)

png("./Figure4.png",width = 15,height = 8.7,units = "cm",res=600,pointsize = 6)
par(bg=NA)
boxplot(res_table_HB[GOI1,1],res_table_CB[GOI1,1],
        res_table_HB[GOI2,1],res_table_CB[GOI2,1],
        res_table_HB[GOI3,1],res_table_CB[GOI3,1],
        res_table_HB[GOI4,1],res_table_CB[GOI4,1],
        border=c("salmon","royalblue"),xaxt="n",
        ylab=expression(log[2](FC)))
abline(h=0,lwd=2,lty=2)
legend("bottomleft",pch=0,col=c("salmon","royalblue"),legend=c("Hot-evolved", "Cold-evolved"),bty = "n")
mtext(c("multicellular organism reproduction\n(GO:0032504)",
        "mitochondrial electron transport,\nNADH to ubiquinone\n(GO:0006120)",
        "cellular response to heat\n(GO:0034605)",
        "antibacterial humoral response\n(GO:0019731)"),
      side = 1,line = c(2,2.5,2,2),at=c(1.5,3.5,5.5,7.5))
dev.off()



