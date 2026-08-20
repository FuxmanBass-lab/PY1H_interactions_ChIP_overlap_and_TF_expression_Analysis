metadata=read.table('metadata.tsv',header = T,row.names = 1)
Total=list()

for (i in 1:nrow(metadata)) {
  temp=Read10X(paste0('Counts/',rownames(metadata)[i]))
  temp=CreateSeuratObject(temp)
  temp$sample=rownames(metadata)[i]
  temp$donor=(metadata$donor)[i]
  temp$organ_tissue=(metadata$organ_tissue)[i]
  temp$gender=(metadata$gender)[i]
  Total[[i]]=temp
  print(i)
  rm(i,temp)
} 

Total=merge(Total[[1]],Total[-1])

Total[['percent.mt']]=PercentageFeatureSet(Total,'^MT-')
pdf('Results/QC.pdf',width = 25)
VlnPlot(Total,'nFeature_RNA',group.by = 'sample',pt.size = 0)+NoLegend()
VlnPlot(Total,'nCount_RNA',group.by = 'sample',pt.size = 0)+NoLegend()
VlnPlot(Total,'percent.mt',group.by = 'sample',pt.size = 0)+NoLegend()
dev.off()

saveRDS(Total,'Unfiltered.rds')

Total=subset(Total,subset=nFeature_RNA>500&nFeature_RNA<7500&nCount_RNA<10000&percent.mt<25)

saveRDS(Total,'Filtered.rds')



files=list.files('Results/',pattern = '.rds')


for (i in files) {
  temp=readRDS(paste0('Results/',i))
  DefaultAssay(temp)='RNA'
  temp=NormalizeData(temp)
  temp.markers=FindAllMarkers(temp,only.pos = T)
  library(dplyr)
  
  top10 = temp.markers %>% group_by(cluster) %>% top_n(wt=avg_log2FC,n=30)
  
  temp=ScaleData(temp,assay = 'RNA',features = top10$gene)
  
  pdf(paste0('Results/',gsub('.rds','',i),'.Celltype.Heatmap.pdf'),width = 10,height = 30)
  plot(DoHeatmap(temp,features = top10$gene,assay = 'RNA',label = F))
  dev.off()
  a=data.frame(table(temp$Celltype))
  colnames(a)=c('Celltype','Cells')
  pdf(paste0('Results/',gsub('.rds','',i),'.Celltype.composition.pdf'),width = 10,height = 10)
  plot(ggplot(a,aes(y=Celltype,x=Cells,fill=Celltype))+geom_bar(stat = 'identity'))
  dev.off()
}




