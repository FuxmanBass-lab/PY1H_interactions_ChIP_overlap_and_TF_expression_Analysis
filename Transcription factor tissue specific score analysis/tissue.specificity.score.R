TF.CPM=read.table('Transcription.factors.Tissue_Celltype.CPM.fixed.log2.transformed.tsv',header = T,row.names = 1,sep = '\t')

TF.specificity=data.frame(genes=rownames(TF.CPM))
TF.specificity$Tissue.Celltype.specific.score=0
TF.specificity$number.of.tissue.celltype.expressed=0

for (i in 1:nrow(TF.specificity)) {
  TF.specificity$number.of.tissue.celltype.expressed[i]=sum(as.numeric(TF.CPM[TF.specificity$genes[i],])!=0)
  average.expression=mean(as.numeric(TF.CPM[TF.specificity$genes[i],]))
  Total.expression=sum(as.numeric(TF.CPM[TF.specificity$genes[i],]))
  pi=as.numeric(TF.CPM[TF.specificity$genes[i],])/Total.expression
  specificities=c()
  for (j in 1:length(pi)) {
    
    specificity=pi[j]*log2(pi[j]/mean(pi))
    if( !is.nan(specificity)) {
      specificities=c(specificities,specificity)
    }
  }
  TF.specificity$Tissue.Celltype.specific.score[i]=sum(specificities)
  
}
write.table(TF.specificity,'Tabula.sapiens.TF.tissue.celltype.specificity.score.tsv',quote = F,sep = '\t',row.names = F)



TF.CPM=read.table('Cytokine.Tissue_Celltype.CPM.fixed.tsv',header = T,row.names = 1,sep = '\t')
TF.CPM=log(TF.CPM+1)
TF.specificity=data.frame(genes=rownames(TF.CPM))
TF.specificity$number.of.tissue.celltype.expressed=0

for (i in 1:nrow(TF.specificity)) {
  TF.specificity$number.of.tissue.celltype.expressed[i]=sum(as.numeric(TF.CPM[TF.specificity$genes[i],])!=0)
  average.expression=mean(as.numeric(TF.CPM[TF.specificity$genes[i],]))
  Total.expression=sum(as.numeric(TF.CPM[TF.specificity$genes[i],]))
  pi=as.numeric(TF.CPM[TF.specificity$genes[i],])/Total.expression
  specificities=c()
  for (j in 1:length(pi)) {
    
    specificity=pi[j]*log2(pi[j]/mean(pi))
    if( !is.nan(specificity)) {
      specificities=c(specificities,specificity)
    }
  }
  TF.specificity$Tissue.Celltype.specific.score[i]=sum(specificities)
  
}



write.table(TF.specificity,'Tabula.sapiens.Cytokine.tissue.celltype.specificity.score.tsv',quote = F,sep = '\t',row.names = F)



