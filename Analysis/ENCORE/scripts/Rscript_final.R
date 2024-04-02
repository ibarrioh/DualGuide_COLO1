
###This is the code for the adhoc figure generation, paths for the figure location as well as input must be modified
###All input needed for this code is coming from the 'PILOT' sub folder system or from the 'ENCORE/input', 
###both paths should we modified as needed when running locally.
###RColorBrewer together with ineq are needed to run these lines
###This code is tailored to the analysis that we are doing here and meant to reproduce the results from each of the figures
###
####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for FIGURE1
####### ####### ####### ####### ####### ####### ####### ####### 

####### FIG1G
## path_stat_plasmid -> path where the stats for the COLO1 plasmid are located, in the github rep thwy are under ANALYSIS/ENCORE/input

#load table and modify 
stat=as.matrix(read.delim('input/lib-COLO-1.stats',sep=' '))
temp=stat[,1]
stat=apply(stat[,2:4],2,as.numeric)
rownames(stat)=temp

##barplot
pdf('/pathFortheFigure/Fig1G_COLO1_plasmid_reads.pdf',width=5,height=5)

vect=c(stat['MATCH','lib.COLO.1.2']*100,
       stat['SWAP','lib.COLO.1.2']*100)
x=barplot(vect,ylim=c(0,100),ylab='GI library reads (% of reads)',
          names=c('gRNA1+gRNA2','SWAP'),las=3,border=NA)
text(x[,1],vect,round(vect,2),pos=4,srt=90)

dev.off()

###### FIG1-E
#Counts for the pilot plasmid
tab_exact=as.matrix(read.delim('PILOT/MERGED/PILOT_EXACT_COUNTS.txt'))

counts=apply(tab_exact[,9:ncol(tab_exact)],2,as.numeric)
rownames(counts)=tab_exact[,'ID']
gini_plm=counts[,'PILOT_Plasmid']

#package to calculate GINI index
library("ineq")

ineq::Gini(log(counts[,'PILOT_Plasmid']+1))

###lines for the curve

plm=cbind(counts[counts[,'PILOT_Plasmid']>0,'PILOT_Plasmid'],0,0)
plm=plm[order(plm[,1]),]
colnames(plm)=c('ori','read','sgRNA')

plm[,'read']=cumsum(plm[,'ori'])/sum(plm[,'ori'])
plm[,'sgRNA']=1:nrow(plm)
plm[,'sgRNA']=plm[,'sgRNA']/max(plm[,'sgRNA'])

pdf('/pathFortheFigure/Fig1E_plasmid1_reads.pdf',width=5,height=5)

plot(plm[,'sgRNA'],plm[,'read'],type='l',col=rgb(0.2,0.7,0.4),lwd=2,
     xlab='Fraction of total sgRNAs',ylab='Fraction of total reads')

points(c(0,1),c(0,1),type='l')
dev.off()


####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for FIGURE2
####### ####### ####### ####### ####### ####### ####### ####### 

###### FIG2-A

anot=as.matrix(read.delim('PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']

anot[anot[,'Gene2']=='','Gene2']=anot[anot[,'Gene2']=='','sgRNA2_Class']
anot[anot[,'Gene1']=='','Gene1']=anot[anot[,'Gene1']=='','sgRNA1_Class']

lfc_plm=readRDS('PILOT/MERGED/PILOC_EXACT_LFC_D3.rds')

#####Average replicas

mat=cbind(apply(lfc_plm[,c('PILOT_100X_D14_R1','PILOT_100X_D14_R2','PILOT_100X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_500X_D14_R1','PILOT_500X_D14_R2','PILOT_500X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_PCR500X_D14_R1','PILOT_PCR500X_D14_R2','PILOT_PCR500X_D14_R3')],1,mean))
colnames(mat)=c('PILOT_100X_D14','PILOT_500X_D14','PILOT_PCR500X_D14')

pal1=c(rgb(108/255,116/255,134/255),rgb(105/255,156/255,201/255),rgb(205/255,232/255,236/255))

pos=c('essential + non-targeting','non-targeting + essential','intergenic + essential','essential + intergenic')
neg=c('non-essential + non-essential','non-targeting + non-targeting','intergenic + intergenic')

pdf('pathFortheFigure/Fig2A_boxplot_500PX_D3.pdf',width=4,height=5)
boxplot(mat[anot[,'Notes']!='DistanceCut',2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%neg,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%neg,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%neg,2],
        col=alpha(pal1[c(1,2,3,2,3,2,3)],1),xlab='sgRNA Scaffold',ylab='log2FC',las=3)
dev.off()

###### FIG2-C

####Aggregate gen level
sel1=anot[,'Scaffold']=='Wildtype'
sel2=anot[,'Scaffold']=='Modified6'
sel3=anot[,'Scaffold']=='Modified7'

temp=data.frame(G1=anot[sel1,'Gene1'],G2=anot[sel1,'Gene2'],lfc=mat[sel1,'PILOT_500X_D14'])
temp <- aggregate(. ~ G1 + G2, data=temp, mean)
temp=as.matrix(temp)
wt=temp
temp=data.frame(G1=anot[sel2,'Gene1'],G2=anot[sel2,'Gene2'],lfc=mat[sel2,'PILOT_500X_D14'])
temp <- aggregate(. ~ G1 + G2 , data=temp, mean)
temp=as.matrix(temp)
m6=temp

temp=data.frame(G1=anot[sel3,'Gene1'],G2=anot[sel3,'Gene2'],lfc=mat[sel3,'PILOT_500X_D14'])
temp <- aggregate(. ~ G1 + G2 , data=temp, mean)
temp=as.matrix(temp)
m7=temp

rownames(wt)=paste(wt[,1],wt[,2],sep='_')
temp=wt[,'lfc']
names(temp)=paste(wt[,2],wt[,1],sep='_')
wt=wt[rownames(wt)%in%names(temp),]
wt=cbind(wt,temp[rownames(wt)])

rownames(m6)=paste(m6[,1],m6[,2],sep='_')
temp=m6[,'lfc']
names(temp)=paste(m6[,2],m6[,1],sep='_')
m6=m6[rownames(m6)%in%names(temp),]
m6=cbind(m6,temp[rownames(m6)])

rownames(m7)=paste(m7[,1],m7[,2],sep='_')
temp=m7[,'lfc']
names(temp)=paste(m7[,2],m7[,1],sep='_')
m7=m7[rownames(m7)%in%names(temp),]
m7=cbind(m7,temp[rownames(m7)])

wt=cbind(as.numeric(wt[,3]),as.numeric(wt[,4]))
colnames(wt)=c('sgRNA1','sgRNA2')
m6=cbind(as.numeric(m6[,3]),as.numeric(m6[,4]))
colnames(m6)=c('sgRNA1','sgRNA2')
m7=cbind(as.numeric(m7[,3]),as.numeric(m7[,4]))
colnames(m7)=c('sgRNA1','sgRNA2')

###WT
data <- data.frame(sgRNA1 = wt[,1], sgRNA2 = wt[,2])
# Fit linear model
reg<-lm(sgRNA1 ~ sgRNA2, data = data)
# Scatter plot

pdf('pathFortheFigure/Fig2C_WT_500PX_D3.pdf',width=15,height=5)
par(mfrow=c(1,3))

data <- data.frame(sgRNA1 = wt[,1], sgRNA2 = wt[,2])
reg<-lm(sgRNA1 ~ sgRNA2, data = data)

plot(data,xlim=c(-4,2),ylim=c(-4,2),
     main='Wildtype (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m6[,1], sgRNA2 = m6[,2])
reg<-lm(sgRNA1 ~ sgRNA2, data = data)

plot(data,xlim=c(-4,2),ylim=c(-4,2),
     main='M6 (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m7[,1], sgRNA2 = m7[,2])
reg<-lm(sgRNA1 ~ sgRNA2, data = data)

plot(data,xlim=c(-4,2),ylim=c(-4,2),
     main='M& (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))
dev.off()

###### FIG2-D

######Selection of guides and colors

selection=c('essential + non-targeting','non-targeting + essential','intergenic + essential','essential + intergenic',
            'non-essential + non-essential','non-targeting + non-targeting','intergenic + intergenic')
colores=c("#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#3182bd", "#393b79","#756bb1")

###vectors to store values

AUCS_wt=rep(0,length(selection))
names(AUCS_wt)=selection
AUCS_m6=rep(0,length(selection))
names(AUCS_m6)=selection
AUCS_m7=rep(0,length(selection))
names(AUCS_m7)=selection
AUCS_ze=rep(0,length(selection))
names(AUCS_ze)=selection

n_wt=AUCS_wt
n_m6=AUCS_m6
n_m7=AUCS_m7

AUCS_max=rep(0,3)
names(AUCS_max)=c('wt_500X','m6_500X','m7_500X')
n_max=AUCS_max

####RECALL CURVES

### WT
pdf('pathFortheFigure/Fig2D_recall_500PX_D3.pdf',width=12,height=4)
par(mfrow=c(1,3))

stype_df=sort(mat[scaf1,'PILOT_500X_D14'])

plot(0,type='n',ylim=c(0,1),xlim=c(0,1),xlab=c("Percent-rank of vectors"),ylab=c('Cumulative fraction'),
     main='PILOT_500X_D14_WT')

for (i in 1:length(selection)){
  
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  y=cumsum(binary)/sum(binary)
  x=(1:length(stype_df))/length(stype_df)
  
  AUCS_wt[i]=integrate(approxfun(x,y),lower=min(x),upper=1,subdivisions = 2000)$value
  n_wt[i]=sum(binary==1)
  
  points(x,y,type='l',col=colores[i]) 
  
}
binary=rep(0,length(stype_df))
binary[1:sum(names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID'])]=1
y=cumsum(binary)/sum(binary)
x=(1:length(stype_df))/length(stype_df)
points(x,y,type='l',col=rgb(0.5,0.5,0.5),lty=2) 

AUCS_max['wt_500X']=integrate(approxfun(x,y),lower=min(x),upper=1,subdivisions = 2000)$value
n_max['wt_500X']=sum(binary==1)
###M6
stype_df=sort(mat[scaf2,'PILOT_500X_D14'])

plot(0,type='n',ylim=c(0,1),xlim=c(0,1),xlab=c("Percent-rank of vectors"),ylab=c('Cumulative fraction'),
     main='PILOT_500X_D14_M6')

for (i in 1:length(selection)){
  
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  y=cumsum(binary)/sum(binary)
  x=(1:length(stype_df))/length(stype_df)
  
  AUCS_m6[i]=integrate(approxfun(x,y),lower=min(x),upper=1,subdivisions = 2000)$value
  n_m6[i]=sum(binary==1)
  
  points(x,y,type='l',col=colores[i]) 
  
}
binary=rep(0,length(stype_df))
binary[1:sum(names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID'])]=1
y=cumsum(binary)/sum(binary)
x=(1:length(stype_df))/length(stype_df)
points(x,y,type='l',col=rgb(0.5,0.5,0.5),lty=2) 

AUCS_max['m6_500X']=integrate(approxfun(x,y),lower=min(x),upper=1,subdivisions = 2000)$value
n_max['m6_500X']=sum(binary==1)
##M7
stype_df=sort(mat[scaf3,'PILOT_500X_D14'])

plot(0,type='n',ylim=c(0,1),xlim=c(0,1),xlab=c("Percent-rank of vectors"),ylab=c('Cumulative fraction'),
     main='PILOT_500X_D14_M7')

for (i in 1:length(selection)){
  
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  y=cumsum(binary)/sum(binary)
  x=(1:length(stype_df))/length(stype_df)
  
  AUCS_m7[i]=integrate(approxfun(x,y),lower=min(x),upper=1,subdivisions = 2000)$value
  n_m7[i]=sum(binary==1)
  
  points(x,y,type='l',col=colores[i]) 
  
}
binary=rep(0,length(stype_df))
binary[1:sum(names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID'])]=1
y=cumsum(binary)/sum(binary)
x=(1:length(stype_df))/length(stype_df)
points(x,y,type='l',col=rgb(0.5,0.5,0.5),lty=2) 

AUCS_max['m7_500X']=integrate(approxfun(x,y),lower=min(x),upper=1,subdivisions = 2000)$value
n_max['m7_500X']=sum(binary==1)
dev.off()

####### FIG2B

##########SCORE comparaisson

score_ag=readRDS('/input/SCORE_median_lfc.rds')

lfc=readRDS('/PILOT/MERGED/PILOC_EXACT_LFC_D3.rds')
anot=as.matrix(read.delim('PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']

mat=cbind(apply(lfc[,c('PILOT_100X_D14_R1','PILOT_100X_D14_R2','PILOT_100X_D14_R3')],1,mean),
          apply(lfc[,c('PILOT_500X_D14_R1','PILOT_500X_D14_R2','PILOT_500X_D14_R3')],1,mean),
          apply(lfc[,c('PILOT_PCR500X_D14_R1','PILOT_PCR500X_D14_R2','PILOT_PCR500X_D14_R3')],1,mean))
colnames(mat)=c('PILOT_100X_D14','PILOT_500X_D14','PILOT_PCR500X_D14')

###Preparamos los genes

anot=cbind(anot,'')
colnames(anot)[ncol(anot)]='GEN'

sel=(anot[,'sgRNA1_Class']=='essential' & anot[,'sgRNA2_Class']%in%c('intergenic','non-targeting')) |
  (anot[,'sgRNA1_Class']=='non-essential' & anot[,'sgRNA2_Class']%in%c('intergenic','non-targeting')) |
  (anot[,'sgRNA1_Class']=='unclassified' & anot[,'sgRNA2_Class']%in%c('intergenic','non-targeting'))

anot[sel,'GEN']=anot[sel,'Gene1']

sel=(anot[,'sgRNA2_Class']=='essential' & anot[,'sgRNA1_Class']%in%c('intergenic','non-targeting')) |
  (anot[,'sgRNA2_Class']=='non-essential' & anot[,'sgRNA1_Class']%in%c('intergenic','non-targeting')) |
  (anot[,'sgRNA2_Class']=='unclassified' & anot[,'sgRNA1_Class']%in%c('intergenic','non-targeting'))

anot[sel,'GEN']=anot[sel,'Gene2']

#######
sel=anot[,'Scaffold']=='Wildtype' & anot[,'Scaffold']!='DistanceCut' & anot[,'GEN']!=''

scaff_wt=data.frame(gen=anot[sel,'GEN'],
                    P_100X_D14=mat[sel,1],
                    P_500X_D14=mat[sel,2],
                    P_PCR500X_D14=mat[sel,3])
scaff_wt=as.matrix(aggregate(. ~ gen , data=scaff_wt, mean))

sel=anot[,'Scaffold']=='Modified6' & anot[,'Scaffold']!='DistanceCut' & anot[,'GEN']!=''
scaff_m6=data.frame(gen=anot[sel,'GEN'],
                    P_100X_D14=mat[sel,1],
                    P_500X_D14=mat[sel,2],
                    P_PCR500X_D14=mat[sel,3])
scaff_m6=as.matrix(aggregate(. ~ gen , data=scaff_m6, mean))

sel=anot[,'Scaffold']=='Modified7' & anot[,'Scaffold']!='DistanceCut' & anot[,'GEN']!=''
scaff_m7=data.frame(gen=anot[sel,'GEN'],
                    P_100X_D14=mat[sel,1],
                    P_500X_D14=mat[sel,2],
                    P_PCR500X_D14=mat[sel,3])
scaff_m7=as.matrix(aggregate(. ~ gen , data=scaff_m7, mean))

temp=scaff_wt[,1]
scaff_wt=apply(scaff_wt[,2:4],2,as.numeric)
rownames(scaff_wt)=temp

temp=scaff_m6[,1]
scaff_m6=apply(scaff_m6[,2:4],2,as.numeric)
rownames(scaff_m6)=temp

temp=scaff_m7[,1]
scaff_m7=apply(scaff_m7[,2:4],2,as.numeric)
rownames(scaff_m7)=temp

scaff_wt=cbind(scaff_wt,score_ag[rownames(scaff_wt)])
scaff_m6=cbind(scaff_m6,score_ag[rownames(scaff_m6)])
scaff_m7=cbind(scaff_m7,score_ag[rownames(scaff_m7)])

scaff_wt=na.omit(scaff_wt)
scaff_m6=na.omit(scaff_m6)
scaff_m7=na.omit(scaff_m7)

##500X Wildtype scaffold
data <- data.frame(DUAL = scaff_wt[,'P_500X_D14'], SCO = scaff_wt[,4])
# Fit linear model
reg<-lm(SCO ~ DUAL, data = data)

# Scatter plot

pdf('pathFortheFigure/Fig2B_WT_500PX_D3.pdf',width=5,height=5)

plot(data,xlim=c(-6,1),ylim=c(-6,1),
     main='Gene level correlation single KO',xlab='DualKO (log2FC)',ylab='SingleKO-SCORE (log2FC)',
     col=rgb(0.2,0.5,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
text(-6,-1,paste('cor:',round(cor(data[,1],data[,2]),2)),pos=4)
text(-6,-1.5,paste('reg:',round(coef(reg)[2],2)),pos=4)
text(-6,-2,paste('n:',nrow(data)),pos=4)

dev.off()

####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for FIGURE3
####### ####### ####### ####### ####### ####### ####### ####### 
############## Fig3A
##numeric values
colo1=readRDS('/input/colo1_plm.rds')
##character values
IDs=readRDS('/input/IDs_colo1.rds')
colo1=colo1[,grep('136',colnames(colo1))]
colo1=apply(colo1[,1:3],1,mean)
##Selecting conditions
sel1=c('AnchorCombinations','LibraryCombinations','GIControlsCombinations')
sel2=c('AnchorSingletons','LibrarySingletons','GIControlsSingletons')
sel3=c('ESSENTIAL-INTERGENIC','ESSENTIAL-NONTARGET',
       'INTERGENIC-INTERGENIC','NONESSENTIAL-NONESSENTIAL','NONTARGET-NONTARGET')

pdf('pathFortheFigure/Fig3A_boxplot_ENCORE_mean_v2.pdf',width=4,height=5)
boxplot(colo1[IDs[,'MyNote']%in%sel1],colo1[IDs[,'MyNote']%in%sel2],colo1[IDs[,'MyNote']==sel3[1]],
        colo1[IDs[,'MyNote']==sel3[2]],colo1[IDs[,'MyNote']==sel3[3]],colo1[IDs[,'MyNote']==sel3[4]],
        colo1[IDs[,'MyNote']==sel3[5]],names=c('DoubleKO','SingleKO',sel3),las=3,ylab='ENCORE-HT19 log(FC)',col='white')
dev.off()

############# Fig3B 

IDs_c1=readRDS('/input/IDs_colo1.rds')
c1=as.matrix(readRDS('/input/colo1_plm_GEN.rds'))

row_mean_na_rm <- function(x) {
  mean(x, na.rm = TRUE)
}

cell=unique(unlist(strsplit(colnames(c1),'_'))[c(T,F)])
temp=matrix(0,nrow(c1),length(cell))
rownames(temp)=rownames(c1)
colnames(temp)=cell

for (i in 1:ncol(temp)){
  temp[,i]=apply(c1[,grep(colnames(temp)[i],colnames(c1))],1,row_mean_na_rm)
}
c1=temp

####Genetic interactions 
map_mat=c1[grep('GIControls',rownames(c1)),]
map_ids=IDs_c1[IDs_c1[,'Note']=='GIControls',]
map_ids=cbind(map_ids,paste(map_ids[,'Note'],map_ids[,'Gene1'],map_ids[,'Gene2'],sep='_'))
colnames(map_ids)[ncol(map_ids)]='pair'
comb=map_mat[rownames(map_mat)%in%map_ids[map_ids[,'MyNote']=='GIControlsCombinations','pair'],]
sing=map_mat[rownames(map_mat)%in%map_ids[map_ids[,'MyNote']=='GIControlsSingletons','pair'],]
ctrl=c('ADAD1','CYLC2','KLK12')
sum(map_ids[map_ids[,'MyNote']=='GIControlsSingletons','Gene1']%in%ctrl)
anch=sing[rownames(sing)%in%map_ids[map_ids[,'MyNote']=='GIControlsSingletons' & map_ids[,'Gene1']%in%ctrl,'pair'],]
libr=sing[rownames(sing)%in%map_ids[map_ids[,'MyNote']=='GIControlsSingletons' & map_ids[,'Gene2']%in%ctrl,'pair'],]

combi_c1=rbind(c1[c(grep('LibraryCombinations',rownames(c1)),grep('AnchorCombinations',rownames(c1))),],comb)
a_c1=rbind(c1[grep('AnchorSingletons',rownames(c1)),],anch)
l_c1=rbind(c1[grep('LibrarySingletons',rownames(c1)),],libr)

IDs_combi_c1=cbind(unlist(strsplit(rownames(combi_c1),'_'))[c(F,T,F)],
                   unlist(strsplit(rownames(combi_c1),'_'))[c(F,F,T)])

colnames(IDs_combi_c1)=c('lib','anc')

temp=data.frame(a_c1,gen=unlist(strsplit(rownames(a_c1),'_'))[c(F,F,T)])
temp=aggregate(. ~ gen , data=temp, median, na.action = NULL, na.rm=T)
a_c1=apply(as.matrix(temp[,2:ncol(temp)]),2,as.numeric)
rownames(a_c1)=temp$gen

temp=data.frame(l_c1,gen=unlist(strsplit(rownames(l_c1),'_'))[c(F,T,F)])
temp=aggregate(. ~ gen , data=temp, median, na.action = NULL, na.rm=T)
l_c1=apply(as.matrix(temp[,2:ncol(temp)]),2,as.numeric)
rownames(l_c1)=temp$gen

#######PARALOGS from Ryan et al https://doi.org/10.1016/j.cels.2021.08.006

rownames(IDs_combi_c1)=paste(IDs_combi_c1[,1],IDs_combi_c1[,2],sep='_')

m_c1=c(a_c1[,grep('136',colnames(a_c1))],l_c1[,grep('136',colnames(a_c1))])

combi_c1=combi_c1[,grep('136',colnames(combi_c1))]

combi_c1=cbind(combi_c1,
               m_c1[IDs_combi_c1[,'lib']],
               m_c1[IDs_combi_c1[,'anc']])

rownames(combi_c1)=paste(unlist(strsplit(rownames(combi_c1),'_'))[c(F,T,F)],
                         unlist(strsplit(rownames(combi_c1),'_'))[c(F,F,T)],sep='_')
colnames(combi_c1)=c('pair','lib','anc')

temp=rep(0,nrow(combi_c1))

for (i in 1:nrow(combi_c1)){
  
  temp[i]=sum(is.na(combi_c1[i,]))
  
}

combi_c1=combi_c1[temp==0,]

combi_c1=cbind(combi_c1,combi_c1[,'pair']-(combi_c1[,'anc']+combi_c1[,'lib']))
colnames(combi_c1)[ncol(combi_c1)]='sco'

combi_c1=combi_c1[order(combi_c1[,'sco']),]
combi_c1[combi_c1[,'pair']<=-1 & combi_c1[,'sco']<=-1 ,]

###uploading table
para=as.matrix(read.delim('/inpout/paralog_scores_ryan.csv',sep=','))
rownames(para)=para[,'sorted_gene_pair']
para=para[,c('prediction_rank','prediction_percentile','prediction_score','sorted_gene_pair','A1','A2' )]

#####BOXPLOT

sel1=as.numeric(para[,'prediction_percentile'])<=10 & as.numeric(para[,'prediction_score'])<0.2
sel2=as.numeric(para[,'prediction_score'])>=0.2

temp1=c(paste(para[sel1,'A1'],para[sel1,'A2'],sep='_'),
        paste(para[sel1,'A2'],para[sel1,'A1'],sep='_'))
temp2=c(paste(para[sel2,'A1'],para[sel2,'A2'],sep='_'),
        paste(para[sel2,'A2'],para[sel2,'A1'],sep='_'))


sort_bliss=cbind(c(1:nrow(combi_c1)),
                 combi_c1[,'sco'])

x1=sort_bliss[!rownames(combi_c1)%in%c(temp1,temp2),2]
x2=sort_bliss[rownames(combi_c1)%in%temp1,2]
x3=sort_bliss[rownames(combi_c1)%in%temp2,2]

pdf('/pathFortheFigure/fig3B_bliss.pdf',width=3,height = 5)

box=boxplot(x1,c(x2,x3),ylim=c(-2,1.5),col='white',
            ylab='bliss score (observed-expected)',main='Bliss score based on logFC (HT29)')
points(jitter(rep(2,length(x3)),amount=0.2),x3,pch=16,col=rgb(0.7,0.2,0.2))
text(1.5,1,paste('p.val:',round(t.test(x1,c(x2,x3),exact=T)$p.value,4)))

dev.off()


#############

###########Figure 3C 

colo1=readRDS('/input/colo1_plm_scale.rds')
IDs=readRDS('/input/IDs_colo1.rds')
colo1=colo1[,grep('136',colnames(colo1))]
colo1=apply(colo1[,1:3],1,mean)


table(IDs[IDs[,'Gene1']=='MAPK1' & 
            IDs[,'Gene2']=='MAPK3'  ,'MyNote'])

table(IDs[IDs[,'Gene1']=='CNOT7','MyNote'])
table(IDs[IDs[,'Gene2']=='MAPK3','MyNote'])


IDs[IDs[,'Gene1']=='HDAC1' & IDs[,'Gene2']!='HDAC2' & IDs[,'Note']=='GIControls',]

h1=colo1[IDs[IDs[,'Gene1']=='HDAC1' & IDs[,'Note']=='LibrarySingletons','ID']]
h2=colo1[IDs[IDs[,'Gene1']=='HDAC1' & IDs[,'Gene2']=='HDAC2','ID']]
h3=colo1[IDs[IDs[,'Gene2']=='HDAC2' & IDs[,'MyNote']=='GIControlsSingletons','ID']]
a1=colo1[IDs[IDs[,'Gene1']=='ASF1A' & IDs[,'Note']=='LibrarySingletons','ID']]
a2=colo1[IDs[IDs[,'Gene1']=='ASF1A' & IDs[,'Gene2']=='ASF1B','ID']]
a3=colo1[IDs[IDs[,'Gene2']=='ASF1B' & IDs[,'MyNote']=='GIControlsSingletons','ID']]
c1=colo1[IDs[IDs[,'Gene1']=='CNOT7' & IDs[,'MyNote']=='GIControlsSingletons','ID']]
c2=colo1[IDs[IDs[,'Gene1']=='CNOT7' & IDs[,'Gene2']=='CNOT8','ID']]
c3=colo1[IDs[IDs[,'Gene2']=='CNOT8' & IDs[,'MyNote']=='GIControlsSingletons','ID']]
m1=colo1[IDs[IDs[,'Gene1']=='MAPK1' & IDs[,'MyNote']=='GIControlsSingletons','ID']]
m2=colo1[IDs[IDs[,'Gene1']=='MAPK1' & IDs[,'Gene2']=='MAPK3','ID']]
m3=colo1[IDs[IDs[,'Gene2']=='MAPK3' & IDs[,'MyNote']=='GIControlsSingletons','ID']]

pal=c(rgb(170/255,222/255,135/255),
      rgb(255/255,230/255,128/255),
      rgb(170/255,222/255,135/255))

pdf('pathFortheFigure/Fig3C_boxplot_examples.pdf',width=5,height=6)

boxplot(h1,h2,h3,a1,a2,a3,c1,c2,c3,m1,m2,m3,
        names=c('HDAC1','HDAC1+HDAC2','HDAC2',
                'ASF1A','ASF1A+ASF1B','ASF1B',
                'CNOT7','CNOT7+CNOT8','CNOT8',
                'MAPK1','MAPK1+MAPK3','MAPK3'),col='white',ylim=c(-3,2),outline=F,las=3,
        ylab='logFC(HT29-sgRNA) known paralogs')

points(jitter(rep(1,length(h1)),amount=0.3),h1,col=pal[1],pch=16)
points(jitter(rep(2,length(h2)),amount=0.3),h2,col=pal[2],pch=16)
points(jitter(rep(3,length(h3)),amount=0.3),h3,col=pal[3],pch=16)
points(jitter(rep(4,length(a1)),amount=0.3),a1,col=pal[1],pch=16)
points(jitter(rep(5,length(a2)),amount=0.3),a2,col=pal[2],pch=16)
points(jitter(rep(6,length(a3)),amount=0.3),a3,col=pal[3],pch=16)
points(jitter(rep(7,length(c1)),amount=0.3),c1,col=pal[1],pch=16)
points(jitter(rep(8,length(c2)),amount=0.3),c2,col=pal[2],pch=16)
points(jitter(rep(9,length(c3)),amount=0.3),c3,col=pal[3],pch=16)
points(jitter(rep(10,length(m1)),amount=0.3),m1,col=pal[1],pch=16)
points(jitter(rep(11,length(m2)),amount=0.3),m2,col=pal[2],pch=16)
points(jitter(rep(12,length(m3)),amount=0.3),m3,col=pal[3],pch=16)

t_h=c(round(t.test(h1,h2)$p.value,4),round(t.test(h3,h2)$p.value,4))
t_a=c(round(t.test(a1,a2)$p.value,4),round(t.test(a3,a2)$p.value,4))
t_c=c(round(t.test(c1,c2)$p.value,4),round(t.test(c3,c2)$p.value,4))
t_m=c(round(t.test(m1,m2)$p.value,4),round(t.test(m3,m2)$p.value,4))

text(c(1,3),c(1,1),paste('pval:',t_h),pos=4,srt=90)
text(c(4,6),c(1,1),paste('pval:',t_a),pos=4,srt=90)
text(c(7,9),c(1,1),paste('pval:',t_c),pos=4,srt=90)
text(c(10,12),c(1,1),paste('pval:',t_m),pos=4,srt=90)

dev.off()



####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for SUPP FIGURE1
####### ####### ####### ####### ####### ####### ####### ####### 

##### Foir pannel B, plasmid mappings

###Columns to keep plus naming
pos=c(0,10,11,12,13,14,16,17,19,20,22,23,25,26,29,30,31,32,33)+1
nam=c('Readname','VectorID1','V1edit1R1','V1direction1R1','V1edit2R1','V1direction2R1','Scaffold','Scaffold_edit','Linker','Linker_edit','tRNA','tRNA_edit','Improved','Improved_edit','VectorID2','V2edit1R1','V2direction1R1','V2edit2R2','V2direction2R2')
##REsults table
tab_final=matrix(0,4,3)

colnames(tab_final)=c('wt','m6','m7')
rownames(tab_final)=c('sgrna1_exists','sgrna2_exists',
                      'sgrna_match','swap')
###plasmid lib
x=as.matrix(read.table("PILOT/MAPPING/WALK_280/lib_S5_walk280.mapping.out", sep="\t", skip=3, header=FALSE))
x=x[,pos]
colnames(x)=nam

for (j in 1:ncol(x)){
  x[,j]=gsub(' ','',x[,j])
}

####Los NA nos estan jodiendo, vamso a intentar otra estratagema

for (j in 1:ncol(x)){
  x[is.na( x[,j]),j]='No_Hay'
}

##scaffolds
wt=x[x[,'Scaffold']=='Wildtype',]
m6=x[x[,'Scaffold']=='Modified6',]
m7=x[x[,'Scaffold']=='Modified7',]

#WT

a=wt

##Misma estrategia, a para limpia de NAs, selectiones y b para counts
total_number_of_reads=nrow(a)

###########

sel=(a[,'VectorID1']!="OEM2" & a[,'VectorID1']!="OEM2_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM2" & a[,'VectorID2']!="OEM2_WO" & a[,'VectorID2']!="NOMATCH")

b = sum(as.numeric(a[sel,'V1edit1R1'])==0 & 
          as.numeric(a[sel,'V2edit1R1'])==0 & 
          a[sel,'V1direction1R1']=="FALSE" & 
          a[sel,'V2direction1R1']=="FALSE")

sgRNA1_exists=(b/total_number_of_reads)*100

###########

sel=(a[,'VectorID1']!="OEM1" & a[,'VectorID1']!="OEM1_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM1" & a[,'VectorID2']!="OEM1_WO" & a[,'VectorID2']!="NOMATCH")

b= sum(a[sel,'V1edit2R1']==0 & 
         a[sel,'V2edit2R2']==0 & 
         a[sel,'V1direction2R1']=="FALSE" & 
         a[sel,'V2direction2R2']=="R")

sgRNA2_exists=(b/total_number_of_reads)*100

#################

b=sum(  a[,'VectorID1']==a[,'VectorID2'] & 
          a[,'V1edit2R1']==0 & 
          a[,'V2edit1R1']==0 &
          a[,'V1direction1R1']=="FALSE" & 
          a[,'V1direction2R1']=="FALSE" &
          a[,'V2edit1R1']==0 &
          a[,'V2edit2R2']==0 &
          a[,'V2direction1R1']=="FALSE" & 
          a[,'V2direction2R2']=="R")

sgRNA_match = (b/total_number_of_reads)*100

#################

b =sum(a[,'VectorID1']==a[,'VectorID2'] & 
         a[,'VectorID1']=="SWAP" & 
         a[,'V1edit2R1']==0 & 
         a[,'V2edit1R1']==0 &
         a[,'V2edit1R1']==0 & 
         a[,'V2edit2R2']==0)
swap = (b/total_number_of_reads)*100

########
tab_final[,'wt']=c(sgRNA1_exists,sgRNA2_exists,sgRNA_match,swap)

#m6

a=m6

##Misma estrategia, a para limpia de NAs, selectiones y b para counts
total_number_of_reads=nrow(a)

###########

sel=(a[,'VectorID1']!="OEM2" & a[,'VectorID1']!="OEM2_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM2" & a[,'VectorID2']!="OEM2_WO" & a[,'VectorID2']!="NOMATCH")

b = sum(as.numeric(a[sel,'V1edit1R1'])==0 & 
          as.numeric(a[sel,'V2edit1R1'])==0 & 
          a[sel,'V1direction1R1']=="FALSE" & 
          a[sel,'V2direction1R1']=="FALSE")

sgRNA1_exists=(b/total_number_of_reads)*100

###########

sel=(a[,'VectorID1']!="OEM1" & a[,'VectorID1']!="OEM1_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM1" & a[,'VectorID2']!="OEM1_WO" & a[,'VectorID2']!="NOMATCH")

b= sum(a[sel,'V1edit2R1']==0 & 
         a[sel,'V2edit2R2']==0 & 
         a[sel,'V1direction2R1']=="FALSE" & 
         a[sel,'V2direction2R2']=="R")

sgRNA2_exists=(b/total_number_of_reads)*100

#################

b=sum(  a[,'VectorID1']==a[,'VectorID2'] & 
          a[,'V1edit2R1']==0 & 
          a[,'V2edit1R1']==0 &
          a[,'V1direction1R1']=="FALSE" & 
          a[,'V1direction2R1']=="FALSE" &
          a[,'V2edit1R1']==0 &
          a[,'V2edit2R2']==0 &
          a[,'V2direction1R1']=="FALSE" & 
          a[,'V2direction2R2']=="R")

sgRNA_match = (b/total_number_of_reads)*100

#################

b =sum(a[,'VectorID1']==a[,'VectorID2'] & 
         a[,'VectorID1']=="SWAP" & 
         a[,'V1edit2R1']==0 & 
         a[,'V2edit1R1']==0 &
         a[,'V2edit1R1']==0 & 
         a[,'V2edit2R2']==0)
swap = (b/total_number_of_reads)*100

########
tab_final[,'m6']=c(sgRNA1_exists,sgRNA2_exists,sgRNA_match,swap)

#m7

a=m7

##Misma estrategia, a para limpia de NAs, selectiones y b para counts
total_number_of_reads=nrow(a)

###########

sel=(a[,'VectorID1']!="OEM2" & a[,'VectorID1']!="OEM2_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM2" & a[,'VectorID2']!="OEM2_WO" & a[,'VectorID2']!="NOMATCH")

b = sum(as.numeric(a[sel,'V1edit1R1'])==0 & 
          as.numeric(a[sel,'V2edit1R1'])==0 & 
          a[sel,'V1direction1R1']=="FALSE" & 
          a[sel,'V2direction1R1']=="FALSE")

sgRNA1_exists=(b/total_number_of_reads)*100

###########

sel=(a[,'VectorID1']!="OEM1" & a[,'VectorID1']!="OEM1_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM1" & a[,'VectorID2']!="OEM1_WO" & a[,'VectorID2']!="NOMATCH")

b= sum(a[sel,'V1edit2R1']==0 & 
         a[sel,'V2edit2R2']==0 & 
         a[sel,'V1direction2R1']=="FALSE" & 
         a[sel,'V2direction2R2']=="R")

sgRNA2_exists=(b/total_number_of_reads)*100

#################

b=sum(  a[,'VectorID1']==a[,'VectorID2'] & 
          a[,'V1edit2R1']==0 & 
          a[,'V2edit1R1']==0 &
          a[,'V1direction1R1']=="FALSE" & 
          a[,'V1direction2R1']=="FALSE" &
          a[,'V2edit1R1']==0 &
          a[,'V2edit2R2']==0 &
          a[,'V2direction1R1']=="FALSE" & 
          a[,'V2direction2R2']=="R")

sgRNA_match = (b/total_number_of_reads)*100

#################

b =sum(a[,'VectorID1']==a[,'VectorID2'] & 
         a[,'VectorID1']=="SWAP" & 
         a[,'V1edit2R1']==0 & 
         a[,'V2edit1R1']==0 &
         a[,'V2edit1R1']==0 & 
         a[,'V2edit2R2']==0)
swap = (b/total_number_of_reads)*100

########
tab_final[,'m7']=c(sgRNA1_exists,sgRNA2_exists,sgRNA_match,swap)

pdf('pathFortheFigure/Fig1suppB_overall_stats.pdf',width=7,height=5)
x=barplot(t(tab_final),beside=T,ylim=c(0,100),legend=c('WT','M6','M7'),col=c(rgb(233/255,164/255,49/255),
                                                                             rgb(205/255,228/255,241/255),
                                                                             rgb(67/255,82/255,114/255)),border=NA,
          main='Overall stats (Pilot library plasmid)',ylab='percentage')

text(x[,1],tab_final[1,]+2,round(tab_final[1,],1))
text(x[,2],tab_final[2,]+2,round(tab_final[2,],1))
text(x[,3],tab_final[3,]+2,round(tab_final[3,],1))
text(x[,4],tab_final[4,]+2,round(tab_final[4,],1))
dev.off()

####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for SUPP FIGURE2
####### ####### ####### ####### ####### ####### ####### ####### 

######### Supp figure 2A


anot=as.matrix(read.delim('/PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']
lfc_plm=readRDS('/PILOT/MERGED/PILOC_EXACT_LFC_D3.rds')

scaf1=anot[,'Scaffold']=='Wildtype'& anot[,'Notes']!='DistanceCut'
scaf2=anot[,'Scaffold']=='Modified6'& anot[,'Notes']!='DistanceCut'
scaf3=anot[,'Scaffold']=='Modified7'& anot[,'Notes']!='DistanceCut'

##We cut scaf WT to be same pairs as the others 
guides_ID=paste(anot[,'sgRNA1_ID'],anot[,'sgRNA2_ID'],sep='')
scaf1=anot[,'Scaffold']=='Wildtype'& anot[,'Notes']!='DistanceCut' & guides_ID%in%guides_ID[scaf2]

#####HAcemos la media para los 3  

mat=cbind(apply(lfc_plm[,c('PILOT_100X_D14_R1','PILOT_100X_D14_R2','PILOT_100X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_500X_D14_R1','PILOT_500X_D14_R2','PILOT_500X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_PCR500X_D14_R1','PILOT_PCR500X_D14_R2','PILOT_PCR500X_D14_R3')],1,mean))
colnames(mat)=c('PILOT_100X_D14','PILOT_500X_D14','PILOT_PCR500X_D14')

pal1=c(rgb(108/255,116/255,134/255),rgb(105/255,156/255,201/255),rgb(205/255,232/255,236/255))
colores=pal1[c(1,1,1,2,2,2,3,3,3,2,2,2,3,3,3,2,2,2,3,3,3)]
pos=c('essential + non-targeting','non-targeting + essential','intergenic + essential','essential + intergenic')
neg=c('non-essential + non-essential','non-targeting + non-targeting','intergenic + intergenic')

pdf('/pathFortheFigure/Fig2ASupp_boxplot_D3.pdf',width=4,height=5)
boxplot(mat[anot[,'Notes']!='DistanceCut',1],
        mat[anot[,'Notes']!='DistanceCut',2],
        mat[anot[,'Notes']!='DistanceCut',3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%pos,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%pos,3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%neg,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%neg,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%neg,3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%pos,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%pos,3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%neg,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%neg,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%neg,3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%pos,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%pos,3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%neg,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%neg,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%neg,3],
        col=alpha(colores,1),xlab='sgRNA Scaffold',ylab='log2FC',las=3)
dev.off()

############################################################################################################################
############################################################################################################################
######### Supp figure 2B


scaf1=anot[,'Scaffold']=='Wildtype'& anot[,'Notes']!='DistanceCut'
scaf2=anot[,'Scaffold']=='Modified6'& anot[,'Notes']!='DistanceCut'
scaf3=anot[,'Scaffold']=='Modified7'& anot[,'Notes']!='DistanceCut'

breaksList = seq(0.5, 1, length.out = 100)

colnames(lfc_plm)=c("100X_R1","100X_R2","100X_R3",
                    "500X_R1","500X_R2","500X_R3",
                    "PCR500X_R1","PCR500X_R2","PCR500X_R3")

pdf('/pathFortheFigure/Fig2Bsup_corr_D3.pdf',width=5,height=5)

pheatmap(cor(lfc_plm),col=colorRampPalette(brewer.pal(5,'Reds'))(100),
         fontsize = 7,cellwidth = 7, cellheight = 7,
         breaks = breaksList)

dev.off()


#####Positional effect Supp fig 2C

anot=as.matrix(read.delim('/PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']

anot[anot[,'Gene2']=='','Gene2']=anot[anot[,'Gene2']=='','sgRNA2_Class']
anot[anot[,'Gene1']=='','Gene1']=anot[anot[,'Gene1']=='','sgRNA1_Class']

lfc_plm=readRDS('/PILOT/MERGED/PILOC_EXACT_LFC_D3.rds')

#####HAcemos la media para los 3  (aunque usaremos 500X)

mat=cbind(apply(lfc_plm[,c('PILOT_100X_D14_R1','PILOT_100X_D14_R2','PILOT_100X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_500X_D14_R1','PILOT_500X_D14_R2','PILOT_500X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_PCR500X_D14_R1','PILOT_PCR500X_D14_R2','PILOT_PCR500X_D14_R3')],1,mean))
colnames(mat)=c('PILOT_100X_D14','PILOT_500X_D14','PILOT_PCR500X_D14')

####Agregamos con todos juntos

####Vamos a agregar

sel1=anot[,'Scaffold']=='Wildtype'
sel2=anot[,'Scaffold']=='Modified6'
sel3=anot[,'Scaffold']=='Modified7'

temp=data.frame(G1=anot[sel1,'Gene1'],G2=anot[sel1,'Gene2'],lfc_500x=mat[sel1,'PILOT_500X_D14'],lfc_100x=mat[sel1,'PILOT_100X_D14'],lfc_PCR500x=mat[sel1,'PILOT_PCR500X_D14'])
temp <- aggregate(. ~ G1 + G2, data=temp, mean)
temp=as.matrix(temp)
wt=temp
temp=data.frame(G1=anot[sel2,'Gene1'],G2=anot[sel2,'Gene2'],lfc_500x=mat[sel2,'PILOT_500X_D14'],lfc_100x=mat[sel2,'PILOT_100X_D14'],lfc_PCR500x=mat[sel2,'PILOT_PCR500X_D14'])
temp <- aggregate(. ~ G1 + G2 , data=temp, mean)
temp=as.matrix(temp)
m6=temp

temp=data.frame(G1=anot[sel3,'Gene1'],G2=anot[sel3,'Gene2'],lfc_500x=mat[sel3,'PILOT_500X_D14'],lfc_100x=mat[sel3,'PILOT_100X_D14'],lfc_PCR500x=mat[sel3,'PILOT_PCR500X_D14'])
temp <- aggregate(. ~ G1 + G2 , data=temp, mean)
temp=as.matrix(temp)
m7=temp


rownames(wt)=paste(wt[,1],wt[,2],sep='_')
temp=wt[,c("lfc_500x","lfc_100x","lfc_PCR500x")]
rownames(temp)=paste(wt[,2],wt[,1],sep='_')
wt=wt[rownames(wt)%in%rownames(temp),]
wt=cbind(wt,temp[rownames(wt),])
wt=wt[,3:8]
colnames(wt)=c('lfc_500x_1','lfc_100x_1','lfc_PCR500x_1',
               'lfc_500x_2','lfc_100x_2','lfc_PCR500x_2')
temp=rownames(wt)
wt=apply(wt,2,as.numeric)
rownames(wt)=temp

rownames(m6)=paste(m6[,1],m6[,2],sep='_')
temp=m6[,c("lfc_500x","lfc_100x","lfc_PCR500x")]
rownames(temp)=paste(m6[,2],m6[,1],sep='_')
m6=m6[rownames(m6)%in%rownames(temp),]
m6=cbind(m6,temp[rownames(m6),])
m6=m6[,3:8]
colnames(m6)=c('lfc_500x_1','lfc_100x_1','lfc_PCR500x_1',
               'lfc_500x_2','lfc_100x_2','lfc_PCR500x_2')
temp=rownames(m6)
m6=apply(m6,2,as.numeric)
rownames(m6)=temp

rownames(m7)=paste(m7[,1],m7[,2],sep='_')
temp=m7[,c("lfc_500x","lfc_100x","lfc_PCR500x")]
rownames(temp)=paste(m7[,2],m7[,1],sep='_')
m7=m7[rownames(m7)%in%rownames(temp),]
m7=cbind(m7,temp[rownames(m7),])
m7=m7[,3:8]
colnames(m7)=c('lfc_500x_1','lfc_100x_1','lfc_PCR500x_1',
               'lfc_500x_2','lfc_100x_2','lfc_PCR500x_2')
temp=rownames(m7)
m7=apply(m7,2,as.numeric)
rownames(m7)=temp

pdf('pathFortheFigure/FigSupp2C_all_D3.pdf',width=10,height=10)
par(mfrow=c(3,3))

###WT
data <- data.frame(sgRNA1 = wt[,'lfc_500x_1'], sgRNA2 = wt[,'lfc_500x_2'])
# Fit linear model
reg<-lm(sgRNA1 ~ sgRNA2, data = data)
# Scatter plot
plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='Wildtype (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = wt[,'lfc_100x_1'], sgRNA2 = wt[,'lfc_100x_2'])
reg<-lm(sgRNA1 ~ sgRNA2, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='Wildtype (X100-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = wt[,'lfc_PCR500x_1'], sgRNA2 = wt[,'lfc_PCR500x_2'])
reg<-lm(sgRNA1 ~ sgRNA2, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='Wildtype (PCR500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

###M6
data <- data.frame(sgRNA1 = m6[,'lfc_500x_1'], sgRNA2 = m6[,'lfc_500x_2'])
# Fit linear model
reg<-lm(sgRNA1 ~ sgRNA2, data = data)
# Scatter plot
plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M6 (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m6[,'lfc_100x_1'], sgRNA2 = m6[,'lfc_100x_2'])
reg<-lm(sgRNA1 ~ sgRNA2, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M6 (X100-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m6[,'lfc_PCR500x_1'], sgRNA2 = m6[,'lfc_PCR500x_2'])
reg<-lm(sgRNA1 ~ sgRNA2, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M6 (PCR500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))


###M7
data <- data.frame(sgRNA1 = m7[,'lfc_500x_1'], sgRNA2 = m7[,'lfc_500x_2'])
# Fit linear model
reg<-lm(sgRNA1 ~ sgRNA2, data = data)
# Scatter plot
plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M7 (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m7[,'lfc_100x_1'], sgRNA2 = m7[,'lfc_100x_2'])
reg<-lm(sgRNA1 ~ sgRNA2, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M7 (X100-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m7[,'lfc_PCR500x_1'], sgRNA2 = m7[,'lfc_PCR500x_2'])
reg<-lm(sgRNA1 ~ sgRNA2, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M7 (PCR500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.5),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

dev.off()


####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for SUPP FIGURE3
####### ####### ####### ####### ####### ####### ####### ####### 

#######  Supp fig 3A A Gini index

colo1=as.matrix(read.delim('/input/COLO1_RUNMERGED_EXACT_ANNOTATED.txt'))

selec=c("lib.COLO.1",
        "SIDM00136_CPID2437","SIDM00136_CPID2440","SIDM00136_CPID2443",
        "SIDM00136_CPID1020","SIDM00136_CPID1023","SIDM00136_CPID1026")

gini_plm=counts[,'lib.COLO.1']

library("ineq")

##curve only plasmid
plm=cbind(counts[counts[,'lib.COLO.1']>0,'lib.COLO.1'],0,0)
plm=plm[order(plm[,1]),]
colnames(plm)=c('ori','read','sgRNA')

plm[,'read']=cumsum(plm[,'ori'])/sum(plm[,'ori'])
plm[,'sgRNA']=1:nrow(plm)
plm[,'sgRNA']=plm[,'sgRNA']/max(plm[,'sgRNA'])

pdf('/pathFortheFigure/FigSupp3A_plasmidCOLO_reads.pdf',width=5,height=5)

plot(plm[,'sgRNA'],plm[,'read'],type='l',col=rgb(0.2,0.7,0.4),lwd=2,
     xlab='Fraction of total sgRNAs',ylab='Fraction of total reads')
legend('topleft',paste('GINI:',round(ineq::Gini(counts[counts[,'lib.COLO.1']>0,'lib.COLO.1']),2)))
points(c(0,1),c(0,1),type='l')
dev.off()

##### Pannel B stats ll
stat=as.matrix(read.delim('/Users/ih6/projects/Basset/2Gpaper/ENCORE_data/colo1/All.stats'))
plmd=as.matrix(read.delim('/Users/ih6/projects/Basset/2Gpaper/ENCORE_data/lib-COLO-1.stats',sep=' '))

colnames(plmd)=c('fet','plasmid','plasmid_total','plasmid_rel')

length(grep('.1',colnames(stat)))

colnames(stat)[grep('\\.1',colnames(stat))]=
  
  temp=colnames(stat)[grep('\\.1',colnames(stat))][1:5]

colnames(stat)=gsub('\\.1','_total',colnames(stat))
colnames(stat)=gsub('\\.2','_rel',colnames(stat))


temp=stat[,1]
stat=cbind(apply(plmd[,2:4],2,as.numeric),
           apply(stat[,grep("SIDM00136_CPID2437",colnames(stat))],2,as.numeric),
           apply(stat[,grep("SIDM00136_CPID2440",colnames(stat))],2,as.numeric),
           apply(stat[,grep("SIDM00136_CPID2443",colnames(stat))],2,as.numeric),
           apply(stat[,grep("SIDM00136_CPID1020",colnames(stat))],2,as.numeric),
           apply(stat[,grep("SIDM00136_CPID1023",colnames(stat))],2,as.numeric),
           apply(stat[,grep("SIDM00136_CPID1026",colnames(stat))],2,as.numeric))
rownames(stat)=temp

vect=stat[c('MATCH','SWAP'),grep('_rel',colnames(stat))]

pdf('/pathFortheFigure/FigSupp3B_plasmidCOLO_reads.pdf',width=5,height=5)

x=barplot(vect*100,beside=T,ylim=c(0,100),las=3,legend=c('gRNA1+gRNA2','SWAP'),
          names=c('plasmid',
                  'cas9_neg (r1)','cas9_neg (r2)','cas9_neg (r3)',
                  'HT29 (r1)','HT29 (r2)','HT29 (r3)'),ylab='HT29 reads (% total reads)')
text(x[1,],vect[1,]*100,round(vect[1,]*100,2),pos=4,srt=90)
text(x[2,],vect[2,]*100,round(vect[2,]*100,2),pos=4,srt=90)
dev.off()

############
#### Pannel C

score_ag=readRDS('/inpout/SCORE_median_lfc.rds')

colo1=readRDS('/input/colo1_c91.rds')

IDs_colo1=readRDS('/input/IDs_colo1.rds')
sel=IDs_colo1[,'Note']%in%c('AnchorSingletons','LibrarySingletons','NegativeControls','PositiveControls')
IDs_colo1=IDs_colo1[sel,]
colo1=colo1[sel,]

#####Agregamos ahora nuestros datos 

colo1_ag=data.frame(gene=c(IDs_colo1[IDs_colo1[,'Note']%in%c('AnchorSingletons'),'Gene2'],
                           IDs_colo1[IDs_colo1[,'Note']%in%c('LibrarySingletons'),'Gene1']),
                    lfc1=c(colo1[IDs_colo1[,'Note']%in%c('AnchorSingletons'),1],
                           colo1[IDs_colo1[,'Note']%in%c('LibrarySingletons'),1]),
                    lfc2=c(colo1[IDs_colo1[,'Note']%in%c('AnchorSingletons'),2],
                           colo1[IDs_colo1[,'Note']%in%c('LibrarySingletons'),2]),
                    lfc3=c(colo1[IDs_colo1[,'Note']%in%c('AnchorSingletons'),3],
                           colo1[IDs_colo1[,'Note']%in%c('LibrarySingletons'),3]),row.names=NULL)
colo1_ag= aggregate(. ~ gene, data=colo1_ag, median)
colo1_ag=as.matrix(colo1_ag)
temp=colo1_ag[,1]
colo1_ag=apply(colo1_ag[,2:4],2,as.numeric)
rownames(colo1_ag)=temp

colo1_ag=colo1_ag[rownames(colo1_ag)%in%names(score_ag),]
score_ag=score_ag[names(score_ag)%in%c(rownames(colo1_ag))]

colo1_ag=apply(colo1_ag,1,median)

colo1_ag=cbind(colo1_ag,score_ag[names(colo1_ag)])

x=c(colo1_ag[,1])
y=c(colo1_ag[,2])

##500X Wildtype scaffold
data <- data.frame(X = x, Y = y)
# Fit linear model
reg<-lm(Y ~ X, data = data)

pdf("/pathFortheFigure/FigSupp3C_COLOvsSCORE.pdf",height=5, width=5)

plot(data,xlim=c(-4,2),ylim=c(-4,2),
     main='Gene level correlation single KO',xlab='logFC(ENCORE-HT29)',ylab='logFC(Score-HT29)',
     col=rgb(0.2,0.5,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))
dev.off()















