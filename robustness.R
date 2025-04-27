setwd("E:/Users/LEE/Desktop/network/rubustness/bacrubustnss")

##read otu table
otutab<-read.table("bacrubustness.txt",header = T,row.names=1,sep="\t")
otutab[is.na(otutab)]<-0
##keep 12 otus
counts<-rowSums(otutab>0)
otutab<-otutab[counts>=12,]

comm<-t(otutab)
sp.ra<-colMeans(comm)/32958   #relative abundance of each species

###### there are two choices to get the correlation matrix #######
###### choice 1 (slow): calculate correlation matrix from OTU table
cormatrix=matrix(0,ncol(comm),ncol(comm))
for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi <- sapply(1:nrow(comm), function(k){
      ifelse(comm[k,i]>0, comm[k,i], ifelse(comm[k,j]>0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm), function(k){
      ifelse(comm[k,j]>0, comm[k,j], ifelse(comm[k,i]>0, 0.01, NA))
    })
    # 检查方差是否为零
    if(var(log(speciesi)[!is.na(speciesi)]) > 0 & var(log(speciesj)[!is.na(speciesj)]) > 0){
      corij <- cor(log(speciesi)[!is.na(speciesi)], log(speciesj)[!is.na(speciesj)])
    } else {
      corij <- NA # 如果方差为零，则将相关性设置为NA
    }
    
    cormatrix[i,j] <- cormatrix[j,i] <- corij
  }
}

row.names(cormatrix)<-colnames(cormatrix)<-colnames(comm) # if processed using MENAP, OTU order should match in the original OTU table and the correlation matrix downloaded from MENAP.
#--选择相关系数大于0.8的连线
cormatrix2<-cormatrix*(abs(cormatrix)>=0.80) 
#存在某些情况计算不出来相关系数，定义相关为0
cormatrix2[is.na(cormatrix2)]<-0
#-去除自相关的点
diag(cormatrix2)<-0  
#-查看网络边的数量
sum(abs(cormatrix2)>0)/2
#网络中节点的数量
sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.
#去除没有任何相关的节点.
network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
#对应的删除otu表格otu
sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
sum(row.names(network.raw)==names(sp.ra2))  #check if matched

## 鲁棒性评估robustness simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

rm.p.list = seq(0.05,0.2,by=0.05)

rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  #-随机挑选出一定百分比的OTU
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw
  #这些节点和其他节点连接全部去除
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    #网络矩阵乘以物种平均丰度，改变相关性值的大小
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  # 每一个节点的平均链接数
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  remain.percent
}
install.packages("tidyverse")
library(tidyverse)

rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 year=rep(2014,40),treat=rep("Warmed",40))

currentdat<-dat1

write.csv(currentdat,"random_removal_result_Y14_W.csv")

##plot
library(ggplot2)

currentdat$year

p <- ggplot(currentdat[currentdat$weighted=="weighted",], aes(x=Proportion.removed, y=remain.mean, group=treat, color=treat)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()+
  facet_wrap(~year, ncol=3)
ggsave(paste("test1.pdf", sep=""), p, width = 6, height = 5.1)

p1 <- ggplot(currentdat[currentdat$weighted=="unweighted",], aes(x=Proportion.removed, y=remain.mean, group=treat, color=treat)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()
ggsave(paste("test2.pdf", sep=""), p1, width = 6, height = 5.1)

