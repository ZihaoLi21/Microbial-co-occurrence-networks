install.packages("ggpmisc")
install.packages("conquer")
install.packages("ggplot2")
install.packages("tidyverse")
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")
install.packages("jsonlite", type = "source")
library(devtools)
#devtools::install_github('zdk123/SpiecEasi')
library(SpiecEasi)
library(igraph)
remove__ <- function(s, patten='__') {
if (grepl(patten, s)) {
s <- strsplit(s, patten)[[1]][2]
}
return(s)
}
setwd("D:/360MoveData/Users/Administrator/Desktop/zxmx")
setwd("D:/360MoveData/Users/Administrator/Desktop/zxmx")
setwd("D:/360MoveData/Users/Administrator/Desktop")
D=data(cars)
View(cars)
View(cars)
View(cars)
View(cars)
D=data(mycars）
D=data(mytcars）
D=data(mtcars）
D=data("mtcars"）
data(mtcars）
library(ggplot2)
data(mtcars）
d=data(mtcars）
# Bar chart example
p <- ggplot(mtcars, aes(factor(cyl)))
# Default plotting
p + geom_bar()
data(mtcars）
data(mtcars）
dd=data(mtcars）
dd=data("mtcars"）
# Bar chart example
ggplot(mtcars, aes(x=cyl,y=count)+geom_point()
# Bar chart example
ggplot(mtcars, aes(x=cyl,y=count))+geom_point()
# Bar chart example
ggplot(mtcars, aes(x="cyl",y="count"))+geom_point()
View(cars)
# Bar chart example
ggplot(cars, aes(x="speed",y="dist"))+geom_point()
# Default plotting
p + geom_bar()
library(ggplot2)
# Bar chart example
ggplot(cars, aes(x="speed",y="dist"))+geom_point()
# Bar chart example
ggplot(cars, aes(x=speed,y=dist))+geom_point()
# Bar chart example
ggplot(cars, aes(x=speed,y=dist,color="red"))+geom_point()
# Bar chart example
ggplot(cars, aes(x=speed,y=dist,color="red",size="8"))+geom_point()
# Bar chart example
ggplot(cars, aes(x=speed,y=dist,color="red",size="8"))+geom_point()+gemo_smooth()
# Bar chart example
ggplot(cars, aes(x=speed,y=dist,color="red",size="8"))+geom_point()+geom_smooth()
# Bar chart example
ggplot(cars, aes(x=speed,y=dist,color="red",size="8"))+geom_point()+
geom_smooth()+theme_bw()
# Bar chart example
ggplot(cars, aes(x=speed,y=dist,color="red",size="8"))+geom_point()+
geom_smooth()+theme_bw()+theme(legend.position = "none")
# Bar chart example
ggplot(cars, aes(x=speed,y=dist,color="red",size="8"))+geom_point()+
geom_smooth()+theme_bw()+theme(legend.position = "up")
# Bar chart example
ggplot(cars, aes(x=speed,y=dist,color="red",size="8"))+geom_point()+
geom_smooth()+theme_bw()+theme(legend.position = "none")
# Bar chart example
ggplot(cars, aes(x=speed,y=dist,color="red",size="dist"))+geom_point()+
geom_smooth()+theme_bw()+theme(legend.position = "none")+
# Default plotting
p + geom_bar()
# Bar chart example
ggplot(cars, aes(x=speed,y=dist,color="red",size=dist))+geom_point()+
geom_smooth()+theme_bw()+theme(legend.position = "none")+
# Default plotting
p + geom_bar()
library(devtools)
devtools::install_github('zdk123/SpiecEasi')
library(SpiecEasi)
library(igraph)
remove__ <- function(s, patten='__') {
if (grepl(patten, s)) {
s <- strsplit(s, patten)[[1]][2]
}
return(s)
}
#设置路径
setwd("E:/desktop/anosim")
#加载vegan
library(vegan)
#导入数据
otu = read.table("Fungi.txt", header=T, row.names= 1, sep="\t")
#加载分组文件
design = read.table("group.txt", header=T, row.names= 1, sep="\t")
#构建检验函数
pairwise.anosim <-function(x,factors)
{
library(vegan)
co = as.matrix(combn(unique(factors),2))
pairs = c()
R = c()
p.value = c()
for(elem in 1:ncol(co)){
ad = anosim(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ,
factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] );
pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
R = c(R,ad$statistic);
p.value = c(p.value,ad$signif)}
pairw.res = data.frame(pairs,R,p.value)
return(pairw.res)
}
#检验
anosim_result <- pairwise.anosim(t(otu), design$Group)
#导入数据
otu = read.table("Fun2.txt", header=T, row.names= 1, sep="\t")
#加载分组文件
design = read.table("group.txt", header=T, row.names= 1, sep="\t")
#构建检验函数
pairwise.anosim <-function(x,factors)
{
library(vegan)
co = as.matrix(combn(unique(factors),2))
pairs = c()
R = c()
p.value = c()
for(elem in 1:ncol(co)){
ad = anosim(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ,
factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] );
pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
R = c(R,ad$statistic);
p.value = c(p.value,ad$signif)}
pairw.res = data.frame(pairs,R,p.value)
return(pairw.res)
}
#检验
anosim_result <- pairwise.anosim(t(otu), design$Group)
#结果书写当前文件夹
write.table(anosim_result,"anosim_result.txt", quote = F, row.names = T, col.names = T, sep = "\t")
#作图部分
nmds <- metaMDS(t(otu), distance = 'bray', k = 2)
nmds.point <- data.frame(nmds$point)
#匹配分组
group <- design[match(rownames(t(otu)),rownames(design)),]
sample_site <- nmds.point[1:2]
sample_site$names <- rownames(sample_site)
colnames(sample_site)[1:2] <- c('NMDS1', 'NMDS2')
#重组数据框
sample_site <- cbind(sample_site,group)
#加载包
library(ggplot2)
#加颜色
p <- ggplot(data = sample_site, aes(NMDS1, NMDS2,color = group )) +
geom_point(size = 5, alpha = 0.8) +
#scale_shape_manual(values = c(17, 16)) +
scale_color_manual(values = c("#e67e22","#e74c3c","#2980b9","#27ae60","grey","black")) + theme_bw() + stat_ellipse(level = 0.68)
#生成图片
#点击回车
p
#回车
#改颜色
c("#e67e22","#e74c3c","#2980b9","#27ae60","grey","black")) ；如，"#FFB6C1";
#形状加颜色（大组套小组）
p <- ggplot(data = sample_site, aes(NMDS1, NMDS2,color = group , shape = group)) +
geom_point(size = 5, alpha = 0.8) +
scale_shape_manual(values = c(15,15,16,16,17,17)) +
scale_color_manual(values = c("#2980b9","#e74c3c","#2980b9","#e74c3c","#2980b9","#e74c3c")) + theme_bw() + stat_ellipse(level = 0.68)
setwd("E:/desktop/anosim")
library(ggplot2)
setwd("E:/desktop/anosim")
library(ggplot2)
library(vegan)
# 读取NMDS数据文件
df = read.delim("fun2.txt",# 这里读取了网络上的demo数据，将此处换成你自己电脑里的文件
header = T,    # 指定第一行是列名
row.names = 1  # 指定第一列是行名
)
# 读取样本分组数据文件
dfGroup = read.delim("group.txt",row.names = 1)
# NMDS计算
dfNmds<-metaMDS(df,distance="bray",k = 2)
# 绘图前的数据整理
data = data.frame(dfNmds$points)
data$group = dfGroup$Group
install.packages("ggforce")
library(ggplot2)
library(ggforce)
library(vegan)
# 读取NMDS数据文件
df = read.delim("fun2.txt", header = T, row.names = 1)
df = t(df)
# 读取样本分组数据文件
dfGroup = read.delim("group.txt", row.names = 1)
# NMDS计算
dfNmds <- metaMDS(df, distance = "bray", k = 2)
# 绘图前的数据整理
data = data.frame(dfNmds$points)
data$group = dfGroup$Group
# 绘图
ggplot(data, aes(x = MDS1, y = MDS2, color = group, group = group, fill = group)) +
geom_point(size = 2) +
theme_classic() +
geom_ellipse(aes(group = group), alpha = 0.3, level = 0.95) +  # 添加置信区间
geom_text(aes(label = rownames(data)), vjust = 1.5, size = 2, color = "black") +
labs(subtitle = paste("stress=", round(dfNmds$stress, 3), sep = ""))
ggplot(data, aes(x = MDS1, y = MDS2, color = group, group = group, fill = group)) +
geom_point(size = 2) +
theme_classic() +
stat_ellipse(aes(group = group), geom = "polygon", level = 0.95, alpha = 0.3) +  # 使用stat_ellipse添加置信区间
geom_text(aes(label = rownames(data)), vjust = 1.5, size = 2, color = "black") +
labs(subtitle = paste("stress=", round(dfNmds$stress, 3), sep = ""))
