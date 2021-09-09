suppressMessages(require(tidyverse))
suppressMessages(require(dplyr))
suppressMessages(require(Matrix))
library(ggplot2)

namechange = function(geph){
  a = apply(as.matrix(colnames(geph)),1,grepl, pattern = "Schulte")
  colnames(geph)[a] = gsub(
    ".*Fu_","",colnames(geph)[a]
  )
  colnames(geph)[a] = gsub(
    ".raw.PG.Quantity","",colnames(geph)[a]
  )
  geph
}

ms = read.delim("data/data.xls.txt", stringsAsFactors = FALSE)

an = namechange(ms)
an1 = an[,c(1, which(apply(as.matrix(colnames(an)),1,grepl, pattern = "R1")))]
an1=t(t(an1))

#the following three lines always togheter!!! never alone
a=as.double(substring(colnames(an1)[2:41],4))
an1[,2:41] = an1[,(order(a)+1)]
colnames(an1)[2:41] = colnames(an1)[(order(a)+1)]

for(j in 1:dim(an1)[1]){
  an1[j,which(an1[j,]=="Filtered")] = "0" 
}

class(an1) <- "numeric"

an1 = cbind(an1,matrix(apply(an1[,2:41],1,sum),ncol=1)) 
colnames(an1)[42]="summe"

an1 =an1[order(an1[,42]),]

an1.df = as.data.frame(an1)
an1.df$ds = apply(an1[,2:41],1,mean)
an1.df$Geph= FALSE
an1.df$Geph[an1.df$Index==450] = TRUE

pgm = ggplot(an1.df,aes((R1_39)^(1/4),(R1_40)^(1/4), color= Geph))+
  geom_point()
pgm

an.c = an[,c(1,
             which(apply(as.matrix(colnames(an)),1,grepl, pattern = "R1")),
             which(apply(as.matrix(colnames(an)),1,grepl, pattern = "R2"))
             )]
#alles zusammen
b= as.double(substring(colnames(an.c)[2:81],4))
an.c[,2:81] = an.c[,(order(b)+1)]
colnames(an.c)[2:81] = colnames(an.c)[(order(b)+1)]

for(j in 1:dim(an.c)[1]){
  an.c[j,which(an.c[j,]=="Filtered")] = "0" 
}
for(i in 2:dim(an.c)[2]){
  class(an.c[,i]) <- "numeric"
}

an.c[,-1] = (an.c[,-1])^(1/4)

set.seed(1337)
for(j in 1:dim(an.c)[1]){
  an.c[j,which(an.c[j,]==0)] = rnorm(sum(an.c[j,]==0),5.25,0.5) 
}

an.c$m = apply(an.c[2:81],1,mean)

pgm = ggplot(an.c,aes(m,R1_40))+
  geom_point()+
  geom_smooth(method = "lm")

pgm

scal = function(an.c,stand){
  slopes=data.frame(Index= 1:dim(an.c)[1], "1" = lm(formula = m~R1_1, data = an.c)$coefficients[2])
  for(b in 3: 81){
    e = data.frame(A=an.c$m, B= an.c[,b])
    slopes$A= lm(formula = A~B, data = e)$coefficients[2]
    colnames(slopes)[b]=paste0(b-1,"")
  }
  an.c[,2:81]=an.c[,2:81]*(slopes[,stand]/slopes[,2:81])
  an.c
}
an.c=scal(an.c,2)










