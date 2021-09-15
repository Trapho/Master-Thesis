suppressMessages(require(tidyverse))
suppressMessages(require(dplyr))
suppressMessages(require(Matrix))
library(ggplot2)
library(EnvStats)

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

an.c[,-1] = log2(an.c[,-1])

set.seed(1337)
for(j in 1:dim(an.c)[1]){
  an.c[j,which(an.c[j,]==-Inf)] = rnorm(sum(an.c[j,]==-Inf),12,0.5) 
}

an.c$m = apply(an.c[2:41],1,mean)
an.c$m2 = apply(an.c[42:81],1,mean)

pgm = ggplot(an.c,aes(m,R1_40))+
  geom_point()+
  geom_smooth(method = "lm")

pgm

#scal = function(an.c,stand){
#  slopes=data.frame(Index= 1:dim(an.c)[1], "1" = lm(formula = m~R1_1, data = an.c)$coefficients[2])
#  for(b in 3: 81){
#    e = data.frame(A=an.c$m, B= an.c[,b])
#    slopes$A= lm(formula = A~B, data = e)$coefficients[2]
#    colnames(slopes)[b]=paste0(b-1,"")
#  }
#  an.c[,2:81]=an.c[,2:81]*(slopes[,stand]/slopes[,2:81])
#  an.c
#}
#an.c=scal(an.c,2)
#we dont have to scal
###Now the real challange begins
## calculate the gene specific fluctuations
## according to tusher et al. 
## First we will do a test run by comparing the pos ctrl(40) to the neg ctrl(39)

outliers = function(x, control, an.c, outlier=TRUE) {  # x = c(1:40) druch apply, control = 82 (immer)

  an.c =an.c
  
fluc = data.frame("A"=an.c[,control], "B"=an.c[,control+1], "C" = an.c[,x*2], "D" = an.c[,x*2+1])    

rediff = function(a,b,d,e,variation=0,opt=FALSE){     # das ist die function aus tusher et al
  xi=mean(c(a,b))
  xu=mean(c(d,e))
  
  m=sum((c(a,b)-xi)^2)
  n=sum((c(d,e)-xu)^2)
  s=sqrt(0.5*(m+n))
  
  di=(xi-xu)/(s+variation)
  
  ifelse(opt==FALSE,di,s)
}
#wichtig di40 ist einfach di (war nur zu faul das zu ändern)

di40= apply(fluc,1,function(x) rediff(x[3],x[4],x[1],x[2],variation = 1))
#si40= apply(fluc,1,function(x) rediff(x[1],x[2],x[3],x[4],variation = 1, opt=TRUE))
#si40 =apply(fluc,1,function(x) rediff(x[1],x[2],x[3],x[4],variation=5,opt=TRUE))

#plot(si40,di40)
#plot(1:3475,di40)
boxplot(di40, main =paste0(21,""))

out <- boxplot.stats(di40)$out
out_ind <- which(di40 %in% c(out))
#out_ind
test <- rosnerTest(di40,
                   k = length(out_ind)
)
#test
print(paste0(x, "done"))

ol = test[["all.stats"]][["Obs.Num"]][which(test[["all.stats"]][["Outlier"]] ==TRUE & test[["all.stats"]][["Value"]] > mean(di40))]

#ifelse(outlier== TRUE, ol, di40)
ol
} #end of function outlier


outlier = lapply(1:40,function(x) outliers(x=x, control = 82, an.c = an.c, outlier = TRUE))
outlier_order = unique(unlist(outlier)[order(unlist(outlier))])
di = lapply(1:40,function(x) outliers(x=x, control = 82, an.c = an.c, outlier = TRUE))  # achtung funktion muss manuel geändert werden

di.df = as.data.frame(di)
colnames(di.df) = 1:40
di.df$PG.Gene = ms$PG.Genes
di.df$PG.Prot = ms$PG.ProteinNames
di.df.filter= di.df[outlier_order,]


#sd(di40)/mean(di40)

#calc_s0 = function(fluc){
  
  calc= function(fluc,x){
    apply(fluc,1,function(x) rediff(x[1],x[2],x[3],x[4],variation = c[x]))
  }
  i=0
  j=3
  fin=FALSE
  
  while(fin == FALSE){
    s0=c(i,i+j,i+j*2,i+j*3)
    di=apply(s0,2,function(x)calc(fluc=fluc))
  
    if(which(min(di))==4){
      i=i+9
    }else{
      if(which(min(di))==3)
    }
  }
} # for one comparison between + and - it can be done by hand.. but later we will need a fuction, like this one

#di40c= apply(fluc,1,function(x) rediff(x[1],x[3],x[2],x[4],variation = 1))
#si40c= apply(fluc,1,function(x) rediff(x[1],x[3],x[2],x[4],variation = 2.5, opt = TRUE))

#plot(1:3475,di40c) 
#sd(di40c)/mean(di40c)

plot(si40,di40)
plot(si40c,di40c)

#jede disc wird jetzt 25 mal permutiert mit den fluctuations von di40c (und ich glaube, bin aber nicht sicher dass dich dafür sd8di40c nehemn muss)

#permu= function(fluc,stadev,durch,repeats){
  fluc2=data.frame(A= fluc[,1]+rnorm(1,durch,stadev),
                   B= fluc[,2]+rnorm(1,durch,stadev),
                   C= fluc[,3]+rnorm(1,durch,stadev),
                   D=fluc[,3]+rnorm(1,durch,stadev))
  di= apply(fluc2,1,function(x) rediff(x[1],x[2],x[3],x[4],variation = 2.5))
  di=data.frame(A=di)
  
  for(t in 2:repeats){
    fluc2=data.frame(A= fluc[,1]+rnorm(1,durch,stadev),
                     B= fluc[,2]+rnorm(1,durch,stadev),
                     C= fluc[,3]+rnorm(1,durch,stadev),
                     D=fluc[,3]+rnorm(1,durch,stadev))
    di2= apply(fluc2,1,function(x) rediff(x[1],x[2],x[3],x[4],variation = 2.5))
    di2=data.frame(A=di2)
    di=cbind(di,di2)
  }
  di
}

#permutations = permu(fluc,sd(si40c),mean(di40c),25)
#permutations$m = apply(permutations,1,mean)

plot(permutations$m,di40)

fluc$index= 1:3475
permutations$index= 1:3475

pgm_det= data.frame(di40=di40, Index=1:3475)
pgm_det$id= FALSE
pgm_det$id[pgm_det$di40<=-2.5] = ms$PG.Genes[pgm_det$di40<=-2.5]
pgm_det$geph =1
pgm_det$geph[450]=5


ggplot(pgm_det, aes(Index,di40,color=id,size=geph))+
  geom_point()





####lets test vulcano

fluc.vulc = fluc
fluc.vulc = data.frame(A = fluc$A^4, B = fluc$B^4, C = fluc$C^4, D = fluc$D^4) 

fold.change = function(fluc.vulc){
  a = log2(mean(c(fluc.vulc[1], fluc.vulc[2])))
  b = log2(mean(c(fluc.vulc[3], fluc.vulc[4])))
  b-a
}
fluc.vulc$fold2 = apply(fluc.vulc,1,fold.change)

y.axis = function(x){
  -log10(t.test(c(x[1],x[2]),c(x[3],x[4]))$p.value)  
}
fluc.vulc$t = apply(fluc.vulc,1,y.axis)

plot(fluc.vulc$fold2, fluc.vulc$t)


