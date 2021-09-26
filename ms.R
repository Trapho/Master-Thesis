suppressMessages(require(tidyverse))
suppressMessages(require(dplyr))
suppressMessages(require(Matrix))
library(ggplot2)
library(EnvStats)
library(cowplot)
theme_set(theme_cowplot())

##load the data
ms = read.delim("data/data.xls.txt", stringsAsFactors = FALSE)

##the names need to be changed
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
an = namechange(ms)

##we only need the important cols
an.c = an[,c(1,
             which(apply(as.matrix(colnames(an)),1,grepl, pattern = "R1")),
             which(apply(as.matrix(colnames(an)),1,grepl, pattern = "R2"))
)]

#bring cols in right order
b= as.double(substring(colnames(an.c)[2:81],4))
an.c[,2:81] = an.c[,(order(b)+1)]
colnames(an.c)[2:81] = colnames(an.c)[(order(b)+1)]

#First set all "filtered" to "0". this is important because right now all numbers are saved as strings and we want to change that.


for(j in 1:dim(an.c)[1]){
  an.c[j,which(an.c[j,]=="Filtered")] = "0" 
}
for(i in 2:dim(an.c)[2]){
  class(an.c[,i]) <- "numeric"
}

# now transfer everything with the log2 expect of the index... 0(which are the filtereds) will become inf
an.c[,-1] = log2(an.c[,-1])

#now the filtered which are now -Inf need numbers to simulate low background. 
set.seed(1337)
for(j in 1:dim(an.c)[1]){
  an.c[j,which(an.c[j,]==-Inf)] = rnorm(sum(an.c[j,]==-Inf),12,0.5) 
}


## if u want to delete columns make a backup of an.c and delete the wanted cols

an.c_save = an.c
an.c = an.c_save[,c(-(7*2),-(7*2+1) )]  


#calculated the mean for R1 and R2 seperately.
an.c$m = apply(an.c[seq(2,78,2)],1,mean)
an.c$m2 = apply(an.c[seq(3,79,2)],1,mean)


# mark outlier red with this function. This function can be used for fig1

figure1 = function(compare ,color1 = "black", color2 = "red", x = 60, file, replicate){

  rot1 =compare-an.c$m
  
  if(replicate == 1){
    rot1 = data.frame(A = compare, B = an.c$m, C = rot1)
  }else{
    rot1 = data.frame(A = compare, B = an.c$m2, C = rot1)
  }
  rot1$D = "black"
  rot1$D[head(order(-rot1$C),x)] = "red"

  pgm = ggplot(rot1,aes(B,A))+
    geom_point(aes(color = D))+
    geom_smooth(method = "lm")+
    scale_color_manual(values = c(
      "black" =color1,
      "red" =color2
    ))+
    theme(legend.position = "none")+
    xlab("log2(mean)")+ylab(paste0("log2"))
  
  svg(paste0("imgs/",file, ".svg"),width=10,height=10)
  pgm
  dev.off()
  pgm
}

#####
figure1(compare = an.c$R1_8,x =65, file = "R1_8", replicate = 1)  
# x is the number of red marked
# type for compare which disc u want to compare
#####

##function to detect every outlier for each disc seperately

outliers = function(x, control, an.c, outlier=TRUE, capt = FALSE, rosner=TRUE) {  
  an.c =an.c
  
  fluc = data.frame("A"=an.c[,control], "B"=an.c[,control+1], "C" = an.c[,x*2], "D" = an.c[,x*2+1])    

  rediff = function(a,b,d,e,variation=0,opt=FALSE){     # this function calculates the relative difference
    xi=mean(c(a,b))
    xu=mean(c(d,e))
  
    m=sum((c(a,b)-xi)^2)
    n=sum((c(d,e)-xu)^2)
    s=sqrt(0.5*(m+n))
  
    di=(xi-xu)/(s+variation)
  
    ifelse(opt==FALSE,di,s)
  }


  di= apply(fluc,1,function(x) rediff(x[3],x[4],x[1],x[2],variation = 1))
  #si = apply(fluc,1,function(x) rediff(x[3],x[4],x[1],x[2],variation = 1, opt=TRUE))

  #if(capt == TRUE){
    #svg(paste0("imgs/","22",".svg"),width=10,height=10)
    #boxplot(di, main =paste0(22,""))
    #dev.off()
  #}


  out <- boxplot.stats(di)$out
  out_ind <- which(di %in% c(out))
  #out_ind
  if(rosner==TRUE){
    test <- rosnerTest(di,
                     k = length(out_ind)
    )
    #test
    print(paste0(x, "done"))

    ol = test[["all.stats"]][["Obs.Num"]][which(test[["all.stats"]][["Outlier"]] ==TRUE & test[["all.stats"]][["Value"]] > mean(di))]
  }else{
    all_ol = data.frame(A = out_ind, B = di[out_ind])
    ol= all_ol$A[all_ol$B > mean(di)]
  }

  if(outlier == TRUE){
    ol
  }else{
    di
  }
} #end of function outlier

### Here the function is applied: Things you can change: 
      #control is 82 - the amount of cols u deleted
      #rosner: set to FALSE if u want a less aggresiv filter, but the results wont be significant anymore
      #capt: could be turned to true in order to generate one figure. but then it must be defined in the function code which fig.
outlier = lapply(1:((dim(an.c)[2]-3)/2),function(x) outliers(x=x, control = 80, an.c = an.c, outlier = TRUE, rosner = TRUE))
outlier_order = unique(unlist(outlier)[order(unlist(outlier))])
di = lapply(1:((dim(an.c)[2]-3)/2),function(x) outliers(x=x, control = 80, an.c = an.c, outlier = FALSE, rosner = TRUE))  # achtung funktion muss manuel geändert werden

##Here the di data is set to the right format
di.df = as.data.frame(di)
colnames(di.df) = 1:((dim(an.c)[2]-3)/2)
di.df$PG.Gene = ms$PG.Genes
di.df$PG.Prot = ms$PG.ProteinNames

##then the di are filtered, dependent on the outlier
di.df.filter= di.df[outlier_order,] 


###this function is for filtering into the other dimension

heaty = function(x, artefacts, capt= FALSE){

  fml = as.double(x[1:((dim(an.c)[2]-3)/2)])

  #if(capt == TRUE){
  #svg("imgs/gephoutlier.svg",width=10,height=10)
  #boxplot(fml)
  #dev.off()
  #}


  out2 <- boxplot.stats(fml)$out
  out2_ind <- which(fml %in% c(out2)&fml>0)
  out2_ind


  if(length(out2_ind) == 0){
    fml=rep(0,((dim(an.c)[2]-3)/2))
  }else{
  fml[-(out2_ind)]=0
  }

  if(artefacts==TRUE){

  begin = function(fml){
    if(fml[t+1]==0){
      0
    }else{
      fml[t]
    }
  }
  end = function(fml){
    if(fml[t-1]==0){
      0
    }else{
      fml[t]
    }
  }

  mid = function(fml){
    if(fml[t-1]==0&fml[t+1]==0){
      0
    }else{
      fml[t]
    }
  }

for(t in 1:(((dim(an.c)[2]-3)/2)-2)){    ## Dont delete the positive control, for this reason they are not involved. for this raason -2
  if(fml[t] != 0){
    if(t==1){
      fml[t]=begin(fml)
    }else{
      if(t==((dim(an.c)[2]-3)/2)-2){
        fml[t]=end(fml)
      }else{
        fml[t]=mid(fml)
      }
    }
  }else{
    fml[t]=0
  }
}

}
  fml
}#end of heaty

#Do u want to delete all signals without neighbor signals: set artefacts = TRUE
di.df.filter2 = as.data.frame(t(apply(di.df.filter,1,function(x)heaty(x=x, artefacts = TRUE))))

di.df.filter2$s = apply(di.df.filter2,1,sum)

di.df.filter3 = di.df.filter2[di.df.filter2$s > 0,]

di.df.filter3.norm = as.data.frame(t(apply(di.df.filter3[,1:((dim(an.c)[2]-3)/2)],1,function(x) x/max(x))))

data = expand.grid(X = 1:((dim(an.c)[2]-3)/2), Y = ms$PG.Genes[as.double(row.names(di.df.filter3))] )
data$Z = Z = as.double(unlist(t(di.df.filter3.norm)))



svg("imgs/heatmap_man_without7_with_artefacts.svg",width=10,height=15)
ggplot(data, aes(X,Y, fill=Z))+
  geom_tile()
dev.off()





###If the data is not normalized: implement the function shown below, directly after the log2 trans and imputation of missing values

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













