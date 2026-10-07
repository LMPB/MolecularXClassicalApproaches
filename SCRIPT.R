

############dbRDA##############
library(vegan)
library(ggplot2)

#Import data
tab=read.csv(file ="1M/mdbRDA_s.csv", header=T, row.names="id", sep=";")

#environmental data
amb <-read.csv(file="dbRDA/BroaMO_envdbrda.csv", header=T, row.names=1, sep=";") 

#standard environmental variables
denv=decostand(amb[,1:15], na.rm=T, method="standardize")
denv=as.data.frame(denv)

tab=t(tab)

# all environmental variables
rdasenv <- rda(tab~ WATER_TEMPERATURE + zm + secchi + zeu + zm_zeu + IZM + DISSOLVED_OXYGEN + do_sat + pH + DOC + CHLOROPHYLL + tss + OSS + inorg_tss + TP, data=denv) 

# testing the collinearity
vif.cca(rdasenv)

#Parsimonious subsets of explanatory variables (based on #forward selection)
names(denv) 
denv <- denv[, c(1,6,7,8,9,10,11,13,15)] 

names(denv)

dbRDA=capscale(tab ~ WATER_TEMPERATURE +IZM + DISSOLVED_OXYGEN+do_sat+pH+DOC+CHLOROPHYLL+OSS+TP, 
               denv, dist="bray")

##
permutest(dbRDA, permutations=999)
##

#extract axis' importance and significance 
anova(dbRDA)
valor = as.data.frame(dbRDA$CCA$eig)
valor2 = (dbRDA$CCA$eig/sum(dbRDA$CCA$eig))*100
head(valor2)

anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
m=anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

RsquareAdj(dbRDA) 

########## graphic
{
plot(dbRDA, type="n", xlim = c(-2, 2), ylim = c(-3, 3), xlab="DIM1 (49.51%)", ylab="DIM2 (18.61%)")
points(scores(dbRDA, display="sites")[1:16:nrow(scores(dbRDA, display="sites")), ], pch=19, col="black", bg="black", cex=0.8)
text(dbRDA, dis="cn", col="red3", cex=1, font=3)
text("p=0.001", x=-3, y=-2)
}

######## identifying parts of the dbRDA and separating the dry and rainy seasons
{
  s=summary(dbRDA)
  sp=as.data.frame(s$species[,1:2])
  si=as.data.frame(s$sites[,1:2])
  bp=as.data.frame(s$biplot[,1:2])
  
  nv <- as.data.frame (si$Season)
  
  si$Season <- c("cheia", "cheia", "seca", "seca", "seca", "seca", "seca", "seca", "cheia", "cheia", "cheia", "cheia", "cheia", "cheia", "seca", "seca", "seca", "seca", "seca", "seca", "cheia", "cheia", "cheia", "cheia", "cheia", "cheia", "seca", "seca", "seca", "seca", "seca", "seca", "cheia", "cheia", "cheia", "cheia", "cheia", "cheia", "seca", "seca", "seca", "seca", "seca", "seca", "cheia", "cheia", "cheia", "cheia", "cheia", "cheia", "seca", "seca", "seca", "seca", "seca", "seca", "cheia", "cheia", "cheia", "cheia")
  
  cols <- c("chuva"="blue","seca"="limegreen")
}

######## graphic with rainy seasons
{
plot(dbRDA, type="n", xlim = c(-2, 2), ylim = c(-3, 3), xlab="DIM1 (49.51%)", ylab="DIM2 (18.61%)")
points(dbRDA, pch=16, col=cols, bg="black", cex=0.8)
text(dbRDA, dis="bp", col="red", cex=0.9, font = 3)
#text(dbRDA,"sites", col="black", cex=0.5, pos=3) #aqui adiciona nome aos pontos
text("p=0.001", x=-3, y=-2)
}

######## 
{
#Calculation of variance explained by dbRDA axes
# Calculate the eigenvalues of the dbRDA
valor = as.data.frame(dbRDA$CCA$eig)
# Calculate the percentage contribution of each axis
valor2 = (dbRDA$CCA$eig / sum(dbRDA$CCA$eig)) * 100
# Calculate the cumulative percentage contribution of the axes
valor_acumuladoM = cumsum(valor2)
# Display the first cumulative values
head(valor_acumuladoM)
}

#############diversity index ##################
library(BiodiversityR)

 #Import data
  t=read.csv(file ="1M/mdbRDAnew_g.csv",
             header=T,sep=";",row.names="id")
  
  #apply diversity index
  x1=diversityresult(x, y=NULL, level = NULL, factor = NULL,index = "richness", method = "each site", sortit = FALSE, digits = 5)
  ## "each site" for alfa-diversity; "pooled" for gamma-diversity
  x2=diversityresult(x, y=NULL, level = NULL, factor = NULL,index = "Shannon", method = "each site", sortit = FALSE, digits = 5)
  x3=diversityresult(x, y=NULL, level = NULL, factor = NULL,index = "Jevenness", method = "each site", sortit = FALSE, digits = 5)
  x4=diversityresult(x, y=NULL, level = NULL, factor = NULL,index = "Simpson", method = "each site", sortit = FALSE, digits = 5)
  x5=diversityresult(x, y=NULL, level = NULL, factor = NULL,index = "Eevenness", method = "each site", sortit = FALSE, digits = 5)
  
  resultM=cbind(x1,x2,x3,x4,x5)
  resultM[,1]=row.names(x1)
  resultM[,2:6]=cbind(x1,x2,x3,x4,x5)
  colnames(resultM)=c("site","richness M","Shannon M","Jevenness M","Simpson M","Eevenness M")

  #boxplot graphic
library(ggplot2)
##### SHANNON
{
  hlinha= cbind.data.frame(graphic$Shannon.M, graphic$Shannon.16S, 
                           graphic$Shannon.18S, graphic$site)
  
  colnames(hlinha)=c("Microscopy", "16S","18S","site")
  library(tidyr)
  #transformando a tabela em formato longo
  hlinha=pivot_longer(hlinha, names_to = "index",values_to = "Shannon_Diversity_Index", -site)
  
  ggplot(hlinha, aes(x=index,y=Shannon_Diversity_Index))+
    geom_boxplot(fill="slateblue", alpha=0.7)+
    theme_minimal ()+
    theme(axis.line = element_line(colour = "grey60", linewidth = 0.3))
}
##### SIMPSON
{
  d=cbind.data.frame(graphic$Simpson.M, graphic$Simpson.16S, 
                     graphic$Simpson.18S, graphic$site)
  
  colnames(d)=c("Microscopy", "16S","18S","site")
  
  #library(tidyr)
  #transformando a tabela em formato longo
  d=pivot_longer(d, names_to = "index",values_to = "Simpson_Diversity_Index", -site)
  
  ggplot(d, aes(x=index,y=Simpson_Diversity_Index))+
    geom_boxplot(fill="slateblue", alpha=0.7) +
    theme_minimal ()+
    theme(axis.line = element_line(colour = "grey60", linewidth = 0.3)) 
  #xlab("cyl")
  #geom_jitter(width = 0.5, aes(color = site))
}
##### RICHNESS 
{
  r=cbind.data.frame(graphic$richness.M, graphic$richness.16S, 
                     graphic$richness.18S, graphic$site)
  
  colnames(r)=c("Microscopy", "16S","18S","site")
  
  #library(tidyr)
  #transformando a tabela em formato longo
  r=pivot_longer(r, names_to = "index",values_to = "Richness_Index", -site)
  
  ggplot(r, aes(x=index,y=Richness_Index))+
    geom_boxplot(fill="slateblue", alpha=0.7) +
    theme_minimal ()+
    theme(axis.line = element_line(colour = "grey60", linewidth = 0.3)) 
  #xlab("cyl")
  #geom_jitter(width = 0.5, aes(color = site))
}
##### EEVENNES
{
  e=cbind.data.frame(graphic$Eevenness.M, graphic$Eevenness.16S, 
                     graphic$Eevenness.18S, graphic$site)
  
  colnames(e)=c("Microscopy", "16S","18S","site")
  
  #library(tidyr)
  #transformando a tabela em formato longo
  e=pivot_longer(e, names_to = "index",values_to = "Evenness_Index", -site)
  
  ggplot(e, aes(x=index,y=Evenness_Index))+
    geom_boxplot(fill="slateblue", alpha=0.7)+
    theme_minimal ()+
    theme(axis.line = element_line(colour = "grey60", linewidth = 0.3))  
  #xlab("cyl")
  #geom_jitter(width = 0.5, aes(color = site))
}
  
#############Statistical comparison of alpha diversity indices ##################
  #dataframe 
  alfa = cbind.data.frame(result16$richness, result16$Shannon, result16$Simpson, result16$Eevenness,
                          result18$richness, result18$Shannon, result18$Simpson, result18$Eevennes,
                          resultM$richness, resultM$Shannon, resultM$Simpson, resultM$Eevenness,
                          result16$site)
  
  colnames(alfa)=c("richness_16", "shannon_16","simpson_16", "evenness_16", "richness_18", "shannon_18", "simpson_18","evenness_18", "richness_m", "shannon_m", "simpson_m","evennesse_m","site")
  
  library(tidyr)
  #long format table
  alfa=pivot_longer(alfa, names_to = "index",values_to = "values", -site)
  
  #  Friedman test
  resultado <- friedman.test(values ~ index | site, data=alfa)
  print(resultado)
  
  # wilcoxon test
  a=pairwise.wilcox.test(alfa$values,alfa$index,paired=T,p.adjust.method = "fdr")
  View(a)