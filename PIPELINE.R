
############dbRDA##############
library(ggplot2)
library(vegan)
library(tidyr)
library(dplyr)
library(BiodiversityR)
library(ggVennDiagram)
library(pheatmap)

#Import data
# to apply for 16S or 18S replace the table with the data 16S or 18S
tab=read.csv(file = "microscopy.csv",header=T,row.names="id",sep=";")
tab=tab[,1:59] 
tab<-t(tab)
tab <- as.data.frame(tab/rowSums(tab))
rowSums(tab)

rownames(tab) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40","41", "42", "43", "44", "45", "46", "47", "48", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60")

tab=t(tab)

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
  
  si$Season <- c("cheia", "cheia", "seca", "seca", "seca", "seca", "seca", "seca", "cheia", "cheia", "cheia", "cheia", "cheia", "cheia", "seca", "seca", "seca", "seca", "seca", "seca", "cheia", "cheia", "cheia", "cheia", "cheia", "cheia", "seca", "seca", "seca", "seca", "seca", "seca", "cheia", "cheia", "cheia", "cheia", "cheia", "cheia", "seca", "seca", "seca", "seca", "seca", "seca", "cheia", "cheia", "cheia", "cheia", "cheia", "seca", "seca", "seca", "seca", "seca", "seca", "cheia", "cheia", "cheia", "cheia")
  
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
### 
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
 # to apply for 16S or 18S replace the table with the data 16S or 18S
 #Import data
  t=read.csv(file ="microscopy.csv",
             header=T,sep=";",row.names="id")
  t=t[,1:59] 
  x=t(t)

  #apply diversity index in Microscopy
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
##### SHANNON
graphic=read.csv(file ="diversity_index.csv.csv",header=T,sep=";")
  
{
  hlinha= cbind.data.frame(graphic$Shannon.M, graphic$Shannon.16S, 
                           graphic$Shannon.18S, graphic$site)
  
  colnames(hlinha)=c("Microscopy", "16S","18S","site")

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
  
  r=pivot_longer(r, names_to = "index",values_to = "Richness_Index", -site)
  
  ggplot(r, aes(x=index,y=Richness_Index))+
    geom_boxplot(fill="slateblue", alpha=0.7) +
    theme_minimal ()+
    theme(axis.line = element_line(colour = "grey60", linewidth = 0.3)) 
  #xlab("cyl")
  #geom_jitter(width = 0.5, aes(color = site))
}
#####EEVENNES
{
  e=cbind.data.frame(graphic$Eevenness.M, graphic$Eevenness.16S, 
                     graphic$Eevenness.18S, graphic$site)
  
  colnames(e)=c("Microscopy", "16S","18S","site")
  
  e=pivot_longer(e, names_to = "index",values_to = "Evenness_Index", -site)
  
  ggplot(e, aes(x=index,y=Evenness_Index))+
    geom_boxplot(fill="slateblue", alpha=0.7)+
    theme_minimal ()+
    theme(axis.line = element_line(colour = "grey60", linewidth = 0.3))  
  #xlab("cyl")
  #geom_jitter(width = 0.5, aes(color = site))
}
  
########Statistical comparison of diversity indices#########
  #dataframe 
  alfa = cbind.data.frame(result16$richness, result16$Shannon, result16$Simpson, result16$Eevenness,
                          result18$richness, result18$Shannon, result18$Simpson, result18$Eevennes,
                          resultM$richness, resultM$Shannon, resultM$Simpson, resultM$Eevenness,
                          result16$site)
  
  colnames(alfa)=c("richness_16", "shannon_16","simpson_16", "evenness_16", "richness_18", "shannon_18", "simpson_18","evenness_18", "richness_m", "shannon_m", "simpson_m","evennesse_m","site")
  
  #long format table
  alfa=pivot_longer(alfa, names_to = "index",values_to = "values", -site)
  
  #  Friedman test
  resultado <- friedman.test(values ~ index | site, data=alfa)
  print(resultado)
  
  # wilcoxon test
  a=pairwise.wilcox.test(alfa$values,alfa$index,paired=T,p.adjust.method = "fdr")
  View(a)
  
  #### Veen Diagram ####
  #tab of items
  Vg=read.csv("TCC/venn.csv", sep=";")
  
  process_region_data(Venn(Vg))
  
  # 3D Venn diagram
  ggVennDiagram(Vg,label_color="black", set_color="black")+
    scale_fill_gradient(low = "#FFF3E0", high = "#F57C00")
  
  scale_fill_gradient(low = "#F5F5F5", high = "#9E9E9E")
  
  scale_fill_gradient(low = "#FFE0B2", high = "#E65100")
  scale_fill_gradient(low = "#E0F7FA", high = "#006064")
  
  exclusive_Vgm <- setdiff(Vg$genero_m, union(Vg$genero_16, Vg$genero_18))
  print(exclusive_Vgm)
  
  exclusive_V18 <- setdiff(Vg$genero_18, union(Vg$genero_16, Vg$genero_m))
  print(exclusive_V18)
  
  exclusive_V16 <- setdiff(Vg$genero_16, union(Vg$genero_18, Vg$genero_m))
  print(exclusive_V16)
  
  intersecçãom16 = intersect(Vg$genero_m, Vg$genero_16)
  print(intersecçãom16)
  
  intersecçãom18 = intersect(Vg$genero_m, Vg$genero_18)
  print(intersecçãom18)
  
  intersecção1618 = intersect(Vg$genero_16, Vg$genero_18)
  print(intersecção1618)
  
  intersecçãom1618 =  intersect(intersect(Vg$genero_16, Vg$genero_18), Vg$genero_m)
  print(intersecçãom1618)
  
  #### z-scrore ####
  # Table by method
  metodo=read.csv("z_score_method.csv", sep=";", row.names="id")
  
  # Global table
  global=read.csv("z-score_global.csv", sep=";")
  
  metodo <- metodo %>%
    mutate(across(
      where(is.character),
      ~ na_if(., "n.d.a")         
    )) %>%
    mutate(across(
      matches("media|sd"),         
      ~ as.numeric(gsub(",", ".", .))  
    ))
  
  dados_z <- metodo %>%
    mutate(
      z_16S = (media_16 - media_global) / sd_global,
      z_18S = (media_18 - media_global) / sd_global,
      Microscopy = (media_m - media_global) / sd_global
    )
  
  dados=dados_z[,9:11]
  dados= t(dados)
  
  dados=read.csv("dados_Z.csv", sep=";", row.names="id")
  
  { 
    #heatmap
    pheatmap(dados,
             cluster_rows = FALSE,     
             cluster_cols = TRUE,        
             main = "Heatmap",
             color = colorRampPalette(c("blue", "white", "red"))(50),
             fontsize_row = 8,
             fontsize_col = 10,
             cellwidth = 7,
             cellheight = 30,
             border_color = NA)
  }
  
  dados_long <- as.data.frame(dados) %>%    
    mutate(Metodo = rownames(dados)) %>%     
    pivot_longer(
      cols = -Metodo,                        
      names_to = "Genero", 
      values_to = "Valor"
    )
  ordem_generos <- dados_long %>%
    group_by(Genero) %>%
    summarise(media_valor = mean(Valor, na.rm = TRUE)) %>%
    arrange(desc(media_valor)) %>%
    pull(Genero)
  
  dados_long$Genero <- factor(dados_long$Genero, levels = ordem_generos)
  
  ggplot(dados_long, aes(x = Genero, y = Metodo, fill = Valor)) +
    geom_tile(color = NA, width = 1, height = 1) + 
    scale_fill_gradient2(
      low = "#4169E1", mid = "white", high = "red",
      na.value = "grey85", midpoint = 0, name = "Valor"
    ) +
    coord_fixed(ratio = 5) +                       
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = "Heatmap")
  
  
