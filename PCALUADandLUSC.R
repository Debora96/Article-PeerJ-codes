install.packages("factoextra")

library(tidyverse) 
library(stats)
library(factoextra) 


#LUAD

GENESLUADSELEC <- read_csv("../data/GENESLUADSELEC.csv")
LUAD_clin <- read_csv("../data/LUAD_clin.csv")

datatime <- subset(LUAD_clin, select= c(submitter_id,time))
GENESLUADSELEC <- full_join(GENESLUADSELEC, datatime, by = "submitter_id")
GENESLUADSELEC <- GENESLUADSELEC[complete.cases(GENESLUADSELEC[ , c('time')]), ] 

GENESLUADSELEC <- na.omit(GENESLUADSELEC)

GENESLUADSELEC[GENESLUADSELEC==-1] <- 0

dados_pca <- GENESLUADSELEC[,c(2:35)]
sapply(dados_pca, sd)

# PCA with correlation matrix (standard variables)
pca_corr <- prcomp(dados_pca, center = TRUE, scale = TRUE)
summary(pca_corr)

fviz_eig(pca_corr)

summary(pca_corr)$rotation

summary(pca_corr)$x

# individuals graph 
fviz_pca_ind(pca_corr,
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Texto não sobreposto
             legend.title = "Representation"
)


# Variables graph
fviz_pca_var(pca_corr,
             col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     
             legend.title = "Contribution"
)

# Graph of variables and individuals
fviz_pca_biplot(pca_corr, repel = TRUE,
                col.var = "#2E9FDF", 
                col.ind = "#696969"  
)

GENESLUADSELEC <- GENESLUADSELEC %>% mutate(
  vital_status = case_when(vital_status == 0 ~ "DEAD",
                          TRUE ~ "ALIVE"),
)

fviz_pca_ind(pca_corr,
             col.ind = GENESLUADSELEC$vital_status, 
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, 
             ellipse.type = "confidence",
             legend.title = "Vital status",
             repel = TRUE
)

# Gráfico das variáveis e indivíduos
fviz_pca_biplot(pca_corr, repel = TRUE,
                col.var = "black",
                col.ind = as.factor(GENESLUADSELEC$vital_status),  
                addEllipses = TRUE,
                legend.title = "Vital status"
)

#LUSC

GENESLUSCSELEC <- read_csv("../data/GENESLUSCSELEC.csv")
LUSC_clin <- read_csv("../data/LUSC_clin.csv")

datatime <- subset(LUSC_clin, select= c(submitter_id,time))
GENESLUSCSELEC <- full_join(GENESLUSCSELEC, datatime, by = "submitter_id")
GENESLUSCSELEC <- GENESLUSCSELEC[complete.cases(GENESLUSCSELEC[ , c('time')]), ] 

GENESLUSCSELEC <- na.omit(GENESLUSCSELEC)

GENESLUSCSELEC[GENESLUSCSELEC==-1] <- 0

dados_pca <- GENESLUSCSELEC[,c(2:33)]
sapply(dados_pca, sd)

# PCA with correlation matrix (standard variables)
pca_corr <- prcomp(dados_pca, center = TRUE, scale = TRUE)
summary(pca_corr)

fviz_eig(pca_corr)

summary(pca_corr)$rotation

summary(pca_corr)$x

# individuals graph 
fviz_pca_ind(pca_corr,
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Texto não sobreposto
             legend.title = "Representation"
)


# Variables graph
fviz_pca_var(pca_corr,
             col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     
             legend.title = "Contribution"
)

# Graph of variables and individuals
fviz_pca_biplot(pca_corr, repel = TRUE,
                col.var = "#2E9FDF", 
                col.ind = "#696969"  
)

GENESLUSCSELEC <- GENESLUSCSELEC %>% mutate(
  vital_status = case_when(vital_status == 0 ~ "DEAD",
                          TRUE ~ "ALIVE"),
)

fviz_pca_ind(pca_corr,
             col.ind = GENESLUSCSELEC$vital_status, 
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, 
             ellipse.type = "confidence",
             legend.title = "Vital status",
             repel = TRUE
)

# Gráfico das variáveis e indivíduos
fviz_pca_biplot(pca_corr, repel = TRUE,
                col.var = "black", # cor das variáveis
                col.ind = as.factor(GENESLUSCSELEC$vital_status),  
                addEllipses = TRUE,
                legend.title = "Vital status"
)



