title: "R Notebook - Supplementary Figures: Collinearity, Correlation and Variable Ranking"
output: html_notebook
---


```{r}
if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("survival")){install.packages("survival")}
if(!require("glmnet")){install.packages("glmnet")}
if(!require("rms")){install.packages("rms")} # vif - One possible source for a `vif`-function
```

#  Collinearity Analysis with Variance Inflation Factor

```{r}
library(caret)
library(rms)

GENESLUADSELEC <- read_csv("../data/GENESLUADSELEC.csv")
LUAD_clin <- read_csv("../data/LUAD_clin.csv")

datatime <- subset(LUAD_clin, select= c(submitter_id,time))
datagenes <- full_join(GENESLUADSELEC, datatime, by = "submitter_id")
datagenes <- datagenes[complete.cases(datagenes[ , c('time')]), ] 

datagenes <- na.omit(datagenes)

datagenes[datagenes==-1] <- 0


# Selecting the feature
Sign <- c("SI","FAT4","CDH10","PTPRD","RP1L1","ADGRG4","DNAH9","PAPPA2","TNR","KEAP1","DAMTS12","RYR3","APOB","MUC17","PCDH15","PCLO","ANK2","CSMD1","ZNF536","COL11A1","NAV3","FAT3","SPTA1","FLG","XIRP2","KRAS","ZFHX4","USH2A","LRP1B","RYR2","CSMD3","MUC16","TP53","TTN")

genes_vif <- setNames(data.frame(matrix(nrow=0, ncol=3)), c("sign", "gene", "vif"))

sign <- intersect(colnames(datagenes), sign)

luad_sign <- datagenes[, c("time", "vital_status", sign) ]
  
surv_obj <- Surv(time = luad_sign$time, event = luad_sign$vital_status)

mdl1 <- coxph(as.formula(paste0("surv_obj ~ ", paste(sign, collapse= " + "))), data = luad_sign)

  
cvif <- vif(mdl1) 
df_genes_vif <- as.data.frame(cvif)
df_genes_vif$gene  <- rownames(df_genes_vif)
rownames(df_genes_vif) <- NULL
genes_vif <- rbind(genes_vif, df_genes_vif)
print(cvif)
```


```{r}
#Figure S7
pdf("figsupl7.pdf", width = 12, height = 7)
ggplot(data=genes_vif, aes(x=gene, y=cvif)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=round(cvif, digits = 2)), vjust = -0.6, hjust = - 0.3, angle = 45, size=4) +
  theme_minimal() +
   theme(plot.title=element_text(hjust=0.5, size=15), 
        axis.text.x = element_text(angle = 45, hjust=1, size=15),
        axis.text.y = element_text(angle = 0, hjust=1, size=13)) +
    ylim(0, 2.5) +
    xlab("Genes") +
    ylab("Variance Inflation Factors") +
    theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.8), angle = 00))
dev.off()

#rm(list = setdiff(ls(),c("sign", "datagenesLUSC")))
  
```

#  Collinearity Analysis with Correlation

```{r}
if(!require("corrplot")){install.packages("corrplot")}
```


```{r}
dim(luad_sign)

luad_sign %>%
  select_if(is.numeric) %>%
  cor(.) %>%
  as.data.frame() %>%
  mutate(var1 = rownames(.)) %>%
  gather(var2, value, -var1) %>%
  filter(var1 != var2 )  %>%
  arrange(desc(value))
```
```{r}

luad_cor <- cor(luad_sign)

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(luad_sign)
head(p.mat[, 1:5])

```

```{r}
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

pdf("figsupl9.pdf", width = 20, height = 20)
corrplot(luad_cor, method="color", col=col(200),  
         type="lower", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         tl.cex =0.8,
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
         )
dev.off()

#rm(list = setdiff(ls(),c("sign", "data")))

```

# Variable Ranking and Feature Selection

```{r}
if(!require("varrank")){BiocManager::install("varrank")}
```

https://cran.r-project.org/web/packages/varrank/vignettes/varrank.html

Estimates the variables ranks based on mutual information. 

The filter approach based on mutual information is the Minimum Redundancy Maximum Relevance (mRMRe) algorithm.

For the first search, it is advised to use either ''peng'' (faster) or ''estevez'' (more reliable but much slower for large datasets) and, in case the number of variables is large (>100), restrict the ''forward'' search to ''n.var = 100.''

```{r} 
library(varrank)

varrank_df <- datagenesLUSC[, c("vital_status", sign)]

varrank <- varrank(data.df = varrank_df, 
                   #method = "estevez", # more reliable but much slower 
                   method = "peng", # faster
                   variable.important = "vital_status", 
                   discretization.method = "sturges", 
                   n.var = 36,
                   algorithm = "forward", 
                   scheme="mid", 
                   verbose = T)

#summary(varrank)
pdf("../deboravclima/figsupl11.pdf", width = 7, height = 6)
plot(x = varrank, colsep = F, rowsep = F,cellnote = F, labelscex =1)
dev.off()
```
