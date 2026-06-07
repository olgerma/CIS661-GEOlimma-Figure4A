rm(list = ls())

library(GEOquery)
library(limma)
library(hgu133plus2.db)
library(AnnotationDbi)
library(pROC)
library(ggplot2)

# Load GSE8052
gse <- getGEO("GSE8052", GSEMatrix = TRUE, getGPL = FALSE)
gset <- gse[[1]]

expr <- exprs(gset)
pheno <- pData(gset)

group <- factor(
  pheno$`DDAST:ch1`,
  levels = c("CONTROL", "CASE")
)

# Convert GPL570 probe IDs to Entrez Gene IDs
probe_ids <- rownames(expr)

entrez_ids <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)

keep <- !is.na(entrez_ids)

expr_entrez <- expr[keep, ]
entrez_ids <- entrez_ids[keep]

expr_entrez <- avereps(
  expr_entrez,
  ID = entrez_ids
)

# Load GEOlimma prior probabilities
load("CIS 661/GEOlimma_prob.rda")
rownames(prop) <- prop$Entrez

# Load GEOlimma source code
source("CIS 661/12859_2020_3932_MOESM8_ESM.r")

# Full-data Limma analysis used to define reference genes
design_full <- model.matrix(~group)

fit_full <- lmFit(expr_entrez, design_full)
fit_full_limma <- eBayes(fit_full)

truth_limma <- topTable(
  fit_full_limma,
  coef = 2,
  number = Inf,
  sort.by = "B"
)

true_genes <- rownames(truth_limma)[1:100]

# Figure 4A replication
set.seed(661)

n_sizes <- c(6, 9, 12, 15)
n_trials_each <- 40

fig4a_all <- data.frame()
trial_counter <- 1

for(n in n_sizes){
  
  for(i in 1:n_trials_each){
    
    case_idx <- sample(which(group == "CASE"), n)
    control_idx <- sample(which(group == "CONTROL"), n)
    
    subset_idx <- c(case_idx, control_idx)
    
    subset_expr <- expr_entrez[, subset_idx]
    subset_group <- group[subset_idx]
    
    design_sub <- model.matrix(~subset_group)
    
    fit_sub <- lmFit(subset_expr, design_sub)
    
    fit_limma <- eBayes(fit_sub)
    
    fit_geo <- eBayesGEO(
      fit_sub,
      proportion_vector = prop[, 1, drop = FALSE]
    )
    
    limma_sub <- topTable(
      fit_limma,
      coef = 2,
      number = Inf,
      sort.by = "B"
    )
    
    geo_sub <- topTable(
      fit_geo,
      coef = 2,
      number = Inf,
      sort.by = "B"
    )
    
    limma_sub$truth <- ifelse(
      rownames(limma_sub) %in% true_genes,
      1, 0
    )
    
    geo_sub$truth <- ifelse(
      rownames(geo_sub) %in% true_genes,
      1, 0
    )
    
    limma_auc <- auc(
      roc(
        limma_sub$truth,
        limma_sub$B,
        quiet = TRUE
      )
    )
    
    geo_auc <- auc(
      roc(
        geo_sub$truth,
        geo_sub$B,
        quiet = TRUE
      )
    )
    
    fig4a_all <- rbind(
      fig4a_all,
      data.frame(
        Trial = trial_counter,
        SubsetSize = paste0("GEOlimma n=", n),
        Limma_AUC = as.numeric(limma_auc),
        GEOlimma_AUC = as.numeric(geo_auc),
        Improvement = as.numeric(geo_auc - limma_auc)
      )
    )
    
    trial_counter <- trial_counter + 1
  }
}

# Summary results


write.csv(
  fig4a_all,
  "Figure4A_results.csv",
  row.names = FALSE
)

# Plot
fig4a_all$SubsetSize <- factor(
  fig4a_all$SubsetSize,
  levels = c(
    "GEOlimma n=6",
    "GEOlimma n=9",
    "GEOlimma n=12",
    "GEOlimma n=15"
  )
)

ggplot(
  fig4a_all,
  aes(
    x = Trial,
    y = Improvement,
    color = SubsetSize
  )
) +
  geom_point(size = 2.3) +
  geom_hline(
    yintercept = 0,
    color = "red",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = c(
      "GEOlimma n=6" = "red",
      "GEOlimma n=9" = "orange",
      "GEOlimma n=12" = "green",
      "GEOlimma n=15" = "blue"
    )
  ) +
  theme_classic() +
  labs(
    x = "trial",
    y = "AUC improvement",
    color = ""
  ) +
  theme(
    legend.position = c(0.18, 0.15),
    legend.background = element_rect(
      fill = "white",
      color = "black"
    ),
    legend.key = element_blank()
  )

ggsave(
  "Figure4A_full_replication.png",
  width = 6,
  height = 5,
  dpi = 300
)

sessionInfo()