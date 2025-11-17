---
title: 'Analysis with `RegEnrich`'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- 

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- 

::::::::::::::::::::::::::::::::::::::::::::::::

## Explore the data with `RegEnrich`

The TFs included in the package are just for humans, so yopu will need to download the mouse ones here: https://tools.sschmeier.com/tcof/home/

*Note: cam create a challenge that will load `data(TFs)` and see why it doesn't work*


``` r
mouseTFs <- read.csv('data/BrowseTF  TcoF-DB.csv')

logcounts <- filteredcounts[,4:15]
rownames(logcounts) <- filteredcounts$ENTREZID
logcounts <- cpm(logcounts,log=TRUE)

# Define design and contrast
design = model.matrix(~ factordata$CellTypeStatus)
contrast = c(-1, 1,0,0,0,0) # Example contrast

# Initialise a RegenrichSet object
object = RegenrichSet(expr = logcounts,
                      colData = factordata,
                      reg = unique(mouseTFs$GeneID), # regulators
                      method = "limma", # differential expression analysis method
                      design = design, # design model matrix
                      contrast = contrast, # contrast
                      networkConstruction = "COEN", # network inference method
                      enrichTest = "FET") # enrichment analysis method

print(object)


# Perform RegEnrich analysis
set.seed(123)


# This step takes a while
# this step should be last or the results should be given to participants beforehand
object = regenrich_diffExpr(object) %>%
  regenrich_network() %>%
  regenrich_enrich() %>%
  regenrich_rankScore()

# Obtain results
res = results_score(object)
print(res)

# Visualise regulator-target expression
plotRegTarExpr(object, reg = "71371")
```


::::::::::::::::::::::::::::::::::::: keypoints 

-

::::::::::::::::::::::::::::::::::::::::::::::::

