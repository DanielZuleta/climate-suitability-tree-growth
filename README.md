# codes-climate-suitability-tree-growth
These files contain the code to perform analyses and generate the figures in the article:

Bernal-Escobar, M., Zuleta, D. and Feeley, K.J. (2022), Changes in the climate suitability and growth rates of trees in eastern North America. Ecography e06298. https://doi.org/10.1111/ecog.06298

**Steps & file description:** 

(1) Download most updated data: This article is entirely based on published data that are available for free and public download; no new data was used in this study. See: Data availability statement section in the paper. 

(2) run "get_distribution_models.R" to fit distribution models for each species .

(3) run "get_growth_and_suitability.R" to clean ring data, estimate growth, and merge dataset with the predicted climate suitability.

(4) run "get_models_and_figures.R" to perform statistical analyses and generate figures in the paper.
