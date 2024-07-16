# Supplemental Code for "Spatial distribution and determinants of childhood vaccination refusal in the United States"

Authors: Bokgyeong Kang, Sandra Goldlust, Elizabeth Lee, John Hughes, Shweta Bansal, and Murali Haran

We provide instructions for implementing non-spatial and spatial zero-inflated negative binomial (ZINB) regression models.

## Required packages:

The code has been tested with R version 4.2.2, "Innocent and Trusting." The following R packages must be installed before the code will run successfully:

-   `tidyverse`

-   `ngspatial`

-   `foreach`

-   `MASS`

-   `egg`

-   `nimble`

## Non-spatial ZINB regression model

Before running any code, ensure the required R packages have been installed.

### Simulate data

`/nszinb/data.R`

-   Generate a $m \times m$ lattice from a non-spatial ZINB regression model

-   All components are saved in `/nszinb/data/sim.RData`

```         
<img src="/nszinb/fig/simData.png" width="700">
```

### Fit the model and summarize the results
