# The *bbreg* package

The *bbreg* package deals with regression models with response variables being continuous
and bounded. It currently provides implementation of two regression models: beta
regression <https://www.tandfonline.com/doi/abs/10.1080/0266476042000214501> and bessel regression <https://arxiv.org/abs/2003.05157>. For both of these models, the estimation is 
done with the EM-algorithm. The EM-algorithm approach for beta regression
was developed in <https://www.tandfonline.com/doi/abs/10.1080/00949655.2017.1350679>.

## Vignette

The *bbreg* package has a vignette. Check the most recent version of the vignette at <https://rpubs.com/alexandrebsimas/intro-bbreg>

## Installing the *bbreg* package from this repository:

To install the *bbreg* package from this repository, just run the command:

```{r}
devtools::install_github("vpnsctl/bbreg/bbreg")
```

To install the *bbreg* package from this repository **with vignette**, run the command:
```{r}
install_github("vpnsctl/bbreg/bbreg", build_vignettes = TRUE)
```

This repository will always contain the most recent version of the *bbreg* package.
