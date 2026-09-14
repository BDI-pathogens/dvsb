# dvsb

dvsb stands for Disease Variability in Subpopulations from Biomarkers: that's what we estimate.
We use proxy measurements of biomarkers to estimate the prevalence of disease in a population and its subpopulations.

We designed dvsb thinking more specifically about seroprevalence: the proportion of people who are seropositive for some pathogen, having the associated antibody, suggesting previous exposure and some level of current immunity.
The kind of proxy measurements we were thinking of are optical densities / OD values from ELISA experiments.
However dvsb applies more generally to the prevalence of some condition, i.e. the proportion of people who have it currently.
If that more general case describes your application, substitute all of our mentions of: seroprevalence for prevalence, antibody level for biomarker level, serological assay for biomarker assay, seropositive for positive (i.e. having the condition), seronegative for negative.

### Installation

Install dvsb by running the following commands in an R session:

```r
install.packages("pak") # if not already installed
pak::pak("BDI-pathogens/mastiff")
pak::pak("BDI-pathogens/dvsb")
```

To interface to its Stan code, dvsb can use cmdstanr or cmdstan instead of rstan, if you desire.
This requires that you have cmdstan installed.
Do that by running
```r
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
```
(after which you might need to restart your R session).

### What is the statistical model?

See the article

### Usage

See the [vignette](https://BDI-pathogens.github.io/dvsb/articles/estimation_on_simulated_data_simple.html).
