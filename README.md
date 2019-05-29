# uplifteval

## Description
`uplifteval` provides a variety of plots and metrics to evaluate the predictions of uplift models for both binary and continuous outcomes.  Naming conventions are somewhat inconsistent due to the mix of packages that have been ported into this package.

Included with the package are:

| Function      | Description   |
|:------------- |:------------- |
| plot_uplift_guelman  | A copy of the R `uplift` package's Qini chart and metrics. Binary outcomes only. |
| pl_plot | A port of the python pylift package's Qini chart and metrics. Binary outcomes only. |
| plot_uplift | A display of absolute cumulative differential response for bootstrap samples of balanced or unbalanced treatment and control groups, ordered by descending score.  Also displays a histogram of treatment and control groups to confirm if their distributions are balanced and identical.  Both binary and continuous outcomes.|


## Examples
Below is an example using the grf package to estimate uplift via a causal forest.  We then plot the results using three methods: the uplift package, the port of the pylift plot, and an alternative uplift plot. 
```
set.seed(123)

rl <- function(x){
  round(1/(1+exp(-x)))
}
n = 2000; p = 10
X = matrix(rnorm(n*p), n, p)
W = rbinom(n, 1, 0.2)
Y = rl(rl(X[,1]) * W - rl(X[,3]) * W + rnorm(n))
tau.forest = causal_forest(X, Y, W)
tau.hat = predict(tau.forest, X)

# the `uplift` plot
plot_uplift_guelman(tau.hat$predictions, W, Y)

# the port of pylift plot
plue = plUpliftEval(W.test, Y.test, tau.hat$predictions)
pl_plot(plue,
        show_practical_max = TRUE,
        show_theoretical_max = TRUE,
        show_no_dogs = TRUE,
        n_bins=20)

# the alternative uplift plot
plot_uplift(tau.hat$predictions, W, Y)

```

## Contributions
Feedback and contributions are welcome.  Please submit an issue.

## To-do
- when comparing two models, score ranges might be different, consider offering depth option
- add cost of targeting option
