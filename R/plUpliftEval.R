EPS = 1e-20
PLOT_TYPES = c('qini', 'aqini', 'cgains', 'cuplift', 'uplift', 'balance')


#' A non-S3 "constructor" function that returns a list representing a pylift uplift
#' eval object, PlUpliftEval.  This object contains metrics and can be used to generate plots.
#'
#' @param treatment numeric vector of treatment identifiers
#' @param outcome numeric vector of outcomes
#' @param prediction numeric vector of uplift predictions
#' @param p optional "infer", numeric, numeric vector representing treatment
#'   propensities
#' @param n_bins integer number of bins on x-axis; default 20
#' @return a list representing a pylift uplift eval object
#'
#' @export
new_PlUpliftEval <- function(treatment = integer(),
                             outcome = integer(),
                             prediction = numeric(),
                             p = "infer",
                             n_bins = 20){
  # Counts.
  counts = get_counts(treatment, outcome, p)
  tc_counts = get_tc_counts(counts$Nt1o1, counts$Nt0o1, counts$Nt1o0, counts$Nt0o0)
  counts = merge(counts, tc_counts)

  # Calculate maximal curves.
  count_functions = c('get_overfit_counts', 'get_no_sure_thing_counts', 'get_no_sleeping_dog_counts')
  count_names = c('max', 'pmax', 'nosdmax')
  curve_functions = c('maximal_qini_curve', 'maximal_uplift_curve', 'maximal_cuplift_curve')
  curve_names = c('qini', 'uplift', 'cuplift')

  curve_count_xys <- list()

  for (ico in seq(1:length(count_functions))){
    for(icu in seq(1:length(curve_functions))){
      cuf <- match.fun(curve_functions[icu])
      cof <- match.fun(count_functions[ico])
      xy <- cuf(cof, counts$Nt1o1, counts$Nt0o1, counts$Nt1o0, counts$Nt0o0)
      curve_count_xys[[paste0(curve_names[icu],'_',count_names[ico])]] <- xy
    }
  }

  # Calculate Q, q1, q2 scores and other metrics.
  scores = get_scores(treatment, outcome, prediction, p)

  list(
    treatment = treatment,
    outcome = outcome,
    prediction = prediction,
    counts = counts,
    p = p,
    n_bins = n_bins,
    curve_count_xys = curve_count_xys,
    scores = scores
  )
}


#' A helper for the new_PlUpliftEval function that validates the treatment,
#' outcome, prediction, p, and n_bins arguments.
#'
#' @param treatment numeric vector of treatment identifiers
#' @param outcome numeric vector of outcomes
#' @param prediction numeric vector of uplift predictions
#' @param p optional "infer", numeric, numeric vector representing treatment
#'   propensities
#' @param n_bins integer number of bins on x-axis; default 20
#' @return a list representing a pylift uplift eval object
#'
#' @examples
#'
#' set.seed(0)
#' rl <- function(x){
#'   round(1/(1+exp(-x)))
#' }
#' n <- 2000; p <- 3
#' beta <- -0.5
#' X <- matrix(rnorm(n*p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y <- rl(pmax(beta+X[,1], 0) * W + X[,2])
#' p1 <- 1/(1+exp(-(beta+X[,1])))
#' plUpliftEval(W, Y, p1)
#'
#' \donttest{
#' library(grf)
#' set.seed(123)
#'
#' rl <- function(x){
#'   round(1/(1+exp(-x)))
#' }
#' n <- 2000; p <- 10
#' X <- matrix(rnorm(n*p), n, p)
#' W <- rbinom(n, 1, 0.2)
#' Y <- rl(rl(X[,1]) * W - rl(X[,3]) * W + rnorm(n))
#' tau.forest <- causal_forest(X, Y, W)
#' tau.hat <- predict(tau.forest, X)
#' plue <- plUpliftEval(W, Y, tau.hat$predictions)
#' plue
#' }
#'
#' @export
plUpliftEval <- function(treatment, outcome, prediction, p = "infer", n_bins = 20){

  stopifnot(length(treatment)==length(outcome) & length(treatment)==length(prediction))
  stopifnot(is.numeric(treatment))
  stopifnot(is.numeric(outcome))
  stopifnot(is.numeric(prediction))
  stopifnot(n_bins%%1==0)

  # Deal with `p`, in case float or None.
  if(is.character(p)){
    if(p=="infer"){
      p <- rep(1,length(prediction))*length(treatment[treatment==1])/length(treatment)
    }
  } else if (length(p)==1 & is.numeric(p)){
    p <- rep(1,length(prediction))*p
  } else{
    stopifnot(length(p)==length(prediction) & is.numeric(p))
  }

  new_PlUpliftEval(treatment, outcome, prediction, p)

}


pl_calc <- function(self, plot_type, n_bins=20){

  # Create bins.
  bin_range = seq(0, length(self$treatment), length.out=n_bins+1)

  # Sort `self.prediction`, descending, then get the indices in the test set that
  # these correspond to.
  prob_index = order(self$prediction,decreasing=TRUE)

  # Define whether the curve uses all data up to the percentile, or the data within that percentile.
  noncumulative_subset_func <- function(i){
    seq(1:length(self$treatment)) %in% prob_index[bin_range[i]:bin_range[i+1]]
  }
  cumulative_subset_func <- function(i){
    seq(1:length(self$treatment)) %in% prob_index[1:bin_range[i+1]]
  }

  subsetting_functions = list(
    'qini' = cumulative_subset_func,
    'aqini' = cumulative_subset_func,
    'cgains' = cumulative_subset_func,
    'cuplift' = cumulative_subset_func,
    'balance' = noncumulative_subset_func,
    'uplift' = noncumulative_subset_func
  )

  # Define the function that is calculated within the above bins.
  y_calculating_functions = list (
    'qini' = function(nt1o1, nt0o1, nt1, nt0){nt1o1/self$counts$N_treat - nt0o1/self$counts$N_contr},
    'aqini' = function(nt1o1, nt0o1, nt1, nt0){nt1o1/self$counts$N_treat - nt0o1*nt1/(nt0*self$counts$N_treat + EPS)},
    'cgains' = function(nt1o1, nt0o1, nt1, nt0){(nt1o1/(nt1+EPS)- nt0o1/(nt0+EPS))*(nt1+nt0)/self$counts$N},
    'cuplift' = function(nt1o1, nt0o1, nt1, nt0){nt1o1/(nt1+EPS) - nt0o1/(nt0+EPS)},
    'uplift' = function(nt1o1, nt0o1, nt1, nt0){ nt1o1/(nt1+EPS) - nt0o1/(nt0+EPS)},
    'balance' = function(nt1o1, nt0o1, nt1, nt0){nt1/(nt0+nt1+EPS)}
  )

  # Initialize output lists.
  x = c()
  y = c()

  # Calculate qini curve points for each bin.
  for (i in 1:n_bins){
    current_subset = subsetting_functions[[plot_type]](i)
    # Get the values of outcome in this subset for test and control.
    treated_subset = (self$treatment==1) & current_subset
    resp_treated = self$outcome[treated_subset]
    untreated_subset = (self$treatment==0) & current_subset
    resp_untreated = self$outcome[untreated_subset]
    # Get the policy for each of these as well.
    p_treated = self$p[treated_subset]
    p_untreated = self$p[untreated_subset]

    # Count the number of correct values (i.e. y==1) within each of these
    # sections as a fraction of total ads shown.
    nt1o1 = sum(resp_treated*0.5/p_treated)
    nt0o1 = sum(resp_untreated*0.5/(1-p_untreated))
    nt1 = sum(0.5/p_treated)
    nt0 = sum(0.5/(1-p_untreated))
    y = c(y, y_calculating_functions[[plot_type]](nt1o1, nt0o1, nt1, nt0))
    x = c(x,nt1+nt0)
  }

  # For non-cumulative functions, we need to do a cumulative sum of the x
  # values, because the sums in the loop only captured the counts within
  # the non-cumulative bins.
  if (plot_type %in% c('balance', 'uplift')){
    x = cumsum(x)
  }

  # Rescale x so it's between 0 and 1.
  percentile = x/max(x)

  if (! plot_type %in% c('balance', 'uplift', 'cuplift')){
    percentile = c(0,percentile)
    y = c(0,y)
  }

  list(
    x=percentile,
    y=y
  )
}

#' A port of pylift's plot function (https://github.com/wayfair/pylift) as of commit:
#' https://github.com/wayfair/pylift/tree/bb69692388b1fe085001c3ba7edf6dd81d888353
#
#' pylift: Plots the different kinds of percentage-targeted curves.
#'
#' @param plue the result of a call to the plUpliftEval constructor
#' @param plot_type string, optional Either 'qini', 'aqini', 'uplift', 'cuplift', or 'balance'.
#'     'aqini' refers to an adjusted qini plot, 'cuplift' gives a
#'     cumulative uplift plot. 'balance' gives the test-control balance
#'     for each of the bins. All others are self-explanatory.
#' @param n_bins integer, number of population bins; default 20
#' @param show_theoretical_max boolean, optional
#'     Toggle theoretical maximal qini curve, if overfitting to
#'     treatment/control. Only works for Qini-style curves.
#' @param show_practical_max boolean, optional
#'     Toggle theoretical maximal qini curve, if not overfitting to
#'     treatment/control. Only works for Qini-style curves.
#' @param show_no_dogs boolean, optional
#'     Toggle theoretical maximal qini curve, if you believe there are no
#'     sleeping dogs. Only works for Qini-style curves.
#' @param show_random_selection boolean, optional
#'     Toggle straight line indicating a random ordering. Only works for
#'     Qini-style curves.
#' @param ... additional arguments
#' @return a pylift plot
#'
#' @examples
#'
#' set.seed(0)
#' rl <- function(x){
#'   round(1/(1+exp(-x)))
#' }
#' n <- 2000; p <- 3
#' beta <- -0.5
#' X <- matrix(rnorm(n*p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y <- rl(pmax(beta+X[,1], 0) * W + X[,2])
#' p1 <- 1/(1+exp(-(beta+X[,1])))
#' plue <- plUpliftEval(W, Y, p1)
#' pl_plot(plue,
#'         show_practical_max = TRUE,
#'         show_theoretical_max = TRUE,
#'         show_no_dogs = TRUE,
#'         n_bins=20)
#'
#' \donttest{
#' library(grf)
#' set.seed(123)
#'
#' rl <- function(x){
#'   round(1/(1+exp(-x)))
#' }
#' n <- 2000; p <- 10
#' X <- matrix(rnorm(n*p), n, p)
#' W <- rbinom(n, 1, 0.2)
#' Y <- rl(rl(X[,1]) * W - rl(X[,3]) * W + rnorm(n))
#' tau.forest <- causal_forest(X, Y, W)
#' tau.hat <- predict(tau.forest, X)
#' plue <- plUpliftEval(W, Y, tau.hat$predictions)
#' plue
#' pl_plot(plue,
#'         show_practical_max = TRUE,
#'         show_theoretical_max = TRUE,
#'         show_no_dogs = TRUE,
#'         n_bins=20)
#' }
#'
#' @export

pl_plot <- function(plue,
                    plot_type='cgains',
                    n_bins=20,
                    show_theoretical_max=FALSE,
                    show_practical_max=FALSE,
                    show_random_selection=TRUE,
                    show_no_dogs=FALSE,
                    ...){
  print(paste0("plotting: ", plot_type))

  stopifnot(plot_type %in% PLOT_TYPES)

  # Calculate curve (qini, uplift, cumulative uplift).
  xy <- pl_calc(plue, plot_type=plot_type, n_bins=n_bins)

  titles <- list(
    qini = 'Qini curve',
    aqini = 'Adjusted Qini curve',
    cgains = 'Cumulative gain chart',
    cuplift = 'Cumulative uplift curve',
    uplift = 'Uplift curve',
    balance = 'Treatment balance curve'
  )

  ylabels <- list(
    qini = 'Uplift gain',
    aqini = 'Uplift gain',
    cgains = 'Uplift gain',
    cuplift = 'Cumulative lift',
    uplift = 'Lift',
    balance = 'Treatment size / (treatment size + control size)'
  )

  if ((plot_type=='aqini') | (plot_type=='cgains') | (plot_type=='qini')){
    max_plot_type='qini'
  }
  else {
    max_plot_type=plot_type
  }

  plot <- ggplot()+
    geom_line(data=as.data.frame(xy),aes(x,y)) +
    xlab("Fraction of data") +
    ylab(ylabels[[plot_type]]) +
    labs(
      title=titles[[plot_type]]
    ) +
    theme(plot.title = element_text(size=8))

  if(show_random_selection){
    plot <- plot + geom_line(data=data.frame(x=c(0,1),y=c(0,xy$y[length(xy$y)])),aes(x,y), color="grey", linetype="dotted")
  }

  if (plot_type!='balance'){
    if (show_theoretical_max){
      xys <- plue$curve_count_xys[[paste0(max_plot_type,'_max')]]
      plot <- plot +
        geom_line(data=data.frame(x=xys$x,y=xys$y),aes(x,y), color="red", linetype="dotted")
    }
    if (show_practical_max){
      xys <- plue$curve_count_xys[[paste0(max_plot_type,'_pmax')]]
      plot <- plot +
        geom_line(data=data.frame(x=xys$x,y=xys$y),aes(x,y), color="blue", linetype="dotted")
    }
    if (show_no_dogs){
      xys <- plue$curve_count_xys[[paste0(max_plot_type,'_nosdmax')]]
      plot <- plot +
        geom_line(data=data.frame(x=xys$x,y=xys$y),aes(x,y), color="green", linetype="dotted")
    }
  }

  plot

}
