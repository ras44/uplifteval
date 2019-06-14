#' Creates an uplift plot of cumulative differential treatment/control outcomes
#' versus model score.  Also provides a selection of metrics: max uplift as pct
#' of total control outcome, optimum users targeted and optimum score targeting
#' range.
#'
#' @param p1 numeric vector of uplift predictions; can also be predicted outcomes
#'   for treated case (in this case p0 should contain predicted outcomes for
#'   the control case)
#' @param W binary vector {1,0} of treatment assignments
#' @param Y numeric vector of responses
#' @param ns integer number of samples per bootstrap iteration; default min(table(W))
#' @param n_bs integer number of bootstrap iterations
#' @param W_label optional labels for the treatment options (default W)
#' @param p0 optional numeric vector of predicted outcomes for control case
#' @param balanced optional boolean whether to sample equal proportions from
#'   treatment and control cases; default TRUE
#' @param replace optional boolean whether to use replacement when sampling;
#'   default TRUE
#' @param x_interval optional numeric the interval with which to split the
#' @param ... additional arguments (unused)
#'   x-axis
#'
#' @import ggplot2 graphics
#' @importFrom dplyr %>%
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
#' plot_uplift(p1, W, Y, n_bs=20, x_interval = 0.05, balanced = TRUE)
#'
#'
#' set.seed(0)
#' n <- 2000; p <- 3
#' beta <- -0.5
#' X <- matrix(rnorm(n*p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y <- pmax(beta+X[,1], 0) * W + X[,2]
#' p1 <- 1/(1+exp(-(beta+X[,1])))
#' plot_uplift(p1, W, Y, n_bs=20, x_interval = 0.05, balanced = TRUE)
#'
#'
#' \donttest{
#' library(grf)
#' set.seed(123)
#'
#' rl <- function(x){
#'   round(1/(1+exp(-x)))
#' }
#' n = 2000; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.2)
#' Y = rl(rl(X[,1]) * W - rl(X[,3]) * W + rnorm(n))
#' tau.forest = causal_forest(X, Y, W)
#' tau.hat = predict(tau.forest, X)
#' plot_uplift(tau.hat$predictions, W, Y, n_bs=20, x_interval = 0.05, balanced = FALSE)
#' plot_uplift(tau.hat$predictions, W, Y, n_bs=20, x_interval = 0.05, balanced = TRUE)
#'
#' }
#' @export

plot_uplift <- function(p1,
                        W,
                        Y,
                        ns=min(table(W)),
                        n_bs = 1,
                        W_label=W,
                        p0=rep(0,length(p1)),
                        balanced = TRUE,
                        replace = TRUE,
                        x_interval = 0.1,
                        ...
                        ){

  q_pred_t <- c()
  q_pred_c <- c()
  q_resp <- c()
  q_treat <- c()
  q_treat_label <- c()
  q_iter <- c()

  if(balanced){
    for(i in c(1:n_bs)){
      inds_t0 <- sample(seq_along(W)[W==0],ns, replace = replace)
      inds_t1 <- sample(seq_along(W)[W==1],ns, replace = replace)
      q_pred_t <- c(q_pred_t, c(p1[inds_t0],p1[inds_t1]))
      q_pred_c <- c(q_pred_c, c(p0[inds_t0],p0[inds_t1]))
      q_resp <- c(q_resp, c(Y[inds_t0],Y[inds_t1]))
      q_treat <- c(q_treat, c(W[inds_t0],W[inds_t1]))
      q_treat_label <- c(q_treat_label, c(W_label[inds_t0],W_label[inds_t1]))
      q_iter <- c(q_iter, c(rep(i,ns),rep(i,ns)))
    }
  }
  else {
    for(i in c(1:n_bs)){
      print(paste0("boostrap iter: ", i))
      # subsample must be defined for a particular model
      # it must append to the above q_ variables in the calling environment via <<-
      inds <- sample(length(p1), ns, replace = replace)
      q_pred_t <- c(q_pred_t, p1[inds])
      q_pred_c <- c(q_pred_c, p0[inds])
      q_resp <- c(q_resp, Y[inds])
      q_treat <- c(q_treat, W[inds])
      q_treat_label <- c(q_treat_label, W_label[inds])
      q_iter <- c(q_iter, rep(i,ns))
    }
  }


  dif.pred <- q_pred_t - q_pred_c

  if(all(c(0,1) %in% unique(q_treat))){
    q_treat <- 2*q_treat-1
  }

  mm <- cbind(dif.pred = dif.pred, y = q_resp, ct = q_treat, ctl = q_treat_label, dif.pred_r = rank(-dif.pred), i = q_iter)

  stopifnot(all(c(-1,1) %in% unique(q_treat)))

  # Cumulatively sum the outcome by treatment and control groups
  mmo <- mm[order(-dif.pred)[],]
  mmo <- cbind(mmo, cdr = cumsum(mmo[,'y']*mmo[,'ct']))
  # show cumsum(dr) for individual bootstrap iterations
  mmo_df <- as.data.frame(mmo)

  mmo_df <- mmo_df %>%
    dplyr::group_by(i) %>%
    dplyr::arrange(-dif.pred) %>% dplyr::mutate(cdri = cumsum(y*ct))

  # Plot mean with min/max error bars, quantizing scoring (q = 10^x)
  q <- 1/x_interval
  mmo_dfs <- mmo_df %>% dplyr::group_by(dif.pred = floor(dif.pred*q)/q) %>%
    dplyr::summarize(min_cdri = min(cdri),
                     sd_cdri = sd(cdri),
                     max_cdri = max(cdri),
                     mean_cdri = mean(cdri))
  mmo_dfs <- cbind(mmo_dfs, group=rev(seq_along(mmo_dfs$mean_cdri))) # note: group is reversed
  mmo_dfs$sd_cdri[is.na(mmo_dfs$sd_cdri)] <- 0

  max_mean_uplift <- max(mmo_dfs$mean_cdri)
  max_mean_uplift

  max_group <- min(mmo_dfs$group[mmo_dfs$mean_cdri==max_mean_uplift])
  max_group
  dr_opt_min <- min(mmo_dfs$dif.pred[mmo_dfs$group<=max_group])
  dr_opt_min
  dr_opt_max <- max(mmo_dfs$dif.pred[mmo_dfs$group<=max_group])
  dr_opt_max

  mean_total_response_tc <- sum(mmo_df$y[mmo_df$ct==-1])/n_bs
  mean_total_response_tc

  diff_uplift_pct <- 100 * max_mean_uplift / mean_total_response_tc
  diff_uplift_pct

  mean_opt_users_targeted <- length(mmo_df$y[mmo_df$dif.pred>=dr_opt_min & mmo_df$ct==-1])/n_bs
  mean_opt_users_targeted
  mean_total_sampled_users <- length(mmo_df$y[mmo_df$ct==-1])/n_bs
  mean_total_sampled_users

  if(max_mean_uplift == 0){
    dr_opt_min <- dr_opt_max
    opt_users_targeted <- 0
  }

  config <- list(
    ct=0,
    tt=1,
    model_type="cf"
  )

  wr <- function(s,c){
    whisker::whisker.render(s,c)
  }

  s <- sign(mean_total_response_tc)

  p1 <- ggplot(mmo_dfs, aes(x = dif.pred, y = mean_cdri/mean_total_response_tc*s)) +
    geom_point(size = 2) +
    #geom_errorbar(aes(ymin = min_cdri, ymax = max_cdri)) +
    geom_errorbar(aes(ymin = (mean_cdri-sd_cdri)/mean_total_response_tc*s, ymax = (mean_cdri+sd_cdri)/mean_total_response_tc*s)) +
    scale_x_reverse() +
    xlab("score") +
    ylab("CDR") +
    labs(
      title=wr("Cumulative Differential Response T{{{tt}}}-T{{{ct}}}\n  model_type: {{{model_type}}}; {{{n_bs}}} subsample iterations of size {{{mean_total_sampled_users}}}\n  max mean uplift: {{{diff_uplift_pct}}}% of T{{{ct}}} result\n  optimum: {{{opt_pct_users_targeted}}}% of subjects scoring: [{{{dr_opt_min}}},{{{dr_opt_max}}}]; uplift/subject: {{{uplift_per_subject}}}", list(tt=config$tt, ct=config$ct, model_type=config$model_type, n_bs=n_bs, max_mean_uplift=round(max_mean_uplift,2), diff_uplift_pct=round(diff_uplift_pct,2), mean_opt_users_targeted=round(mean_opt_users_targeted,0), mean_total_sampled_users = mean_total_sampled_users, opt_pct_users_targeted = round(100*mean_opt_users_targeted/mean_total_sampled_users,3), dr_opt_min=round(dr_opt_min,3), dr_opt_max=round(dr_opt_max,3), uplift_per_subject = round(max_mean_uplift/mean_opt_users_targeted,3)))
    ) +
    theme(plot.title = element_text(size=8))

  p2 <- ggplot(data.frame(score=mmo_df$dif.pred, treatment=as.factor(mmo_df$ctl)), aes(x=score, fill=treatment)) +
    geom_histogram(position="dodge", bins = 100, linetype="blank") +
    scale_x_reverse() +
    theme(legend.position = c(0.9, 0.8))

  gridExtra::grid.arrange(p1, p2, nrow=2)

}
