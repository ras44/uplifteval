get_counts <- function(treatment,outcome,p){
  Nt1o1 = 0.5*sum(1/p[(treatment == 1) & (outcome > 0)])
  Nt0o1 = 0.5*sum(1/(1-p[(treatment == 0) & (outcome > 0)]))
  Nt1o0 = 0.5*sum(1/p[(treatment == 1) & (outcome == 0)])
  Nt0o0 = 0.5*sum(1/(1-p[(treatment == 0) & (outcome == 0)]))

  list(
    Nt1o1=Nt1o1,
    Nt0o1=Nt0o1,
    Nt1o0=Nt1o0,
    Nt0o0=Nt0o0
  )
}

get_tc_counts <- function(Nt1o1, Nt0o1, Nt1o0, Nt0o0){
  Nt1 = Nt1o0 + Nt1o1
  Nt0 = Nt0o0 + Nt0o1
  N = Nt0 + Nt1
  list(
    Nt1=Nt1,
    Nt0=Nt0,
    N=N
  )
}

get_no_sure_thing_counts <- function(Nt1o1, Nt0o1, Nt1o0, Nt0o0, func=sum){
  counts = get_tc_counts(Nt1o1, Nt0o1, Nt1o0, Nt0o0)
  Nt1=counts$Nt1
  Nt0=counts$Nt0
  N=counts$N
  sure_things = c(0, 0)
  persuadables = c(Nt1o1*Nt0/Nt1, Nt1o1)  # Add the persuadables in control.
  dogs = c(Nt0o1, Nt0o1*Nt1/Nt0)
  lost_causes = c(Nt0-dogs[0]-sure_things[0]-persuadables[0],
                 Nt1-dogs[1]-sure_things[1]-persuadables[1])
  list(
    persuadables = func(persuadables),
    dogs = func(dogs),
    sure_things = func(sure_things),
    lost_causes = func(lost_causes))
}

get_no_sleeping_dog_counts <- function(Nt1o1, Nt0o1, Nt1o0, Nt0o0, func=sum){
  counts = get_tc_counts(Nt1o1, Nt0o1, Nt1o0, Nt0o0)
  Nt1 = counts$Nt1
  Nt0 = counts$Nt0
  N = counts$N
  dogs = c(0,0)
  sure_things = c(Nt0o1, Nt0o1*Nt1/Nt0)
  lost_causes = c(Nt1o0*Nt0/Nt1, Nt1o0)
  persuadables = c(Nt0-dogs[0]-sure_things[0]-lost_causes[0],
                  Nt1-dogs[1]-sure_things[1]-lost_causes[1])
  list(
    persuadables = func(persuadables),
    dogs = func(dogs),
    sure_things = func(sure_things),
    lost_causes = func(lost_causes))
}

get_overfit_counts <- function(Nt1o1, Nt0o1, Nt1o0, Nt0o0, func=sum){
  persuadables = c(0, Nt1o1)
  dogs = c(Nt0o1, 0)
  sure_things = c(0,0)
  lost_causes = c(Nt0o0, Nt1o0)
  list(
    persuadables = func(persuadables),
    dogs = func(dogs),
    sure_things = func(sure_things),
    lost_causes = func(lost_causes))
}


maximal_qini_curve <- function(func, Nt1o1, Nt0o1, Nt1o0, Nt0o0){

  counts = get_tc_counts(Nt1o1, Nt0o1, Nt1o0, Nt0o0)
  Nt1 = counts$Nt1
  Nt0 = counts$Nt0
  N = counts$N
  pdsl = func(Nt1o1, Nt0o1, Nt1o0, Nt0o0)
  persuadables = pdsl$persuadables
  dogs = pdsl$dogs
  sure_things = pdsl$sure_things
  lost_causes = pdsl$lost_causes

  # For the overfit case, this is simply the entire treated group.
  if (identical(func, match.fun('get_overfit_counts'))){
    slope = 2
  } else {
    slope = 1
  }
  x = c(0, persuadables/N, 1-dogs/N, 1)
  # Deal with edge case where number of persuadables is greater than sleeping
  # dogs (common if this is not a treatment/control experiment, but an
  # experiment between two treatments).
  if (x[2]>x[3]){
    new_val = (x[2]+x[3])/2
    x[2] = new_val
    x[3] = new_val
  }
  y = c(0, x[2]*slope, x[2]*slope, (Nt1o1/Nt1-Nt0o1/Nt0))

  list(x=x,y=y)
}


maximal_uplift_curve <- function(func, Nt1o1, Nt0o1, Nt1o0, Nt0o0){
  counts = get_tc_counts(Nt1o1, Nt0o1, Nt1o0, Nt0o0)
  Nt1 = counts$Nt1
  Nt0 = counts$Nt0
  N = counts$N
  pdsl = func(Nt1o1, Nt0o1, Nt1o0, Nt0o0)
  persuadables = pdsl$persuadables
  dogs = pdsl$dogs
  sure_things = pdsl$sure_things
  lost_causes = pdsl$lost_causes

  y = c(1, 1, 0, 0, -1, -1)
  x = c(0, persuadables/N, persuadables/N, 1-dogs/N, 1-dogs/N, 1)
  list(x=x, y=y)
}

maximal_cuplift_curve <- function(func, Nt1o1, Nt0o1, Nt1o0, Nt0o0){

  counts = get_tc_counts(Nt1o1, Nt0o1, Nt1o0, Nt0o0)
  Nt1 = counts$Nt1
  Nt0 = counts$Nt0
  N = counts$N
  pdsl = func(Nt1o1, Nt0o1, Nt1o0, Nt0o0)
  persuadables = pdsl$persuadables
  dogs = pdsl$dogs
  sure_things = pdsl$sure_things
  lost_causes = pdsl$lost_causes

  kink_y = (persuadables[2] + sure_things[2])/(persuadables[2] + sure_things[2] + lost_causes[2])
  - (sure_things[1])/(persuadables[1] + sure_things[1] + lost_causes[1])
  y = c(1, 1, kink_y, (Nt1o1/Nt1-Nt0o1/Nt0))
  x = c(0, sum(persuadables)/N, 1-sum(dogs)/N, 1)
  list(x=x, y=y)
}

get_scores <- function(treatment, outcome, prediction, p, scoring_range=c(0,1), plot_type='all'){

  counts = get_counts(treatment, outcome, p)
  Nt1o1 = counts$Nt1o1
  Nt0o1 = counts$Nt0o1
  Nt1o0 = counts$Nt1o0
  Nt0o0 = counts$Nt0o0
  counts = get_tc_counts(Nt1o1, Nt0o1, Nt1o0, Nt0o0)
  Nt1 = counts$Nt1
  Nt0 = counts$Nt0
  N = counts$N

  riemann <- function(xy){
    x <- xy$x
    y <- xy$y
    stopifnot(length(x)==length(y))
    xi <- utils::head(x, length(x)-1)
    xip1 <- utils::tail(x,length(x)-1)
    yi <- utils::head(y,length(y)-1)
    yip1 <- utils::tail(y, length(y)-1)
    avgy <- (yi+yip1)/2
    dx <- (xip1-xi)
    sum(avgy*dx)
  }

  qini_riemann = riemann(maximal_qini_curve(get_overfit_counts, Nt1o1, Nt0o1, Nt1o0, Nt0o0))
  practical_qini_riemann = riemann(maximal_qini_curve(get_no_sure_thing_counts, Nt1o1, Nt0o1, Nt1o0, Nt0o0))

  overall_lift = (Nt1o1/Nt1-Nt0o1/Nt0)
  qini_max = qini_riemann - 0.5*overall_lift
  practical_qini_max = practical_qini_riemann - 0.5*overall_lift

  # The predicted Qini curve.
  # First we need to reorder the y values and y_pred based on this reordering
  # We calculate TOT roughly here so we have a way of distinguishing those that (ordered, treated) and those that (ordered, untreated).
  y = (2*treatment - 1)*outcome


  sortbyprediction <- function(vec){
    vec[order(prediction,decreasing = TRUE)]
  }

  y_ordered = sortbyprediction(y)
  tr_ordered = sortbyprediction(treatment)
  p_ordered = sortbyprediction(p)

  auc <- function(method='qini'){
    # Calculate the area.
    uplift_last = 0
    nt1o1 = 0
    nt0o1 = 0
    nt1 = EPS
    nt0 = EPS
    pred_riemann = 0
    uplifts = c()

    for(i in seq(round(scoring_range[1]*length(treatment))+1, round(scoring_range[2]*length(treatment)))){

      if(y_ordered[i] > 0){
        nt1o1 = nt1o1 + 0.5*(1/p_ordered[i])
      } else if (y_ordered[i] < 0){
          nt0o1 = nt0o1 + 0.5*(1/(1-p_ordered[i]))
      }

      if (tr_ordered[i] == 1){
        nt1 = nt1 + 0.5*(1/p_ordered[i])
      }
      else{
        nt0 = nt0 + 0.5*(1/(1-p_ordered[i]))
      }

      if(method=='qini'){
        uplift_next = nt1o1/Nt1-nt0o1/Nt0
      } else if (method=='cgains'){
        uplift_next = (nt1o1/nt1-nt0o1/nt0)*(nt1+nt0)/N
      } else if (method=='aqini'){
        uplift_next = nt1o1/Nt1-nt0o1*nt1/(nt0*Nt1 + EPS)
      }

      uplifts = c(uplifts, uplift_next)
      # each point corresponds to an x delta of 1/N
      pred_riemann = pred_riemann + 1/2*(uplift_next+uplift_last)/N
      uplift_last = uplift_next
    }

    AUC = pred_riemann - 0.5*overall_lift*(scoring_range[2]**2 - scoring_range[1]**2)
    maxgain = max(uplifts)
    list(
      AUC = AUC,
      maxgain = maxgain)
  }

  # Dictionary to store all scores.
  scores = list()
  # Raw max scores.
  scores$Q_max = qini_max
  scores$overall_lift = overall_lift
  scores$Q_practical_max = practical_qini_max
  if ((plot_type=='qini') | (plot_type=='all')){
    # Qini curve scores.
    auc = auc(method='qini')
    scores$Q_qini = auc$AUC
    scores$max_qini = auc$maxgain
    scores$q1_qini = scores$Q_qini/scores$Q_max
    scores$q2_qini = scores$Q_qini/scores$Q_practical_max
  } else if ((plot_type=='cgains') | (plot_type=='all')){
    # Scores for cumulative gains curve.
    auc = auc(method='cgains')
    scores$Q_cgains = auc$AUC
    scores$max_cgains = auc$maxgain
    scores$q1_cgains = scores$Q_cgains/scores$Q_max
    scores$q2_cgains = scores$Q_cgains/scores$Q_practical_max
  } else if ((plot_type=='aqini') | (plot_type=='all')){
    # Scores for adjusted qini curve.
    auc = auc(method='aqini')
    scores$Q_aqini = auc$AUC
    scores$max_aqini = auc$maxgain
    scores$q1_aqini = scores$Q_aqini/scores$Q_max
    scores$q2_aqini = scores$Q_aqini/scores$Q_practical_max
  }
  scores
}

