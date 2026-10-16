#' Wrapper function for the implementation of lagged WQS.
#'
#' @param data Data frame containing observations in long format.
#' @param timevar Enquoted variable name identifying the repeated measure / time variable.
#' @param wqs_parms A list containing parameters to be passed to the WQS algorithm. See gWQS package for details.
#' @param outcome An enquoted variable name identifying the outcome measure.
#' @param ID An enquoted variable name identifying the subject identifier.
#' @param rDLM_parms (optional). A list containing parameters to be passed to the GAM algorithm. See gamm4 package for details. Parameters wqs, time, by, and id are created by the lwqs function and passed to the gamm4 function automatically.
#' @param ci_level (optional). Confidence level for the time-varying effect band. Defaults to 0.95.
#'
#' @return The \code{lwqs} function returns a list containing final model output and time-specific model parameters.
#'
#' \item{parameters}{This list contains several objects summarizing different stages of the lagged ensemble model. The first
#' object, \code{res}, contains output from the gWQS algorithm applied to each discrete repeated measure in the
#' overall model; see package gWQS for details. The second output, \code{wqstime}, provides the mixture index,
#' identified as "wqs", estimated for each subject at each discrete time point. The third item, \code{weightstime}, provides
#' the weights estimated for each predictor at each discrete time point.}
#' \item{plot}{A combined plot (as a grob) summarizing the output of the lwqs algorithm.}
#'
#' @examples
#' \dontrun{ # identify predictor variables used in mixture
#' mixvars=names(lwqs_data)[5:9]
#'
#' # run model.
#' model=lwqs(data=lwqs_data,
#'            timevar="time",
#'            wqs_parms=list(formula=out ~ wqs,
#'               data = lwqs_data,
#'               mix_name=mixvars,
#'               b_constr = TRUE,
#'               b1_pos=TRUE,
#'               b = 5,
#'               q = 5,
#'               validation = 0,
#'               family = "gaussian",
#'               seed = 1),
#'               outcome="out",
#'               ID="ID")}
#'
#' @import ggplot2
#' @import gWQS
#' @import plyr
#' @import gamm4
#' @import gridExtra
#' @import data.table
#' @importFrom stats predict qnorm
#'
#' @export lwqs

lwqs <- function(data,
                 timevar,
                 wqs_parms,
                 outcome,
                 ID,
                 rDLM_parms = list(formula = wqs ~ s(time, by = y, bs = "cr"),
                                   random = ~ (1 | id)),
                 ci_level = 0.95) {

  ## ---- input validation -------------------------------------------------
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  if (!is.character(timevar) || length(timevar) != 1L)
    stop("'timevar' must be a single (quoted) column name.", call. = FALSE)
  if (!is.character(outcome) || length(outcome) != 1L)
    stop("'outcome' must be a single (quoted) column name.", call. = FALSE)
  if (!is.character(ID) || length(ID) != 1L)
    stop("'ID' must be a single (quoted) column name.", call. = FALSE)

  missing_cols <- setdiff(c(timevar, outcome, ID), names(data))
  if (length(missing_cols))
    stop("Column(s) not found in 'data': ",
         paste(missing_cols, collapse = ", "), call. = FALSE)

  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      ci_level <= 0 || ci_level >= 1)
    stop("'ci_level' must be a single number strictly between 0 and 1.",
         call. = FALSE)

  ## bind names used only inside aes() to silence R CMD check NOTEs
  Beta <- lower <- upper <- time <- mean_weight <- mix_name <- y <- NULL

  ## ---- per-time-point WQS ------------------------------------------------
  wqs_parameters <- function(x) {
    parms <- wqs_parms
    parms$solve_dir_issue <- "average"
    parms$data <- x
    model <- do.call(gwqs, parms)
    ## attach subject id and outcome to the per-subject wqs index.
    ## Build the data frame by explicit columns so the result does not
    ## depend on column ordering or accidental name collisions.
    model$wqs <- data.frame(wqs = model$wqs,
                            ID  = x[[ID]],
                            y   = x[[outcome]],
                            stringsAsFactors = FALSE)
    model
  }

  res         <- plyr::dlply(data, timevar, wqs_parameters)
  wqstime     <- plyr::ldply(res, function(x) rbind(x$wqs))
  weightstime <- plyr::ldply(res, function(x) rbind(x$final_weights))

  results <- list(res = res, wqstime = wqstime, weightstime = weightstime)

  ## ---- reverse distributed lag model (gamm4) -----------------------------
  rDLM_parms$data       <- wqstime
  rDLM_parms$data$id    <- rDLM_parms$data$ID
  rDLM_parms$data$time  <- rDLM_parms$data[[timevar]]
  fit <- do.call(gamm4, rDLM_parms)

  ## Extract the time-varying smooth term and its standard error through the
  ## documented predict() interface rather than scraping plot.gam()'s
  ## invisibly-returned (and undocumented) list. This avoids both the
  ## positional-index fragility and the side-effect of drawing a plot.
  newd <- data.frame(time = sort(unique(rDLM_parms$data$time)))
  ## predict.gam needs every variable used in the smooth; 'y' enters via
  ## by=y, so hold it at 1 to read off the per-unit (by-variable) effect.
  newd$y  <- 1
  newd$id <- rDLM_parms$data$id[1]

  pr <- predict(fit$gam, newdata = newd, type = "terms", se.fit = TRUE)

  ## locate the column corresponding to the s(time, by = y) term
  term_col <- grep("time", colnames(pr$fit), fixed = TRUE)
  if (length(term_col) < 1L)
    stop("Could not locate the time-varying smooth term in the fitted GAM.",
         call. = FALSE)
  term_col <- term_col[1L]

  crit  <- stats::qnorm(1 - (1 - ci_level) / 2)
  fplot <- data.frame(time = newd$time,
                      Beta = as.numeric(pr$fit[, term_col]),
                      SE   = as.numeric(pr$se.fit[, term_col]))
  fplot$upper <- fplot$Beta + crit * fplot$SE
  fplot$lower <- fplot$Beta - crit * fplot$SE

  ## ---- plots -------------------------------------------------------------
  bplot <- ggplot(data = fplot, aes(x = time, y = Beta)) +
    geom_line(linewidth = 1.5) +
    theme_bw() +
    geom_hline(yintercept = 0, linewidth = 0.6, lty = 2) +
    xlab("Time") + ylab("Time Varying Effect of Mixture")

  bplot <- bplot +
    geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 2, alpha = 0.35)

  weightstime$time <- weightstime[[timevar]]
  wplot <- ggplot(weightstime,
                  aes(x = time, y = mean_weight,
                      group = mix_name, colour = mix_name)) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x,
                linewidth = 1.5) +
    theme_bw() + xlab("Time")

  plot <- gridExtra::arrangeGrob(bplot, wplot, nrow = 1)

  list(parameters = results, plot = plot)
}

#' Function to extract time-varying mixture (wqs) index from lWQS object.
#' @param lobj An object returned from the lwqs function.
#' @return Data frame containing the time index, wqs index estimated at each repeated measure, subject ID, and the outcome variable.
#' @examples
#' \dontrun{
#' # identify predictor variables used in mixture
#' mixvars=names(lwqs_data)[5:9]
#'
#' # run model. Note for example run-time only a small number of bootstraps
#' # (b=5) is used; set b to be >50 in practice.
#' model=lwqs(data=lwqs_data,
#'            timevar="time",
#'            wqs_parms=list(formula=out ~ wqs,
#'               data = lwqs_data,
#'               mix_name=mixvars,
#'               b_constr = TRUE,
#'               b1_pos=FALSE,
#'               b = 5,
#'               q = 5,
#'               validation = 0,
#'               family = "gaussian",
#'               seed = 1),
#'               outcome="out",
#'               ID="ID")
#'
#' # use extract_mixture to access time-varying wqs index
#' mixtime=extract_mixture(model)}
#'
#' @export extract_mixture
extract_mixture <- function(lobj) {
  data.frame(lobj$parameters$wqstime)
}

#' Function to extract time-varying weights from lWQS object.
#' @param lobj An object returned from the lwqs function.
#' @return A (long-form) data frame containing the time index and corresponding variable weights estimated in an lWQS.
#' @examples
#' \dontrun{
#' # identify predictor variables used in mixture
#' mixvars=names(lwqs_data)[5:9]
#'
#' # run model
#' model=lwqs(data=lwqs_data,
#'            timevar="time",
#'            wqs_parms=list(formula=out ~ wqs,
#'               data = lwqs_data,
#'               mix_name=mixvars,
#'               b_constr = TRUE,
#'               b1_pos=TRUE,
#'               b = 5,
#'               q = 5,
#'               validation = 0,
#'               family = "gaussian",
#'               seed = 1),
#'               outcome="out",
#'               ID="ID")
#'
#' # use extract_weights to access time-varying predictor weights
#' timeweights=extract_weights(model)}
#'
#' @export extract_weights
extract_weights <- function(lobj) {
  data.frame(lobj$parameters$weightstime)
}


#' Simulated dataset for accompanying vignette.
#' @name lwqs_data
#' @docType data
#' @return A data frame containing simulated data to explore the \code{lwqs} algorithm. Variables included are as follows:
#' \item{ID}{Variable identifying each simulated subject. Data reflect 30 successive measures per subject.}
#' \item{Sex}{A simulated binary covariate, either 1 or 0.}
#' \item{time}{Variable identifying the successive timing of each repeated measure.}
#' \item{out}{Simulated outcome on standardized scale.}
#' \item{pred1}{First simulated time-varying predictor. This has a large positive association with "out" from times 11-20.}
#' \item{pred2}{Second simulated time-varying predictor. This has a moderate positive association with "out" from times 11-20.}
#' \item{pred3}{Third simulated time-varying predictor. This has a moderate negative association with "out" from times 1-10.}
#' \item{pred4}{Fourth simulated time-varying predictor. This has a strong negative association with "out" from times 1-10.}
#' \item{pred5}{Fifth simulated time-varying predictor. This has no significant association with "out".}
#' @usage data(lwqs_data)
NULL
