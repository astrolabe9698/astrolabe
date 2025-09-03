#' Generate a Data Frame from Inline R Formulas
#'
#' Evaluates a set of R statements (one per line) inside a temporary environment
#' that contains \code{n} and the current session scope, then returns all created
#' objects (except \code{n}) as a data frame.
#'
#' @param formulas Character vector of code lines, or a single multi-line string.
#'   Lines starting with \code{"#"} are ignored.
#' @param n Integer. Sample size available inside the temporary environment (bound to \code{n}).
#'
#' @return A data frame with one column per symbol created by the code (excluding \code{n}).
#'
#' @examples
#' \dontrun{
#' code <- "
#' X <- rnorm(n)
#' Y <- X^3 + rnorm(n, 0, 0.1)
#' S <- rnorm(n)
#' W <- rnorm(n)
#' Z <- X^2 + log(abs(Y)) + W^2 + rnorm(n, 0, 0.1)
#' "
#' df <- dataframe_generation(code, n = 500)
#' str(df)
#'
#' code1 <- "
#' X <- rnorm(n)
#' Y <- X^2 + rnorm(n, 0, 0.1)
#' W <- rnorm(n)
#' Z <- log(abs(W) + 1) + rnorm(n, 0, 0.1)
#' "
#' df1 <- dataframe_generation(code1, n = 500)
#'
#' code2 <- "
#' X <- runif(n, -pi, pi)
#' Y <- sin(X) + rnorm(n, 0, 0.05)
#' Z <- cos(Y) + rnorm(n, 0, 0.05)
#' "
#' df2 <- dataframe_generation(code2, n = 500)
#'
#' code3 <- "
#' L <- rnorm(n)
#' W <- rnorm(n)
#' Z <- rnorm(n)
#' X <- L^4 + rnorm(n, 0 , 0.1)
#' Y <- W^5 + rnorm(n, 0, 0.1)
#' "
#' df3 <- dataframe_generation(code3, n = 500)
#'
#' code4 <- "
#' X <- rnorm(n)
#' Y <- X^2 + rnorm(n, 0, 0.1)
#' Z <- rnorm(n, 0, 0.1)
#' S <- exp(abs(Z)) + rnorm(n, 0, 0.1)
#' "
#' df4 <- dataframe_generation(code4, n = 500)
#'
#' code5 <- "
#' X <- rnorm(n)
#' Y <- 2*X + rnorm(n, 0, 0.1)
#' Z <- 1/(X+1)^2 + rnorm(n, 0, 0.1)
#' S <- atan(Z) + rnorm(n, 0, 0.1)
#' "
#' df5 <- dataframe_generation(code5, n = 500)
#' }
#'
#' @importFrom stats rnorm
#' @export
dataframe_generation <- function(formulas, n = 300) {
  # Temporary environment with access to all loaded R functions
  temp_env <- new.env(parent = globalenv())
  assign("n", n, envir = temp_env)

  # Support for multi-line string input
  if (is.character(formulas) && length(formulas) == 1 && grepl("\n", formulas)) {
    formulas <- unlist(strsplit(formulas, "\n"))
  }

  # Evaluate each line
  for (line in formulas) {
    line <- trimws(line)
    if (nzchar(line) && !grepl("^#", line)) {
      tryCatch(
        eval(parse(text = line), envir = temp_env),
        error = function(e) {
          warning(sprintf("Error in line: '%s'\nMessage: %s", line, e$message))
        }
      )
    }
  }

  # Build dataframe
  vars <- setdiff(ls(temp_env), "n")
  df <- as.data.frame(mget(vars, envir = temp_env), stringsAsFactors = FALSE)
  return(df)
}

#' Simulate Multivariate Short Panel Time Series from Symbolic Formulas
#'
#' Generates a panel dataset with \code{n} independent series over \code{T_points}
#' time points. Each variable is defined by a symbolic formula \code{lhs ~ rhs}
#' where \code{rhs} can depend on current values, a per-time deterministic
#' \code{trend}, Gaussian noise term \code{error}, and individual lags via
#' \code{lag(var)}.
#'
#' @param n Integer. Number of individuals (panel units).
#' @param T_points Integer. Number of time points per individual.
#' @param formulas List of parsed formulas (e.g., created with \code{as.formula}),
#'   each with a single \code{lhs ~ rhs}.
#' @param trend_fun_list Named list of functions; for each \code{lhs} name, a function
#'   of the time index \code{i} returning a deterministic trend contribution.
#' @param sd Numeric. Standard deviation of the Gaussian noise \code{error}.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return A data frame with columns \code{id}, \code{t}, and one column per variable
#'   defined in \code{formulas}. Initial state at \code{t=1} is standard normal.
#'
#' @details
#' Within each formula's environment at time \code{i}:
#' \itemize{
#'   \item \code{trend} is provided from \code{trend_fun_list[[lhs]](i)} if present, else 0;
#'   \item \code{error} is \code{rnorm(1, 0, sd)};
#'   \item \code{lag(var)} returns the previous-time value of \code{var} for the same \code{id}.
#' }
#'
#' @examples
#' \dontrun{
#' fmls <- list(
#'   as.formula(X ~ 0.7 * lag(X) + 0.2 * Y + trend + error),
#'   as.formula(Y ~ 0.5 * lag(Y) + error)
#' )
#' trends <- list(X = function(i) 0.05 * i)
#' panel <- generate_multivariate_time_series(n = 200, T_points = 10,
#'                                            formulas = fmls,
#'                                            trend_fun_list = trends,
#'                                            sd = 0.1, seed = 123)
#' head(panel)
#' }
#'
#' @importFrom stats rnorm
#' @export
generate_multivariate_time_series <- function(
    n = 1000,
    T_points = 4,
    formulas = list(),
    trend_fun_list = list(),
    sd = 0.1,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  var_names <- sapply(formulas, function(f) as.character(f[[2]]))
  id <- rep(1:n, each = T_points)
  t <- rep(1:T_points, times = n)
  data <- data.frame(id = id, t = t)

  # Initialize columns
  for (v in var_names) {
    data[[v]] <- NA
  }

  # Initial values at t = 1
  for (v in var_names) {
    data[data$t == 1, v] <- rnorm(n, 0, 1)
  }

  # Fill time series
  for (j in 1:n) {
    idx <- which(data$id == j)
    for (i in 2:T_points) {
      row <- idx[i]
      for (f in formulas) {
        lhs <- as.character(f[[2]])
        rhs <- f[[3]]

        env <- new.env(parent = baseenv())
        env$trend <- if (lhs %in% names(trend_fun_list)) {
          trend_fun_list[[lhs]](i)
        } else 0
        env$error <- rnorm(1, 0, sd)

        # Provide access to all variable values
        for (v in var_names) env[[v]] <- data[row, v]

        # Simple lag function per individual
        env$lag <- function(var) {
          var_name <- deparse(substitute(var))
          lags <- ave(data[[var_name]], data$id, FUN = function(x) c(NA, head(x, -1)))
          return(lags[row])
        }

        data[row, lhs] <- eval(rhs, env)
      }
    }
  }

  return(data)
}
