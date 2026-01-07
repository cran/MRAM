#' Variable Selection via the Multivariate Regression Association Measure
#'
#' @param y_data A \eqn{n \times d} matrix of responses, where \eqn{n} is the sample size.
#' @param x_data A \eqn{n \times p} matrix of predictors.
#' @description Select a subset of \eqn{\bf X} which can be used to predict \eqn{\bf Y} based on \eqn{T_n}.
#' @details \code{vs_mram} performs forward stepwise variable selection based on the multivariate regression association measure proposed in Shih and Chen (2025). At each step, it selects the predictor with the highest conditional predictability for the response given the previously selected predictors. The algorithm is modified from the FOCI algorithm from Azadkia and Chatterjee (2021).
#'
#' @return The vector containing the indices of the selected predictors in the order they were chosen.
#'
#' @references Azadkia and Chatterjee (2021) A simple measure of conditional dependence, Annals of Statistics, 46(6): 3070-3102.
#' @references Shih and Chen (2026) Measuring multivariate regression association via spatial sign, Computational Statistics & Data Analysis, 215, 108288.
#' @seealso \code{\link{mram}}
#'
#' @export
#'
#' @examples
#' library(MRAM)
#'
#' n = 200
#' p = 10
#'
#' set.seed(1)
#' x_data = matrix(rnorm(p*n),n,p)
#' colnames(x_data) = paste0(rep("X",p),seq(1,p))
#'
#' y_data = x_data[,1]*x_data[,2]+x_data[,1]-x_data[,3]+rnorm(n)
#' colnames(x_data)[vs_mram(y_data,x_data)] # selected variables
#'
#' \dontrun{
#'
#' n = 500
#' p = 10
#'
#' set.seed(1)
#' x_data = matrix(rnorm(p*n),n,p)
#' colnames(x_data) = paste0(rep("X",p),seq(1,p))
#'
#' # Linear
#' y_data = matrix(0,n,2)
#' y_data[,1] = x_data[,1]*x_data[,2]+x_data[,1]-x_data[,3]+rnorm(n)
#' y_data[,2] = x_data[,2]*x_data[,4]+x_data[,2]-x_data[,5]+rnorm(n)
#' colnames(x_data)[vs_mram(y_data,x_data)] # selected variables
#'
#' # Nonlinear
#' y_data = matrix(0,n,2)
#' y_data[,1] = x_data[,1]*x_data[,2]+sin(x_data[,1]*x_data[,3])+0.3*rnorm(n)
#' y_data[,2] = cos(x_data[,2]*x_data[,4])+x_data[,3]-x_data[,4]+0.3*rnorm(n)
#' colnames(x_data)[vs_mram(y_data,x_data)] # selected variables
#'
#' # Non-additive error
#' y_data = matrix(0,n,2)
#' y_data[,1] = abs(x_data[,1]+runif(n))^(sin(x_data[,2])-cos(x_data[,3]))
#' y_data[,2] = abs(x_data[,2]-runif(n))^(sin(x_data[,3])-cos(x_data[,4]))
#' colnames(x_data)[vs_mram(y_data,x_data)] # selected variables
#' }

vs_mram = function(y_data,
                   x_data) {

  y_data = as.matrix(y_data)
  x_data = as.matrix(x_data)

  p = dim(x_data)[2]
  p_seq = c(1:p)

  var_select = 0
  T_temp = 0

  i = 1
  repeat {

    T_step_j = numeric(p)
    for (k in p_seq) {

      T_step_j[k] = mram(y_data,x_data[,c(k,var_select)])$T_est

    }

    T_cond = (T_step_j-T_temp)/(1-T_temp)

    if (max(T_cond) < 0) {break}

    var_select[i] = which.max(T_cond)
    p_seq = setdiff(p_seq,var_select)
    T_temp = max(T_step_j)


    i = i+1

  }

  return(var_select)

}





