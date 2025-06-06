#' Variable Selection via the Multivariate Regression Association Measure
#'
#' @param y_data A \eqn{n \times d} matrix of responses.
#' @param x_data A \eqn{n \times p} matrix of predictors.
#' @description Perform variable selection via the multivariate regression association measure proposed in Shih and Chen (2025).
#' @details \code{vs_mram} is a forward and stepwise variable selection algorithm which utilizes the multivariate regression association measure proposed in Shih and Chen (2025). The Algorithm is modified from the feature ordering by conditional independence (FOCI) algorithm from Azadkia and Chatterjee (2021).
#'
#' @return The vector containing the indices of the selected predictors in the order they were chosen.
#'
#' @references Azadkia and Chatterjee (2021) A simple measure of conditional dependence, Annals of Statistics, 46(6): 3070-3102.
#' @references Shih and Chen (2025) Measuring multivariate regression association via spatial sign (in revision, Computational Statistics & Data Analysis)
#' @seealso \code{\link{mram}}
#'
#' @export
#'
#' @examples
#' n = 200
#' p = 10
#'
#' x_data = matrix(rnorm(p*n),n,p)
#' y_data = x_data[,1]*x_data[,2]+x_data[,1]-x_data[,3]+rnorm(n)
#' colnames(x_data) = paste0(rep("X",p),seq(1,p))
#'
#' library(MRAM)
#' mram_res = vs_mram(y_data,x_data)

vs_mram = function(y_data,
                   x_data) {

  if (is.matrix(y_data)) {

    n = dim(y_data)[1]

  } else {

    n = length(y_data)

  }

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





