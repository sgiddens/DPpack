# This file contains code derived from a Python project originally licensed
# under the Apache License, Version 2.0; you may not use this file except in
# compliance with the License. You may obtain a copy of the License at:
# http://www.apache.org/licenses/LICENSE-2.0.
#
# The original Python project can be found here:
# https://github.com/BorjaBalle/analytic-gaussian-mechanism/blob/master/agm-example.py.
# The code for this is based on original work on the analytic Gaussian mechanism
# described here: https://proceedings.mlr.press/v80/balle18a.
# Code found in this file is essentially the same as the original Python
# project, but has been translated to R.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Calibrate Analytic Gaussian Mechanism
#'
#' Calibrate a Gaussian perturbation for differential privacy using the analytic
#' Gaussian mechanism \insertCite{Balle2018}{DPpack}.
#'
#' @param epsilon Positive real number defining the epsilon privacy parameter.
#' @param delta Positive real number defining the delta privacy parameter.
#' @param sensitivity Real number corresponding to the l2-global sensitivity.
#' @param tol Error tolerance for binary search.
#' @return Standard deviation of Gaussian noise needed to achieve
#'   \code{(epsilon, delta)}-DP for given global sensitivity.
#'
#' @references \insertRef{Balle2018}{DPpack}
#'
#' @keywords internal
#'
#' @export
calibrateAnalyticGaussianMechanism <- function(epsilon, delta,
                                               sensitivity, tol=1e-12){
  caseA <- function(epsilon, s){
    stats::pnorm(sqrt(epsilon*s)) - exp(epsilon)*stats::pnorm(-sqrt(epsilon*(s+2.0)))
  }

  caseB <- function(epsilon, s){
    stats::pnorm(-sqrt(epsilon*s)) - exp(epsilon)*stats::pnorm(-sqrt(epsilon*(s+2.0)))
  }

  doubling_trick <- function(predicate_stop, s_inf, s_sup){
    while(!predicate_stop(s_sup)){
      s_inf <- s_sup
      s_sup <- 2.0*s_inf
    }
    c(s_inf, s_sup)
  }

  binary_search <- function(predicate_stop, predicate_left, s_inf, s_sup){
    s_mid <- s_inf + (s_sup-s_inf)/2.0
    while(!predicate_stop(s_mid)){
      if (predicate_left(s_mid)) s_sup <- s_mid
      else s_inf <-  s_mid
      s_mid <- s_inf + (s_sup-s_inf)/2.0
    }
    s_mid
  }

  delta_thr <- caseA(epsilon, 0.0)

  if (delta == delta_thr) alpha <- 1.0
  else{
    if (delta > delta_thr){
      predicate_stop_DT <- function(s) caseA(epsilon, s) >= delta
      function_s_to_delta <- function(s) caseA(epsilon, s)
      predicate_left_BS <- function(s) function_s_to_delta(s) > delta
      function_s_to_alpha <- function(s) sqrt(1.0 + s/2.0) - sqrt(s/2.0)
    } else{
      predicate_stop_DT <- function(s) caseB(epsilon, s) <= delta
      function_s_to_delta <- function(s) caseB(epsilon, s)
      predicate_left_BS <- function(s) function_s_to_delta(s) < delta
      function_s_to_alpha <- function(s) sqrt(1.0 + s/2.0) + sqrt(s/2.0)
    }

    predicate_stop_BS <- function(s) abs(function_s_to_delta(s) - delta) <= tol

    s_out <- doubling_trick(predicate_stop_DT, 0.0, 1.0)
    s_inf <- s_out[1]
    s_sup <- s_out[2]
    s_final <- binary_search(predicate_stop_BS, predicate_left_BS, s_inf, s_sup)
    alpha <- function_s_to_alpha(s_final)
  }

  alpha*sensitivity/sqrt(2.0*epsilon)
}
