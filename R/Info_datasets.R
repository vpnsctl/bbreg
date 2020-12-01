#############################################
#' Stress/Axiety data set
#' @description Stress and anxiety scores among nonclinical women in Townsville - Queensland, Australia.
#' @docType data
#' @usage data(SA)
#' @format Data frame containing 166 observations on 2 variables.
#' \describe{
#'   \item{stress}{score, linearly transformed to the open unit interval.}
#'   \item{anxiety}{score, linearly transformed to the open unit interval.}
#' }
#' @references
#' arXiv:2003.05157 (\href{https://arxiv.org/abs/2003.05157}{Barreto-Souza, Mayrink and Simas; 2020})
#'
#' DOI:10.1037/1082-989X.11.1.54 (\href{https://psycnet.apa.org/record/2006-03820-004?doi=1}{Smithson and Verkuilen (2006)})
#'
#' @source Data can be obtained from the supplementary materials of \emph{Smithson and Verkuilen (2006)}. See also \emph{Barreto-Souza, Mayrink and Simas (2020)} for details.
#' @examples data(SA)
"SA"

#############################################
#' Weather Task data set
#' @description Weather task data set.
#' @docType data
#' @usage data(WT)
#' @format Data frame containing 345 observations on 3 variables.
#' \describe{
#'   \item{agreement}{probability or the average between two probabilities indicated by each individual.}
#'   \item{priming}{categorical covariate (0 = two-fold, 1 = seven-fold).}
#'   \item{eliciting}{categorical covariate (0 = precise, 1 = imprecise).}
#' }
#' @references
#' arXiv:2003.05157 (\href{https://arxiv.org/abs/2003.05157}{Barreto-Souza, Mayrink and Simas; 2020})
#'
#' DOI:10.1080/15598608.2009.10411918 (\href{https://www.tandfonline.com/doi/abs/10.1080/15598608.2009.10411918?tab=permissions&scroll=top}{Smithson and Verkuilen; 2009})
#'
#' DOI:10.3102/1076998610396893 (\href{https://journals.sagepub.com/doi/abs/10.3102/1076998610396893?journalCode=jebb}{Smithson et al.; 2011})
#'
#' @source Data can be obtained from supplementary materials of \emph{Smithson et al. (2011)}. See also \emph{Barreto-Souza, Mayrink and Simas (2020)} for details.
#' @examples data(WT)
"WT"

#############################################
#' Body Fat data set
#' @description Penrose body fat data set.
#' Response variable is the percentage of body fat and covariates represent several physiologic measurements related to 252 men.
#' All covariates were rescaled dividing their original value by 100.
#' @docType data
#' @usage data(BF)
#' @format Data frame containing 252 observations on 14 variables.
#' \describe{
#'   \item{bodyfat}{percentage of body fat obtained through underwater weighting.}
#'   \item{age}{age in years/100.}
#'   \item{weight}{weight in lbs/100.}
#'   \item{height}{height in inches/100.}
#'   \item{neck}{neck circumference in cm/100.}
#'   \item{chest}{chest circumference in cm/100.}
#'   \item{abdomen}{abdomen circumference in cm/100.}
#'   \item{hip}{hip circumference in cm/100.}
#'   \item{thigh}{thigh circumference in cm/100.}
#'   \item{knee}{knee circumference in cm/100.}
#'   \item{ankle}{ankle circumference in cm/100.}
#'   \item{biceps}{biceps circumference in cm/100.}
#'   \item{forearm}{forearm circumference in cm/100.}
#'   \item{wrist}{wrist circumference in cm/100.}
#' }
#' @references
#' arXiv:2003.05157 (\href{https://arxiv.org/abs/2003.05157}{Barreto-Souza, Mayrink and Simas; 2020})
#'
#' DOI:10.1249/00005768-198504000-00037 (\href{http://dx.doi.org/10.1249/00005768-198504000-00037}{Penrose et al.; 1985})
#'
#' DOI:10.4236/ojs.2016.61010 (\href{https://www.scirp.org/journal/paperinforcitation.aspx?paperid=63650}{Brimacombe; 2016})
#'
#' @source Data is freely available from \emph{Penrose et al. (1985)}. See also \emph{Brimacombe (2016)} and \emph{Barreto-Souza, Mayrink and Simas (2020)} for details.
#' @examples data(BF)
"BF"
