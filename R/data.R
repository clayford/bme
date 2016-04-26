#' Antibody data (2 x 2 matrix)
#' 
#' Portion of data from a cohort study conducted in Bangladesh which investigated whether antibodies present 
#' in breast milk protect infants from diarrhea due to cholera. Table 4.3, page 97 (Newman)
#' 
#' @format a 2 x 2 matrix:
#' \describe{
#'  \item{diarrhea}{yes, no}
#'  \item{antibody}{yes, no}
#' }
#' @source Glass, et al. Protection against cholera in breast-fed children by antibiotics in breast milk. \emph{NEJM} 308, 1389-1392.
"antibody"


#' Breast cancer data (2 x 2 matrix)
#' 
#' A random sample of 199 female breast cancer patients who registered with the Northern Alberta Breast 
#' Cancer Registry in 1985. Entry into the cohort was restricted to women with either stage I, II, or III 
#' disease. Of the 199 subjects in the cohort, seven died of a cause other than breast cancer. These 
#' individuals were dropped from the analysis. Data are provided in different subsets and formats 
#' throughout Newman: Table 4.5(a), page 99 (\code{\link{breast}}); 
#' Table 5.3, page 126 (\code{\link{breast.receptor}}); 
#' Table 5.10, page 140 (\code{\link{breast.stage}}); 
#' Table 9.1, page 175 (\code{\link{breast.survival}}) 
#' 
#' @format 2 x 2 matrix
#' \describe{
#'  \item{survival}{dead, alive}
#'  \item{receptor.level}{amount of estrogen receptor that is present in breast tissue: low, high}
#' }
#' @source Newman (2001)
"breast"

#' Breast cancer data, stratified by receptor level (2 x 3 x 2 array)
#' 
#' A random sample of 199 female breast cancer patients who registered with the Northern Alberta Breast 
#' Cancer Registry in 1985. Entry into the cohort was restricted to women with either stage I, II, or III 
#' disease. Of the 199 subjects in the cohort, seven died of a cause other than breast cancer. These 
#' individuals were dropped from the analysis. Data are provided in different subsets and formats 
#' throughout Newman: Table 4.5(a), page 99 (\code{\link{breast}}); 
#' Table 5.3, page 126 (\code{\link{breast.receptor}}); 
#' Table 5.10, page 140 (\code{\link{breast.stage}}); 
#' Table 9.1, page 175 (\code{\link{breast.survival}}) 
#' 
#' @format \code{breast} 2 x 3 x 2 array
#' \describe{
#'  \item{survival}{dead, alive}
#'  \item{stage}{stage of disease: I, II, III}
#'  \item{receptor.level}{amount of estrogen receptor that is present in breast tissue: low, high}
#' }
#' @source Newman (2001)
"breast.receptor"


#' Breast cancer data, stratified by disease stage (2 x 2 x 3 array)
#' 
#' A random sample of 199 female breast cancer patients who registered with the Northern Alberta Breast 
#' Cancer Registry in 1985. Entry into the cohort was restricted to women with either stage I, II, or III 
#' disease. Of the 199 subjects in the cohort, seven died of a cause other than breast cancer. These 
#' individuals were dropped from the analysis. Data are provided in different subsets and formats 
#' throughout Newman: Table 4.5(a), page 99 (\code{\link{breast}}); 
#' Table 5.3, page 126 (\code{\link{breast.receptor}}); 
#' Table 5.10, page 140 (\code{\link{breast.stage}}); 
#' Table 9.1, page 175 (\code{\link{breast.survival}}) 
#' 
#' @format \code{breast} 2 x 2 x 3 array
#' \describe{
#'  \item{survival}{dead, alive}
#'  \item{receptor.level}{amount of estrogen receptor that is present in breast tissue: low, high}
#'  \item{stage}{stage of disease: I, II, III}
#' }
#' @source Newman (2001)
"breast.stage"

#' Breast cancer survival data (data frame)
#' 
#' A random sample of 199 female breast cancer patients who registered with the Northern Alberta Breast 
#' Cancer Registry in 1985. Entry into the cohort was restricted to women with either stage I, II, or III 
#' disease. This version of the breast cancer data includes the seven subjects who died of a cause other than
#' breast cancer.
#' Data are provided in different subsets and formats throughout Newman: 
#' Table 4.5(a), page 99 (\code{\link{breast}}); 
#' Table 5.3, page 126 (\code{\link{breast.receptor}}); 
#' Table 5.10, page 140 (\code{\link{breast.stage}}); 
#' Table 9.1, page 175 (\code{\link{breast.survival}}) 
#' 
#' @format data frame with 199 observations and 4 variables:
#' \describe{
#'  \item{time}{survival time in months}
#'  \item{status}{censoring variable; 0 = censored; 1 = died of breast cancer}
#'  \item{stage}{stage of disease: I, II, III}
#'  \item{receptor.level}{amount of estrogen receptor that is present in breast tissue: low, high}
#' }
#' @source Newman (2001)
"breast.survival"

#' Ovarian cancer data
#' 
#' Data from a cohort study of women with stage II or IIIA ovarian cancer, where the endpoint is progression of
#' disease.
#' 
#' @format data frame with 35 obervations and 3 variables:
#' \describe{
#'  \item{grade}{indicator of the malignant potential of the tumor: High, Low}
#'  \item{time}{survival time measured in days}
#'  \item{status}{censoring variable; 0 = censored; 1 = died of ovarian cancer}
#' }
#' @source Newman (2001), Table 9.6; Fleming, et al. Modified Kolmogorov-Smirnov test procedures with application to arbitrarily 
#' right-censored data. \emph{Biometrics} 36, 607-625.
"ovarian.cancer"

#' Oral contraceptive data (2 x 2 matrix)
#' 
#' Data from a case-control study investigating oral contraceptives as a risk
#' factor for myocardial infarction.
#' 
#' @format 2 x 2 matrix \describe{ 
#'  \item{Myocardial infarction}{case, control} 
#'  \item{Oral Contraceptive}{yes, no} }
#' @source Newman (2001), Table 11.2; Shapiro, et al. (1979) Oral-contraceptive use in relation to myocardial
#' infarction. \emph{The Lancet} April 7, 743-746
"oral"