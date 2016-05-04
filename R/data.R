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

#' Oral contraceptive data, stratified by age group (2 x 2 x 3 array)
#' 
#' Data from a case-control study investigating oral contraceptives as a risk
#' factor for myocardial infarction, stratified by three age groups.
#' 
#' @format 2 x 2 x 3 array \describe{ 
#'  \item{disease (Myocardial infarction)}{case, control} 
#'  \item{OC (Oral Contraceptive)}{yes, no} 
#'  \item{age group}{25-34, 35-55, 45-49}
#'  }
#' @source Newman (2001), Table 11.5(a); Shapiro, et al. (1979) Oral-contraceptive use in relation to myocardial
#' infarction. \emph{The Lancet} April 7, 743-746
"oral.age"

#' Estrogen data (2 x 2 matrix)
#' 
#' Data from a matched-pairs case-control study investigating estrogen use as a risk factor
#' for endometrial cancer.
#' 
#' @format 2 x 2 matrix
#' \describe{
#' \item{Case}{exposed, unexposed}
#' \item{Control}{exposed, unexposed}
#' }
#' @source Newman (2001), Table 11.10; Antunes, et al. (1979) Endometrial cancer and estrogen use; 
#' Report of a large case-control study. \emph{NEJM} 300, 9-13.
"estrogen"

#' Estrogen data (2 x 5 matrix)
#' 
#' Data from a (1:4) matched-pairs case-control study investigating estrogen use as a risk factor
#' for endometrial cancer.
#' 
#' @format 2 x 5 matrix
#' \describe{
#' \item{Case}{exposed, unexposed}
#' \item{Number of exposed controls}{0,1,2,3,4}
#' }
#' @source Newman (2001), Table 11.14; Mack, et al. (1976). Estrogens and endometrial cancer in a retirement
#' community. \emph{NEJM} 294, 1262-1267.
"estrogen2"


#' Schizophrenia data
#' 
#' Data from a cohort study of mortality in 2122 males who received treatment 
#' for schizophrenia in the province of Alberta, Canada at some time during 1976 -
#' 1985. This is an example of a retrospective cohort study.
#' 
#' @format data frame with 8 rows and 5 variables
#' \describe{
#' \item{age.group}{the age group}
#' \item{cohort.deaths}{number of deaths (of any cause) among schizophrenia cohort}
#' \item{cohort.py}{cohort person years}
#' \item{alberta.deaths}{number of deaths for Alberta males in 1981}
#' \item{alberta.pop}{population of Alberta males in 1981}
#' }
#' @source Newman (2001), Table 12.2(a); Newman and Bland (1991). Mortality in a cohort of patients with 
#' schizophrenia: a record linkage study. \emph{Canadian Journal of Psychiatry} 36, 239 - 245.
"schizophrenia"

#' Female death rates in Canada, all causes
#' 
#' Age-specific death rates for all causes of death in the Candian female population for selected age groups
#' and selected years. These data were taken from official Statistics Canada publications.
#' 
#' @format 5 x 5 matrix
#' \describe{
#'  \item{Year}{1950 - 1990}
#'  \item{Age group}{30-34, 40-44, 50-54, 60-64, 70-74}
#' }
#' @source Newman (2001), Table 12.3.
#' @examples 
#' ## Fig 12.1(a)
#' year <- as.integer(rownames(females))
#' plot(year, females[,1], type="b", ylim = c(0,50), 
#'      ylab = "Death Rate (per 100,000)", xlab="Year")
#' lines(year, females[,2], type = "b", pch = 2)
#' lines(year, females[,3], type = "b", pch = 3)
#' lines(year, females[,4], type = "b", pch = 4)
#' lines(year, females[,5], type = "b", pch = 5)
#' legend("topright", legend = colnames(females), pch = 1:5, title = "Age group")
#'      
#' ## Fig 12.2(b)
#' ag <- unclass(factor(colnames(females)))
#' plot(ag, females[1,], type="b", ylim = c(0,50),
#'      ylab = "Death Rate (per 100,000)", xlab="Age group", axes=F)
#' lines(ag, females[2,], type = "b", pch = 2)
#' lines(ag, females[3,], type = "b", pch = 3)
#' lines(ag, females[4,], type = "b", pch = 4)
#' lines(ag, females[5,], type = "b", pch = 5)
#' legend("top", legend = rownames(females), pch = 1:5, title = "Time period")
#' axis(1, at=1:5, labels = colnames(females))
#' axis(2, at=seq(0,50,10), labels = seq(0,50,10))
#' 
#' ## Fig 12.2(c)
#' ## function to extract cohorts (ie, diagonals)
#' getDiag <- function(m,i,j){
#'   x <- vector(mode = "list", length = length(i))
#'   y <- vector(mode = "list", length = length(i))
#'   k <- pmin(nrow(m) - i, ncol(m) - j)
#'   for(n in seq_along(i)){
#'     y[[n]] <- diag(m[i[n]:(i[n] + k[n]),j[n]:(j[n] + k[n])])    
#'     x[[n]] <- j[n]:(j[n] + k[n])    
#'   }
#'   c(x,y)
#' }
#' ## row and column coordinates for getDiag function
#' i <- c(3,2,1,1,1)
#' j <- c(1,1,1,2,3)
#' dr <- getDiag(females,i,j)
#' plot(dr[[1]], dr[[6]], type="b", xlim = c(1,5), ylim = c(0,50), 
#'      ylab = "Death Rate (per 100,000)", xlab="Age group", axes=F)
#' lines(dr[[2]], dr[[7]], type="b", pch=2)
#' lines(dr[[3]], dr[[8]], type="b", pch=3)
#' lines(dr[[4]], dr[[9]], type="b", pch=4)
#' lines(dr[[5]], dr[[10]], type="b", pch=5)
#' axis(1, at=1:5, labels = colnames(females))
#' axis(2, at=seq(0,50,10), labels = seq(0,50,10))
#' legend("top", legend = seq(1930,1970,10), pch = 1:5, title = "Birth Cohort")
#' rm(i,j,dr, getDiag)
"females"