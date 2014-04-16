\name{coancestry}
\alias{coancestry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate relatedness and inbreeding coefficients}
\description{
Implements Jinliang Wang's code for Coancestry, which allows relatedness to be estimated from codominant genetic data using any of seven estimators, and includes options for considering inbreeding and genotyping errors.
}
\usage{
coancestry(genotype.data, error.rates = 0, allele.freqs = NULL, trioml = 0L, 
wang = 0L, lynchli = 0L, lynchrd = 0L, ritland = 0L, quellergt = 0L, dyadml = 0L, 
ci95.num.bootstrap = 100L, trioml.num.reference = 100L, allow.inbreeding = FALSE, 
rng.seed = NULL, working.directory = tempdir(), output.file = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype.data}{
A data frame containing the genotype data, preferably generated using our \code{\link{readgenotypedata}} function.
}
  \item{error.rates}{
Optional.  If one error rate across all loci is assummed, use that number.  If each locus has a different error rate, create a vector containing the error rate values for each locus, and refer to that vector here.
}
  \item{allele.freqs}{
Optional. If data were read into R using our \code{\link{readgenotypedata}} function, then the allele frequency object will be detected automatically.  If not, then the object you created should be referred to here.
}
  \item{trioml}{
Optional.  The triadic likelihood relatedness estimate (Wang 2007).  If point estimates using this estimator are desired, enter "1", if point estimates and 95\% confidence intervals are desired, enter "2". 
}
  \item{wang}{
Optional.  The relatedness estimate described in Wang (2002).  If point estimates using this estimator are desired, enter "1", if point estimates and 95\% confidence intervals are desired, enter "2". 
}
  \item{lynchli}{
Optional.  The relatedness estimate described in Li et al. (1993).  If point estimates using this estimator are desired, enter "1", if point estimates and 95\% confidence intervals are desired, enter "2".
}
  \item{lynchrd}{
Optional.  The relatedness estimate described in Lynch and Ritland (1999).  If point estimates using this estimator are desired, enter "1", if point estimates and 95\% confidence intervals are desired, enter "2".
}
  \item{ritland}{
Optional.  The relatedness estimate described in Ritland (1996).  If point estimates using this estimator are desired, enter "1", if point estimates and 95\% confidence intervals are desired, enter "2".
}
  \item{quellergt}{
Optional.  The relatedness estimate described in Queller and Goodnight (1989).  If point estimates using this estimator are desired, enter "1", if point estimates and 95\% confidence intervals are desired, enter "2".
}
  \item{dyadml}{
Optional.  The dyadic likelihood estimator, described in Milligan (2003).  If point estimates using this estimator are desired, enter "1", if point estimates and 95\% confidence intervals are desired, enter "2".
}
  \item{ci95.num.bootstrap}{
Optional.  The number of bootstrap iterations to perform to calculate 95\% confidence intervals (default = 100).
}
  \item{trioml.num.reference}{
Optional.  The triadic likelihood estimator requires that you specify the number of reference individuals to use for estimating relatedness.  Enter that number here.  Default = 100.
}
  \item{allow.inbreeding}{
Optional. A logical where inbreeding should, or should not (default), be considered when estimating relatedness.
}
  \item{rng.seed}{
Optional. Can manually set the see of the random number generator.
}
  \item{working.directory}{
Optional.  Can indicate what directory files are in, if not in the current directory.
}
  \item{output.file}{
Optional. Can specify name of the output file prefix (many files will be generated - see below), but can also do this by directing output into an object using standard R commands.
}
}
\value{
\item{relatedness}{A data frame containing all pairwise estimates of relatedness.  This will always have 11 columns: (1) an integer for the pair number; (2) the ID for individual #1; (3) the ID for individual #2; (4) the group assignment (see section 3.5 of accompanying vignette); and (5 - 11) for the 7 relatedness estimators - contain values of 0 for estimators not chosen}

\item{delta7}{A data frame that contains the delta7 estimates for the relatedness estimators that use it (trioml, wang, lynchrd, dyadml). This data frame contains one row for each pairwise comparison, and 8 columns: (1) an integer for the pair number; (2) the ID for individual #1; (3) the ID for individual #2; (4) the group assignment (see section 3.5 of accompanying vignette); and (5 - 8) estimates of delta7 for the 4 relevant estimators, with values of 0 for estimators not chosen.}

\item{delta8}{A data frame that contains the delta8 estimates for the relatedness estimators that use it (trioml, wang, lynchrd, dyadml). This data frame contains one row for each pairwise comparison, and 8 columns: (1) an integer for the pair number; (2) the ID for individual #1; (3) the ID for individual #2; (4) the group assignment (see section 3.5 of accompanying vignette); and (5 - 8) estimates of delta8 for the 4 relevant estimators, with values of 0 for estimators not chosen.} 

\item{inbreeding}{A data frame that contains the inbreeding estimates for each individual, as used in the relatedness estimates.  Only four of the relatedness estimators can account for inbreeding: dyadml, lynchrd, ritland, trioml.  This data frame contains one row for each individual, and 5 columns: (1) individual ID; (2-5) inbreeding estimates for the 4 relatedness estimators.  Estimators not used will have a 0 in the corresponding column.} 

\item{relatedness.ci95}{If confidence intervals are calculated. A data frame containing the low and high cut-off values for the 95\% confidence interval associated with each chosen estimator.  This will always have 18 columns: (1) an integer for the pair number; (2) the ID for individual #1; (3) the ID for individual #2; (4) the group assignment (see section 3.5 of accompanying vignette); (5 - 18) for the high and low values associated with each of the 7 relatedness estimators---contain values of 0 for estimators not chosen.}

\item{delta7.ci95}{If confidence intervals are calculated. A data frame that contains the low and high cut-off values for the 95\% confidence interval for the delta7 estimates associated with each chosen estimator that use it (trioml, wang, lynchrd, dyadml).  This will always have 12 columns: (1) an integer for the pair number; (2) the ID for individual #1; (3) the ID for individual #2; (4) the group assignment (see section 3.5 of accompanying vignette); (5 - 12) for the high and low values associated with each of the 7 relatedness estimators---contain values of 0 for estimators not chosen.}

\item{delta8.ci95}{If confidence intervals are calculated. A data frame that contains the low and high cut-off values for the 95\% confidence interval for the delta8 estimates associated with each chosen estimator that use it (trioml, wang, lynchrd, dyadml).  This will always have 12 columns: (1) an integer for the pair number; (2) the ID for individual #1; (3) the ID for individual #2; (4) the group assignment (see section 3.5 of accompanying vignette); (5 - 12) for the high and low values associated with each of the 7 relatedness estimators---contain values of 0 for estimators not chosen.}

\item{inbreeding.ci95}{If confidence intervals are calculated. A data frame that contains the low and high cut-off values for the 95\% confidence interval for the inbreeding estimates for each individual, as used in the relatedness estimators.  Only four of the relatedness estimators can account for inbreeding: dyadml, lynchrd, ritland, trioml.  This data frame contains one row for each individual, and 9 columns: (1) individual ID; (2-9) inbreeding estimates for the four relatedness estimators.  Estimators not used will have a zero (0) in the corresponding column.}
}
\references{
Li CC, Weeks DE, Chakravarti A (1993) Similarity of DNA fingerprints due to chance and relatedness. \emph{Human Heredity} 43: 45-52.

Lynch M, Ritland K (1999) Estimation of pairwise relatedness with molecular markers. \emph{Genetics} 152: 1753-1766.

Milligan BG (2003) Maximum-likelihood estimation of relatedness. \emph{Genetics} 163: 1153-1167.

Queller DC, Goodnight KF (1989) Estimating relatedness using molecular markers. \emph{Evolution} 43: 258-275.

Ritland (1996) Estimators of pairwise relatedness and inbreeding coefficients. \emph{Genetical Research} 67: 175-186.

Wang J (2002) An estimator of pairwise relatedness using molecular markers. \emph{Genetics} 160: 1203-1215.

Wang J (2007) Triadic IBD coefficients and applications to estimating pairwise relatedness. \emph{Genetical Research} 89: 135-153.

Wang J (2011) COANCESTRY: a program for simulating, estimating and analysing relatedness and inbreeding coefficients. \emph{Molecular Ecology Resources} 11: 141-145.
}

\seealso{
\code{\link{readgenotypedata}}
}
\examples{
	\dontrun{
		#---Read data into R---#
		data(GenotypeData)
		input <- readgenotypedata(GenotypeData)

		#---Calculate Relatedness---#
		output <- coancestry(input$gdata, lynchrd=2, quellergt=2, wang=2)

		#---View Point Estimates---#
		output$relatedness

		#---View 95% Confidence Intervals---#
		output$relatedness.ci95
	}	
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{relatedness}
\keyword{inbreeding}

