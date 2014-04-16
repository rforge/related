\name{familysim}
\alias{familysim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate related individuals}
\description{
Uses observed allele frequencies to generated a user-defined number of pairs of individuals of known relatedness.
}
\usage{
familysim(freqs, ninds = 100L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{freqs}{
An object containing allele frequency data.  Can be generated with the \code{\link{readgenotypedata}} function.
}
  \item{ninds}{
The number of pairs of individuals to simulate for each relatedness values.  For example, entering 100 would generate 100 parent-offspring pairs, 100 full-sib pairs, 100 half-sib pairs, and 100 unrelated pairs.
}
}
\value{
A data frame containing the genotypes of simulated individuals.
}
\author{
Tim Frasier <timothy.frasier@smu.ca>
}

\seealso{
\code{\link{readgenotypedata}}
}
\examples{
	\dontrun{
		#---Load data set---#
		data(GenotypeData)
		input <- readgenotypedata(GenotypeData)

		#---Generate simulated individuals---#
		sim <- familysim(input$freqs, 100)
	}	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
