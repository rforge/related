\name{readgenotypedata}
\alias{readgenotypedata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read a file of genotype data into a data.frame}
\description{
This function will import a genotype file into a proper data frame in R, estimate allele frequencies and store them as an appropriate object for subsequent analyses.
}
\usage{
readgenotypedata(genotype.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype.data}{
The file containing the genotype data to be analyzed.  The file will need to be in R's working directory, and have the following characteristics: (1) It should be a text file (not and Excel file); (2) It should be space- or tab-delimited; (3) Missing data must be represented as zeros (0); and (4) There should not be a header row containing column names.  Column 1 should contain individual identifiers, columns 2 and 3 should contain alleles 1 and 2 for locus 1, columns 4 & 5 should contain alleles 1 and 2 for locus 2, and so on.  Thus, the total number of columns should be 2 x the number of loci + 1.
}
}

\value{
\item{gdata}{The data fram containing the genotype data.  The first column is character data, and the remaining columns are all integers.}
\item{nloci}{An integer containing the number of loci used}
\item{nalleles}{A series of integers specifying the number of alleles at each locus}
\item{ninds}{An integer containing the number of individuals in the genotype file}
\item{freqs}{An object containing the allele frequency data for each locus, which is needed for subsequent analyses}
}

\author{
Jack Pew
}

\seealso{
\code{\link{coancestry}}
}
\examples{
data(GenotypeData)
input <- readgenotypedata(GenotypeData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
