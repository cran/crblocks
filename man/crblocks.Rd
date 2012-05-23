\name{crblocks}
\alias{crblocks}
\alias{catrandstat}
\alias{catrandpvalue}
\alias{catrandpvaluepermute}
\title{Categorical Randomized Block Data Analysis}
\description{
 Implements a statistical test for comparing barplots or histograms of 
 categorical data derived from a randomized block repeated measures layout.
}
\usage{
 catrandstat(rawdata)
 catrandpvalue(datafilename,Nrepeats)
 catrandpvaluepermute(datafilename,Nrepeats)
}
\arguments{
 \item{rawdata}{the data to analyse.}
 \item{datafilename}{a character string giving the name of the data file to analyse.}
 \item{Nrepeats}{the number of Monte Carlo simulated data sets to use in computing the p-value (10000+ recommended).}
}
\details{
 This package implements the statistical test for comparing barplots or histograms of
 categorical data derived from a randomized block repeated measures design described in
 the paper "A Statistical Test for Categorical Randomized Block Sensory Evaluation Data"
 by DJ Best, JCW Rayner and David Allingham (submitted, 2012).  The main functions are
 \code{catrandpvalue} and \code{catrandpvaluepermute}.  They read a dataset from a
 plain-text file can return a p-value, as well as other values of interest, using Monte
 Carlo simulations and permutations, respectively. The function which computes the statistic
 can be called directly if desired.

 \bold{Data format:}

 Using one line of data per judge, each line of the input file contains the category into
 which each product was placed by that judge, with one column for each product. Each judge
 must categorise every product.

 Comments (starting with \#) are allowed (both on their own lines and at the end of lines
 of data). The file should not contain a header of column names: use a comment to include
 such descriptions.

 There are no error checks on the format.  Users should examine the values of
 Njudges and Nproducts in the output to ensure that they are as expected.
}
\value{
 For the \code{catrandstat} function:
 \item{$Njudges}{the number of judges in the data file (number of data lines).}
 \item{$Nproducts}{the number of products tested (number of data columns).}
 \item{$rawdata}{a matrix containing the data that was read from the input file (categories for each product by each judge).}
 \item{$categories}{a vector containing a list of the categories present in the data.}
 \item{$Ncategories}{the number of different categories present in the data (length of $categories).}
 \item{$catCounts}{a matrix containing the number of times each product was placed in each category.}
 \item{$judgeCatCounts}{a matrix containing the number of times each judge used each category.}
 \item{$comparisonStatistic}{the test statistic computed from the ranked data.}

 For the \code{catrandpvalue} function:
 \item{$rawdata}{a matrix containing the data that was read from the input file (categories for each product by each judge).}
 \item{$montecarlodata}{a matrix containing one of the Monte Carlo data sets that was generate.}
 \item{$Ngenerated}{the total number of Monte Carlo data sets generated to give Nrepeats sets with no ties (where a judge places all products into the same category).}
 \item{$Sdata}{the S statistic value for the data read from the file.}
 \item{$Smontecarlo}{the S statistic values for the Monte Carlo data sets.}
 \item{$pvalue}{the p-value for the null hypothesis that there exist no pairwise differences between products.}

 For the \code{catrandpvaluepermute} function:
 \item{$rawdata}{a matrix containing the data that was read from the input file (categories for each product by each judge).}
 \item{$Sdata}{the S statistic value for the data read from the file.}
 \item{$Spermute}{the S statistic values for the Monte Carlo data sets.}
 \item{$pvalue}{the p-value for the null hypothesis that there exist no pairwise differences between products.}
}
\references{
 ``A Statistical Test for Categorical Randomized Block Sensory Evaluation Data'',
 Best, D.J., Rayner, J.C.W. and Allingham, David. Journal of Sensory Studies, submitted, 2012.
}
\author{
 Allingham, David \email{David.Allingham@newcastle.edu.au}

 Best, D.J. \email{John.Best@newcastle.edu.au}
}
\examples{
 ### Analyse the sample dataset provided with this package:
 # Load the data from the file and compute its test statistic
 inputfile = system.file('extdata', 'omahony.txt', package='crblocks')
 omahonydata=read.table(file(inputfile,'r'))
 output = catrandstat(omahonydata)
 output$comparisonStatistic

 # Load the data from the file and compute the p-value for its test statistic using Monte Carlo simulation:
 output = catrandpvalue(inputfile,1000)
 output$pvalue

 # Load the data from the file and compute the p-value for its test statistic using Monte Carlo simulation:
 output = catrandpvaluepermute(inputfile,1000)
 output$pvalue
}
\keyword{barplot comparison}
\keyword{computer-intensive p-values}
\keyword{histogram comparison}
\keyword{repeated measures data}
\keyword{nominal data}
