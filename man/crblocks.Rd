\name{crblocks}
\alias{crblocks}
\alias{catrandstat}
\alias{catrandpvalue}
\alias{catrandpvaluepermute}
\alias{print.crblocks_output}
\title{Categorical Randomized Block Data Analysis}
\description{
 Implements a statistical test for comparing barplots or histograms of 
 categorical data derived from a randomized block repeated measures layout.
}
\usage{
 catrandstat(rawdata)
 catrandpvalue(datafilename,Nrepeats)
 catrandpvaluepermute(datafilename,Nrepeats)
 \method{print}{crblocks_output}(x,\dots)
}
\arguments{
 \item{rawdata}{the data to analyse.}
 \item{datafilename}{a character string giving the name of the data file to analyse.}
 \item{Nrepeats}{the number of Monte Carlo simulated data sets to use in computing the p-value (10000+ recommended).}
 \item{x}{output from catrandstat, catrandpvalue or catrandpvaluepermute}
 \item{\dots}{not used}
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
 \item{$Sstatistic}{the S statistic computed for the data.}
 \item{$Mstatistic}{the M statistic computed for the data.}
 \item{$L2statistic}{the L^2 statistic computed for the data.}
 \item{$Schi2pvalue}{the chi^2 p-value of the S statistic for the data.}
 \item{$Mchi2pvalue}{the chi^2 p-value of the M statistic for the data.}
 \item{$L2chi2pvalue}{the chi^2 p-value of the L^2 statistic for the data.}

 For the \code{catrandpvalue} function:
 \item{$rawdata}{a matrix containing the data that was read from the input file (categories for each product by each judge).}
 \item{$Nproducts}{the number of products tested (number of data columns).}
 \item{$Ncategories}{the number of different categories present in the data (length of $categories).}
 \item{$Njudges}{the number of judges in the data file (number of data lines).}
 \item{$Ngenerated}{the number of Monte Carlo data sets generated in total to produce Nrepeats data sets with no ties (where a judge places all products into the same category).}
 \item{$Sdata}{the S statistic computed for the data.}
 \item{$Mdata}{the M statistic computed for the data.}
 \item{$L2data}{the L^2 statistic computed for the data.}
 \item{$Smontecarlo}{a vector containing the S statistic values computed for each Monte Carlo data set.}
 \item{$Mmontecarlo}{a vector containing the M statistic values computed for each Monte Carlo data set.}
 \item{$L2montecarlo}{a vector containing the L^2 statistic values computed for each Monte Carlo data set.}
 \item{$Spvalue}{the Monte Carlo p-value for the null hypothesis that there exist no pairwise differences between products based on the S statistic.}
 \item{$Mpvalue}{the Monte Carlo p-value for the null hypothesis that there exist no pairwise differences between products based on the M statistic.}
 \item{$L2pvalue}{the Monte Carlo p-value for the null hypothesis that there exist no pairwise differences between products based on the L^2 statistic.}
 \item{$Schi2pvalue}{the chi^2 p-value of the S statistic for the data.}
 \item{$Mchi2pvalue}{the chi^2 p-value of the M statistic for the data.}
 \item{$L2chi2pvalue}{the chi^2 p-value of the L^2 statistic for the data.}

 For the \code{catrandpvaluepermute} function:
 \item{$rawdata}{a matrix containing the data that was read from the input file (categories for each product by each judge).}
 \item{$Nproducts}{the number of products tested (number of data columns).}
 \item{$Ncategories}{the number of different categories present in the data (length of $categories).}
 \item{$Njudges}{the number of judges in the data file (number of data lines).}
 \item{$Sdata}{the S statistic computed for the data.}
 \item{$Mdata}{the M statistic computed for the data.}
 \item{$L2data}{the L^2 statistic computed for the data.}
 \item{$Spermute}{a vector containing the S statistic values computed for each permuted data set.}
 \item{$Mpermute}{a vector containing the M statistic values computed for each permuted data set.}
 \item{$L2permute}{a vector containing the L^2 statistic values computed for each permuted data set.}
 \item{$Spvalue}{the permutation p-value for the null hypothesis that there exist no pairwise differences between products based on the S statistic.}
 \item{$Mpvalue}{the permutation p-value for the null hypothesis that there exist no pairwise differences between products based on the M statistic.}
 \item{$L2pvalue}{the permutation p-value for the null hypothesis that there exist no pairwise differences between products based on the L^2 statistic.}
 \item{$Schi2pvalue}{the chi^2 p-value of the S statistic for the data.}
 \item{$Mchi2pvalue}{the chi^2 p-value of the M statistic for the data.}
 \item{$L2chi2pvalue}{the chi^2 p-value of the L^2 statistic for the data.}
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
  closeAllConnections()
  catrandstat(omahonydata)
### OUTPUT:
#
# Statistic   dof   data value   chi^2 p-value
#  S           6     13.16        0.04058
#  M           3     11.42        0.003311
#  L^2         1     6.671        0.009799
#

 # Load the data from the file and compute the p-value for
 # its test statistic using Monte Carlo simulation:
  catrandpvalue(inputfile,500)
### SAMPLE OUTPUT:
#
# Statistic   dof   data value   chi^2 p-value   Simulated p-value
#  S           6     13.16        0.04058         0.018
#  M           3     11.42        0.003311        0.002
#  L^2         1     6.671        0.009799        0.008
#

 # Load the data from the file, compute the p-value for
 # its test statistic using Monte Carlo simulation, and
 # store the output variables in X:
  Nrepeats = 500
  X = catrandpvalue(inputfile,Nrepeats)
 # This will be a number greater than Nrepeats:
  X$Ngenerated
### SAMPLE OUTPUT:
#
# [1] 6651

 # Load the data from the file and compute the p-value for
 # its test statistic using Monte Carlo simulation:
 catrandpvaluepermute(inputfile,500)
### SAMPLE OUTPUT:
#
# Statistic   dof   data value   chi^2 p-value   Simulated p-value
#  S           6     13.16        0.04058         0.032
#  M           3     11.42        0.003311        0.004
#  L^2         1     6.671        0.009799        0.006
#

}
\keyword{barplot comparison}
\keyword{computer-intensive p-values}
\keyword{histogram comparison}
\keyword{repeated measures data}
\keyword{nominal data}
