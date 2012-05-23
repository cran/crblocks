catrandstat <- function(rawdata){
#######################################################################
# Function to compute the statistic                                   #
#######################################################################

### Extract the meta variables from the data:
 ### Njudges is the number of rows:
  Njudges <- nrow(rawdata)
 ### Nproducts is the number of columns:
  Nproducts <- ncol(rawdata)

### Extract the category labels that are present:
 if (is.matrix(rawdata)){
  ### for the Monte Carlo data (which is a matrix):
   categories <- sort(unique(sort(rawdata)))
 } else {
  ### for the data read from file (which is a data frame):
   categories <- sort(unique(stack(rawdata)[,1]))
 }
 ### Count how many we have:
 Ncategories <- length(categories)

 ### Construct the product-X-categories counts:
  catCounts <- matrix(NaN,Nproducts,Ncategories) # initialise
  for (i in 1:Nproducts){
   for (j in 1:Ncategories){
    catCounts[i,j] <- length(which(rawdata[,i]==categories[j]))
   }
  }
 catColSums <- apply(catCounts,2,sum)

 ### Construct the d vectors (as columns of a matrix:):
  d <- matrix(NaN,Ncategories-1,Nproducts) # initialise
  for (i in 1:Nproducts){
   d[,i]=(Nproducts*catCounts[i,])[1:(Ncategories-1)]-catColSums[1:(Ncategories-1)]
 }

 ### Construct the judges-X-categories counts (judgeCatCounts):
  U <- matrix(NaN,Njudges,Ncategories) # initialise
  for (i in 1:Njudges){
   for (j in 1:Ncategories){
    U[i,j] <- length(which(rawdata[i,]==categories[j]))
   }
  }
  judgeCatCountsFull=U
  U=U[,1:Ncategories-1] # for the statistic we only want the first Ncategories-1 columns

 ### Construct the covariance matrix, V:
  V <- matrix(NaN,Ncategories-1,Ncategories-1) # initialise
  for (i in 1:Ncategories-1){
   for (j in 1:Ncategories-1){
    if (i==j){
     V[i,j] <- Nproducts*sum(U[,i])-crossprod(U[,i])
    } else {
     V[i,j] <- -crossprod(U[,i],U[,j])
    }
   }
  }
 ### find the inverse of V: (first test if it is singular (ie. det(V)==0) (or almost))
  if (det(V)>1e-2){
   ### not singular? okay, compute the inverse:
    Vinverse <- qr.solve(V)

   ### Compute the statistic:
    S <- 0  # initialise
    for (i in 1:Nproducts){
     S <- S + (Nproducts-1)/Nproducts*crossprod(t(crossprod(d[,i],Vinverse)),d[,i])
    }
  } else {
   ### flag the singular matrix by returning NaN
    S <- NaN
    Vinverse <- NaN
  }


#######################################################################
# Return some variables of interest                                   #
#######################################################################
 list(
  Njudges=Njudges,
  Nproducts=Nproducts,
  rawdata=rawdata,
  categories=categories,
  Ncategories=Ncategories,
  catCounts=catCounts,
  judgeCatCounts=judgeCatCountsFull,
  comparisonStatistic=S
 )

} ### end of function











catrandpvalue <- function(datafilename,Nrepeats){
#######################################################################
# Function to compute the p-value for the data                        #
#######################################################################

 ### Check if the data file exists:
  if(!file.exists(datafilename))
   stop('File \'',datafilename,'\' not found')

 ### Open the file read-only:
  inputfile  <-  file(datafilename,"r")
 ### Read the data:
 ### (note that comments (starting with #) are allowed in the data file)
  rawdata  <-  read.table(inputfile, header=FALSE)
 ### Close the file:
  close(inputfile)

 #######################################################################
 # Compute the statistic for the data                                  #
 #######################################################################
  dataoutput <- catrandstat(rawdata)
  Sdata <- dataoutput$comparisonStatistic

 #######################################################################
 # Run the Monte Carlo simulations to compute a p-value                #
 #######################################################################
  Smontecarlo=matrix(NaN,Nrepeats,1) # initialise
  Ngenerated=0
  montecarlodata=array(NaN,c(dataoutput$Njudges,dataoutput$Nproducts,Nrepeats))
  for (i in 1:dataoutput$Njudges){
   indx = 1:Nrepeats # initialise
   while (length(indx)){
    Ngenerated=Ngenerated+length(indx) # count how many we actully generate
#wrong:    montecarlodata[i,,indx]=dataoutput$categories[which(rmultinom(length(indx)*dataoutput$Nproducts,1,dataoutput$judgeCatCounts[i,])==1)%%dataoutput$Ncategories + 1]
    montecarlodata[i,,indx]=dataoutput$categories[apply(rmultinom(length(indx)*dataoutput$Nproducts,1,dataoutput$judgeCatCounts[i,])==1,2,which)]
    ### we have the data, now loop over them and see which are tied, so that we can replace them (those in 'indx'):
    indx=indx[which(apply(apply(apply(montecarlodata[i,,indx],2,sort),2,diff),2,sum)==0)]
    if (length(indx)==1){ # if there is only one element in indx, add another since R doesn't handle the multidimensional matrix indexing well with only one element
     if (indx[1]==1){
      indx[2]=2
     } else {
      indx[2]=1
     }
    }
   }
  }

 ### Compute the statistic for each Monte Carlo data set:
  for (n in 1:Nrepeats){
   Smontecarlo[n]=catrandstat(montecarlodata[,,n])$comparisonStatistic
  }

pvalue=length(Smontecarlo[Smontecarlo>c(Sdata)])/length(Smontecarlo)

### Truncate the Monte Carlo data so that we can return a sample data set
# montecarlodata=montecarlodata[,,1]

#######################################################################
# Return some variables of interest                                   #
#######################################################################
 list(
  rawdata=rawdata,
  montecarlodata=montecarlodata,
  Ngenerated=Ngenerated,
  Sdata=Sdata,
  Smontecarlo=Smontecarlo,
  pvalue=pvalue
 )
} ### end of function








catrandpvaluepermute <- function(datafilename,Nrepeats){
#######################################################################
# Function to compute the p-value for the data using permutation test #
#######################################################################

 ### Check if the data file exists:
  if(!file.exists(datafilename))
   stop('File \'',datafilename,'\' not found')

 ### Open the file read-only:
  inputfile  <-  file(datafilename,"r")
 ### Read the data:
 ### (note that comments (starting with #) are allowed in the data file)
  rawdata  <-  read.table(inputfile, header=FALSE)
 ### Close the file:
  close(inputfile)

 #######################################################################
 # Compute the statistic for the data                                  #
 #######################################################################
  dataoutput <- catrandstat(rawdata)
  Sdata <- dataoutput$comparisonStatistic

 ### Compute the statistic for each permuted data set:
  Spermute=matrix(NaN,Nrepeats,1) # initialise
  for (n in 1:Nrepeats){
   permdata=t(apply(rawdata,1,sample)) ### this permutes each row separately (but transposes too)
   Spermute[n]=catrandstat(permdata)$comparisonStatistic
  }

 pvalue=length(Spermute[Spermute>c(Sdata)])/length(Spermute)

#######################################################################
# Return some variables of interest                                   #
#######################################################################
 list(
  rawdata=rawdata,
  Sdata=Sdata,
  Spermute=Spermute,
  pvalue=pvalue
 )
} ### end of function
