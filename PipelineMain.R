# While the functions are stored in the KeyFunction file, this is where #
# they are actually run. Below is an example of how to run this. #
# However, you can also set the input value to a file path #
source(KeyFunctions.R)

inputdata = 
# Makes spearman, Euclidean, and combined plots #
CombinedPar.R(inputdata)

LV <- function .Last.value
 
# Make list of clusters, move into a function #
MakeClusterList()
