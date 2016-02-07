##########################################################################################
# Dissimilarity measurements for a set of vectors with differing numbers of observations #
##########################################################################################
# Input can be a list of vectors (with no NA's), a matrix padded with NA's or, an mts matrix padded with NA's.
# Current available dissimilarities are 
# 1) DTWARP, 
# 2) CID (without multiplication by Euclidean distance matrix), 
# 3) CIDDTW (CID is multiplied by DTW matrix), 
# 4) HD (Hausdorff Distance).
# Output is the dissimilarity matrix as a dist object.

diss.UNEQUAL <- function(SERIES, METHOD) {
	
	# If input is a matrix padded with NA's, convert to a list of vectors with NA's removed
	if(class(SERIES)=="matrix" || class(SERIES)=="mts"){		
		matrix2list1 <- function(matrix1){
			N <- length(matrix1[1,])
			alist <- vector("list", N)
			vec1 <- rep(0,N)
				for(i in 1:N) {
					vec1 <- matrix1[,i] 
					w <- which(is.na(vec1)==TRUE)
					if(length(w) > 0) {vec1 <- vec1[-w]}
					alist[[i]] <- vec1
				}
		return(alist)
		}

	SERIES <- matrix2list1(SERIES) 
	}

	if(METHOD == "DTWARP" || METHOD == "CIDDTW") {require(dtw)}
	if(METHOD == "HD") {require(pracma)}	
	if(METHOD == "CID" || METHOD == "CIDDTW") {	
		# Define the Complexity Invariant Distance (CID) function
		cid <- function(x,y) {
			cesx <- sqrt(sum(diff(x)^2))
			cesy <- sqrt(sum(diff(y)^2))
			denom <- min(cesx, cesy)
			if(denom == 0) {stop("Cannot divide by zero: A series exists that has complexity zero.")}
			else if(denom > 0){cid1 <- max(cesx, cesy)/denom}
			return(cid1)
		}
	}

	# Get the length of the list
	N <- length(SERIES)

	# Preallocate the outer loop matrix to hold the results from the inner loop
	OuterMatrix <- matrix(0, N, N)
	
	# Preallocate the inner vector list to hold the results from the inner loop		
	InnerVector <- rep(0, N)

	# Start the outer part of the nested loop
	for (i in 1:N) {
						
		# Start the inner loop
		for(j in 1:i) {

			# Calculate the dissimilarity distance (gives lower triangular of the matrix)
			if(j < i){
				InnerVector[j] <- switch(METHOD,
							     HD = hausdorff_dist(as.matrix(SERIES[[i]]), as.matrix(SERIES[[j]])),
							     CID = cid(SERIES[[i]], SERIES[[j]]),
							     CIDDTW = cid(SERIES[[i]], SERIES[[j]])*dtw(SERIES[[i]], SERIES[[j]], distance.only=TRUE)$distance, 
							     DTWARP = dtw(SERIES[[i]], SERIES[[j]], distance.only=TRUE)$distance)
			}
		}

		# Save the output from the inner part of the nested loop
		OuterMatrix[i,] <- InnerVector
	}
	
	return(as.dist(OuterMatrix))
}