#######################################################################################################
# Parallelized Dissimilarity measurements for a set of vectors with differing numbers of observations #
#######################################################################################################
# Input can be a list of vectors (with no NA's), a matrix padded with NA's or, an mts matrix padded with NA's.
# n_cores is the number of processor cores to utilize in the calculation.
# Current available dissimilarities are 
# 1) DTWARP, 
# 2) CID (without multiplication by Euclidean distance matrix), 
# 3) CIDDTW (CID is multiplied by DTW matrix), 
# 4) HD (Hausdorff Distance).
# Output is the dissimilarity matrix as a dist object.

diss.UNEQUAL.Parallel <- function(SERIES, METHOD, n_cores) {
	require(parallel)

	# If input is a matrix padded with NA's, convert to a list of vectors with NA's removed
	if(class(SERIES)=="matrix" || class(SERIES)=="mts"){		

		matrix2list1 <- function(matrix1){
			N <- ncol(matrix1)
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

	# Define the database dtw function
	if(METHOD == "DTWARP" || METHOD == "CIDDTW") {
		require(dtw)
		ParallelDTWUnequal <- function(alist, n_cores) {
		cl <- makeCluster(n_cores)
			N <- length(alist)
			mat1 <- matrix(0,N,N)
			x1 <- foreach(i=1:N) %:%
			foreach(j=1:i) %dopar% {
				mat1[i,j] <- dtw(alist[[i]], alist[[j]])$distance
			}
		return(as.dist(mat1))
		}
	}

	# Define the database Hausdorff distance function
	if(METHOD == "HD") {
		ParallelHDUnequal <- function(alist, n_cores) {
		require(pracma)
		cl <- makeCluster(n_cores)
		N <- length(alist)
		mat1 <- matrix(0,N,N)
		x1 <- foreach(i=1:N) %:%
			foreach(j=1:i) %dopar% {
				mat1[i,j] <- hausdorff_dist(alist[[i]], alist[[j]])
			}
		return(as.dist(mat1))
		}
	}	

	# Define the complexity invariant distance function and the database complexity invariant distance function
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

		# Define the database complexity invariant distance function
		ParallelCIDUnequal <- function(alist, n_cores) {
			N <- length(alist)
			cl <- makeCluster(n_cores)
			mat1 <- matrix(0,N,N)
			x1 <- foreach(i=1:N) %:%
			foreach(j=1:i) %dopar% {
				mat1[i,j] <- cid(alist[[i]], alist[[j]])
			}
		return(as.dist(mat1))
		}
	}

	# Run the selected method
	out1 <- switch(METHOD,
				HD = ParallelHDUnequal(SERIES, n_cores),
				CID = ParallelCIDUnequal(SERIES, n_cores),
				CIDDTW = ParallelDTWUnequal(SERIES, n_cores)*ParallelCIDUnequal(SERIES, n_cores),
				DTWARP = ParallelDTWUnequal(SERIES, n_cores))			
	return(out1)
}