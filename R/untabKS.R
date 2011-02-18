untabKS <-
function(data, freqvar = NULL){
	dset = data.frame(data)
	if("Freq" %in% names(dset) & is.null(freqvar)){ freqvar = "Freq" }
	dset = data.frame(dset)
	stopifnot(freqvar %in% names(dset))
	ind = which(names(dset) != freqvar)
	fi = which(names(dset) == freqvar)
	names(dset)[fi] = "Freq"
	n = ncol(dset)
	m = nrow(dset)
	X = data.frame(matrix(ncol=n-1,nrow=0))
	zind = which(dset$Freq > 1)
	zero = which(dset$Freq == 0)

	X = sapply(zind,function(x) spread(dset[x,ind],nrow = dset[x,fi]))
	X = do.call("rbind",X)
	vn = names(dset)[ind]
	names(X)=vn
	X = rbind(X,as.matrix(dset[c(1:m)[-c(zero,zind)],ind]))
	return(suppressWarnings(data.frame(X)))
}

spread = function(M,ncol = 1,nrow = 1){
	# 	>>> this function turns each matrix entry into a nrow x ncol - matrix
	M = as.matrix(M)
	M2 = apply(M,2,function(x) rep(x,each = nrow)  )
	if( (nrow(M) == 1)&(nrow == 1) ){
		M2 = matrix(M2,ncol = ncol(M))
	#dim(M2) = dim(M)
	}
	
	M3 = t(apply(t(M2),2,function(x) rep(x,each = ncol)  ))
	if( (ncol(M) == 1)&(ncol == 1) ){
		M3 = matrix(M3,nrow = nrow(M2))
		#dim(M3) = dim(M2)
	}
	return(M3)
}