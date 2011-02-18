#_____________________________________________________________________#
#-------------------------SUMMARY FOR KS------------------------------#

summary.KS=function(object,npattern=5,...){
	KS=object
	ans=list()
	class(ans)="summary.KS"
	
	ans$domain=KSdomain(KS)
	
	if (!is.null(attributes(KS)$memberships)){
		KSdf=convert(KS,return.dataframe=TRUE)
		ans$Freq=colSums(KSdf[,1:(ncol(KSdf)-1)]*KSdf[,ncol(KSdf)])
		names(ans$Freq) = ans$domain
		ans$N = sum(KSdf[,ncol(KSdf)])
		ord=order(KSdf$Freq,decreasing=TRUE)
		KSdf=KSdf[ord,]
		rownames(KSdf)=NULL
		ans$pattern = KSdf[1:npattern,]
	}else{
		KSdf=convert(KS,return.dataframe=TRUE)
		ans$Freq=colSums(KSdf)
		names(ans$Freq) = ans$domain
	}
	
	KS = as.set(KS)
	
	ans$notions=KSnotions(KS)
	ans$Nstates=length(KS)
	ans$atoms=KSatoms(KS)
	ans$KSpace=is.KSpace(KS)
	
	if (ans$KSpace){
		ans$base=KSbase(KS)
	}
	
	return(ans)
}

#_____________________________________________________________________#
#-------------------------IS.KS FUNCTIONS-----------------------------#

is.KS = function(x, domain=NULL,...){
	UseMethod("is.KS")
}

is.KS.set = function(x,domain=NULL,...){
	x = asCharSet(x)
	if (is.null(domain)){
		domain = as.set(unlist(x))
	}
	domain = asCharSet(set(domain))
	items = as.set(unlist(x))
	if(!set_is_subset(items,as.set(unlist(domain)))){
		x = KStrace(x, domain)
	}
	return((set(set()) %in% x) & (domain %in% x))
}

is.KS.gset = function(x,domain=NULL,...){
	return(is.KS.set(as.set(x),domain))
#	if (is.null(domain)){
#		domain = as.set(unlist(x))
#	}
#	items = as.set(unlist(x))
#	stopifnot(set_is_subset(items,domain))
#	return((set(set()) %in% x) & (set(domain) %in% x))
}

is.KS.matrix = function(x,domain=NULL,...){
	return(is.KS.set(convert(x),domain))
#	if (is.null(domain)){
#		domain = as.set(unlist(x))
#	}
#	items = as.set(unlist(x))
#	stopifnot(set_is_subset(items,domain))
#	return((set(set()) %in% x) & (set(domain) %in% x))
}

#_____________________________________________________________________#
#--------------------------AS.KS FUNCTIONS----------------------------#

as.KS = function(x,type = "set", space= TRUE,...){
	UseMethod("as.KS")
}

as.KS.set = function(x,type = "set",...){
	if( !(type %in% c("set","matrix")) ){
		stop("wrong type argument")	
	}
	if(type == "set"){
		ks = asCharSet(x)
	}
	if(type == "matrix"){
		ks = convert(x)
	}
	class(ks)=unique(c("KS",class(ks)))
	return(ks)
}
as.KS.gset = function(x,type = "set",...){

	if( !(type %in% c("set","matrix")) ){
		stop("wrong type argument")	
	}
	if(type == "set"){
		ks = asCharSet(x)
	}
	if(type == "matrix"){
		ks = convert(x)
	}
	class(ks)=unique(c("KS",class(ks)))
	return(ks)
}	
as.KS.relation = function(x,type = "set",space = TRUE,...){
	
	RI = relation_incidence(x)
	diag(RI) = 1
#	rownames(RI)=strsplit(rownames(RI),split="L")
	
    ks = apply(RI,2, function(y){
				as.set(rownames(RI)[ which( y == 1 ) ])
				})
	names(ks) = NULL
	ks = as.set(as.set(ks) + set(set()) + set(as.set(unlist(ks))))
	if(space){
		ks = closure(ks, operation = "union")
	}
			
	if(type == "matrix"){
		ks = convert(ks)
	}
	class(ks)=unique(c("KS",class(ks)))
	return(ks)
}
as.KS.matrix = function(x,type = "set",...){
	# maybe a check for rownames here?
	if(is.null(rownames(x))){
		rownames(x) = as.vector(sapply(c("",1:(floor(dim(x)[1]/26))),function(y) paste(letters,y,sep="")))[1:dim(x)[1]]
	}
	ks = (x > 0) + 0
	if( !(type %in% c("set","matrix")) ){
		stop("wrong type argument")	
	}
	if(type == "set"){
		ks = convert(ks)
	}
	if(type == "matrix"){
		ks = ks
	}
	class(ks)=unique(c("KS",class(ks)))
	return(ks)
}

#_____________________________________________________________________#
#---------------------------KS -> RELATION----------------------------#

as.relation.KS = function(x, empirical = FALSE, v = 1,...){
	if (inherits(x,"set")){
		empirical = FALSE
	}
	if(empirical == FALSE){
		if(inherits(x,"gset")){
		x = convert(x)
		}
		mat = apply(x, 1, function(y) {
						apply(as.matrix(x[,which(as.logical(y))],nrow=nrow(x)), 1, prod)})
		rel = endorelation(incidence = mat)
	}else{
		if(inherits(x,"gset") & !inherits(x,"set")){
			dat = untabKS(convert(x,return.dataframe=TRUE))
			names(dat) = unlist(KSdomain(x))
			graph=iita(dat,v=v)$implications
			rel0 = endorelation(domain=as.list(1:length(KSdomain(x))),graph=graph)
			rel0 = relabel.relation(rel0,data=dat)
			mat=relation_incidence(rel0)
			diag(mat)=1
#			rownames(mat)=strsplit(rownames(mat),split="L")
#			colnames(mat)=strsplit(colnames(mat),split="L")
			rel = endorelation(incidence=mat)
		}
	}
	return (rel)
}

#_____________________________________________________________________#
#----------------------------- NOTIONS -------------------------------#

KSnotions = function(KS){
	stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
#	stopifnot(is.KS(KS))
	if( inherits(KS,"set") | inherits(KS,"gset") ){
		KS = convert(KS)
	}
	
		cont = apply(KS,1,function(x){
			apply(t(KS) - x, 2, function(y){
				return(max(y) < 1) # TRUE if x-row contains row 	
			})
		})
		notionmat = cont*t(cont)
		mode(notionmat) = "logical"
		notions = unique(apply(notionmat,2,function(x){
							#return( do.call(paste,c(as.list(rownames(KS)[which(x)]),sep="")))
							return( as.list(rownames(KS)[which(x)]))
					}))
	return(lapply(notions,as.set))
}

#_____________________________________________________________________#
#----------------------------- ATOMS ---------------------------------#

KSatoms = function(KS){
#	stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	stopifnot(is.KS(KS))
	if( inherits(KS,"matrix") ){
		KS = convert(KS)
	}
	if( inherits(KS,"gset") & !inherits(KS,"set")){
		KS = as.set(KS)
	}
	items = KSdomain(KS)
	KSlist = as.list(KS - set(set()))
	contained = sapply(items,function(x) sapply(KSlist, function(y) x %in% y))
	
	subs = sapply(KSlist,function(x) sapply(KSlist, function(y) set_is_subset(x,y) ) )
	atom_cont = apply(contained,2,function(x){
			return(colSums(x*t(subs)) == 1)
	})
#set_is_subset macht evtl probleme		
	atomlist = apply(atom_cont,2,function(x){
			as.set(lapply(which(as.logical(x)), function(y) KSlist[[y]]))
		})
	names(atomlist) = items
	return(atomlist)
}

#_____________________________________________________________________#
#-----------------------------  TRACE --------------------------------#

KStrace = function(KS,subs){
	stopifnot(inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix"))
	if( inherits(KS,"matrix") ){
		type = "matrix"
		KS = convert(KS)
	}
	if( inherits(KS,"set") ){
		type = "set"
	}
	KS = asCharSet(KS)
	domain=KSdomain(KS)
	subs = as.set(as.character(unlist(subs)))
#	domain = as.set(as.character(unlist(domain)))
#	stopifnot(subs<=domain)
	if( subs>=domain){
		return(KS)
	}
	if (inherits(KS,"gset") & (!inherits(KS,"set"))){
		KSdf = convert(KS, return.dataframe = TRUE)
		dims=dim(KSdf)
		colnames(KSdf)[1:(dims[2]-1)] = domain
		colnames(KSdf)[dims[2]] = "Freq"
		tracedf = subKS(KSdf,cols=which(colnames(KSdf) %in% subs), keep.zero = TRUE)
		#tracedf = as.data.frame(xtabs(KSdf$Freq~.,data = KSdf[,which(colnames(KSdf) %in% subs)]))
		dims2 = dim(tracedf)
		colnames(tracedf)[1:(dims2[2]-1)] = subs
		colnames(tracedf)[dims2[2]] = "Freq"
		ret = as.gset(tracedf)
		return(ret)
	}
	KS = convert(KS)
	if (inherits(KS,"matrix")){
		tracemat = KS[which(rownames(KS) %in% subs),]
		if (length(subs) == 1){
			ret = matrix(unique(tracemat),ncol=2,nrow=1)
			rownames(ret) = unlist(subs)
		}else{
			ret = t(unique(t(tracemat)))
		}
		if (type == "matrix"){
			class(ret) = c("KS",class(ret))
			return(ret)
		}else{
			if (type == "set"){
				ret = convert(ret)
				return(ret)
			}
		}
	}
}

#_____________________________________________________________________#
#------------------------------ BASE ---------------------------------#

KSbase = function(KS){
#	stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	stopifnot(is.KS(KS))
	if ((inherits(KS,"gset"))&(!inherits(KS,"set"))){
		KS = as.set(KS)
	}
	stopifnot(is.KSpace(KS))
	base = do.call(c,KSatoms(KS))
	names(base) = NULL
	return(base)
}

#_____________________________________________________________________#
#------------------------------ DOMAIN -------------------------------#

KSdomain = function(KS){
 	stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
#	stopifnot(is.KS(KS))
	if( inherits(KS,"set") | inherits(KS,"gset") ){
		dm = unique(unlist(KS))
	}else{
		dm = rownames(KS)
	}
	return(as.set(dm))
}


#_____________________________________________________________________#
#---------------------- DISCRIMINATIVE REDUCTION ---------------------#

discred = function(KS,concat = ":"){
	stopifnot((inherits(KS,"matrix"))|(inherits(KS,"gset")))
	if (inherits(KS,"gset")){
		if (inherits(KS,"set")){
			KS=as.set(KS)
			type="set"
		}else{
			memberships=attributes(KS)$memberships
			KS=as.set(KS)
			type="gset"
		}
	}else{
		KS = convert(KS)
		type="matrix"
	}
	
	nlist = KSnotions(KS)
	
	MKS = unique(convert(KS))
	rownames(MKS) = sapply(nlist, function(x) do.call(paste, c(as.list(x), sep=concat)) )
	
	#for( i in 1:nrow(MKS) ){
	#	print( sapply(nlist,function(x) rownames(MKS)[i] %in% x ) )
		#dimnames(MKS)[[1]][i] =  do.call(  paste, c(   nlist[[ which( sapply(nlist,function(x) which(rownames(MKS)[i] %in% x)))]] , sep=""  )  )
	#}
	ret=MKS
	if (type=="set"){
		ret=convert(MKS)
	}
	if (type=="gset"){
		mat=rbind(MKS,memberships)
		dat=as.data.frame(t(mat))
		names(dat)[dim(dat)[2]]="Freq"
		ret=as.gset(dat)
	}
	return(ret)
} 
#_____________________________________________________________________#
#--------------------------KS CONVERSION------------------------------#

convert = function(KS, domain = NULL, return.dataframe=FALSE){
	if (is.matrix(KS)){
		
		if(is.null(row.names(KS))){
			items = letters[c(1:nrow(KS))]#should that be changed?
		}else{
			items = row.names(KS)
		}		
		tmp = apply(KS,2,function(x) as.set(items[which(as.logical(x))]) )
		names(tmp) = rep("", length(tmp))
		
		return( as.KS(as.set( tmp ))  )
	}
	if ((!is.set(KS)) & (is.gset(KS))){
		Freq=attributes(KS)$memberships
		KS=as.set(KS)
	}else{
		Freq=matrix(nrow=length(KS),ncol=0)
	}
	if (is.set(KS)){
		# also accepts class kstructure
		#print("Converting KS set to matrix design...")
		n = length(KS)
		items = unique(unlist(KS))
		items = sort(items)
		if(!is.null(domain)){
				stopifnot( min(items %in% domain) > 0  )
				items = domain
				if( n == 0 ){ 
					mat = matrix(nrow=length(items),ncol=0)
					dimnames(mat)[[1]] = items
					return(mat)
				}
		}else{
			stopifnot( !(n==0) )
		}
		m = length(items)
		#sapply(KS,function(x) items %in% x)
		mat = matrix(sapply(as.list(KS),function(x){
				 						tmp = rep(0,m)
										#ind = which(letters %in% x)# alphabetische reihenfolge
				 						ind = which(items %in% x)
										tmp[ind] = 1
										return(tmp)
										}
				),ncol = n, nrow = length(items))
		if(length(items) > 1){ dimnames(mat)[[1]] = items }
		
		if (return.dataframe){
			dat=data.frame(t(mat),Freq)
			names(dat)[1:dim(mat)[1]]=rownames(mat)
			return(dat)
		}
		return( mat	)
	}
		print("KS is neither matrix nor set! No conversion will be applied.")
		return(KS)
}

#_____________________________________________________________________#
#-----------------------DATA FRAME -> GSET----------------------------#

as.gset.data.frame=function(x){
	if (!("Freq" %in% names(x))){
		x=data.frame(ftable(x))
		ind=which(names(x)!="Freq")
	}else{
		ind=which(names(x)!="Freq")
		x[,ind] = data.frame(sapply(x[,ind],as.factor))
	}
	x=x[x$"Freq"!=0,]
	
	ord1 = do.call(order,c(x,decreasing=T))
	x[,ind]=sapply(x[,ind],as.integer)-1
	x=x[ord1,]
	dims = dim(x)
	if (dims[2] == 2){
		ord2 = order(rowSums(matrix(x[,ind],nrow=dims[1],ncol=dims[2]-1)))
		x=x[ord2,]
		mat = matrix(x[,ind],nrow=dims[1],ncol=dims[2]-1)
		colnames(mat) = colnames(x)[colnames(x)!="Freq"]
	}
	else{
		ord2 = order(rowSums(x[,ind]))
		x=x[ord2,]
		mat = x[,ind]
	}
	return(gset(convert(t(mat)),memberships=x$"Freq"))
}

#_____________________________________________________________________#
#------------------- RELABELLING DATA IN A RELATION ------------------#

relabel.relation = function(er,domain = NULL, data = NULL){
	stopifnot( !is.null(domain) | !is.null(data) )
	stopifnot(inherits(er,"relation"))
	er = endorelation(incidence=relation_incidence(er))
	if(!is.null(data)){
		stopifnot(inherits(data,"data.frame"))
		domain = names(data)
	}
	er$".Data"$domain[1]$X = domain
	er$".Data"$domain[2]$X = domain
	return(er)
}
