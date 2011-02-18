#_____________________________________________________________________#
#-------------REPLACE ALL ITEM NAMES WITH CHARACTER NAMES-------------#

asCharSet=function(KS){
	if (!inherits(KS,"set") & !inherits(KS,"gset")){
		stop("KS must be of class set or gset")
	}
	
	if (!is.null(attributes(KS)$memberships)){
		memb=attributes(KS)$memberships
		KS=as.set(lapply(KS,function(x) as.set(as.character(x))))
		KS=gset(support=KS,memberships=memb)
	}else{
		KS=as.set(lapply(KS,function(x) as.set(as.character(x))))
	}
	return(KS)
}

#_____________________________________________________________________#
#----------------------DUAL KNOWLEDGE STRUCTURE-----------------------#

dualKS=function(KS){
	stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
#	stopifnot(is.KS(KS))
	S=KS
	if (inherits(S,"set") | inherits(S,"gset")){
		KS=convert(S)
	}
	dualmat=1-KS
	#Dualmat=dualmat
	if (inherits(S,"set") | inherits(S,"gset")){
		dualmat=convert(dualmat)#D
	}
	return (dualmat)#D
}

#_____________________________________________________________________#
#--------------SKILL ASSIGNMENT -> KNOWLEDGE STRUCTURE----------------#

probfun=function(SA,model="disjunctive",alldata=FALSE){
	if (is.vector(SA)){
		KS=cbind(0,SA)
	}else{
		SA=SA[which(rowSums(SA)!=0),which(colSums(SA)!=0)]
		pot=t(expand.grid( lapply(data.frame(replicate(dim(SA)[2],c(0,1))),I)  ))
		if (!is.null(colnames(SA))){
			rownames(pot) = colnames(SA)
		}
		else{
			rownames(pot) = as.vector(sapply(c("",1:(floor(dim(SA)[2]/26))),function(x) paste(letters,x,sep="")))[1:dim(SA)[2]]
		}
		mult=SA %*% pot
		if( model=="disjunctive" ){
			if (dim(SA)[1] != length(unique(rownames(SA)))){
				stop("More than one set of skills assigned to one item -> disjunctive model not suitable")
				return()
			}else{
				KS=t(unique(data.frame(t(mult>0))))
				row.names(KS)=row.names(SA)
				KS=as.KS(KS+0,type="matrix")
				attr(KS,"model")=model
				KSall=t(data.frame(t(mult>0)))
				row.names(KSall)=row.names(SA)
				KSall=KSall+0
				if (alldata){
					attr(KS,"KSall")=as.KS(KSall,type="matrix")
					attr(KS,"pot")=pot
				}
				return(KS)
			}
		}else{
			if( model=="conjunctive" ){
				if (dim(SA)[1] != length(unique(rownames(SA)))){
					stop("More than one set of skills assigned to one item -> conjunctive model not suitable")
					return()
				}else{
					matsum=apply(SA,1,sum)
					shift=matsum-mult==0
					KS=t(unique(data.frame(t(shift))))
					row.names(KS)=row.names(shift)
					KS=as.KS(KS+0,type="matrix")
					attr(KS,"model")=model
					KSall=t(data.frame(t(shift)))
					row.names(KSall)=row.names(shift)
					KSall=KSall+0
					if (alldata){
						attr(KS,"KSall")=as.KS(KSall,type="matrix")
						attr(KS,"pot")=pot
					}
					return(KS)
				}
			}else{
				if( model=="skill multimap" ){
					matsum=apply(SA,1,sum)
					shift=matsum-mult
					shift2=apply(shift,2,function(x) tapply(x,row.names(shift),min))==0
					KS=t(unique(data.frame(t(shift2))))
					KS=as.KS(KS+0,type="matrix")
					attr(KS,"model")=model
					row.names(KS)=row.names(shift2)
					KSall=t(data.frame(t(shift2)))
					row.names(KSall)=row.names(shift2)
					KSall=KSall+0
					if (alldata){
						attr(KS,"KSall")=as.KS(KSall,type="matrix")
						attr(KS,"pot")=pot
					}
					return(KS)
				}else{
					stop("model neither disjunctive nor conjunctive nor skill multimap")
				}
			}
		}
	}	
}

#_____________________________________________________________________#
#-----------------MINIMALIZATION OF SKILL ASSIGNMENTS-----------------#

minSA=function(SA,model){
	stopifnot(inherits(SA,"matrix"))
	wholeKS=probfun(SA,model)
	smaller=sapply(1:dim(SA)[2],function(x) {
				SA1=SA[,which(1:dim(SA)[2] != x)]
				KS1=probfun(SA1,model)
				return(!(dim(wholeKS)[2]>dim(KS1)[2]))
	})
	if (any(smaller)){
		smallerSA=SA[,which(1:dim(SA)[2] != min(which(smaller)))]
		SA=Recall(smallerSA,model)
	}
	return(SA)			
}

#_____________________________________________________________________#
#---------------------------KNOWLEDGE SPACE---------------------------#

#test:closed under union
is.KSpace=function(KS,domain=NULL) {
	#stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	stopifnot(is.KS(KS))
	if (inherits(KS,"matrix") ){
		KS=convert(KS)
	}
	if(is.null(domain)){
		domain=KSdomain(KS)
	}
	cl=closure(KS,"union")
	missings=cl-KS
	ret=missings=={}
	if (!ret){
		attr(ret,"missings")=missings
	}
	return(ret)
}

#join-irreducible elements
#__________________________________________________________________________
JoinIrrEl=function(KS){
	stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
#	stopifnot(is.KS(KS))
	if (inherits(KS,"gset")){
		KS=as.set(KS)
		type = "set"
	}
	if (inherits(KS,"matrix") ){
		type="matrix"
		KS=convert(KS)
	}
	ret = reduction(KS,"union")-set(set())
	if (type =="matrix"){
		ret=convert(ret)
		}
	return(ret)
}

#skill assignment (disjunctive) D&G
KSkillAss=function(KS){
	#stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	stopifnot(is.KS(KS))
	if (inherits(KS,"matrix") )
		KS=convert(KS)
	if (is.KSpace(KS)){
		model="disjunctive"
		SA=convert(JoinIrrEl(KS)-set(set()))
		#colnames(SA)=letters[1:dim(SA)[2]]
		colnames(SA) = as.vector(sapply(c("",1:(floor(dim(SA)[2]/26))),function(x) paste(letters,x,sep="")))[1:dim(SA)[2]]
		return (list(SA=SA,model=model))
	}else{
		print("KS ist not a knowledge space as required!")
		return(set())	# default rueckgabetyp set???
	}
}

#non-minimal-skill assignment D&F
nonminKSkillAss=function(KS){
	#stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	stopifnot(is.KS(KS))
	if (inherits(KS,"set") | inherits(KS,"gset"))
		KS=convert(KS)
	if (is.KSpace(KS)){
		model="disjunctive"
		SA=KS[,which(colSums(KS)!=0)]
		#colnames(SA)=letters[1:dim(SA)[2]]
		colnames(SA) = as.vector(sapply(c("",1:(floor(dim(SA)[2]/26))),function(x) paste(letters,x,sep="")))[1:dim(SA)[2]]
        return (list(SA=SA,model=model))
	}else{
		print("KS ist not a knowledge space as required!")
		return(set())	# default rueckgabetyp set???
	}
}	

#_____________________________________________________________________#
#--------------------------CLOSURE SPACE------------------------------#

#test:closed under intersection
is.CSpace=function(KS,domain=NULL) {
	#stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	stopifnot(is.KS(KS))
	if (inherits(KS,"matrix") ){
		KS=convert(KS)
	}
	if(is.null(domain)){
		domain=KSdomain(KS)
	}
	cl=closure(KS,"intersection")
	missings=cl-KS
	ret=missings=={}
	if (!ret){
		attr(ret,"missings")=missings
	}
	return(ret)
}

#meet-irreducible elements
#___________________________________________________________________
MeetIrrEl=function(KS){
	stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
#	stopifnot(is.KS(KS))
	if (inherits(KS,"gset")){
		KS=as.set(KS)
		type = "set"
	}
	if (inherits(KS,"matrix") ){
		KS=convert(KS)
		type="matrix"
	}
	ret=reduction(KS,"intersection")-set(KSdomain(KS))
	if (type=="matrix"){
		ret=convert(ret)
	}
	return(ret)
}

#skill assignment (conjunctive) D&G
ClSkillAss=function(KS){
	#stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	stopifnot(is.KS(KS))
	if (inherits(KS,"matrix") )
		KS=convert(KS)
	if (is.CSpace(KS)){
		model="conjunctive"
		SA=1-convert(MeetIrrEl(KS))
		SA=SA[,which(colSums(SA)!=0)]
		#colnames(SA)=letters[1:dim(SA)[2]]
		colnames(SA) = as.vector(sapply(c("",1:(floor(dim(SA)[2]/26))),function(x) paste(letters,x,sep="")))[1:dim(SA)[2]]
        return (list(SA=SA,model=model))
	}else{
		print("KS ist not a closure space as required!")
		return(set())	# default rueckgabetyp set???
	}
}

#non-minimal-skill assignment D&F
nonminClSkillAss=function(KS){
	#stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	stopifnot(is.KS(KS))
	if (inherits(KS,"set") | inherits(KS,"gset"))
		KS=convert(KS)
	if (is.CSpace(KS)){
		model="conjunctive"
		SA=dualKS(KS)[,which(colSums(dualKS(KS))!=0)]
		#colnames(SA)=letters[1:dim(SA)[2]]
		colnames(SA) = as.vector(sapply(c("",1:(floor(dim(SA)[2]/26))),function(x) paste(letters,x,sep="")))[1:dim(SA)[2]]
        return (list(SA=SA,model=model))
	}else{
		print("KS ist not a closure space as required!")
		return(set())	# default rueckgabetyp set???
	}
}

#_____________________________________________________________________#
#----------------------TEST FOR QUASI-ORDINAL SPACE--------------------#

is.quordSpace=function(KS,domain=NULL){
	stopifnot(is.KS(KS))
	if (inherits(KS,"matrix") ){
		KS=convert(KS)
	}
	if(is.null(domain)){
		domain=KSdomain(KS)
	}
	ret=(is.CSpace(KS,domain) & is.KSpace(KS,domain))
	if (!ret){
		clos = closure(closure(KS,"intersection"),"union")
		missings = clos-KS
		attr(ret,"missings")=missings
	}
	return(ret)	
}

#_____________________________________________________________________#
#------------------MINIMIZE KS TO INCOMPARABLE STATES-----------------#

incomparable=function(KS){
	stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	#stopifnot(is.KS(KS))
	if (inherits(KS,"set") | inherits(KS,"gset")){
		M=convert(KS-set(set()))
	}else{
		if (inherits(KS,"matrix") ){
			M=KS[,which(colSums(KS)!=0)]
		}
	}
	MAT=t(M) %*% M
	diag(MAT)=-1
	comp=apply(M,2,sum)-apply(MAT,2,max)
	if (prod(comp)==0){
		incomp=F
	}else{
		incomp=T
	}
	return(incomp)
}

#_____________________________________________________________________#
#------------------DÜNTSCH & GEDIGA: SKILL MULTIMAP-------------------#

SkillAssDGgen=function(KS){
	#stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	stopifnot(is.KS(KS))
	if (inherits(KS,"set") | inherits(KS,"gset")){
		KS=convert(KS)
	}
	Q=rownames(KS)
	K=KS[,which(colSums(KS)!=0)]
	
	mats=lapply(1:dim(K)[2],function(x){
				Nrow=colSums(K)[x]
				r=rep(0,dim(K)[2])
				r[x]=1
				mat=matrix(rep(r,Nrow),nrow=Nrow,byrow=T)
				rownames(mat)=rownames(K)[which(K[,x]==1)]
				return(mat)
	})
	SA=do.call(rbind,mats)
	
	len=dim(K)[2]
	pot=t(expand.grid( lapply(data.frame(replicate(len,c(0,1))),I)  ))
	colnames(pot)=1:2^len
	myunion=apply(pot,2,function(x) {
				v=rep(1,dim(K)[1])
				v[which(rowSums(as.matrix(K[,which(x==1)]))==0)]=0
				return(v)
	})
	if (dim(KS)[1] %in% colSums(KS)){
		myunionRest=myunion[,which(colSums(myunion)!=dim(myunion)[1])]
	}else{
		myunionRest=myunion
	}
	lis=lapply(data.frame(myunionRest),I)
	Kis=lapply(data.frame(K),I)
	cond2bmat=myunionRest[,which((lis %in% Kis)==F)]
	cond2b=colnames(cond2bmat)      #states, delineated by an element of the power set, but not contained in KS
	
	if (is.vector(cond2bmat)){
		#print("Condition 2b is never satisfied -> no competencies necessary")
	}else{
		cond2a=sapply(as.integer(colnames(cond2bmat)[2:dim(cond2bmat)[2]]),function(x) {
						return(incomparable(as.matrix(K[,which(pot[,x]==1)])))
		})        #are the states, which are associated to the above computed elements of the power set, incomparable
		cond2amat=cbind(cond2bmat[,1],cond2bmat[,1+which(cond2a)])
		colnames(cond2amat)=c(1,colnames(cond2bmat)[1+which(cond2a)])
		rownames(cond2amat)=rownames(KS)
#		if (dim(KS)[1] %in% colSums(KS)){
#			if (dim(cond2amat)[2]==2){
#				x=pot[,as.integer(colnames(cond2amat)[2])]
#				Nrow=dim(cond2amat)[1]-sum(cond2amat[,2])
#				mat=matrix(rep(x,Nrow),nrow=Nrow,byrow=T)
#				rownames(mat)=rownames(K)[which(cond2amat[,2]==0)]
#				news=mat
#			}else{
#				ne=lapply(2:dim(cond2amat)[2],function(x){
#						v=cond2amat[,x]
#						y=pot[,as.integer(colnames(cond2amat)[x])]
#						Nrow=dim(cond2amat)[1]-sum(v)
#						mat=matrix(rep(y,Nrow),nrow=Nrow,byrow=T)
#						rownames(mat)=rownames(K)[which(v==0)]
#						return(mat)
#				})
#				news=do.call(rbind,ne)
#			}
#			SA=rbind(SA,news)
#		}else{
#			print("State Q does not exist!")
#			print("Skill assignment is not inverse to problem function!")
			MAT=t(cond2amat[,2:dim(cond2amat)[2]]) %*% K
			rownames(MAT)=colnames(cond2amat)[2:dim(cond2amat)[2]]
			if (dim(cond2amat)[2]==2){
				if (sum(MAT>=sum(cond2amat[,2]))==0){
					print("union is not contained in any other state A")
				}else{
					A=matrix(K[,as.integer(which(MAT==sum(cond2amat[,2])))],nrow=nrow(K))
					rownames(A)=rownames(K)
					uni=ifelse(rowSums(A)>0,1,0)
					Nrow=sum(uni-cond2amat[,2])
					x=pot[,as.integer(colnames(cond2amat)[2])]
					mat=matrix(rep(x,Nrow),nrow=Nrow,byrow=T)
					rownames(mat)=rownames(K)[which(uni != cond2amat[,2])]
					news=mat
				}
			}else{
				ne=lapply(2:dim(cond2amat)[2],function(x){
						v=cond2amat[,x]
						y=pot[,as.integer(colnames(cond2amat)[x])]
						A=matrix(K[,which(MAT[x-1,]==sum(cond2amat[,x]))],nrow=nrow(K))
						rownames(A)=rownames(K)
						if (dim(A)[2]==1){
							name=A-v
						}else{
							p = t(A) %*% A - colSums(A)
							diag(p) = 1
							name = rowSums(matrix(A[,which(colSums(p==0)==0)],nrow=nrow(A)))-v
						}
						Nrow=sum(name!=0)
						mat=matrix(rep(y,Nrow),nrow=Nrow,byrow=TRUE)
						rownames(mat)=rownames(A)[name!=0]
#						if (length(which(MAT[x-1,]==sum(cond2amat[,x])))==0){
#							mat=NULL
#						}else{
#							if (length(which(MAT[x-1,]==sum(cond2amat[,x])))==1){
#								uni=A
#							}else{
#								uni=ifelse(rowSums(A)>0,1,0)
#							}
#							Nrow=sum(uni-cond2amat[,x])
#							mat=matrix(rep(y,Nrow),nrow=Nrow,byrow=T)
#							rownames(mat)=rownames(K)[which(uni != cond2amat[,x])]
#						}
						return(mat)
				})
				news=do.call(rbind,ne)
			}
			SA=rbind(SA,news)
#		}
	}
	SA=SA[order(rownames(SA)),]
	SA=SA[,1:(ncol(SA)-1)]   #the state Q is not necessary as skill
	SA=SA[which(rowSums(SA)!=0),]
	#colnames(SA)=letters[1:dim(SA)[2]]
	colnames(SA) = as.vector(sapply(c("",1:(floor(dim(SA)[2]/26))),function(x) paste(letters,x,sep="")))[1:dim(SA)[2]]
    model="skill multimap"
	return(list(SA=SA,model=model))
}
#_____________________________________________________________________#
#-----------------DOIGNON & FALMAGNE: SKILL ASSIGNMENT----------------#
SkillAssDFgen=function(KS){
	#stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	stopifnot(is.KS(KS))
	if (inherits(KS,"set") | inherits(KS,"gset"))
		KS=convert(KS)
	SA=matrix(1,nrow=sum(KS),ncol=dim(KS)[2])
	#dimnames(KS)[[1]]=1:dim(KS)[1]
	k=1
#	for (i in 1:dim(KS)[1]){
#		for (j in 1:dim(KS)[2]){
#			if (KS[i,j]==1){
#				SA[k,j]=0
#				k=k+1
#			}
#		}	
#	}

	m = sum(KS)
	cs = cumsum(t(KS))*t(KS)
	SA = apply(cs,1,function(x){
		ret = rep(1,m)
		ret[x] = 0
		return(ret)
	})

	rownames(SA)=rep(row.names(KS),apply(KS,1,sum))
	#colnames(SA)=letters[1:dim(SA)[2]]
	colnames(SA) = as.vector(sapply(c("",1:(floor(dim(SA)[2]/26))),function(x) paste(letters,x,sep="")))[1:dim(SA)[2]]
    model="skill multimap"
	return (list(SA=SA,model=model))
}

#_____________________________________________________________________#
#-------------THE ALLMIGTHY SKILL ASSIGNMENT FUNCTION-----------------#

SkillAss=function(KS,method="DG",minimize=F,model="disjunctive"){
	#stopifnot( inherits(KS,"set") | inherits(KS,"gset") | inherits(KS,"matrix") )
	stopifnot(is.KS(KS))
	if (model=="disjunctive"){
		if (method=="DG"){
			if (is.KSpace(KS)){
				ret=KSkillAss(KS)$SA
				attr(ret,"method")="DG"
				attr(ret,"model")="disjunctive"
				return(ret)
			}else{
				stop("disjunctive model requires a knowledge space")				
				#ret=SkillAssDGgen(KS)$SA
				#attr(ret,"method")="DG"
				#attr(ret,"model")="skill multimap"
				#print("KS not closed under union -> skill multimap")
				#return(ret)
			}
		}
		if (method=="DF"){
			if (is.KSpace(KS)){
				if (minimize){
					ret=minSA(nonminKSkillAss(KS)$SA,"disjunctive")
				}else{
					ret=nonminKSkillAss(KS)$SA
				}
				attr(ret,"method")="DF"
				attr(ret,"model")="disjunctive"
				return(ret)
			}else{
				stop("disjunctive model requires a knowledge space")
				#ret=SkillAssDFgen(KS)$SA
				#attr(ret,"method")="DF"
				#attr(ret,"model")="skill multimap"
				#print("KS not closed under union -> skill multimap")
				#return(ret)
			}
		}
	}
	if (model=="conjunctive"){
		if (method=="DG"){
			if (is.CSpace(KS)){
				ret=ClSkillAss(KS)$SA
				attr(ret,"method")="DG"
				attr(ret,"model")="conjunctive"
				return(ret)
			}else{
				stop("conjunctive model requires a closure space")
				#SkillAssDGgen(KS)$SA
				#ret=SkillAssDGgen(KS)$SA
				#attr(ret,"method")="DG"
				#attr(ret,"model")="skill multimap"
				#print("KS not closed under intersection -> skill multimap")
				#return(ret)
			}
		}
		if (method=="DF"){
			if (is.CSpace(KS)){
				if (minimize){
					ret=minSA(nonminClSkillAss(KS)$SA,"conjunctive")
				}else{
					ret=nonminClSkillAss(KS)$SA
				}
				attr(ret,"method")="DF"
				attr(ret,"model")="conjunctive"
				return(ret)
			}else{
				stop("conjunctive model requires a closure space")
				#ret=SkillAssDFgen(KS)$SA
				#attr(ret,"method")="DF"
				#attr(ret,"model")="skill multimap"
				#print("KS not closed under intersection -> skill multimap")
				#return(ret)
			}
		}
	}
	if (model=="skill multimap"){
		if (method=="DG"){
			ret=SkillAssDGgen(KS)$SA
			attr(ret,"method")="DG"
			attr(ret,"model")="skill multimap"
			return(ret)
		}
		if (method=="DF"){
			ret=SkillAssDFgen(KS)$SA
			attr(ret,"method")="DF"
			attr(ret,"model")="skill multimap"
			return(ret)
		}
	}
	stop("inaccurate model or method specification")
}








