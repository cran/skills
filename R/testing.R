#_____________________________________________________________________#
#------------------------SYMMETRIC DIFFERENCE-------------------------#

hammingDist=function(s1,s2){
	if(is.vector(s1)&is.vector(s2)){
		m=as.list(convert(cbind(s1,s2)))
		s1=m[[1]]
		s2=m[[2]]
	}
	return(length(set_symdiff(s1,s2)))
}

#_____________________________________________________________________#
#-------------------------RELATION MAPPINGS---------------------------#

#necessity operator
NecOp=function(Rel,X){
	out=set()
	for (i in 1:(dim(Rel)[2])){
		if (set_is_subset(as.set(unlist(convert(as.matrix(Rel[,i])))),as.set(as.character(unlist(X)))))
			out=as.set(out)+set(colnames(Rel)[i])
	}
	out=as.set(out)
	return (as.set(out))
}

#possibility operator
PosOp=function(Rel,X){
	out=set()
	for (i in 1:(dim(Rel)[2])){
		if (set_intersection(as.set(unlist(convert(as.matrix(Rel[,i])))),as.set(as.character(unlist(X))))!=set())
			out=as.set(out)+set(colnames(Rel)[i])
	}
	out=as.set(out)
	return (as.set(out))
}

#sufficiency operator
SufOp=function(Rel,X){
	out=set()
	for (i in 1:(dim(Rel)[2])){
		if (set_is_subset(as.set(as.character(unlist(X))),as.set(unlist(convert(as.matrix(Rel[,i]))))))
			out=as.set(out)+set(colnames(Rel)[i])
	}
	return (as.set(out))
}

#dual sufficiency operator
DSufOp=function(Rel,X){
	out=set()
	for (i in 1:(dim(Rel)[2])){
		if (!set_is_equal(set_union(as.set(unlist(convert(as.matrix(Rel[,i])))),as.set(as.character(unlist(X)))),rownames(Rel)))
			out=as.set(out)+set(colnames(Rel)[i])
	}
	return (as.set(out))
}

#_____________________________________________________________________#
#----------------------LOWER BOUND / UPPER BOUND----------------------#

lb_eKS=function(eKS,tKS,model){
	KS=asCharSet(as.set(eKS))
	tKS=asCharSet(tKS)
	SA=SkillAss(tKS,model=model)
	lb=vector(length=length(KS),mode="list")
	attr(lb,"model")=model
	if (model=="disjunctive"){
		for (i in 1:length(KS)){
			newset=PosOp(t(SA),NecOp(SA,as.list(KS)[[i]]))
			lb[[i]]=newset
		}
		return(lb)
	}else{
		if (model=="conjunctive"){
			for (i in 1:length(KS)){
				if (set_contains_element(tKS,as.list(KS)[[i]])){
					lb[[i]]=as.list(KS)[[i]]
				}else{
					partset=set(set())
					for (j in 1:length(tKS)){
						if (set_is_subset(as.list(tKS)[[j]],as.list(KS)[[i]])){
							partset=as.set(partset)+set(as.list(tKS)[[j]])
						}
					}
					partset=as.set(partset)
					M=convert(partset-set(set()))
					MAT=t(M) %*%M
					diag(MAT)=-1
					comp=apply(M,2,sum)-apply(MAT,2,max)
					if (sum(comp!=0)==1){
						lb[[i]]=convert(as.matrix(M[,which(comp!=0)]))
					}else{
						intsec=apply(M[,which(comp!=0)],1,prod)
						lb[[i]]=as.set(unlist(convert(as.matrix(intsec))))
					}
					
					lb[[i]]=as.set(unlist(lb[[i]]))
				}
			}
			return(lb)
		}else{
			stop("model is neither disjunctive nor conjunctive")
		}
	}
}

ub_eKS=function(eKS,tKS,model){
	KS=asCharSet(as.set(eKS))
	tKS=asCharSet(tKS)
	SA=SkillAss(tKS,model=model)
	ub=vector(length=length(KS),mode="list")
	attr(ub,"model")=model
	if (model=="disjunctive"){
		for (i in 1:length(KS)){
				if (set_contains_element(tKS,as.list(KS)[[i]])){
					ub[[i]]=as.list(KS)[[i]]
				}else{
					partset=set(set())
					for (j in 1:length(tKS)){
						if (set_is_subset(as.list(KS)[[i]],as.list(tKS)[[j]])){
							partset=as.set(partset)+set(as.list(tKS)[[j]])
						}
					}
					if (length(partset-set(set()))==1){
						ub[[i]]=as.set(unlist(partset-set(set())))
					}else{
						M=convert(partset-set(set()))
						Mdual=dualKS(M)
						MAT=t(Mdual) %*%Mdual
						diag(MAT)=-1
						comp=apply(Mdual,2,sum)-apply(MAT,2,max)
						if (sum(comp!=0)==1){
							ub[[i]]=convert(as.matrix(M[,which(comp!=0)]))
						}else{
							intsec=apply(M[,which(comp!=0)],1,prod)
							ub[[i]]=as.set(unlist(convert(as.matrix(intsec))))
						}
							ub[[i]]=as.set(unlist(ub[[i]]))
					}	
				}
			}
		return(ub)
	}else{
		if (model=="conjunctive"){
			for (i in 1:length(KS)){
				newset=NecOp(t(SA),PosOp(SA,as.list(KS)[[i]]))
			ub[[i]]=newset
			}
			return(ub)
		}else{
			stop("model is neither disjunctive nor conjunctive")
		}
	}
}

#_____________________________________________________________________#
#----------------------------TEST PARAMETERS--------------------------#

#consistency index
cind=function(eKS,tKS){
	eKS=asCharSet(as.set(eKS))
	tKS=asCharSet(tKS)
	if (typeof(tKS)=="double")
		tKS=convert(tKS)
	#if (set_union(as.set(unlist(as.set(eKS))))<=set_union(as.set(unlist(tKS)))){
		ind=sum(eKS %in% tKS)/length(eKS)
		return(ind)
	#}else{
	#	print("Empirical knowledge structure consists of more/other items than theoretical knowledge structure")
	#}
}

#consistency index weighted by persons
w_cind=function(eKS,tKS){
	eKSnew=asCharSet(as.set(eKS))
	tKS=asCharSet(tKS)
	if (typeof(tKS)=="double")
		tKS=convert(tKS)
	#if (set_union(as.set(unlist(eKSnew)))<=set_union(as.set(unlist(tKS)))){
		ind=sum((eKSnew %in% tKS)*(attributes(eKS)$memberships))/length(eKS)
		return(ind)
	#}else{
	#	print("Empirical knowledge structure consists of more/other items than theoretical knowledge structure")
	#}
}

#_____________________________________________________________________#
#permutation test: random skill assignment of problems to skill sets?

#compute a matrix of all permutations of a given vector (only advisible for length <= 7)
swap=function(v,i,j){
	tmp=v[i]
	v[i]=v[j]
	v[j]=tmp
	return(v)
}

rotateLeft=function(v,beg,last){
	if (beg==1){
		if (last==length(v)){
			v=v[c(2:last,1)]
		}else{
			v=v[c(2:last,1,(last+1):length(v))]
		}
	}else{
		if (last==length(v)){
			v=v[c(1:(beg-1),(beg+1):last,beg)]
		}else{
			v=v[c(1:(beg-1),(beg+1):last,beg,(last+1):length(v))]
		}
	}
	return(v)
}

permutations=function (v, beg = 1, last = length(v), env = .GlobalEnv){
    if (beg <= last) {
        i = last - 1
        while (i >= beg) {
            for (j in (i + 1):last) {
                v = swap(v, i, j)
                env$mat = cbind(env$mat, v)
                v = Recall(v, i + 1, last, env)
            }
            v = rotateLeft(v, i, last)
            i = i - 1
        }
    }
    return(v)
}

permute=function (v, beg = 1, last = length(v)){
    if (is.vector(v)) {
        e1 = new.env()
        e1$mat = v
        permutations(v, beg, last, e1)
        return(e1$mat)
    }
    else {
        print("v is no vector -> v is returned unchanged")
        return(v)
    }
}

#compute test parameter
############## FUNKTIONIERT NUR FÜR CONJUNCTIVE / DISJUNCTIVE #####################

permtest=function(eKS,tKS,model){
	eKS=asCharSet(eKS)
	tKS=asCharSet(tKS)
	SA=SkillAss(tKS,model=model)
	if (typeof(tKS)=="list"){
		tKSmat=convert(tKS)
	}else{
		tKSmat=tKS
	}
	Gsigma=apply(permute(rownames(SA)),2,function(x){
												rownames(tKSmat)=x
												tKSsigma=convert(tKSmat)
												w_cind(eKS,tKSsigma)
												return(w_cind(eKS,tKSsigma))
												}
			)
	pvalue=sum((Gsigma>=Gsigma[1]))/factorial(dim(SA)[1])
	return(pvalue)
}

#_____________________________________________________________________#
#expected value
ExpectedValue=function(eKS,tKS,model){
	eKS=asCharSet(eKS)
	tKS=asCharSet(tKS)
	SA=SkillAss(tKS,model=model)
	if (typeof(tKS)=="list"){
		tKSmat=convert(tKS)
	}else{
		tKSmat=tKS
	}
	Gsigma=apply(permute(rownames(SA)),2,function(x){
												rownames(tKSmat)=x
												tKSsigma=convert(tKSmat)
												w_cind(eKS,tKSsigma)
												return(w_cind(eKS,tKSsigma))
												}
			)
	expv=sum(Gsigma)/factorial(dim(SA)[1])
	return(expv)
}


#adjusted consistency measure
acm=function(eKS,tKS,model){
	kappa=(w_cind(eKS,tKS)-ExpectedValue(eKS,tKS,model))/(1-ExpectedValue(eKS,tKS,model))
	return(kappa)
}
#_____________________________________________________________________#
#----------------------------ITEM CONSISTENCY-------------------------#
############## FUNKTIONIERT NUR FÜR CONJUNCTIVE / DISJUNCTIVE #####################
item_cind=function(eKS,tKS,model){
	eKS=asCharSet(eKS)
	tKS=asCharSet(tKS)
	SA=SkillAss(tKS,model=model)
	tKSmat=convert(tKS)
	if (model=="skill multimap"){
		print("No item consistency will be computed because problem function is neither disjunctive nor conjunctive")
	}else{
		eKSmat=rbind(convert(as.set(eKS)),attributes(eKS)$memberships)
		rownames(eKSmat)=c(rownames(eKSmat)[1:(length(rownames(eKSmat))-1)],"Freq")
		gammaqs=sapply((1:(length(rownames(eKSmat))-1)),function(x){
							if (x==1){
								eKSmatq=t(subKS(as.data.frame(t(eKSmat)),c(2:(length(rownames(eKSmat))-1)),keep.zero=FALSE,allfactor=FALSE))
								rownames(eKSmatq)=c(rownames(eKSmat)[2:(length(rownames(eKSmat))-1)],"Freq")
								tKSmatq=tKSmat[c(2:(length(rownames(eKSmat))-1)),]
							}else{
								if (x==(length(rownames(eKSmat))-1)){
									eKSmatq=t(subKS(as.data.frame(t(eKSmat)),c(1:(length(rownames(eKSmat))-2)),keep.zero=FALSE,allfactor=FALSE))
									rownames(eKSmatq)=c(rownames(eKSmat)[1:(length(rownames(eKSmat))-2)],"Freq")
									tKSmatq=tKSmat[c(1:(length(rownames(eKSmat))-2)),]
								}else{
									eKSmatq=t(subKS(as.data.frame(t(eKSmat)),c(1:(x-1),(x+1):(length(rownames(eKSmat))-1)),keep.zero=FALSE,allfactor=FALSE))
									rownames(eKSmatq)=c(rownames(eKSmat)[c(1:(x-1),(x+1):(length(rownames(eKSmat))-1))],"Freq")
									tKSmatq=tKSmat[c(1:(x-1),(x+1):(length(rownames(eKSmat))-1)),]
								}
							}
							eKSq=as.gset(as.data.frame(t(eKSmatq)))
							tKSq=convert(tKSmatq)
							return(w_cind(eKSq,tKSq))
		})
		mat=c(w_cind(eKS,tKS),gammaqs)
		names(mat)=c("none",rownames(eKSmat)[1:(length(rownames(eKSmat))-1)])
		return(mat)
	}
}

#_____________________________________________________________________#
#---------------------------SKILL CONSISTENCY-------------------------#
########### FUNKTIONIERT NUR FÜR CONJUNCTIVE / DISJUNCTIVE ############
skill_cind=function(eKS,tKS,model,allsubs=FALSE){
	eKS=asCharSet(eKS)
	tKS=asCharSet(tKS)
	SA=SkillAss(tKS,model=model)
	if (allsubs == TRUE){
		pot=t(expand.grid( lapply(data.frame(replicate(dim(SA)[2],c(0,1))),I)  ))
		ord=order(colSums(pot))
		pot=pot[,ord]
		rownames(pot)=colnames(SA)
		tKSmat=convert(tKS)
		if (model=="skill multimap"){
			print("No item consistency will be computed because problem function is neither disjunctive nor conjunctive")
		}else{
			gammaT=sapply(1:(2^(dim(SA)[2])-1),function(x){
							SAT=SA[,which((pot[,x])==0)]
							tKST=convert(probfun(SAT,model))
							return(w_cind(eKS,tKST))
							})
			names(gammaT)=c("none",sapply(2:((2^(dim(SA)[2]))-1),function(x) do.call(paste,as.list(rownames(pot)[as.logical(pot[,x])]))))
			return(gammaT)
		}
	}
	else{
		size=ncol(SA)
		pot=matrix(0,ncol=size,nrow=size)
		diag(pot)=1
		rownames(pot)=colnames(SA)
		tKSmat=convert(tKS)
		if (model=="skill multimap"){
			print("No item consistency will be computed because problem function is neither disjunctive nor conjunctive")
		}else{
			gammaT=sapply(1:size,function(x){
							SAT=SA[,which((pot[,x])==0)]
							tKST=convert(probfun(SAT,model))
							return(w_cind(eKS,tKST))
							})
			gammaT=c(w_cind(eKS,tKS),gammaT)
			names(gammaT)[1]="none"
			names(gammaT)[2:length(gammaT)]=colnames(SA)
			return(gammaT)
		}
	}
}

#_____________________________________________________________________#
#-----------------------DISTANCE TO BOUNDARIES------------------------#

lowerDist=function(eKS,tKS,model){
	listeks = as.list(asCharSet(eKS))
	lb = lb_eKS(asCharSet(eKS),asCharSet(tKS),model)
	ham=lapply(1:length(as.set(eKS)),function(x) hammingDist(listeks[[x]],lb[[x]]))
	attr(ham,"model")=attr(lb,"model")
	return(ham)
}

upperDist=function(eKS,tKS,model){
	listeks = as.list(asCharSet(eKS))
	ub = ub_eKS(asCharSet(eKS),asCharSet(tKS),model)
	ham=lapply(1:length(as.set(eKS)),function(x) hammingDist(listeks[[x]],ub[[x]]))
	attr(ham,"model")=attr(ub,"model")
	return(ham)
}

lowerRep=function(eKS,tKS,model){
	REP=1-(do.call(sum,lowerDist(eKS,tKS,model)))/(length(eKS)*length(as.set(unlist(eKS))))
	attr(REP,"model")=attr(lowerDist(eKS,tKS,model),"model")
	return(REP)
}

upperRep=function(eKS,tKS,model){
	REP=1-(do.call(sum,upperDist(eKS,tKS,model)))/(length(eKS)*length(as.set(unlist(eKS))))
	attr(REP,"model")=attr(upperDist(eKS,tKS,model),"model")
	return(REP)
}
