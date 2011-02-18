plot.KS = function (x, y = NULL, ord = NULL, max.scale = "rel", lab.cex = 1.5, 
    border = "black", bg = c("grey70", "skyblue", "tomato"), 
    fill = c("tomato", "dodgerblue"), shape = "rect", ...) 
{
    KS = x
    KS2 = y
    if (!inherits(KS, "KS") & !(is.matrix(KS))) {
        stop(print("KS is not of class KS or matrix"))
    }
    if (length(bg) < 3) {
        bg = c(bg[1], "skyblue", "tomato")
    }
    if (length(fill) < 2) {
        fill = c("tomato", "dodgerblue")
    }
    anygset = FALSE
    stopifnot(max.scale %in% c("abs", "rel") | (max.scale > 0 & 
        max.scale <= 1))
    if (!is.null(KS2)) {
        if (!inherits(KS2, "KS") & !(is.matrix(KS2))) {
            stop(print("KS2 is not of class KS or matrix"))
        }
            
		dm1 = KSdomain(KS)
		dm2 = KSdomain(KS2)
		dm = set_union(dm1,dm2)
		if (is.matrix(KS)) {
			KS = convert(convert(KS),domain = dm)
		}	
		if (is.matrix(KS2)) {
			KS2 = convert(convert(KS2),domain = dm)
		}
        if (is.matrix(KS)) {
            Freq <- rep(1, ncol(KS))
            KSdf = data.frame(t(KS), Freq)
			EE = matrix(0,ncol=length(dm),nrow=length(dm))  # Einheitsmatrix
			diag(EE) = 1
			KSdf = rbind(as.matrix(KSdf),cbind(EE,rep(0,length(dm))))# mit Freq 0 anhaengen
            KSdf = subKS(KSdf, 1:nrow(KS), keep.zero = TRUE)
			#KSdf = as.data.frame(xtabs(KSdf$Freq~., data = KSdf[,1:nrow(KS)]))
			
            ks.gset = 0
        }
        else {
            if (!is.null(attributes(KS)$memberships)) {
                anygset = TRUE
                KSdf = convert(KS, return.dataframe = T, domain = dm)
                KSdf = subKS(KSdf, 1:(ncol(KSdf) - 1), keep.zero = TRUE)
				#KSdf = as.data.frame(xtabs(KSdf$Freq~., data = KSdf[,1:(ncol(KSdf)-1)]))
				
                ks.gset = 1
            }
            else {
                Freq <- rep(1, length(KS))
                KSdf = data.frame(convert(KS, return.dataframe = TRUE, domain = dm), 
                  Freq)
				EE = matrix(0,ncol=length(dm),nrow=length(dm))  # Einheitsmatrix
				diag(EE) = 1
				KSdf = rbind(as.matrix(KSdf),cbind(EE,rep(0,length(dm))))# mit Freq 0 anhaengen
                KSdf = subKS(KSdf, 1:(ncol(KSdf) - 1), keep.zero = TRUE)
				#KSdf = as.data.frame(xtabs(KSdf$Freq~., data = KSdf[,1:(ncol(KSdf)-1)]))
                ks.gset = 0
            }
        }
		
        if (is.matrix(KS2)) {
			Freq <- rep(1, ncol(KS2))
            KS2df = data.frame(t(KS2), Freq)
			EE = matrix(0,ncol=length(dm),nrow=length(dm))  # Einheitsmatrix
			diag(EE) = 1
			KS2df = rbind(as.matrix(KS2df),cbind(EE,rep(0,length(dm))))# mit Freq 0 anhaengen
                
            KS2df = subKS(KS2df, 1:nrow(KS2), keep.zero = TRUE)
			#KS2df = as.data.frame(xtabs(KS2df$Freq~., data = KS2df[,1:nrow(KS2)]))
            ks2.gset = 0
        }
        else {
            if (!is.null(attributes(KS2)$memberships)) {
                anygset = TRUE
                KS2df = convert(KS2, return.dataframe = T)
                KS2df = subKS(KS2df, 1:(ncol(KS2df) - 1), keep.zero = TRUE)
				#KS2df = as.data.frame(xtabs(KS2df$Freq~., data = KS2df[,1:(ncol(KS2df)-1)]))
                ks2.gset = 1
            }
            else {
                Freq <- rep(1, length(KS2))
                KS2df = data.frame(convert(KS2, domain = dm, return.dataframe = T), 
                  Freq)
				EE = matrix(0,ncol=length(dm),nrow=length(dm))  # Einheitsmatrix
				diag(EE) = 1
				KS2df = rbind(as.matrix(KS2df),cbind(EE,rep(0,length(dm))))# mit Freq 0 anhaengen
                
                KS2df = subKS(KS2df, 1:(ncol(KS2df) - 1), keep.zero = TRUE)
				#KS2df = as.data.frame(xtabs(KS2df$Freq~., data = KS2df[,1:(ncol(KS2df)-1)]))
                ks2.gset = 0
            }
        }
        if (!anygset) {
            max.scale = "abs"
        }
        KSdf = data.frame(sapply(KSdf, as.integer))
        ind1 = which(KSdf$Freq > 0 & KS2df$Freq > 0)
        ind2 = which(KSdf$Freq > 0 & KS2df$Freq == 0)
        ind3 = which(KSdf$Freq == 0 & KS2df$Freq > 0)
        f1 = KSdf$Freq[c(ind1, ind2, ind3)] * ks.gset
        f2 = KS2df$Freq[c(ind1, ind2, ind3)] * ks2.gset
        if (max.scale != "abs") {
            f1 = f1/sum(KSdf$Freq)
            f2 = f2/sum(KS2df$Freq)
        }
        KS = t(KSdf[c(ind1, ind2, ind3), 1:(ncol(KSdf) - 1)]) - 
            1
        rownames(KS) = dm
        grp = c(rep(1, length(ind1)), rep(2, length(ind2)), rep(3, 
            length(ind3)))
        nmax = max(max.scale, f1, f2)
        if (max.scale == "rel") {
            nmax = max(c(f1, f2))
        }
        if (max.scale == "abs") {
            nmax = max(c(f1, f2, 1))
        }
        ratio1 = (f1)/nmax
        diff1 = (f1 - f2)
        abs1 = abs(diff1)/nmax
        sign1 = sign(diff1) >= 0
    }
    else {
        if (inherits(KS, "set")) {
            KS = convert(KS)
        }
        if (inherits(KS, "gset") & !(inherits(KS, "set"))) {
            dm = KSdomain(KS)
            anygset = TRUE
            KSdf = convert(KS, return.dataframe = T)
            KS = t(KSdf[, 1:(ncol(KSdf) - 1)])
            nmax = max(KSdf$Freq)
            ratio1 = (KSdf$Freq)/nmax
            abs1 = ratio1
            sign1 = rep(T, length(ratio1))
            rownames(KS) = dm
        }
        grp = rep(1, ncol(KS))
    }
    if (!is.null(ord)) {
        KS = KS[ord, ]
    }
    ord1 = do.call(order, c(data.frame(t(KS)), decreasing = T))
    KS = KS[, ord1]
    grp = grp[ord1]
    if (anygset) {
        abs1 = abs1[ord1]
        sign1 = sign1[ord1]
        ratio1 = ratio1[ord1]
    }
    ord2 = order(colSums(KS))
    KS = KS[, ord2]
    grp = grp[ord2]
    if (anygset) {
        abs1 = abs1[ord2]
        sign1 = sign1[ord2]
        ratio1 = ratio1[ord2]
    }
    y = colSums(KS)
    yn = max(y)
    tt = table(y)
    xn = max(tt)
    x = do.call(c, lapply(tt, function(x) {
        return(seq(from = -(x - 1), to = (x - 1), by = 2))
    }))
    n = ncol(KS)
    incmat1 = apply(KS, 2, function(x) {
        tmp = KS - x
        return(apply(tmp, 2, function(y) {
            return(min(y) >= 0)
        }))
    })
    diag(incmat1) = FALSE
    if (!is.null(KS2)) {
        incmat1[which(grp == 2), ] = F
        incmat1[, which(grp == 2)] = F
    }
    incmat2 = apply(incmat1, 2, function(x) {
        if (any(x)) {
            tmp = matrix(incmat1[, which(x)] * x, ncol = sum(x))
            return(rowSums(tmp) > 0)
        }
        else {
            return(rep(0, length(x)))
        }
    })
    incmat = t(incmat1 - incmat2)
    if (!is.null(KS2)) {
        incmat1 = apply(KS, 2, function(x) {
            tmp = KS - x
            return(apply(tmp, 2, function(y) {
                return(min(y) >= 0)
            }))
        })
        diag(incmat1) = FALSE
        incmat1[which(grp == 3), ] = F
        incmat1[, which(grp == 3)] = F
        incmat2 = apply(incmat1, 2, function(x) {
            if (any(x)) {
                tmp = matrix(incmat1[, which(x)] * x, ncol = sum(x))
                return(rowSums(tmp) > 0)
            }
            else {
                return(rep(0, length(x)))
            }
        })
        incmat = (incmat + t(incmat1 - incmat2)) > 0
    }
    x1 = matrix(rep(x, times = n), ncol = n)
    x2 = t(x1)
    y1 = matrix(rep(y, times = n), ncol = n)
    y2 = t(y1)
    dx = x2 - x1
    dy = y2 - y1
    coords = mapply(function(x, y) c(x, y), x = x, y = y, SIMPLIFY = F)
    dy[dy == 0] = 1
    dy[incmat == 0] = 1
    dz = dx/dy
    cmat = matrix(0, ncol = n, nrow = n)
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            if (dy[i, j] > 1) {
                tmp = mapply(function(x, y) c(x, y), x = 1:(dy[i, 
                  j] - 1) * dz[i, j] + x1[i, j], y = y1[i, j] + 
                  1:(dy[i, j] - 1), SIMPLIFY = F)
                if (any(tmp %in% coords)) {
                  cmat[i, j] = ifelse(sign(dz[i, j]) == 0, (-1)^(y1[i, 
                    j]%%2), sign(dz[i, j]))
                }
            }
        }
    }
    x = x + max(x)
    xm = max(x)
    if (xm == 0) {
        x = x + 0.5
    }
    else {
        x = x/xm
    }
    y = y - min(y)
    ym = max(1, y)
    yxmax = max(which(tt == max(tt)))/ym
    y = y/ym
    if (anygset) {
        W = 0.8/xn
        shape = "rect"
    }
    else {
        W = sapply(y, function(x) min((x + 1)/xn/(yxmax + 1) * 
            0.8, 1))
    }
    dev.new(width = 1024, height = 768)
    vp0 <- viewport(x = 0, y = 0, w = 1, h = 1, just = c("left", 
        "bottom"), name = "vp0")
    pushViewport(vp0)
    vp1 <- viewport(x = 0.5, y = 0.5, w = 0.98 - max(W)/2/(1 + 
        max(W)), h = 0.98 - 0.3/ym/2, just = "centre", name = "vp1")
    pushViewport(vp1)
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            if (incmat[i, j] == 1) {
                grid.curve(x1 = x[i], x2 = x[j], y1 = y[i], y2 = y[j], 
                  curvature = cmat[i, j]/(8^(1/(dy[i, j] + xn))), 
                  square = F, vp = vp1, gp = gpar(col = bg[max(grp[i], 
                    grp[j])], lwd = 3))
            }
        }
    }
    if (anygset) {
        bg = rep(bg[1], 3)
    }
    if (shape == "rect") {
        grid.rect(x = x, y = y, width = W, height = 0.3/ym, just = "centre", 
            default.units = "npc", name = NULL, gp = gpar(col = border, 
                fill = sapply(grp, function(x) bg[x]), alpha = 1), 
            draw = TRUE, vp = vp1)
    }
    if (shape == "roundrect") {
        mapply(function(a, b, d, e) {
            grid.roundrect(x = a, y = b, width = d, height = e, 
                r = unit(0.1, "snpc"), just = "centre", default.units = "npc", 
                name = NULL, gp = gpar(col = border, fill = sapply(grp, 
                  function(x) bg[x]), alpha = 1), vp = vp1)
        }, a = x, b = y, d = W, e = 1/ym/3)
    }
    if (shape == "ellipsis") {
        mapply(function(a, b, d, e) {
            grid.roundrect(x = a, y = b, width = d, height = e, 
                r = unit(0.5, "snpc"), just = "centre", default.units = "npc", 
                name = NULL, gp = gpar(col = border, fill = sapply(grp, 
                  function(x) bg[x]), alpha = 1), vp = vp1)
        }, a = x, b = y, d = W, e = 1/ym/3)
    }
    if (anygset) {
        w2 = abs1 * W
        x2 = x - W/2 + ratio1 * W
        y2 = y
        h2 = 0.3/ym
        grid.rect(x = x2[which(sign1)], y = y2[which(sign1)], 
            width = w2[which(sign1)], height = h2, just = "right", 
            default.units = "npc", name = NULL, gp = gpar(col = border, 
                fill = fill[2], alpha = 1), draw = TRUE, vp = vp1)
        if (any(!sign1)) {
            grid.rect(x = x2[which(!sign1)], y = y2[which(!sign1)], 
                width = w2[which(!sign1)], height = 0.3/ym, just = "left", 
                default.units = "npc", name = NULL, gp = gpar(col = border, 
                  fill = fill[1], alpha = 1), draw = TRUE, vp = vp1)
        }
        mapply(function(x, y) {
            grid.lines(x = c(x, x), y = c(y - h2/2, y + h2/2), 
                default.units = "npc", arrow = NULL, name = NULL, 
                gp = gpar(lwd = 3), draw = TRUE, vp = vp1)
        }, y = y2, x = x2)
    }
    grid.text(label = apply(KS, 2, function(x) {
        if (any(x == 1)) {
            return(paste(" ", do.call(paste, as.list(rownames(KS)[x == 
                1])), " "))
        }
        else {
            return("{ }")
        }
    }), y = y, x = x, just = "centre", rot = 0, gp = gpar(cex = lab.cex), 
        vp = vp1)
    if (anygset) {
        grid.text(label = paste("max.scale = ", round(nmax, digits = 3)), 
            x = min(x) + 0.3/xn, y = max(y) - 0.3/yn)
        popViewport(1)
    }
}



plothasse = function(KS, empirical =FALSE,v=1,concat=":",lab.cex=2){
	red = discred(as.KS(as.relation(as.KS(KS),empirical=empirical,v=v)),concat=concat)
	item = relation_incidence(as.relation(red))
	y = colSums(item)
	ord = order(y)
	item = item[ord,ord]
	y = rep(1,ncol(item))
	M = item
	diag(M) = 0	
	for( i in 1:ncol(item) ){
		y = apply(M*y,2,max)+1
	}
	tt = table(y)
	xn = max(tt)
	x=do.call(c,lapply(tt,function(x){
			return(seq(from=-(x-1),to=(x-1),by=2))
	}))
	
	n = ncol(item)
	
	incmat1 = apply(item,2,function(z){
				tmp =item-z
				return(apply(tmp,2,function(y){
							return(min(y)>=0)
				}))
	})

	diag(incmat1) = FALSE
	
	incmat2 = apply(incmat1,2,function(z){
			if(any(z)){
				tmp = matrix(incmat1[,which(z)]*z,ncol=sum(z))
				return(rowSums(tmp)>0)
			}else{
				return(rep(0,length(z)))	
			}
	})
	incmat = t(incmat1 - incmat2)
	
	x1 = matrix(rep(x,times = n),ncol = n)
	x2 = t(x1) #matrix(rep(x,times = n),ncol = n,byrow = T)
	y1 = matrix(rep(y,times = n),ncol = n)
	y2 = t(y1) #matrix(rep(y,times = n),ncol = n,byrow = T)

	dx = x2-x1
	dy = y2-y1
	
	coords = mapply(function(x,y) c(x,y),x=x,y=y,SIMPLIFY = F)
	
	dy[dy == 0] = 1
	dy[incmat == 0] = 1
	dz = dx/dy
	
	cmat = matrix(0,ncol=n,nrow=n)
	for( i in 1:(n-1) ){
		for(j in (i+1):n){
			if(dy[i,j] > 1){
				tmp = mapply(function(x,y) c(x,y), x = 1:(dy[i,j]-1)*dz[i,j] + x1[i,j], y = y1[i,j] + 1:(dy[i,j]-1),SIMPLIFY = F)

				if( any( tmp %in% coords ) ){
					cmat[i,j] = ifelse(sign(dz[i,j])==0, (-1)^(y1[i,j]%%2),sign(dz[i,j])) 	
				}

			}
		}
	}
	
	x = x+max(x)
	xm = max(x)
	if(xm == 0){
		x = x + 0.5
	}else{
		x = x/xm
	}
	y = y-min(y)
	ym = max(1,y)
	yxmax = max(which(tt==max(tt)))/ym
	y = y/ym
	
	W = sapply(y,function(x) min( (x+1)/xn/(yxmax+1)*0.8 ,1 ))
	
	dev.new(width = 1024, height = 768)
	vp0 <- viewport(x = 0, y = 0, w = 1, h = 1, just = c("left","bottom"), name = "vp0")
    pushViewport(vp0)
	vp1 <- viewport(x = 0.5, y = 0.5, w = 0.98-max(W)/2/(1+max(W)), h = 0.98-0.3/ym/2, just = "centre", name = "vp1")
    pushViewport(vp1)	
	
	for(i in 1:(n-1)){
		for(j in (i+1):n){
			if(incmat[i,j] == 1){
				grid.curve(x1=x[i],x2=x[j],y1=y[i],y2=y[j],curvature = cmat[i,j]/(8^(1/(dy[i,j]+xn))),
				 square = F,vp=vp1)
			}
		}		
	}
	
#		grid.rect(x = x, y = y, width = W, height = 0.3/ym, just = "centre", 
#        default.units = "npc", name = NULL, draw = TRUE, vp = vp1)

		grid.rect(x = x, y = y, width = 0.3/ym, height = 0.3/ym, just = "centre", 
        default.units = "npc", name = NULL, draw = TRUE, vp = vp1)


	grid.text(label = as.list(rownames(item))
				,y = y, x = x, just = "centre", rot = 0, gp = gpar(cex = lab.cex),vp=vp1)

#	if (anygset){
#	grid.rect(x = min(x) +0.05, y = max(y) -0.05, 
#            width = 0.06, height = 0.06, just = "centre", 
#            default.units = "npc", name = NULL, gp = gpar(col = border, 
#                fill = fill[1], alpha = 1), draw = TRUE)
#	grid.text(label = "x",x = min(x) +0.1,y=max(y) -0.05)
# 
#	grid.rect(x = min(x) +0.05, y = max(y) -0.12, 
#            width = 0.06, height = 0.06, just = "centre", 
#            default.units = "npc", name = NULL, gp = gpar(col = border, 
#                fill = fill[2], alpha = 1), draw = TRUE)
#	grid.text(label = "y",x = min(x) +0.1,y=max(y) -0.12)
#}



#	grid.text(label = apply(item,2,function(z){
#				if(any(z==1)){
#				 	return( paste( " ", do.call(paste,as.list(rownames(item)[z==1]))," ") )
#				}else{
#					return("{ }")
#				}
#				})
#				,y = y, x = x, just = "centre", rot = 0, gp = gpar(cex = max(1.3,3^(1/2/sqrt(xn)))),vp=vp1)
				
	popViewport(1)
	
}
