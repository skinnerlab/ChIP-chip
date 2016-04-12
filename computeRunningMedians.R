require("Ringo")
computeRM <- 
function (object, 
		  eSetA=NULL,
		  probeAnno, 
		  design = rep(1,time=ncol(object)),
		  groupings = 1:factor(sampleNames(object)), 
		  allChr = c(1:19, "X", "Y"), 
		  winHalfSize = 400, 
		  min.probes = 5, 
		  quant = 0.5, 
		  combineReplicates = FALSE, 
		  checkUnique = FALSE, 
		  uniqueCodes = c(1),
		  verbose = TRUE) 
{
	stopifnot(all(is.character(allChr)), 
			is.numeric(quant), (quant >= 0) & (quant <= 1), length(quant) == 1)
	 exprs <- t(apply(object$M,1,function(x) x*design))
	 rownames(exprs) <- object$genes$PROBE_ID
	 eSetA <- object$A
	 rownames(eSetA) <- object$genes$PROBE_ID
	 if (combineReplicates){ 
	        grouping <- as.factor(groupings)
	 } else { grouping <- factor(colnames(exprs)) }
     newExprs <- matrix(NA, nrow = nrow(exprs), ncol = nlevels(grouping))
	 colnames(newExprs) <- grouping
	 rownames(newExprs) <- object$genes$PROBE_ID

	 if (!is.null(eSetA)) newExprsA <- newExprs

	 counts <- numeric(length=nrow(exprs))
	 
	 for (chr in allChr) {
		# chr <- allChr[2]
        if (verbose) 
            cat("\nChromosome", chr, "...")
		idx <- which(probeAnno$CHROMOSOME == chr)
		chridx <- probeAnno[idx,"PROBE_ID"]
		chrsta <- probeAnno[idx,"POSITION"]
        chrend <- chrsta + probeAnno[idx,"LENGTH"]
        chrmid <- round((chrsta + chrend)/2)
        if (checkUnique) {
            chruni <- probeAnno[idx,"COUNT"]
            stopifnot(length(chruni) == length(chridx))
            chridx <- chridx[chruni %in% uniqueCodes]
            chrmid <- chrmid[chruni %in% uniqueCodes]
        }
		#ind <- which(chridx %in% intersect(chridx,rownames(newExprs)))
		ind <- match(chridx,rownames(newExprs))
		cat(length(ind)," probes\n")
		for (i in 1:nlevels(grouping)) {
			# i <- 1
            modSamples <- which(grouping == levels(grouping)[i])
            if (verbose) 
                cat(colnames(exprs)[modSamples], "... ")
			combined.pos <- rep(chrmid, each = length(modSamples))
			combined.dat <- as.vector(t(exprs[ind, modSamples, drop = FALSE]))
            slidingRes <- sliding.quantile(positions = combined.pos, 
                scores = combined.dat, half.width = winHalfSize, 
                prob = quant, return.counts = TRUE)
            slidingRes <- slidingRes[seq(1, nrow(slidingRes) + 
                1 - length(modSamples), by = length(modSamples)), 
                , drop = FALSE]
			chrrm <- slidingRes[, "quantile"]
			slidingRes[, "count"] <- slidingRes[, "count"]/length(modSamples)
			
			areBelow <- slidingRes[, "count"] < min.probes
			if (any(areBelow)){
				chrrm[areBelow] <- NA
			}
			stopifnot(length(chrrm) == length(chrmid))
			newExprs[ind, i] <- chrrm
			if (!is.null("eSetA")){ 
				combined.A <- as.vector(t(eSetA[ind, modSamples, drop = FALSE]))
				slidingA <- sliding.quantile(positions = combined.pos, 
						scores = combined.A, half.width = winHalfSize, 
						prob = quant, return.counts = TRUE)
				slidingA <- slidingA[seq(1, nrow(slidingA) + 
						1 - length(modSamples), by = length(modSamples)), 
						, drop = FALSE]
				chrA <- slidingA[,"quantile"]
				if (any(areBelow)){
					chrA[areBelow] <- NA
				}
				stopifnot(length(chrA) == length(chrmid))
				newExprsA[ind, i] <- chrA
				
			}
        }
		table(slidingRes[,"count"])
		counts[ind] <- slidingRes[, "count"]
    }
	cat("\n")
    sample.labels <- vector("character", nlevels(grouping))
    for (i in 1:nlevels(grouping)) sample.labels[i] <- as.character(levels(grouping)[i])
    object$M <- newExprs
	object$A <- newExprsA
	return(object)
}

