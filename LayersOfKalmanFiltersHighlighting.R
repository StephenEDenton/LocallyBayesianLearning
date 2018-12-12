# LayersOfKalmanFiltersHighlighting.R

# Stephen Denton

rm(list=ls(all=TRUE))
graphics.off()
library(rootSolve)

# Default values and toggles
nSubj = 1
fitParams = F
plotOutput = F
permuteTrials = FALSE
dynChange = FALSE
outputMu = T
outputCov = F
trialOutput = T
outSubjRespProb = T
nDigS = 4 # Significant digits for output
# Set seed
seed = 2009
set.seed(seed)

# Kalman Filter Initial Parameters
vOut = 0.5 # output noise; should be non-zero. Bigger => slower mean shift.
vHid = 10 # hidden noise
vAll = 1
initVarOut = 100 # Variance for output connections
initVarHid = 100 # Variance for hidden lateral connections
initVarAll = 100
initVarS = 1/100 # Small variance for direct connections
attWt = 1 # Initial attention weight for direct connections
outWt = 1 # Initial output weight values for set connections
uncertAdd = 0  # Additive uncertianty for dynamic change (U)
uncertMult = 1  # Multiplicative uncertainty for dynamic change (D)
chDec = 1 # Choice decisiveness parameter

params2fit = c("vAll")

initParams = NULL
for (i in 1:length(params2fit) ) {
    initParams = c(initParams, eval(parse(text=params2fit[i])) )
}

# Notify user that program is invoked
cat('\n---- Local Bayes (layers of KFs) Called ----\n')

## EXPERIMENT STRUCTURE --------------------------------------------------------
nCues = 3
nOut = 2
nExp = 2
nPh = 1
trnPats = vector("list", nExp)
nBlks4Ph = matrix(0, nExp, nPh)
nPats4Blk = nBlks4Ph

# Highlighting
expIdx = 1
trnPats[[expIdx]] = matrix( c(
   rep( c( 1,1,0,  1 ) , 10 ),
   rep( c( 0,1,1,  2 ) , 10 )
   ) , ncol=4, byrow=T )
colnames( trnPats[[expIdx]] ) = c("PE","I","PL","E/L")

# Ambiguous Cue
expIdx = 2
trnPats[[expIdx]] = matrix( c(
   rep( c( 1,1,0,  1 ) , 1 ),
   rep( c( 0,1,1,  2 ) , 1 ),
   rep( c( 1,1,0,  1 ) , 1 ),
   rep( c( 0,1,1,  2 ) , 1 ),
   rep( c( 1,1,0,  1 ) , 1 ),
   rep( c( 0,1,1,  2 ) , 1 ),
   rep( c( 1,1,0,  1 ) , 1 ),
   rep( c( 0,1,1,  2 ) , 1 ),
   rep( c( 1,1,0,  1 ) , 1 ),
   rep( c( 0,1,1,  2 ) , 1 ),
   rep( c( 1,1,0,  1 ) , 1 ),
   rep( c( 0,1,1,  2 ) , 1 ),
   rep( c( 1,1,0,  1 ) , 1 ),
   rep( c( 0,1,1,  2 ) , 1 ),
   rep( c( 1,1,0,  1 ) , 1 ),
   rep( c( 0,1,1,  2 ) , 1 ),
   rep( c( 1,1,0,  1 ) , 1 ),
   rep( c( 0,1,1,  2 ) , 1 ),
   rep( c( 1,1,0,  1 ) , 1 ),
   rep( c( 0,1,1,  2 ) , 1 )
   ) , ncol=4, byrow=T )
colnames( trnPats[[expIdx]] ) = c("PE","I","PL","E/L")

# Set up special test items for specified experiment structure
nTest = 6
test = matrix( c(
   rep( c( 1,1,0 ) , 1 ),
   rep( c( 0,1,1 ) , 1 ),
   rep( c( 0,1,0 ) , 1 ),
   rep( c( 1,0,1 ) , 1 ),
   rep( c( 1,0,0 ) , 1 ),
   rep( c( 0,0,1 ) , 1 )
   ) , ncol=3, byrow=T )

# Set up test recording matrix and labels
out4Test = matrix(0, nrow=nTest, ncol=nOut)
cueLabels = colnames(trnPats[[1]])[-(nCues+1)]
outLabel = colnames(trnPats[[1]])[nCues+1]
testLabels = character(nTest)
for (tIdx in 1:nTest ) {
    testLabels[tIdx] = paste(cueLabels[which(test[tIdx,]==1)],collapse=".")
}
#respLabels = character(nOut)
#for (rIdx in 1:nOut ) { respLabels[rIdx] = paste(outLabel, rIdx, sep="") }
respLabels = c("E", "L")
dimnames(out4Test) = list(testLabels, respLabels)

# Set up KF matrices
CHid = array(0,c(nCues,nCues,nCues))
COut = array(0,c(nCues,nCues,nOut))

# Set up recording matrices
resp4Subj = array(0,c(nTest,nOut,nSubj,nExp))
rownames(resp4Subj) = dimnames(out4Test)[[1]]
colnames(resp4Subj) = dimnames(out4Test)[[2]]
respPer = array(0, c(nTest,nOut,nExp))
rownames(respPer) = dimnames(out4Test)[[1]]
colnames(respPer) = dimnames(out4Test)[[2]]

## FUNCTIONS -------------------------------------------------------------------

updateCov = function( covMat , inVec , outVar ) {
    # Make sure that inVec is a column vector:
    dim( inVec ) = c( length( inVec ) , 1 )
    CNextNumer = covMat %*% inVec %*% t(inVec) %*% covMat
    CNextDenom = as.numeric( outVar + t(inVec) %*% covMat %*% inVec )
    CNext = covMat - CNextNumer / CNextDenom
    return( CNext )
}

updateMu = function( covMat , inVec , outVar , teach , muOld ) {
    # Make sure that inVec is a column vector:
    dim( inVec ) = c( length( inVec ) , 1 )
    denom = as.numeric( outVar + t(inVec) %*% covMat %*% inVec )
    muNext = muOld + ( covMat %*% inVec %*% ( teach - t(inVec) %*% muOld ) / denom )
    return( muNext )
}

# Function to find the hidden target that maximizes the p(out=teach)
findHid = function( mu, covMats, teachers, outVar, initPars ) {
    initPar = initPars[initPars!=0]
    parLoc = which(initPars!=0)
    hidAtMaxP = multiroot(zfun2, initPar, mu=mu, covMats=covMats, teachers=teachers,
        outVar=outVar, parLoc=parLoc) 
    hidAtMaxPars = initPars
    hidAtMaxPars[parLoc] = hidAtMaxP$root
    return( hidAtMaxPars )
}

# Function that will be minimized by to hidden teacher
zfun = function( x, mu, covMats, teachers, outVar, parLoc ) {
    fit = 1
    for ( outIdx in 1:length(teachers) ) {
        fit = fit * dnorm(teachers[outIdx], mean=(sum(mu[outIdx,parLoc] * x)),
            sd=sqrt( outVar + x %*% covMats[parLoc,parLoc,outIdx] %*% x ) )
    }
    return( fit )
}

# Derivative of above function which can be set to zero to find roots
zfun2 = function( x, mu, covMats, teachers, outVar, parLoc ) {
    fun = rep(0,length(x))
    for ( outIdx in 1:length(teachers) ) {
        variance = c(outVar + x %*% covMats[parLoc,parLoc,outIdx] %*% x)
        predDiff = c(teachers[outIdx] - mu[outIdx,parLoc]%*%x)
        Cx = covMats[parLoc,parLoc,outIdx] %*% x
        fun = fun + variance*predDiff*mu[outIdx,parLoc] + Cx*predDiff^2 - Cx*variance
    }
    return(fun)
}

plotKF = function( mu, covMat ) {
    wtPlotLim = 2
    wtPlotLen = 10*4+1
    
    # 1D marginals
    maxDens = 0
    for ( cueIdx in 1:nCues ) {
        cueMean = mu[cueIdx]
        cueSD = sqrt( covMat[cueIdx,cueIdx] )
        maxDens = max( maxDens , dnorm( cueMean , mean=cueMean , sd=cueSD ) )
    }
    for ( cueIdx in 1:nCues ) {
        wValComb = seq(-wtPlotLim,wtPlotLim,length=wtPlotLen)
        cueMean = mu[cueIdx]
        cueSD = sqrt( covMat[cueIdx,cueIdx] )
        pWcue = dnorm( wValComb , mean=cueMean , sd=cueSD )
        par(pty="s")
        par(mar=c(2,2,0.5,0.3))
        plot( wValComb , pWcue ,  type="l" , ylim=c(0,1), ann=FALSE, cex.axis=.7 )
    }
    # 2D marginals
    colComb = seq(-wtPlotLim,wtPlotLim,length=wtPlotLen)
    rowComb = seq(-wtPlotLim,wtPlotLim,length=wtPlotLen)
    z = matrix( 0 , nrow=length(rowComb) , ncol=length(colComb) )
    for ( colIdx in 1:(nCues-1) ) {
        for ( rowIdx in (colIdx+1):nCues ) {
            for ( rIdx in 1:dim(z)[1] ) {
                z[rIdx,] = dmvnorm(
                    cbind(rep(rowComb[rIdx],dim(z)[2]),colComb) ,
                    mean = c( mu[rowIdx] , mu[colIdx] ) ,
                    sigma = covMat[ (c(rowIdx,colIdx)) , (c(rowIdx,colIdx)) ] )
            }
            par(pty="s")
            par(mar=c(2,2,0.5,0.3))
            contour( colComb , rowComb , t(z) , drawlabels=F, cex.axis=.7 )
        }
    }
}

## Experiment function
experiment = function( params ) {

  # Reset seed each time experiment function is called
  set.seed(seed)
  # If there are negative params set return a NaN as the fit
  if ( any(params<0) ) {
      return( NaN )
  }
  cat("\n###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat("Parameters used:\n  <", paste(params2fit, collapse = ", "), ">\n")
  cat(params)
  
  for ( i in 1:length(params2fit) ) {
      assign(params2fit[i], params[i] )
  }

  if ("vAll" %in% params2fit) {
      vOut = vAll
      vHid = vAll
  }
  if ("initVarAll" %in% params2fit) {
      initVarHid = initVarAll
      initVarOut = initVarAll
  }
  if ( uncertMult > 1.05  ) {
      return( NaN )
  }
  
  # Set up KF matrices
  U = diag( uncertAdd, nCues ) # Uncertainty added every trial; could be zero.
  D = diag( uncertMult, nCues ) # Dynamic shift; will assume to be identity matrix.

  #-----------------------------------------------------------------------------
  # Initialize experiment
  
  for ( expIdx in 1:nExp ) {
    
    cat('\nExperiment', expIdx, '========================================\n')

    nPats = dim(trnPats[[expIdx]])[1]

    #---------------------------------------------------------------------------
    # Initialize each simulated subject
    for ( subjIdx in 1:nSubj ) {
      cat('S', subjIdx, ' ', sep='')
      if ( subjIdx %% 10 == 0 ) { cat('\n') }

      # Set up permutation index
      permuteIdx = NULL
      cumTrls = 0
      if ( permuteTrials ) {
        for ( phIdx in 1:dim(nBlks4Ph)[2] ) {
          for ( blkIdx in 1:nBlks4Ph[expIdx,phIdx] ) {
            permuteIdx = c(permuteIdx, sample((cumTrls+1):(cumTrls+nPats4Blk[expIdx,phIdx])) )
            cumTrls = cumTrls + nPats4Blk[expIdx,phIdx]
          }
        }
      } else {
        permuteIdx = 1:nPats
      }

      muHid = diag(nCues) * attWt
      for ( hidIdx in 1:nCues ) {
          CHid[,,hidIdx] = diag( initVarHid, nCues )
          CHid[hidIdx,hidIdx,hidIdx] = initVarS
      }
      muOut = matrix( 0, nrow=nOut, ncol=nCues )
      for ( outIdx in 1:nOut ) {
          COut[,,outIdx] = diag( initVarOut, nCues )
      }
      

      #-------------------------------------------------------------------------
      # Do the training/updating.

      # Trial loop, goes to nTrials +1 bcz initial state counts as cycle thru the loop.
      for ( trialIdx in 1:(nPats+1) ) {
      
          if ( plotOutput ) {

              # Plot current beliefs for the Output Kalman Filter
              windows(height=4, width=4*nOut, xpos=200, ypos=0)
              layoutm = NULL
              plotNum = 1
              for (outIdx in 1:nOut ) {
                  ilayoutm =  diag( plotNum:(plotNum+nCues-1) )
                  plotNum = plotNum + nCues - 1
                  ilayoutm[ lower.tri( ilayoutm ) ] = plotNum + 1:( nCues * (nCues-1) / 2 )
                  plotNum = plotNum + ( nCues * (nCues-1) / 2 ) + 1
                  layoutm = cbind(layoutm, ilayoutm)
                  # layoutm[1,2] = max(layoutm)+1
              }
              layout( layoutm )
              for (outIdx in 1:nOut ) {
                  plotKF( muOut[outIdx,], COut[,,outIdx] )
              }

              # Plot current beliefs for the Hidden Kalman Filters
              windows(height=4, width=4*nCues, xpos=40, ypos=600)
              layoutm = NULL
              plotNum = 1
              for ( hidIdx in 1:nCues ) {
                  ilayoutm = diag( plotNum:(plotNum+nCues-1) )
                  plotNum = plotNum + nCues - 1
                  ilayoutm[ lower.tri( ilayoutm ) ] = plotNum + 1:( nCues * (nCues-1) / 2 )
                  plotNum = plotNum + ( nCues * (nCues-1) / 2 ) + 1
                  layoutm = cbind(layoutm, ilayoutm)
              }
              layout( layoutm )
              for (hidIdx in 1:nCues) {
                  plotKF( muHid[hidIdx,], CHid[,,hidIdx] )
              }
          }

          # Print some output
          if (outputMu) {
              cat("\nmuOut (rows are seperate KFs) \n"); print(round(muOut,nDigS))
              cat("muHid (rows are seperate KFs) \n"); print(round(muHid,nDigS))
          }
          if (outputCov) {
              cat("COut\n"); print(COut)
              cat("CHid\n"); print(CHid)
          }

          # Exit loop before update after the last trial
          if (trialIdx == (nPats+1)) { break }

          # Apply dynamic change
          if ( dynChange ) {
            for ( hidIdx in 1:nCues ) {
              muHid[hidIdx,] = D %*% muHid[hidIdx,]
              CHid[,,hidIdx] = D %*% CHid[,,hidIdx] %*% t(D) + U
            }
            for ( outIdx in 1:nOut ) {
              muOut[outIdx,] = D %*% muOut[outIdx,]
              COut[,,outIdx] = D %*% COut[,,outIdx] %*% t(D) + U
            }
          }

          # Get next trial datum from trnPats array
          avec = trnPats[[expIdx]][permuteIdx[trialIdx],-(nCues+1)]
          teachVal = trnPats[[expIdx]][permuteIdx[trialIdx],nCues+1]
          tvals = rep(0,nOut); tvals[teachVal] = 1

          # Compute Attention Activations
          hidAtt = muHid %*% avec
          # Hidden activation is the product of inputs times att wts just computed
          hidAct = hidAtt * avec
          outAct = muOut %*% hidAct
          outInpt = hidAct
          
          # Update output weights
          for ( outIdx in 1:nOut ) {
            muOut[outIdx,] = updateMu( COut[,,outIdx], outInpt, vOut, tvals[outIdx], muOut[outIdx,] )
            COut[,,outIdx] = updateCov( COut[,,outIdx], outInpt, vOut )
          }
          # Find the best hidden teacher
          hidTeach = findHid( muOut, COut, tvals, vOut, hidAct )

          # Update hidden weights
          for ( hidIdx in 1:nCues ) {
            muHid[hidIdx,] = updateMu( CHid[,,hidIdx], avec, vHid, hidTeach[hidIdx], muHid[hidIdx,] )
            CHid[,,hidIdx] = updateCov( CHid[,,hidIdx], avec , vHid )
          }

          # Show trial by trial output
          if ( trialOutput ) {
            cat('\nTrial', trialIdx, '========================================\nCues:'
                , avec, ', Teacher:', tvals, '\n')
            cat('Output Act:', round(outAct,nDigS), '\n')
            cat('Hidden Att:', round(hidAtt,nDigS), '\n')
            cat('Hidden Act:', round(hidAct,nDigS), '\n')
            cat('Hidden Act for Max P:', round(hidTeach,nDigS), '\n')
          }
      } # end training loop

      # Test the trained network
      for (tIdx in 1:nTest) {
          # hidAtt = muHid %*% test[tIdx,]
          # out4Test[tIdx,] = muOut %*% (hidAtt * test[tIdx,])
          out4Test[tIdx,] = muOut %*% ((muHid %*% test[tIdx,]) * test[tIdx,])
      }

      # Convert to choice probabilities (using exponential rule
      choiceProb = out4Test
      choiceProb = exp(chDec*choiceProb) / rowSums(exp(chDec*choiceProb))

      if ( outSubjRespProb ) {
        cat('\nExpected Output for Cue Combinations:\n')
        print(out4Test)
        cat("\nRespone prob. using exp. choice rule, chDec =", chDec, "\n")
        print(choiceProb)
        cat('\n')
      }
      # Record subject response percentages
      resp4Subj[,,subjIdx,expIdx] = choiceProb
    }

    # Collapse response probabilities across subjects
    if ( nSubj > 1) {
      respProb = rowSums(resp4Subj[,,,expIdx], dims=2)/nSubj
    } else {
      respProb = resp4Subj[,,,expIdx]
    }
    
    respPer[,,expIdx] = respProb
  }

  return(respPer)
}

respPer = experiment( initParams )

print(respPer)

# End Program