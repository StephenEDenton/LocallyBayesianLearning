# GlobalKalmanFilter.R
# Published in Kruschke & Denton (2010)
# Stephen Denton

rm(list=ls(all=TRUE))
graphics.off()

# Default values and toggles
nSubj = 40
fitParams = F
saving = FALSE
permuteTrials = TRUE
dynChange = FALSE
outputMu = FALSE
outputCov = FALSE
trialOutput = FALSE
outSubjRespProb = FALSE
nDigS = 4 # Significant digits for output
# Set seed
seed = 2009
set.seed(seed)

# Kalman Filter Initial Parameters
vOut = .1 # output noise; should be non-zero. Bigger => slower mean shift.
initVarOut = 100 # Variance for context cue output connections
initVarDirOut = 1/10000  # Small variance for direct output "specialist" connections 
outWt = 1 # Initial output weight values for set connections
uncertAdd = 0  # Additive uncertianty for dynamic change (U)
uncertMult = 1  # Multiplicative uncertainty for dynamic change (D)
chDec = 1 # Choice decisiveness parameter

params2fit = c("vOut")
params2fit = c("vOut", "initVarOut", "initVarDirOut")
params2fit = c("initVarOut", "initVarDirOut")

initParams = NULL
for (i in 1:length(params2fit) ) {
  initParams = c(initParams, eval(parse(text=params2fit[i])) )
}
# initParams = c(0.06305312)  # --> SSD = 470.9497 


# Notify user that program is invoked
cat('\n---- Global Kalman Filter Called ----\n')
# Create filename for saving
filename = paste0("GKFfit",nSubj,"Subj", format(Sys.time(),"%y-%m-%d-%H-%M-%S"),".Rdata")

## EXPERIMENT STRUCTURE --------------------------------------------------------
nCues = 11
nOut = 2
nExp = 3
nPh = 2
trnPats = vector("list", nExp)
nBlks4Ph = matrix(0,nExp,nPh)
nPats4Blk = nBlks4Ph
## Backward Blocking
expIdx = 1
nBlks4Ph[expIdx,] = c( 3, 2 ) # Number of blocks in phase 1 and 2
trnPatsB1 = matrix( c(
    rep( c( 1,1,0,0,0,1,0,0,1,0,0,  1 ) , 1 ),
    rep( c( 1,1,0,0,0,0,1,1,0,0,0,  2 ) , 1 ),
    rep( c( 0,0,1,1,0,0,1,1,0,0,0,  1 ) , 1 ),
    rep( c( 0,0,1,1,0,1,0,0,1,0,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,0,1,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,1,0,0,1,  2 ) , 1 ),
    rep( c( 1,1,0,0,0,1,0,0,0,0,1,  1 ) , 1 ),
    rep( c( 1,1,0,0,0,0,1,0,0,1,0,  2 ) , 1 ),
    rep( c( 0,0,1,1,0,0,0,1,0,0,1,  1 ) , 1 ),
    rep( c( 0,0,1,1,0,0,0,0,1,1,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,1,0,0,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,1,0,0,0,0,1,  2 ) , 1 )
    ), ncol=12, byrow=T )
trnPatsB2 = matrix( c(
    rep( c( 1,0,0,0,0,1,0,0,1,0,0,  1 ) , 1 ),
    rep( c( 1,0,0,0,0,0,1,1,0,0,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,0,1,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,1,0,0,1,  2 ) , 1 ),
    rep( c( 1,0,0,0,0,1,0,0,0,0,1,  1 ) , 1 ),
    rep( c( 1,0,0,0,0,0,1,0,0,1,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,1,0,0,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,1,0,0,0,0,1,  2 ) , 1 )
    ), ncol=12, byrow=T )
# Set up additional variable to be used later to permute trials
nPats4Blk[expIdx,] = c( dim(trnPatsB1)[1], dim(trnPatsB2)[1] )
# Combine all all blocks into one training structure.
trnPats[[expIdx]] = rbind(
  matrix(1,nBlks4Ph[expIdx,][1],1) %x% trnPatsB1 ,
  matrix(1,nBlks4Ph[expIdx,][2],1) %x% trnPatsB2 )
colnames( trnPats[[expIdx]] ) = c("A","B","C","D","E","S1f","S1j","S2f","S2j","S3f","S3j","R")

## Forward Blocking
expIdx = 2
nBlks4Ph[expIdx,] = c( 4, 2 ) # Number of blocks in phase 1 and 2
trnPatsB1 = matrix( c(
    rep( c( 1,0,0,0,0,1,0,0,1,0,0,  1 ) , 1 ),
    rep( c( 1,0,0,0,0,0,1,1,0,0,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,0,1,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,1,0,0,1,  2 ) , 1 ),
    rep( c( 1,0,0,0,0,1,0,0,0,0,1,  1 ) , 1 ),
    rep( c( 1,0,0,0,0,0,1,0,0,1,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,1,0,0,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,1,0,0,0,0,1,  2 ) , 1 )
    ), ncol=12, byrow=T )
trnPatsB2 = matrix( c(
    rep( c( 1,1,0,0,0,1,0,0,1,0,0,  1 ) , 1 ),
    rep( c( 1,1,0,0,0,0,1,1,0,0,0,  2 ) , 1 ),
    rep( c( 0,0,1,1,0,0,1,1,0,0,0,  1 ) , 1 ),
    rep( c( 0,0,1,1,0,1,0,0,1,0,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,0,1,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,1,0,0,1,  2 ) , 1 ),
    rep( c( 1,1,0,0,0,1,0,0,0,0,1,  1 ) , 1 ),
    rep( c( 1,1,0,0,0,0,1,0,0,1,0,  2 ) , 1 ),
    rep( c( 0,0,1,1,0,0,0,1,0,0,1,  1 ) , 1 ),
    rep( c( 0,0,1,1,0,0,0,0,1,1,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,1,0,0,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,1,0,0,0,0,1,  2 ) , 1 )
    ), ncol=12, byrow=T )
# Set up additional variable to be used later to permute trials
nPats4Blk[expIdx,] = c( dim(trnPatsB1)[1], dim(trnPatsB2)[1] )
# Combine all blocks into one training structure.
trnPats[[expIdx]] = rbind(
    matrix(1,nBlks4Ph[expIdx,][1],1) %x% trnPatsB1 ,
    matrix(1,nBlks4Ph[expIdx,][2],1) %x% trnPatsB2 )
colnames( trnPats[[expIdx]] ) = colnames( trnPats[[1]] )

## Backward Blocking
expIdx = 3
nBlks4Ph[expIdx,] = c( 3, 2 ) # Number of blocks in phase 1 and 2
trnPatsB1 = matrix( c(
    rep( c( 1,1,0,0,0,1,0,0,1,0,0,  1 ) , 1 ),
    rep( c( 1,1,0,0,0,0,1,1,0,0,0,  2 ) , 1 ),
    rep( c( 0,0,1,1,0,0,1,1,0,0,0,  1 ) , 1 ),
    rep( c( 0,0,1,1,0,1,0,0,1,0,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,0,1,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,1,0,0,1,  2 ) , 1 ),
    rep( c( 1,1,0,0,0,1,0,0,0,0,1,  1 ) , 1 ),
    rep( c( 1,1,0,0,0,0,1,0,0,1,0,  2 ) , 1 ),
    rep( c( 0,0,1,1,0,0,0,1,0,0,1,  1 ) , 1 ),
    rep( c( 0,0,1,1,0,0,0,0,1,1,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,1,0,0,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,1,0,0,0,0,1,  2 ) , 1 )
    ), ncol=12, byrow=T )
trnPatsB2 = matrix( c(
    rep( c( 1,0,0,0,0,1,0,0,1,0,0,  1 ) , 1 ),
    rep( c( 1,0,0,0,0,0,1,1,0,0,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,0,1,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,0,0,1,0,0,1,  2 ) , 1 ),
    rep( c( 1,0,0,0,0,1,0,0,0,0,1,  1 ) , 1 ),
    rep( c( 1,0,0,0,0,0,1,0,0,1,0,  2 ) , 1 ),
    rep( c( 0,0,0,0,1,0,1,0,0,1,0,  1 ) , 1 ),
    rep( c( 0,0,0,0,1,1,0,0,0,0,1,  2 ) , 1 )
    ), ncol=12, byrow=T )
# Set up additional variable to be used later to permute trials
nPats4Blk[expIdx,] = c( dim(trnPatsB1)[1], dim(trnPatsB2)[1] )
# Combine all all blocks into one training structure.
trnPats[[expIdx]] = rbind(
    matrix(1,nBlks4Ph[expIdx,][1],1) %x% trnPatsB1 ,
    matrix(1,nBlks4Ph[expIdx,][2],1) %x% trnPatsB2 )
colnames( trnPats[[expIdx]] ) = colnames( trnPats[[1]] )

# Set up special test items for specified experiment structure
nTest = 8
test = matrix( c(
   rep( c( 1,0,1,0,0,1,0,0,1,0,0 ) , 1 ),
   rep( c( 1,0,0,1,0,1,0,0,1,0,0 ) , 1 ),
   rep( c( 0,1,1,0,0,1,0,0,1,0,0 ) , 1 ),
   rep( c( 0,1,0,1,0,1,0,0,1,0,0 ) , 1 ),
   rep( c( 1,0,1,0,0,0,1,1,0,0,0 ) , 1 ),
   rep( c( 1,0,0,1,0,0,1,1,0,0,0 ) , 1 ),
   rep( c( 0,1,1,0,0,0,1,1,0,0,0 ) , 1 ),
   rep( c( 0,1,0,1,0,0,1,1,0,0,0 ) , 1 )
   ) , ncol=11, byrow=T )

# Set up test recording matrix and labels
out4Test = matrix(0, nrow=nTest, ncol=nOut)
cueLabels = colnames(trnPats[[1]])[-(nCues+1)]
outLabel = colnames(trnPats[[1]])[nCues+1]
testLabels = character(nTest)
respLabels = character(nOut)
for (tIdx in 1:nTest ) {
    testLabels[tIdx] = paste(cueLabels[which(test[tIdx,]==1)],collapse=".")
}
for (rIdx in 1:nOut ) { respLabels[rIdx] = paste(outLabel, rIdx, sep="") }
dimnames(out4Test) = list(testLabels, respLabels)

# Set up recording matrices
resp4Subj = array(0,c(nTest,nOut,nSubj,nExp))
rownames(resp4Subj) = dimnames(out4Test)[[1]]
colnames(resp4Subj) = dimnames(out4Test)[[2]]
respPer = array(0, c(2,nOut,nExp))
rownames(respPer) = c("A.C   ", "B.D   ")
colnames(respPer) = c("A/B", "C/D")

## DATA TO FIT -----------------------------------------------------------------

humanRespPer = respPer
humanRespPer[,,1] = matrix( c(
    c( 54.5, 45.5 ),
    c( 42.1, 57.9 )
    ) , ncol=2, byrow=T )
humanRespPer[,,2] = matrix( c(
    c( 59.8, 40.2 ),
    c( 44.5, 55.5 )
    ) , ncol=2, byrow=T )
humanRespPer[,,3] = matrix( c(
    c( 52.3, 47.8 ),
    c( 45.3, 54.7 )
    ) , ncol=2, byrow=T )

## FUNCTIONS -------------------------------------------------------------------

updateMu = function( covMat, inVec, outVar, teach, muOld ) {
    # Make sure that inVec is a column vector:
    dim( inVec ) = c( length( inVec ), 1 )
    denom = as.numeric( outVar + t(inVec) %*% covMat %*% inVec )
    muNext = muOld + ( covMat %*% inVec %*% ( teach - t(inVec) %*% muOld ) / denom )
    return( muNext )
}

updateCov = function( covMat, inVec, outVar ) {
    # Make sure that inVec is a column vector:
    dim( inVec ) = c( length( inVec ) , 1 )
    CNextNumer = covMat %*% inVec %*% t(inVec) %*% covMat
    CNextDenom = as.numeric( outVar + t(inVec) %*% covMat %*% inVec )
    CNext = covMat - CNextNumer / CNextDenom
    return( CNext )
}

# Experiment function
experiment = function( params ) {
  
  # Reset seed each time experiment function is called
  set.seed(seed)
  # If there are negative params set return a NaN as the fit
  if ( any(params<0) ) {
    return( NaN )
  }
  cat("\n###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat("\"", filename, "\"\n", sep='')
  cat("Parameters used:\n  <", paste(params2fit, collapse = ", "), ">\n")
  cat(params)
  
  for ( i in 1:length(params2fit) ) {
    assign(params2fit[i], params[i] )
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

      # Set up KF matrices
      COut = array(0,c(nCues,nCues,nOut))

      muOut = matrix( 0, nrow=nOut, ncol=nCues )
      for ( outIdx in 1:nOut ) {
          COut[,,outIdx] = diag( initVarOut, nCues )
          diag(COut[,,outIdx])[6:11] = initVarDirOut
      }
      
      # Set direct connections between 'specialists' and outcomes
      muOut[1,seq(6,11,2)] = outWt
      muOut[2,seq(7,11,2)] = outWt

      #-------------------------------------------------------------------------
      # Do the training/updating.

      # Trial loop, goes to nTrials +1 bcz initial state counts as cycle thru the loop.
      for ( trialIdx in 1:(nPats+1) ) {

          # Print some output
          if (outputMu) {
              cat("\nmuOut (rows are seperate KFs) \n"); print(round(muOut,5))
          }
          if (outputCov) {
              cat("COut\n"); print(COut[,,1])
          }

          # Exit loop before update after the last trial
          if (trialIdx == (nPats+1)) { break }

          # Apply dynamic change
          if ( dynChange ) {
            for ( outIdx in 1:nOut ) {
              muOut[outIdx,] = D %*% muOut[outIdx,]
              COut[,,outIdx] = D %*% COut[,,outIdx] %*% t(D) + U
            }
          }

          # Get next trial datum from trnPats array
          avec = trnPats[[expIdx]][permuteIdx[trialIdx],-(nCues+1)]
          teachVal = trnPats[[expIdx]][permuteIdx[trialIdx],nCues+1]
          tvals = rep(0,nOut); tvals[teachVal] = 1

          # Compute Output Activations
          outAct = muOut %*% avec
          outInpt = avec

          # Update output weights
          for ( outIdx in 1:nOut) {
            muOut[outIdx,] = updateMu( COut[,,outIdx], outInpt, vOut, tvals[outIdx], muOut[outIdx,] )
            COut[,,outIdx] = updateCov( COut[,,outIdx], outInpt, vOut )
          }

          # Show trial by trial output
          if ( trialOutput ) {
            cat('\nTrial', trialIdx, '========================================\nCues:'
                , avec, ', Teacher:', tvals, '\n')
            cat('Output Act:', round(outAct, nDigS), '\n')
          }

      } # end training loop

      # Test the trained network
      for (tIdx in 1:nTest) {
          out4Test[tIdx,] = muOut %*% test[tIdx,]
      }
      
      # Convert to choice probabilities (using exponential softmax)
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

    # Collapse choice probabilites to match table in paper (hard coded for ease)
    respPer[1,,expIdx] = colMeans(rbind(respProb[1:2,],respProb[5:6,c(2,1)])) * 100
    respPer[2,,expIdx] = colMeans(rbind(respProb[3:4,],respProb[7:8,c(2,1)])) * 100
  }

  cat("\nCollapsed Resp. Percentages (1->BB, 2->FB, 3->FB) \n")
  print(respPer)

  SSD = sum((respPer-humanRespPer)^2)
  cat("Sum squared deviation:", SSD, "\n")
  return(SSD)
}

if ( fitParams ) {
  # bestFit = optimize(experiment, c(0,1))
  bestFit = optim(initParams, experiment, control=list(trace=1,reltol=1e-16,maxit=10000))
  print(bestFit)
} else {
  experiment( initParams )
}

print(humanRespPer)

if (saving) { save.image(file = filename) }

# End Program