#pragma rtGlobals=1		// Use modern global access method.

#include ":ProcessRawBeam"

Constant m = 3.817543e-26		//Na
//Constant m = 			//K
//Constant m = 			//Rb

Function ProcessDiffraction(fileName)
	String fileName
	
	Variable runTimer = 1
	If(runTimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Make/O/WAVE/N=6 WaveRefs
	// WaveRefs[0] = counts
	// WaveRefs[1] = position
	// WaveRefs[2] = counts error
	// WaveRefs[3] = counts resampled
	// WaveRefs[4] = position resampled
	// WaveRefs[5] = counts error resampled
	
	String path = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:091021:diffscans:"
	
	LoadAndDisplayScan2(path, fileName)
	Wave countsWave = WaveRefs[0]; Wave positionWave = WaveRefs[1]
	
	Variable maxPosition = 1.4e3		//microns
	RezeroScan2(positionWave, countsWave)
	
	CropScan(positionWave, countsWave, maxPosition)
	
	// Computer the error wave
	Wave errWave = GenErrWave(countsWave)
	WaveRefs[2] = errWave
	
	if(1)
		// Resample the position and counts waves and then recompute the error wave
		Variable resampleRate = .15
		resampleWaves(positionWave, countsWave, errWave, resampleRate)
		Wave resampledCounts = WaveRefs[3]; Wave resampledPosition = WaveRefs[4]
		Wave resampledErr = GenErrWave(resampledCounts)
		WaveRefs[5] = resampledErr
	
		DisplayScan(ResampledCounts, resampledPosition, resampledErr)
		GoFitDiffraction(ResampledCounts, resampledPosition, resampledErr)
		CleanUpGraph(resampledCounts, resampledErr)
	else
		// process entire data set
		DisplayScan(countsWave, positionWave, ErrWave)
		GoFitDiffraction(countsWave, positionWave, ErrWave)
		CleanUpGraph(countsWave, ErrWave)
	endif
	
	If(runTimer)
		Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
	endIf
End


//Constant nSigmas = 4.5
//Constant velPoints = 100

Function GenVelDist2(nSigmas,velPoints)
	Variable nSigmas,velPoints
	velPoints=5001
	
	Wave InitialGuesses = InitialGuesses1
	
	Variable FlowV = InitialGuesses[0]
	Variable SigV = InitialGuesses[1]
	
	Variable/G minVel = FlowV*(1-nSigmas/15); minVel = 0
	Variable/G maxVel = FlowV*(1+nSigmas/15); maxVel = 5000
	Variable/G velStep = (maxVel-minVel)/velPoints
	
	Make/D/O/N=(velPoints) velwave = minVel+p*velStep
	Duplicate/O velwave velprobs2
	SetScale/I x minVel,maxVel,"m/s", velprobs2
	Velprobs2 = exp(-(x-FlowV)^2/(2*SigV^2))*x^3
	Variable velsum = sum(velprobs2,-inf,inf)
	Velprobs2/=velsum
	
	if(0)
	Wavestats/Q velprobs
	print "Most probable velocity: ", V_maxloc
	print "velprobs[0]: ", velprobs2[0]
	print "velprobs[max]: ", velprobs2[(maxVel-minVel)/velStep]
	endif
End


Function DisplayVelocityDist()
	Wave velprobs, velwave

	Display/k=1 velprobs vs velwave
	Label left "Probability"
	Label bottom "Velocity (m/s)"
End




Function GoFitDiffraction(ywave, xwave, errwave)
	Wave ywave, xwave, errwave
	
	GenVelDist2(8,1000)
	
	Wave InitialGuesses = InitialGuesses1
	Wave ConstraintsConv
	
	Variable DiffordersToFit = 6
	Variable OtherParams = 10
	Duplicate/O InitialGuesses w_coef
	//DeletePoints (OtherParams+DiffordersToFit), (numpnts(w_coef)-DiffordersToFit), w_coef
	
	String hold = ProcessHoldAndConstraintsWave2()

//ywave[429]=nan
//ywave[412]=nan
	FuncFit/TBOX=768/L=400/H=hold FitDiffractionScan W_coef  ywave /X=xwave /R /D /W=errwave /I=1 /C=ConstraintsConv		
	//Add /N to turn off updates
	//Add /O to only graph initial guesses
	
	Variable reducedChiSqrd = (V_chisq/(V_npnts)
	printf "Reduced Chi-Squared: %g\r", reducedChiSqrd
	
	Variable vratio = w_coef[0]/w_coef[1]
	printf "vratio: %2.2f\r", vratio
End



Function TestFitDiffraction()
	Variable testpnts = 401
	Variable halfScanWidth = 2000
	Variable resolution = 2*halfScanWidth/testpnts
	Make/D/O/N=(testpnts) testYdif, testXdif
	SetScale/I x -halfScanWidth, halfScanWidth, "um", testXdif, testYdif
	testXdif = x
	
	GenVelDist2(5,100)
	
	Wave InitialGuesses1
	Duplicate/O InitialGuesses1 w_coef

	FitDiffractionScan(w_coef, testYdif, testXdif) 
End



Function FitDiffractionScan(pw, yw, xw) : FitFunc
	Wave pw, yw, xw
	
	// initialize fittable variables
	Variable v0	= pw[0]			// flow velocity
	Variable sv = pw[1]				// velocity distribution
	Variable cmol = pw[2]			// molecule fraction
	Variable amp = pw[3]			// normalization amplitude
	Variable x0 = pw[4]			// x-offset
	Variable background = pw[5]	// background count rate
	Variable s1 = pw[6]			// 1s width
	Variable s2 = pw[7]			// 2s width
	Variable GaussAmp = pw[8]		// raw beam wing amplitude
	Variable GaussWidth = pw[9]	// raw beam wing width
									// order weights are not explicitly assigned
							
	Variable OtherParams = 10	// used in assigning order weights
	
	// regenerate normalized velocity distribution for new flow velocity and velocity distribution
	// use this to find normalization factor
	Wave velprobs2
	Velprobs2 = exp(-(x-v0)^2/(2*sv^2))*x^3
	Variable velsum = area(velprobs2,0,inf)
	Velprobs2/=velsum
	
	
	
	// set geometrical parameters
	Variable L1d = 2.3974				// dist to detector in m
	Variable Lgd = 96.5*25.4/1000; L1d=Lgd	// override 1g-det dist with magic grating distance
	
	//Variable m = 22.9898*1.66054e-27	// Na mass	// now set at top of procedure
	Variable a = 100e-9						// grating period in microns
	Variable h = 6.62607E-34				// h in kg*m^2/s^2
	Variable xPreFactor = L1d*h/m/a		// = 0.41613 m^2/s for sodium
//print xPreFactor



	// generate diffraction pattern for an infinitely thin beam and detector
	
	Variable dX = 1 // spatial resolution of model in microns		

	Variable diffScanHalfwidth = 2000	//microns		// must be large enough that all orders fit entirely inside -halfwidth to +halfwidth!
	Variable diffPnts = diffScanHalfwidth*2/dX +1
	Make/D/O/N=(diffPnts) FofX=0, SingleOrder = 0, SingleOrderNeg = 0, SingleOrderMol = 0, SingleOrderMolNeg = 0
	SetScale/P x -diffScanHalfwidth, dX, "um", FofX, SingleOrder, SingleOrderMol, SingleOrderNeg, SingleOrderMolNeg
	//SetScale/P x diffScanHalfwidth, -dX, "um", SingleOrderNeg, SingleOrderMolNeg
	
	//Here's the most important lines
	// First we ask: What velocity would be needed to put a diffracted atom at a position x? the answer is v = ( L*h / m*a ) * n/x = xPreFactor*n/x
	// Then we lookup what the probability for an atom to have that velocity, given our v0 and sig, P(v(x)) 
	// This probability is just the diffraction amplitude FofX at that location x.   
	// Note that x is in microns but VelProbs scaling is in m/s, so x must be converted to meters first.
		
	Variable totalOrders = numpnts(pw)-OtherParams
	variable n=1
	Variable cn = pw[n+OtherParams-1] 
	Variable DiffOrderNorm
	Variable DiffOrderNormMol
	do
		// Depending on the scaling of FofX and VelProbs the diffraction orders will not all have integrated unit probability automatically
		// so we compute each diffraction order separately to find the proper scaling.
//Variable timerRefNum = Startmstimer
		cn = pw[n+OtherParams-1] 
		
		if(0)
			SingleOrder = velprobs2(xPreFactor*n/x*1e6)		//this operation is one of the slower ones so we only do it once and use the Reverse command
			DiffOrderNorm = sum(SingleOrder)
			if(DiffOrderNorm < .01)	// if the model window isn't large enough to include most of the highest order then the normalization will be very wrong
				DiffOrderNorm=1		// this prevents the model from blowing up, although it's not very accurate
				//print diffordernorm
			endif
			SingleOrder/=DiffOrderNorm
			Reverse/P singleorder /d=singleorderneg			
		
			SingleOrderMol = velprobs2(xPreFactor*n/x/2*1e6)
			DiffOrderNormMol = sum(SingleOrderMol)
			SingleOrderMol/=DiffOrderNormMol
			Reverse/P singleordermol /d=singleordermolneg
			
			//FofX += cn* (   1/DiffOrderNorm*(velprobs(xPreFactor/(x*1e-6/n))+velprobs(-xPreFactor/(x*1e-6/n))) + cmol/DiffOrderNormMol* (velprobs(xPreFactor/(x*1e-6/n*2))+velprobs(-xPreFactor/(x*1e-6/n*2)))  )
			//FofX += cn* (  (SingleOrder+SingleOrderNeg) + cmol* (SingleOrderMol+SingleOrderMolNeg)  )
			MatrixOP/O/S FofX = FofX+ cn* (  (SingleOrder+SingleOrderNeg) + cmol* (SingleOrderMol+SingleOrderMolNeg)  )
		endif
		
		SingleOrder = 1/velsum*(xPrefactor)^4*n^4/(x)^5*exp(-(xPrefactor*n/x*1e6-v0)^2/(2*sv^2))*dX*1e24
		Reverse/P SingleOrder /D=SingleOrderNeg	
		//FofX +=  cn*1/velsum*(xPrefactor)^4*n^4/(x)^5*exp(-(xPrefactor*n/x*1e6-v0)^2/(2*sv^2))*dX*1e24 
		//FofX +=  cn*1/velsum*(xPrefactor)^4*n^4/(-x)^5*exp(-(-xPrefactor*n/x*1e6-v0)^2/(2*sv^2))*dX*1e24 
		
		SingleOrderMol = 1/velsum*(xPrefactor/2)^4*n^4/(x)^5*exp(-(xPrefactor/2*n/x*1e6-v0)^2/(2*sv^2))*dX*1e24
		Reverse/P SingleOrderMol /D=SingleOrderMolNeg
		//FofX +=  cn*cmol*1/velsum*(xPrefactor/2)^4*n^4/(x)^5*exp(-(xPrefactor/2*n/x*1e6-v0)^2/(2*sv^2))*dX*1e24 
		//FofX +=  cn*cmol*1/velsum*(xPrefactor/2)^4*n^4/(-x)^5*exp(-(-xPrefactor/2*n/x*1e6-v0)^2/(2*sv^2))*dX*1e24 
		
		MatrixOP/O/S FofX = FofX+ cn* (  (SingleOrder+SingleOrderNeg) + cmol* (SingleOrderMol+SingleOrderMolNeg)  )

//Print "time: ", StopMSTimer(timerRefNum)*10^-6
		n+=1
	while(n<=totalOrders || n<=1)
	
	//Add the 0th order atoms and molecules and scale the whole pattern
	FofX[x2pnt(FofX,0)] = 1+cmol

	
	//Generate the raw beam waves
	Variable rawScanHalfWidth = 800
	Variable rawBeamPnts = rawScanHalfWidth*2/dX +1
	Make/D/O/N=(rawBeamPnts) rawY, rawX
	SetScale/P x -rawScanHalfWidth, dX, "um", rawY, rawX
	rawX = x
	
	// Choose beam shape
	if(0)	// pure convolution of 1s, 2s, and det
		Make/D/O/N=1 rawbeamparams={1, s1,s2,0,0}
		FitRawBeam(rawbeamparams, rawY, rawX)
	elseif(1) // convolution plus gaussian
		Variable GaussRatio = 1/40
		Make/D/O/N=1 rawbeamparams={1, s1, s2, 0, 0, GaussAmp, GaussWidth}
		FitRawBeamGauss(rawbeamparams, rawY, rawX)
	elseif(0) // gaussian only
		rawY = GaussAmp*Gauss(x, 0, GaussWidth)
	endif
	
	Variable rawYmax = wavemax(rawY)
	Variable rawYarea = area(rawY)
	rawY /= rawYarea
	
	FofX*=amp//*rawYarea
	
	//convolve the infinitely thin beam/detector diffraction pattern with the beam shape
	Duplicate/O FofX FofXConvolved	
	Convolve/A rawY FofXConvolved
	
	//assign output
	yw = FofXConvolved(xw[p]-x0)+background	
End







Function/S ProcessHoldAndConstraintsWave2()
	Variable n
	Wave HoldWave = HoldConv
	Wave/T Constraints = ConstraintsConv, Constraints1= ConstraintsConv1, Constraints2 = ConstraintsConv2
	Duplicate/O/T Constraints1 Constraints1Dup
	Duplicate/O/T Constraints2 Constraints2Dup
	
	String hold = ""
	For(n=0;n<numpnts(HoldWave);n+=1)
		If(HoldWave[n]==1)
			hold += "1"
		Else
			hold += "0"
		EndIf
	EndFor
	
	Variable/G numHolds=0
	For(n=0;n<numpnts(HoldWave);n+=1)
		If(HoldWave[n]==1)
			Constraints1Dup[n]=""; Constraints2Dup[n]=""
			numHolds += 1
		EndIf
	EndFor
	
	Concatenate/O/NP/T {Constraints1Dup, Constraints2Dup}, $nameofwave(Constraints)
	Variable conLength = numpnts(Constraints)
	For(n=0;n<conLength;n+=1)
		If (stringmatch(Constraints[n],"")==1)
			DeletePoints n,1, Constraints
			n -= 1
		Endif
	EndFor
	
	Return hold
end










Function VelocityComb()
	Wave velwave, velprobs
	
	Variable totalOrders = 3
	
	Duplicate/O velwave xPos1 
	
	Variable L1d = 2.3974				// dist to detector in m
	Variable m = 22.9898*1.66054e-27	// Na mass
	Variable a = 100e-9						// grating period in microns
	Variable h = 6.62607E-34				// h in kg*m^2/s^2
	
	xPos1 = (L1d*h/m/a)/velwave*1e6 // in microns
	
	Variable VelPnts = numpnts(velwave)
	Make/O/N=(velPnts*totalOrders) xPosTot, FofX
	variable n=1
	do
		xPosTot[(n-1)*VelPnts,n*VelPnts] = xPos1
		n+=1
	while(n<=totalOrders)
	
	FofX = velprobs
End



Function VelocityComb2()
	Wave velwave, velprobs
	
	Wave cn
	Variable totalOrders = numpnts(cn)-1
	
	Variable L1d = 2.3974				// dist to detector in m
	Variable m = 22.9898*1.66054e-27	// Na mass
	Variable a = 100e-9						// grating period in microns
	Variable h = 6.62607E-34				// h in kg*m^2/s^2
	Variable xPreFactor = L1d*h/m/a; print xprefactor	//same as NaFactor from old fitting routine

	Variable xPnts = 3000
	Variable ScanHalfwidth = 1500	//microns
	Variable ScanResolution = ScanHalfwidth*2/xPnts; //print scanresolution
	Make/O/N=(xPnts) FofX=0
	SetScale/P x -ScanHalfwidth, ScanResolution, "um", FofX
	
	variable n=1
	do
		FofX += cn[n]*(velprobs(xPreFactor/(x*1e-6/n))+velprobs(-xPreFactor/(x*1e-6/n)))
		n+=1
	while(n<=totalOrders)
	
	Variable cmol= 0.05	//molecule fraction
	n=1
	do
		FofX += cmol*cn[n]*(velprobs(xPreFactor/(x*1e-6/n*2))+velprobs(-xPreFactor/(x*1e-6/n*2)))
		n+=1
	while(n<=totalOrders)
	
	FofX[x2pnt(FofX,0)] = cn[0]+cmol
End



Function QuickTest()
	Wave FofX, testY
	
	Duplicate/O testY testYnorm
	Variable testYsum = sum(testYnorm)
	testYnorm /= testYsum
	
	Duplicate/O FofX FofXConvolved
	
	Convolve/A testYnorm FofXConvolved
End





Window DiffParamsTable() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(141,44,1099,477) KWave1,ParamList1,InitialGuesses1,w_coef,w_sigma,HoldConv,ConstraintsConv1
	AppendToTable ConstraintsConv2,ConstraintsConv
	ModifyTable format(Point)=1,width(KWave1)=54,width(ParamList1)=104,width(InitialGuesses1)=84
	ModifyTable width(ConstraintsConv1)=90,width(ConstraintsConv2)=98,width(ConstraintsConv)=92
EndMacro





Function GenPofX()
	Make/D/O/N=700 PofX
	SetScale/P x 1, 1, "um", PofX
	
	Variable L1d = 2.3974				// dist to detector in m
	Variable Lgd = 96.5*25.4/1000; L1d=Lgd	// override 1g-det dist with magic grating distance
	
	//Variable m = 22.9898*1.66054e-27	// Na mass	// now set at top of procedure
	Variable a = 100e-9						// grating period in microns
	Variable h = 6.62607E-34				// h in kg*m^2/s^2
	Variable xPreFactor = L1d*h/m/a		// = 0.41613 m^2/s for sodium
	
	Variable v0=2000
	Variable sig=100
	
	Make/D/O/N=5000 PofV
	SetScale/P x, 0,1, "m/s", PofV
	PofV = x^3*exp(-(x-v0)^2/(2*sig^2))
	Variable PofVsum = sum(pofV)
	print PofVsum
	PofV/=PofVsum
	
	Variable n = 2
	PofX = 1/PofVsum*(xPrefactor)^3*(n/x*1e6)^3*exp(-(xPrefactor*n/x*1e6-v0)^2/(2*sig^2))*xPrefactor*n/x^2*1e6
	
	Variable PofXsum = sum(PofX,0,inf)
	print PofXsum
	PofX/=PofXsum
End