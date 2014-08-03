#pragma rtGlobals=1		// Use modern global access method.

#include ":ProcessRawBeam"
#include ":Constants"


///////////////////////////////////////////////////////////////////////////////
//
//	Functions for processing diffraction data to determine beam velocity.
//	The code is poorly commented, but I'm not going to improve it since we're not using it actively.
//	It was written in 2009 by Will Holmgren. This intro was written in 2013, also by Will.
//
//	The very first thing you should do is run:
//		DiffParamsTable()
//	This acts as a control panel of sorts.
//	After setting your constraints, you need to run ProcessHoldAndConstraintsWave2()
//
//	ProcessDiffraction(fileName) is the master function. 
//		Takes a single data file name as its argument.
//		Before running this function (or any that call it), make a global string named 'atom'
//		and set it to 'Na', 'K', etc. This selects the appropriate mass.
//		Loads the data, rezeros the scan (this can sometimes go wrong), and crops it to a max +/- position.
//		Next, optionally resample the data, display the scan, fit it, and then clean up the graph.
//		
//	FitMultipleScans(series, start, stop) will process a diffraction data series.
//	BatchProcess() was used at one point to help control the initial guesses for different files (I think)
//
//	GoFitDiffraction(ywave, xwave, errwave) is the prep function for fitting to the diffraction model.
//	FitDiffractionScanIsotopes(pw, yw, xw) is the primary diffraction model.
//
//	TestFitDiffractionIsotopes() is useful for testing your initial guesses.
//
//	VelSigMolGraph(series) and VelSigMolChiGraph(series) graph the results of your series fits.
//
//	KillFileResults(series, fileNumber) eliminates the offending point from the graphs.
//
//	DisplayVelocityDist() shows the velocity distribution determined by your fit.
//
//	MinusScansStats(wavein, waveout) and PlusScansStats(wavein, waveout) can help when
//	trying to see systematic difference between scan directions.
//
///////////////////////////////////////////////////////////////////////////////



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
	
	// Make a global string variable called dateFolder with the name of the data folder
	// or just set the data folder in the first half of the if-else statement
	String basePath = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:" 
	String path
	String dateFolderLoc = StrVarOrDefault("dateFolder", "")
	if(StringMatch(dateFolderLoc,""))
		dateFolderLoc = "090720"
		path = basePath + dateFolderLoc + ":"
	else
		path = basePath + dateFolderLoc + ":"
	endif
	
	Variable/G m
	String atomLoc = StrVarOrDefault("atom", "Na")
	if(StringMatch(atomLoc, "Na"))
		m = 3.817543e-26
		Variable/G numIsotopes = 1
		Make/D/O massWave = {Na23Mass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {1}
	elseif(StringMatch(atomLoc, "K"))
		m = 6.49243e-26
		Variable/G numIsotopes = 1
		Make/D/O massWave = {KAvgMass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {1}
	elseif(StringMatch(atomLoc, "KIso"))
		Variable/G numIsotopes = 3
		Make/D/O massWave = {K39mass, K40mass, K41mass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {K39frac, K40frac, K41frac}
	elseif(StringMatch(atomLoc, "Rb"))
		m = 1.41923e-25
		Variable/G numIsotopes = 1
		Make/D/O massWave = {RbAvgMass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {1}
	elseif(StringMatch(atomLoc, "RbIso"))
		Variable/G numIsotopes = 2
		Make/D/O massWave = {Rb85mass, Rb87mass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {Rb85frac, Rb87frac}
	elseif(StringMatch(atomLoc, "Cs"))
		Variable/G numIsotopes = 1
		Make/D/O massWave = {CsAvgMass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {1}
	else
		DoAlert 0,"Uncertain atom type"
		Abort
	endif
	
	LoadAndDisplayScan2(path, fileName, 0)
	Wave countsWave = WaveRefs[0]; Wave positionWave = WaveRefs[1]
	
	Variable maxPosition = 1.2e3		//microns
	RezeroScan2(positionWave, countsWave)
	
	CropScan(positionWave, countsWave, maxPosition)
	
	// Computer the error wave
	Wave errWave = GenErrWave(countsWave)
	WaveRefs[2] = errWave
	
	if(1)
		// Resample the position and counts waves and then recompute the error wave
		Variable resampleRate = .1
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
	
	//Variable DiffordersToFit = 6
	Variable OtherParams = 11
	Duplicate/O InitialGuesses w_coef, eps
	eps = 1e-6; eps[8]=10; eps[1]=1
	
	String hold = ProcessHoldAndConstraintsWave2()
	//Variable V_fitOptions=2
	Variable V_FitTol = .00001

	//FuncFit/M=2/L=400/H=hold/ODR=0 FitDiffractionScan W_coef  ywave /X=xwave /R /D /I=1 /C=ConstraintsConv /W=errwave 	/E=eps
	FuncFit/M=2/L=400/H=hold/ODR=0 FitDiffractionScanIsotopes W_coef  ywave /X=xwave /R /D /I=1 /C=ConstraintsConv /W=errwave 	/E=eps
	//Add /N to turn off updates
	//Add /O to only graph initial guesses
	
	Wave M_Covar
	Duplicate/O M_Covar, CorMat	 
	CorMat = M_Covar[p][q]/sqrt(M_Covar[p][p]*M_Covar[q][q])

	AppendResults(ParamList1)

	NVAR numHolds
	Variable/G reducedChiSqrd = V_chisq/(V_npnts-(numpnts(W_coef)-numHolds))
	printf "Reduced Chi-Squared: %g\r", reducedChiSqrd
	TextBox/C/N=text1/Z=1/X=85.00/Y=30.00 "\\F'Symbol'c\\F'Geneva'\\S2\\M/dof = "+num2str(reducedChiSqrd)
	
	Variable vratio = w_coef[0]/w_coef[1]
	printf "vratio: %2.2f\r", vratio
	
	String w_coefName = nameofwave(ywave) + "_w_coef"
	String w_sigmaName = nameofwave(ywave) + "_w_sigma"
	
	Duplicate/O w_coef $w_coefName; Wave w_coefSaved = $w_coefName
	Duplicate/O w_sigma $w_sigmaName; Wave w_sigmaSaved = $w_sigmaName
End



Function TestFitDiffraction()
	Variable testpnts = 2001
	Variable halfScanWidth = 2000
	Variable resolution = 2*halfScanWidth/testpnts
	Make/D/O/N=(testpnts) testYdif, testXdif
	SetScale/I x -halfScanWidth, halfScanWidth, "um", testXdif, testYdif
	testXdif = x
	
	GenVelDist2(5,100)
	
	Wave InitialGuesses1, w_coef
	Duplicate/O InitialGuesses1 w_coef
	
	Variable/G m
	String atomLoc = StrVarOrDefault("atom", "Na")
	if(StringMatch(atomLoc, "Na"))
		m = 3.817543e-26
	elseif(StringMatch(atomLoc, "He"))
		m = 6.64648e-27
	elseif(StringMatch(atomLoc, "K"))
		m = 6.49243e-26
	elseif(StringMatch(atomLoc, "Rb"))
		m = 1.41923e-25
	elseif(StringMatch(atomLoc, "Rb8587"))
		Variable/G numIsotopes = 2
		Make/D/O massWave = {Rb85mass, Rb87mass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {Rb85frac, Rb87frac}
	else
		DoAlert 0,"Uncertain atom type"
		Abort
	endif

	FitDiffractionScan(w_coef, testYdif, testXdif) 
End



Function TestFitDiffractionIsotopes()
	Variable testpnts = 2001
	Variable halfScanWidth = 2000
	Variable resolution = 2*halfScanWidth/testpnts
	Make/D/O/N=(testpnts) testYdif, testXdif
	SetScale/I x -halfScanWidth, halfScanWidth, "um", testXdif, testYdif
	testXdif = x
	
	GenVelDist2(5,100)
	
	Wave InitialGuesses1, w_coef
	Duplicate/O InitialGuesses1 w_coef
	
	Variable/G m
	String atomLoc = StrVarOrDefault("atom", "Na")
	if(StringMatch(atomLoc, "Na"))
		m = 3.817543e-26
	elseif(StringMatch(atomLoc, "K"))
		m = 6.49243e-26
	elseif(StringMatch(atomLoc, "Rb"))
		m = 1.41923e-25
	elseif(StringMatch(atomLoc, "Rb8587"))
		Variable/G numIsotopes = 2
		Make/D/O massWave = {Rb85mass, Rb87mass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {Rb85frac, Rb87frac}
	elseif(StringMatch(atomLoc, "Cs"))
		Variable/G numIsotopes = 1
		Make/D/O massWave = {CsAvgMass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {1}
	else
		DoAlert 0,"Uncertain atom type"
		Abort
	endif

	FitDiffractionScanIsotopes(w_coef, testYdif, testXdif) 
End





Function FitDiffractionScan(pw, yw, xw) : FitFunc
	Wave pw, yw, xw
	
	// initialize fittable variables
	Variable v0	= pw[0]			// flow velocity
	Variable sigv = pw[1]				// velocity distribution
	Variable cmol = pw[2]/(100-pw[2])			// molecule ratio
	Variable amp = pw[3]			// normalization amplitude
	Variable x0 = pw[4]			// x-offset
	Variable background = pw[5]	// background count rate
	Variable s1 = pw[6]			// 1s width
	Variable s2 = pw[7]			// 2s width
	Variable s3 = pw[8]
	Variable GaussAmpRatio = pw[9]		// raw beam wing amplitude
	Variable GaussWidth = pw[10]	// raw beam wing width
									// order ratios are not explicitly assigned
							
	Variable OtherParams = 11	// # of parameters other than order weights. used when looping over order weights
	
	// regenerate normalized velocity distribution for new flow velocity and velocity distribution
	// use this to find normalization factor 
	Wave velprobs2
	Variable sv = 2*sigv^2
	MultiThread Velprobs2 = exp(-(x-v0)^2/sv)*x^3
	Variable VelNorm// = area(velprobs2,0,inf); print VelNorm
	VelNorm = sqrt(2*pi)*sigv*v0*(v0^2+3*sigv^2)//; print VelNorm
	Velprobs2/=VelNorm
	
	




	// generate diffraction pattern for an infinitely thin beam and detector
	
	Variable dX = 1 // spatial resolution of model in microns		

	Variable diffScanHalfwidth = 2000	//microns		// must be large enough that all orders fit entirely inside -halfwidth to +halfwidth!
	Variable diffPnts = diffScanHalfwidth*2/dX +1
	Make/D/O/N=(diffPnts) FofX=0, SingleOrder = 0, SingleOrderNeg = 0, SingleOrderMol = 0, SingleOrderMolNeg = 0
	SetScale/P x -diffScanHalfwidth, dX, "um", FofX, SingleOrder, SingleOrderMol, SingleOrderNeg, SingleOrderMolNeg
	//SetScale/P x diffScanHalfwidth, -dX, "um", SingleOrderNeg, SingleOrderMolNeg
	
	//Here's the most important lines
	// First we ask: What velocity would be needed to put a diffracted atom at a position x? the answer is v = ( L*h / m*a ) * n/x = xPreFactor*n/x
	// Then we calculate what the probability for an atom to have that velocity is, given our v0 and sig, P(v(x)) 
	// The diffraction amplitude FofX at that location x is given by: .   
	// Note that x is in microns but VelProbs scaling is in m/s, so x must be converted to meters first.
		
	Variable totalOrders = numpnts(pw)-OtherParams
	variable n=1
	Variable cn = pw[n+OtherParams-1] 
	//Variable DiffOrderNorm
	//Variable DiffOrderNormMol
	
		// set geometrical parameters
	Variable L1d =	2.3724		// 2.3974				// dist from 1G to detector in m
	//Variable L1g2g = .93927; L1d=L1g2g
	//Variable Lgd = 96.5*25.4/1000; L1d=Lgd	// override 1g-det dist with magic grating distance
	
	Variable a = 100e-9						// grating period in microns
	Variable h = 6.62607E-34				// h in kg*m^2/s^2
	NVAR m
	Variable xPreFactor = L1d*h/(m*a) * 1e6		// = 416114 um*m/s for sodium
	//print xPreFactor
	
	Variable NormConst = xPrefactor^4/VelNorm*dX		// fewer unnecessary computations

	do
		// Depending on the scaling of FofX and VelProbs the diffraction orders will not all have integrated unit probability automatically
		// so we compute each diffraction order separately to find the proper scaling.
//Variable timerRefNum = Startmstimer
		cn = pw[n+OtherParams-1]// ; print n, cn
		
		Multithread		SingleOrder = n^4/(x)^5*exp(-(xPrefactor*n/x-v0)^2/sv)
		//SingleOrder = 1/VelNorm*(xPrefactor)^4*n^4/(x)^5*exp(-(xPrefactor*n/(x*1e-6)-v0)^2/(2*sigv^2))*dX*1e24
		//SingleOrder = abs(1/VelNorm*xPrefactor*n/(x)^2*exp(-(xPrefactor*n/(x*1e-6)-v0)^2/(2*sigv^2))*dX*1e6)
		Reverse/P SingleOrder /D=SingleOrderNeg	
		//FofX +=  cn*1/VelNorm*(xPrefactor)^4*n^4/(x)^5*exp(-(xPrefactor*n/x*1e6-v0)^2/(2*sigv^2))*dX*1e24 
		//FofX +=  cn*1/VelNorm*(xPrefactor)^4*n^4/(-x)^5*exp(-(-xPrefactor*n/x*1e6-v0)^2/(2*sigv^2))*dX*1e24 
		
		Multithread		SingleOrderMol = n^4/(x)^5*exp(-((xPrefactor/2)*n/x-v0)^2/sv)/16
		//SingleOrderMol = 1/VelNorm*(xPrefactor/2)^4*n^4/(x)^5*exp(-((xPrefactor/2)*n/x*1e6-v0)^2/(2*sigv^2))*dX*1e24
		//SingleOrderMol = abs(1/VelNorm*xPrefactor/2*n/(x)^2*exp(-((xPrefactor/2)*n/(x*1e-6)-v0)^2/(2*sigv^2))*dX*1e6)
		Reverse/P SingleOrderMol /D=SingleOrderMolNeg
		//FofX +=  cn*cmol*1/VelNorm*(xPrefactor/2)^4*n^4/(x)^5*exp(-(xPrefactor/2*n/x*1e6-v0)^2/(2*sigv^2))*dX*1e24 
		//FofX +=  cn*cmol*1/VelNorm*(xPrefactor/2)^4*n^4/(-x)^5*exp(-(-xPrefactor/2*n/x*1e6-v0)^2/(2*sigv^2))*dX*1e24 
		
		MatrixOP/O/S FofX = FofX + cn* (  (SingleOrder+SingleOrderNeg) + cmol* (SingleOrderMol+SingleOrderMolNeg)  )
		
//Print "time: ", StopMSTimer(timerRefNum)*10^-6
		n+=1
	while(n<=totalOrders)
	
	MatrixOP/O/S FofX=FofX*NormConst
	
	//Add the 0th order atoms and molecules and scale the whole pattern
	FofX[x2pnt(FofX,0)] = 1+cmol

	
	//Generate the raw beam waves
	Variable rawScanHalfWidth = 1600
	Variable rawBeamPnts = rawScanHalfWidth*2/dX +1
	Make/D/O/N=(rawBeamPnts) rawY, rawX
	SetScale/P x -rawScanHalfWidth, dX, "um", rawY, rawX
	rawX = x
	
	// Choose beam shape here by turning on (1) one if/elseif and turning off (0) the rest
	if(0)	// pure convolution of 1s, 2s, and det
		Make/D/O/N=1 rawbeamparams={1, s1,s2,0,0}
		FitRawBeam(rawbeamparams, rawY, rawX)
	elseif(1) // convolution plus gaussian
		Make/D/O/N=1 rawbeamparams={1, s1, s2, 0, 0, GaussAmpRatio, GaussWidth}
		FitRawBeamGauss(rawbeamparams, rawY, rawX)
	elseif(0) // convolution plus gaussian
		Make/D/O/N=1 rawbeamparams={1, s1, s2, 0, 0, GaussAmpRatio, GaussWidth, s3}
		FitRawBeam3sGauss(rawbeamparams, rawY, rawX)
	elseif(0) // two gaussians
		rawY = sqrt(2*pi)*GaussWidth*(Gauss(x, 0, GaussWidth)+GaussAmpRatio*Gauss(x,0,GaussWidth*50))
	endif
	
	Variable rawYmax = wavemax(rawY)
	Variable rawYarea = area(rawY)
	rawY /= rawYarea
	
	MultiThread FofX*=amp//*rawYarea
	
	//convolve the infinitely thin beam/detector diffraction pattern with the beam shape
	Duplicate/O FofX FofXConvolved	
	Convolve/A rawY FofXConvolved
	
	//assign output
	yw = FofXConvolved(xw[p]-x0)+background	
End




Function FitDiffractionScanIsotopes(pw, yw, xw) : FitFunc
	Wave pw, yw, xw
	
	// initialize fittable variables
	Variable v0	= pw[0]			// flow velocity
	Variable sigv = pw[1]				// velocity distribution
	Variable cmol = pw[2]/(100-pw[2])			// molecule ratio
	Variable amp = pw[3]			// normalization amplitude
	Variable x0 = pw[4]			// x-offset
	Variable background = pw[5]	// background count rate
	Variable s1 = pw[6]			// 1s width
	Variable s2 = pw[7]			// 2s width
	Variable s3 = pw[8]
	Variable GaussAmpRatio = pw[9]		// raw beam wing amplitude
	Variable GaussWidth = pw[10]	// raw beam wing width
									// order ratios are not explicitly assigned
							
	Variable OtherParams = 11	// # of parameters other than order weights. used when looping over order weights
	
	// regenerate normalized velocity distribution for new flow velocity and velocity distribution
	// use this to find normalization factor 
	Wave velprobs2
	Variable sv = 2*sigv^2
	MultiThread Velprobs2 = exp(-(x-v0)^2/sv)*x^3
	Variable VelNorm// = area(velprobs2,0,inf); print VelNorm
	VelNorm = sqrt(2*pi)*sigv*v0*(v0^2+3*sigv^2)//; print VelNorm
	Velprobs2/=VelNorm


	// generate diffraction pattern for an infinitely thin beam and detector
	
	Variable dX = 1 // spatial resolution of model in microns		

	Variable diffScanHalfwidth = 2000	//microns		// must be large enough that all orders fit entirely inside -halfwidth to +halfwidth!
	Variable diffPnts = diffScanHalfwidth*2/dX +1
	Make/D/O/N=(diffPnts) FofX=0, SingleOrder = 0, SingleOrderNeg = 0, SingleOrderMol = 0, SingleOrderMolNeg = 0, FofXsum=0
	SetScale/P x -diffScanHalfwidth, dX, "um", FofX, SingleOrder, SingleOrderMol, SingleOrderNeg, SingleOrderMolNeg
	//SetScale/P x diffScanHalfwidth, -dX, "um", SingleOrderNeg, SingleOrderMolNeg
	
	//Here's the most important lines
	// First we ask: What velocity would be needed to put a diffracted atom at a position x? the answer is v = ( L*h / m*a ) * n/x = xPreFactor*n/x
	// Then we calculate what the probability for an atom to have that velocity is, given our v0 and sig, P(v(x)) 
	// The diffraction amplitude FofX at that location x is given by: .   
	// Note that x is in microns but VelProbs scaling is in m/s, so x must be converted to meters first.
		
	Variable totalOrders = numpnts(pw)-OtherParams
	variable n=1
	Variable cn = pw[n+OtherParams-1] 
	//Variable DiffOrderNorm
	//Variable DiffOrderNormMol
	
		// set geometrical parameters
	Variable L1d =	2.3724		// 2.3974				// dist from 1G to detector in m
	//Variable L1g2g = .93927; L1d=L1g2g
	//Variable Lgd = 96.5*25.4/1000; L1d=Lgd	// override 1g-det dist with magic grating distance
	
	Variable a = 100e-9						// grating period in microns
	Variable h = 6.62607E-34				// h in kg*m^2/s^2
	
	Variable, xPreFactor, NormConst, m, isotopeFraction
	Wave massWave, isotopeFracs
	String FofXname

	NVAR numIsotopes
	Make/O/WAVE/N=(numIsotopes) FofXRefs
	
	Variable i = 0
	for(i=0; i<numIsotopes; i+=1)
	
		m = massWave[i]
		isotopeFraction = IsotopeFracs[i]
	
		xPreFactor = L1d*h/(m*a) * 1e6		// = 416114 um*m/s for sodium
		//print xPreFactor
	
		NormConst = xPrefactor^4/VelNorm*dX		// fewer unnecessary computations 
	
		n=1
		FofX=0
	
		do
			// Depending on the scaling of FofX and VelProbs the diffraction orders will not all have integrated unit probability automatically
			// so we compute each diffraction order separately to find the proper scaling.
			//Variable timerRefNum = Startmstimer
			cn = pw[n+OtherParams-1]// ; print n, cn
		
			Multithread		SingleOrder = n^4/(x)^5*exp(-(xPrefactor*n/x-v0)^2/sv)
			//SingleOrder = 1/VelNorm*(xPrefactor)^4*n^4/(x)^5*exp(-(xPrefactor*n/(x*1e-6)-v0)^2/(2*sigv^2))*dX*1e24
			//SingleOrder = abs(1/VelNorm*xPrefactor*n/(x)^2*exp(-(xPrefactor*n/(x*1e-6)-v0)^2/(2*sigv^2))*dX*1e6)
			Reverse/P SingleOrder /D=SingleOrderNeg	
			//FofX +=  cn*1/VelNorm*(xPrefactor)^4*n^4/(x)^5*exp(-(xPrefactor*n/x*1e6-v0)^2/(2*sigv^2))*dX*1e24 
			//FofX +=  cn*1/VelNorm*(xPrefactor)^4*n^4/(-x)^5*exp(-(-xPrefactor*n/x*1e6-v0)^2/(2*sigv^2))*dX*1e24 
		
			Multithread		SingleOrderMol = n^4/(x)^5*exp(-((xPrefactor/2)*n/x-v0)^2/sv)/16
			//SingleOrderMol = 1/VelNorm*(xPrefactor/2)^4*n^4/(x)^5*exp(-((xPrefactor/2)*n/x*1e6-v0)^2/(2*sigv^2))*dX*1e24
			//SingleOrderMol = abs(1/VelNorm*xPrefactor/2*n/(x)^2*exp(-((xPrefactor/2)*n/(x*1e-6)-v0)^2/(2*sigv^2))*dX*1e6)
			Reverse/P SingleOrderMol /D=SingleOrderMolNeg
			//FofX +=  cn*cmol*1/VelNorm*(xPrefactor/2)^4*n^4/(x)^5*exp(-(xPrefactor/2*n/x*1e6-v0)^2/(2*sigv^2))*dX*1e24 
			//FofX +=  cn*cmol*1/VelNorm*(xPrefactor/2)^4*n^4/(-x)^5*exp(-(-xPrefactor/2*n/x*1e6-v0)^2/(2*sigv^2))*dX*1e24 
		
			MatrixOP/O/S FofX = FofX + cn* (  (SingleOrder+SingleOrderNeg) + cmol* (SingleOrderMol+SingleOrderMolNeg)  )
		
			//Print "time: ", StopMSTimer(timerRefNum)*10^-6
			n+=1
		while(n<=totalOrders)
		//	abort
		
		FofX*=isotopeFraction
	
		FofXname = "FofX"+num2str(i)
		Duplicate/O FofX $FofXname
	
		FofXsum += FofX
	
	endfor
	
	
	FofX = FofXsum
	
	MatrixOP/O/S FofX=FofX*NormConst
	
	//Add the 0th order atoms and molecules and scale the whole pattern
	FofX[x2pnt(FofX,0)] = 1+cmol

	
	//Generate the raw beam waves
	Variable rawScanHalfWidth = 1600
	Variable rawBeamPnts = rawScanHalfWidth*2/dX +1
	Make/D/O/N=(rawBeamPnts) rawY, rawX
	SetScale/P x -rawScanHalfWidth, dX, "um", rawY, rawX
	rawX = x
	
	// Choose beam shape here by turning on (1) one if/elseif and turning off (0) the rest
	if(0)	// pure convolution of 1s, 2s, and det
		Make/D/O/N=1 rawbeamparams={1, s1,s2,0,0}
		FitRawBeam(rawbeamparams, rawY, rawX)
	elseif(1) // convolution plus gaussian
		Make/D/O/N=1 rawbeamparams={1, s1, s2, 0, 0, GaussAmpRatio, GaussWidth}
		FitRawBeamGauss(rawbeamparams, rawY, rawX)
	elseif(0) // convolution plus gaussian
		Make/D/O/N=1 rawbeamparams={1, s1, s2, 0, 0, GaussAmpRatio, GaussWidth, s3}
		FitRawBeam3sGauss(rawbeamparams, rawY, rawX)
	elseif(0) // two gaussians
		rawY = sqrt(2*pi)*GaussWidth*(Gauss(x, 0, GaussWidth)+GaussAmpRatio*Gauss(x,0,GaussWidth*50))
	endif
	
	Variable rawYmax = wavemax(rawY)
	Variable rawYarea = area(rawY)
	rawY /= rawYarea
	
	MultiThread FofX*=amp//*rawYarea
	
	//convolve the infinitely thin beam/detector diffraction pattern with the beam shape
	Duplicate/O FofX FofXConvolved	
	Convolve/A rawY FofXConvolved
	
	//assign output
	yw = FofXConvolved(xw[p]-x0)+background	
End






Function FitDiffractionScanK(pw, yw, xw) : FitFunc
	Wave pw, yw, xw
	
	// initialize fittable variables
	Variable v0	= pw[0]			// flow velocity
	Variable sv = pw[1]				// velocity distribution
	Variable cmol = pw[2]/(1-pw[2])			// molecule fraction
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
	Variable VelNorm// = area(velprobs2,0,inf); print VelNorm
	VelNorm = sqrt(2*pi)*sv*v0*(v0^2+3*sv^2)//; print VelNorm
	Velprobs2/=VelNorm
	
	
	// set geometrical parameters
	Variable L1d = 2.3974				// dist to detector in m
	//Variable Lgd = 96.5*25.4/1000; L1d=Lgd	// override 1g-det dist with magic grating distance
	
	Variable m 
	Variable a = 100e-9						// grating period in microns
	Variable h = 6.62607E-34				// h in kg*m^2/s^2
	Variable xPreFactor 		// = 0.41613 m^2/s for sodium
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
	Variable frac
	Variable j=1
	Variable totalOrders = numpnts(pw)-OtherParams
	variable n=1
	Variable cn = pw[n+OtherParams-1] 
	do
	if(j==1)
	m=K39mass*amu2kg
	frac=K39frac
	elseif(j==2)
	m=K39mass*amu2kg
	frac=K40frac
	elseif(j==3)
	m=K41mass*amu2kg
	frac=K41frac
	else
	m=K39mass
	frac=0
	endif
	xPreFactor = L1d*h/m/a
	n=1
	do
		// Depending on the scaling of FofX and VelProbs the diffraction orders will not all have integrated unit probability automatically
		// so we compute each diffraction order separately to find the proper scaling.
//Variable timerRefNum = Startmstimer
		cn = pw[n+OtherParams-1]// ; print n, cn
		
		SingleOrder = 1/VelNorm*(xPrefactor)^4*n^4/(x)^5*exp(-(xPrefactor*n/(x*1e-6)-v0)^2/(2*sv^2))*dX*1e24
		//SingleOrder = abs(1/VelNorm*xPrefactor*n/(x)^2*exp(-(xPrefactor*n/(x*1e-6)-v0)^2/(2*sv^2))*dX*1e6)
		Reverse/P SingleOrder /D=SingleOrderNeg	
		//FofX +=  cn*1/VelNorm*(xPrefactor)^4*n^4/(x)^5*exp(-(xPrefactor*n/x*1e6-v0)^2/(2*sv^2))*dX*1e24 
		//FofX +=  cn*1/VelNorm*(xPrefactor)^4*n^4/(-x)^5*exp(-(-xPrefactor*n/x*1e6-v0)^2/(2*sv^2))*dX*1e24 
		
		SingleOrderMol = 1/VelNorm*(xPrefactor/2)^4*n^4/(x)^5*exp(-((xPrefactor/2)*n/x*1e6-v0)^2/(2*sv^2))*dX*1e24
		//SingleOrderMol = abs(1/VelNorm*xPrefactor/2*n/(x)^2*exp(-((xPrefactor/2)*n/(x*1e-6)-v0)^2/(2*sv^2))*dX*1e6)
		Reverse/P SingleOrderMol /D=SingleOrderMolNeg
		//FofX +=  cn*cmol*1/VelNorm*(xPrefactor/2)^4*n^4/(x)^5*exp(-(xPrefactor/2*n/x*1e6-v0)^2/(2*sv^2))*dX*1e24 
		//FofX +=  cn*cmol*1/VelNorm*(xPrefactor/2)^4*n^4/(-x)^5*exp(-(-xPrefactor/2*n/x*1e6-v0)^2/(2*sv^2))*dX*1e24 
		
		MatrixOP/O/S FofX = FofX+ frac* cn* (  (SingleOrder+SingleOrderNeg) + cmol* (SingleOrderMol+SingleOrderMolNeg)  )
//Print "time: ", StopMSTimer(timerRefNum)*10^-6
		n+=1
	while(n<=totalOrders)
	//print m, frac
	j+=1
	while(j<=3)
	
	//Add the 0th order atoms and molecules and scale the whole pattern
	FofX[x2pnt(FofX,0)] = 1+cmol

	
	//Generate the raw beam waves
	Variable rawScanHalfWidth = 1600
	Variable rawBeamPnts = rawScanHalfWidth*2/dX +1
	Make/D/O/N=(rawBeamPnts) rawY, rawX
	SetScale/P x -rawScanHalfWidth, dX, "um", rawY, rawX
	rawX = x
	
	// Choose beam shape
	if(0)	// pure convolution of 1s, 2s, and det
		Make/D/O/N=1 rawbeamparams={1, s1,s2,0,0}
		FitRawBeam(rawbeamparams, rawY, rawX)
	elseif(1) // convolution plus gaussian
		//Variable GaussRatio = 1/40
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











Window DiffParamsTable() : Table
	Make/T/O/N=20   KWave1={"K0","K1","K2","K3","K4","K5","K6","K7","K8","K9","K10","K11","K12","K13","K14","K15","K16","K17","K18","K19","K20","K21", "K22"}
	Make/T/O/N=20  ParamList1= {"flow velocity","sigma","molecule fraction","amplitude","center offset","background","1s width","2s width", "3s param", "Gauss amp ratio","Gauss width","1st order","2nd order","3rd order","4th order","5th order", "6th order","7th order","8th order","9th order","10th order","11th order","12th order"}
	Make/T/O/N=20     ConstraintsConv1={"","K1>K0/30","K2>.00001","","","","K6>5","K7>5","K8>5","K9>.000001","K10>10","K11>.1","","K13>.001","K14>.001","K15>.001","K16>.0001","K17>.0001","K18>.0001","K19>.0001","K20>.0001", "K21>0.0001","K22>0.0001"}
	Make/T/O/N=20  ConstraintsConv2={"","","K2<20","","","","K6<200","K7<200","K8<200","K9<.05","K10<1000","","","K13<.2","K14<.15","K15<.1","K16<.1","K17<.1","K18<.1","K19<.1","K20<.1","K21<.1","K22<.1"}
	Make/T/O/N=20   ConstraintsConv= {"K1>K0/30","K2>.0001","K6>20","K7>20","K8>.001","K9>10","K10>.1","K12>.001","K13>.001","K14>.001","K15>.0001","K16>.0001","K2<.2","K6<200","K7<200","K8<.2","K9<500","K12<.2","K13<.1","K14<.1","K15<.05","K16<.05"}
	Make/D/O/N=20  InitialGuesses1= {3000,150,0,35000,0,0.35,25,25,50,0.001,200,0.6,0.1,0.04,0.04,0.02,0.02,0.015,0.01,0.01,0.01,0.01,0.01}
	Make/D/O/N=20  W_coef= {1962.53,120,0,31421.8,-0.747397,0.35,55.6519,18.9738,57.9848,5.62764e-05,769.341,0.609617,0.0938087,0.0409004,0.0403856,0.0142732,0.0161595,0.0106368,0.00870566,0.00860006,0.00548227,0.01,0.01}
	Make/D/O/N=20 W_sigma={1.48604,0,0,929.777,0.0773076,0,16.719,8.64061,2.55465,1.01211e-05,404.541,0.00239707,0.00175734,0.00156781,0.00188425,0.0027069,0.00380549,0.00465392,0.00536845,0.00547534,0.00630387,0,0}
	Make/O   HoldConv= {0,0,1,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1}
	
	PauseUpdate; Silent 1		// building window...
	Edit/W=(141,44,1099,477) KWave1,ParamList1,InitialGuesses1,w_coef,w_sigma,HoldConv,ConstraintsConv1,ConstraintsConv2,ConstraintsConv
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
	NVAR m
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





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



Function FitMultipleScans(series, start, stop)
	String series
	Variable start
	Variable stop

	Make/O/N=(stop-start+1) velocities, sigmas, velunc, sigunc, MolFracs, MolFracsUnc, RchiSquared
	SetScale/P x 1,1, "file number", velocities, sigmas, velunc, sigunc, MolFracs, MolFracsUnc, RchiSquared
	
	Wave W_coef, W_sigma, velocities, sigmas, velunc, sigunc, MolFracs, MolFracsUnc
	
	Variable i
	For(i = start; i <= stop; i += 1)
	If(1)
		ProcessDiffraction(series + num2str(i))
		velocities[i-1]=W_coef[0]; sigmas[i-1]=W_coef[1]; velunc[i-1]=W_sigma[0]; sigunc[i-1]=W_sigma[1]; MolFracs[i-1]=W_coef[2]; MolFracsUnc[i-1]=W_sigma[2];
		NVAR reducedChiSqrd
		RchiSquared[i-1]= reducedChiSqrd
	EndIf
	EndFor
	
	DuplicateSeriesResults(series,1)
	
	VelSigMolChiGraph(series)
End




Function DuplicateSeriesResults(series, displayToggle)
	String series
	Variable displayToggle
	
	String name
	name = "velocities_" + series; Duplicate/O velocities $name
	name = "velunc_" + series; Duplicate/O velunc $name
	name = "sigmas_" + series; Duplicate/O sigmas $name
	name = "sigunc_" + series; Duplicate/O sigunc $name
	name = "MolFrac_"+series; Duplicate/O MolFracs $name
	name = "MolFracUnc_"+series; Duplicate/O MolFracsUnc $name
	name = "RchiSquared_"+series; Duplicate/O RchiSquared $name
End



//Function SigmasGraph(sigmasw, siguncw) : Graph
//	Wave sigmasw, siguncw
//	String sigName = NameOfWave(sigmasw); String sigUncName = NameOfWave(siguncw)
//	
//	Wavestats/Q sigmasw
//	Variable AvgSig = V_avg
//	Variable StdSig = V_sdev
//	Variable PercentBars = 5 // in units of percentage
//	Make/O/N=(numpnts(sigmas)) $sigName+"Plus1" = AvgSig*(1+PercentBars/100), $sigName+"Minus1" = AvgSig*(1-PercentBars/100)
//		
//	PauseUpdate; Silent 1		// building window...
//	Display/K=1 sigmasw,$sigName+"Plus1",$sigName+"Minus1"
//	ModifyGraph mode($sigName)=4
//	ModifyGraph marker($sigName)=8
//	ModifyGraph lStyle($sigName+"Plus1")=3,lStyle($sigName+"Minus1")=3
//	ModifyGraph rgb($sigName+"Plus1")=(0,0,0),rgb($sigName+"Minus1")=(0,0,0)
//	ModifyGraph grid(bottom)=1
//	ModifyGraph nticks(bottom)=8
//	ModifyGraph sep(bottom)=4
//	Label left "Velocity Distribution (m/s)"
//	Label bottom "File index"
//	ErrorBars $sigName Y,wave=($sigUncName,$sigUncName)
//	TextBox/C/N=text1 /A=RT "Sig_avg = " + num2str(AvgSig) + " ± " + num2str(StdSig)
//	TextBox/C/N=text2 /A=LT "Dashed Lines: ± " + num2str(PercentBars) + "%"
//	SetAxis bottom 1,numpnts(sigmasw)-1
//	
//	WindowNamer(sigName)
//	
//	MinusScansStats(sigmasw)
//End
//
//
//
//Function VelocitiesGraph(velocitiesw, veluncw) : Graph
//	Wave velocitiesw, veluncw
//	String velName = NameOfWave(velocitiesw); String velUncName = NameOfWave(veluncw)
//	
//	Wavestats/Q velocitiesw
//	Variable AvgVel = V_avg
//	Variable StdVel = V_sdev
//	Variable PercentBars = 0.5 // in units of percentage
//	Make/O/N=(numpnts(velocities)) $velName+"Plus1" = AvgVel*(1+PercentBars/100), $velName+"Minus1" = AvgVel*(1-PercentBars/100)
//	
//	PauseUpdate; Silent 1		// building window...
//	Display/K=1 velocitiesw,$velName+"Plus1",$velName+"Minus1"
//	ModifyGraph mode($velName)=4
//	ModifyGraph marker($velName)=8
//	ModifyGraph lStyle($velName+"Plus1")=3,lStyle($velName+"Minus1")=3
//	ModifyGraph rgb($velName+"Plus1")=(0,0,0),rgb($velName+"Minus1")=(0,0,0)
//	ModifyGraph grid(bottom)=1
//	ModifyGraph nticks(bottom)=8
//	ModifyGraph sep(bottom)=4
//	Label left "Velocity (m/s)"
//	Label bottom "File index"
//	ErrorBars $velName Y,wave=($velUncName,$velUncName)
//	TextBox/C/N=text1 /A=RT "V_avg = " + num2str(AvgVel) + " ± " + num2str(StdVel)
//	TextBox/C/N=text2 /A=LT "Dashed Lines: ± " + num2str(PercentBars) + "%"
//	SetAxis bottom 1,numpnts(velocitiesw)-1
//	
//	WindowNamer(velName)
//	
//	MinusScansStats(velocitiesw)
//End
//
//
//Function MolFracGraph(MolFracw, MolFracuncw) : Graph
//	Wave MolFracw, MolFracuncw
//	String velName = NameOfWave(MolFracw); String velUncName = NameOfWave(MolFracuncw)
//	
//	Wavestats/Q MolFracw
//	Variable AvgVel = V_avg
//	Variable StdVel = V_sdev
//	Variable PercentBars = 10 // in units of percentage
//	Make/O/N=(numpnts(MolFracs)) $velName+"Plus1" = AvgVel*(1+PercentBars/100), $velName+"Minus1" = AvgVel*(1-PercentBars/100)
//	
//	PauseUpdate; Silent 1		// building window...
//	Display/K=1 MolFracw,$velName+"Plus1",$velName+"Minus1"
//	ModifyGraph mode($velName)=4
//	ModifyGraph marker($velName)=8
//	ModifyGraph lStyle($velName+"Plus1")=3,lStyle($velName+"Minus1")=3
//	ModifyGraph rgb($velName+"Plus1")=(0,0,0),rgb($velName+"Minus1")=(0,0,0)
//	ModifyGraph grid(bottom)=1
//	ModifyGraph nticks(bottom)=8
//	ModifyGraph sep(bottom)=4
//	Label left "Molecule Fraction"
//	Label bottom "File index"
//	ErrorBars $velName Y,wave=($velUncName,$velUncName)
//	TextBox/C/N=text1 /A=RT "V_avg = " + num2str(AvgVel) + " ± " + num2str(StdVel)
//	TextBox/C/N=text2 /A=LT "Dashed Lines: ± " + num2str(PercentBars) + "%"
//	SetAxis bottom 1,numpnts(MolFracw)-1
//	
//	WindowNamer(velName)
//	
//	MinusScansStats(MolFracw)
//End


Function VelSigMolGraph(series) : Graph
	String series
	
	String velname = "velocities_" + series; Wave velocities =$velname
	String veluncname = "velunc_" + series; Wave velunc =$veluncname
	String signame = "sigmas_" + series; Wave sigmas =$signame
	String sigUncname = "sigunc_" + series; Wave sigunc =$sigUncname
	String molFracname = "MolFrac_"+series; Wave MolFrac =$molFracname
	String molFracUncname = "MolFracUnc_"+series; Wave MolFracUnc =$molFracUncname
	
	Wavestats/Q velocities
	Variable AvgVel = V_avg
	Variable StdVel = V_sdev/sqrt(V_npnts)
	Variable PercentBars = 0.5 // in units of percentage
	Make/O/N=(numpnts(velocities)) $velName+"Plus1" = AvgVel*(1+PercentBars/100), $velName+"Minus1" = AvgVel*(1-PercentBars/100)
	SetScale/P x 1, 1, "file number", $velName+"Plus1", $velName+"Minus1"
	
	PauseUpdate; Silent 1		// building window...
	Display/K=1/W=(100,100,550,600) velocities,$velName+"Plus1",$velName+"Minus1"
	ModifyGraph mode($velName)=4
	ModifyGraph marker($velName)=8
	ModifyGraph lStyle($velName+"Plus1")=3,lStyle($velName+"Minus1")=3
	ModifyGraph rgb($velName+"Plus1")=(0,0,0),rgb($velName+"Minus1")=(0,0,0)
	ModifyGraph grid(bottom)=1
	ModifyGraph nticks(bottom)=(numpnts(velocities))
	ModifyGraph sep(bottom)=4
	Label left "Velocity (m/s)"
	ErrorBars $velName Y,wave=($velUncName,$velUncName)
	
	String tboxText
	sprintf tboxText, "V avg = %.1f ± %.1f", AvgVel, StdVel
	TextBox/C/N=velAvg/Z=0/X=0/Y=5.00 /A=RT tboxText
	TextBox/C/N=velDashes/X=4/Y=5.00 /A=LT "Dashed Lines: ± " + num2str(PercentBars) + "%"
	
	Make/O/N=2 results
	MinusScansStats(velocities, results)
	Variable MinusAverage = results[0]
	Variable MinusStdDev = results[1]
	sprintf tboxText, "Minus vel avg = %.1f ± %.1f", MinusAverage, MinusStdDev
	TextBox/C/N=velAvgMinus/Z=0/X=0/Y=25.00 /A=RT tboxText
	
	
	Wavestats/Q sigmas
	Variable AvgSig = V_avg
	Variable StdSig = V_sdev/sqrt(V_npnts)
	PercentBars = 5 // in units of percentage
	Make/O/N=(numpnts(sigmas)) $sigName+"Plus1" = AvgSig*(1+PercentBars/100), $sigName+"Minus1" = AvgSig*(1-PercentBars/100)
	SetScale/P x 1, 1, "file number", $sigName+"Plus1", $sigName+"Minus1"
		
	AppendToGraph/L=sigAxis sigmas,$sigName+"Plus1",$sigName+"Minus1"
	ModifyGraph mode($sigName)=4
	ModifyGraph marker($sigName)=8
	ModifyGraph lStyle($sigName+"Plus1")=3,lStyle($sigName+"Minus1")=3
	ModifyGraph rgb($sigName+"Plus1")=(0,0,0),rgb($sigName+"Minus1")=(0,0,0)
	Label sigAxis "Velocity Distribution (m/s)"
	ErrorBars $sigName Y,wave=($sigUncName,$sigUncName)
	
	sprintf tboxText, "Sig avg = %.1f ± %.1f", AvgSig, StdSig
	TextBox/C/N=sigAvg/Z=0/X=0/Y=42.00 /A=RT tboxText
	TextBox/C/N=sigDashes/Z=0/X=4/Y=42.00 /A=LT "Dashed Lines: ± " + num2str(PercentBars) + "%"
	
	MinusScansStats(sigmas, results)
	 MinusAverage = results[0]
	 MinusStdDev = results[1]
	 sprintf tboxText, "Minus sig avg = %.1f ± %.1f", MinusAverage, MinusStdDev
	TextBox/C/N=sigAvgMinus/Z=0/X=0/Y=60.00 /A=RT tboxText
	
	
	Wavestats/Q MolFrac
	Variable AvgMol = V_avg
	Variable StdMol = V_sdev/sqrt(V_npnts)
	PercentBars = 5 // in units of percentage
	Make/O/N=(numpnts(MolFracs)) $molFracName+"Plus1" = AvgMol*(1+PercentBars/100), $molFracName+"Minus1" = AvgMol*(1-PercentBars/100)
	SetScale/P x 1, 1, "file number", $molFracName+"Plus1", $molFracName+"Minus1"
	
	AppendToGraph/L=molAxis MolFrac,$molFracName+"Plus1",$molFracName+"Minus1"
	ModifyGraph mode($molFracName)=4
	ModifyGraph marker($molFracName)=8
	ModifyGraph lStyle($molFracName+"Plus1")=3,lStyle($molFracName+"Minus1")=3
	ModifyGraph rgb($molFracName+"Plus1")=(0,0,0),rgb($molFracName+"Minus1")=(0,0,0)
	Label molAxis "Molecule Fraction (%)"
	ErrorBars $molFracName Y,wave=($MolFracUncName,$MolFracUncName)
	
	sprintf tboxText, "Mol avg = %.1f ± %.1f", AvgMol, StdMol
	TextBox/C/N=molAvg/Z=0/X=0/Y=72.00 /A=RT tboxText
	TextBox/C/N=molDashes/Z=0/X=4/Y=72.00 /A=LT  "Dashed Lines: ± " + num2str(PercentBars) + "%"
	
	MinusScansStats(molFrac, results)
	 MinusAverage = results[0]
	 MinusStdDev = results[1]
	  sprintf tboxText, "Minus mol avg = %.1f ± %.1f", MinusAverage, MinusStdDev
	TextBox/C/N=molAvgMinus/Z=0/X=0/Y=95.00 /A=RT tboxText
	
	
	
	ModifyGraph axisEnab(left)={0.7,1},axisEnab(sigAxis)={0.35,0.65},axisEnab(molAxis)={0,0.3}
	ModifyGraph freePos(sigAxis)={0,bottom}, freePos(molAxis)={0,bottom}
	ModifyGraph lblPosMode(sigAxis)=1, lblPosMode(molAxis)=1
	SetAxis bottom 0.8,*
	
	ModifyGraph gfSize=10
	
	WindowNamer("Diff Scan Stats "+series)
End



Function VelSigMolChiGraph(series) : Graph
	String series
	
	String velname = "velocities_" + series; Wave velocities =$velname
	String veluncname = "velunc_" + series; Wave velunc =$veluncname
	String signame = "sigmas_" + series; Wave sigmas =$signame
	String sigUncname = "sigunc_" + series; Wave sigunc =$sigUncname
	String molFracname = "MolFrac_"+series; Wave MolFrac =$molFracname
	String molFracUncname = "MolFracUnc_"+series; Wave MolFracUnc =$molFracUncname
	String RchiSquaredname = "RchiSquared_"+series; Wave RchiSquared = $RchiSquaredname
	
	Wavestats/Q velocities
	Variable AvgVel = V_avg
	Variable StdVel = V_sdev/sqrt(V_npnts)
	Variable PercentBars = 0.5 // in units of percentage
	Make/O/N=(numpnts(velocities)) $velName+"Plus1" = AvgVel*(1+PercentBars/100), $velName+"Minus1" = AvgVel*(1-PercentBars/100)
	SetScale/P x 1, 1, "file number", $velName+"Plus1", $velName+"Minus1"
	
	PauseUpdate; Silent 1		// building window...
	Display/K=1/W=(100,100,550,700) velocities,$velName+"Plus1",$velName+"Minus1"
	ModifyGraph mode($velName)=4
	ModifyGraph marker($velName)=8
	ModifyGraph lStyle($velName+"Plus1")=3,lStyle($velName+"Minus1")=3
	ModifyGraph rgb($velName+"Plus1")=(0,0,0),rgb($velName+"Minus1")=(0,0,0)
	ModifyGraph grid(bottom)=1
	ModifyGraph nticks(bottom)=(numpnts(velocities))
	ModifyGraph sep(bottom)=4
	Label left "Velocity (m/s)"
	ErrorBars $velName Y,wave=($velUncName,$velUncName)
	
	String tboxText
	sprintf tboxText, "V avg = %.1f ± %.1f", AvgVel, StdVel
	TextBox/C/N=velAvg/Z=0/X=0/Y=25.00 /A=RT tboxText
	TextBox/C/N=velDashes/X=4/Y=25.00 /A=LT "Dashed Lines: ± " + num2str(PercentBars) + "%"
	
	Make/O/N=2 results
	MinusScansStats(velocities, results)
	Variable MinusAverage = results[0]
	Variable MinusStdDev = results[1]
	sprintf tboxText, "Minus vel avg = %.1f ± %.1f", MinusAverage, MinusStdDev
	TextBox/C/N=velAvgMinus/Z=0/X=0/Y=40.00 /A=RT tboxText
	
	PlusScansStats(velocities, results)
	Variable PlusAverage = results[0]
	Variable PlusStdDev = results[1]
	sprintf tboxText, "Plus vel avg = %.1f ± %.1f", PlusAverage, PlusStdDev
	TextBox/C/N=velAvgPlus/Z=0/X=0/Y=30.00 /A=RT tboxText	
	
	
	Wavestats/Q rChiSquared
	Variable AvgrChiSquared = V_avg
	Variable StdrChiSquared = V_sdev/sqrt(V_npnts)
	//PercentBars = 5 // in units of percentage
	//Make/O/N=(numpnts(sigmas)) $sigName+"Plus1" = AvgSig*(1+PercentBars/100), $sigName+"Minus1" = AvgSig*(1-PercentBars/100)
	//SetScale/P x 1, 1, "file number", $sigName+"Plus1", $sigName+"Minus1"
		
	AppendToGraph/L=chiAxis rChiSquared//,$sigName+"Plus1",$sigName+"Minus1"
	ModifyGraph mode($rChiSquaredName)=4
	ModifyGraph marker($rChiSquaredName)=8
	//ModifyGraph lStyle($sigName+"Plus1")=3,lStyle($sigName+"Minus1")=3
	//ModifyGraph rgb($sigName+"Plus1")=(0,0,0),rgb($sigName+"Minus1")=(0,0,0)
	Label chiAxis "Red. Chi-Squared"
	//ErrorBars $sigName Y,wave=($sigUncName,$sigUncName)
	
	sprintf tboxText, "Red Chi^2 avg = %.3f ± %.3f", AvgrChiSquared, StdrChiSquared
	TextBox/C/N=chiAvg/Z=0/X=0/Y=4.00 /A=RT tboxText
	//TextBox/C/N=sigDashes/Z=0/X=4/Y=42.00 /A=LT "Dashed Lines: ± " + num2str(PercentBars) + "%"
	
	MinusScansStats(rChiSquared, results)
	 MinusAverage = results[0]
	 MinusStdDev = results[1]
	 sprintf tboxText, "Minus Red Chi^2 avg = %.3f ± %.3f", MinusAverage, MinusStdDev
	TextBox/C/N=chiAvgMinus/Z=0/X=0/Y=20.00 /A=RT tboxText

	PlusScansStats(rChiSquared, results)
	 PlusAverage = results[0]
	 PlusStdDev = results[1]
	 sprintf tboxText, "Plus Red Chi^2 avg = %.3f ± %.3f", PlusAverage, PlusStdDev
	TextBox/C/N=chiAvgPlus/Z=0/X=0/Y=9.00 /A=RT tboxText
	
	
	
	Wavestats/Q sigmas
	Variable AvgSig = V_avg
	Variable StdSig = V_sdev/sqrt(V_npnts)
	PercentBars = 5 // in units of percentage
	Make/O/N=(numpnts(sigmas)) $sigName+"Plus1" = AvgSig*(1+PercentBars/100), $sigName+"Minus1" = AvgSig*(1-PercentBars/100)
	SetScale/P x 1, 1, "file number", $sigName+"Plus1", $sigName+"Minus1"
		
	AppendToGraph/L=sigAxis sigmas,$sigName+"Plus1",$sigName+"Minus1"
	ModifyGraph mode($sigName)=4
	ModifyGraph marker($sigName)=8
	ModifyGraph lStyle($sigName+"Plus1")=3,lStyle($sigName+"Minus1")=3
	ModifyGraph rgb($sigName+"Plus1")=(0,0,0),rgb($sigName+"Minus1")=(0,0,0)
	Label sigAxis "Velocity Distribution (m/s)"
	ErrorBars $sigName Y,wave=($sigUncName,$sigUncName)
	
	sprintf tboxText, "Sig avg = %.1f ± %.1f", AvgSig, StdSig
	TextBox/C/N=sigAvg/Z=0/X=0/Y=52.00 /A=RT tboxText
	TextBox/C/N=sigDashes/Z=0/X=4/Y=52.00 /A=LT "Dashed Lines: ± " + num2str(PercentBars) + "%"
	
	MinusScansStats(sigmas, results)
	 MinusAverage = results[0]
	 MinusStdDev = results[1]
	 sprintf tboxText, "Minus sig avg = %.1f ± %.1f", MinusAverage, MinusStdDev
	TextBox/C/N=sigAvgMinus/Z=0/X=0/Y=70.00 /A=RT tboxText

	PlusScansStats(sigmas, results)
	 PlusAverage = results[0]
	 PlusStdDev = results[1]
	 sprintf tboxText, "Plus sig avg = %.1f ± %.1f", PlusAverage, PlusStdDev
	TextBox/C/N=sigAvgPlus/Z=0/X=0/Y=57.00 /A=RT tboxText
	
	
	Wavestats/Q MolFrac
	Variable AvgMol = V_avg
	Variable StdMol = V_sdev/sqrt(V_npnts)
	PercentBars = 5 // in units of percentage
	Make/O/N=(numpnts(MolFracs)) $molFracName+"Plus1" = AvgMol*(1+PercentBars/100), $molFracName+"Minus1" = AvgMol*(1-PercentBars/100)
	SetScale/P x 1, 1, "file number", $molFracName+"Plus1", $molFracName+"Minus1"
	
	AppendToGraph/L=molAxis MolFrac,$molFracName+"Plus1",$molFracName+"Minus1"
	ModifyGraph mode($molFracName)=4
	ModifyGraph marker($molFracName)=8
	ModifyGraph lStyle($molFracName+"Plus1")=3,lStyle($molFracName+"Minus1")=3
	ModifyGraph rgb($molFracName+"Plus1")=(0,0,0),rgb($molFracName+"Minus1")=(0,0,0)
	Label molAxis "Molecule Fraction (%)"
	ErrorBars $molFracName Y,wave=($MolFracUncName,$MolFracUncName)
	
	sprintf tboxText, "Mol avg = %.1f ± %.1f", AvgMol, StdMol
	TextBox/C/N=molAvg/Z=0/X=0/Y=76.00 /A=RT tboxText
	TextBox/C/N=molDashes/Z=0/X=4/Y=76.00 /A=LT  "Dashed Lines: ± " + num2str(PercentBars) + "%"
	
	MinusScansStats(molFrac, results)
	 MinusAverage = results[0]
	 MinusStdDev = results[1]
	  sprintf tboxText, "Minus mol avg = %.1f ± %.1f", MinusAverage, MinusStdDev
	TextBox/C/N=molAvgMinus/Z=0/X=0/Y=95.00 /A=RT tboxText
	
	PlusScansStats(molFrac, results)
	 PlusAverage = results[0]
	 PlusStdDev = results[1]
	  sprintf tboxText, "Plus mol avg = %.1f ± %.1f", PlusAverage, PlusStdDev
	TextBox/C/N=molAvgPlus/Z=0/X=0/Y=81.00 /A=RT tboxText
	
	
	
	ModifyGraph axisEnab(left)={0.51,0.74},axisEnab(sigAxis)={0.26,0.49},axisEnab(molAxis)={0,0.24}, axisEnab(chiAxis)={0.76,1}
	ModifyGraph freePos(sigAxis)={0,bottom}, freePos(molAxis)={0,bottom}, freePos(chiAxis)={0,bottom}
	ModifyGraph lblPosMode(sigAxis)=1, lblPosMode(molAxis)=1, lblPosMode(chiAxis)=1
	SetAxis bottom 0.8,*
	
	ModifyGraph gfSize=10
	
	WindowNamer("Diff Scan Stats "+series)
End





Function MinusScansStats(wavein, waveout)
		Wave wavein, waveout
		
		Duplicate/O wavein tempWave
		
		tempWave/=(Mod(x,2)==0)
		
		WaveStats/Q tempWave
		
		waveout[0] = V_avg
		waveout[1] = V_sdev/sqrt(V_npnts)
End


Function PlusScansStats(wavein, waveout)
		Wave wavein, waveout
		
		Duplicate/O wavein tempWave
		
		tempWave/=(Mod(x,2)!=0)
		
		WaveStats/Q tempWave
		
		waveout[0] = V_avg
		waveout[1] = V_sdev/sqrt(V_npnts)
End





Function KillFileResults(series, fileNumber)
	String series
	Variable fileNumber
	
	String velname = "velocities_" + series; Wave velocities =$velname
	String veluncname = "velunc_" + series; Wave velunc =$veluncname
	String signame = "sigmas_" + series; Wave sigmas =$signame
	String sigUncname = "sigunc_" + series; Wave sigunc =$sigUncname
	String molFracname = "MolFrac_"+series; Wave MolFrac =$molFracname
	String molFracUncname = "MolFracUnc_"+series; Wave MolFracUnc =$molFracUncname
	String chiName = "rChiSquared_"+series; Wave chi = $chiName
	
	velocities[fileNumber-1]=NaN; sigmas[fileNumber-1]=NaN; MolFrac[fileNumber-1]=NaN; chi[fileNumber-1]=NaN
	
	VelSigMolChiGraph(series)
End





Function AppendResults(params)
	Wave/T params
	Wave w_coef, w_sigma
	
	String text = "", textTotal = "\\Zr100"
	Variable n = 0
	do
		sprintf text, "  %s = %g ± %g  ", params[n], w_coef[n], w_sigma[n]
		print text
		textTotal += text
		if(n<numpnts(w_coef)-1)
			textTotal+="\r"
		endif
		n+=1
	while(n<numpnts(w_coef))
	//print textTotal
	TextBox/C/N=ResultsTBox /A=RT textTotal
End






Function BatchProcess()
	Wave/T SeriesName
	Wave BackgroundWave, SigmaWave
	
	Wave InitialGuesses1
	
	Variable numSeries = numpnts(SeriesName)
	Variable i=0
	
	For(i=0; i<numpnts(SeriesName); i+=1)
		InitialGuesses1[5]=BackgroundWave[i]
		InitialGuesses1[1]=SigmaWave[i]
		FitMultipleScans(SeriesName[i],1,8)
	EndFor
	
End


Function BatchProcess2()
	Wave/T SeriesName
	Wave BackgroundWave, SigmaWave
	
	Wave InitialGuesses1
	
	Variable numSeries = numpnts(SeriesName)
	Variable i=0
	
	For(i=0; i<numpnts(SeriesName); i+=1)
		InitialGuesses1[5]=BackgroundWave[i]
		InitialGuesses1[1]=SigmaWave[i]
		FitMultipleScans(SeriesName[i],1,8)
	EndFor
	
End