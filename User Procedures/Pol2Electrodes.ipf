#pragma rtGlobals=1		// Use modern global access method.

#include ":PhysicalConstants"
#include ":LabConstants"
#include ":WindowNamer"
#include ":RemoveNaNs"

/////////////////////////////////////////////////////////////////////////////////////////
//
//	PredictPhase2Electrodes(alpha, xoffset, phaseOffset)
//		Wrapper function for calculating a predicted (or best fit) phase shift vs prediction.
//
//	PhaseFit2Electrodes(SeriesName, OneKTwoK, holdOffsets)
//		Automated function for running the fit routine and displaying nice results.
//		SeriesName (string) sets the series to analyze		
//		OneKTwoK (string) controls if you should process the 1k2k or 2k2k data
//		holdOffsets (string) controls holding of the x offset or phi offset parameters in the fit.
//
//	ExtractPhaseContrast2(SeriesName, oneKtwoK, makegraphsYesNo)
//		Does the messy work of extracting useful phase and contrast data from the FringeFit results and plotting it.	
//
//	PolParams2()
//		Displays prompts to set the relavent parameters (velocity, atom, etc) for polarizabilty fits for this data folder.
//	
//	Phi2Electrodes(pw, output, y)
// 		Does the real work of calculating phase shift as a function of position
//	
//	PhiMonochromatic(positions)
//		Calculates phase shift / polarizability for v0. 
//		After multiplying by input or fit polarizability, this is used for unwrapping the more sophisticated phase shift calculation.
//
//	AppendFittedPolWave2electrodes(PositionWave, PhaseWave)
//		Append fitted wave to graph if it's not already there.
//	
//	BatchGraph(appendedName)
//		Takes the average of a data series and makes a nice plot of results.	
//
//	BatchDuplicate(appendedName)
//		Useful for backing up or saving the results of many data sets before recalculating using different parameters.
//
//	FindKillablesPolv2(series_name, sigmas)
//		Finds points to ignore based on residual size.
//
//	AutoFitAllDF()
//		Run fitting routine on all data series in data folders formatted as e.g. '130726a', '130726b', etc...
//
//	CalculateConstants2Electrodes()
//		Calculates some of the parameters needed in the polarizability calculation.
//
//	SagnacAndGravityPhase()
//		Calculate Sagnac and gravitation phase shifts as a function of velocity
//
//	NormalizeVelocitiesAn(velocity, sigma, n)
//		Prepare the desired velocity distribution.
//
//	PolPhaseAn()
//		Calculate the position and velocity dependent phase shift / pol.
//		This only needs to be done once and then the fit routine can multiply the results by polarizability on each iteration.
//		
//	CleanupPolWaves()
//		Clean up all the big unneeded waves after the calculation is done.
//
/////////////////////////////////////////////////////////////////////////////////////////


Static Constant number_v = 100	// points for velocity distribution integral
Static Constant n = 3				// n= 3 for M-B velocity distribution or 0 for gaussian


// Calculate a predicted (or best fit) phase shift vs prediction.
Function PredictPhase2Electrodes(alpha, xoffset, phaseOffset)
	Variable alpha
	Variable xoffset
	Variable phaseOffset
	
	NVAR velocity, vratio
	Variable sigma = velocity/vratio
	
	Variable timerStart = StopMStimer(-1)
	
	Make/D/O paramWave = {alpha, xoffset, phaseOffset}
	
	Variable positionOffset = xoffset//.0001	//in microns
	Make/D/O/N=80 position = -2000+50 * x+positionOffset, predictedPhase	//in microns
	
	
	Variable/G num_pos = numpnts(position)
	MatrixOP/O position_shifted =  position					// account for y offset
	Make/D/o/n=(num_pos, number_v) yy    						
	yy = position_shifted[p]	

	Make/o/D/n=(num_pos) IFMphase, IFMprob, Contrast, yshifted	

	CalculateConstants2Electrodes()
	
	NormalizeVelocitiesAn(velocity, sigma, n)
	
	PhiMonochromatic(position_shifted)

	PolPhaseAn()

	SagnacAndGravityPhase()
	
	Phi2Electrodes(paramWave, predictedPhase, position)
	
	Duplicate/O Contrast, predictedContrast
	
	CleanupPolWaves()
	
	//Print "New: ", (timerstart-StopMSTimer(-1))*10^-6
End


// Calculates phase shift / polarizability for v0. This is used for unwrapping the phase shift.
Function PhiMonochromatic(positions)
	Wave positions
	
	NVAR velocity, d, s
	
	Duplicate/O positions phaseOverPolMono
	
	phaseOverPolMono = -1/velocity*(1/(positions^2-d^2)-1/((positions-s/velocity)^2-d^2))
End





Function FindKillablesPol(series_name, sigmas)
	String series_name
	Variable sigmas
	
	String traces = TraceNameList("", ";", 1)							// Get all traces on the graph
	String resTraceName = StringFromList(0,GrepList(traces, "Res"))	// Find the name of the trace that starts with "Res"
	
	Wave resWave = $resTraceName
	
	wavestats/q resWave
	
	Duplicate/O resWave killables
	
	killables = abs(resWave) > sigmas*v_sdev
	
	String/G autoKill = ""
	String/G autoKillc = ""
	String/G autoKillfiles = ""
	String/G autoKillfilesSemi=""
	
	String filesWaveName = "file_index_"+series_name+"_red"
	Wave files = $filesWaveName
	
	Variable i
	for(i=0; i<numpnts(killables); i+=1)
		if(killables[i] == 1)
			//print i
			autoKill += num2str(i)+";"
			autoKillc += num2str(i)+","
			autoKillfiles += num2str(files[i])+"," 
			autoKillfilesSemi += num2str(files[i])+";" 
		endif
	endfor
	
	//print autokill
	print "points: ", autokillc
	print "files: ", autokillfiles
	print "files: ", autokillfilesSemi
End



// Finds points to ignore based on residual size.
Function FindKillablesPolv2(series_name, sigmas)
	String series_name
	Variable sigmas
	
	String traces = TraceNameList("", ";", 1)							// Get all traces on the graph
	String resTraceName = StringFromList(0,GrepList(traces, "Res"))	// Find the name of the trace that starts with "Res"
	
	Wave resWave = $resTraceName
	
	resWave = resWave == 0 ? NaN : resWave
	
	wavestats/q resWave
	
	Duplicate/O resWave killables
	
	killables = abs(resWave) > sigmas*v_sdev
	
	String/G autoKill = ""
	String/G autoKillc = ""
	String/G autoKillfiles = ""
	String/G autoKillfilesSemi=""
	
	String filesWaveName = "file_index_"+series_name+"_red"
	Wave files = $filesWaveName
	
	Variable i
	for(i=0; i<numpnts(killables); i+=1)
		if(killables[i] == 1)
			//print i
			autoKill += num2str(i)+";"
			autoKillc += num2str(i)+","
			autoKillfiles += num2str(files[i])+"," 
			autoKillfilesSemi += num2str(files[i])+";" 
		endif
	endfor
	
	if(waveexists(AutoFilesToKill))
		Duplicate/O AutoFilesToKill AutoFilesToKillOld
	else
		Make/O/N=0 AutoFilesToKill, AutoFilesToKillOld
	endif
	
	Variable NumAutoFilesToKillNew = ItemsInList(autoKillfilesSemi)					//take a string of numbers seperated by ; and count how many there are
	if(NumAutoFilesToKillNew != 0)
		Make/O/N=(NumAutoFilesToKillNew) AutoFilesToKillNew
		AutoFilesToKillNew = Str2Num(StringFromList(p, autoKillfilesSemi))
		Concatenate/NP/O {AutoFilesToKillNew, AutoFilesToKillOld}, AutoFilesToKill
		//Variable/G AutoKillFilesYesNo = 1
	else
		//Make/O/N=0 AutoFilesToKill = NaN
		//Variable/G AutoKillFilesYesNo = 0
	EndIf
	
	//print autokill
	print "points: ", autokillc
	print "files: ", autokillfiles
	print "files: ", autokillfilesSemi
	
	print autofilestokill
End




//		Automated function for running the fit routine and displaying nice results.
//		SeriesName (string) sets the series to analyze		
//		OneKTwoK (string) controls if you should process the 1k2k or 2k2k data
//		holdOffsets (string) controls holding of the x offset or phi offset parameters in the fit.
Function/Wave PhaseFit2Electrodes(SeriesName, OneKTwoK, holdOffsets)
	String SeriesName
	Variable OnekTwok
	String holdOffsets
	
	String onektwokStr = ""
	if(OnekTwok==1)
		onektwokStr = "1k2k"
	elseif(OnekTwok == 2)
		onektwokStr = "2k2k"
	endif
	
	
	String posWaveName = "pos"+"_"+SeriesName+"_polfit_red"; Wave PositionWave = $posWaveName
	String phaseWaveName = "phase"+onektwokStr+"_"+SeriesName+"_polfit_red"; Wave PhaseWave = $phaseWaveName
	String phaseErrorWaveName = "phase_error"+onektwokStr+"_"+SeriesName+"_polfit_red"; Wave ErrorWave = $phaseErrorWaveName
	
	NVAR velocity, vratio
	variable sigma = velocity/vratio
	
	Variable starttime = StopMStimer(-2)
	
	Variable/G num_pos = numpnts(PositionWave)
	Make/D/o/n=(num_pos, number_v) yy = positionWave[p]
	Make/D/o/n=(num_pos) IFMphase, IFMprob, Contrast		// Waves for: total phase shift of monochrom. beam; phase shift of polychromatic beam
	Make/D/O/n=(num_pos) yshifted
	
	CalculateConstants2Electrodes()
	
	NormalizeVelocitiesAn(velocity, sigma, n)

	SagnacAndGravityPhase()
	
	
	String resWaveName = "Res_"+NameOfWave(PhaseWave)
	if(WaveExists($resWaveName))
		Wave resWave = $resWaveName
		resWave = NaN
	endif
	
	String holdStr = "0" + holdOffsets
	
	NVAR initpol
	Variable LastFitPointLoc = NumVarOrDefault("LastFitPoint", 1000)
	Make/D/O W_coef={initpol+.5, 0, 0} 					// Initial Guess
	Make/D/O Epsilonwave={1e-3, 1e-6, 1e-6}			// Amount to vary the fit parameter by in each iteration	
	Make/T/O Constraints = {"K0>20", "K0<70"}
	FuncFit/Q/N/M=2/TBOX=768/H=(holdStr) Phi2Electrodes W_coef PhaseWave[0,LastFitPointLoc]  /X=PositionWave /W=ErrorWave /I=1 /R /C=Constraints /E=EpsilonWave /M=mask //// 
	// add /O to calculate initial guesses only
	
	Wave w_sigma
	
	Variable reducedChiSqrd =  V_chisq/(V_npnts-(numpnts(W_coef))
	String rChiSqrdStr; sprintf rChiSqrdStr, "%.2f", reducedChiSqrd; printf "Reduced Chi-Squared: %.2f\r", reducedChiSqrd
	TextBox/C/N=ChiBox/Z=0/X=82.00/Y=30.00 "\\F'Symbol'c\\F'Geneva'\\S2\\M/dof = "+rChiSqrdStr

	Wave M_Covar
	Duplicate/O M_Covar, CorMat	 
	CorMat = M_Covar[p][q]/sqrt(M_Covar[p][p]*M_Covar[q][q])
	
	//String contrastWaveName=ReplaceString("phase", NameOfWave(PhaseWave), "contrast")
	//Wavestats/Q $contrastWaveName
	//Duplicate/O Contrast $contrastWaveName; Wave contrastWave = $contrastWaveName
	
	AppendFittedPolWave2electrodes(PositionWave, PhaseWave)				

	CleanupPolWaves()
	
	Variable/G bestFitPol = w_coef[0]
	Variable/G bestFitx0 = w_coef[1]
	Variable/G bestFitphi0 = w_coef[2]
	
	Variable/G bestFitPolSig = w_sigma[0]
	Variable/G bestFitx0Sig = w_sigma[1]
	Variable/G bestFitphi0Sig = w_sigma[2]
	
	Variable/G chisqrdDOF = reducedChiSqrd

	//Print "New fit: ", (StopMSTimer(-2)-startTime)*10^-6
	
	Make/O/D bestFitParams =  {bestFitPol, bestFitPolSig, bestFitx0, bestFitx0Sig, bestFitphi0, bestFitphi0Sig, chisqrdDOF}
	
	printf "alpha = %.3f ± %.3f A^3\r", bestFitPol, bestFitPolSig
	printf "x0 = %.1f ± %.1f um\r", bestFitx0, bestFitx0Sig
	printf "phi0 = %.0f ± %.0f mrad\r", bestFitPhi0*1000, bestFitPhi0Sig*1000
	
	Return bestFitParams
End



// Append fitted wave to graph if it's not already there.
Function AppendFittedPolWave2electrodes(PositionWave, PhaseWave)
	Wave PositionWave, PhaseWave
	Wave W_coef
	
	PredictPhase2electrodes(W_coef[0], w_coef[1], w_coef[2])			// Run the predict phase routine with the fitted polarizability and fixed offset
	
	Wave position, predictedPhase, predictedContrast
		
	String fitPosWaveName= NameOfWave(positionWave)+"_fit"			// make some descriptive wave names
	String fitPhaseWaveName = NameOfWave(phaseWave) +"_fit"
	String fitContrastWaveName = ReplaceString("phase", fitPhaseWaveName, "contrast")//NameOfWave(contrastWave) + "_fit"
	
	Duplicate/O position $fitPosWaveName; Wave fitPosWave = $fitPosWaveName					// assign the calculated position and phase waves to the right wave names
	Duplicate/O predictedPhase $fitPhaseWaveName; Wave fitPhaseWave = $fitPhaseWaveName
	Duplicate/O predictedContrast $fitContrastWaveName; Wave fitContrastWave = $fitContrastWaveName
	
	//fitposwave*=1e6	// calculations were done in meters, but this will make them display in microns
	
	If(StringMatch(TraceInfo("", fitPhaseWaveName,0),""))		// if the trace is not already on the graph then put it there
		AppendToGraph fitPhaseWave vs fitPosWave
		ModifyGraph rgb($fitPhaseWaveName)=(1,16019,65535)
		ModifyGraph zero(Res_Left)=1
		SetAxis/A Res_Left
		
		AppendToGraph/l=contrast fitContrastWave vs fitPosWave
	EndIf
	
//	If(StringMatch(TraceInfo("", fitcontrastWaveName,0),""))		// if the trace is not already on the graph then put it there
//		AppendToGraph/l=contrast fitcontrastWave vs fitPosWave
//		ModifyGraph rgb($fitcontrastWaveName)=(1,16019,65535)
//	EndIf
End



// Calculates phase as a function of position
Function Phi2Electrodes(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
	
	Variable alpha = pw[0]
	Variable xoffset = pw[1]
	Variable phaseoffset = pw[2]
	
	Wave yy, yshifted
	yshifted = (y - xoffset )*1e-6 //convert microns to meters
	yy = yshifted[p]
//	MatrixOP/O yy = y + xoffset	
	PolPhaseAn()									// calculate phase/polarizability
	PhiMonochromatic(yshifted)
	
	Wave IFMprob, density2D, Contrast, IFMphase_v_overpolM, IFMphase_v_overpolP, IFMphase_v_overpolMmol, IFMphase_v_overpolPmol, phaseOverPolMono, IFMphase_v_overpolM2, IFMphase_v_overpolP2
	NVAR molFrac, molPol, AnalyticPrefactor
	
	Variable molFracLoc = molFrac/100
	
	Variable pol = alpha*AnalyticPrefactor// pw[0] is alpha in cgs units but without the factor of 10^-24
	Variable polmol = molPol*AnalyticPrefactor
	//print molpol
	
	Wave phiSagGravLayered2D
	Variable sagPhaseSign = 1
	MatrixOP/O IFMphase_vM = pol * IFMphase_v_overpolM+sagPhaseSign*phiSagGravLayered2D			// now we have the actual phase shift to take the real and imaginary parts of
	MatrixOP/O IFMphase_vP = pol * IFMphase_v_overpolP+sagPhaseSign*phiSagGravLayered2D
	MatrixOP/O IFMphase_vMmol =  polmol * IFMphase_v_overpolMmol+sagPhaseSign*phiSagGravLayered2D
	MatrixOP/O IFMphase_vPmol = polmol * IFMphase_v_overpolPmol+sagPhaseSign*phiSagGravLayered2D
	MatrixOP/O IFMphase_vM2 = pol * IFMphase_v_overpolM2+sagPhaseSign*phiSagGravLayered2D
	MatrixOP/O IFMphase_vP2 = pol * IFMphase_v_overpolP2+sagPhaseSign*phiSagGravLayered2D
	
	Variable sndOrder = 0 //weight of 2nd order IFMs
	
	//0.5 to account for two IFMs										
	MatrixOP/O/s imagpart = 0.5  * density2D * ( (1-sndOrder) *	((1-molFracLoc)*(sin(IFMphase_vM) + sin(IFMphase_vP))    +    molFracLoc*(sin(IFMphase_vMmol) + sin(IFMphase_vPmol))) + sndOrder*(sin(IFMphase_vM2) + sin(IFMphase_vP2))	 )			// weight imaginary part by prob(v) 
	MatrixOP/O/s realpart =  0.5  * density2D * ( (1-sndOrder) *	((1-molFracLoc)*(cos(IFMphase_vM) + cos(IFMphase_vP))    +    molFracLoc*(cos(IFMphase_vMmol) + cos(IFMphase_vPmol))) + sndOrder*(cos(IFMphase_vM2) + cos(IFMphase_vP2))	 )				// weight real part by prob(v)
	Integrate/DIM=1 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][number_v-1]) + magsqr(realpart[p][number_v-1]) )
	IFMprob = atan2(imagpart[p][number_v-1], realpart[p][number_v-1])
	MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
//	Variable n = NumVarOrDefault("k", 0)
//	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
//	if(IFMprob[0]-2*pi*n<0)
//		IFMprob+=2*pi//*(k+1)		// Correction necessary if Sag. phase is large, and thus the "0" point is in the wrong domain.
//	endif
	
	NVAR AvgSagGravPhi
	output = IFMprob - sagPhaseSign*AvgSagGravPhi
	
//	FindLevel/Q y 0
//	Variable zeroPhase = output[V_LevelX]
//	Variable piWraps = round(zeroPhase / (2*pi))
//	output -= piWraps*(2*pi)
	
	phaseOverPolMono*=pol

	unwrap 2*pi, output
	output -= 2*pi*round((output[p]-phaseOverPolMono[p])/(2*pi))
	
	//output-=2*pi*floor((output[p]-phaseOverPolMono[p])/(2*pi))
	//output-=2*pi*round((phaseOverPolMono[p]-output[p])/(2*pi))
	
	output += phaseoffset
	
	//SetDataFolder saveDFR
End



Function CalculateConstants2Electrodes()
//	Variable/G Vwire = (6035 * HVmon)	 +0.7							// Applied electrode voltage [V]
	
	NVAR HVp, HVm, mass
	
	Variable gapLoc = 1.905e-3//1.905e-3
	Variable rLoc = 25.4/2/2*1e-3
	
	//Variable/G Vwire = GetHVKeithley(HVmon, ZeroMon)*1000/2//; Vwire=9000
	Variable/G Vwire = (HVp+HVm)/2
	Variable/G d = gapLoc*sqrt(1+2*rLoc/gapLoc)//sqrt((R+gapLoc)^2 - R^2)								// Distance from ground plane to the "infinite wire" [m]
	Variable/G lambda = 2*pi*eps0*vwire*(ln((gapLoc+rLoc+d)/(gapLoc+rLoc-d)))^-1 //Vwire*pi*eps0 / asinh(d/R)				// Charge density of the "infinite wire"
	Variable/G E0 = -lambda / (pi*eps0)		
	Variable/G s = (LoneGint*hbar*2*pi)/(mass*GratingPeriod)								// beam separation at center of grad e region * velocity
	Variable/G AnalyticPrefactor=lambda^2*d/(hbar*pi*eps0^2)*4*pi*eps0*1E-30
	//Variable/G AnalyticPrefactor=4*pi*Vwire^2*d/hbar*((ln((gap+R+d)/(gap+R-d)))^-2)*4*pi*eps0*1E-30
	//Variable lambda2 = Vwire*2*pi*eps0/ln((gap+R+d)/(gap+R-d))
	//print lambda
	//print lambda2
	//print lambda2*2
	
	//print 43*AnalyticPrefactor/(d^2-.0017^2)/1800
	//print 43*AnalyticPrefactor/(d^2-.0017^2)/1800- 43*AnalyticPrefactor/(d^2-(.0017+50e-6)^2)/1800
End







//	Set the relavent parameters for polarizabilty fits for this data folder.
Function PolParams2()
	String CurrentDataFolder = GetDataFolder(0)
	
	//     General parameters
	Variable velLoc = NumVarOrDefault("velocity", 3000)
	Prompt velLoc, "Velocity (m/s): "
	Variable vratioLoc = NumVarOrDefault("vratio", 18)
	Prompt vratioLoc, "vratio (v/sig): "
	Variable molFracLoc = NumVarOrDefault("molFrac",0)
	Prompt molFracLoc, "Molecule Fraction (%): "
	Variable HVpLoc = NumVarOrDefault("HVp", 6000)
	Prompt HVpLoc, "HV+ (V): "
	Variable HVmLoc = NumVarOrDefault("HVm", 6000)
	Prompt HVmLoc, "HV- (V): "
//	Variable ZeroMonLoc = NumVarOrDefault("ZeroMon", .013)
//	Prompt ZeroMonLoc, "ZeroMon (V): "
//	Variable OffsetLoc = NumVarOrDefault("Offset", 1E-5)
//	Prompt OffsetLoc, "Offset (m): "
	String AtomLoc = StrVarOrDefault("Atom", "Potassium Avg Mass")
//	Prompt AtomLoc, "Atom Type: ", popup "Sodium;Potassium;Potassium 39;Rubidium"
	Prompt AtomLoc, "Atom Type: ", popup "Lithium Avg Mass;Lithium weighted;Sodium;Potassium Avg Mass;Potassium weighted;Potassium Avg test;Rubidium Avg Mass;Rubidium avg test;Rubidium weighted;Cesium"
	
	Variable FirstFileLoc = NumVarOrDefault("FirstFile", 1)
	Prompt FirstFileLoc, "First File: "
	Variable LastFileLoc = NumVarOrDefault("LastFile", 0)
	Prompt LastFileLoc, "Last File: "
	
	String FilesToKillStrLoc  = StrVarOrDefault("FilesToKillStr", "")
	Prompt FilesToKillStrLoc, "List files to kill separated by ;"
	
	Variable C0Loc = NumVarOrDefault("c0", 1)
	Prompt C0Loc, "Initial contrast: "
	
	Variable firstPosLGLoc = NumVarOrDefault("firstPosLG", 0)
	Prompt firstPosLGLoc, "First phi position (LG): "
	
	Variable lastPosLGLoc = NumVarOrDefault("lastPosLG", 0)
	Prompt lastPosLGLoc, "Last phi position (LG): "
	
	DoPrompt CurrentDataFolder+ " General Parameters", velLoc, vratioLoc, HVpLoc, HVmLoc, AtomLoc, MolFracLoc, FirstFileLoc, LastFileLoc, FilesToKillStrLoc, C0Loc
	
	DoPrompt CurrentDataFolder+ " Eclipse Parameters", firstPosLGLoc, lastPosLGloc
	
	Variable/G velocity = velLoc
	Variable/G vratio = vratioLoc
	Variable/G HVp = HVpLoc
	Variable/G HVm = HVmLoc
//	Variable/G ZeroMon = ZeroMonLoc
//	Variable/G offset = OffsetLoc
	String/G Atom = AtomLoc
	Variable/G molFrac = molFracLoc
	Variable/G c0 = c0Loc
	
	Variable/G firstPosLG = firstPosLGLoc
	Variable/G lastPosLG = lastPosLGLoc
	
	Variable/G mass
	Variable/G molPol			// molPol from Tarnovsky et al (1993)
	Variable/G initPol
	If(stringmatch(Atom, "Sodium"))
		mass = Na23mass*amu2kg	
		molPol = 40
		initPol = 24
	elseif(stringmatch(atom, "Lithium Avg Mass"))
		mass = Liavgmass*amu2kg
		molPol = 40
		initPol = 24
	elseif(stringmatch(atom, "Lithium weighted"))
		Variable/G numIsotopes = 2
		Make/D/O massWave = {Li6mass, Li7mass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {Li6frac, Li7frac}
		molPol = 40
		initPol = 24	
	elseif(stringmatch(atom, "Potassium Avg Mass"))
		mass = Kavgmass*amu2kg
		molPol = 77
		initPol = 43
	elseif(stringmatch(atom, "Potassium 39"))
		mass =  K39mass*amu2kg
		molPol = 77
		initPol = 43	
	elseif(stringmatch(atom, "Potassium weighted"))
		Variable/G numIsotopes = 3
		Make/D/O massWave = {K39mass, K40mass, K41mass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {K39frac, K40frac, K41frac}
		molPol = 77
		initPol = 43	
	elseif(stringmatch(atom, "Potassium avg test"))
		Variable/G numIsotopes = 1
		Make/D/O massWave = {Kavgmass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {1}
		molPol = 77
		initPol = 43
	elseif(stringmatch(atom, "Rubidium Avg Mass"))
		mass =   RbAvgmass*amu2kg	
		molPol = 79
		initPol = 47
	elseif(stringmatch(atom, "Rubidium avg test"))
		Variable/G numIsotopes = 1
		Make/D/O massWave = {Rbavgmass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {1}
		molPol = 79
		initPol = 47
	elseif(stringmatch(atom, "Rubidium weighted"))
		Variable/G numIsotopes = 2
		Make/D/O massWave = {Rb85mass, Rb87mass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {Rb85frac, Rb87frac}
		molPol = 79
		initPol = 47
	elseIf(stringmatch(Atom, "Cesium"))
		mass = CsAvgMass*amu2kg	
		molPol = 104
		initPol = 59
	EndIf
	
	Variable/G FirstFile = FirstFileLoc
	Variable/G LastFile = LastFileLoc

	Variable i = 0
	
	
	// parse FilesToKill string. end result is a numeric wave of file numbers to kill
	String/G FilesToKillStr = FilesToKillStrLoc
	String ranges = GrepList(FilesToKillStr, "-")//; print "ranges: ", ranges; 
	String oneRange
	Variable minRange, maxRange, pntsToInsert
	Make/O/N=0 rangesToKill						//Initialize rangesToKill wave
	for(i=0; i<ItemsInList(ranges); i+=1)
		oneRange = StringFromList(i, ranges)//; print i, " oneRange: ", oneRange
		minRange = str2num(StringFromList(0, oneRange, "-"))//; print minRange
		maxRange = str2num(StringFromList(1, oneRange, "-"))//; print maxRange
		pntsToInsert = maxRange-minRange+1
		InsertPoints 0, pntsToInsert, rangesToKill
		rangesToKill[0, pntsToInsert-1] = minRange+p
	endfor
	Variable NumFilesToKill = ItemsInList(FilesToKillStr)					//take a string of numbers seperated by ; and count how many there are
	if(NumFilesToKill != 0)
		Make/O/N=(NumFilesToKill) SingleFilesToKill
		SingleFilesToKill = Str2Num(StringFromList(p, FilesToKillStr))
		Concatenate/NP/O {rangesToKill, SingleFilestoKill}, FilesToKill
		Variable/G KillFiles = 1
	else
		Make/O/N=0 FilesToKill = NaN
		Variable/G KillFiles = 0
	EndIf

	
	//          HV phase unwrapping parameters
	Variable PlusMinusLoc = NumVarOrDefault("plusminus", 1)
	Prompt plusminusloc, "Add 2pi (1) or subtract 2pi (-1): "
	Variable UnwrapInitEndLoc = NumVarOrDefault("UnwrapInitEnd", 0)
	Prompt UnwrapInitEndLoc, "Last file for initial unwrapping: "
	Variable UnwrapInitEnd2Loc = NumVarOrDefault("UnwrapInitEnd2", 0)
	Prompt UnwrapInitEnd2Loc, "Last file for initial unwrapping 2: "
	Variable UnwrapInitEnd3Loc = NumVarOrDefault("UnwrapInitEnd3", 0)
	Prompt UnwrapInitEnd3Loc, "Last file for initial unwrapping 3: "
	Variable PlusMinusFinalLoc = NumVarOrDefault("plusminusFinal", 1)
	Prompt plusminusFinalloc, "Add 2pi (1) or subtract 2pi (-1) final: "
	Variable UnwrapFinalStartLoc = NumVarOrDefault("UnwrapFinalStart", 0)
	Prompt UnwrapFinalStartLoc, "First file for final unwrapping: "
	Variable UnwrapFinalStart2Loc = NumVarOrDefault("UnwrapFinalStart2", 0)
	Prompt UnwrapFinalStart2Loc, "First file for final unwrapping 2: "
	Variable UnwrapFinalStart3Loc = NumVarOrDefault("UnwrapFinalStart3", 0)
	Prompt UnwrapFinalStart3Loc, "First file for final unwrapping 3: "	
	
	DoPrompt CurrentDataFolder+" HV unwrapping", plusminusloc, UnwrapInitEndLoc, UnwrapInitEnd2Loc, UnwrapInitEnd3Loc, plusminusFinalloc, UnwrapFinalStartLoc, UnwrapFinalStart2Loc, UnwrapFinalStart3Loc//, Unwrap3EndLoc, 
	
	Variable/G plusminus = plusminusloc
	Variable/G plusminusFinal = plusminusFinalLoc
	
	if(UnwrapInitEndLoc != 0)
		Variable/G UnwrapInitEndYesNo = 1
		Variable/G UnwrapInitEnd = UnwrapInitEndLoc
	else
		Variable/G UnwrapInitEndYesNo = 0
		Variable/G UnwrapInitEnd = 0
	endif
	
	if(UnwrapInitEnd2Loc != 0)
		Variable/G UnwrapInitEnd2 = UnwrapInitEnd2Loc
	else
		Variable/G UnwrapInitEnd2 = 0
	endif
	
	if(UnwrapInitEnd3Loc != 0)
		Variable/G UnwrapInitEnd3 = UnwrapInitEnd3Loc
	else
		Variable/G UnwrapInitEnd3 = 0
	endif
	
	if(UnwrapFinalStartLoc != 0)
		Variable/G UnwrapFinalStartYesNo = 1
		Variable/G UnwrapFinalStart = UnwrapFinalStartLoc
	else
		Variable/G UnwrapFinalStartYesNo = 0
		Variable/G UnwrapFinalStart = 0
	endif
	
	if(UnwrapFinalStart2Loc != 0)
		Variable/G UnwrapFinalStart2 = UnwrapFinalStart2Loc
	else
		Variable/G UnwrapFinalStart2 = 0
	endif
	
	if(UnwrapFinalStart3Loc != 0)
		Variable/G UnwrapFinalStart3 = UnwrapFinalStart3Loc
	else
		Variable/G UnwrapFinalStart3 = 0
	endif

	
	
	//    Reference phase unwrapping
	Variable PlusMinusRefLoc = NumVarOrDefault("plusminusref", 1)
	Prompt plusminusrefloc, "Add 2pi (1) or subtract 2pi (-1): "
	Variable UnwrapRefStartLoc = NumVarOrDefault("UnwrapRefStart", 0)
	Prompt UnwrapRefStartLoc, "First file for unwrapping: "
//	String RefFilesToKillStrLoc  = StrVarOrDefault("RefFilesToKillStr", "")
//	Prompt RefFilesToKillStrLoc, "List Ref files to kill separated by ; "
	Variable InvertPhaseLoc = NumVarOrDefault("InvertPhase",0)
	Prompt InvertPhaseLoc, "Multiply phase shift by -1? yes (1) no (0)"
	
	DoPrompt CurrentDataFolder+" Reference phase unwrapping", plusminusrefloc, UnwrapRefStartLoc, InvertPhaseLoc
	
	Variable/G plusminusref = plusminusrefloc
	Variable/G InvertPhase = InvertPhaseLoc
	
	if(UnwrapRefStartLoc != 0)
		Variable/G UnwrapRefStartYesNo = 1
		Variable/G UnwrapRefStart = UnwrapRefStartLoc
	else
		Variable/G UnwrapRefStartYesNo = 0
		Variable/G UnwrapRefStart = 0
	endif
	
End







// Does the messy work of extracting useful phase and contrast data from the FringeFit results and plotting it.	
Function ExtractPhaseContrast2(SeriesName, oneKtwoK, makegraphsYesNo)
	String SeriesName		// ex: "e"
	Variable OnekTwok
	Variable makegraphsYesNo
	
	Variable autofilterRefContrast =3		// number of times to run the autofilter
	Variable sigmasRefContrast = 2			// sigmas for autofilter
	
	Variable autofilterRefPhase =3			// number of times to run the autofilter
	Variable sigmasRefPhase = 2			// sigmas for autofilter
	
	Variable refPolyTerms = 3				// number of terms to include in reference phase fit. 2 for a line, 3 for a quadratic, 4 for a cubic, etc.
	
	String onektwokStr = ""
	if(OnekTwok==1)
		onektwokStr = "1k2k"
	elseif(OnekTwok == 2)
		onektwokStr = "2k2k"
	endif
	
	String initposWaveName="pos_"+SeriesName
	String initphaseWaveName = "phase" + onektwokStr + "_"+SeriesName
	String initphase_errorWaveName = "phase_error" + onektwokStr + "_"+SeriesName
	String initcontrastWaveName = "contrast" + onektwokStr + "_"+SeriesName
	String initcontrast_errorWaveName = "contrast_error" + onektwokStr + "_"+SeriesName
	String inittimeWaveName = "time_"+SeriesName
	String initfileindexWaveName = "file_index_"+SeriesName
	String inithvpWaveName = "hvp_"+SeriesName		//hv postive voltage wave name
	String inithvmWaveName = "hvm_"+SeriesName
	
	Wave poswave = $initposWaveName
	Wave phase = $initphaseWaveName
	Wave phase_error = $initphase_errorWaveName
	Wave contrast = $initcontrastWaveName
	Wave contrast_error = $initcontrast_errorWaveName
	Wave timeWave = $inittimeWaveName
	Wave FileNumbers = $initfileindexWaveName
	Wave hvpWave = $inithvpWaveName
	Wave hvmWave = $inithvmWaveName
	
	Variable pntsInSeries = numpnts(phase)

	NVAR FirstFile, LastFile
	Variable minFile = FirstFile
	Variable maxFile = LastFile	!= 0 ? LastFile : numpnts(phase)			// "0" can be used for include all files
	
	Variable n=0; Variable i=0									
		
	String PosWaveName = NameOfWave(poswave) + "_polfit"
	String PhaseWaveName = NameOfWave(phase) + "_polfit"
	String PhaseErrorWaveName = NameOfWave(phase_error) + "_polfit"
	String contrastWaveName = NameOfWave(contrast) + "_polfit"
	String contrastErrorWaveName = NameOfWave(contrast_error) + "_fit"
	String hvpWaveName = NameOfWave(hvpWave) + "_polfit"
	String hvmWaveName = NameOfWave(hvmWave) + "_polfit"
	String filenumbersWaveName = NameOfWave(filenumbers) + "_polfit"
	String timeWaveName = NameOfWave(timewave) + "_polfit"


	Make/O/N=(maxFile-minFile+1) $PosWaveName
	Make/O/N=(maxFile-minFile+1) $PhaseWaveName
	Make/O/N=(maxFile-minFile+1) $PhaseErrorWaveName
	Make/O/N=(maxFile-minFile+1) $contrastWaveName
	Make/O/N=(maxFile-minFile+1) $contrastErrorWaveName
	Make/O/N=(maxFile-minFile+1) $timeWaveName
	Make/O/N=(maxFile-minFile+1) $fileNumbersWaveName
	Make/O/N=(maxFile-minFile+1) $hvpWaveName
	Make/O/N=(maxFile-minFile+1) $hvmWaveName
	
	Wave PolFitPos = $PosWaveName
	Wave PolFitPhase = $PhaseWaveName
	Wave PolFitPhaseError = $PhaseErrorWaveName
	Wave PolFitcontrast = $contrastWaveName
	Wave PolFitcontrastError = $contrastErrorWaveName
	Wave PolFitTimes = $timeWaveName
	Wave PolFitFileNumbers = $fileNumbersWaveName
	Wave PolFithvp = $hvpWaveName
	Wave PolFithvm = $hvmWaveName

	// assign 
	PolFitPos = poswave[p+minFile-1]
	PolFitPhase = phase[p+minFile-1]
	PolFitPhaseError = phase_error[p+minFile-1]
	PolFitcontrast = contrast[p+minFile-1]
	PolFitcontrastError = contrast_error[p+minFile-1]
	PolFitTimes = timeWave[p+minFile-1]
	PolFitFileNumbers = filenumbers[p+minFile-1]
	PolFitHvp = hvpWave[p+minFile-1]
	PolFitHvm = hvmWave[p+minFile-1]
	
	
	// Generate a manual tick wave for file numbers
	String timeTicksWaveName = NameOfWave(timewave) + "_polfit_ticks"
	String fileNumbersTicksWaveName = NameOfWave(filenumbers) + "_polfit_ticks"
	Extract/O PolFitFileNumbers, filenumbersTemp, mod(p,5)==0
	Extract/O PolFitTimes, timeticksTemp, mod(p,5)==0
	Duplicate/O timeTicksTemp $timeTicksWaveName
	Make/O/T/N=(numpnts(filenumbersTemp)) $fileNumbersTicksWaveName = num2str(filenumbersTemp)
	Wave PolFitFileNumbersTicks = $fileNumbersTicksWaveName
	Wave PolFitTimeTicks = $timeTicksWaveName
	
	
	String RefPhaseWaveName = NameOfWave(PolFitPhase)+"_ref"
	String RefPhaseErrorWaveName = NameOfWave(PolFitPhaseError)+"_ref"
	String RefPhaseResidWN = "Res_"+NameOfWave(PolFitPhase)+"_ref"
	String PhaseFitTextBoxName = "CF_"+NameOfWave(PolFitPhase)+"_ref"
	String RefcontrastWaveName = NameOfWave(PolFitcontrast)+"_ref"
	String RefcontrastErrorWaveName = NameOfWave(PolFitcontrastError)+"_ref"
	String RefcontrastResidWN = "Res_"+NameOfWave(PolFitcontrast)+"_ref"
	String ContrastFitTextBoxName = "CF_"+NameOfWave(PolFitcontrast)+"_ref"
	
	//Kill bad points
	NVAR KillFiles
	Wave FilesToKill
	Variable killindex
	If(KillFiles)
		for(i = 0; i<numpnts(FilesToKill); i+=1)
			Killindex = FilesToKill[i]-minFile
			PolFitPhase[Killindex] = NaN				//NaN the points to be killed, including the reference phase points
			PolFitPhaseError[Killindex] = NaN
			PolFitContrast[Killindex] = NaN
			PolFitContrastError[Killindex] = NaN
		endfor
	EndIf
	
	// Kill files hit by the auto filter
	//NVAR AutoKillFilesYesNo
	Wave/Z AutoFilesToKill
	If(numpnts(AutoFilesToKill)!=0)
		for(i = 0; i<numpnts(AutoFilesToKill); i+=1)
			Killindex = AutoFilesToKill[i]-minFile
			PolFitPhase[Killindex] = NaN				
			PolFitPhaseError[Killindex] = NaN
			PolFitContrast[Killindex] = NaN
			PolFitContrastError[Killindex] = NaN
		endfor
	EndIf
	
	
	Duplicate/O PolFitPhase $RefPhaseWaveName				//copy the phase data into a wave that will eventually only hold reference phase data
	Duplicate/O PolFitPhaseError $RefPhaseErrorWaveName
	Duplicate/O PolFitcontrast $RefcontrastWaveName				//copy the contrast data into a wave that will eventually only hold reference contrast data
	Duplicate/O PolFitcontrastError $RefcontrastErrorWaveName
	
	Wave RefPhase = $RefPhaseWaveName
	Wave RefPhaseError = $RefPhaseErrorWaveName
	Wave Refcontrast = $RefcontrastWaveName
	Wave RefcontrastError = $RefcontrastErrorWaveName


	RefPhase = PolFitHvp==0 && PolFitHvm==0 ? PolFitphase : NaN						//	Assign reference phase wave
	RefPhaseError = PolFitHvp==0 && PolFitHvm==0? PolFitphaseError : NaN			//	Assign reference phase error wave
	RefContrast = PolFitHvp==0 && PolFitHvm==0? PolFitcontrast : NaN					//	Assign reference contrast wave
	RefContrastError = PolFitHvp==0 && PolFitHvm==0 ? PolFitcontrastError : NaN	//	Assign reference contrast error wave
	
	PolFitPhase = PolFitHvp!=0 || PolFitHvm!=0 ? PolFitphase : NaN					//   Assign phase shift wave
	PolFitPhaseError = PolFitHvp!=0 || PolFitHvm!=0  ? PolFitphaseError : NaN		//   Assign phase shift wave
	PolFitContrast = PolFitHvp!=0 || PolFitHvm!=0  ? PolFitContrast : NaN				//   Assign contrast shift wave
	PolFitContrastError = PolFitHvp!=0 || PolFitHvm!=0  ? PolFitcontrastError : NaN	//   Assign contrast shift wave
	
	
	String DataFolderNameInit = GetDataFolder(0)
	String DataFolderName = ReplaceString("'", DataFolderNameInit, "")
	
	
	//kill first file from each set
	if(0)
		RefPhase = mod(p,5)==0 ? NaN : RefPhase						
		RefPhaseError = mod(p,5)==0 ? NaN : RefPhaseError			
		RefContrast = mod(p,5)==0 ? NaN : RefContrast					
		RefContrastError = mod(p,5)==0 ? NaN : RefContrastError
	
		PolFitPhase = mod(p,5)==0 ? NaN : PolFitPhase					
		PolFitPhaseError =mod(p,5)==0 ? NaN : PolFitPhaseError	
		PolFitContrast = mod(p,5)==0 ? NaN : PolFitContrast			
		PolFitContrastError =mod(p,5)==0 ? NaN : PolFitContrastError	
	endif


	Variable maxUnwrapIndex, NumPointsToUnwrap, UnwrapIndex
	Variable startUnwrapIndex

	NVAR PlusMinusRef, UnwrapRefStart, UnwrapRefStartYesNo
	
	if(UnwrapRefStartYesNo)
		//Unwrap
		i=0
		startUnwrapIndex = UnwrapRefStart-minFile
		do
			RefPhase[startUnwrapIndex+i] += plusminusRef*2*pi
			i+=1
		while(i<=numpnts(refphase))
	endif
	

	NVAR plusminus, plusminusFinal, UnwrapInitEndYesNo, UnwrapInitEnd, UnwrapFinalStartYesNo, UnwrapFinalStart, UnwrapInitEnd2, UnwrapInitEnd3, UnwrapFinalStart2, UnwrapFinalStart3


	if(UnwrapInitEndYesNo)
		//Unwrap
		i=0
		maxUnwrapIndex = UnwrapInitEnd-minFile
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<=maxUnwrapIndex)
	endif


	if(UnwrapInitEnd2)
		//Unwrap
		i=0
		maxUnwrapIndex = UnwrapInitEnd2-minFile
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<=maxUnwrapIndex)
	endif
	
	
	if(UnwrapInitEnd3)
		//Unwrap
		i=0
		maxUnwrapIndex = UnwrapInitEnd3-minFile
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<=maxUnwrapIndex)
	endif

	
	
	Variable startUnwrapIndex2
	if(UnwrapFinalStartYesNo)
		//Unwrap again
		i=0
		startUnwrapIndex2 = UnwrapFinalStart-minFile
		do
			PolFitPhase[startUnwrapIndex2+i] += plusminusFinal*2*pi
			i+=1
		while(i<numpnts(polfitphase)-startUnwrapIndex2)
	endif


	if(UnwrapFinalStart2)
		//Unwrap again
		i=0
		startUnwrapIndex2 = UnwrapFinalStart2-minFile
		do
			PolFitPhase[startUnwrapIndex2+i] += plusminusFinal*2*pi
			i+=1
		while(i<numpnts(polfitphase)-startUnwrapIndex2)
	endif
	
	
	if(UnwrapFinalStart3)
		//Unwrap again
		i=0
		startUnwrapIndex2 = UnwrapFinalStart3-minFile
		do
			PolFitPhase[startUnwrapIndex2+i] += plusminusFinal*2*pi
			i+=1
		while(i<numpnts(polfitphase)-startUnwrapIndex2)
	endif



	if(makegraphsYesNo)
	Display/K=1/W=(350,150,1250,650) PolFitcontrast vs FileNumbers
	ErrorBars $contrastWaveName Y,wave=(PolFitcontrastError,PolFitcontrastError)
	Label left "Contrast"; Label bottom "File number"
	AppendToGraph/L='Refcontrast' Refcontrast vs FileNumbers
	ModifyGraph mode=3,marker=8
	ErrorBars $RefcontrastWaveName Y,wave=(RefcontrastError,RefcontrastError)
	Label Refcontrast "Ref contrast"
	ModifyGraph axisEnab(left)={.6,1}
	ModifyGraph axisEnab(Refcontrast)={0,0.57}
	ModifyGraph freePos(Refcontrast)={0,bottom}
	ModifyGraph lblPosMode(refcontrast)=1
	TextBox/C/N=text0/A=MC DataFolderName
	TextBox/C/N=text0/X=9.58/Y=45.65
	ModifyGraph grid(bottom)=1,nticks(bottom)=20
	ModifyGraph gfSize=11
	WindowNamer(DataFolderName+" Contrast and Ref contrast " + onektwokStr)
	endif

	if(1)
		//Fit the reference contrast to a function and adjust the contrast
		Make/O/N=2 w_coef
		CurveFit/Q/NTHR=1/TBOX=768 line,  Refcontrast /W=RefcontrastError /X=PolFitFileNumbers /I=1 /D /R
		Wave RefContrastResid = $RefcontrastResidWN
		RefContrastResid = RefContrastResid == 0 ? NaN : RefContrastResid
		
		// autofilter on the reference contrast, and NaN the corresponding points in the reference phase
		for(i=1; i<=autofilterRefContrast; i+=1)		
			wavestats/q RefContrastResid
			Duplicate/O RefContrastResid killables
			killables = abs(RefContrastResid) > sigmasRefContrast*v_sdev
			RefContrast = killables == 1 ? NaN : RefContrast
			RefContrastResid = killables == 1 || numtype(RefContrast)==2? NaN : RefContrastResid
			RefPhase = killables == 1 ? NaN : RefPhase
			print "ref contrast files killed: ", sum(killables)
			CurveFit/Q/NTHR=1/TBOX=768 line,  Refcontrast /W=RefcontrastError /X=PolFitFileNumbers /I=1 /D /R
		endfor
		
		if(1)
			PolFitcontrast /= w_coef[0] + w_coef[1]*FileNumbers[p] 
			PolFitcontrastError /= w_coef[0] + w_coef[1]*FileNumbers[p] 
		else
			//Variable AvgRefContrastLoc = numVarOrDefault("AvgRefContrast",0.2)
			NVAR c0
			PolFitContrast /= c0
			PolFitContrastError /= c0
		endif
		
		if(makegraphsYesNo)
		TextBox/C/N=$contrastFitTextBoxName/X=10.7/Y=56
		ModifyGraph zero(Res_Refcontrast)=1
		endif
	endif
	


	if(0)
		Display/K=1/W=(350,150,1250,650) PolFitPhase vs PolFitFileNumbers
		AppendToGraph/L='RefPhase'/T RefPhase vs FileNumbers
		ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
		ErrorBars $RefPhaseWaveName Y,wave=(RefPhaseError,RefPhaseError)
		Label left "Phase"; Label bottom "File number";	Label RefPhase "RefPhase"
		ModifyGraph mode=3,marker=8
		ModifyGraph axisEnab(left)={0,.4}
		ModifyGraph axisEnab(RefPhase)={.43,1}
		ModifyGraph freePos(RefPhase)={0,bottom}
		ModifyGraph lblPosMode(refphase)=1
		ModifyGraph rgb($PhaseWaveName)=(0,0,65535)
		TextBox/C/N=text0/A=MC DataFolderName
		TextBox/C/N=text0/X=9.58/Y=45.65
		ModifyGraph gfSize=11
		WindowNamer(DataFolderName+" Phase and Ref Phase " + onektwokStr + " vs File")
	endif

	if(makegraphsYesNo)
	Display/K=1/W=(350,150,1250,650) PolFitPhase vs PolFitTimes
	AppendToGraph/L='RefPhase'/T RefPhase vs PolFitTimes
	ModifyGraph userticks(top)={PolFitTimeTicks,PolFitFileNumbersTicks}
	ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
	ErrorBars $RefPhaseWaveName Y,wave=(RefPhaseError,RefPhaseError)
	Label left "Phase"; Label bottom "Time";	Label RefPhase "RefPhase"
	ModifyGraph mode=3,marker=8
	ModifyGraph axisEnab(left)={0,.4}
	ModifyGraph axisEnab(RefPhase)={.43,1}
	ModifyGraph freePos(RefPhase)={0,bottom}
	ModifyGraph lblPosMode(refphase)=1
	ModifyGraph rgb($PhaseWaveName)=(0,0,65535)
	TextBox/C/N=text0/A=MC DataFolderName
	TextBox/C/N=text0/X=9.58/Y=45.65
	ModifyGraph gfSize=11
	ModifyGraph grid(top)=1
	WindowNamer(DataFolderName+" Phase and Ref Phase " + onektwokStr + " vs Time")
	endif
	
	Variable fitvsTime = 1		// 1 for fit vs time. 0 for fit vs file index
	
	if(1)
		//Fit the reference phase to a quadratic polynomial and adjust the phase
		Variable k
		Make/O w_coef
		if(fitvsTime)
			if(refPolyTerms==2)
				CurveFit/Q/TBOX=768 line  RefPhase /W=RefPhaseError /X=PolFitTimes /I=1 /D /R
			else
				CurveFit/Q/TBOX=768 poly refPolyTerms,  RefPhase /W=RefPhaseError /X=PolFitTimes /I=1 /D /R
			endif		
			Wave RefPhaseResid = $RefPhaseResidWN
			RefPhaseResid = RefPhaseResid == 0 ? NaN : RefPhaseResid

			// the autofilter
			for(i=1; i<=autofilterRefPhase; i+=1)		
				wavestats/q RefPhaseResid
				Duplicate/O RefPhaseResid killables
				killables = abs(RefPhaseResid) > sigmasRefPhase*v_sdev
				RefPhase = killables == 1 ? NaN : RefPhase
				RefPhaseResid = killables == 1 || numtype(refPhase)==2 ? NaN : RefPhaseResid
				print "ref phase files killed: ", sum(killables)
				if(refPolyTerms==2)
					CurveFit/Q/TBOX=768 line  RefPhase /W=RefPhaseError /X=PolFitTimes /I=1 /D /R
				else
					CurveFit/Q/TBOX=768 poly refPolyTerms,  RefPhase /W=RefPhaseError /X=PolFitTimes /I=1 /D /R
				endif	
			endfor

			for(k=0; k<refPolyTerms; k+=1)
				PolFitPhase -= w_coef[k]*TimeWave[p]^k
				//RefPhase -= w_coef[k]*TimeWave[p]^k
			endfor
			
		else
			CurveFit/Q/NTHR=1/TBOX=768 poly refPolyTerms,  RefPhase /W=RefPhaseError /X=PolFitFileNumbers /I=1 /D /R
			for(k=0; k<refPolyterms; k+=1)
				PolFitPhase -= w_coef[k]*FileNumbers[p]^k
				//RefPhase -= w_coef[k]*FileNumbers[p]^k
			endfor
		endif
	
		if(makegraphsYesNo)
		TextBox/C/N=$PhaseFitTextBoxName/X=10.7/Y=56
		ModifyGraph zero(Res_RefPhase)=1
		endif
	endif
	
	
	NVAR InvertPhase

	if(InvertPhase)							//Make the phase shift positive
		PolfitPhase*=-1
	endif

	// normalize contrast to c0
//	NVAR c0
//	PolFitContrast/= c0
//	PolFitContrastError/=c0
	
	// Make new wave names
	String PositionReducedWaveName = NameOfWave(PolFitpos) +"_red"
	String PhaseReducedWaveName= NameOfWave(PolFitphase) + "_red"
	String ErrorReducedWaveName= NameOfWave(PolFitPhaseerror) + "_red"
	String contrastReducedWaveName= NameOfWave(PolFitcontrast) + "_red"
	String ContrastErrorReducedWaveName= NameOfWave(PolFitcontrasterror) + "_red"
	String FileNumbersReducedWaveName = NameOfWave(FileNumbers)+"_red"
	
	// The reduced waves are initially just copies of the raw waves
	Duplicate/O polfitpos $PositionReducedWaveName; Wave PositionReduced = $PositionReducedWaveName
	Duplicate/O PolFitphase $PhaseReducedWaveName; Wave PhaseReduced = $PhaseReducedWaveName
	Duplicate/O PolFitphaseerror $ErrorReducedWaveName; Wave ErrorReduced = $ErrorReducedWaveName
	Duplicate/O PolFitcontrast $contrastReducedWaveName; Wave contrastReduced = $contrastReducedWaveName
	Duplicate/O PolFitContrasterror $ContrastErrorReducedWaveName; Wave ContrastErrorReduced = $ContrastErrorReducedWaveName
	Duplicate/O PolFitFileNumbers $FileNumbersReducedWaveName; Wave FileNumbersReduced = $FileNumbersReducedWaveName
	
	// If a phase, position, or error point is a NaN then we remove it and the corresponding points in the other waves
	RemoveNaNsXYZABC(PhaseReduced, PositionReduced, ErrorReduced, ContrastReduced, ContrastErrorReduced, FileNumbersReduced)
	
	// Sort the waves so that point 0 is closest to the ground plane (makes unwrapping easier)
	Sort PositionReduced,PositionReduced,PhaseReduced,ErrorReduced, ContrastReduced, ContrastErrorReduced, FileNumbersReduced
	
	if(makegraphsYesNo || 0)
	Display/K=1/W=(40,100,940,700) PhaseReduced vs PositionReduced
	AppendToGraph/L=contrast ContrastReduced vs PositionReduced
	ModifyGraph axisEnab(left)={.40,1.00}
	ModifyGraph axisEnab(contrast)={0,.35}
	ModifyGraph freePos(contrast)=0
	ModifyGraph lblPosMode(left)=1,lblPosMode(contrast)=1,lblMargin(left)=10;DelayUpdate
	ModifyGraph lblMargin(contrast)=10;DelayUpdate
	Label contrast "Contrast"
	ModifyGraph mode=3,marker=8
	ErrorBars $PhaseReducedWaveName Y,wave=(ErrorReduced,ErrorReduced)
	ErrorBars $ContrastReducedWaveName Y,wave=(ContrastErrorReduced,ContrastErrorReduced)
	Label left "Phase"; Label bottom "Pillar position (um)"
	ModifyGraph zero(left)=1
	//SetAxis left 0,*
	//SetAxis bottom 0,*
	//SetAxis contrast 0,*
	ModifyGraph nticks(contrast)=10;DelayUpdate
	ModifyGraph gmSize=4
	if(1)	// set to true for finer ticks and measure distance in mm
		ModifyGraph nticks(bottom)=10,minor(bottom)=1
		Label bottom "Pillar position (um)"
		ModifyGraph minor(left)=1
		ModifyGraph sep(left)=4
	endif
	ModifyGraph gfSize=11
	TextBox/C/N=text0/A=MC DataFolderName
	TextBox/C/N=text0/X=-16.00/Y=21.00
	WindowNamer(DataFolderName+" Pillars Data " + onektwokStr)
	endif
	
	if(0)
		Print PositionReducedWaveName
		Print PhaseReducedWaveName
		Print ErrorReducedWaveName
		Print ContrastReducedWaveName
		Print ContrastErrorReducedWaveName
		Print FileNumbersReducedWaveName
	endif
	
	if(numVarOrDefault("startPosAbs",0) && 1)
		PositionReduced += numVarOrDefault("startPosAbs",0)
		ModifyGraph zero(bottom)=1
	endif
End



// run fitting routine on all data series in data folders formatted as e.g. '130726a', '130726b', etc...
Function AutoFitAllDF()
	setdatafolder root:
	
	String holdStr = "11"
	
	DFREF rootDF = root:
	
	Variable numDFs = CountObjectsDFR(rootDF,4)
	//numDFs-=4
	
	String currentDFStr
	String dataSetStr
	String windowName
	String windowList = WinList("*",";","WIN:1")
	String series_name 
	Variable thedate
	
	make/D/o/n=7 bestFitParamsLoc
	make/D/o/n=(numDFs) bestFitPolWave, bestFitPolSigWave, bestFitx0Wave, bestFitx0SigWave, bestFitphi0Wave, bestFitphi0SigWave, chisqrdDOFWave
	make/O/T/n=(numDFs) dataSeriesNames
	
	Variable polcutsigma = 2.5
	Variable polcutIterations = 3
	
	DFREF currentDF
	Variable i, j
	For(i=0; i<(numDFs); i+=1)
		dataSetStr = GetIndexedObjNameDFR(rootDF,4,i)
		currentDFStr = "root:'"+dataSetStr+"'"
		SetDataFolder $currentDFStr
		
		windowName = StringFromList(0, GrepList(windowList,"WN"+datasetStr,0,";"))
		DoWindow/F $windowName
		
		sscanf dataSetStr, "%d%s", thedate, series_name
		dataSeriesNames[i]= dataSetStr
		
		printf "\r%.0f%s\r\r", thedate, series_name
		
		Make/O/N=0 AutoFilesToKill, AutoFilesToKillOld
		
		ExtractPhaseContrast2(series_name, 1,0)
		Wave bestFitOutputWave = PhaseFit2Electrodes(series_name, 1, holdStr)
		
		for(j=1; j<=polcutIterations; j+=1)
			FindKillablesPolv2(series_name, polcutsigma)
			ExtractPhaseContrast2(series_name, 1,0)
			Wave bestFitOutputWave = PhaseFit2Electrodes(series_name, 1, holdStr)
		endfor
		
		bestFitPolWave[i]=bestFitOutputWave[0]
		bestFitPolSigWave[i]=bestFitOutputWave[1]
		bestFitx0Wave[i]=bestFitOutputWave[2]
		bestFitx0SigWave[i]=bestFitOutputWave[3]
		bestFitphi0Wave[i]=bestFitOutputWave[4]
		bestFitphi0SigWave[i]=bestFitOutputWave[5]
		chisqrdDOFWave[i]=bestFitOutputWave[6]
	EndFor
	
	setdatafolder root:
	
	wavestats/q bestfitpolwave
	
	printf "\r\r AutoFitAllDF() results:\r alpha = %.3f pm %.3f\r", v_avg, v_sem
	
	//edit dataseriesNames, bestFitPolWave, bestFitPolSigWave, bestFitx0Wave, bestFitx0SigWave, bestFitphi0Wave, bestFitphi0SigWave, chisqrdDOFWave
End




Function BrowseDF()
	setdatafolder root:
	
	DFREF rootDF = root:
	
	Variable numDFs = CountObjectsDFR(rootDF,4)
	
	String currentDFStr
	String dataSetStr
	String windowName
	String windowList = WinList("*",";","WIN:1")
	String series_name 
	Variable thedate
	

	
	DFREF currentDF
	Variable i
	For(i=0; i<numDFs; i+=1)
		dataSetStr = GetIndexedObjNameDFR(rootDF,4,i)
		currentDFStr = "root:'"+dataSetStr+"'"
		SetDataFolder $currentDFStr
		
		windowName = StringFromList(0, GrepList(windowList,"WN"+datasetStr,0,";"))
		DoWindow/F $windowName
		
		sscanf dataSetStr, "%d%s", thedate, series_name
		printf "\r%.0f%s\r\r", thedate, series_name
		
		Wave FilesToKill
		print filestokill
		
		//Wave autofilestokill
		//print autofilestokill
		
		polparams2()
	EndFor
	
	setdatafolder root:
End



// Useful for backing up or saving the results of many data sets before recalculating using different parameters
Function BatchDuplicate(appendedName)
	String appendedName
	
	String dataSeriesNamesStr = "dataSeriesNames" + appendedName
	String bestFitPolWaveStr = "bestFitPolWave" + appendedName
	String bestFitPolSigWaveStr = "bestFitPolSigWave" + appendedName
	String bestFitx0WaveStr = "bestFitx0Wave" + appendedName
	String bestFitx0SigWaveStr = "bestFitx0SigWave" + appendedName
	String bestFitphi0WaveStr = "bestFitphi0Wave" + appendedName
	String bestFitphi0SigWaveStr = "bestFitphi0SigWave" + appendedName
	String chisqrdDOFWaveStr = "chisqrdDOFWave" + appendedName
	
	Duplicate/O dataSeriesNames $dataSeriesNamesStr; 
	Duplicate/O bestFitPolWave $bestFitPolWaveStr; 
	Duplicate/O bestFitPolSigWave $bestFitPolSigWaveStr; 
	Duplicate/O bestFitx0Wave $bestFitx0WaveStr; 
	Duplicate/O bestFitx0SigWave $bestFitx0SigWaveStr; 
	Duplicate/O bestFitphi0Wave $bestFitphi0WaveStr; 
	Duplicate/O bestFitphi0SigWave $bestFitphi0SigWaveStr; 
	Duplicate/O chisqrdDOFWave $chisqrdDOFWaveStr; 
	
	BatchGraph(appendedName)
end



//  Takes the average of a data series and makes a nice plot of results.
Function BatchGraph(appendedName)
	String appendedName
	
	// get references to the appropriate waves
	String dataSeriesNamesStr = "dataSeriesNames" + appendedName
	String bestFitPolWaveStr = "bestFitPolWave" + appendedName
	String bestFitPolSigWaveStr = "bestFitPolSigWave" + appendedName
	String bestFitx0WaveStr = "bestFitx0Wave" + appendedName
	String bestFitx0SigWaveStr = "bestFitx0SigWave" + appendedName
	String bestFitphi0WaveStr = "bestFitphi0Wave" + appendedName
	String bestFitphi0SigWaveStr = "bestFitphi0SigWave" + appendedName
	String chisqrdDOFWaveStr = "chisqrdDOFWave" + appendedName
	
	Wave/T dataSeriesNamesNew = $dataSeriesNamesStr 
	Wave bestFitPolWaveNew = $bestFitPolWaveStr 
	Wave bestFitPolSigWaveNew = $bestFitPolSigWaveStr 
	Wave bestFitx0WaveNew = $bestFitx0WaveStr 
	Wave bestFitx0SigWaveNew = $bestFitx0SigWaveStr 
	Wave bestFitphi0WaveNew = $bestFitphi0WaveStr 
	Wave bestFitphi0SigWaveNew = $bestFitphi0SigWaveStr 
	Wave chisqrdDOFWaveNew = $chisqrdDOFWaveStr 
	
	// add an entry for the average in the data series names
	Redimension/N=(numpnts(bestFitPolWaveNew)+1) dataSeriesNamesNew
	dataSeriesNamesNew[numpnts(bestFitPolWaveNew)] = "Avg"
	
	String dataSeriesPlotNumStr = "dataSeriesPlotNum"+appendedName
	make/o/n=(numpnts(bestFitPolWaveNew)+1) $dataSeriesPlotNumStr = p
	Wave dataSeriesPlotNum = $dataSeriesPlotNumStr
	
	String avgPlotNumStr = "avgPlotNum"+appendedName
	make/o/n=1 $avgPlotNumStr = numpnts(bestFitPolWaveNew)
	Wave avgPlotNum = $avgPlotNumStr

	//Prepare waves to hold the average, std dev, and std err. use a 2nd average wave to display both the std dev and std err	
	String dataSeriesAvgStr = "dataSeriesAvg"+appendedName
	String dataSeriesAvg2Str = "dataSeriesAvg2"+appendedName
	String dataSeriesAvgSigStr = "dataSeriesAvgSig"+appendedName
	String dataSeriesAvgSemStr = "dataSeriesAvgSem"+appendedName
	Make/o/n=1 $dataSeriesAvgStr, $dataSeriesAvg2Str, $dataSeriesAvgSigStr, $dataSeriesAvgSemStr
	Wave dataSeriesAvg = $dataSeriesAvgStr
	Wave dataSeriesAvg2 = $dataSeriesAvg2Str
	Wave dataSeriesAvgSig = $dataSeriesAvgSigStr
	Wave dataSeriesAvgSem = $dataSeriesAvgSemStr
	
	Variable i
	
	//Prepare a wave to hold only the data series we want to average
	String bestFitPolWaveToAvgStr = "bestFitPolWaveToAvg"+appendedName
	Duplicate/O bestFitPolWaveNew $bestFitPolWaveToAvgStr
	Wave bestFitPolWaveToAvg = $bestFitPolWaveToAvgStr
	
	//Perpare waves to overlay X's on top of ignored data points
	String dataSeriesNixedPolStr = "dataSeriesNixedPol" + appendedName
	Make/o/n=0 $dataSeriesNixedPolStr
	Wave dataSeriesNixedPol = $dataSeriesNixedPolStr
	String dataSeriesNixedPntsStr = "dataSeriesNixedPnts" + appendedName
	Make/o/n=0 $dataSeriesNixedPntsStr
	Wave dataSeriesNixedPnts = $dataSeriesNixedPntsStr
	
	//Find the kill series wave. If it exists, match up the specified series with the corresponding data point.
	String dataSeriesList 
	Variable pntNum
	String killSeriesWaveStr = "killseries"+appendedName
	Wave/T/Z killSeries = $killSeriesWaveStr

	if(waveexists(killSeries))
		Make/o/n=(numpnts(killSeries)) killSeriesPnts
		For(i=0; i<numpnts(killSeries); i+=1)
			dataSeriesList = TextWaveToStringList(dataSeriesNamesNew)
			pntNum = WhichListItem(killSeries[i], dataSeriesList)
			killSeriesPnts[i]= pntNum
			bestFitPolWaveToAvg[pntNum] = NaN
			Redimension/N=(numpnts(dataSeriesNixedPol)+1) dataSeriesNixedPol, dataSeriesNixedPnts
			dataSeriesNixedPol[i] = bestFitPolWaveNew[pntNum]
			dataSeriesNixedPnts[i] = pntNum
		EndFor
	else

	endif
	
	// calculate and print the average, std dev, and std err (igorning data series that we don't want)	
	wavestats/Q bestFitPolWaveToAvg
	dataSeriesAvg = V_avg
	dataSeriesAvg2 = V_avg
	dataSeriesAvgSig = V_sdev
	dataSeriesAvgSem = V_sem
	
	printf "alpha = %2.3f\r" v_avg
	printf "std dev = %.3f\r" v_sdev
	printf "std err = %.3f\r"  v_sem
	
	//Make a wave for an average line
	String dataSeriesAvgLineStr = "dataSeriesAvgLine"+appendedName
	String dataSeriesAvgLinePntsStr = "dataSeriesAvgLinePnts"+appendedName
	Make/O/N=2 $dataSeriesAvgLineStr, $dataSeriesAvgLinePntsStr
	Wave dataSeriesAvgLine = $dataSeriesAvgLineStr
	Wave dataSeriesAvgLinePnts = $dataSeriesAvgLinePntsStr
	dataSeriesAvgLine = v_avg
	dataSeriesAvgLinePnts[0] = -1
	dataSeriesAvgLinePnts[1] = dataSeriesPlotNum[numpnts(dataSeriesPlotNum)-1]
	
	//make a nice graph
	display/k=1/W=(30,50,700,400) bestFitPolWaveNew vs dataSeriesPlotNum
	appendToGraph dataSeriesAvg vs avgPlotNum
	appendToGraph dataSeriesAvg2 vs avgPlotNum
	AppendToGraph dataSeriesNixedPol vs dataSeriesNixedPnts
	AppendToGraph dataSeriesAvgLine vs dataSeriesAvgLinePnts
	ModifyGraph mode=3,marker=8,rgb=(1,4,52428)
	ModifyGraph mode($dataSeriesAvgLineStr)=0
	ModifyGraph marker($dataSeriesAvgStr)=19
	ErrorBars $bestFitPolWaveStr Y,wave=($bestFitPolSigWaveStr, $bestFitPolSigWaveStr)
	ErrorBars $dataSeriesAvgStr Y,wave=($dataSeriesAvgSigStr, $dataSeriesAvgSigStr)
	ErrorBars $dataSeriesAvg2Str Y,wave=($dataSeriesAvgSemStr, $dataSeriesAvgSemStr)
	ModifyGraph marker($dataSeriesNixedPolStr)=1
	ModifyGraph msize($dataSeriesNixedPolStr)=8
	ModifyGraph userticks(bottom)={dataSeriesPlotNum,dataSeriesNamesNew}
	ModifyGraph tkLblRot(bottom)=90
	ModifyGraph standoff(left)=0
	SetAxis bottom -0.5,*
	ModifyGraph minor(left)=1
	ShowInfo
	
	String tboxStr
	sprintf tboxStr,  "\\Zr150\\F'Symbol'a = %2.3f(%2.0f)\rs = %0.3f\rN = %.0f", v_avg, v_sem*1000, v_sdev, v_npnts
	TextBox/C/N=alpha tboxStr
	
	WindowNamer("Pol Analysis Summary " + appendedName)
End







Function/S TextWaveToStringList(aTextWave)
	Wave/T aTextWave
	
	String listOut =""
	
	Variable i
	For(i=0; i<numpnts(aTextWave);  i+=1)
		listOut += aTextWave[i]+";"
	EndFor
	
	return listOut
End





Function AddSeriesToList()
	Variable pnt = pcsr(a)
	
	Wave seriesNames
	
	String listOfNames = TextWaveToStringList(seriesNames)
End



// Calculate Sagnac and gravitation phase shifts as a function of velocity
Function SagnacAndGravityPhase()
	Wave vel, dens, vrecip
	
	Variable verbose=0	// set = 1 if you want phase shift and contrast printed to history
	
	Variable GratingTilt = 0e-3	// in radians         // 0.5*pi/180
	Variable L1g2g = 0.94		// meters
	Variable g = TucsonGravity				
	Variable gFactor = g*sin(GratingTilt)*L1g2g^2/GratingPeriod*2*pi
	MatrixOP/O phiGrav = gFactor * powR(vel,-2)

	Variable Latitude = 32.2	//degrees
	Variable OmegaEarth = 2*pi/(24*3600)*(366.25/365.25)*cos((90-Latitude)*pi/180)    //; print OmegaEarth	// earth rotation rate for a given latitude
	//omegaEarth =0
	Variable SagFactor = OmegaEarth*L1g2g^2*4*pi/1e-7
	//print sagfactor/3000
	MatrixOP/O phiSag = SagFactor *rec( vel)
	
	MatrixOP/O phiSagGrav = phiGrav+phiSag
	MatrixOP/O phiSagGravWeighted = dens*phiSagGrav
	Variable phiSagGravAvg = area(phiSagGravWeighted); 
	
	// Calculate how much we'll need to unwrap the phase when we take the atan2 of the integrated e^iphi
	Variable/G k = -1
	do
		k+=1
	while(phiSagGravAvg-k*2*pi > pi)
	//print k
		
	MatrixOP/O phiSagGravImag = dens * sin(phiSagGrav)
	MatrixOP/O phiSagGravReal = dens * cos(phiSagGrav)
	Variable phiSagGravImagInt = area(phiSagGravImag)
	Variable phiSagGravRealInt = area(phiSagGravReal)
	
	Variable/G AvgSagGravPhi = atan2(phiSagGravImagInt, phiSagGravRealInt) + 2*pi*k
	Variable/G AvgSagGravCon = sqrt(phiSagGravImagInt^2+phiSagGravRealInt^2)
	
	if(verbose)
		printf "Improper velocity averaged Sagnac and gravity phase shift: %2.7f\r", phiSagGravAvg
		printf "Velocity averaged Sagnac and gravity phase shift: %2.7f\r", avgSagGravphi
		printf "Velocity averaged Sagnac and gravity contrast: %1.3f\r", avgSagGravCon
	endif
	
	Duplicate/O density2D phiSagGravLayered2D
	phiSagGravLayered2D = phiSagGrav[y]
End



// Calculate the position and velocity dependent phase shift / pol.
// This only needs to be done once and then the fit routine can multiply the results by polarizability on each iteration.
Function PolPhaseAn()
	Wave vrecip,yy
	NVAR AnalyticPrefactor, d,s
	Variable dsqr=d^2
	
	MatrixOP/O IFMphase_v_overpolM = vrecip*(rec(powr(yy,2)-dsqr)-rec(powr(yy-s*vrecip,2)-dsqr))*(-1)
	MatrixOP/O IFMphase_v_overpolP = vrecip*(rec(powr(yy+s*vrecip,2)-dsqr)-rec(powr(yy,2)-dsqr))*(-1)
	MatrixOP/O IFMphase_v_overpolMmol = vrecip*(rec(powr(yy,2)-dsqr)-rec(powr(yy-s*vrecip/2,2)-dsqr))*(-1)
	MatrixOP/O IFMphase_v_overpolPmol = vrecip*(rec(powr(yy+s*vrecip/2,2)-dsqr)-rec(powr(yy,2)-dsqr))*(-1)

	MatrixOP/O IFMphase_v_overpolM2 = vrecip*(rec(powr(yy-s*vrecip,2)-dsqr)-rec(powr(yy-2*s*vrecip,2)-dsqr))*(-1)		// second order IFM. phi_(-1)-phi_(-2)
	MatrixOP/O IFMphase_v_overpolP2 = vrecip*(rec(powr(yy+2*s*vrecip,2)-dsqr)-rec(powr(yy+s*vrecip,2)-dsqr))*(-1)
End



// Prepare the desired velocity distribution.
Function NormalizeVelocitiesAn(velocity, sigma, n)			//Find the normalization for the density profile. 
	Variable velocity, sigma, n

	variable flowv = velocity
	Variable nSigmas = 7
	Variable minVel = FlowV*(1-nSigmas/15)//; minVel = 0
	Variable maxVel = FlowV*(1+nSigmas/15)//; maxVel = 5000
	Variable velStep = (maxVel-minVel)/number_v


	NVAR num_pos
	Make/D/o/n=(number_v) vel, dens
	SetScale/I x minVel,maxVel,"m/s", vel//, dens
	vel = x //velocity + (x * v_step - 0.5 * number_v * v_step)			// make a wave for velocities centered on "velocity"
	dens = vel^n* exp (-(vel-velocity)^2/(2*sigma^2)) 	// prob(v)
	Variable vel_normalization = area(dens)							// used to normalize Int(prob(v)) to 1. builds in xscaling 
	dens /= vel_normalization
	
	Variable maxprob = WaveMax(dens)
	if(dens[0]/maxprob > .001 || dens[numpnts(dens)] > .001)
		print "VELOCITY DISTRIBUTION IS NOT WIDE ENOUGH"
		//DoAlert 0, "Velocity distribution is not wide enough. Adjust number_v or v_step and rerun."
	endif
	
	if(0)
		Display/K=1 dens vs vel
		Label left "Probability"
		Label bottom "Velocity (m/s)"
	endif
	
	Make/D/o/n=(num_pos, number_v) v, density, vrecip, density2D
	SetScale/I y minVel,maxVel,"m/s", v//, density, vrecip, density2D
	v =y//( (velocity) - (.5 * number_v * v_step ))+(v_step * q) // all columns have the same velocity
	
	MatrixOP/O/S density = powR(v,n)* exp (-magsqr(v-(velocity))/(2*magsqr(sigma))) / vel_normalization	//Probability density distribution. same layering as v
	MatrixOP/O/S vrecip = rec(v)	

	density2D = dens[q]
end




// Clean up all the big unneeded waves after the calculation is done.
Function CleanupPolWaves()
	KillWaves/Z dens, density, DminusYsep, DminusYY, dphi, dphi_s, DplusYsep, DplusYY, Etot, EtotMagSqr, Etot_s, Etot_sMagSqr, v
	KillWaves/Z imagpart, realpart, xx, xxsqr, yy
	KillWaves/Z dphi_sM, dphi_sP, EtotMagSqr0, EtotMagSqrM, EtotMagSqrP, IFMphase, IFMphase_v, IFMphase_vM, IFMPhase_vP, IFMphase_v_overpol, IFMphase_v_overpolM, IFMphase_v_overpolP, vel, vrecip
	KillWaves/Z EtotMagSqrMmol, EtotMagSqrPmol, dphiMmol, dphiPmol, IFMphase_v_overpolMmol, IFMphase_v_overpolPmol, IFMphase_vMmol, IFMphase_vPmol
	KillWaves/z phiSag, phiSagAvg, phiSagLayered
	KillWaves/Z density2D, phiGrav, phiSagGrav, phiSagGravImag, phiSagGravLayered2D, phiSagGravReal, phiSagGravWeighted, phiSagImag, phiSagLayered2D, phiSagReal
	KillWaves/Z EtotMagSqrM0, EtotMagSqrMmol0, EtotMagSqrP0, EtotMagSqrPmol0, EtotMagSqrM1, EtotMagSqrMmol1, EtotMagSqrP1, EtotMagSqrPmol1,EtotMagSqrM2, EtotMagSqrMmol2, EtotMagSqrP2, EtotMagSqrPmol2
	KillWaves/Z  IFMphase_v_overpolM0, IFMphase_v_overpolP0,IFMphase_v_overpolMmol0, IFMphase_v_overpolPmol0, IFMphase_v_overpolM1, IFMphase_v_overpolP1,IFMphase_v_overpolMmol1, IFMphase_v_overpolPmol1, IFMphase_v_overpolM2, IFMphase_v_overpolP2,IFMphase_v_overpolMmol2, IFMphase_v_overpolPmol2
	KillWaves/Z EsqrdSepWaveRefs, IFMphase_v_overpol_Refs
	KillWaves/Z dphi_sM2, dphi_sP2, dphiM2mol, dphiP2mol, IFMphase_v_overpolM2, IFMphase_v_overpolP2, IFMphase_v_overpolM2mol, IFMphase_v_overpolP2mol
	KillWaves/Z EtotMagSqrM2, EtotMagSqrP2, EtotMagSqrM2mol, EtotMagSqrP2mol
	KillWaves/Z IFMphase_vM2, IFMphase_vP2, IFMphase_v_overpolM2, IFMphase_v_overpolP2
End