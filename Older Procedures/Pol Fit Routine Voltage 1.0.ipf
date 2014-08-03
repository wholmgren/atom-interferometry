#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.0		// 1.0 = 1 IFM model, 0 position = 50% counts, high polarizabilities

#include ":Pol Fit Routine"


Function PolParamsVoltage()
	String CurrentDataFolder = GetDataFolder(0)
	
	//     General parameters
	Variable velLoc = NumVarOrDefault("velocity", 2000)
	Prompt velLoc, "Velocity: "
	Variable sigLoc = NumVarOrDefault("sigma", 100)
	Prompt sigLoc, "Sigma: "
	Variable InitialGradEMotorPositionLoc = NumVarOrDefault("InitialGradEMotorPosition", .002)
	Prompt InitialGradEMotorPositionLoc, "InitialGradEMotorPosition (m): "
	Variable FinalGradEMotorPositionLoc = NumVarOrDefault("FinalGradEMotorPosition", .0017)
	Prompt FinalGradEMotorPositionLoc, "FinalGradEMotorPosition (m): "
	Variable PositionLoc = NumVarOrDefault("Position", .00175)
	Prompt PositionLoc, "Position (m): "
	String AtomLoc = StrVarOrDefault("Atom", "Sodium")
//	Prompt AtomLoc, "Atom Type: ", popup "Sodium;Potassium;Potassium 39;Rubidium"
	Prompt AtomLoc, "Atom Type: ", popup "Sodium;Potassium Avg Mass;Rubidium Avg Mass"
	
	Variable LastEclipseLoc = NumVarOrDefault("LastEclipse", 5)
	Prompt LastEclipseLoc, "Last Eclipse File: "
	Variable LastFileLoc = NumVarOrDefault("LastFile", 80)
	Prompt LastFileLoc, "Last File: "
//	Variable KillEclipseFilesLoc = NumVarOrDefault("KillEclipseFiles", 0)
//	Prompt KillEclipseFilesLoc, "Kill Eclipse Files: yes (1) no (0) "
	String EclipseFilesToKillStrLoc  = StrVarOrDefault("EclipseFilesToKillStr", "")
	Prompt EclipseFilesToKillStrLoc, "List eclipse files to kill sep. by ; "
	
	
	DoPrompt CurrentDataFolder+ " General Parameters", velLoc, sigLoc, InitialGradEMotorPositionLoc, FinalGradEMotorPositionLoc, PositionLoc, AtomLoc, LastEclipseLoc, LastFileLoc, EclipseFilesToKillStrLoc
	
	
	Variable/G velocity = velLoc
	Variable/G sigma = sigLoc
	Variable/G InitialGradEMotorPosition = InitialGradEMotorPositionLoc
	Variable/G FinalGradEMotorPosition = FinalGradEMotorPositionLoc
	Variable/G Position = PositionLoc
	String/G Atom = AtomLoc
	
	Variable/G mass
	If(stringmatch(Atom, "Sodium"))
		mass = 3.81754e-26	
	elseif(stringmatch(atom, "Potassium Avg Mass"))
		mass = 6.49242e-26	
	elseif(stringmatch(atom, "Potassium 39"))
		mass =  6.4761e-26	
	elseif(stringmatch(atom, "Rubidium Avg Mass"))
		mass =   1.41923e-25	
	EndIf
	
	Variable/G LastEclipse = LastEclipseLoc
	Variable/G LastFile = LastFileLoc
	//Variable/G KillEclipseFiles = KillEclipseFilesLoc
	
	String/G EclipseFilesToKillStr = EclipseFilesToKillStrLoc
	Variable i = 0
	Variable NumEclipseFilesToKill = ItemsInList(EclipseFilesToKillStr)
	if(NumEclipseFilesToKill != 0)
		Make/O/N=(NumEclipseFilesToKill) EclipseFilesToKill
		do
			EclipseFilesToKill[i] = Str2Num(StringFromList(i, EclipseFilesToKillStr))
			i+=1
		while(i < NumEclipseFilesToKill)
		Variable/G KillEclipseFiles = 1
	else
		If(WaveExists(EclipseFilesToKill))
			EclipseFilesToKill = NaN
		EndIf
		Variable/G KillEclipseFiles = 0
	EndIf
	
	
	
	//          HV phase unwrapping
	Variable PlusMinusLoc = NumVarOrDefault("plusminus", 1)
	Prompt plusminusloc, "Add 2pi (1) or subtract 2pi (-1): "
	Variable Unwrap1Loc = NumVarOrDefault("Unwrap1", 0)
//	Prompt Unwrap1Loc, "Turn on first unwrap? 1 or 0"
	Variable Unwrap1EndLoc = NumVarOrDefault("Unwrap1End", 0)
	Prompt Unwrap1EndLoc, "First file for first unwrapping: "
	Variable Unwrap2Loc = NumVarOrDefault("Unwrap2", 0)
//	Prompt Unwrap2Loc, "Turn on second unwrap? 1 or 0"
	Variable Unwrap2EndLoc = NumVarOrDefault("Unwrap2End", 0)
	Prompt Unwrap2EndLoc, "First file for second unwrapping: "
	Variable Unwrap3Loc = NumVarOrDefault("Unwrap3", 0)
//	Prompt Unwrap3Loc, "Turn on second unwrap? 1 or 0"
	Variable Unwrap3EndLoc = NumVarOrDefault("Unwrap3End", 0)
	Prompt Unwrap3EndLoc, "First file for third unwrapping: "
	String HVFilesToKillStrLoc  = StrVarOrDefault("HVFilesToKillStr", "")
	Prompt HVFilesToKillStrLoc, "List HV files to kill separated by ; "
	Variable AddSubAllLoc = NumVarOrDefault("AddSubAll", 0)
	Prompt AddSubAllLoc, "Add(1), Subtract(-1), Nothing(0) 2pi from all phase shift data: "
		
//	DoPrompt CurrentDataFolder, plusminusloc, Unwrap1Loc, Unwrap1EndLoc, Unwrap2loc, Unwrap2EndLoc, Unwrap3Loc, Unwrap3EndLoc
	DoPrompt CurrentDataFolder+" HV unwrapping", plusminusloc, Unwrap1EndLoc, Unwrap2EndLoc, Unwrap3EndLoc, HVFilesToKillStrLoc, AddSubAllLoc
	
	Variable/G plusminus = plusminusloc
	Variable/G AddSubAll = AddSubAllLoc
	
	if(Unwrap1EndLoc != 0)
		Variable/G Unwrap1 = 1
		Variable/G Unwrap1End = Unwrap1EndLoc
	else
		Variable/G Unwrap1 = 0
		Variable/G Unwrap1End = 0
	endif
	
	if(Unwrap2EndLoc != 0)
		Variable/G Unwrap2 = 1
		Variable/G Unwrap2End = Unwrap2EndLoc
	else
		Variable/G Unwrap2 = 0
		Variable/G Unwrap2End = 0
	endif
	
	if(Unwrap3EndLoc != 0)
		Variable/G Unwrap3 = 1
		Variable/G Unwrap3End = Unwrap3EndLoc
	else
		Variable/G Unwrap3 = 0
		Variable/G Unwrap3End = 0
	endif
	
	String/G HVFilesToKillStr = HVFilesToKillStrLoc
	i = 0
	Variable NumHVFilesToKill = ItemsInList(HVFilesToKillStr)
	if(NumHVFilesToKill != 0)
		Make/O/N=(NumHVFilesToKill) HVFilesToKill
		do
			HVFilesToKill[i] = Str2Num(StringFromList(i, HVFilesToKillStr))
			i+=1
		while(i < NumHVFilesToKill)
		Variable/G KillHVFiles = 1
	else
		If(WaveExists(HVFilesToKill))
			HVFilesToKill = NaN
		EndIf
		Variable/G KillHVFiles = 0
	EndIf
	
	
	//    Reference phase unwrapping
	Variable PlusMinusRefLoc = NumVarOrDefault("plusminusref", 1)
	Prompt plusminusrefloc, "Add 2pi (1) or subtract 2pi (-1): "
	Variable UnwrapRef1Loc = NumVarOrDefault("UnwrapRef1", 0)
	Variable UnwrapRef1EndLoc = NumVarOrDefault("UnwrapRef1End", 0)
	Prompt UnwrapRef1EndLoc, "Last file for unwrapping: "
	String RefFilesToKillStrLoc  = StrVarOrDefault("RefFilesToKillStr", "")
	Prompt RefFilesToKillStrLoc, "List Ref files to kill separated by ; "
	Variable InvertPhaseLoc = NumVarOrDefault("InvertPhase",1)
	Prompt InvertPhaseLoc, "Multiply phase shift by -1? yes (1) no (0)"
	
	DoPrompt CurrentDataFolder+" Reference phase unwrapping", plusminusrefloc, UnwrapRef1EndLoc, RefFilesToKillStrLoc, InvertPhaseLoc
	
	Variable/G plusminusref = plusminusrefloc
	Variable/G InvertPhase = InvertPhaseLoc
	
	if(UnwrapRef1EndLoc != 0)
		Variable/G UnwrapRef1 = 1
		Variable/G UnwrapRef1End = UnwrapRef1EndLoc
	else
		Variable/G UnwrapRef1 = 0
		Variable/G UnwrapRef1End = 0
	endif
	
	String/G RefFilesToKillStr = RefFilesToKillStrLoc
	i = 0
	Variable NumRefFilesToKill = ItemsInList(RefFilesToKillStr)
	if(NumRefFilesToKill != 0)
		Make/O/N=(NumRefFilesToKill) RefFilesToKill
		do
			RefFilesToKill[i] = Str2Num(StringFromList(i, RefFilesToKillStr))
			i+=1
		while(i < NumRefFilesToKill)
		Variable/G KillRefFiles = 1
	else
		If(WaveExists(RefFilesToKill))
			RefFilesToKill = NaN
		EndIf
		Variable/G KillRefFiles = 0
	EndIf
	
End



Function PredictPhaseVoltage(alpha)
	Variable alpha
	//Variable position
	
	NVAR velocity, sigma, position
	
	Variable runtimer = 1
	If(runtimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Make/O/N=1 paramWave = {alpha}
	
	Make/O/N=20 HVmonWave = .1 * p, predictedPhase
	Duplicate/O HVmonWave Vwire1d, VwireSqrd
	Variable ZeroMon = .0137
	Vwire1d = GetHVKeithley(HVmonWave[p], ZeroMon)*1000
	VwireSqrd = Vwire1d^2
	
	Variable/G num_HV = numpnts(HVmonWave)
	Make/o/n=(num_HV,num_x, number_v) Vwire3d xx   						// yy: along a column we increase in y position with each row. Along a row we step along x in each column.
	xx[][][] = ((q * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_HV) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam

	Vwire3d = Vwire1d[p]

	CalculateConstantsVoltage()
	
	NormalizeVelocitiesVoltage(velocity, sigma, n)
	
	//Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	//MatrixOP/O position_shifted = offset + position					// account for y offset
	Calc_EfieldVoltage(position)							// calculate E-field at each data point
	
	CalcPhiOverPol()										// calculate phase/polarizability
	
	Phi(paramWave, predictedPhase, Vwire1d)
	
	//CleanupPolWavesVoltage()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End


// Vwire1d should be the actual voltage on the wire in units of volts (not HVmon). Use GetHVKeithley*1000 to convert HVmon into Vwire in volts. 
Function PhaseFitVoltage(VwireWave, PhaseWave, ErrorWave)
	Wave VwireWave, PhaseWave, ErrorWave
	
	NVAR velocity, sigma, position
	
	Variable runTimer = 0
	If(runTimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	//Make/O/N=20 HVmonWave = .1 * p, predictedPhase
	//Duplicate/O HVmonWave Vwire1d, VwireSqrd
	//Variable ZeroMon = .0137
	//Vwire1d = HVwave//GetHVKeithley(HVmonWave[p], ZeroMon)*1000
	Duplicate/O VwireWave VwireSqrd
	VwireSqrd = VwireWave^2
	
	Variable/G num_HV = numpnts(VwireWave)
	Make/o/n=(num_HV,num_x, number_v) Vwire3d xx   						// yy: along a column we increase in y position with each row. Along a row we step along x in each column.
	xx[][][] = ((q * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_HV) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam

	Vwire3d = VwireWave[p]

	CalculateConstantsVoltage()
	
	Variable SwitchIFM = 1
	if(SwitchIFM)
		NVAR s; 	s*=-1
	endif
	
	NormalizeVelocitiesVoltage(velocity, sigma, n)
	
	//Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	//MatrixOP/O position_shifted = offset + position					// account for y offset
	Calc_EfieldVoltage(position)							// calculate E-field at each data point
	
	CalcPhiOverPol()										// calculate phase/polarizability
	
	if(SwitchIFM)
		Wave IFMphase_v_overpol; IFMphase_v_overpol*=-1
	endif
	
	Variable LastFitPointLoc = NumVarOrDefault("LastFitPoint", 1000)		// works so long as we're not fitting more than 1000 data points (the * option does not work with NumVarOrDefault() )

	Make/D/O W_coef={43} 				// Initial Guess
	Make/D/O Epsilonwave={1e-1}			// Amount to vary the fit parameter by in each iteration	
	FuncFit/N/M=2/TBOX=768 Phi W_coef PhaseWave[0,LastFitPointLoc]  /X=VwireWave /W=ErrorWave /I=1 /R /E=EpsilonWave// /M=mask_g
	
	AppendFittedPolWaveVoltage(VwireWave, PhaseWave)				

	CleanupPolWavesVoltage()

	if(runTimer)
		Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
	endIf
End




Function CalculateConstantsVoltage()
	NVAR mass
	
	//Variable/G Vwire = GetHVKeithley(HVmon, ZeroMon)*1000
	Variable/G d = sqrt((R+gap)^2 - R^2)								// Distance from ground plane to the "infinite wire" [m]
	Variable/G lambdaOverVwire = pi*eps0 / asinh(d/R)				// Charge density of the "infinite wire"/Vwire
	Variable/G E0overVwire = -lambdaOverVwire / (pi*eps0)		
	Variable/G s = (L*hbar*2*pi)/(mass*a)								// beam separation at center of grad e region
End


Function NormalizeVelocitiesVoltage(velocity, sigma, n)			//Find the normalization for the density profile. 
	Variable velocity, sigma, n

	NVAR num_HV
	Make/o/n=(number_v) vel dens
	vel = velocity + (x * v_step - 0.5 * number_v * v_step)			//make a wave for velocities centered on "velocity"
	dens = v_step * (vel)^n* exp (-(vel-velocity)^2/(2*sigma^2)) 	// prob(v)* v_step to find area under curve
	Variable vel_normalization = area(dens)							// used to normalize Int(prob(v)) to 1
	dens /= vel_normalization
	
	Variable maxprob = WaveMax(dens)
	if(dens[0]/maxprob > .001 || dens[numpnts(dens)] > .001)
		print "VELOCITY DISTRIBUTION IS NOT WIDE ENOUGH"
		DoAlert 0, "Velocity distribution is not wide enough. Adjust number_v or v_step and rerun."
	endif
	
	if(0)
		Display/K=1 dens vs vel
		Label left "Probability"
		Label bottom "Velocity (m/s)"
	endif
		
	
	Make/o/n=(num_HV,num_x, number_v) v
	v[][][] =( (velocity) - (.5 * number_v * v_step ))+(v_step * z) // all rows and columns of a given layer have the same velocity
	MatrixOP/O density = powR(v,n)* exp (-magsqr(v-(velocity))/(2*magsqr(sigma))) / vel_normalization	//Probability density distribution. same layering as v
	MatrixOP/O vrecip = rec(v)
end


Function Calc_EFieldVoltage(position)			// Calculates the electric field of the gradient region and makes waves for phase, prob, and contrast
	Variable position						// positions at which to calculate the E-field
	//Wave HVmonWave
	Nvar num_HV, d, E0overVwire, s
	
	Wave Vwire3d, xx, xxsqr, vrecip
	//HV3d[][][] = HVmonWave[p] 			// populate the HV3d wave with the HVmons. Each row contains a different HV[p], but the HV is the same across columns and layers.
	
	// some auxiliary waves that make the E field calculations more efficient
	Duplicate/O vrecip DplusYY, DminusYY; DplusYY=d+position; DminusYY=d-position
//	MatrixOP/O DplusYY = d+position			
//	MatrixOP/O DminusYY = -1*position + d
	MatrixOP/O DplusYsep = -1*s*vrecip+d+position//d + position - (s * vrecip)
	MatrixOP/O DminusYsep =  s*vrecip - position + d 
	
	// The electric field magnitude of the primary and diffracted beams
	MatrixOP/O Etot= (E0overVwire) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O Etot_s= (E0overVwire)*sqrt(	magsqr(DplusYsep/(magsqr(DplusYsep)+xxsqr) + DminusYSep/(magsqr(DminusYSep)+xxsqr)) + magsqr(xx/(magsqr(DplusYsep)+xxsqr) - xx/(magsqr(DminusYsep)+xxsqr))		)

	// Need the magnitude squared.
	MatrixOP/O EtotMagSqr = magsqr(Etot*Vwire3d)
	MatrixOP/O Etot_sMagSqr = magsqr(Etot_s*Vwire3d)
End


//
//Function CalcPhiOverPolVoltage()
//	Wave EtotMagSqr, Etot_sMagSqr, vrecip, IFMphase 
//	// Find the phase shift for each velocity component of the main and diffracted beam
//	MatrixOP/O dphi = EtotMagSqr * vrecip							//phase shift of first path:   alpha*E^2/v	
//	MatrixOP/O dphi_s = Etot_sMagSqr * vrecip					//phase shift of second path a distance "sep" away
//	Integrate/dim=1 dphi, dphi_s 
//	MatrixOP/O IFMphase_v_overpol = dphi - dphi_s				// phase shift divided by polarizability
//		
//	// Just for curiousity (and debugging): the phase shift of the average velocity path:
//	//IFMphase = IFMphase_v_overpol[p][num_x-1][.5*number_v]		
//End
//

//
//Function PhiVoltage(pw, output, y) :FitFunc
//	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
//	
//	Wave Etot, Etot_s, IFMprob, density, Contrast, IFMphase_v_overpol
//	
//	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
//	
//	MatrixOP/O IFMphase_v = pol * IFMphase_v_overpol			// now we have the actual phase shift to take the real and imaginary parts of
//	
//	// We need to find the total velocity weighted phase shift
//	// phi_total = ArcTan( Int( Prob(v) * exp(i * delta_phi) ) ) 
//	//			 = ArcTan( Int(Prob(v)*sin(phi)) / Int(Prob(v)*cos(phi)) )
//	// the code to implement this follows:
//	
//	Variable loc_v_step = v_step											
//	MatrixOP/O imagpart = loc_v_step * density * sin(IFMphase_v)			// weight imaginary part by prob(v) 
//	MatrixOP/O realpart = loc_v_step * density * cos(IFMphase_v)				// weight real part by prob(v)
//	Integrate/DIM=2 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
//	
//	Contrast = sqrt( magsqr(imagpart[p][num_x-1][number_v-1]) + magsqr(realpart[p][num_x-1][number_v-1]) )
//	IFMprob = atan2(imagpart[p][num_x-1][number_v-1],realpart[p][num_x-1][number_v-1])
//	MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
//	
//	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
//	if(IFMprob[0]<0)
//		IFMprob+=2*pi
//	endif
//	
//	output = IFMprob
//	
//// the following is wrong because we would integrate just the phase itself over velocity (rather than exp(i phi)):
////	MatrixOP/O ProbWeightedPhase = loc_v_step * density * IFMphase_v	// weight each calculation by the appropriate amount
////	Integrate/DIM=2 ProbWeightedPhase									// must integrate along dv to find total phase shift
////	output = ProbWeightedPhase[p][num_x-1][number_v-1]				// select the column corresponding to the total integral over x and v
//End




Function CleanupPolWavesVoltage()
	KillWaves/Z dens, density, DminusYsep, DminusYY, dphi, dphi_s, DplusYsep, DplusYY, Etot, EtotMagSqr, Etot_s, Etot_sMagSqr, IFMphase_v, IFMphase_v_overpol, IFMphase, IFMprob
	KillWaves/Z imagpart, realpart, xx, xxsqr, yy, Vwire3d, v, vrecip
End




Function AppendFittedPolWaveVoltage(PositionWave, PhaseWave)
	Wave PositionWave, PhaseWave
	Wave W_coef
	

	PredictPhaseVoltage(W_coef[0])

		
	Wave Vwire1d, predictedPhase
		
	String fitPosWaveName= NameOfWave(positionWave)+"_fit"			// make some descriptive wave names
	String fitPhaseWaveName = NameOfWave(phaseWave) +"_fit"
	
	Duplicate/O Vwire1d $fitPosWaveName; Wave fitPosWave = $fitPosWaveName					// assign the calculated position and phase waves to the right wave names
	Duplicate/O predictedPhase $fitPhaseWaveName; Wave fitPhaseWave = $fitPhaseWaveName
	
	If(StringMatch(TraceInfo("", fitPhaseWaveName,0),""))		// if the trace is not already on the graph then put it there
		AppendToGraph fitPhaseWave vs fitPosWave
		ModifyGraph rgb($fitPhaseWaveName)=(1,16019,65535)
		ModifyGraph zero(Res_Left)=1
		SetAxis left 0,*
		SetAxis/A Res_Left
	EndIf
End





Function GaussTest()
	Wave predictedphasenoisy, predictedphase, vwire1d, predictedphaseerror, w_coef
	
	Variable numtests=10
	Make/o/n=(numtests) fittedpol
	
	Variable i=0
	do
	PredictPhaseVoltage(43)
	predictedPhaseNoisy=predictedphase*(1+gnoise(.02))
	PhaseFitVoltage(Vwire1d, PredictedPhaseNoisy, PredictedPhaseError)
	
	fittedpol[i]=w_coef[0]
	i+=1
	while(i<numtests)
End	




Function ExtractPhaseDataVoltage(poswave, phase, phase_error)
	Wave poswave, phase, phase_error
	
	NVAR LastEclipse, LastFile
	Variable minFile = LastEclipse+1
	Variable maxFile = LastFile		
	
	String series
	
	Variable maxIndex
	if(maxFile==0)
		maxFile = numpnts(phase)
	endif
	maxIndex = maxFile/5
	
	Variable FirstHVPosIndex = (minFile-1)/5
	String PosWaveName = NameOfWave(poswave) + "_polfit"
	Make/O/N=(maxIndex-FirstHVPosIndex) $PosWaveName
	Wave Pos = $PosWaveName
	Pos = poswave[x+FirstHVPosIndex]
	
	String PosPadWaveName = NameOfWave(poswave) + "_polfit_pad"
	Make/O/N=(maxFile-minFile+1) $PosPadWaveName
	Wave PosPad = $PosPadWaveName
	Variable n=0; Variable i=0
	do	
		PosPad[i] = pos[n]
		if(mod(i+1, 5)==0 && i>0)
			n += 1
		endif
		i+=1
	while(i<maxFile)										
		
	String PhaseWaveName = NameOfWave(phase) + "_polfit"
	String PhaseErrorWaveName = NameOfWave(phase_error) + "_polfit"

	Make/O/N=(maxFile-minFile+1) $PhaseWaveName
	Make/O/N=(maxFile-minFile+1) $PhaseErrorWaveName
	
	Wave PolFitPhase = $PhaseWaveName
	Wave PolFitPhaseError = $PhaseErrorWaveName

	PolFitPhase = phase[x+minFile-1]
	PolFitPhaseError = phase_error[x+minFile-1]
	
	String RefPhaseWaveName = NameOfWave(PolFitPhase)+"_ref"
	String RefPhaseErrorWaveName = NameOfWave(PolFitPhaseError)+"_ref"
	String RefPhaseResidWN = "Res_"+NameOfWave(PolFitPhase)+"_ref"
	
	Duplicate/O PolFitPhase $RefPhaseWaveName				//copy the phase data into a wave that will eventually only hold reference phase data
	Duplicate/O PolFitPhaseError $RefPhaseErrorWaveName
	
	Wave RefPhase = $RefPhaseWaveName
	Wave RefPhaseError = $RefPhaseErrorWaveName
	
		
		
	//	extract ref phase
	if(1)
		n=0; i=0
		do
			RefPhase[i]=NaN				//This is a HV phase, so erase it from the ref phase wave
			RefPhaseError[i]=NaN
			
			if(n==9)						// if n=9 we've gone through 10 HV phase measurements, so the next 5 will be reference phases
				n=0						
				i+=5						// so we'll skip the eraser over them in the next iteration
			else
				n+=1
			endif
			
			i+=1
		while(i<numpnts(RefPhase))
		
		i=0
		do										//remove reference phase measurements from phase data
			if(PolFitPhase[i]==RefPhase[i])	// returns true if the point is a reference phase measurement
				PolFitPhase[i]=NaN
				PolFitPhaseError[i]=NaN
			endif
			i+=1
		while(i<numpnts(PolFitPhase))
	endif



	// Insert last eclipse data points as a ref phase measurement, pad 
	InsertPoints 0, 5, $PhaseWaveName, $PhaseErrorWaveName, $RefPhaseWaveName, $RefPhaseErrorWaveName, $PosPadWaveName
	i=0
	do
		PolFitPhase[i]=NaN
		PolFitPhaseError[i]=NaN
		RefPhase[i]=phase[i+minFile-6]
		RefPhaseError[i]=phase_error[i+minFile-6]
		PosPad[i]=PosWave[FirstHVPosIndex-1]
		i+=1
	while(i<5)

	//Display/K=1/W=(40,100,550,600) 
	Display/K=1/W=(140,150,850,600) PolFitPhase
	ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
	Label left "Phase"; Label bottom "File number minus 1"
	AppendToGraph/L='RefPhase' RefPhase
	ModifyGraph mode=3,marker=8
	ErrorBars $RefPhaseWaveName Y,wave=(RefPhaseError,RefPhaseError)
	Label RefPhase "RefPhase"
	ModifyGraph axisEnab(left)={.7,1}
	ModifyGraph axisEnab(RefPhase)={0,0.67}
	ModifyGraph freePos(RefPhase)={0,bottom}
	ModifyGraph lblPosMode(refphase)=1
	String DataFolderNameInit = GetDataFolder(0)
	String DataFolderName = ReplaceString("'", DataFolderNameInit, "")
	TextBox/C/N=text0/A=MC DataFolderName
	TextBox/C/N=text0/X=9.58/Y=45.65
	
	
	
//	Keep4(PosPad)
//	Keep4(PolFitPhase)
//	Keep4(PolFitPhaseError)
//	Keep4(RefPhase)
//	Keep4(RefPhaseError)
	
	Keep3(PosPad)
	Keep3(PolFitPhase)
	Keep3(PolFitPhaseError)
	Keep3(RefPhase)
	Keep3(RefPhaseError)
	

	Variable maxUnwrapIndex, NumPointsToKill, Killindex, NumPointsToUnwrap, UnwrapIndex

	NVAR PlusMinusRef, UnwrapRef1, UnwrapRef1End, KillRefFiles//, UnwrapRefSpecial
	Wave RefFilesToKill, RefFilesToUnwrap
	
	if(UnwrapRef1)
		//Unwrap
		i=0
		maxUnwrapIndex = UnwrapRef1End-minFile+5+1
		do
			RefPhase[i] += plusminusRef*2*pi
			i+=1
		while(i<maxUnwrapIndex)
	endif
	
	
	If(KillRefFiles)
		NumPointsToKill = numpnts(RefFilesToKill)
		i = 0
		do
			Killindex = RefFilesToKill[i]-minFile+5
			RefPhase[Killindex] = NaN
			RefPhaseError[Killindex] = NaN
			i+=1
		while(i<NumPointsToKill)
	EndIf
	
//	If(UnwrapRefSpecial)
//		NumPointsToUnwrap = numpnts(RefFilesToUnwrap)
//		i = 0
//		do
//			UnwrapIndex = RefFilesToUnwrap[i]-minFile+5
//			RefPhase[UnwrapIndex] = NaN
//			i+=1
//		while(i<NumPointsToUnwrap)
//	EndIf
	

	NVAR plusminus, UnwrapAll, Unwrap1, Unwrap1End, Unwrap2, Unwrap2End, Unwrap3, Unwrap3End, AddSubAll

	If(UnwrapAll)
		PolFitPhase+=2*pi
	EndIf


	if(Unwrap1)
		//Unwrap
		//i=0
		maxUnwrapIndex = Unwrap1End-minFile+5+1
		PolFitPhase[maxUnwrapIndex,*] += plusminus*2*pi
		//do
		//	PolFitPhase[i] += plusminus*2*pi
		//	i+=1
		//while(i<maxUnwrapIndex)
	endif
	
	
	if(Unwrap2)
		//Unwrap again
		//i=0
		Variable maxUnwrapIndex2 = Unwrap2End-minFile+5+1
		PolFitPhase[maxUnwrapIndex2,*] += plusminus*2*pi
		//do
		//	PolFitPhase[i] += plusminus*2*pi
		//	i+=1
		//while(i<maxUnwrapIndex2)
	endif
	
	
	if(Unwrap3)
		//Unwrap again
		//i=0
		Variable maxUnwrapIndex3 = Unwrap3End-minFile+5+1
		PolFitPhase[maxUnwrapIndex3,*] += plusminus*2*pi
		//do
		//	PolFitPhase[i] += plusminus*2*pi
		//	i+=1
		//while(i<maxUnwrapIndex2)
	endif
	
	
	if(AddSubAll!=0)
		PolFitPhase+=AddSubAll*2*pi
	endif
	
	
	NVAR KillHVFiles
	Wave HVFilesToKill
	
	if(KillHVFiles)
		NumPointsToKill = numpnts(HVFilesToKill)
		i = 0
		do
			Killindex = HVFilesToKill[i]-minFile+5
			PolFitPhase[Killindex] = NaN
			i+=1
		while(i<NumPointsToKill)
	endif

	
	if(1)
		//Unwrap special
		PolFitPhase[103] += 2*pi
	endif




/// fix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!???????????	
	if(1)
		//Fit the reference phase to a quadratic polynomial and adjust the phase shift data via interpolation
		Make/O/N=3 w_coef
		CurveFit/NTHR=1/TBOX=768 poly 3,  RefPhase /W=RefPhaseError /I=1 /D /R
		//TextBox/C/N=CF_phase_polfit_ref/X=7.95/Y=44.71
		PolFitPhase -= w_coef[0] + w_coef[1]*x + w_coef[2]*x^2
		if(1)
		Wave RefPhaseResid = $RefPhaseResidWN
		RefPhaseResid/=RefPhase
		RefPhaseResid*=RefPhase
//		Duplicate/O RefPhaseResid RefPhaseResidTemp
//		RefPhaseResid /= RefPhaseResid
//		RefPhaseResid *= RefPhaseResidTemp
		endif
		ModifyGraph zero(Res_RefPhase)=1
		//TextBox/C/N=CF_phase_polfit_ref/X=5.40/Y=52.00
	endif

	if(0)
		Display/K=1/W=(40,100,440,325) PolFitPhase vs PosPad
		ModifyGraph mode=3,marker=8
		ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
		Label left "Phase"; Label bottom "Distance to Ground Plane (m)"
	endif
	
	
	
	NVAR InvertPhase

	if(InvertPhase)							//Make the phase shift positive
		PolfitPhase*=-1
	endif


	ReducePolWavesVoltage(PosPad,PolFitPhase, PolFitPhaseError)
End




// Reduces the raw position, phase, and error waves to just the important parts
Function ReducePolWavesVoltage(position,phase,error)			//phase OR contrast - doesn't matter
	Wave position, phase, error
	
	// Make new wave names
	String PositionReducedWaveName = NameOfWave(position) +"_red"
	String PhaseReducedWaveName= NameOfWave(phase) + "_red"
	String ErrorReducedWaveName= NameOfWave(error) + "_red"
	
	// The reduced waves are initially just copies of the raw waves
	Duplicate/O position $PositionReducedWaveName; Wave PositionReduced = $PositionReducedWaveName
	Duplicate/O phase $PhaseReducedWaveName; Wave PhaseReduced = $PhaseReducedWaveName
	Duplicate/O error $ErrorReducedWaveName; Wave ErrorReduced = $ErrorReducedWaveName
	
	// If a phase, position, or error point is a NaN then we remove it and the corresponding points in the other waves
	RemoveNaNsXYZ(PhaseReduced, PositionReduced, ErrorReduced)
	
	// Sort the waves so that point 0 is closest to the ground plane (makes unwrapping easier)
	Sort PositionReduced,PositionReduced,PhaseReduced,ErrorReduced
	
	Display/K=1/W=(40,100,440,325) PhaseReduced vs PositionReduced
	ModifyGraph mode=3,marker=8
	ErrorBars $PhaseReducedWaveName Y,wave=(ErrorReduced,ErrorReduced)
	Label left "Phase shift (rad)"; Label bottom "Vwire (V)"
	SetAxis left 0,*
	SetAxis bottom 0,*
	if(1)	// set to true for finer ticks and measure distance in mm
		ModifyGraph nticks(bottom)=10,minor(bottom)=1,prescaleExp(bottom)=-3
		Label bottom "Vwire (kV)"
		ModifyGraph minor(left)=1
		ModifyGraph sep(left)=4
	endif
	
	String DataFolderNameInit = GetDataFolder(0)
	String DataFolderName = ReplaceString("'", DataFolderNameInit, "")
	TextBox/C/N=text0/A=MC DataFolderName
	String WindowName = "PP"+DataFolderName
	DoWindow/T $WindowName, DataFolderName+": Phase vs Vwire"
	
	Print PositionReducedWaveName
	Print PhaseReducedWaveName
	Print ErrorReducedWaveName
End




Function ExtractEclipseVoltage()

	String DataFolderNameInit = GetDataFolder(0)
	String DataFolderName = ReplaceString("'", DataFolderNameInit, "")

	String datay = "EclipseCountRate"+DataFolderName	// counts
	String datax = "EclipseMotorCounts"+DataFolderName	// motor position

	LoadWave/G/D/W/N/O 
	Duplicate/O wave0 $datay; Duplicate/o wave1 $datax 		// copy counts and motor position into waves named $datay and $datax
	Wave positionWave = $datax; Wave countsWave = $datay;	// create references to the datay and datax waves
	
	positionWave-=(positionWave>1e6)*  16777216			// if the counter has "wrapped" around then correct the record
	
	positionWave /= 1e7										// convert counts from 100nm to meters

	If(1)
		Display/k=1/W=(40,100,650,400) countsWave vs positionWave	
		//ModifyGraph mode=4,marker=8,msize=2
		TextBox/C/N=text0/F=0/A=LT DataFolderName
		Label bottom "Position (m)"
		Label left "Atom Count Rate"
		ModifyGraph rgb($datay)=(0,0,65535)
	EndIf
	
	Make/O/D/N=4  W_coef={10,17000,-0.0004,1}
	FuncFit/NTHR=0/TBOX=768 errfitneg W_coef  countsWave /X=PositionWave /D 
	
	Wave w_sigma
	
	Print "offset = " + num2str(w_coef[2]*1e6)  + "±" + num2str(w_sigma[2]*1e6)+" microns"
	
	NVAR InitialGradEMotorPosition, FinalGradEMotorPosition
	
	Variable ZeroPosition = InitialGradEMotorPosition+w_coef[2]
	
	Variable/G Position = FinalGradEMotorPosition - ZeroPosition + gap
	print "position = " + num2str(position) + " meters"
	TextBox/C/N=text1/F=0/A=LT/X=(60)/Y=50 "HV position = "+num2str(position)
End