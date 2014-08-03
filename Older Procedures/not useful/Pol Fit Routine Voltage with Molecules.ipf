#pragma rtGlobals=1		// Use modern global access method.

#include ":Pol Fit Routine Voltage"
#include ":Pol Fit Routine with Molecules"


Function PredictPhaseVoltage2IFMmol(alpha)
	Variable alpha
	//Variable position
	
	NVAR velocity, sigma, position
	
	Variable runtimer = 1
	If(runtimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Make/O/N=1 paramWave = {alpha}
	
	Make/O/N=20 HVmonWave = .06 * p, predictedPhase
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
	
	Calc_EfieldVoltage2IFMmol(position)							// calculate E-field at each data point
	
	CalcPhiOverPol2IFMmol()										// calculate phase/polarizability
	
	Phi2IFMmol(paramWave, predictedPhase, Vwire1d)
	
	CleanupPolWaves()
	CleanupPolWavesVoltage()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End





//Function PhaseFitVoltage2IFMmol(VwireWave, PhaseWave, ErrorWave)
//	Wave VwireWave, PhaseWave, ErrorWave
Function PhaseFitVoltage2IFMmol(SeriesName)
	String SeriesName
	
	String posWaveName = "Vwire_"+SeriesName+"_polfit_pad_red"; Wave VwireWave = $posWaveName
	String phaseWaveName = "phase_"+SeriesName+"_polfit_red"; Wave PhaseWave = $phaseWaveName
	String phaseErrorWaveName = "phase_error_"+SeriesName+"_polfit_red"; Wave ErrorWave = $phaseErrorWaveName

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
	Make/o/n=(num_HV) IFMphase IFMprob Contrast					//Total phase shift of monochrom. beam; phase shift of polychromatic beam

	Vwire3d = VwireWave[p]

	CalculateConstantsVoltage()
	
	NormalizeVelocitiesVoltage(velocity, sigma, n)
	
	Calc_EfieldVoltage2IFMmol(position)							// calculate E-field at each data point
	
	CalcPhiOverPol2IFMmol()										// calculate phase/polarizability
	
	Variable LastFitPointLoc = NumVarOrDefault("LastFitPoint", 1000)		// works so long as we're not fitting more than 1000 data points (the * option does not work with NumVarOrDefault() )

	Make/D/O W_coef={47} 				// Initial Guess
	Make/D/O Epsilonwave={1e-1}			// Amount to vary the fit parameter by in each iteration	
	FuncFit/N/M=2/TBOX=768 Phi2IFMmol W_coef PhaseWave[0,LastFitPointLoc]  /X=VwireWave /W=ErrorWave /I=1 /R /E=EpsilonWave// /M=mask_g
	
	AppendFittedPolWaveVolt2IFMmol(VwireWave, PhaseWave)				

	CleanupPolWaves()
	CleanupPolWavesVoltage()

	if(runTimer)
		Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
	endIf
End




Function Calc_EFieldVoltage2IFMmol(position)			// Calculates the electric field of the gradient region and makes waves for phase, prob, and contrast
	Variable position						// positions at which to calculate the E-field
	//Wave HVmonWave
	Nvar num_HV, d, E0overVwire, s		// d = distance from ground plane to infinitly thin effective wire
											// s = (beam separation)*velocity at center of GradE region. so s*vrecip = beam separation
	
	Wave Vwire3d, xx, xxsqr, vrecip
	//HV3d[][][] = HVmonWave[p] 			// populate the HV3d wave with the HVmons. Each row contains a different HV[p], but the HV is the same across columns and layers.
	
	// zeroth order
	Duplicate/O vrecip DplusYY, DminusYY; 
	DplusYY=d+position
	DminusYY=d-position
	MatrixOP/O Etot= (E0overVwire) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqr0 = magsqr(Etot*Vwire3d)
	
	// minus ifm
	MatrixOP/O DplusYY = -1*s*vrecip+d+position//     d + position - (|s| * vrecip)
	MatrixOP/O DminusYY =  s*vrecip - position + d 
	MatrixOP/O Etot= (E0overVwire) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrM = magsqr(Etot*Vwire3d)
	
	// plus ifm
	s*=-1
	MatrixOP/O DplusYY = -1*s*vrecip+d+position//d + position + (|s| * vrecip)
	MatrixOP/O DminusYY =  s*vrecip - position + d 	// d - position - |s| * vrecip
	MatrixOP/O Etot= (E0overVwire) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrP = magsqr(Etot*Vwire3d)
	
	// plus ifm molecules
	s*=0.5
	MatrixOP/O DplusYY = -1*s*vrecip+d+position//d + position + (|s| * vrecip)
	MatrixOP/O DminusYY =  s*vrecip - position + d 	// d - position - |s| * vrecip
	MatrixOP/O Etot= (E0overVwire) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrPmol = magsqr(Etot*Vwire3d)
	
	// minus ifm molecules
	s*=-1
	MatrixOP/O DplusYY = -1*s*vrecip+d+position//d + position + (|s| * vrecip)
	MatrixOP/O DminusYY =  s*vrecip - position + d 	// d - position - |s| * vrecip
	MatrixOP/O Etot= (E0overVwire) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrMmol = magsqr(Etot*Vwire3d)
End










Function AppendFittedPolWaveVolt2IFMmol(PositionWave, PhaseWave)
	Wave PositionWave, PhaseWave
	Wave W_coef
	

	PredictPhaseVoltage2IFMmol(W_coef[0])

		
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