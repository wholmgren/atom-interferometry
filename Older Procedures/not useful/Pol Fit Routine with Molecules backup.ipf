#pragma rtGlobals=1		// Use modern global access method.

#include ":Pol Fit Routine"

Function PredictPhase2IFMmol(alpha, offset)
	Variable alpha
	Variable offset
	
	NVAR velocity, sigma
	
	Variable runtimer = 0
	If(runtimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Make/O/N=1 paramWave = {alpha}
	
	Variable positionOffset = 0//.0001
	Make/O/N=40 position = .05 * x/1000+positionOffset, predictedPhase
	
	Variable/G num_pos = numpnts(position)
	Make/o/n=(num_pos,num_x, number_v) yy xx   						// yy: along a column we increase in y position with each row. Along a row we step along x in each column.
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam

	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
										// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = offset + position					// account for y offset
	Calc_Efield2IFMmol(position_shifted)							// calculate E-field at each data point
	
	CalcPhiOverPol2IFMmol()										// calculate phase/polarizability
	 SagPhase()
	
	Phi2IFMmol(paramWave, predictedPhase, position)
	
	Duplicate/O Contrast, predictedContrast
	
	CleanupPolWaves()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End



Function PhaseFit2IFMmol(SeriesName)
	String SeriesName
	
	String posWaveName = "pos_"+SeriesName+"_polfit_pad_red"; Wave PositionWave = $posWaveName
	String phaseWaveName = "phase_"+SeriesName+"_polfit_red"; Wave PhaseWave = $phaseWaveName
	String phaseErrorWaveName = "phase_error_"+SeriesName+"_polfit_red"; Wave ErrorWave = $phaseErrorWaveName
	
	NVAR velocity, sigma, offset
	
	Variable runTimer = 1
	If(runTimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Variable/G num_pos = numpnts(PositionWave)
	Make/D/o/n=(num_pos,num_x, number_v) yy xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)						// The x distance squared
	Make/D/o/n=(num_pos) IFMphase IFMprob Contrast		// Waves for: total phase shift of monochrom. beam; phase shift of polychromatic beam
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = PositionWave + shift		// account for y offset
	Calc_Efield2IFMmol(position_shifted)							// calculate E-field at each data point

	CalcPhiOverPol2IFMmol()										// calculate phase/polarizability
	 SagPhase()
	
	Variable LastFitPointLoc = NumVarOrDefault("LastFitPoint", 1000)
	
	String resWaveName = "Res_"+NameOfWave(PhaseWave)
	if(WaveExists($resWaveName))
		Wave resWave = $resWaveName
		resWave = NaN
	endif
	
	Make/D/O W_coef={35} 				// Initial Guess
	Make/D/O Epsilonwave={1e-1}			// Amount to vary the fit parameter by in each iteration	
	FuncFit/N/M=2/TBOX=768 Phi2IFMmol W_coef PhaseWave[0,LastFitPointLoc]  /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave// /M=mask_g
	
	Variable reducedChiSqrd =  V_chisq/(V_npnts-(numpnts(W_coef))
	String rChiSqrdStr; sprintf rChiSqrdStr, "%.2f", reducedChiSqrd; printf "Reduced Chi-Squared: %.2f\r", reducedChiSqrd
	TextBox/C/N=ChiBox/Z=0/X=82.00/Y=30.00 "\\F'Symbol'c\\F'Geneva'\\S2\\M/dof = "+rChiSqrdStr
	
//	String contrastWaveName=ReplaceString("phase", NameOfWave(PhaseWave), "contrast")
//	Duplicate/O Contrast $contrastWaveName; Wave contrastWave = $contrastWaveName
	AppendFittedPolWave2IFMmol(PositionWave, PhaseWave)				

	CleanupPolWaves()

	if(runTimer)
		Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
	endIf
End


Function PhaseFit2IFMmolVratioFit(SeriesName)
	String SeriesName
	
	String posWaveName = "pos_"+SeriesName+"_polfit_pad_red"; Wave PositionWave = $posWaveName
	String phaseWaveName = "phase_"+SeriesName+"_polfit_red"; Wave PhaseWave = $phaseWaveName
	String phaseErrorWaveName = "phase_error_"+SeriesName+"_polfit_red"; Wave ErrorWave = $phaseErrorWaveName
	
	NVAR velocity, sigma, offset
	
	Variable runTimer = 0
	If(runTimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Variable/G num_pos = numpnts(PositionWave)
	Make/D/o/n=(num_pos,num_x, number_v) yy xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)						// The x distance squared
	Make/D/o/n=(num_pos) IFMphase IFMprob Contrast		// Waves for: total phase shift of monochrom. beam; phase shift of polychromatic beam
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = PositionWave + shift		// account for y offset
	Calc_Efield2IFMmol(position_shifted)							// calculate E-field at each data point

	CalcPhiOverPol2IFMmol()										// calculate phase/polarizability
	// SagPhase()
	
	Variable LastFitPointLoc = NumVarOrDefault("LastFitPoint", 1000)
	
	String resWaveName = "Res_"+NameOfWave(PhaseWave)
	if(WaveExists($resWaveName))
		Wave resWave = $resWaveName
		resWave = NaN
	endif
	
	Make/D/O W_coef={35,120} 				// Initial Guess
	Make/D/O Epsilonwave={1e-1,1}			// Amount to vary the fit parameter by in each iteration	
	Make/O/T Constraints = {"K1>80"}
	FuncFit/N/M=2/TBOX=768 Phi2IFMmolVRatio W_coef PhaseWave[0,LastFitPointLoc]  /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave /C=Constraints // /M=mask_g
	
	Variable reducedChiSqrd =  V_chisq/(V_npnts-(numpnts(W_coef))
	String rChiSqrdStr; sprintf rChiSqrdStr, "%.2f", reducedChiSqrd; printf "Reduced Chi-Squared: %.2f\r", reducedChiSqrd
	TextBox/C/N=ChiBox/Z=0/X=82.00/Y=30.00 "\\F'Symbol'c\\F'Geneva'\\S2\\M/dof = "+rChiSqrdStr
	
//	String contrastWaveName=ReplaceString("phase", NameOfWave(PhaseWave), "contrast")
//	Duplicate/O Contrast $contrastWaveName; Wave contrastWave = $contrastWaveName
	AppendFittedPolWave2IFMmol(PositionWave, PhaseWave)				

	CleanupPolWaves()

	if(runTimer)
		Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
	endIf
End



Function PhaseFit2IFMmolFit(SeriesName)
	String SeriesName
	
	String posWaveName = "pos_"+SeriesName+"_polfit_pad_red"; Wave PositionWave = $posWaveName
	String phaseWaveName = "phase_"+SeriesName+"_polfit_red"; Wave PhaseWave = $phaseWaveName
	String phaseErrorWaveName = "phase_error_"+SeriesName+"_polfit_red"; Wave ErrorWave = $phaseErrorWaveName
	
	NVAR velocity, sigma, offset
	
	Variable runTimer = 0
	If(runTimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Variable/G num_pos = numpnts(PositionWave)
	Make/o/n=(num_pos,num_x, number_v) yy xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)						// The x distance squared
	Make/o/n=(num_pos) IFMphase IFMprob Contrast		// Waves for: total phase shift of monochrom. beam; phase shift of polychromatic beam
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = PositionWave + shift		// account for y offset
	Calc_Efield2IFMmol(position_shifted)							// calculate E-field at each data point

	CalcPhiOverPol2IFMmol()										// calculate phase/polarizability
	
	Variable LastFitPointLoc = NumVarOrDefault("LastFitPoint", 1000)
	
	Make/D/O W_coef={24,50} 				// Initial Guess
	Make/D/O Epsilonwave={1e-1, 1}			// Amount to vary the fit parameter by in each iteration	
	FuncFit/N/M=2/TBOX=768 Phi2IFMmol2 W_coef PhaseWave[0,LastFitPointLoc]  /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave// /M=mask_g
	
	AppendFittedPolWave2IFMmol(PositionWave, PhaseWave)				

	CleanupPolWaves()

	if(runTimer)
		Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
	endIf
End



Function Calc_EField2IFMmol(positionwave)			// Calculates the electric field of the gradient region and makes waves for phase, prob, and contrast
	Wave positionwave						// positions at which to calculate the E-field
	Nvar num_pos, Vwire, d, lambda, E0, s             //; E0=0
	
	Variable sLoc = s
	
	Wave yy, xx, xxsqr, vrecip
	yy[][][] = positionwave[p]				// populate the yy wave with the positions. Each row contains a different position[p], but initially the position is the same across columns and layers.
											// vrecip is the same across rows and columns, but different across layers
											// the beam propagation direction is along a single row
	
	// zeroth order
	MatrixOP/O DplusYY = d+yy			
	MatrixOP/O DminusYY = -1*yy + d
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqr0 = magsqr(Etot)
	
	// minus ifm
	MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
	MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrM = magsqr(Etot)

	// plus ifm
	sLoc*=-1
	MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
	MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrP = magsqr(Etot)
	
	// plus ifm molecules
	sLoc*=0.5
	MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
	MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrPmol = magsqr(Etot)
	
	// minus ifm molecules
	sLoc*=-1
	MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
	MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrMmol = magsqr(Etot)
End



Function CalcPhiOverPol2IFMmol()
	Wave EtotMagSqr0, EtotMagSqrM, EtotMagSqrP, vrecip, IFMphase, EtotMagSqrMmol, EtotMagSqrPmol
	// Find the phase shift for each velocity component of the main and diffracted beam
	MatrixOP/O dphi = EtotMagSqr0 * vrecip							//phase shift of first path:   alpha*E^2/v	// this is the +1 diffraction order
	MatrixOP/O dphi_sM = EtotMagSqrM * vrecip					//phase shift of second path a distance "sep" away		// this is the 0th order
	MatrixOP/O dphi_sP = EtotMagSqrP * vrecip				//phase shift of third path a distance "sep" away		// this is the -1 order
	MatrixOP/O dphiMmol = EtotMagSqrMmol * vrecip
	MatrixOP/O dphiPmol = EtotMagSqrPmol * vrecip
	
	Integrate/dim=1 dphi, dphi_sM, dphi_sP, dphiMmol, dphiPmol		// integrate along beam path to find total acquired phase
	MatrixOP/O IFMphase_v_overpolM = dphi - dphi_sM				// phase shift divided by polarizability for minus ifm
	MatrixOP/O IFMphase_v_overpolP = dphi_sP - dphi			// phase shift divided by polarizability for plus ifm
	MatrixOP/O IFMphase_v_overpolMmol = dphi - dphiMmol				// phase shift divided by polarizability for minus ifm
	MatrixOP/O IFMphase_v_overpolPmol = dphiPmol - dphi			// phase shift divided by polarizability for plus ifm
	
End


Function SagPhase()
	Wave vel, dens, vrecip
	Variable loc_v_step = v_step
	
	Variable Latitude = 32.2	//degrees
	Variable OmegaEarth = 2*pi/(24*3600)*cos((90-Latitude)*pi/180)	// earth rotation rate for a given latitude
	Variable L1g2g = 0.94		// meters
	Variable SagFactor = OmegaEarth*L1g2g^2*4*pi/1e-7
	MatrixOP/O phiSag = SagFactor *rec( vel)
//	MatrixOP/O phiSagAvg = phiSag*dens; Variable/G AvgSagPhi = area(phiSagAvg); printf "Velocity averaged Sagnac phase shift: %2.7f\r", avgsagphi
	MatrixOP/O phiSagImag = dens * sin(phiSag)
	MatrixOP/O phiSagReal = dens * cos(phiSag)
	Variable phiSagImagInt = area(phiSagImag)
	Variable phiSagRealInt = area(phiSagReal)
	
	Variable/G AvgSagPhi = atan2(phiSagImagInt, phiSagRealInt)
	Variable/G AvgSagCon = sqrt(phiSagImagInt^2+phiSagRealInt^2)
	//printf "Velocity averaged Sagnac phase shift: %2.7f\r", avgsagphi; printf "Velocity averaged Sagnac contrast: %1.3f\r", avgsagcon
	
	Duplicate/O vrecip phiSagLayered
	phiSagLayered = phiSag[z]
End




Function Phi2IFMmol(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
	
	Wave IFMprob, density, Contrast, IFMphase_v_overpolM, IFMphase_v_overpolP, IFMphase_v_overpolMmol, IFMphase_v_overpolPmol
	NVAR molFrac, molPol
	
	Variable molFracLoc = molFrac/100
	
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	Variable polmol = molPol*4*pi*eps0*1E-30* x_step / (2 * hbar)
	//print molpol
	
	if(0)
		// phi(v) for plus and minus interferometers for atoms and molecules
		MatrixOP/O IFMphase_vM = pol * IFMphase_v_overpolM			// now we have the actual phase shift to take the real and imaginary parts of
		MatrixOP/O IFMphase_vP = pol * IFMphase_v_overpolP
		MatrixOP/O IFMphase_vMmol =  polmol * IFMphase_v_overpolMmol
		MatrixOP/O IFMphase_vPmol = polmol * IFMphase_v_overpolPmol
	else
		Wave phiSagLayered
		MatrixOP/O IFMphase_vM = pol * IFMphase_v_overpolM+phiSagLayered			// now we have the actual phase shift to take the real and imaginary parts of
		MatrixOP/O IFMphase_vP = pol * IFMphase_v_overpolP+phiSagLayered
		MatrixOP/O IFMphase_vMmol =  polmol * IFMphase_v_overpolMmol+phiSagLayered
		MatrixOP/O IFMphase_vPmol = polmol * IFMphase_v_overpolPmol+phiSagLayered
	endif
	
	//Variable loc_v_step = v_step	
	//0.5 to account for two IFMs										
	MatrixOP/O imagpart = 0.5  * density * (	(1-molFracLoc)*(sin(IFMphase_vM) + sin(IFMphase_vP))    +    molFracLoc*(sin(IFMphase_vMmol) + sin(IFMphase_vPmol))	 )			// weight imaginary part by prob(v) 
	MatrixOP/O realpart = 0.5  * density * (		(1-molFracLoc)*(cos(IFMphase_vM) + cos(IFMphase_vP))    +    molFracLoc*(cos(IFMphase_vMmol) + cos(IFMphase_vMmol)) 	 )				// weight real part by prob(v)
	Integrate/DIM=2 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][num_x-1][number_v-1]) + magsqr(realpart[p][num_x-1][number_v-1]) )
	IFMprob = atan2(imagpart[p][num_x-1][number_v-1], realpart[p][num_x-1][number_v-1])
	MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
	if(IFMprob[0]<0)
		IFMprob+=2*pi
	endif
	
	NVAR AvgSagPhi
	output = IFMprob - AvgSagPhi
End




Function Phi2IFMmolVratio(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
	
	NVAR velocity, sigma;  sigma= pw[1]
	NormalizeVelocities(velocity, sigma, 3)
	CalcPhiOverPol2IFMmol()										// calculate phase/polarizability
	SagPhase()
	
	Wave IFMprob, density, Contrast, IFMphase_v_overpolM, IFMphase_v_overpolP, IFMphase_v_overpolMmol, IFMphase_v_overpolPmol
	NVAR molFrac, molPol
	
	Variable molFracLoc = molFrac/100
	
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	Variable polmol = molPol*4*pi*eps0*1E-30* x_step / (2 * hbar)
	//print molpol
	
	if(0)
		// phi(v) for plus and minus interferometers for atoms and molecules
		MatrixOP/O IFMphase_vM = pol * IFMphase_v_overpolM			// now we have the actual phase shift to take the real and imaginary parts of
		MatrixOP/O IFMphase_vP = pol * IFMphase_v_overpolP
		MatrixOP/O IFMphase_vMmol =  polmol * IFMphase_v_overpolMmol
		MatrixOP/O IFMphase_vPmol = polmol * IFMphase_v_overpolPmol
	else
		Wave phiSagLayered
		MatrixOP/O IFMphase_vM = pol * IFMphase_v_overpolM+phiSagLayered			// now we have the actual phase shift to take the real and imaginary parts of
		MatrixOP/O IFMphase_vP = pol * IFMphase_v_overpolP+phiSagLayered
		MatrixOP/O IFMphase_vMmol =  polmol * IFMphase_v_overpolMmol+phiSagLayered
		MatrixOP/O IFMphase_vPmol = polmol * IFMphase_v_overpolPmol+phiSagLayered
	endif

	Variable loc_v_step = v_step	
	//0.5 to account for two IFMs										
	MatrixOP/O imagpart = 0.5 * loc_v_step * density * (	(1-molFracLoc)*(sin(IFMphase_vM) + sin(IFMphase_vP))    +    molFracLoc*(sin(IFMphase_vMmol) + sin(IFMphase_vPmol))	 )			// weight imaginary part by prob(v) 
	MatrixOP/O realpart = 0.5 * loc_v_step * density * (		(1-molFracLoc)*(cos(IFMphase_vM) + cos(IFMphase_vP))    +    molFracLoc*(cos(IFMphase_vMmol) + cos(IFMphase_vPmol)) 	 )				// weight real part by prob(v)
	Integrate/DIM=2 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][num_x-1][number_v-1]) + magsqr(realpart[p][num_x-1][number_v-1]) )
	IFMprob = atan2(imagpart[p][num_x-1][number_v-1], realpart[p][num_x-1][number_v-1])
	MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
	if(IFMprob[0]<0)
		IFMprob+=2*pi
	endif
	
	output = IFMprob
End




//fits both alpha atom and alpha mol:
Function Phi2IFMmol2(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
	
	Wave IFMprob, density, Contrast, IFMphase_v_overpolM, IFMphase_v_overpolP, IFMphase_v_overpolMmol, IFMphase_v_overpolPmol
	NVAR molFrac, molPol
	
	Variable molFracLoc = molFrac/100
	
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	Variable polmol = pw[1]*4*pi*eps0*1E-30* x_step / (2 * hbar)
	//print molpol
	
	// phi(v) for plus and minus interferometers for atoms and molecules
	MatrixOP/O IFMphase_vM = pol * IFMphase_v_overpolM			// now we have the actual phase shift to take the real and imaginary parts of
	MatrixOP/O IFMphase_vP = pol * IFMphase_v_overpolP
	MatrixOP/O IFMphase_vMmol =  polmol * IFMphase_v_overpolMmol
	MatrixOP/O IFMphase_vPmol = polmol * IFMphase_v_overpolPmol

	Variable loc_v_step = v_step											
	MatrixOP/O imagpart = loc_v_step * density * (	(1-molFracLoc)*(sin(IFMphase_vM) + sin(IFMphase_vP))    +    molFracLoc*(sin(IFMphase_vMmol) + sin(IFMphase_vPmol))	 )			// weight imaginary part by prob(v) 
	MatrixOP/O realpart = loc_v_step * density * (		(1-molFracLoc)*(cos(IFMphase_vM) + cos(IFMphase_vP))    +    molFracLoc*(cos(IFMphase_vMmol) + cos(IFMphase_vMmol)) 	 )				// weight real part by prob(v)
	Integrate/DIM=2 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][num_x-1][number_v-1]) + magsqr(realpart[p][num_x-1][number_v-1]) )
	IFMprob = atan2(imagpart[p][num_x-1][number_v-1], realpart[p][num_x-1][number_v-1])
	MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
	if(IFMprob[0]<0)
		IFMprob+=2*pi
	endif
	
	output = IFMprob
End



Function AppendFittedPolWave2IFMmol(PositionWave, PhaseWave)
	Wave PositionWave, PhaseWave
	Wave W_coef
	
	NVAR offset
	
	PredictPhase2IFMmol(W_coef[0], offset)			// Run the predict phase routine with the fitted polarizability and fixed offset
	
	Wave position, predictedPhase, predictedContrast
		
	String fitPosWaveName= NameOfWave(positionWave)+"_fit"			// make some descriptive wave names
	String fitPhaseWaveName = NameOfWave(phaseWave) +"_fit"
	String fitContrastWaveName = ReplaceString("phase", fitPhaseWaveName, "contrast")//NameOfWave(contrastWave) + "_fit"
	
	Duplicate/O position $fitPosWaveName; Wave fitPosWave = $fitPosWaveName					// assign the calculated position and phase waves to the right wave names
	Duplicate/O predictedPhase $fitPhaseWaveName; Wave fitPhaseWave = $fitPhaseWaveName
	Duplicate/O predictedContrast $fitContrastWaveName; Wave fitContrastWave = $fitContrastWaveName
	
	If(StringMatch(TraceInfo("", fitPhaseWaveName,0),""))		// if the trace is not already on the graph then put it there
		AppendToGraph fitPhaseWave vs fitPosWave
		ModifyGraph rgb($fitPhaseWaveName)=(1,16019,65535)
		ModifyGraph zero(Res_Left)=1
		SetAxis left 0,*
		SetAxis/A Res_Left
	EndIf
	
	If(StringMatch(TraceInfo("", fitcontrastWaveName,0),""))		// if the trace is not already on the graph then put it there
		AppendToGraph/l=contrast fitcontrastWave vs fitPosWave
		ModifyGraph rgb($fitcontrastWaveName)=(1,16019,65535)
	EndIf
End