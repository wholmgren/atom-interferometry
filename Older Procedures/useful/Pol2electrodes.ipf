#pragma rtGlobals=1		// Use modern global access method.

#include ":PolAnalytic"


Function PredictPhase2Electrodes(alpha, offset)
	Variable alpha
	Variable offset
	
	NVAR velocity, sigma
	
	Variable runtimer = 0
	If(runtimer)
		Variable timerRefNum = Startmstimer
	EndIf

	
	Make/D/O/N=1 paramWave = {alpha}
	
	Variable positionOffset = 0//.0001
	Make/D/O/N=80 position = -.002+.05 * x/1000+positionOffset, predictedPhase

	
	Variable/G num_pos = numpnts(position)
	MatrixOP/O position_shifted = offset + position					// account for y offset
	Make/D/o/n=(num_pos, number_v) yy    						// yy: along a column we increase in y position with each row. Along a row we step along x in each column.
	yy[][] = position_shifted[p]	
	//xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	//MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/D/n=(num_pos) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam


	CalculateConstants2Electrodes()
	
	NormalizeVelocitiesAn(velocity, sigma, n)

	PolPhaseAn()

	SagnacAndGravityPhase()
	
	Phi2Electrodes(paramWave, predictedPhase, position)
	
	Duplicate/O Contrast, predictedContrast
	
	CleanupPolWaves()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End




Function Phi2Electrodes(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam

	//DFREF saveDFR = GetDataFolderDFR(); SetDataFolder root:'090625g'	// Uncomment when performing a GlobalFit with ContrastVratio()
	
	Wave IFMprob, density2D, Contrast, IFMphase_v_overpolM, IFMphase_v_overpolP, IFMphase_v_overpolMmol, IFMphase_v_overpolPmol
	NVAR molFrac, molPol, AnalyticPrefactor
	
	Variable molFracLoc = molFrac/100
	
	Variable pol = pw[0]*AnalyticPrefactor// pw[0] is alpha in cgs units but without the factor of 10^-24
	Variable polmol = molPol*AnalyticPrefactor
	//print molpol
	
	Wave phiSagGravLayered2D
	MatrixOP/O IFMphase_vM = pol * IFMphase_v_overpolM+phiSagGravLayered2D			// now we have the actual phase shift to take the real and imaginary parts of
	MatrixOP/O IFMphase_vP = pol * IFMphase_v_overpolP+phiSagGravLayered2D
	MatrixOP/O IFMphase_vMmol =  polmol * IFMphase_v_overpolMmol+phiSagGravLayered2D
	MatrixOP/O IFMphase_vPmol = polmol * IFMphase_v_overpolPmol+phiSagGravLayered2D
	
	
	//0.5 to account for two IFMs										
	MatrixOP/O imagpart = 0.5  * density2D * (	(1-molFracLoc)*(sin(IFMphase_vM) + sin(IFMphase_vP))    +    molFracLoc*(sin(IFMphase_vMmol) + sin(IFMphase_vPmol))	 )			// weight imaginary part by prob(v) 
	MatrixOP/O realpart =  0.5  * density2D * (	(1-molFracLoc)*(cos(IFMphase_vM) + cos(IFMphase_vP))    +    molFracLoc*(cos(IFMphase_vMmol) + cos(IFMphase_vPmol)) 	 )				// weight real part by prob(v)
	Integrate/DIM=1 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][number_v-1]) + magsqr(realpart[p][number_v-1]) )
	IFMprob = atan2(imagpart[p][number_v-1], realpart[p][number_v-1])
	MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	Variable n = NumVarOrDefault("k", 0)
	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
	if(IFMprob[0]-2*pi*n<0)
		IFMprob+=2*pi//*(k+1)		// Correction necessary if Sag. phase is large, and thus the "0" point is in the wrong domain.
	endif
	
	NVAR AvgSagGravPhi
	output = IFMprob - AvgSagGravPhi
	
	FindLevel/Q y 0
	Variable zeroPhase = output[V_LevelX]
	Variable piWraps = round(zeroPhase / (2*pi))
	output -= piWraps*(2*pi)
	
	//SetDataFolder saveDFR
End



Function CalculateConstants2Electrodes()
//	Variable/G Vwire = (6035 * HVmon)	 +0.7							// Applied electrode voltage [V]
	
	NVAR HVmon, ZeroMon, mass
	
	Variable gapLoc = 2e-3
	Variable rLoc = 25.4/2/2*1e-3
	
	Variable/G Vwire = GetHVKeithley(HVmon, ZeroMon)*1000/2//; Vwire=9000
	Variable/G d = gapLoc*sqrt(1+2*rLoc/gapLoc)//sqrt((R+gapLoc)^2 - R^2)								// Distance from ground plane to the "infinite wire" [m]
	Variable/G lambda = 2*pi*eps0*vwire*(ln((gapLoc+rLoc+d)/(gapLoc+rLoc-d)))^-1 //Vwire*pi*eps0 / asinh(d/R)				// Charge density of the "infinite wire"
	Variable/G E0 = -lambda / (pi*eps0)		
	Variable/G s = (L*hbar*2*pi)/(mass*a)								// beam separation at center of grad e region * velocity
	Variable/G AnalyticPrefactor=lambda^2*d/(hbar*pi*eps0^2)*4*pi*eps0*1E-30
	//Variable/G AnalyticPrefactor=4*pi*Vwire^2*d/hbar*((ln((gap+R+d)/(gap+R-d)))^-2)*4*pi*eps0*1E-30
	//Variable lambda2 = Vwire*2*pi*eps0/ln((gap+R+d)/(gap+R-d))
	//print lambda
	//print lambda2
	//print lambda2*2
	
	//print 43*AnalyticPrefactor/(d^2-.0017^2)/1800
	//print 43*AnalyticPrefactor/(d^2-.0017^2)/1800- 43*AnalyticPrefactor/(d^2-(.0017+50e-6)^2)/1800
End