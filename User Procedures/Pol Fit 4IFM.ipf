#pragma rtGlobals=1		// Use modern global access method.

#include ":Pol Fit Routine with Molecules"

Function PredictPhase4IFM(alpha, offset)
	Variable alpha
	Variable offset
	
	NVAR velocity, sigma
	
	Variable runtimer = 0
	If(runtimer)
		Variable timerRefNum = Startmstimer
	EndIf

	
	Make/D/O/N=1 paramWave = {alpha}
	
	Variable positionOffset = 0//.0001
	Make/D/O/N=40 position = .05 * x/1000+positionOffset, predictedPhase
	
	Variable/G num_pos = numpnts(position)
	Make/D/o/n=(num_pos,num_x, number_v) yy xx   						// yy: along a column we increase in y position with each row. Along a row we step along x in each column.
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/D/n=(num_pos) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam

	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	MatrixOP/O position_shifted = offset + position					// account for y offset

	Calc_Efield4IFM(position_shifted)							// calculate E-field at each data point

	CalcPhiOverPol4IFM()										// calculate phase/polarizability
	//SagPhase()
	//GravityPhase()
	SagnacAndGravityPhase()
	
	Phi4IFM(paramWave, predictedPhase, position)
	
	Duplicate/O Contrast, predictedContrast
	
	CleanupPolWaves()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End



Function Calc_EField4IFM(positionwave)			// Calculates the electric field of the gradient region and makes waves for phase, prob, and contrast
	Wave positionwave						// positions at which to calculate the E-field
	Nvar num_pos, Vwire, d, lambda, E0, s             //; E0=0
	
	Variable sLoc = s
	
	Wave yy, xx, xxsqr, vrecip
	yy[][][] = positionwave[p]				// populate the yy wave with the positions. Each row contains a different position[p], but initially the position is the same across columns and layers.
											// vrecip is the same across rows and columns, but different across layers
											// the beam propagation direction is along a single row
	
	//print "s: ", s
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
	
	
	sLoc=2*s
	// minus ifm
	MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
	MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrM2 = magsqr(Etot)

	// plus ifm
	sLoc*=-1
	MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
	MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrP2 = magsqr(Etot)
	
	// plus ifm molecules
	sLoc*=0.5
	MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
	MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrP2mol = magsqr(Etot)
	
	// minus ifm molecules
	sLoc*=-1
	MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
	MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrM2mol = magsqr(Etot)
End





Function CalcPhiOverPol4IFM()
	Wave vrecip, IFMphase, EtotMagSqr0, EtotMagSqrM, EtotMagSqrP, EtotMagSqrMmol, EtotMagSqrPmol, EtotMagSqrM2, EtotMagSqrP2, EtotMagSqrM2mol, EtotMagSqrP2mol
	// Find the phase shift for each velocity component of the main and diffracted beam
	MatrixOP/O dphi = EtotMagSqr0 * vrecip							//phase shift of first path:   alpha*E^2/v	// this is the 0th diffraction order
	MatrixOP/O dphi_sM = EtotMagSqrM * vrecip					//phase shift of second path a distance "sep" away		// this is the -1st order
	MatrixOP/O dphi_sP = EtotMagSqrP * vrecip				//phase shift of third path a distance "sep" away		// this is the +1 order
	MatrixOP/O dphiMmol = EtotMagSqrMmol * vrecip
	MatrixOP/O dphiPmol = EtotMagSqrPmol * vrecip
	MatrixOP/O dphi_sM2 = EtotMagSqrM2 * vrecip					//2nd order IFMs
	MatrixOP/O dphi_sP2 = EtotMagSqrP2 * vrecip				
	MatrixOP/O dphiM2mol = EtotMagSqrM2mol * vrecip
	MatrixOP/O dphiP2mol = EtotMagSqrP2mol * vrecip
	
	Integrate/dim=1 dphi, dphi_sM, dphi_sP, dphiMmol, dphiPmol, dphi_sM2, dphi_sP2, dphiM2mol, dphiP2mol		// integrate along beam path to find total acquired phase
	
	NVAR num_pos
	Make/D/O/N=(num_pos, number_v) IFMphase_v_overpolM, IFMphase_v_overpolP, IFMphase_v_overpolMmol, IFMphase_v_overpolPmol, IFMphase_v_overpolM2, IFMphase_v_overpolP2, IFMphase_v_overpolM2mol, IFMphase_v_overpolP2mol
	
	// these waves are proportional to phi(x[row], v[column])
	// later we'll take the real and imaginary parts of them and integrate across columns
	IFMphase_v_overpolM = dphi[p][100][q] - dphi_sM[p][100][q]
	IFMphase_v_overpolP = dphi_sP[p][100][q] - dphi[p][100][q]
	IFMphase_v_overpolMmol = dphi[p][100][q] - dphiMmol[p][100][q]
	IFMphase_v_overpolPmol = dphiPmol[p][100][q] - dphi[p][100][q]
	IFMphase_v_overpolM2 = dphi_sM[p][100][q] - dphi_sM2[p][100][q]
	IFMphase_v_overpolP2 = dphi_sP2[p][100][q] - dphi_sP[p][100][q]
	IFMphase_v_overpolM2mol = dphiMmol[p][100][q] - dphiM2mol[p][100][q]
	IFMphase_v_overpolP2mol = dphiP2mol[p][100][q] - dphiPmol[p][100][q]
	
//	MatrixOP/O IFMphase_v_overpolM = dphi - dphi_sM				// phase shift divided by polarizability for minus ifm
//	MatrixOP/O IFMphase_v_overpolP = dphi_sP - dphi			// phase shift divided by polarizability for plus ifm
//	MatrixOP/O IFMphase_v_overpolMmol = dphi - dphiMmol				// phase shift divided by polarizability for minus ifm
//	MatrixOP/O IFMphase_v_overpolPmol = dphiPmol - dphi			// phase shift divided by polarizability for plus ifm
	
End





Function Phi4IFM(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam

	//DFREF saveDFR = GetDataFolderDFR(); SetDataFolder root:'090625g'	// Uncomment when performing a GlobalFit with ContrastVratio()
	
	Wave IFMprob, density, density2D, Contrast, IFMphase_v_overpolM, IFMphase_v_overpolP, IFMphase_v_overpolMmol, IFMphase_v_overpolPmol, IFMphase_v_overpolM2, IFMphase_v_overpolP2, IFMphase_v_overpolM2mol, IFMphase_v_overpolP2mol
	NVAR molFrac, molPol
	
	Variable molFracLoc = molFrac/100
	
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	Variable polmol = molPol*4*pi*eps0*1E-30* x_step / (2 * hbar)
	//print molpol
	
		// phi(v) for plus and minus interferometers for atoms and molecules
		
	Wave phiSagGravLayered2D
	MatrixOP/O IFMphase_vM = pol * IFMphase_v_overpolM+phiSagGravLayered2D			// now we have the actual phase shift to take the real and imaginary parts of
	MatrixOP/O IFMphase_vP = pol * IFMphase_v_overpolP+phiSagGravLayered2D
	MatrixOP/O IFMphase_vMmol =  polmol * IFMphase_v_overpolMmol+phiSagGravLayered2D
	MatrixOP/O IFMphase_vPmol = polmol * IFMphase_v_overpolPmol+phiSagGravLayered2D
	MatrixOP/O IFMphase_vM2 = pol * IFMphase_v_overpolM2+phiSagGravLayered2D			// now we have the actual phase shift to take the real and imaginary parts of
	MatrixOP/O IFMphase_vP2 = pol * IFMphase_v_overpolP2+phiSagGravLayered2D
	MatrixOP/O IFMphase_vM2mol =  polmol * IFMphase_v_overpolM2mol+phiSagGravLayered2D
	MatrixOP/O IFMphase_vP2mol = polmol * IFMphase_v_overpolP2mol+phiSagGravLayered2D
	
	
	//0.5 to account for two IFMs
	Variable SecondOrderIFMweight = 0.1										
	MatrixOP/O imagpart = 0.5  * density2D * ( (1-SecondOrderIFMweight)* (	(1-molFracLoc)*(sin(IFMphase_vM) + sin(IFMphase_vP))    +    molFracLoc*(sin(IFMphase_vMmol) + sin(IFMphase_vPmol)))	+ SecondOrderIFMweight* ( 	(1-molFracLoc)*(sin(IFMphase_vM2) + sin(IFMphase_vP2))    +    molFracLoc*(sin(IFMphase_vM2mol) + sin(IFMphase_vP2mol)) )	)		// weight imaginary part by prob(v) 
	MatrixOP/O realpart =  0.5  * density2D * ( (1-SecondOrderIFMweight)* ( (1-molFracLoc)*(cos(IFMphase_vM) + cos(IFMphase_vP))    +    molFracLoc*(cos(IFMphase_vMmol) + cos(IFMphase_vPmol)))	+ SecondOrderIFMweight* (		(1-molFracLoc)*(cos(IFMphase_vM2) + cos(IFMphase_vP2))    +    molFracLoc*(cos(IFMphase_vM2mol) + cos(IFMphase_vP2mol)) )	)			// weight real part by prob(v)
	Integrate/DIM=1 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][number_v-1]) + magsqr(realpart[p][number_v-1]) )
	IFMprob = atan2(imagpart[p][number_v-1], realpart[p][number_v-1])
	MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	NVAR k
	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
	if(IFMprob[0]-2*pi*k<0)
		IFMprob+=2*pi//*(k+1)
	endif
	
	NVAR AvgSagGravPhi
	output = IFMprob - AvgSagGravPhi
	
	//SetDataFolder saveDFR
End
