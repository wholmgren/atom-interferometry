#pragma rtGlobals=1		// Use modern global access method.

#include ":Pol Fit Routine with Molecules"


Function PredictPhase2IFMmol_iso(alpha, offset)
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

	CalculateConstants_Isotopes()
	
	NormalizeVelocities(velocity, sigma, n)
	
	MatrixOP/O position_shifted = offset + position					// account for y offset

	Calc_Efield2IFMmol_iso(position_shifted)							// calculate E-field at each data point

	CalcPhiOverPol2IFMmol_iso()										// calculate phase/polarizability
	//SagPhase()
	//GravityPhase()
	SagnacAndGravityPhase()
	
	Phi2IFMmol_iso(paramWave, predictedPhase, position)
	
	Duplicate/O Contrast, predictedContrast
	
	CleanupPolWaves()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End



Function Calc_EField2IFMmol_iso(positionwave)			// Calculates the electric field of the gradient region and makes waves for phase, prob, and contrast
	Wave positionwave						// positions at which to calculate the E-field
	Nvar num_pos, Vwire, d, lambda, E0, s             //; E0=0
	
	Wave SeparationWave
	Variable sLoc
	
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
	
	
	String EsqrdWaveName
	Make/O/WAVE/N=(numpnts(SeparationWave),4) EsqrdSepWaveRefs
	Variable i = 0, j=0
	do
		// For the diffracted beam now:
		sLoc = SeparationWave[i]
		// minus ifm
		MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
		MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
		MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
		MatrixOP/O Etot_sMagSqr = magsqr(Etot)
		EsqrdWaveName = "EtotMagSqrM" + num2str(i)
		Duplicate/O Etot_sMagSqr $EsqrdWaveName
		EsqrdSepWaveRefs[i][0] = $EsqrdWaveName
	
		// plus ifm
		sLoc*=-1
		MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
		MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
		MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
		MatrixOP/O Etot_sMagSqr = magsqr(Etot)
		EsqrdWaveName = "EtotMagSqrP" + num2str(i)
		Duplicate/O Etot_sMagSqr $EsqrdWaveName
		EsqrdSepWaveRefs[i][1] = $EsqrdWaveName
		
		// plus ifm molecules
		sLoc*=0.5
		MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
		MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
		MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
		MatrixOP/O Etot_sMagSqr = magsqr(Etot)
		EsqrdWaveName = "EtotMagSqrPmol" + num2str(i)
		Duplicate/O Etot_sMagSqr $EsqrdWaveName
		EsqrdSepWaveRefs[i][2] = $EsqrdWaveName
		
		// minus ifm molecules
		sLoc*=-1
		MatrixOP/O DplusYY = d + yy - (sLoc * vrecip)
		MatrixOP/O DminusYY =  sLoc*vrecip - yy + d 
		MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
		MatrixOP/O Etot_sMagSqr = magsqr(Etot)
		EsqrdWaveName = "EtotMagSqrMmol" + num2str(i)
		Duplicate/O Etot_sMagSqr $EsqrdWaveName
		EsqrdSepWaveRefs[i][3] = $EsqrdWaveName
	
		i+=1
	while(i<numpnts(SeparationWave))
	
	KillWaves/Z Etot_sMagSqr
End





Function CalcPhiOverPol2IFMmol_iso()
	Wave EtotMagSqr0, vrecip//, EtotMagSqrM, EtotMagSqrP, IFMphase, EtotMagSqrMmol, EtotMagSqrPmol
	NVAR num_pos
	// Find the phase shift for each velocity component of the main and diffracted beam
	MatrixOP/O dphi = EtotMagSqr0 * vrecip							//phase shift of first path:   alpha*E^2/v	// this is the +1 diffraction order
	Integrate/dim=1 dphi 
	
	Wave isotopeFracs
	Wave/Wave EsqrdSepWaveRefs
	String IFMphase_v_overpolWaveName
	Make/O/WAVE/N=(numpnts(isotopeFracs),4) IFMphase_v_overpol_Refs
	Variable i = 0, j=0
	do
		// For the diffracted beams now:
		Wave LocalEsqrd =  EsqrdSepWaveRefs[i][0]
		IFMphase_v_overpolWaveName = "IFMphase_v_overpolM"+num2str(i)
		MatrixOP/O dphi_s = LocalEsqrd* vrecip
		Integrate/dim=1 dphi_s
		Make/D/O/N=(num_pos,number_v) $IFMphase_v_overpolWaveName; Wave TempIFMphase_v_overpol = $IFMphase_v_overpolWaveName
		TempIFMphase_v_overpol = dphi[p][100][q] - dphi_s[p][100][q]
		IFMphase_v_overpol_Refs[i][0] = $IFMphase_v_overpolWaveName
		
		Wave LocalEsqrd =  EsqrdSepWaveRefs[i][1]
		IFMphase_v_overpolWaveName = "IFMphase_v_overpolP"+num2str(i)
		MatrixOP/O dphi_s = LocalEsqrd* vrecip
		Integrate/dim=1 dphi_s
		Make/D/O/N=(num_pos,number_v) $IFMphase_v_overpolWaveName; Wave TempIFMphase_v_overpol = $IFMphase_v_overpolWaveName
		TempIFMphase_v_overpol = dphi_s[p][100][q] - dphi[p][100][q]
		IFMphase_v_overpol_Refs[i][1] = $IFMphase_v_overpolWaveName

		Wave LocalEsqrd =  EsqrdSepWaveRefs[i][2]
		IFMphase_v_overpolWaveName = "IFMphase_v_overpolMmol"+num2str(i)
		MatrixOP/O dphi_s = LocalEsqrd* vrecip
		Integrate/dim=1 dphi_s
		Make/D/O/N=(num_pos,number_v) $IFMphase_v_overpolWaveName; Wave TempIFMphase_v_overpol = $IFMphase_v_overpolWaveName
		TempIFMphase_v_overpol = dphi[p][100][q] - dphi_s[p][100][q]
		IFMphase_v_overpol_Refs[i][2] = $IFMphase_v_overpolWaveName
		
		Wave LocalEsqrd =  EsqrdSepWaveRefs[i][3]
		IFMphase_v_overpolWaveName = "IFMphase_v_overpolPmol"+num2str(i)
		MatrixOP/O dphi_s = LocalEsqrd* vrecip
		Integrate/dim=1 dphi_s
		Make/D/O/N=(num_pos,number_v) $IFMphase_v_overpolWaveName; Wave TempIFMphase_v_overpol = $IFMphase_v_overpolWaveName
		TempIFMphase_v_overpol = dphi_s[p][100][q] - dphi[p][100][q]
		IFMphase_v_overpol_Refs[i][3] = $IFMphase_v_overpolWaveName
		
		i+=1
	while(i<numpnts(SeparationWave))	
	
	
//	MatrixOP/O dphi_sM = EtotMagSqrM * vrecip					//phase shift of second path a distance "sep" away		// this is the 0th order
//	MatrixOP/O dphi_sP = EtotMagSqrP * vrecip				//phase shift of third path a distance "sep" away		// this is the -1 order
//	MatrixOP/O dphiMmol = EtotMagSqrMmol * vrecip
//	MatrixOP/O dphiPmol = EtotMagSqrPmol * vrecip
//	
//	Integrate/dim=1 dphi, dphi_sM, dphi_sP, dphiMmol, dphiPmol		// integrate along beam path to find total acquired phase
//	
//
//	Make/D/O/N=(num_pos, number_v) IFMphase_v_overpolM, IFMphase_v_overpolP, IFMphase_v_overpolMmol, IFMphase_v_overpolPmol
//	
//	// these waves are proportional to phi(x[row], v[column])
//	// later we'll take the real and imaginary parts of them and integrate across columns
//	IFMphase_v_overpolM = dphi[p][100][q] - dphi_sM[p][100][q]
//	IFMphase_v_overpolP = dphi_sP[p][100][q] - dphi[p][100][q]
//	IFMphase_v_overpolMmol = dphi[p][100][q] - dphiMmol[p][100][q]
//	IFMphase_v_overpolPmol = dphiPmol[p][100][q] - dphi[p][100][q]
	
//	MatrixOP/O IFMphase_v_overpolM = dphi - dphi_sM				// phase shift divided by polarizability for minus ifm
//	MatrixOP/O IFMphase_v_overpolP = dphi_sP - dphi			// phase shift divided by polarizability for plus ifm
//	MatrixOP/O IFMphase_v_overpolMmol = dphi - dphiMmol				// phase shift divided by polarizability for minus ifm
//	MatrixOP/O IFMphase_v_overpolPmol = dphiPmol - dphi			// phase shift divided by polarizability for plus ifm
	
End



Function Phi2IFMmol_iso(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam

	//DFREF saveDFR = GetDataFolderDFR(); SetDataFolder root:'090625g'	// Uncomment when performing a GlobalFit with ContrastVratio()
	
	Wave IFMprob, density, density2D, Contrast //,IFMphase_v_overpolM, IFMphase_v_overpolP, IFMphase_v_overpolMmol, IFMphase_v_overpolPmol
	NVAR molFrac, molPol
	
	Variable molFracLoc = molFrac/100
	
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	Variable polmol = molPol*4*pi*eps0*1E-30* x_step / (2 * hbar)
	//print molpol
	Wave phiSagGravLayered2D

	NVAR num_pos
	Make/O/D/N=(num_pos, number_v) imagpart=0, realpart=0
	Wave/Wave IFMphase_v_overpol_Refs
	Wave IsotopeFracs
	Variable isotopeFrac
	Variable i = 0, j=0
	do
		isotopeFrac = IsotopeFracs[i]
		
		Wave IFMphase_v_overpolTemp=IFMphase_v_overpol_Refs[i][0]; MatrixOP/O IFMphase_vM=pol*IFMphase_v_overpolTemp+phiSagGravLayered2D
		Wave IFMphase_v_overpolTemp=IFMphase_v_overpol_Refs[i][1]; MatrixOP/O IFMphase_vP=pol*IFMphase_v_overpolTemp+phiSagGravLayered2D
		Wave IFMphase_v_overpolTemp=IFMphase_v_overpol_Refs[i][2]; MatrixOP/O IFMphase_vMmol=molpol*IFMphase_v_overpolTemp+phiSagGravLayered2D
		Wave IFMphase_v_overpolTemp=IFMphase_v_overpol_Refs[i][3]; MatrixOP/O IFMphase_vPmol=molpol*IFMphase_v_overpolTemp+phiSagGravLayered2D
		
		MatrixOP/O imagpart = imagpart + isotopeFrac*((1-molFracLoc)*(sin(IFMphase_vM) + sin(IFMphase_vP))    +    molFracLoc*(sin(IFMphase_vMmol) + sin(IFMphase_vPmol))	 )
		MatrixOP/O realpart =  realpart  + isotopeFrac*((1-molFracLoc)*(cos(IFMphase_vM) + cos(IFMphase_vP))    +    molFracLoc*(cos(IFMphase_vMmol) + cos(IFMphase_vPmol)) 	 )
		
		i+=1
	while(i<numpnts(IsotopeFracs))	

	//0.5 to account for two IFMs										
	MatrixOP/O imagpart = 0.5  * density2D * imagpart			// weight imaginary part by prob(v) 
	MatrixOP/O realpart =  0.5  * density2D * realpart			// weight real part by prob(v)
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



