#pragma rtGlobals=1		// Use modern global access method.






Function PredictPhase2IFM(alpha, offset)
	Variable alpha
	Variable offset
	
	NVAR velocity, sigma
	
	Variable timerRefNum = Startmstimer
	
	Make/O/N=1 paramWave = {alpha}
	
	Make/O/N=20 position = .1 * x/1000+.0001, predictedPhase
	//Make/O/N=10 position = .1 * x/1000+.0011, predictedPhase
	
	Variable/G num_pos = numpnts(position)
	Make/o/n=(num_pos,num_x, number_v) yy xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam
	Make/o/n=(num_pos) IFMphase2 IFMprob2 Contrast2

	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = offset + position					// account for y offset
	Calc_Efield2IFM(position_shifted)							// calculate E-field at each data point
	
	CalcPhiOverPol2IFM()										// calculate phase/polarizability
	
	Phi2IFMv2(paramWave, predictedPhase, position)
	
	//Print "New: ", StopMSTimer(timerRefNum)*10^-6
End




// the positions in 2IFMv3 represent the position of the 0th order beam
Function PredictPhase2IFMv3(alpha, offset)
	Variable alpha
	Variable offset
	
	Variable timerRefNum = Startmstimer
	
	Make/O/N=1 paramWave = {alpha}
	
	Make/O/N=20 position = .1 * x/1000, predictedPhase		//position of 0th order
	//Make/O/N=10 position = .1 * x/1000+.0011, predictedPhase
	
	Variable/G num_pos = numpnts(position)
	Make/o/n=(num_pos,num_x, number_v) yy xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam
	Make/o/n=(num_pos) IFMphase2 IFMprob2 Contrast2

	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = offset + position					// account for y offset
	Calc_Efield2IFMv3(position_shifted)							// calculate E-field at each data point
	
	CalcPhiOverPol2IFMv3()										// calculate phase/polarizability
	
	Phi2IFMv3(paramWave, predictedPhase, position)
	
	//Print "New: ", StopMSTimer(timerRefNum)*10^-6
End



Function WideBeam2IFM()
	//Wave predictedPhase
	
	Variable halfwidth=25	//microns
	String phasewavename
	Variable i
	for(i=-halfwidth;i<=halfwidth;i+=5)	// Initialize variables;continue test
		PredictPhase2IFMv3(43, i*1E-6)	
		phasewavename= "predictedPhase"+Num2Str(i)				
		Duplicate/O predictedPhase $phasewavename
		print phasewavename
	endfor										
	
	PredictPhase2IFMv3(43, 0)//just an inefficient way of resetting the x axis to be centered on the beam again.
End






Function PhaseFit2IFM(PositionWave, PhaseWave, ErrorWave)
	Wave PositionWave, PhaseWave, ErrorWave
	
	Variable timerRefNum = Startmstimer
	
	Variable/G num_pos = numpnts(PositionWave)
	Make/o/n=(num_pos,num_x, number_v) yy xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)						// The x distance squared
	Make/o/n=(num_pos) IFMphase IFMprob Contrast		// Waves for: total phase shift of monochrom. beam; phase shift of polychromatic beam
	Make/o/n=(num_pos) IFMphase2 IFMprob2 Contrast2
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = PositionWave + shift		// account for y offset
	Calc_Efield2IFM(position_shifted)							// calculate E-field at each data point

	CalcPhiOverPol2IFM()										// calculate phase/polarizability
	
	Make/D/O W_coef={43} 				// Initial Guess
	Make/D/O Epsilonwave={1e-1}			// Amount to vary the fit parameter by in each iteration	
	Make/T/O PolConstraint = {"K0>40","K0<46"}
//	FuncFit/N/M=2/TBOX=768 Phi2IFMv2 W_coef PhaseWave[20,*]  /X=PositionWave[20,*] /W=ErrorWave[20,*] /I=1 /R /E=EpsilonWave ///M=mask_b
	FuncFit/N/M=2/TBOX=768 Phi2IFMv2 W_coef PhaseWave  /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave ///M=mask_b /C=PolConstraint

	AppendFittedPolWave2IFM(PositionWave, PhaseWave)				

	//Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
End





Function PhaseFit2IFMv3(PositionWave, PhaseWave, ErrorWave)
	Wave PositionWave, PhaseWave, ErrorWave
	
	Variable timerRefNum = Startmstimer
	
	Variable/G num_pos = numpnts(PositionWave)
	Make/o/n=(num_pos,num_x, number_v) yy xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)						// The x distance squared
	Make/o/n=(num_pos) IFMphase IFMprob Contrast		// Waves for: total phase shift of monochrom. beam; phase shift of polychromatic beam
	Make/o/n=(num_pos) IFMphase2 IFMprob2 Contrast2
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = PositionWave + shift		// account for y offset
	Calc_Efield2IFMv3(position_shifted)							// calculate E-field at each data point

	CalcPhiOverPol2IFMv3()										// calculate phase/polarizability
	
	Make/D/O W_coef={43} 				// Initial Guess
	Make/D/O Epsilonwave={1e-1}			// Amount to vary the fit parameter by in each iteration	
	Make/T/O PolConstraint = {"K0>40","K0<46"}
//	FuncFit/N/M=2/TBOX=768 Phi2IFMv2 W_coef PhaseWave[20,*]  /X=PositionWave[20,*] /W=ErrorWave[20,*] /I=1 /R /E=EpsilonWave ///M=mask_b
	FuncFit/N/M=2/TBOX=768 Phi2IFMv3 W_coef PhaseWave  /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave ///M=mask_b /C=PolConstraint

	AppendFittedPolWave2IFMv3(PositionWave, PhaseWave)				

	//Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
End






Function PhaseFitOffset2IFM(PositionWave, PhaseWave, ErrorWave)
	Wave PositionWave, PhaseWave, ErrorWave
	
	Variable timerRefNum = Startmstimer
	
	Variable/G num_pos = numpnts(PositionWave)
	Make/o/n=(num_pos,num_x, number_v) yy xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam
	Make/o/n=(num_pos) IFMphase2 IFMprob2 Contrast2
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Make/D/O W_coef={43, 48.7e-6}	// initial guesses {alpha, shift}
	Make/D/O Epsilonwave={1e-1,1e-6}	// Amount to vary the fit parameters by in each iteration
	FuncFit/N/M=2 PhiOffset2IFM W_coef PhaseWave /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave 
	
	AppendFittedPolWave2IFM(PositionWave, PhaseWave)	
	
	//Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
End



Function PhaseFitOffset2IFMv3(PositionWave, PhaseWave, ErrorWave)
	Wave PositionWave, PhaseWave, ErrorWave
	
	Variable timerRefNum = Startmstimer
	
	Variable/G num_pos = numpnts(PositionWave)
	Make/o/n=(num_pos,num_x, number_v) yy xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam
	Make/o/n=(num_pos) IFMphase2 IFMprob2 Contrast2
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Make/D/O W_coef={43, 48.7e-6}	// initial guesses {alpha, shift}
	Make/D/O Epsilonwave={1e-1,1e-6}	// Amount to vary the fit parameters by in each iteration
	FuncFit/N/M=2 PhiOffset W_coef PhaseWave /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave 
	
	AppendFittedPolWave2IFMv3(PositionWave, PhaseWave)	
	
	//Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
End








Function Calc_EField2IFM(positionwave)			// Calculates the electric field of the gradient region and makes waves for phase, prob, and contrast
	Wave positionwave						// positions at which to calculate the E-field
	Nvar num_pos, Vwire, d, lambda, E0, s
	
	Wave yy, xx, xxsqr, vrecip
	yy[][][] = positionwave[p]				// populate the yy wave with the positions. Each row contains a different position[p], but the position is the same across columns and layers.
	
	// some auxiliary waves that make the E field calculations more efficient
	MatrixOP/O DplusYY = d+yy			
	MatrixOP/O DminusYY = -1*yy + d
	MatrixOP/O DplusYsep = d + yy - (s * vrecip)
	MatrixOP/O DminusYsep =  s*vrecip - yy + d 
	MatrixOP/O DplusYsep2 =  d + yy - 2*(s * vrecip)	//for IFM2
	MatrixOP/O DminusYsep2 = 2*s*vrecip - yy + d 
	
	// The electric field magnitude of the primary and diffracted beams
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O Etot_s= (E0)*sqrt(	magsqr(DplusYsep/(magsqr(DplusYsep)+xxsqr) + DminusYSep/(magsqr(DminusYSep)+xxsqr)) + magsqr(xx/(magsqr(DplusYsep)+xxsqr) - xx/(magsqr(DminusYsep)+xxsqr))		)
	MatrixOP/O Etot_s2 = (E0)*sqrt(	magsqr(DplusYsep2/(magsqr(DplusYsep2)+xxsqr) + DminusYSep2/(magsqr(DminusYSep2)+xxsqr)) + magsqr(xx/(magsqr(DplusYsep2)+xxsqr) - xx/(magsqr(DminusYsep2)+xxsqr))		)
	
	
	// Need the magnitude squared.
	MatrixOP/O EtotMagSqr = magsqr(Etot)
	MatrixOP/O Etot_sMagSqr = magsqr(Etot_s)
	MatrixOP/O Etot_s2MagSqr = magsqr(Etot_s2)
End





Function Calc_EField2IFMv3(positionwave)			// Calculates the electric field of the gradient region and makes waves for phase, prob, and contrast
	Wave positionwave						// positions at which to calculate the E-field
	Nvar num_pos, Vwire, d, lambda, E0, s
	
	Wave yy, xx, xxsqr, vrecip
	yy[][][] = positionwave[p]				// populate the yy wave with the positions. Each row contains a different position[p], but the position is the same across columns and layers.
	
	//	yy is the position of the 0th order beam in the 2IFMv3 program!!!
	
	// some auxiliary waves that make the E field calculations more efficient
	MatrixOP/O DplusYYplus = d+yy	 + (s * vrecip)					//
	MatrixOP/O DminusYYplus = -1*s*vrecip - yy + d					//		the position of the +1 beam with respect to the ideal line charge
	MatrixOP/O DplusYYminus = d + yy - (s * vrecip)		// 
	MatrixOP/O DminusYYminus =  s*vrecip - yy + d 		//		
	
	// The electric field magnitude of the primary and diffracted beams
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYYplus/(magsqr(DplusYYplus)+xxsqr) + DminusYYplus/(magsqr(DminusYYplus)+xxsqr)) + magsqr(xx/(magsqr(DplusYYplus)+xxsqr) - xx/(magsqr(DminusYYplus)+xxsqr))   		)
	MatrixOP/O Etot_s= (E0)*sqrt(	magsqr(DplusYYminus/(magsqr(DplusYYminus)+xxsqr) + DminusYYminus/(magsqr(DminusYYminus)+xxsqr)) + magsqr(xx/(magsqr(DplusYYminus)+xxsqr) - xx/(magsqr(DminusYYminus)+xxsqr))		)

	// Need the magnitude squared.
	MatrixOP/O EtotMagSqr = magsqr(Etot)
	MatrixOP/O Etot_sMagSqr = magsqr(Etot_s)
End






Function CalcPhiOverPol2IFM()
	Wave EtotMagSqr, Etot_sMagSqr, Etot_s2MagSqr, vrecip, IFMphase 
	// Find the phase shift for each velocity component of the main and diffracted beam
	MatrixOP/O dphi = EtotMagSqr * vrecip							//phase shift of first path:   alpha*E^2/v	// this is the +1 diffraction order
	MatrixOP/O dphi_s = Etot_sMagSqr * vrecip					//phase shift of second path a distance "sep" away		// this is the 0th order
	MatrixOP/O dphi_s2 = Etot_s2MagSqr * vrecip				//phase shift of third path a distance "sep" away		// this is the -1 order
	
	Integrate/dim=1 dphi, dphi_s, dphi_s2
	MatrixOP/O IFMphase_v_overpol = dphi - dphi_s				// phase shift divided by polarizability
	MatrixOP/O IFMphase_v_overpol2 = dphi_s - dphi_s2			// phase shift divided by polarizability for path 2
		
	// Just for curiousity (and debugging): the phase shift of the average velocity path:
	//IFMphase = IFMphase_v_overpol[p][num_x-1][.5*number_v]		
End



Function CalcPhiOverPol2IFMv3()
	Wave EtotMagSqr, Etot_sMagSqr, vrecip, IFMphase 
	// Find the phase shift for each velocity component of the main and diffracted beam
	MatrixOP/O dphi = EtotMagSqr * vrecip							//phase shift of first path:   alpha*E^2/v	
	MatrixOP/O dphi_s = Etot_sMagSqr * vrecip					//phase shift of second path a distance "sep" away
	Integrate/dim=1 dphi, dphi_s 
	MatrixOP/O IFMphase_v_overpol = (dphi - dphi_s)/2				// phase shift divided by polarizability
		
	// Just for curiousity (and debugging): the phase shift of the average velocity path:
	//IFMphase = IFMphase_v_overpol[p][num_x-1][.5*number_v]		
End











Function Phi2IFM(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
	
	Wave Etot, Etot_s, Etot_s2, IFMprob, IFMprob2, density, Contrast, Contrast2, IFMphase_v_overpol, IFMphase_v_overpol2
	
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	
	MatrixOP/O IFMphase_v = pol * IFMphase_v_overpol			// now we have the actual phase shift to take the real and imaginary parts of
	MatrixOP/O IFMphase_v2 = pol * IFMphase_v_overpol2
	
	// We need to find the total velocity weighted phase shift
	// phi_total = ArcTan( Int( Prob(v) * exp(i * delta_phi) ) ) 
	//			 = ArcTan( Int(Prob(v)*sin(phi)) / Int(Prob(v)*cos(phi)) )
	// the code to implement this follows:
	
	Variable loc_v_step = v_step											
	MatrixOP/O imagpart = loc_v_step * density * sin(IFMphase_v)			// weight imaginary part by prob(v) 
	MatrixOP/O realpart = loc_v_step * density * cos(IFMphase_v)				// weight real part by prob(v)
	Integrate/DIM=2 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	MatrixOP/O imagpart2 = loc_v_step * density * sin(IFMphase_v2)			// weight imaginary part by prob(v) 
	MatrixOP/O realpart2 = loc_v_step * density * cos(IFMphase_v2)				// weight real part by prob(v)
	Integrate/DIM=2 imagpart2, realpart2			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][num_x-1][number_v-1]) + magsqr(realpart[p][num_x-1][number_v-1]) )
	IFMprob = atan(imagpart[p][num_x-1][number_v-1]/realpart[p][num_x-1][number_v-1])
	//MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	Contrast2 = sqrt( magsqr(imagpart2[p][num_x-1][number_v-1]) + magsqr(realpart2[p][num_x-1][number_v-1]) )
	IFMprob2 = atan(imagpart2[p][num_x-1][number_v-1]/realpart2[p][num_x-1][number_v-1])
	MatrixOP/O PhaseWrapped2 = IFMprob2									// Duplicates IFMprob before unwrap for debugging

	Unwrap Pi, IFMprob		//Correction from taking the arctan.  
	Unwrap Pi, IFMprob2
	
	if(IFMprob[0]<-0.5)
		IFMprob+=pi
	endif
	
	if(IFMprob2[0]<-0.5)
		IFMprob2+=pi
	endif
	
	
	output = (IFMprob+IFMprob2)/2
End




Function Phi2IFMv2(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
	
	Wave Etot, Etot_s, Etot_s2, IFMprob, IFMprob2, density, Contrast, Contrast2, IFMphase_v_overpol, IFMphase_v_overpol2
	
	Variable polcgs = pw[0]
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	
	MatrixOP/O IFMphase_v = pol * IFMphase_v_overpol			// now we have the actual phase shift to take the real and imaginary parts of
	MatrixOP/O IFMphase_v2 = pol * IFMphase_v_overpol2
	
	// We need to find the total velocity weighted phase shift
	// phi_total = ArcTan( Int( Prob(v) * exp(i * delta_phi) ) ) 
	//			 = ArcTan( Int(Prob(v)*sin(phi)) / Int(Prob(v)*cos(phi)) )
	// the code to implement this follows:
	
	Variable loc_v_step = v_step											
	MatrixOP/O imagpart = loc_v_step * density * sin(IFMphase_v)			// weight imaginary part by prob(v) 
	MatrixOP/O realpart = loc_v_step * density * cos(IFMphase_v)				// weight real part by prob(v)
	Integrate/DIM=2 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	MatrixOP/O imagpart2 = loc_v_step * density * sin(IFMphase_v2)			// weight imaginary part by prob(v) 
	MatrixOP/O realpart2 = loc_v_step * density * cos(IFMphase_v2)				// weight real part by prob(v)
	Integrate/DIM=2 imagpart2, realpart2			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][num_x-1][number_v-1]) + magsqr(realpart[p][num_x-1][number_v-1]) )
	IFMprob = atan2(imagpart[p][num_x-1][number_v-1],realpart[p][num_x-1][number_v-1])
	MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	Contrast2 = sqrt( magsqr(imagpart2[p][num_x-1][number_v-1]) + magsqr(realpart2[p][num_x-1][number_v-1]) )
	IFMprob2 = atan2(imagpart2[p][num_x-1][number_v-1],realpart2[p][num_x-1][number_v-1])
	MatrixOP/O PhaseWrapped2 = IFMprob2									// Duplicates IFMprob before unwrap for debugging


	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
	Unwrap 2*Pi, IFMprob2
	
	if(IFMprob[0] < 0)
		IFMprob+=2*pi
	endif
	
	if(IFMprob2[0] < 0)
		IFMprob2+=2*pi
	endif
	
	//IFMprob+=2*pi
	//IFMprob2+=2*pi
	
	output = (IFMprob+IFMprob2)/2
End



Function Phi2IFMv3(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
	
	Wave Etot, Etot_s, IFMprob, density, Contrast, IFMphase_v_overpol
	
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	
	MatrixOP/O IFMphase_v = pol * IFMphase_v_overpol			// now we have the actual phase shift to take the real and imaginary parts of
	
	// We need to find the total velocity weighted phase shift
	// phi_total = ArcTan( Int( Prob(v) * exp(i * delta_phi) ) ) 
	//			 = ArcTan( Int(Prob(v)*sin(phi)) / Int(Prob(v)*cos(phi)) )
	// the code to implement this follows:
	
	Variable loc_v_step = v_step											
	MatrixOP/O imagpart = loc_v_step * density * sin(IFMphase_v)			// weight imaginary part by prob(v) 
	MatrixOP/O realpart = loc_v_step * density * cos(IFMphase_v)				// weight real part by prob(v)
	Integrate/DIM=2 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][num_x-1][number_v-1]) + magsqr(realpart[p][num_x-1][number_v-1]) )
	IFMprob = atan2(imagpart[p][num_x-1][number_v-1],realpart[p][num_x-1][number_v-1]) 
	MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
	if(IFMprob[0]<-.2)
		IFMprob+=2*pi
	endif
	
	//IFMprob/=2
	
	output = IFMprob
	
// the following is wrong because we would integrate just the phase itself over velocity (rather than exp(i phi)):
//	MatrixOP/O ProbWeightedPhase = loc_v_step * density * IFMphase_v	// weight each calculation by the appropriate amount
//	Integrate/DIM=2 ProbWeightedPhase									// must integrate along dv to find total phase shift
//	output = ProbWeightedPhase[p][num_x-1][number_v-1]				// select the column corresponding to the total integral over x and v
End






Function PhiOffset2IFM(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
	
	Wave Etot, Etot_s, Etot_s2, IFMprob, IFMprob2, density, Contrast, Contrast2, IFMphase_v_overpol, IFMphase_v_overpol2	
	
	MatrixOP/O position_shifted = pw[1] + y				// pw[1] is the y offset fit parameter
	Calc_Efield2IFM(position_shifted)							// calculate E-field at each data point
	CalcPhiOverPol2IFM()										// calculate phase/polarizability
	
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	
	MatrixOP/O IFMphase_v = pol * IFMphase_v_overpol
	MatrixOP/O IFMphase_v2 = pol * IFMphase_v_overpol2
	
	// We need to find the total velocity weighted phase shift
	// phi_total = ArcTan( Int( Prob(v) * exp(i * delta_phi) ) ) 
	//			 = ArcTan( Int(Prob(v)*sin(phi)) / Int(Prob(v)*cos(phi)) )
	// the code to implement this follows:
	
	Variable loc_v_step = v_step											
	MatrixOP/O imagpart = loc_v_step * density * sin(IFMphase_v)			// weight imaginary part by prob(v) 
	MatrixOP/O realpart = loc_v_step * density * cos(IFMphase_v)				// weight real part by prob(v)
	Integrate/DIM=2 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	MatrixOP/O imagpart2 = loc_v_step * density * sin(IFMphase_v2)			// weight imaginary part by prob(v) 
	MatrixOP/O realpart2 = loc_v_step * density * cos(IFMphase_v2)				// weight real part by prob(v)
	Integrate/DIM=2 imagpart2, realpart2			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][num_x-1][number_v-1]) + magsqr(realpart[p][num_x-1][number_v-1]) )
	IFMprob = atan2(imagpart[p][num_x-1][number_v-1],realpart[p][num_x-1][number_v-1])
	//MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	Contrast2 = sqrt( magsqr(imagpart2[p][num_x-1][number_v-1]) + magsqr(realpart2[p][num_x-1][number_v-1]) )
	IFMprob2 = atan2(imagpart2[p][num_x-1][number_v-1],realpart2[p][num_x-1][number_v-1])
	//MatrixOP/O PhaseWrapped2 = IFMprob2									// Duplicates IFMprob before unwrap for debugging
	
	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
	if(IFMprob[0]<0)
		IFMprob+=2*pi
	endif
	
	Unwrap 2*Pi, IFMprob2
	if(IFMprob2[0]<0)
		IFMprob2+=2*pi
	endif
	
	output = (IFMprob+IFMprob2)/2
End







Function AppendFittedPolWave2IFM(PositionWave, PhaseWave)
	Wave PositionWave, PhaseWave
	Wave W_coef, position, predictedPhase
	
	If(numpnts(W_coef)==1)						// Determine if the offset was fixed or a fit parameter
		PredictPhase2IFM(W_coef[0], offset)			// Run the predict phase routine with the fitted polarizability and fixed offset
	else
		PredictPhase2IFM(W_coef[0],W_coef[1])		// Run the predict phase routine with the fitted polarizaiblity and offset
	EndIf
	
	String fitPosWaveName= NameOfWave(positionWave)+"_fit"			// make some descriptive wave names
	String fitPhaseWaveName = NameOfWave(phaseWave) +"_fit"
	
	Duplicate/O position $fitPosWaveName; Wave fitPosWave = $fitPosWaveName					// assign the calculated position and phase waves to the right wave names
	Duplicate/O predictedPhase $fitPhaseWaveName; Wave fitPhaseWave = $fitPhaseWaveName
	
	If(StringMatch(TraceInfo("", fitPhaseWaveName,0),""))		// if the trace is not already on the graph then put it there
		AppendToGraph fitPhaseWave vs fitPosWave
		ModifyGraph rgb($fitPhaseWaveName)=(1,16019,65535)
		ModifyGraph zero(Res_Left)=1
		SetAxis left 0,*
		SetAxis/A Res_Left
	EndIf
End


Function AppendFittedPolWave2IFMv3(PositionWave, PhaseWave)
	Wave PositionWave, PhaseWave
	Wave W_coef, position, predictedPhase
	
	If(numpnts(W_coef)==1)						// Determine if the offset was fixed or a fit parameter
		PredictPhase2IFMv3(W_coef[0], offset)			// Run the predict phase routine with the fitted polarizability and fixed offset
	else
		PredictPhase2IFMv3(W_coef[0],W_coef[1])		// Run the predict phase routine with the fitted polarizaiblity and offset
	EndIf
	
	String fitPosWaveName= NameOfWave(positionWave)+"_fit"			// make some descriptive wave names
	String fitPhaseWaveName = NameOfWave(phaseWave) +"_fit"
	
	Duplicate/O position $fitPosWaveName; Wave fitPosWave = $fitPosWaveName					// assign the calculated position and phase waves to the right wave names
	Duplicate/O predictedPhase $fitPhaseWaveName; Wave fitPhaseWave = $fitPhaseWaveName
	
	If(StringMatch(TraceInfo("", fitPhaseWaveName,0),""))		// if the trace is not already on the graph then put it there
		AppendToGraph fitPhaseWave vs fitPosWave
		ModifyGraph rgb($fitPhaseWaveName)=(1,16019,65535)
		ModifyGraph zero(Res_Left)=1
		SetAxis left 0,*
		SetAxis/A Res_Left
	EndIf
End