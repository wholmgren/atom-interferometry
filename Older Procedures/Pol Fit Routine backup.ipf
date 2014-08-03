-#pragma rtGlobals=1		// Use modern global access method.

#include <Remove Points>



//////////////////////////////////////////////////////////////////////////////////////////////////
//																												     //
//																												     //
//				POLARIZABILITY FITTING																		     //
//																												     //
//																												     //
//////////////////////////////////////////////////////////////////////////////////////////////////




// There are three primary programs in this procedure:
// 	1. PredictPhase(alpha, offset): Predicts the phase shift given a polarizability and offset from the ground plane
//	2. PhaseFit(PositionWave, PhaseWave, ErrorWave): Finds the best fit polarizability from a reduced data set. Requires the ground plane offset to be found otherwise and set manually
//	3. PhaseFitOffset(PositionWave, PhaseWave, ErrorWave): Finds the best fit polarizability and ground plane offset from a reduced data set.

// Before running the fit procedures you should reduce your data sets for efficiency. 
// The function ReducePolWaves(position,phase,error) will take a long set of positions, phases, and error waves and pull out just the
// points that are not NaNs. Then it will reorder the points so that point 0 is closest to the ground plane and the last point is closest to the electrode.
// The input waves will not be touched. Instead, the function will create new waves with the same name plus "_red"(uced) appended.

// The beam velocity and standard deviation parameters are set in the constants section.
// The offset parameter for fixed offset fitting is also set in the constants section.

// 	The assumed geometry is described below.
// "x" is parallel to the beam
// "y" is the distance from the ground plane to the diffracted beam. "0" when the diffracted beam is centered on the ground plane and the 0th order is trying to go through it


//							|			      0   +1
//							|			     ^^^  ^^^
//							|				|	|
//							|				|	|	
//							|				|	|	  (	    )
//			Ground			|				|	|	(		) HV electrode
//			Plane			|				|	|	  (	    )
//							|				|	|		
//							|				|	|
//							|
//							|------------->|
//							|		   y



// Fundamental constants
Constant eps0 = 8.8541878e-12
Constant hbar = 1.054572e-34
//Constant m=   3.81753e-26				//Na = 22.9897 * 1e-3 / 6.02214e23 	//mass of element in kg: 85.4678  39.0983 22.9897 	
Constant m =     6.49243e-26			//K natural abundance
//Constant m = 	  6.4761e-26				//K39 only 

// Geometry constants
Constant gap = 1.998e-3					//Gap between the ground and cylindrical electrodes [m]
Constant R =  .0063315						//Radius of the electrode [m]
Constant L = .8044 							//.7813	//Distance from 1g to the grad E region
Constant a = 1e-7							//Grating period

// x-axis integration constants
Constant num_x=70							// position integration points along the beam path
Constant x_step =   0.000999				// = (gap*35)/num_x		//35 mm from gap to edge of ground plane

// velocity integration constants
Constant number_v = 50					// number of individual velocities for integration
Constant v_step = 30						// separation of the individual velocities in m/s

Constant n=3								//n= 3, Maxwellian;    0, gaussian

// velocity parameters
Constant velocity = 1452.5		//g
//Constant velocity = 1450.4		//h
//Constant velocity = 1449.95		//m
//Constant velocity = 1990.84		// b 4/7/09

//Constant sigma = 112.8 		 	//g
Constant sigma = 112.5 		 	//h
//Constant sigma = 112.4 		 	//m
//Constant sigma = 154.4				//b 4/7/09

	
// y offset parameter
Constant offset = 74.1e-6 	//g
//Constant offset = 49.6e-6 	//h
//Constant offset =   54.6e-6 	//h 2IFM
//Constant offset =   43.3e-6 	//h 2IFM counts=70%
//Constant offset = 28.8e-6		//m counts contrast = 50%
//Constant offset = 2.13e-5		//m counts=50%
//Constant offset =   5.35e-05		//m counts = 80%
//Constant offset = -35.5e-6	//q 10/16/08
//Constant offset = 4.39542e-05 // b 4/7/09

Constant HVmon = 0.9165	//g
//Constant HVmon = 1.143		//h
//Constant HVmon = 0.961 		//m
//Constant HVmon = 1.010	//q 10/16/08
//Constant HVmon= 1.57		// b 4/7/09


Function PredictPhase(alpha, offset)
	Variable alpha
	Variable offset
	
	Variable timerRefNum = Startmstimer
	
	Make/O/N=1 paramWave = {alpha}
	
	Make/O/N=20 position = .1 * x/1000+.0001, predictedPhase
	
	Variable/G num_pos = numpnts(position)
	Make/o/n=(num_pos,num_x, number_v) yy xx   						// yy: along a column we increase in y position with each row. Along a row we step along x in each column.
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam

	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = offset + position					// account for y offset
	Calc_Efield(position_shifted)							// calculate E-field at each data point
	
	CalcPhiOverPol()										// calculate phase/polarizability
	
	Phi(paramWave, predictedPhase, position)
	
	//Print "New: ", StopMSTimer(timerRefNum)*10^-6
End





Function PredictPhase2IFM(alpha, offset)
	Variable alpha
	Variable offset
	
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




//PhaseFit(pos_b_polfit_pad_red, phase_b_polfit_red, phase_error_b_polfit_red)
Function PhaseFit(PositionWave, PhaseWave, ErrorWave)
	Wave PositionWave, PhaseWave, ErrorWave
	
	Variable timerRefNum = Startmstimer
	
	Variable/G num_pos = numpnts(PositionWave)
	Make/o/n=(num_pos,num_x, number_v) yy xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)						// The x distance squared
	Make/o/n=(num_pos) IFMphase IFMprob Contrast		// Waves for: total phase shift of monochrom. beam; phase shift of polychromatic beam
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = PositionWave + shift		// account for y offset
	Calc_Efield(position_shifted)							// calculate E-field at each data point

	CalcPhiOverPol()										// calculate phase/polarizability
	
	Make/D/O W_coef={42.5} 				// Initial Guess
	Make/D/O Epsilonwave={1e-1}			// Amount to vary the fit parameter by in each iteration	
	FuncFit/N/M=2/TBOX=768 Phi W_coef PhaseWave  /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave ///M=mask_b
	
	AppendFittedPolWave(PositionWave, PhaseWave)				

	//Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
End






Function PhaseFitOffset(PositionWave, PhaseWave, ErrorWave)
	Wave PositionWave, PhaseWave, ErrorWave
	
	Variable timerRefNum = Startmstimer
	
	Variable/G num_pos = numpnts(PositionWave)
	Make/o/n=(num_pos,num_x, number_v) yy xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Make/D/O W_coef={43, 48.7e-6}	// initial guesses {alpha, shift}
	Make/D/O Epsilonwave={1e-1,1e-6}	// Amount to vary the fit parameters by in each iteration
	FuncFit/N/M=2 PhiOffset W_coef PhaseWave /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave 
	
	AppendFittedPolWave(PositionWave, PhaseWave)	
	
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




Function CalculateConstants()
	Variable/G Vwire = (6035 * HVmon)	 +0.7							// Applied electrode voltage [V]
	Variable/G d = sqrt((R+gap)^2 - R^2)								// Distance from ground plane to the "infinite wire" [m]
	Variable/G lambda = Vwire*pi*eps0 / asinh(d/R)				// Charge density of the "infinite wire"
	Variable/G E0 = -lambda / (pi*eps0)		
	Variable/G s = (L*hbar*2*pi)/(m*a)								// beam separation at center of grad e region
	//Variable lambda2 = Vwire*2*pi*eps0/ln((gap+R+d)/(gap+R-d))
	//print lambda
	//print lambda2
	//print lambda2*2
End



Function NormalizeVelocities(velocity, sigma, n)			//Find the normalization for the density profile. 
	Variable velocity, sigma, n

	NVAR num_pos
	Make/o/n=(number_v) vel dens
	vel = velocity + (x * v_step - 0.5 * number_v * v_step)			//make a wave for velocities centered on "velocity"
	dens = v_step * (vel)^n* exp (-(vel-velocity)^2/(2*sigma^2)) 	// prob(v)* v_step to find area under curve
	Variable vel_normalization = area(dens)							// used to normalize Int(prob(v)) to 1
	
	//Display/K=1 dens vs vel
	
	Make/o/n=(num_pos,num_x, number_v) v
	v[][][] =( (velocity) - (.5 * number_v * v_step ))+(v_step * z) // all rows and columns of a given layer have the same velocity
	MatrixOP/O density = powR(v,n)* exp (-magsqr(v-(velocity))/(2*magsqr(sigma))) / vel_normalization	//Probability density distribution. same layering as v
	MatrixOP/O vrecip = rec(v)
end




Function Calc_EField(positionwave)			// Calculates the electric field of the gradient region and makes waves for phase, prob, and contrast
	Wave positionwave						// positions at which to calculate the E-field
	Nvar num_pos, Vwire, d, lambda, E0, s
	
	Wave yy, xx, xxsqr, vrecip
	yy[][][] = positionwave[p]				// populate the yy wave with the positions. Each row contains a different position[p], but the position is the same across columns and layers.
	
	// some auxiliary waves that make the E field calculations more efficient
	MatrixOP/O DplusYY = d+yy			
	MatrixOP/O DminusYY = -1*yy + d
	MatrixOP/O DplusYsep = d + yy - (s * vrecip)
	MatrixOP/O DminusYsep =  s*vrecip - yy + d 
	
	// The electric field magnitude of the primary and diffracted beams
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O Etot_s= (E0)*sqrt(	magsqr(DplusYsep/(magsqr(DplusYsep)+xxsqr) + DminusYSep/(magsqr(DminusYSep)+xxsqr)) + magsqr(xx/(magsqr(DplusYsep)+xxsqr) - xx/(magsqr(DminusYsep)+xxsqr))		)

	// Need the magnitude squared.
	MatrixOP/O EtotMagSqr = magsqr(Etot)
	MatrixOP/O Etot_sMagSqr = magsqr(Etot_s)
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




Function CalcPhiOverPol()
	Wave EtotMagSqr, Etot_sMagSqr, vrecip, IFMphase 
	// Find the phase shift for each velocity component of the main and diffracted beam
	MatrixOP/O dphi = EtotMagSqr * vrecip							//phase shift of first path:   alpha*E^2/v	
	MatrixOP/O dphi_s = Etot_sMagSqr * vrecip					//phase shift of second path a distance "sep" away
	Integrate/dim=1 dphi, dphi_s 
	MatrixOP/O IFMphase_v_overpol = dphi - dphi_s				// phase shift divided by polarizability
		
	// Just for curiousity (and debugging): the phase shift of the average velocity path:
	//IFMphase = IFMphase_v_overpol[p][num_x-1][.5*number_v]		
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





Function Phi(pw, output, y) :FitFunc
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
	if(IFMprob[0]<0)
		IFMprob+=2*pi
	endif
	
	output = IFMprob
	
// the following is wrong because we would integrate just the phase itself over velocity (rather than exp(i phi)):
//	MatrixOP/O ProbWeightedPhase = loc_v_step * density * IFMphase_v	// weight each calculation by the appropriate amount
//	Integrate/DIM=2 ProbWeightedPhase									// must integrate along dv to find total phase shift
//	output = ProbWeightedPhase[p][num_x-1][number_v-1]				// select the column corresponding to the total integral over x and v
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





Function PhiOffset(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
	
	Wave Etot, Etot_s, IFMprob, density, Contrast, IFMphase_v_overpol
	
	MatrixOP/O position_shifted = pw[1] + y				// pw[1] is the y offset fit parameter
	Calc_Efield(position_shifted)							// calculate E-field at each data point
	CalcPhiOverPol()										// calculate phase/polarizability
	
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	
	MatrixOP/O IFMphase_v = pol * IFMphase_v_overpol
	
	// We need to find the total velocity weighted phase shift
	// phi_total = ArcTan( Int( Prob(v) * exp(i * delta_phi) ) ) 
	//			 = ArcTan( Int(Prob(v)*sin(phi)) / Int(Prob(v)*cos(phi)) )
	// the code to implement this follows:
	
	Variable loc_v_step = v_step											
	MatrixOP/O imagpart = loc_v_step * density * sin(IFMphase_v)			// weight imaginary part by prob(v) 
	MatrixOP/O realpart = loc_v_step * density * cos(IFMphase_v)				// weight real part by prob(v)
	Integrate/DIM=2 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][num_x-1][number_v-1]) + magsqr(realpart[p][num_x-1][number_v-1]) )
	IFMprob = atan(imagpart[p][num_x-1][number_v-1]/realpart[p][num_x-1][number_v-1])
	//MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	Unwrap Pi, IFMprob		//Correction from taking the arctan.  
	if(IFMprob[0]<-1)
		IFMprob+=pi
	endif
	
	
	output = IFMprob
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





// Reduces the raw position, phase, and error waves to just the important parts
Function ReducePolWaves(position,phase,error)			//phase OR contrast - doesn't matter
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
	Label left "Phase"; Label bottom "Distance to Ground Plane (m)"
	
	Print PositionReducedWaveName
	Print PhaseReducedWaveName
	Print ErrorReducedWaveName
End




Function RemoveNaNsXYZ(theXWave, theYWave, theZWave) //adapted from Remove Points.ipf
	Wave theXWave
	Wave theYWave
	Wave theZWave

	Variable p, numPoints, numNaNs
	Variable xval, yval, zval
	
	numNaNs = 0
	p = 0											// the loop index
	numPoints = numpnts(theXWave)			// number of times to loop

	do
		xval = theXWave[p]
		yval = theYWave[p]
		zval = theZWave[p]
		if ((numtype(xval)==2) %| (numtype(yval)==2) %| (numtype(zval)==2))		// either is NaN?
			numNaNs += 1
		else										// if not an outlier
			theXWave[p - numNaNs] = xval		// copy to input wave
			theYWave[p - numNaNs] = yval		// copy to input wave
			theZWave[p - numNaNs] = zval
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numNaNs, numNaNs, theXWave, theYWave, theZWave
End



Function AppendFittedPolWave(PositionWave, PhaseWave)
	Wave PositionWave, PhaseWave
	Wave W_coef, position, predictedPhase
	
	If(numpnts(W_coef)==1)						// Determine if the offset was fixed or a fit parameter
		PredictPhase(W_coef[0], offset)			// Run the predict phase routine with the fitted polarizability and fixed offset
	else
		PredictPhase(W_coef[0],W_coef[1])		// Run the predict phase routine with the fitted polarizaiblity and offset
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



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function errfit(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = amplitude*(Erf(ramp*(x-offset))+1)+background
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = offset
	//CurveFitDialog/ w[1] = amplitude
	//CurveFitDialog/ w[2] = background
	//CurveFitDialog/ w[3] = ramp

	return w[1]*(Erf(w[3]*(x-w[0]))+1)+w[2]
End

Function errfitneg(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = -amplitude*(Erf(ramp*(x-offset))-1)+background
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = amplitude
	//CurveFitDialog/ w[1] = ramp
	//CurveFitDialog/ w[2] = offset
	//CurveFitDialog/ w[3] = background

	return -w[0]*(Erf(w[1]*(x-w[2]))-1)+w[3]
End



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Function ExtractEclipse0(poswave, maxFile)
	Wave poswave				// Master position wave
	Variable maxFile			// Last file with eclipse data
	
	Variable lastEclipseIndex = maxFile/5
	String EclipsePosWaveName = NameOfWave(poswave) + "_eclipse"

	Make/O/N=(lastEclipseIndex) $EclipsePosWaveName
	
	Wave EclipsePos = $EclipsePosWaveName
	EclipsePos = poswave
End


// example ExtractEclipseData(pos_h,contrast_h,contrast_error_h,avg_counts_h,avg_counts_error_h,75,"h")
Function ExtractEclipseData(poswave, contrast, contrast_error, counts, counts_error, maxFile, series)
	Wave poswave, contrast, contrast_error, counts, counts_error
	Variable maxFile			// Last file with eclipse data
	String series
	
	Variable lastEclipseIndex = maxFile/5
	String EclipsePosWaveName = NameOfWave(poswave) + "_eclipse"
	Make/O/N=(lastEclipseIndex) $EclipsePosWaveName
	Wave EclipsePos = $EclipsePosWaveName
	EclipsePos = poswave
	
	String EclipsePosPadWaveName = NameOfWave(poswave) + "_eclipse_pad"
	Make/O/N=(maxFile) $EclipsePosPadWaveName
	Wave EclipsePosPad = $EclipsePosPadWaveName
	Variable n=0; Variable i=0
	do	// Initialize variables;continue test
		EclipsePosPad[i] = poswave[n]
		if(mod(i+1, 5)==0 && i>0)
			n += 1
		endif
		i+=1
	while(i<maxFile)										
		
	String EclipseContrastWaveName = NameOfWave(contrast) + "_eclipse"
	String EclipseContrastErrorWaveName = NameOfWave(contrast_error) + "_eclipse"
	String EclipseCountsWaveName = NameOfWave(counts) + "_eclipse"
	String EclipseCountsErrorWaveName = NameOfWave(counts_error) + "_eclipse"
	String EclipseCCWaveName = "CountsContrast_" + series
	String EclipseCCErrorWaveName = "CountsContrast_error_" + series
	
	Make/O/N=(maxFile) $EclipseContrastWaveName
	Make/O/N=(maxFile) $EclipseContrastErrorWaveName
	Make/O/N=(maxFile) $EclipseCountsWaveName
	Make/O/N=(maxFile) $EclipseCountsErrorWaveName
	Make/O/N=(maxFile) $EclipseCCWaveName
	Make/O/N=(maxFile) $EclipseCCErrorWaveName
	
	Wave EclipseContrast = $EclipseContrastWaveName
	Wave EclipseContrastError = $EclipseContrastErrorWaveName
	Wave EclipseCounts = $EclipseCountsWaveName
	Wave EclipseCountsError = $EclipseCountsErrorWaveName
	Wave EclipseCountsContrast = $EclipseCCWaveName
	Wave EclipseCountsContrastError = $EclipseCCErrorWaveName
	
	EclipseContrast = contrast
	EclipseContrastError = contrast_error
	EclipseCounts = counts
	EclipseCountsError = counts_error
	EclipseCountsContrast = EclipseCounts*EclipseContrast
	EclipseCountsContrastError = EclipseCountsContrast*sqrt((EclipseCountsError/EclipseCounts)^2 + (EclipseContrastError/EclipseContrast)^2)
	
	Keep4(EclipseContrast)
	Keep4(EclipseContrastError)
	Keep4(EclipseCounts)
	Keep4(EclipseCountsError)
	Keep4(EclipseCountsContrast)
	Keep4(EclipseCountsContrastError)
	Keep4(EclipsePosPad)
	
	RemoveNaNsXYZABCD(EclipseContrast, EclipseContrastError,EclipseCounts,EclipseCountsError,EclipseCountsContrast,EclipseCountsContrastError,EclipsePosPad)
	
	Variable pnts = numpnts(EclipseContrast)
	WaveStats/R=[pnts-1,pnts-4] EclipseContrast
	EclipseContrast /= V_avg; EclipseContrastError /= V_avg
	WaveStats/R=[pnts-1,pnts-4] EclipseCounts
	EclipseCounts /= V_avg; EclipseCountsError /= V_avg
	WaveStats/R=[pnts-1,pnts-4] EclipseCountsContrast
	EclipseCountsContrast /= V_avg; EclipseCountsContrastError /= V_avg
	
	Display/K=1/W=(40,100,550,600) EclipseContrast vs EclipsePosPad
	ModifyGraph mode=3,marker=8
	ErrorBars $EclipseContrastWaveName Y,wave=(EclipseContrastError,EclipseContrastError)
	AppendToGraph/R EclipseCounts vs EclipsePosPad
	ModifyGraph mode($EclipseCountsWaveName)=3
	ModifyGraph marker($EclipseCountsWaveName)=7
	ModifyGraph rgb($EclipseCountsWaveName)=(0,0,0)
	ErrorBars $EclipseCountsWaveName Y,wave=(EclipseCountsError,EclipseCountsError)
	Label left "\\K(65535,0,0)Normalized Contrast"
	Label bottom "Position (m)"
	Label right "Normalized Counts"
	ModifyGraph axisOnTop(right)=1
	ModifyGraph axisEnab(left)={0.53,1},axisEnab(right)={0.53,1}
	AppendToGraph/l=countscontrast EclipseCountsContrast vs EclipsePosPad
	ModifyGraph axisEnab(countscontrast)={0,0.47}
	ModifyGraph freePos(countscontrast)={0,bottom}
	ModifyGraph lblPosMode(countscontrast)=1
	ModifyGraph mode($EclipseCCWaveName)=3,marker($EclipseCCWaveName)=8
	ModifyGraph rgb($EclipseCCWaveName)=(1,16019,65535)
	ErrorBars $EclipseCCWaveName Y,wave=(EclipseCountsContrastError,EclipseCountsContrastError)
	Label countscontrast "\\K(1,16019,65535)Counts*Contrast"; Label bottom "Position"
	SetAxis/A/N=1 bottom
	ModifyGraph gfSize=10
	
//	Display/K=1/W=(40,350,440,575) EclipseCountsContrast vs EclipsePosPad
//	ModifyGraph mode=3,marker=8
//	ModifyGraph rgb($EclipseCCWaveName)=(1,16019,65535)
//	ErrorBars $EclipseCCWaveName Y,wave=(EclipseCountsContrastError,EclipseCountsContrastError)
//	Label left "Counts*Contrast"; Label bottom "Position"
	
	Make/O/D/N=4 W_coef
	W_coef[0] = {2,40000,.0019200,.01}
	FuncFit/NTHR=1/TBOX=768 errfitneg W_coef  EclipseCountsContrast /X=EclipsePosPad /W=EclipseCountsContrastError /I=1 /D 
	
	String tboxname = "CF_"+EclipseCCWaveName
	TextBox/C/N=$tboxname/X=(-5)/Y=50
	
	Variable polfitoffset = gap-w_coef[2]
	String offsetboxtext = "Offset: " + num2str(polfitoffset) + " m"
	TextBox/C/N=offsetbox/X=10/Y=70 offsetboxtext
	
	print "pol fit offset: " , polfitoffset
End






// Extracts the phase data from a the standard data sets (10 HV on, 5 Ref...)
// minfile is the first hv file
// Use maxFile = 0 to automatically use the last file
// example: ExtractPhaseData(pos_g,phase_g,phase_error_g,91,0,"g")
Function ExtractPhaseData(poswave, phase, phase_error, minFile, maxFile, series)
	Wave poswave, phase, phase_error
	
	Variable minFile, maxFile			
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
	
	Duplicate/O PolFitPhase $RefPhaseWaveName
	Duplicate/O PolFitPhaseError $RefPhaseErrorWaveName
	
	Wave RefPhase = $RefPhaseWaveName
	Wave RefPhaseError = $RefPhaseErrorWaveName
	
		
		
	//	extract ref phase
	if(1)
		n=0; i=0
		do
			RefPhase[i]=NaN
			RefPhaseError[i]=NaN
			
			if(n==9)
				n=0
				i+=5
			else
				n+=1
			endif
			
			i+=1
		while(i<numpnts(RefPhase))
		
		i=0
		do
			if(PolFitPhase[i]==RefPhase[i])
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


	Display/K=1/W=(240,100,640,325) PolFitPhase
	ModifyGraph mode=3,marker=8
	ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
	Label left "Phase"; Label bottom "Files"
	
	Display/K=1/W=(240,350,640,575) RefPhase
	ModifyGraph mode=3,marker=8
	ErrorBars $RefPhaseWaveName Y,wave=(RefPhaseError,RefPhaseError)
	Label left "RefPhase"; Label bottom "Files"
	
	
	
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
	
	
		
	if(1)
		//Unwrap
		i=0
		Variable maxUnwrapIndex = 46
		Variable plusminus = 1	// plusminus = 1 if phase needs to be adjusted up, plusminus = -1 if phase needs to be adjusted down
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<maxUnwrapIndex)
	endif
	
	if(1)
		//Unwrap again
		i=0
		Variable maxUnwrapIndex2 = 130
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<maxUnwrapIndex2)
	endif

	
	if(0)
		//Unwrap special
		PolFitPhase[118] -= 2*pi
	endif

	
	if(1)
		//Fit the reference phase to a quadratic polynomial and adjust the phase
		Make/O/N=3 w_coef
		CurveFit/NTHR=1/TBOX=768 poly 3,  RefPhase /W=RefPhaseError /I=1 /D 
		PolFitPhase -= w_coef[0] + w_coef[1]*x + w_coef[2]*x^2
	endif

	if(0)
		Display/K=1/W=(40,100,440,325) PolFitPhase vs PosPad
		ModifyGraph mode=3,marker=8
		ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
		Label left "Phase"; Label bottom "Distance to Ground Plane (m)"
	endif
	

	if(0)							//Make the phase shift positive
		PolfitPhase*=-1
	endif
	
	

	ReducePolWaves(PosPad,PolFitPhase, PolFitPhaseError)
End






Function Keep4(mywave)
	Wave mywave
	
	Variable i=0
	do
		if(mod(i,5)==0)
			mywave[i]=NaN
		endif
		i+=1
	while(i<numpnts(mywave))
End




Function Keep3(mywave)
	Wave mywave
	
	Variable i=0
	do
		if(mod(i,5)==0 || mod(i-1,5)==0)
			mywave[i]=NaN
		endif
		i+=1
	while(i<numpnts(mywave))
End







Function RemoveNaNsXYZAB(theXWave, theYWave, theZWave, theAWave, theBWave) //adapted from Remove Points.ipf
	Wave theXWave
	Wave theYWave
	Wave theZWave
	Wave theAWave
	Wave theBWave

	Variable p, numPoints, numNaNs
	Variable xval, yval, zval, aval, bval
	
	numNaNs = 0
	p = 0											// the loop index
	numPoints = numpnts(theXWave)			// number of times to loop

	do
		xval = theXWave[p]
		yval = theYWave[p]
		zval = theZWave[p]
		aval = theAWave[p]
		bval = theBWave[p]
		if ((numtype(xval)==2) %| (numtype(yval)==2) %| (numtype(zval)==2) %| (numtype(aval)==2) %| (numtype(bval)==2) )		// either is NaN?
			numNaNs += 1
		else										// if not an outlier
			theXWave[p - numNaNs] = xval		// copy to input wave
			theYWave[p - numNaNs] = yval		// copy to input wave
			theZWave[p - numNaNs] = zval
			theAWave[p - numNaNs] = aval
			theBWave[p - numNaNs] = bval
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numNaNs, numNaNs, theXWave, theYWave, theZWave, theAWave, theBWave
End





Function RemoveNaNsXYZABCD(theXWave, theYWave, theZWave, theAWave, theBWave, theCWave, theDWave) //adapted from Remove Points.ipf
	Wave theXWave
	Wave theYWave
	Wave theZWave
	Wave theAWave
	Wave theBWave
	Wave theCWave
	Wave theDWave

	Variable p, numPoints, numNaNs
	Variable xval, yval, zval, aval, bval, cval, dval
	
	numNaNs = 0
	p = 0											// the loop index
	numPoints = numpnts(theXWave)			// number of times to loop

	do
		xval = theXWave[p]
		yval = theYWave[p]
		zval = theZWave[p]
		aval = theAWave[p]
		bval = theBWave[p]
		cval = theCWave[p]
		dval = theDWave[p]
		if ((numtype(xval)==2) %| (numtype(yval)==2) %| (numtype(zval)==2) %| (numtype(aval)==2) %| (numtype(bval)==2) %| (numtype(cval)==2) %| (numtype(dval)==2))		// either is NaN?
			numNaNs += 1
		else										// if not an outlier
			theXWave[p - numNaNs] = xval		// copy to input wave
			theYWave[p - numNaNs] = yval		// copy to input wave
			theZWave[p - numNaNs] = zval
			theAWave[p - numNaNs] = aval
			theBWave[p - numNaNs] = bval
			theCWave[p - numNaNs] = cval
			theDWave[p - numNaNs] = dval
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numNaNs, numNaNs, theXWave, theYWave, theZWave, theAWave, theBWave, theCWave, theDWave
End




// 		savefringedata("m",0)
// 		ExtractEclipseData(pos_q,contrast_q,contrast_error_q,avg_counts_q,avg_counts_error_q,75,"q")
// 		ExtractPhaseData(pos_m,phase_m,phase_error_m,76,0,"m")
//		phasefit(pos_m_polfit_pad_red,phase_m_polfit_red,phase_error_m_polfit_red)







Function BeamProfile()
	Variable beamWidth = 30e-6
	Variable FirstOrder = 0.55
	Variable SecondOrder = 0.10
	Variable ThirdOrder = 0.05
	
	Variable profilepoints = 1000
	Variable profiledistance = 200e-6
	Make/O/N=(profilepoints) ProfilePosition = 2*x*profiledistance/profilepoints - profiledistance
	
	
End
	