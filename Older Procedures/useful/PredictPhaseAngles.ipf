#pragma rtGlobals=1		// Use modern global access method.

#include ":Pol Fit Routine with Molecules"

Function PredictPhase2IFMmolThetaZ(alpha, offset)
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

	Calc_Efield2IFMmolThetaZ(position_shifted)							// calculate E-field at each data point

	CalcPhiOverPol2IFMmol()										// calculate phase/polarizability
	//SagPhase()
	//GravityPhase()
	SagnacAndGravityPhase()
	
	Phi2IFMmol(paramWave, predictedPhase, position)
	
	Duplicate/O Contrast, predictedContrast
	
	CleanupPolWaves()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End



Function Calc_EField2IFMmolThetaZ(positionwave)			// Calculates the electric field of the gradient region and makes waves for phase, prob, and contrast
	Wave positionwave						// positions at which to calculate the E-field
	Nvar num_pos, Vwire, d, lambda, E0, s             //; E0=0
	
	Variable sLoc = s
	
	Wave yy, xx, xxsqr, vrecip
	yy[][][] = positionwave[p]				// populate the yy wave with the positions. Each row contains a different position[p], but initially the position is the same across columns and layers.
											// vrecip is the same across rows and columns, but different across layers
											// the beam propagation direction is along a single row
	Variable thetaZ=0.1e-3	// in rad										
	yy[][][] = positionwave[p]+xx*thetaZ										// modify positions along a row
	
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
End
