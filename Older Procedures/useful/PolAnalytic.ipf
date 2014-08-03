#pragma rtGlobals=1		// Use modern global access method.

#include ":Pol Fit Routine with Molecules"


Function PredictPhase2IFMmolAn(alpha, offset)
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
Make/D/O/N=40 position = .05 * x/1000/2+positionOffset, predictedPhase
	
	Variable/G num_pos = numpnts(position)
	MatrixOP/O position_shifted = offset + position					// account for y offset
	Make/D/o/n=(num_pos, number_v) yy    						// yy: along a column we increase in y position with each row. Along a row we step along x in each column.
	yy[][] = position_shifted[p]	
	//xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	//MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/D/n=(num_pos) IFMphase IFMprob Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam


	CalculateConstants2()
	
	NormalizeVelocitiesAn(velocity, sigma, n)

	PolPhaseAn()

	SagnacAndGravityPhase()
	
	Phi2IFMmolAn(paramWave, predictedPhase, position)
	
	Duplicate/O Contrast, predictedContrast
	
	CleanupPolWaves()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End



Function PolPhaseAn()
	Wave vrecip,yy
	NVAR AnalyticPrefactor, d,s
	Variable dsqr=d^2
	
	MatrixOP/O IFMphase_v_overpolM = vrecip*(rec(powr(yy,2)-dsqr)-rec(powr(yy-s*vrecip,2)-dsqr))*(-1)
	MatrixOP/O IFMphase_v_overpolP = vrecip*(rec(powr(yy+s*vrecip,2)-dsqr)-rec(powr(yy,2)-dsqr))*(-1)
	MatrixOP/O IFMphase_v_overpolMmol = vrecip*(rec(powr(yy,2)-dsqr)-rec(powr(yy-s*vrecip/2,2)-dsqr))*(-1)
	MatrixOP/O IFMphase_v_overpolPmol = vrecip*(rec(powr(yy+s*vrecip/2,2)-dsqr)-rec(powr(yy,2)-dsqr))*(-1)
	
End



Function Phi2IFMmolAn(pw, output, y) :FitFunc
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
	
	NVAR k
	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
	if(IFMprob[0]-2*pi*k<0)
		IFMprob+=2*pi//*(k+1)
	endif
	
	NVAR AvgSagGravPhi
	output = IFMprob - AvgSagGravPhi
	
	//SetDataFolder saveDFR
End





Function NormalizeVelocitiesAn(velocity, sigma, n)			//Find the normalization for the density profile. 
	Variable velocity, sigma, n

	NVAR num_pos
	Make/D/o/n=(number_v) vel dens
	vel = velocity + (x * v_step - 0.5 * number_v * v_step)			// make a wave for velocities centered on "velocity"
	dens = (vel)^n* exp (-(vel-velocity)^2/(2*sigma^2)) 	// prob(v)
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
		
	
	Make/D/o/n=(num_pos, number_v) v
	v[][] =( (velocity) - (.5 * number_v * v_step ))+(v_step * q) // all columns have the same velocity
	MatrixOP/O density = powR(v,n)* exp (-magsqr(v-(velocity))/(2*magsqr(sigma))) / vel_normalization	//Probability density distribution. same layering as v
	MatrixOP/O vrecip = rec(v)
	
	Make/D/O/N=(num_pos, number_v) density2D = dens[q]
end




Function CalculateConstants2()
//	Variable/G Vwire = (6035 * HVmon)	 +0.7							// Applied electrode voltage [V]
	
	NVAR HVmon, ZeroMon, mass
	Variable gap = 1e-3
	Variable R =1.55/2 *1e-3
	Variable L = 13*25.4e-3
	
	Variable/G Vwire = GetHVKeithley(HVmon, ZeroMon)*1000//; Vwire=9000
	Variable/G d = gap*sqrt(1+2*R/gap)//sqrt((R+gap)^2 - R^2)								// Distance from ground plane to the "infinite wire" [m]
	Variable/G lambda = 2*pi*eps0*vwire*(ln((gap+r+d)/(gap+r-d)))^-1 //Vwire*pi*eps0 / asinh(d/R)				// Charge density of the "infinite wire"
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