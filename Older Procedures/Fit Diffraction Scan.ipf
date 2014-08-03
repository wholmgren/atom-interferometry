#pragma rtGlobals=1		// Use modern global access method.
#include <Remove Points>

#include ":WindowNamer"
#include ":Weighted Average"

// 			INSTRUCTIONS
// 1. Run the procedure     RecreateDiffTableFromProcedure()    to generate a table with the fit parameters and initial guesses
// 2. Quick Fit the background counts scan to a line and input the initial guess setting.
// 3. Adjust initial guesses for mean velocity and std deviation.
// 4. Run  GenVelDist(nSigmas,velPoints)          //Adjust minVel, maxVel, and velStep for your initial guesses, then run DisplayVelocityDistribution() to insure that you're integrating over an appropriate range
// 5. If the velocity distribution looks ok, go ahead, if not, repeat previous step.
// 6. In the DiffractionAAOoptimized(pw, output, pos) function you can select what type of atom you are using
// 6. Run FitDiffraction("filename.txt")   
// 7. Adjust initial guesses to be roughly equal to the fit parameters, then double check that the velocity range still looks ok
// 8. Adjust parameters in FitMultipleScans() and run it to fit a whole diffraction series
// 10. Use WeightedAverage() to find the average velocity for two diffraction sets



Strconstant atomtype = "na"
//Strconstant atomtype = "k"
//Strconstant atomtype = "rb"
////Strconstant atomtype = "sr"


Function FitDiffraction(fileName)
	String fileName
	
	Variable runTimer = 1
	If(runTimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Make/O/T/N=5 wavenames
	
	//   CHOOSE YOUR DATA PATH HERE BY COMMENTING/UNCOMMENTING THE APPROPRIATE PATH
	//String path = "Macintosh HD:Users:holmgren:Desktop:090407:"
	String path = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:090929:"
	//String path = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:090929:"
	//String path =  "C:Documents and Settings:data:Desktop:DAQ:2009:090203:"
	
	LoadAndDisplayScan(path, fileName)
	Wave countsWave = $wavenames[0]; Wave positionWave = $wavenames[1]	//2nd half of other half of LoadAndDisplayScan()'s "return statement"
	
	KillSpikes(positionWave, countsWave)
	
	Variable maxPosition = 1.2								// Counts = 0 beyond this position (in mm). must hold higher orders in fitting if scan range does not include them
	RezeroScan(positionWave, countsWave, maxPosition)

	KillSpikes(positionWave, countsWave)

	Variable clumpFactor = 10								// Factor by which to reduce the number of data points
	Clump(positionWave, countsWave, clumpFactor, fileName)
	Wave countsWaveClumped = $wavenames[2]; Wave positionWaveClumped = $wavenames[3]	// "returned waves"

	MakeErrorWave(countsWaveClumped, fileName, clumpFactor)
	Wave errorWave = $wavenames[5]		// "returned wave"
	
	GenVelDist(nSigmas,velPoints)
	
	if(1)	
	GoFitAAO(positionWaveClumped, countsWaveClumped, errorWave)
	AppendFittedDiffWave(positionWave, countsWave, maxPosition)
	else
	positionWaveClumped*=1000
//	GoFitDiffraction(countsWaveClumped, positionWaveClumped, ErrorWave)
	endif
//	GoFitDoubleGauss(positionWaveClumped, countsWaveClumped, errorWave)
	

	
	If(runTimer)
		Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
	endIf
End


//Constant minVel = 1500
//Constant maxVel =2500
//Constant velStep = 30

Constant nSigmas = 4.5
Constant velPoints = 100

Function GenVelDist(nSigmas,velPoints)
	Variable nSigmas,velPoints
	
	Wave InitialGuesses
	
	Variable AvgV = InitialGuesses[4]//; Variable Sig = InitialGuesses[5]
	
	Variable/G minVel = AvgV*(1-nSigmas/15)
	Variable/G maxVel = AvgV*(1+nSigmas/15)
	Variable/G velStep = (maxVel-minVel)/velPoints
	
	//DisplayVelocityDistribution()
End



Function DisplayVelocityDistribution()
	Wave InitialGuesses
	
	NVAR minVel, maxVel, velStep
	//minVel=0; maxVel=5000; velStep=1
	
	Make/o/n=((maxVel-minVel)/velStep+1) velwave = minVel+x*velStep
	Duplicate/o velwave velprobs
	Velprobs = exp(-(velwave-InitialGuesses[4])^2/(2*InitialGuesses[5]^2))*velwave^3
	Wavestats/Q velprobs
	velprobs /= V_max
	Display/k=1 velprobs vs velwave
	Label left "Relative probability"
	Label bottom "Velocity (m/s)"
	print "Most probable velocity: ", velwave[V_maxloc]
	print "velprobs[0]: ", velprobs[0]
	print "velprobs[max]: ", velprobs[(maxVel-minVel)/velStep]
End




Function SetVelocityDistribution()
	Wave InitialGuesses
	
	NVAR minVel, maxVel, velStep
	
	Make/o/n=((maxVel-minVel)/velStep+1) velwave = minVel+x*velStep
	Duplicate/o velwave velprobs
	Velprobs = exp(-(velwave-InitialGuesses[4])^2/(2*InitialGuesses[5]^2))*velwave^3
	Wavestats/Q velprobs
	velprobs /= V_max
	
	
End


Function LoadAndDisplayScan(path, fileName)
	String path, filename
	//Variable motorCountsDivider

	String datay = fileName + "y"	// counts
	String datax = fileName + "x" 	// motor position

	LoadWave/G/D/W/N/O path + filename
	Duplicate/O wave0 $datay; Duplicate/o wave1 $datax 		// copy counts and motor position into waves named $datay and $datax
	Wave positionWave = $datax; Wave countsWave = $datay;	// create references to the datay and datax waves
	
	//positionWave /= 1000/.1		// convert motor position from units of 100 nm to 1 mm 
	//positionWave /= (1000/.16)	// convert motor position from units of 160 nm to 1 mm
	positionWave /= 1000/.4		// 400 nm to 1 mm
	//positionWave /= 1000/motorCountsDivider

	If(0)
		Display/k=1 countsWave vs positionWave			// display a nice graph (will not ask for confirmation when closing)
		ModifyGraph mode=4,marker=8,msize=2
		TextBox/C/N=text0/F=0/A=LT fileName
	EndIf
	
	Wave/T wavenames
	wavenames[0] = datay; wavenames[1] = datax;		//this is effectively a return statement. can be updated now that igor 6.1 supports WAVE waves.
End





Function RezeroScan(positionWave, countsWave, maxPosition)
	Wave positionWave, countsWave
	Variable maxPosition
	
	Wavestats/q countsWave;
	Variable xmax = positionWave[V_maxloc]		// the position of the maximum count rate
	positionWave -= xmax							// shift x-axis
	countsWave *= (abs(positionWave) < maxPosition); 	// counts = 0 if abs(position) > maxPosition (killspikes will delete these points later)
End





Function KillSpikes(positionWave, countsWave)
	Wave positionWave, countsWave
	
	Duplicate/o countsWave dd2y; differentiate dd2y; differentiate dd2y; dd2y = abs(dd2y)		// find the absolute value of the 2nd derivative of the counts
	
	Variable i=0, pcut=50 , ncut = 15
	if(1)
		countsWave = (dd2y > ncut)*(-1) +  (dd2y<=ncut)*countsWave
		countsWave = (dd2y > pcut*countsWave/250)*(-1) +  (dd2y<=pcut*countsWave/250)*countsWave
	endif
	
	Wavestats/q countsWave
	Do
		If (countsWave[i]<1)
			DeletePoints i,1, countsWave, positionWave
			i -= 1
		Endif
		i += 1
	While(i < V_npnts)
End





Function Clump(positionWave, countsWave, clumpFactor, fileName)
	Wave positionWave, countsWave
	Variable clumpFactor
	String fileName

	String cdatay = "c" + fileName + "y"; String cdatax = "c" + fileName + "x"
	Wavestats/Q countsWave
	Variable dataLength = trunc(V_npnts/clumpFactor)-1		// new number of data points
	Make/O/N=(dataLength) $cdatay=0, $cdatax=0
	Wave countsWaveClumped = $cdatay; Wave positionWaveClumped = $cdatax
	
	variable i, j, ingrp
	do
		countsWaveClumped[j] = sum(countsWave,i,i+clumpFactor)/(clumpFactor+1)
		positionWaveClumped[j] = sum(positionWave,i,i+clumpFactor)/(clumpFactor+1)
		i += clumpFactor
		j += 1
	while(j <= dataLength)

	if(1)
		display/K=1 countsWaveClumped vs positionWaveClumped
		ModifyGraph mode=3,marker=8,msize=2
		TextBox/C/N=text0/F=0/A=LT cdatay
		ModifyGraph log(left)=1
		ModifyGraph rgb=(0,0,0)
	endif
	
	
	String WindowName = fileName+"_1"
	Variable k=1
	do
		if(StringMatch(WinList(WindowName,";",""),""))
			break
		else
			k+=1
			WindowName = fileName+"_"+num2str(k)
		endif
	while(1)
	
	if(k==1)
		DoWindow/C/T $WindowName, fileName + " diffraction"
	else
		DoWindow/C/T $WindowName, fileName + " diffraction " + num2str(k)
	endif
	
	Wave/T wavenames
	wavenames[2] = cdatay; wavenames[3] = cdatax	// "return statement"
end






Function MakeErrorWave(counts, fileName, clumpFactor)
	Wave counts
	String fileName
	Variable clumpFactor
	
	Wave/T wavenames
	
	String errwave = fileName + "error"
	Duplicate/O counts $errwave
	Wave errorWave = $errwave
	errorWave = sqrt(counts/clumpFactor)		
	
	ErrorBars $wavenames[2] Y, wave=(errorWave, errorWave)

	wavenames[5] = errWave
End



Function GoFitAAO(positionWave, countsWave, errorWave)	// all at once
	Wave positionWave, countsWave, errorWave
	Wave InitialGuesses, HoldWave
	
	Variable timeryesno = 0
	if(timeryesno)
		Variable timer = StartMSTimer
	endif
	
	Duplicate/O InitialGuesses W_coef
	
	Wave/T Constraints
 	String hold = ProcessHoldAndConstraintsWave()
	
	if(1)
		FuncFit/N/H=hold DiffractionAAOoptimized W_coef  countsWave /X=positionWave /I=1 /W=errorWave /C=Constraints /R
//	elseif(stringmatch(atomtype, "na"))
//		FuncFit/Q/N/H=hold DiffractionAAOna W_coef  countsWave /X=positionWave /I=1 /W=errorWave /C=Constraints /R
//	elseif(stringmatch(atomtype, "k"))
//		FuncFit/Q/N/H=hold DiffractionAAOk W_coef  countsWave /X=positionWave /I=1 /W=errorWave /C=Constraints /R
//	elseif(stringmatch(atomtype, "rb"))
//		FuncFit/Q/N/H=hold DiffractionAAOrb W_coef  countsWave /X=positionWave /I=1 /W=errorWave /C=Constraints /R
//	else
		//Abort "Atom type not properly selected."
	endif
	
	if(timeryesno)
		Print "Elapsed Time: ", stopMSTimer(timer)*1E-6, "sec"
	endif
	
	PrintResults()
	
	NVAR numHolds
	print "Chi-Squared: ", V_chisq, "	Chi-Squared per DOF: ", V_chisq/(V_npnts-(numpnts(W_coef)-numHolds))
	print ""
End


//Constant a0tocm = 0.148184	// * 10^-24


		//  x1  =  diffraction order spacing
		//        =   diffn angle * 1g-detector throw
		//	    =  lambda_dB / d_grating * throw
		//	= h/mv / dg * throw
		//	=  6.62607E-34 /  (22.9898 * 1.66054e-27 * 1000 ) / 1e-7 * 2397.4 
		//     =  .41611 mm    [for velocity = 1000 m/s]. 	(previously .4020 mm)
		// c2_vel = 2563 if x1 = .4020; 	 c2_vel = 2579 if x1 = .40449; 	c2_vel diff = 0.6%. 	x1 diff = 0.6%
		// relies on knowing d_grating relatively well. how well do we really know it?
		// 1g-det distance 2397.4 +/- 1.1
		// how about mass_K?
		// mass K_avg = 39.0983
		// mass Rb = 85.4678

Constant LiAvgmass = 6.941
Constant Li6mass = 6.015122795
Constant Li6frac = 0.0759
Constant Li7mass = 7.01600455
Constant Li7frac = 0.9241

Constant Na23mass = 22.9897692809

Constant Kavgmass = 39.0983
Constant K39mass = 38.96370668
Constant K39frac = 0.932581
Constant K40mass = 39.96399848
Constant K40frac = 0.000117
Constant K41mass = 40.96182576
Constant K41frac = 0.067302

Constant RbAvgmass = 85.4678
Constant Rb87mass = 86.909180527
Constant Rb87frac = 0.2783
Constant Rb85mass = 84.911789738
Constant Rb85frac = 0.7217




Constant piFactor = 0.3989422804  //       = 1/(2 * pi)^.5

Constant LiFactor = 1378.24
Constant Li6Factor = 1590.39
Constant Li7Factor = 1363.51

Constant NaFactor = 416.114			//	= x1*1000
//Constant NaFactor = 425.434		// for magic grating

Constant KFactor =   244.675			// 1st order displacement (in microns) for atoms at 1000 m/s
Constant K39Factor = 245.52
Constant K40Factor = 239.375
Constant K41Factor = 233.544

Constant RbFactor =   111.93
Constant Rb85Factor =   112.662
Constant Rb87Factor =  110.073

Constant SrFactor =   109.18

Constant KFactorGradE =     82.0958		// L =   0.8044 for grad E


Function DiffractionAAOoptimizedTesting(pw, output, pos) : FitFunc
	Wave pw, output, pos	// pw = parameter wave, output = calculated counts wave, pos = position wave
	
	Variable sb = 2*pw[2]^2		// diffraction order width^2 * 2
	Variable v0 = pw[4]			// mean velocity
	Variable sv = 2*pw[5]^2		// standard deviation of velocity distribution^2 *2
	Variable mfac = pw[6]			// molecule fraction
	Variable p1 = pw[7]			// weighting of 1st order relative to 0th
	Variable p2 = pw[8]			// weighting of 2nd order 
	Variable p3 = pw[9]			// weighting of 3rd order 
	Variable p4 = pw[10]			// weighting of 4th order 
	Variable p5 = pw[11]			// weighting of 5th order 
	Variable p6 = pw[12]			// weighting of 6th order
	Variable p7 = pw[13]			// weighting of 7th order
	Variable p8 = pw[14]			// weighting of 8th order
	Variable p9 = pw[15]			// weighting of 9th order
	Variable p10 = pw[16]			// weighting of 10th order
	Variable par = pw[17]			// Ped. Amp. Ratio: pedestal is this many times as brighter than peak
	
	Variable v, pv, x1, x1m		// v = velocity, pv = probability of velocity v, x1 = diff. order spacing
	
	output = 0				// initialize output
	
	MatrixOP/O xx = pos - pw[0]	//account for x offset in diffraction
	MatrixOP/O xxx0 = -powR(xx,2)/sb
		
			
	Variable psdr = 50				// Ped. Std. Dev. Ratio: pedestal is this many times wider than peak
	
	Variable sbprof = 2*(pw[2]*psdr)^2
	MatrixOp/O xxProf0 = -powR(xx,2)/sbprof
	
	mfac = mfac/(1-mfac)
	
	NVAR minVel, maxVel, velStep
	//minVel=0; maxVel=5000; velStep=1
	Make/O/N=((maxVel-minVel)/VelStep+1) velprobs, velwave, x1wave
	
	Variable i
	Variable pvsum = 0
	For( v=minVel; v<maxVel; v+=velStep )
		pv =exp(-(v-v0)^2/sv)//*v^3		//*(v/v0)^3	// probability of v from Maxwell-Boltzmann distribution
		velprobs[(v-minVel)/velStep]=pv
		velwave[(v-minVel)/velStep]=v
		
		//Select sodium, potassium...
		x1 = NaFactor/v; 		x1m = NaFactor/(2*v)
//		x1 = KFactor/v;		x1m = KFactor/(2*v)
//		x1 = KFactorGradE/v;	x1m = KFactorGradE/(2*v)
//		x1 = RbFactor/v; 		x1m = RbFactor/(2*v)
//		x1 = SrFactor/v; 		x1m = SrFactor/(2*v)
		x1wave[(v-minVel)/velStep] = x1
		
		MatrixOp/O xxx1mp= -powR(xx+x1m,2)/sb; MatrixOp/O xxx1mm = -powR(xx-x1m,2)/sb
		MatrixOp/O xxx3mp= -powR(xx+3*x1m,2)/sb; MatrixOp/O xxx3mm = -powR(xx-3*x1m,2)/sb
		
		//FindLevel/Q xx, x1
		//output[V_LevelX]=pv
		
		output=output+pv*p1*exp(-(xx-x1)^2/sb)
		//MatrixOp/O output=output+pv*((1+mfac)*exp(xxx0)+(p1+p2*mfac)*(exp(-powR(xx+x1,2)/sb)+exp(-powR(xx-x1,2)/sb))+(p2+p4*mfac)*(exp(-powR(xx+2*x1,2)/sb)+exp(-powR(xx-2*x1,2)/sb))+p3*(exp(-powR(xx+3*x1,2)/sb)+exp(-powR(xx-3*x1,2)/sb))+p4*(exp(-powR(xx+4*x1,2)/sb)+exp(-powR(xx-4*x1,2)/sb))+p5*(exp(-powR(xx+5*x1,2)/sb)+exp(-powR(xx-5*x1,2)/sb))+p6*(exp(-powR(xx+6*x1,2)/sb)+exp(-powR(xx-6*x1,2)/sb)))
		//MatrixOp/O output=output+pv*(p7*(exp(-powR(xx+7*x1,2)/sb)+exp(-powR(xx-7*x1,2)/sb))+p8*(exp(-powR(xx+8*x1,2)/sb)+exp(-powR(xx-8*x1,2)/sb))+p9*(exp(-powR(xx+9*x1,2)/sb)+exp(-powR(xx-9*x1,2)/sb))+p10*(exp(-powR(xx+10*x1,2)/sb)+exp(-powR(xx-10*x1,2)/sb)))
		
		//Odd order molecules:
		//MatrixOp/O output=output+pv*mfac*(p1*(exp(-powR(xx+x1m,2)/sb)+exp(-powR(xx-x1m,2)/sb)) + p3*(exp(-powR(xx+3*x1m,2)/sb)+exp(-powR(xx-3*x1m,2)/sb)) + p5*(exp(-powR(xx+5*x1m,2)/sb)+exp(-powR(xx-5*x1m,2)/sb)) + p6*(exp(-powR(xx+6*x1m,2)/sb)+exp(-powR(xx-6*x1m,2)/sb)) + p7*(exp(-powR(xx+7*x1m,2)/sb)+exp(-powR(xx-7*x1m,2)/sb)) + p8*(exp(-powR(xx+8*x1m,2)/sb)+exp(-powR(xx-8*x1m,2)/sb))  )
		//MatrixOp/O output=output+pv*mfac*(p9*(exp(-powR(xx+9*x1m,2)/sb)+exp(-powR(xx-9*x1m,2)/sb)) + p10*(exp(-powR(xx+10*x1m,2)/sb)+exp(-powR(xx-10*x1m,2)/sb)))
		
		//Beam profile
		//MatrixOp/O output=output+pv*par*(exp(xxProf0)+p1*(exp(-powR(xx+x1,2)/sbprof)+exp(-powR(xx-x1,2)/sbprof))+p2*(exp(-powR(xx+2*x1,2)/sbprof)+exp(-powR(xx-2*x1,2)/sbprof))+p3*(exp(-powR(xx+3*x1,2)/sbprof)+exp(-powR(xx-3*x1,2)/sbprof)))
		
		pvsum += pv
	Endfor
	
	//velprobs/=pvsum
	
	Variable normconst = pw[1]/pvsum
	MatrixOp/O output = output*normconst + pw[3]	//normalize and add background
End


Function DiffractionAAOoptimized(pw, output, pos) : FitFunc
	Wave pw, output, pos	// pw = parameter wave, output = calculated counts wave, pos = position wave
	
	Variable sb = 2*pw[2]^2		// diffraction order width^2 * 2
	Variable v0 = pw[4]			// mean velocity
	Variable sv = 2*pw[5]^2		// standard deviation of velocity distribution^2 *2
	Variable mfac = pw[6]			// molecule fraction
	Variable p1 = pw[7]			// weighting of 1st order relative to 0th
	Variable p2 = pw[8]			// weighting of 2nd order 
	Variable p3 = pw[9]			// weighting of 3rd order 
	Variable p4 = pw[10]			// weighting of 4th order 
	Variable p5 = pw[11]			// weighting of 5th order 
	Variable p6 = pw[12]			// weighting of 6th order
	Variable p7 = pw[13]			// weighting of 7th order
	Variable p8 = pw[14]			// weighting of 8th order
	Variable p9 = pw[15]			// weighting of 9th order
	Variable p10 = pw[16]			// weighting of 10th order
	Variable par = pw[17]			// Ped. Amp. Ratio: pedestal is this many times as brighter than peak
	
	Variable v, pv, x1, x1m		// v = velocity, pv = probability of velocity v, x1 = diff. order spacing
	
	output = 0				// initialize output
	
	MatrixOP/O xx = pos - pw[0]	//account for x offset in diffraction
	MatrixOP/O xxx0 = -powR(xx,2)/sb
		
			
	Variable psdr = 50				// Ped. Std. Dev. Ratio: pedestal is this many times wider than peak
	
	Variable sbprof = 2*(pw[2]*psdr)^2
	MatrixOp/O xxProf0 = -powR(xx,2)/sbprof
	
	mfac = mfac/(1-mfac)
	
	NVAR minVel, maxVel, velStep
	//minVel=0; maxVel=5000; velStep=1
	//Make/O/N=((maxVel-minVel)/VelStep+1) velprobs, velwave
	
	Variable pvsum = 0
	For( v=minVel; v<maxVel; v+=velStep )
		pv =exp(-(v-v0)^2/sv)*v^3		//*(v/v0)^3	// probability of v from Maxwell-Boltzmann distribution
		//velprobs[(v-minVel)/velStep]=pv
		//velwave[(v-minVel)/velStep]=v
		
		//Select sodium, potassium...
		x1 = NaFactor/v; 		x1m = NaFactor/(2*v)
//		x1 = KFactor/v;		x1m = KFactor/(2*v)
//		x1 = KFactorGradE/v;	x1m = KFactorGradE/(2*v)
//		x1 = RbFactor/v; 		x1m = RbFactor/(2*v)
//		x1 = SrFactor/v; 		x1m = SrFactor/(2*v)
		
		MatrixOp/O xxx1mp= -powR(xx+x1m,2)/sb; MatrixOp/O xxx1mm = -powR(xx-x1m,2)/sb
		MatrixOp/O xxx3mp= -powR(xx+3*x1m,2)/sb; MatrixOp/O xxx3mm = -powR(xx-3*x1m,2)/sb
		
		//Atoms and even order molecules:
		MatrixOp/O output=output+pv*((1+mfac)*exp(xxx0)+(p1+p2*mfac)*(exp(-powR(xx+x1,2)/sb)+exp(-powR(xx-x1,2)/sb))+(p2+p4*mfac)*(exp(-powR(xx+2*x1,2)/sb)+exp(-powR(xx-2*x1,2)/sb))+p3*(exp(-powR(xx+3*x1,2)/sb)+exp(-powR(xx-3*x1,2)/sb))+p4*(exp(-powR(xx+4*x1,2)/sb)+exp(-powR(xx-4*x1,2)/sb))+p5*(exp(-powR(xx+5*x1,2)/sb)+exp(-powR(xx-5*x1,2)/sb))+p6*(exp(-powR(xx+6*x1,2)/sb)+exp(-powR(xx-6*x1,2)/sb)))
		MatrixOp/O output=output+pv*(p7*(exp(-powR(xx+7*x1,2)/sb)+exp(-powR(xx-7*x1,2)/sb))+p8*(exp(-powR(xx+8*x1,2)/sb)+exp(-powR(xx-8*x1,2)/sb))+p9*(exp(-powR(xx+9*x1,2)/sb)+exp(-powR(xx-9*x1,2)/sb))+p10*(exp(-powR(xx+10*x1,2)/sb)+exp(-powR(xx-10*x1,2)/sb)))
		
		//Odd order molecules:
		MatrixOp/O output=output+pv*mfac*(p1*(exp(-powR(xx+x1m,2)/sb)+exp(-powR(xx-x1m,2)/sb)) + p3*(exp(-powR(xx+3*x1m,2)/sb)+exp(-powR(xx-3*x1m,2)/sb)) + p5*(exp(-powR(xx+5*x1m,2)/sb)+exp(-powR(xx-5*x1m,2)/sb)) + p6*(exp(-powR(xx+6*x1m,2)/sb)+exp(-powR(xx-6*x1m,2)/sb)) + p7*(exp(-powR(xx+7*x1m,2)/sb)+exp(-powR(xx-7*x1m,2)/sb)) + p8*(exp(-powR(xx+8*x1m,2)/sb)+exp(-powR(xx-8*x1m,2)/sb))  )
		MatrixOp/O output=output+pv*mfac*(p9*(exp(-powR(xx+9*x1m,2)/sb)+exp(-powR(xx-9*x1m,2)/sb)) + p10*(exp(-powR(xx+10*x1m,2)/sb)+exp(-powR(xx-10*x1m,2)/sb)))
		
		//Beam profile
		MatrixOp/O output=output+pv*par*(exp(xxProf0)+p1*(exp(-powR(xx+x1,2)/sbprof)+exp(-powR(xx-x1,2)/sbprof))+p2*(exp(-powR(xx+2*x1,2)/sbprof)+exp(-powR(xx-2*x1,2)/sbprof))+p3*(exp(-powR(xx+3*x1,2)/sbprof)+exp(-powR(xx-3*x1,2)/sbprof)))
		
		pvsum += pv
	Endfor
	
	//velprobs/=pvsum
	
	Variable normconst = pw[1]/pvsum
	MatrixOp/O output = output*normconst + pw[3]	//normalize and add background
End




Function DiffractionAAOna(pw, output, pos) : FitFunc
	Wave pw, output, pos	// pw = parameter wave, output = calculated counts wave, pos = position wave
	
	Variable sb = 2*pw[2]^2		// diffraction order width^2 * 2
	Variable v0 = pw[4]			// mean velocity
	Variable sv = 2*pw[5]^2		// standard deviation of velocity distribution^2 *2
	Variable mfac = pw[6]			// molecule fraction
	Variable p1 = pw[7]			// weighting of 1st order relative to 0th
	Variable p2 = pw[8]			// weighting of 2nd order 
	Variable p3 = pw[9]			// weighting of 3rd order 
	Variable p4 = pw[10]			// weighting of 4th order 
	Variable p5 = pw[11]			// weighting of 5th order 
	Variable p6 = pw[12]			// weighting of 6th order
	Variable p7 = pw[13]			// weighting of 7th order
	Variable p8 = pw[14]			// weighting of 8th order
	Variable p9 = pw[15]			// weighting of 9th order
	Variable p10 = pw[16]			// weighting of 10th order
	
	Variable v, pv, x1, x1m		// v = velocity, pv = probability of velocity v, x1 = diff. order spacing
	
	output = 0				// initialize output
	
	MatrixOP/O xx = pos - pw[0]	//account for x offset in diffraction
	MatrixOP/O xxx0 = -powR(xx,2)/sb
		
	Variable par = pw[15]			// Ped. Amp. Ratio: pedestal is this many times as brighter than peak
	Variable psdr = 50				// Ped. Std. Dev. Ratio: pedestal is this many times wider than peak
	
	Variable sbprof = 2*(pw[2]*psdr)^2
	MatrixOp/O xxProf0 = -powR(xx,2)/sbprof
	
	NVAR minVel, maxVel, velStep
	Variable pvsum = 0
	For( v=minVel; v<maxVel; v+=velStep )
		pv =exp(-(v-v0)^2/sv)*v^3		//*(v/v0)^3	// probability of v from Maxwell-Boltzmann distribution
		
		x1 = NaFactor/v; 		x1m = NaFactor/(2*v)
		
		//Atoms and even order molecules:
		MatrixOp/O output=output+pv*((1+mfac)*exp(xxx0)+(p1+p2*mfac)*(exp(-powR(xx+x1,2)/sb)+exp(-powR(xx-x1,2)/sb))+(p2+p4*mfac)*(exp(-powR(xx+2*x1,2)/sb)+exp(-powR(xx-2*x1,2)/sb))+p3*(exp(-powR(xx+3*x1,2)/sb)+exp(-powR(xx-3*x1,2)/sb))+p4*(exp(-powR(xx+4*x1,2)/sb)+exp(-powR(xx-4*x1,2)/sb))+p5*(exp(-powR(xx+5*x1,2)/sb)+exp(-powR(xx-5*x1,2)/sb))+p6*(exp(-powR(xx+6*x1,2)/sb)+exp(-powR(xx-6*x1,2)/sb)))
		MatrixOp/O output=output+pv*(p7*(exp(-powR(xx+7*x1,2)/sb)+exp(-powR(xx-7*x1,2)/sb))+p8*(exp(-powR(xx+8*x1,2)/sb)+exp(-powR(xx-8*x1,2)/sb))+p9*(exp(-powR(xx+9*x1,2)/sb)+exp(-powR(xx-9*x1,2)/sb))+p10*(exp(-powR(xx+10*x1,2)/sb)+exp(-powR(xx-10*x1,2)/sb)))
		
		//Odd order molecules:
		MatrixOp/O output=output+pv*mfac*(p1*(exp(-powR(xx+x1m,2)/sb)+exp(-powR(xx-x1m,2)/sb)) + p3*(exp(-powR(xx+3*x1m,2)/sb)+exp(-powR(xx-3*x1m,2)/sb)) + p5*(exp(-powR(xx+5*x1m,2)/sb)+exp(-powR(xx-5*x1m,2)/sb)) )
		
		//Beam profile
		MatrixOp/O output=output+pv*par*(exp(xxProf0)+p1*(exp(-powR(xx+x1,2)/sbprof)+exp(-powR(xx-x1,2)/sbprof))+p2*(exp(-powR(xx+2*x1,2)/sbprof)+exp(-powR(xx-2*x1,2)/sbprof))+p3*(exp(-powR(xx+3*x1,2)/sbprof)+exp(-powR(xx-3*x1,2)/sbprof)))
		
		pvsum += pv
	Endfor
	
	Variable normconst = pw[1]/pvsum
	MatrixOp/O output = output*normconst + pw[3]	//normalize and add background
End





Function DiffractionAAOk(pw, output, pos) : FitFunc
	Wave pw, output, pos	// pw = parameter wave, output = calculated counts wave, pos = position wave
	
	Variable sb = 2*pw[2]^2		// diffraction order width^2 * 2
	Variable v0 = pw[4]			// mean velocity
	Variable sv = 2*pw[5]^2		// standard deviation of velocity distribution^2 *2
	Variable mfac = pw[6]			// molecule fraction
	Variable p1 = pw[7]			// weighting of 1st order relative to 0th
	Variable p2 = pw[8]			// weighting of 2nd order 
	Variable p3 = pw[9]			// weighting of 3rd order 
	Variable p4 = pw[10]			// weighting of 4th order 
	Variable p5 = pw[11]			// weighting of 5th order 
	Variable p6 = pw[12]			// weighting of 6th order
	Variable p7 = pw[13]			// weighting of 7th order
	Variable p8 = pw[14]			// weighting of 8th order
	Variable p9 = pw[15]			// weighting of 9th order
	Variable p10 = pw[16]			// weighting of 10th order
	
	Variable v, pv, x1, x1m		// v = velocity, pv = probability of velocity v, x1 = diff. order spacing
	
	output = 0				// initialize output
	
	Duplicate/O output out39, out40, out41
	
	MatrixOP/O xx = pos - pw[0]	//account for x offset in diffraction
	MatrixOP/O xxx0 = -powR(xx,2)/sb
		
	Variable par = pw[15]			// Ped. Amp. Ratio: pedestal is this many times as brighter than peak
	Variable psdr = 50				// Ped. Std. Dev. Ratio: pedestal is this many times wider than peak
	
	Variable sbprof = 2*(pw[2]*psdr)^2
	MatrixOp/O xxProf0 = -powR(xx,2)/sbprof
	
	NVAR minVel, maxVel, velStep
	Variable pvsum = 0
	For( v=minVel; v<maxVel; v+=velStep )
		pv =exp(-(v-v0)^2/sv)*v^3		//*(v/v0)^3	// probability of v from Maxwell-Boltzmann distribution

		x1 = K39Factor/v; 		x1m = K39Factor/(2*v)
		
		//Atoms and even order molecules:
		MatrixOp/O out39=out39+pv*((1+mfac)*exp(xxx0)+(p1+p2*mfac)*(exp(-powR(xx+x1,2)/sb)+exp(-powR(xx-x1,2)/sb))+(p2+p4*mfac)*(exp(-powR(xx+2*x1,2)/sb)+exp(-powR(xx-2*x1,2)/sb))+p3*(exp(-powR(xx+3*x1,2)/sb)+exp(-powR(xx-3*x1,2)/sb))+p4*(exp(-powR(xx+4*x1,2)/sb)+exp(-powR(xx-4*x1,2)/sb))+p5*(exp(-powR(xx+5*x1,2)/sb)+exp(-powR(xx-5*x1,2)/sb))+p6*(exp(-powR(xx+6*x1,2)/sb)+exp(-powR(xx-6*x1,2)/sb)))
		MatrixOp/O out39=out39+pv*(p7*(exp(-powR(xx+7*x1,2)/sb)+exp(-powR(xx-7*x1,2)/sb))+p8*(exp(-powR(xx+8*x1,2)/sb)+exp(-powR(xx-8*x1,2)/sb))+p9*(exp(-powR(xx+9*x1,2)/sb)+exp(-powR(xx-9*x1,2)/sb))+p10*(exp(-powR(xx+10*x1,2)/sb)+exp(-powR(xx-10*x1,2)/sb)))
		
		//Odd order molecules:
		MatrixOp/O out39=out39+pv*mfac*(p1*(exp(-powR(xx+x1m,2)/sb)+exp(-powR(xx-x1m,2)/sb)) + p3*(exp(-powR(xx+3*x1m,2)/sb)+exp(-powR(xx-3*x1m,2)/sb)) + p5*(exp(-powR(xx+5*x1m,2)/sb)+exp(-powR(xx-5*x1m,2)/sb)) )
		
		//Beam profile
		MatrixOp/O out39=out39+pv*par*(exp(xxProf0)+p1*(exp(-powR(xx+x1,2)/sbprof)+exp(-powR(xx-x1,2)/sbprof))+p2*(exp(-powR(xx+2*x1,2)/sbprof)+exp(-powR(xx-2*x1,2)/sbprof))+p3*(exp(-powR(xx+3*x1,2)/sbprof)+exp(-powR(xx-3*x1,2)/sbprof)))
		
		
		x1 = K40Factor/v; 		x1m = K40Factor/(2*v)
		
		
		MatrixOp/O xxx1mp= -powR(xx+x1m,2)/sb; MatrixOp/O xxx1mm = -powR(xx-x1m,2)/sb
		MatrixOp/O xxx3mp= -powR(xx+3*x1m,2)/sb; MatrixOp/O xxx3mm = -powR(xx-3*x1m,2)/sb
		
		//Atoms and even order molecules:
		MatrixOp/O out40=out40+pv*((1+mfac)*exp(xxx0)+(p1+p2*mfac)*(exp(-powR(xx+x1,2)/sb)+exp(-powR(xx-x1,2)/sb))+(p2+p4*mfac)*(exp(-powR(xx+2*x1,2)/sb)+exp(-powR(xx-2*x1,2)/sb))+p3*(exp(-powR(xx+3*x1,2)/sb)+exp(-powR(xx-3*x1,2)/sb))+p4*(exp(-powR(xx+4*x1,2)/sb)+exp(-powR(xx-4*x1,2)/sb))+p5*(exp(-powR(xx+5*x1,2)/sb)+exp(-powR(xx-5*x1,2)/sb))+p6*(exp(-powR(xx+6*x1,2)/sb)+exp(-powR(xx-6*x1,2)/sb)))
		MatrixOp/O out40=out40+pv*(p7*(exp(-powR(xx+7*x1,2)/sb)+exp(-powR(xx-7*x1,2)/sb))+p8*(exp(-powR(xx+8*x1,2)/sb)+exp(-powR(xx-8*x1,2)/sb))+p9*(exp(-powR(xx+9*x1,2)/sb)+exp(-powR(xx-9*x1,2)/sb))+p10*(exp(-powR(xx+10*x1,2)/sb)+exp(-powR(xx-10*x1,2)/sb)))
		
		//Odd order molecules:
		MatrixOp/O out40=out40+pv*mfac*(p1*(exp(-powR(xx+x1m,2)/sb)+exp(-powR(xx-x1m,2)/sb)) + p3*(exp(-powR(xx+3*x1m,2)/sb)+exp(-powR(xx-3*x1m,2)/sb)) + p5*(exp(-powR(xx+5*x1m,2)/sb)+exp(-powR(xx-5*x1m,2)/sb)) )
		
		//Beam profile
		MatrixOp/O out40=out40+pv*par*(exp(xxProf0)+p1*(exp(-powR(xx+x1,2)/sbprof)+exp(-powR(xx-x1,2)/sbprof))+p2*(exp(-powR(xx+2*x1,2)/sbprof)+exp(-powR(xx-2*x1,2)/sbprof))+p3*(exp(-powR(xx+3*x1,2)/sbprof)+exp(-powR(xx-3*x1,2)/sbprof)))
		
		
		
		x1 = K41Factor/v; 		x1m = K41Factor/(2*v)
		
		MatrixOp/O xxx1mp= -powR(xx+x1m,2)/sb; MatrixOp/O xxx1mm = -powR(xx-x1m,2)/sb
		MatrixOp/O xxx3mp= -powR(xx+3*x1m,2)/sb; MatrixOp/O xxx3mm = -powR(xx-3*x1m,2)/sb
		
		//Atoms and even order molecules:
		MatrixOp/O out41=out41+pv*((1+mfac)*exp(xxx0)+(p1+p2*mfac)*(exp(-powR(xx+x1,2)/sb)+exp(-powR(xx-x1,2)/sb))+(p2+p4*mfac)*(exp(-powR(xx+2*x1,2)/sb)+exp(-powR(xx-2*x1,2)/sb))+p3*(exp(-powR(xx+3*x1,2)/sb)+exp(-powR(xx-3*x1,2)/sb))+p4*(exp(-powR(xx+4*x1,2)/sb)+exp(-powR(xx-4*x1,2)/sb))+p5*(exp(-powR(xx+5*x1,2)/sb)+exp(-powR(xx-5*x1,2)/sb))+p6*(exp(-powR(xx+6*x1,2)/sb)+exp(-powR(xx-6*x1,2)/sb)))
		MatrixOp/O out41=out41+pv*(p7*(exp(-powR(xx+7*x1,2)/sb)+exp(-powR(xx-7*x1,2)/sb))+p8*(exp(-powR(xx+8*x1,2)/sb)+exp(-powR(xx-8*x1,2)/sb))+p9*(exp(-powR(xx+9*x1,2)/sb)+exp(-powR(xx-9*x1,2)/sb))+p10*(exp(-powR(xx+10*x1,2)/sb)+exp(-powR(xx-10*x1,2)/sb)))
		
		//Odd order molecules:
		MatrixOp/O out41=out41+pv*mfac*(p1*(exp(-powR(xx+x1m,2)/sb)+exp(-powR(xx-x1m,2)/sb)) + p3*(exp(-powR(xx+3*x1m,2)/sb)+exp(-powR(xx-3*x1m,2)/sb)) + p5*(exp(-powR(xx+5*x1m,2)/sb)+exp(-powR(xx-5*x1m,2)/sb)) )
		
		//Beam profile
		MatrixOp/O out41=out41+pv*par*(exp(xxProf0)+p1*(exp(-powR(xx+x1,2)/sbprof)+exp(-powR(xx-x1,2)/sbprof))+p2*(exp(-powR(xx+2*x1,2)/sbprof)+exp(-powR(xx-2*x1,2)/sbprof))+p3*(exp(-powR(xx+3*x1,2)/sbprof)+exp(-powR(xx-3*x1,2)/sbprof)))
		
		pvsum += pv
	Endfor
	
	output = out39*K39frac+out40*K40frac+out41*K41frac
	
	Variable normconst = pw[1]/pvsum
	MatrixOp/O output = output*normconst + pw[3]	//normalize and add background
End







Function DiffractionAAOrb(pw, output, pos) : FitFunc
	Wave pw, output, pos	// pw = parameter wave, output = calculated counts wave, pos = position wave
	
	Variable sb = 2*pw[2]^2		// diffraction order width^2 * 2
	Variable v0 = pw[4]			// mean velocity
	Variable sv = 2*pw[5]^2		// standard deviation of velocity distribution^2 *2
	Variable mfac = pw[6]			// molecule fraction
	Variable p1 = pw[7]			// weighting of 1st order relative to 0th
	Variable p2 = pw[8]			// weighting of 2nd order 
	Variable p3 = pw[9]			// weighting of 3rd order 
	Variable p4 = pw[10]			// weighting of 4th order 
	Variable p5 = pw[11]			// weighting of 5th order 
	Variable p6 = pw[12]			// weighting of 6th order
	Variable p7 = pw[13]			// weighting of 7th order
	Variable p8 = pw[14]			// weighting of 8th order
	Variable p9 = pw[15]			// weighting of 9th order
	Variable p10 = pw[16]			// weighting of 10th order
	
	Variable v, pv, x1, x1m		// v = velocity, pv = probability of velocity v, x1 = diff. order spacing
	
	output = 0				// initialize output
	
	Duplicate/O output out85 out87
	
	out85 = 0; out87=0
	
	MatrixOP/O xx = pos - pw[0]	//account for x offset in diffraction
	MatrixOP/O xxx0 = -powR(xx,2)/sb
		
	Variable par = pw[15]			// Ped. Amp. Ratio: pedestal is this many times as brighter than peak
	Variable psdr = 50				// Ped. Std. Dev. Ratio: pedestal is this many times wider than peak
	
	Variable sbprof = 2*(pw[2]*psdr)^2
	MatrixOp/O xxProf0 = -powR(xx,2)/sbprof
	
	NVAR minVel, maxVel, velStep
	Variable pvsum = 0
	For( v=minVel; v<maxVel; v+=velStep )
		pv =exp(-(v-v0)^2/sv)*v^3		//*(v/v0)^3	// probability of v from Maxwell-Boltzmann distribution
		
		x1 = Rb85Factor/v; 		x1m = Rb85Factor/(2*v)
		
		//Atoms and even order molecules:
		MatrixOp/O out85=out85+pv*((1+mfac)*exp(xxx0)+(p1+p2*mfac)*(exp(-powR(xx+x1,2)/sb)+exp(-powR(xx-x1,2)/sb))+(p2+p4*mfac)*(exp(-powR(xx+2*x1,2)/sb)+exp(-powR(xx-2*x1,2)/sb))+p3*(exp(-powR(xx+3*x1,2)/sb)+exp(-powR(xx-3*x1,2)/sb))+p4*(exp(-powR(xx+4*x1,2)/sb)+exp(-powR(xx-4*x1,2)/sb))+p5*(exp(-powR(xx+5*x1,2)/sb)+exp(-powR(xx-5*x1,2)/sb))+p6*(exp(-powR(xx+6*x1,2)/sb)+exp(-powR(xx-6*x1,2)/sb)))
		MatrixOp/O out85=out85+pv*(p7*(exp(-powR(xx+7*x1,2)/sb)+exp(-powR(xx-7*x1,2)/sb))+p8*(exp(-powR(xx+8*x1,2)/sb)+exp(-powR(xx-8*x1,2)/sb))+p9*(exp(-powR(xx+9*x1,2)/sb)+exp(-powR(xx-9*x1,2)/sb))+p10*(exp(-powR(xx+10*x1,2)/sb)+exp(-powR(xx-10*x1,2)/sb)))
		
		//Odd order molecules:
		MatrixOp/O out85=out85+pv*mfac*(p1*(exp(-powR(xx+x1m,2)/sb)+exp(-powR(xx-x1m,2)/sb)) + p3*(exp(-powR(xx+3*x1m,2)/sb)+exp(-powR(xx-3*x1m,2)/sb)) + p5*(exp(-powR(xx+5*x1m,2)/sb)+exp(-powR(xx-5*x1m,2)/sb)) )
		
		//Beam profile
		MatrixOp/O out85=out85+pv*par*(exp(xxProf0)+p1*(exp(-powR(xx+x1,2)/sbprof)+exp(-powR(xx-x1,2)/sbprof))+p2*(exp(-powR(xx+2*x1,2)/sbprof)+exp(-powR(xx-2*x1,2)/sbprof))+p3*(exp(-powR(xx+3*x1,2)/sbprof)+exp(-powR(xx-3*x1,2)/sbprof)))
		
		
		x1 = Rb87Factor/v; 		x1m = Rb87Factor/(2*v)
		
		
		MatrixOp/O xxx1mp= -powR(xx+x1m,2)/sb; MatrixOp/O xxx1mm = -powR(xx-x1m,2)/sb
		MatrixOp/O xxx3mp= -powR(xx+3*x1m,2)/sb; MatrixOp/O xxx3mm = -powR(xx-3*x1m,2)/sb
		
		//Atoms and even order molecules:
		MatrixOp/O out87=out87+pv*((1+mfac)*exp(xxx0)+(p1+p2*mfac)*(exp(-powR(xx+x1,2)/sb)+exp(-powR(xx-x1,2)/sb))+(p2+p4*mfac)*(exp(-powR(xx+2*x1,2)/sb)+exp(-powR(xx-2*x1,2)/sb))+p3*(exp(-powR(xx+3*x1,2)/sb)+exp(-powR(xx-3*x1,2)/sb))+p4*(exp(-powR(xx+4*x1,2)/sb)+exp(-powR(xx-4*x1,2)/sb))+p5*(exp(-powR(xx+5*x1,2)/sb)+exp(-powR(xx-5*x1,2)/sb))+p6*(exp(-powR(xx+6*x1,2)/sb)+exp(-powR(xx-6*x1,2)/sb)))
		MatrixOp/O out87=out87+pv*(p7*(exp(-powR(xx+7*x1,2)/sb)+exp(-powR(xx-7*x1,2)/sb))+p8*(exp(-powR(xx+8*x1,2)/sb)+exp(-powR(xx-8*x1,2)/sb))+p9*(exp(-powR(xx+9*x1,2)/sb)+exp(-powR(xx-9*x1,2)/sb))+p10*(exp(-powR(xx+10*x1,2)/sb)+exp(-powR(xx-10*x1,2)/sb)))
		
		//Odd order molecules:
		MatrixOp/O out87=out87+pv*mfac*(p1*(exp(-powR(xx+x1m,2)/sb)+exp(-powR(xx-x1m,2)/sb)) + p3*(exp(-powR(xx+3*x1m,2)/sb)+exp(-powR(xx-3*x1m,2)/sb)) + p5*(exp(-powR(xx+5*x1m,2)/sb)+exp(-powR(xx-5*x1m,2)/sb)) )
		
		//Beam profile
		MatrixOp/O out87=out87+pv*par*(exp(xxProf0)+p1*(exp(-powR(xx+x1,2)/sbprof)+exp(-powR(xx-x1,2)/sbprof))+p2*(exp(-powR(xx+2*x1,2)/sbprof)+exp(-powR(xx-2*x1,2)/sbprof))+p3*(exp(-powR(xx+3*x1,2)/sbprof)+exp(-powR(xx-3*x1,2)/sbprof)))
		
		pvsum += pv
	Endfor
	
	output = out85*Rb85frac+out87*Rb87frac
	
	Variable normconst = pw[1]/pvsum
	MatrixOp/O output = output*normconst + pw[3]	//normalize and add background
End









Function/S ProcessHoldAndConstraintsWave()
	Variable n
	Wave HoldWave
	Wave/T Constraints, Constraints1, Constraints2
	Duplicate/O/T Constraints1 Constraints1Dup
	Duplicate/O/T Constraints2 Constraints2Dup
	
	String hold = ""
	For(n=0;n<numpnts(HoldWave);n+=1)
		If(HoldWave[n]==1)
			hold += "1"
		Else
			hold += "0"
		EndIf
	EndFor
	
	Variable/G numHolds=0
	For(n=0;n<numpnts(HoldWave);n+=1)
		If(HoldWave[n]==1)
			Constraints1Dup[n]=""; Constraints2Dup[n]=""
			numHolds += 1
		EndIf
	EndFor
	
	Concatenate/O/NP/T {Constraints1Dup, Constraints2Dup}, Constraints
	Variable conLength = numpnts(Constraints)
	For(n=0;n<conLength;n+=1)
		If (stringmatch(Constraints[n],"")==1)
			DeletePoints n,1, Constraints
			n -= 1
		Endif
	EndFor
	
	Return hold
end




Function PrintResults()
	Wave results = W_coef; Wave uncertainties = W_sigma
	
	print " "
	print "center position = ", results[0], "±", uncertainties[0]
	print "normalization = ", results[1], "±", uncertainties[1]
	print "diff. order width = ", results[2], "±", uncertainties[2]
	print "background = ", results[3], "±", uncertainties[3]
	//print " "
	print "mean velocity = ", results[4], "±", uncertainties[4]
	print "standard deviation = ", results[5], "±", uncertainties[5]
	print "v/sigma = ", results[4]/results[5]
	//print " "
	print "molecule fraction =  ", results[6], "±", uncertainties[6]
	//print " "
	print "1st order fraction = ", results[7], "±", uncertainties[7]
	print "2nd order fraction = ", results[8], "±", uncertainties[8]
	print "3rd order fraction = ", results[9], "±", uncertainties[9]
	print "4th order fraction = ", results[10], "±", uncertainties[10]
	print "5th order fraction = ", results[11], "±", uncertainties[11]
	print "6th order fraction = ", results[12], "±", uncertainties[12]
	print "7th order fraction = ", results[13], "±", uncertainties[13]
	print "8th order fraction = ", results[14], "±", uncertainties[14]
	print "9th order fraction = ", results[15], "±", uncertainties[15]
	print "10th order fraction = ", results[16], "±", uncertainties[16]
	//print " "
	print "Ped. Amp. Ratio = ", results[17], "±", uncertainties[17]
	//print " "
	

	
	TextBox/C/N=text1 /A=RT "vel = " + num2str(results[4]) + " ± " + num2str(uncertainties[4]) + "\rsig = " + num2str(results[5]) + " ± " + num2str(uncertainties[5]) + "\rmfrac = " + num2str(results[6]) + " ± " + num2str(uncertainties[6])
End



Function AppendResults()
	Wave w_coef, w_sigma
	Wave/T paramList
	
	String text = "", textTotal = "\\Zr060"
	Variable n = 0
	do
		sprintf text, "%s = %g ± %g", ParamList[n], w_coef[n], w_sigma[n]
		print text
		textTotal += text
		if(n<numpnts(w_coef)-1)
			textTotal+="\r"
		endif
		n+=1
	while(n<numpnts(w_coef))
	print textTotal
	TextBox/C/N=text1 /A=RT textTotal
End

//
//	sprintf RbStr, "Rb = %4.2f ± %1.2f (%1.2f%)\r", RbAvg, RbStdErr, RbPercentError
//	sprintf KStr, "K = %4.2f ± %1.2f (%1.2f%)\r", KAvg, KStdErr, KPercentError
//	sprintf NaStr, "Na = %4.2f ± %1.2f (%1.2f%)", NaAvg, NaStdErr, NaPercentError
//			
//	LegendString = "\\s(RbCMol) "+RbStr+"\\s(KcMol) "+KStr+"\\s(NaCmol) "+NaStr
//	
//	Legend/C/N=PolAbsTBox/J "\\s(RbCMol) "+RbStr+"\\s(KcMol) "+KStr+"\\s(NaCmol) "+NaStr



Function AppendFittedDiffWave(positionWave, countsWave, maxPosition)
	Wave positionWave, countsWave
	Variable maxPosition
	
	Variable fitPoints = 400
	String fitPosWaveName= "fit_"+NameOfWave(positionWave)
	String fitCountsWaveName = "fit_"+NameOfWave(countsWave)
	Make/O/N=(fitPoints) $fitPosWaveName; Wave fitPosWave = $fitPosWaveName
	Make/O/N=(fitPoints) $fitCountsWaveName; Wave fitCountsWave = $fitCountsWaveName
	
	fitPosWave = -maxPosition + 2*maxPosition/fitPoints * x
	//DiffractionAAOoptimized(W_coef, fitCountsWave, fitPosWave)
	if(1)
		DiffractionAAOoptimized(W_coef, fitCountsWave, fitPosWave)
	elseif(stringmatch(atomtype, "na"))
		DiffractionAAOna(W_coef, fitCountsWave, fitPosWave)
	elseif(stringmatch(atomtype, "k"))
		DiffractionAAOk(W_coef, fitCountsWave, fitPosWave)
	elseif(stringmatch(atomtype, "rb"))
		DiffractionAAOrb(W_coef, fitCountsWave, fitPosWave)
	else
		Abort "Atom type not properly selected."
	endif
	
	AppendToGraph fitCountsWave vs fitPosWave
	ModifyGraph zero(Res_Left)=1
End





///////////////////////////////////////////////////////////////////////////////////////////////
//																															//
//								OTHER FUNCTIONS																		//
//																															//
///////////////////////////////////////////////////////////////////////////////////////////////




macro gomain()
	String name = "maxPos1"
	String velName = name + "Vel"
	String sigName = name + "Sig"
	
	Make/O/N=5 $velName = 0, $sigName =0
	
	main("a1")
	$velName[0] = W_coef[4]
	$sigName[0] = W_coef[5]
	
	main("a2")
	$velName[1] = W_coef[4]
	$sigName[1] = W_coef[5]
	
	main("a4")
	$velName[2] = W_coef[4]
	$sigName[2] = W_coef[5]
	
	Wavestats/R=[0,2] $velName
	$velName[3] = V_avg
	
	Wavestats/R=[0,2] $sigName
	$sigName[3] = V_avg
	
	$velName[4] = ($velName[3] - baseVel[3]) / baseVel[3] * 100
	$sigName[4] = ($sigName[3] - baseSig[3]) / baseSig[3] * 100
	
	AppendToTable $velName, $sigName
end



macro testmacro()
	String velName = "test1"
	String sigName = "test2"
	
	Make/O/N=4 $velName = 0, $sigName =0
end




Function PredictDiffraction()
	Variable mst, dt

	Variable runTimer = 0
	If(runTimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Variable nSigmas=10, velPoints=400
	GenVelDist(nSigmas,velPoints)
	
	Make/D/O/N=8001 initGuessX = x/2000 - 2
	Duplicate/O initGuessX initGuessY
//	Make/D/O W_coef= {0,50,0.027,1,2500,250,.60,.2,.1,0.045,.01,.02,0.05}
//	Make/O/N=13 testParams = {0,200,.03,3,2500,250,.6,.12,.04,.03,.01,.005,.05}
	
	//Display/K=1 testwaveout vs testwave
	//ModifyGraph mode=4,marker=8,msize=2
	DiffractionAAOoptimized(InitialGuesses,initGuessY,initGuessX)
	
	if(0)
		Display/k=1 initGuessY vs initGuessX
		ModifyGraph log(left)=1
		ModifyGraph nticks(bottom)=10
	endif
	
	If(runTimer)
		Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
	endIf	

	initguessx*=1000
End



Function PredictDiffraction2mass()
	Variable mst, dt
//	mst=startMSTimer
	
	Make/O/N=400 initGuessX2m = x/100 - 2
	Duplicate/O initGuessX2m initGuessY2m
//	Make/D/O W_coef= {0,50,0.027,1,2500,250,.60,.2,.1,0.045,.01,.02,0.05}
//	Make/O/N=13 testParams = {0,200,.03,3,2500,250,.6,.12,.04,.03,.01,.005,.05}
	
	//Display/K=1 testwaveout vs testwave
	//ModifyGraph mode=4,marker=8,msize=2
	//DiffractionAAOoptimized2mass(InitialGuesses,initGuessY2m,initGuessX2m)
	
	Display/k=1 initGuessY2m vs initGuessX2m
	ModifyGraph log(left)=1
	ModifyGraph nticks(bottom)=10
	
//	dt=stopMSTimer(mst)
//	print "AAO: ",dt*1E-6,"sec"
End





Function AppendInitialGuess(range)
	Variable range
	
	Make/O/N=200 initGuessX = -range + 2*range*x/200
	Duplicate/O initGuessX initGuessY
	
	DiffractionAAOoptimized(InitialGuesses,initGuessY,initGuessX)
	
	RemoveFromGraph initGuessY
	AppendToGraph initGuessY vs initGuessX
End






Function makeTestWavesAAO()
	Variable mst, dt
	mst=startMSTimer
	
	Make/O/N=400 testwaveAAO = x/200 - 1
	Duplicate/O testwaveAAO testwaveoutAAO
		//Make/D/O W_coef= {0,50,0.027,1,2500,250,.60,.2,.1,0.045,.01,.02,0.05}
//	Make/O/N=13 testParams = {0,200,.03,3,2500,250,.6,.12,.04,.03,.01,.005,.05}
	
	//Display/K=1 testwaveout vs testwave
	//ModifyGraph mode=4,marker=8,msize=2
	DiffractionAAOoptimized(InitialGuesses,testwaveoutAAO,testwaveAAO)
	
	dt=stopMSTimer(mst)
	print "AAO: ",dt*1E-6,"sec"
End



Function VelMot(wavein)
	Wave wavein
	Variable i
	For(i=0; i < numpnts(wavein); i+=1)
		If(mod(i,2)==1)
			wavein[i]=NaN
		EndIf
	EndFor
End


Function VelMot2(velocities_motor)
	Wave velocities_motor
	Variable velmotLength = numpnts(velocities_motor)
	Variable n
	For(n=0;n<velmotLength;n+=1)
		If (velocities_motor[n]==NaN)
			DeletePoints n,1, velocities_motor
			n -= 1
		Endif
	EndFor
End


Function FitMultipleScans(series, start, stop)
	String series
	Variable start
	Variable stop
	//Variable motorCountsDivider = 0.4

	Make/O/N=(stop-start+2) velocities sigmas velunc sigunc motorCountsDividerWave MolFracs MolFracsUnc
	velocities[0]=NaN; sigmas[0]=NaN; MolFracs[0]=NaN

	
	Wave W_coef, W_sigma, velocities, sigmas, velunc, sigunc, MolFracs, MolFracsUnc
	
	Variable i
	For(i = start; i <= stop; i += 1)
	If(1)
		FitDiffraction(series + num2str(i))//, motorCountsDivider)
		velocities[i]=W_coef[4]; sigmas[i]=W_coef[5]; velunc[i]=W_sigma[4]; sigunc[i]=W_sigma[5]; MolFracs[i]=W_coef[6]; MolFracsUnc[i]=W_sigma[6];
	EndIf
	EndFor
	
	Wavestats/Q velocities; Variable AvgVel = V_avg
	Make/O/N=(numpnts(velocities)) VavgPlus1 = AvgVel*(1.005), VavgMinus1 = AvgVel*(.995)
	
	DuplicateSeriesResults(series,1)
	VelocitiesGraph($"velocities_"+series,$"velunc_"+series)
	SigmasGraph($"sigmas_"+series,$"sigunc_"+series)
	MolFracGraph($"MolFrac_"+series,$"MolFracUnc_"+series)
End




Function DuplicateSeriesResults(series, displayToggle)
	String series
	Variable displayToggle
	
	String name
	name = "velocities_" + series; Duplicate/O velocities $name
	name = "velunc_" + series; Duplicate/O velunc $name
	name = "sigmas_" + series; Duplicate/O sigmas $name
	name = "sigunc_" + series; Duplicate/O sigunc $name
	name = "MolFrac_"+series; Duplicate/O MolFracs $name
	name = "MolFracUnc_"+series; Duplicate/O MolFracsUnc $name
End



Function SigmasGraph(sigmasw, siguncw) : Graph
	Wave sigmasw, siguncw
	String sigName = NameOfWave(sigmasw); String sigUncName = NameOfWave(siguncw)
	
	Wavestats/Q sigmasw
	Variable AvgSig = V_avg
	Variable StdSig = V_sdev
	Variable PercentBars = 5 // in units of percentage
	Make/O/N=(numpnts(sigmas)) $sigName+"Plus1" = AvgSig*(1+PercentBars/100), $sigName+"Minus1" = AvgSig*(1-PercentBars/100)
		
	PauseUpdate; Silent 1		// building window...
	Display/K=1 sigmasw,$sigName+"Plus1",$sigName+"Minus1"
	ModifyGraph mode($sigName)=4
	ModifyGraph marker($sigName)=8
	ModifyGraph lStyle($sigName+"Plus1")=3,lStyle($sigName+"Minus1")=3
	ModifyGraph rgb($sigName+"Plus1")=(0,0,0),rgb($sigName+"Minus1")=(0,0,0)
	ModifyGraph grid(bottom)=1
	ModifyGraph nticks(bottom)=8
	ModifyGraph sep(bottom)=4
	Label left "Velocity Distribution (m/s)"
	Label bottom "File index"
	ErrorBars $sigName Y,wave=($sigUncName,$sigUncName)
	TextBox/C/N=text1 /A=RT "Sig_avg = " + num2str(AvgSig) + " ± " + num2str(StdSig)
	TextBox/C/N=text2 /A=LT "Dashed Lines: ± " + num2str(PercentBars) + "%"
	SetAxis bottom 1,numpnts(sigmasw)-1
	
	WindowNamer(sigName)
	
	MinusScansStats(sigmasw)
End



Function VelocitiesGraph(velocitiesw, veluncw) : Graph
	Wave velocitiesw, veluncw
	String velName = NameOfWave(velocitiesw); String velUncName = NameOfWave(veluncw)
	
	Wavestats/Q velocitiesw
	Variable AvgVel = V_avg
	Variable StdVel = V_sdev
	Variable PercentBars = 0.5 // in units of percentage
	Make/O/N=(numpnts(velocities)) $velName+"Plus1" = AvgVel*(1+PercentBars/100), $velName+"Minus1" = AvgVel*(1-PercentBars/100)
	
	PauseUpdate; Silent 1		// building window...
	Display/K=1 velocitiesw,$velName+"Plus1",$velName+"Minus1"
	ModifyGraph mode($velName)=4
	ModifyGraph marker($velName)=8
	ModifyGraph lStyle($velName+"Plus1")=3,lStyle($velName+"Minus1")=3
	ModifyGraph rgb($velName+"Plus1")=(0,0,0),rgb($velName+"Minus1")=(0,0,0)
	ModifyGraph grid(bottom)=1
	ModifyGraph nticks(bottom)=8
	ModifyGraph sep(bottom)=4
	Label left "Velocity (m/s)"
	Label bottom "File index"
	ErrorBars $velName Y,wave=($velUncName,$velUncName)
	TextBox/C/N=text1 /A=RT "V_avg = " + num2str(AvgVel) + " ± " + num2str(StdVel)
	TextBox/C/N=text2 /A=LT "Dashed Lines: ± " + num2str(PercentBars) + "%"
	SetAxis bottom 1,numpnts(velocitiesw)-1
	
	WindowNamer(velName)
	
	MinusScansStats(velocitiesw)
End


Function MolFracGraph(MolFracw, MolFracuncw) : Graph
	Wave MolFracw, MolFracuncw
	String velName = NameOfWave(MolFracw); String velUncName = NameOfWave(MolFracuncw)
	
	Wavestats/Q MolFracw
	Variable AvgVel = V_avg
	Variable StdVel = V_sdev
	Variable PercentBars = 10 // in units of percentage
	Make/O/N=(numpnts(MolFracs)) $velName+"Plus1" = AvgVel*(1+PercentBars/100), $velName+"Minus1" = AvgVel*(1-PercentBars/100)
	
	PauseUpdate; Silent 1		// building window...
	Display/K=1 MolFracw,$velName+"Plus1",$velName+"Minus1"
	ModifyGraph mode($velName)=4
	ModifyGraph marker($velName)=8
	ModifyGraph lStyle($velName+"Plus1")=3,lStyle($velName+"Minus1")=3
	ModifyGraph rgb($velName+"Plus1")=(0,0,0),rgb($velName+"Minus1")=(0,0,0)
	ModifyGraph grid(bottom)=1
	ModifyGraph nticks(bottom)=8
	ModifyGraph sep(bottom)=4
	Label left "Molecule Fraction"
	Label bottom "File index"
	ErrorBars $velName Y,wave=($velUncName,$velUncName)
	TextBox/C/N=text1 /A=RT "V_avg = " + num2str(AvgVel) + " ± " + num2str(StdVel)
	TextBox/C/N=text2 /A=LT "Dashed Lines: ± " + num2str(PercentBars) + "%"
	SetAxis bottom 1,numpnts(MolFracw)-1
	
	WindowNamer(velName)
	
	MinusScansStats(MolFracw)
End



Function MinusScansStats(wavein)
		Wave wavein
		
		Variable N= (numpnts(wavein)-1)/2, i
		Variable MinusAvg=0
		For(i=1; i<numpnts(wavein); i+=1)
			If(Mod(i,2)==0)
				MinusAvg+= wavein[i]
			EndIf
		Endfor
		MinusAvg /= N
		
		Variable StdDevSum = 0
		For(i=1; i<numpnts(wavein); i+=1)
			If(Mod(i,2)==0)
				StdDevSum += (wavein[i]-MinusAvg)^2
			EndIf
		Endfor
		
		Variable StdDev = Sqrt(1/N * StdDevSum)
		TextBox/C/N=text3 /A=RB "Minus Avg = " + num2str(MinusAvg) + " ± " + num2str(StdDev)
End




Function KillFileResults(fileNumber)
	Variable fileNumber
	Wave velocities, sigmas
	velocities[fileNumber]=NaN; sigmas[fileNumber]=NaN
End





Function manualprocessmultiplescans()
	Wave W_coef, W_sigma, velocities, sigmas, velunc, sigunc
	Make/O/N=3 velocities sigmas velunc sigunc motorCountsDividerWave; velocities[0]=NaN; sigmas[0]=NaN; 
//	FitDiffraction("g6",.4); velocities[1]=W_coef[4]; sigmas[1]=W_coef[5]; velunc[1]=W_sigma[4]; sigunc[1]=W_sigma[5]
//	FitDiffraction("g8",.4); velocities[2]=W_coef[4]; sigmas[2]=W_coef[5]; velunc[2]=W_sigma[4]; sigunc[2]=W_sigma[5]
//	main("e10",.1); velocities[3]=W_coef[4]; sigmas[3]=W_coef[5]; velunc[3]=W_sigma[4]; sigunc[3]=W_sigma[5]
//	main("e12",.1); velocities[4]=W_coef[4]; sigmas[4]=W_coef[5]; velunc[4]=W_sigma[4]; sigunc[4]=W_sigma[5]
//	main("e14",.4); velocities[5]=W_coef[4]; sigmas[5]=W_coef[5]; velunc[5]=W_sigma[4]; sigunc[5]=W_sigma[5]
//	main("e16",.4); velocities[6]=W_coef[4]; sigmas[6]=W_coef[5]; velunc[6]=W_sigma[4]; sigunc[6]=W_sigma[5]
//	main("d22",.1); velocities[7]=W_coef[4]; sigmas[7]=W_coef[5]; velunc[7]=W_sigma[4]; sigunc[7]=W_sigma[5]
//	main("d24",.1); velocities[8]=W_coef[4]; sigmas[8]=W_coef[5]; velunc[8]=W_sigma[4]; sigunc[8]=W_sigma[5]
//	main("b26",.4); velocities[9]=W_coef[4]; sigmas[9]=W_coef[5]; velunc[9]=W_sigma[4]; sigunc[9]=W_sigma[5]
//	main("b28",.4); velocities[10]=W_coef[4]; sigmas[10]=W_coef[5]; velunc[10]=W_sigma[4]; sigunc[10]=W_sigma[5]
//	main("b30",.4); velocities[11]=W_coef[4]; sigmas[11]=W_coef[5]; velunc[11]=W_sigma[4]; sigunc[11]=W_sigma[5]
//	main("b32",.4); velocities[12]=W_coef[4]; sigmas[12]=W_coef[5]; velunc[12]=W_sigma[4]; sigunc[12]=W_sigma[5]
//	main("d13"); velocities[13]=W_coef[4]; sigmas[13]=W_coef[5]; velunc[13]=W_sigma[4]; sigunc[13]=W_sigma[5]
//	main("d14"); velocities[14]=W_coef[4]; sigmas[14]=W_coef[5]; velunc[14]=W_sigma[4]; sigunc[14]=W_sigma[5]
//	main("d15"); velocities[15]=W_coef[4]; sigmas[15]=W_coef[5]; velunc[15]=W_sigma[4]; sigunc[15]=W_sigma[5]
//	main("d16"); velocities[16]=W_coef[4]; sigmas[16]=W_coef[5]; velunc[16]=W_sigma[4]; sigunc[16]=W_sigma[5]
End






///////////////////////////////////////////////////////////////////////////////////////////////
//																															//
//								RAW BEAM FITTING FUNCTIONS															//
//																															//
///////////////////////////////////////////////////////////////////////////////////////////////







Function GoFitDoubleGauss(positionWave, countsWave, errorWave)
	Wave positionWave, countsWave, errorWave
	
	Make/D/O W_coef = {.4,366,0,.008,.04,.8}	//initial guesses
	Make/T/O dgconstraints = {"K1 > 0","K3 > 0","K4 > 0","K5 > 0"}
	
	FuncFit/TBOX=256/H="100000" DoubleGauss W_coef countsWave /X=positionWave /I=1 /D /W=errorWave /C=dgconstraints 
	
	Wavestats/Q countsWave; print "Chi-Squared per DOF: ", V_chisq/(V_npnts-13-1)
	
	PrintProfileResults()
End



Function DoubleGauss(params, output, pos) : FitFunc
	Wave params, output, pos
	
	Variable bkg = params[0] 		// background counts
	Variable scale = params[1]		// scaling factor
	Variable x0 = params[2]			// mean
	Variable peakRatio = params[3]	// ratio of pedestal to peak
	Variable peakSig = params[4]	// std dev of peak gaussian
	Variable pedSig = params[5]		// std dev of pedestal
	
	output = 0																							// set output = 0
	output += scale * Exp (-(pos - x0)^2/(2 * peakSig^2))			// add central peak
	output += scale * peakRatio * Exp (-(pos - x0)^2/(2 * pedSig^2))	// add pedestal
	output += bkg																						// add background
End


Function PrintProfileResults()
	Wave results = W_coef; Wave uncertainties = W_sigma
	
	print " "
	print "background = ", results[0], "±", uncertainties[0]
	print "scale = ", results[1], "±", uncertainties[1]
	print "x0 = ", results[2], "±", uncertainties[2]
	print "Peak Ratio = ", results[3], "±", uncertainties[3]
	print "Peak Sig = ", results[4], "±", uncertainties[4]
	print "Ped. Sig = ", results[5], "±", uncertainties[5]
	print " "
	print "Ped Ratios =  ", results[5]/results[4]
	print " "
End




Function DoubleGaussTest()
	Make/O/N=600 testwave = x/100 - 3
	Duplicate/O testwave testwaveout
	Make/O/N=5 testParams = {1,366,0,.008,.025,.8}
	
	//Display/K=1 testwaveout vs testwave
	//ModifyGraph mode=4,marker=8,msize=2
	DoubleGauss(testParams,testwaveout,testwave)
End







// These three functions recreate a convenient table for predicting/fitting diffraction
Function RecreateDiffTableFromFiles()

	NewPath/O/Q myPath "Macintosh HD:Users:holmgren:Documents:School:UA:Cronin labs:Diffraction Analysis:waves:"
	
	String fileName
	Variable index=0
	
	//algorithm from:		DisplayHelpTopic "Loading All of the Files in a Folder"
	do			// Loop through each file in folder
		fileName = IndexedFile(mypath, index, ".ibw")
		print fileName
		if (strlen(fileName) == 0)			// No more files?
			break									// Break out of loop
		endif
		
		LoadWave/H/O/P=mypath fileName
		
		index += 1
	while (1)
	
	Execute "DiffParamsTable()"
End



Function RecreateDiffTableFromProcedure()
	Make/T/O/N=18   KWave= {"K0","K1","K2","K3","K4","K5","K6","K7","K8","K9","K10","K11","K12","K13","K14","K15","K16","K17"}
	Make/T/O/N=18  ParamList= {"center position","normalization","diff. order width","background","mean velocity","standard deviation","molecule fraction","1st order","2nd order","3rd order","4th order","5th order","6th order","7th order","8th order","9th order","10th order","ped amp ratio"}
	Make/T/O/N=18     Constraints1={"","K1>10","","K3>.01","K4<3300","K5>10"," K6>.0001","K7>.3","K8>.01","K9>.01","K10>.001","K11>.001","K12>.001","K13>.0001","K14>.0001","K15>.0001","K16>.0001","K17>.0001"}
	Make/T/O/N=18  Constraints2= {"","","","","K4>100","K5<600","K6<.5","","","","K10<.2","K11<.1","K12<.1","K13<.05","K14<.05","K15<.03","K16<.03","K17<.01"}
	Make/T/O/N=18   Constraints= {"K1>10","K4<3300","K5>10"," K6>.001","K7>.3","K8>.01","K9>.01","K10>.001","K11>.001","K12>.001","K4>100","K5<600","K6<.5","K10<.1","K11<.05","K12<.05"}
	Make/O/N=18  InitialGuesses= {0,300,0.033,0.5,901,40,0.02,0.55,0.05,0.05,0.02,0.01,0.028,0.01,0.005,0.009,0.005,0.003}
	Make/O/N=18  W_coef= {-0.000697805,425.451,0.0312131,0.15,1449.7,110.792,0.001,0.622123,0.114097,0.0511681,0.0481266,0.020466,0.0231163,0.01,0.003,0.003}
	Duplicate/O W_coef W_sigma
	Make/O/N=18  HoldWave= {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1}
	
	PauseUpdate; Silent 1		// building window...
	Edit/W=(141,44,1099,477) KWave,ParamList,InitialGuesses,w_coef,w_sigma,HoldWave,Constraints1,Constraints2,Constraints
	ModifyTable format(Point)=1,width(KWave)=54,width(ParamList)=104,width(InitialGuesses)=84
	ModifyTable width(Constraints1)=90,width(Constraints2)=98,width(Constraints)=92
End







Function GlobalFitDiffraction(Series)
	String Series
	
	Make/O/T/N=1 FitFuncName = "DiffractionAAOoptimized"
	
	Make/O/T/N=(4,3) DataToFit
	SetDimLabel 1,2,'Weights', DataToFit
	DataToFit[][0]={"ca2y","ca4y","ca6y","ca8y"}	//y waves
	DataToFit[][1]={"ca2x","ca4x","ca6x","ca8x"}	//x waves
	DataToFit[][%Weights]={"a2error","a4error","a6error","a8error"}
	
	Make/O/N=(4,21) DataSetLinkage
								// a matrix wave with a row for each data set and N+2 columns, where N is the maximum number of coefficients
								// used by any of the fit functions. It looks like this for a hypothetical case of two functions and four
								// data sets:
								
								//		|	f()	first	last	N			c0	c1	c2	c3	c4	c5
								//	---|------------------------------------------
								//	ds1	|	0	0			100		5	0	1	2	3	4	-1
								//	ds2	|	0	101		150		5	5	6	2	7	4	-1
								//	ds3	|	1	151		220		6	8	9	2	10	11	12
								//	ds4	|	1	221		300		6	13	14	2	15	16	12

								// In this matrix, I imagine fitting to two functions, one of which takes 5 coefficients, the other 6. 
								// Coefficients 0, 1, and 3 for function f1 are local- they have distinct coefficient array indices 
								// everywhere. Coefficient 2 is global- the same coefficient array index is used for every data set. 
								// Coefficient 4 is "shared local" (group global?)- it is global for ds1 and ds2. The corresponding 
								// coefficient for ds3 and ds4 is local- it probably refers to something entirely different. Function 
								// f1 has no coefficient 5- hence the -1. For f2, coefficient 5 is shared between the data sets (ds3 
								// and ds4) which use f2. The column labelled N is the number of coefficients needed by the fit function.
								// The column labelled f() has an index into the FitFuncNames wave.
								// These columns are set up by the function in its private copy. You can set them to zero:
								// The column labelled first contains the point number where that data set starts in the cumulative waves.
								// The column labelled last contains the point number where the last point of that data set resides in the cumulative waves
								
	DataSetLinkage[][0] = 0	// Links to FitFuncName index
	DataSetLinkage[][1] = 0	//	GlobalFit will overwrite
	DataSetLinkage[][2] = 0	//	GlobalFit will overwrite
	DataSetLinkage[][3] = {16,17,18,19}	//	center position column
	DataSetLinkage[][4] = {20,21,22,23}	//	normalization
	DataSetLinkage[][5] = 0	//	diffraction order width
	DataSetLinkage[][6] = 1	//	background
	DataSetLinkage[][7] = 2	//	flow velocity
	DataSetLinkage[][8] = 3	//	"std deviation"
	DataSetLinkage[][9] = 4	//	molecule fraction
	DataSetLinkage[][10] = 5	//	1st order
	DataSetLinkage[][11] = 6	//	2
	DataSetLinkage[][12] = 7	//	3
	DataSetLinkage[][13] = 8	//	4
	DataSetLinkage[][14] = 9	//	5
	DataSetLinkage[][15] = 10	//	6
	DataSetLinkage[][16] = 11	//	7
	DataSetLinkage[][17] = 12	//	8
	DataSetLinkage[][18] = 13	//	9
	DataSetLinkage[][19] = 14	//	10
	DataSetLinkage[][20] = 15	//	pedestal amplitude ratio
	
	
	//DoNewGlobalFit(FitFuncName, DataToFit, CoefDataSetLinkage, CoefWave, CoefNames, ConstraintWave, Options, FitCurvePoints, DoAlertsOnError, [errorName, errorMessage, maxIters, resultWavePrefix, resultDF])
End	

	
	
	
	
Function MakeInitGuessesForGlobal()
	Variable NumberOfDataSets = 8
	
	Wave w_coef, holdwave
	
	Make/O/N=(NumberOfDataSets*18) InitGuessesForGlobal = w_coef, HoldForGlobal = HoldWave
	
	Variable i = 1
	Variable StartIndex
	Variable EndIndex 
	do
		StartIndex = i*18
		EndIndex = 17+i*18
		InitGuessesForGlobal[StartIndex, EndIndex] = w_coef[p-18*(i)]
		HoldForGlobal[StartIndex, EndIndex] = HoldWave[p-18*(i)]
		i+=1
	while(i<NumberOfDataSets)
	
End