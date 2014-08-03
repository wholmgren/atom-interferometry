#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 2.0

#include <Remove Points>
#include ":HV Calibration"
#include ":RemoveNaNs"
#include ":WindowNamer"
#include ":Weighted Average"
#include ":PhysicalConstants"


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
//Constant eps0 = 8.8541878e-12
//Constant hbar = 1.054572e-34
//Constant m=   3.81753e-26				//Na = 22.9897 * 1e-3 / 6.02214e23 	//mass of element in kg: 85.4678  39.0983 22.9897 	
//Constant m =     6.49243e-26			//K natural abundance
//Constant m = 	  6.4761e-26				//K39 only 
//Constant a0tocm = 0.148184

// Geometry constants
//Constant gap = 1.998e-3					//Gap between the ground and cylindrical electrodes [m]
Static Constant gap = 1.000e-3					//Gap between the ground and cylindrical electrodes [m]
//Constant R =  .0063315						//Radius of the electrode [m]
//Constant R =  .7239e-3							//Radius if diameter is .001" larger
//Constant L = .80257						// .80257 8/8/08-3/1/10	used for 2010 pol paper	//.8044  7/8/08-8/8/08	//.7813 pre 7/8/08	//Distance from 1g to the grad E region
Constant L = .8289				// 2 electrode distance

Constant a = 1e-7							//Grating period

// x-axis integration constants
Constant num_x=100							// position integration points along the beam path
Constant x_step =   0.001				// = (gap*35)/num_x		//35 mm from gap to edge of ground plane

// velocity integration constants
Constant number_v = 100					// number of individual velocities for integration
Constant v_step = 15						// separation of the individual velocities in m/s

Constant n=3								//n= 3, Maxwellian;    0, gaussian


//Constant K39masskg = 6.4761e-26
//Constant K40masskg = 6.4761e-26
//Constant K41masskg = 6.4761e-26


Function PolParams()
	String CurrentDataFolder = GetDataFolder(0)
	
	//     General parameters
	Variable velLoc = NumVarOrDefault("velocity", 2000)
	Prompt velLoc, "Velocity (m/s): "
	Variable sigLoc = NumVarOrDefault("sigma", 100)
	Prompt sigLoc, "Sigma (m/s): "
	Variable molFracLoc = NumVarOrDefault("molFrac",0)
	Prompt molFracLoc, "Molecule Fraction (%): "
	Variable HVmonLoc = NumVarOrDefault("HVmon", 2)
	Prompt HVmonLoc, "HVmon (V): "
	Variable ZeroMonLoc = NumVarOrDefault("ZeroMon", .013)
	Prompt ZeroMonLoc, "ZeroMon (V): "
	Variable OffsetLoc = NumVarOrDefault("Offset", 1E-5)
	Prompt OffsetLoc, "Offset (m): "
	String AtomLoc = StrVarOrDefault("Atom", "Sodium")
//	Prompt AtomLoc, "Atom Type: ", popup "Sodium;Potassium;Potassium 39;Rubidium"
	Prompt AtomLoc, "Atom Type: ", popup "Sodium;Potassium Avg Mass;Potassium weighted;Potassium Avg test;Rubidium Avg Mass;Rubidium avg test;Rubidium weighted"
	
	Variable LastEclipseLoc = NumVarOrDefault("LastEclipse", 60)
	Prompt LastEclipseLoc, "Last Eclipse File: "
	Variable LastFileLoc = NumVarOrDefault("LastFile", 200)
	Prompt LastFileLoc, "Last File: "
//	Variable KillEclipseFilesLoc = NumVarOrDefault("KillEclipseFiles", 0)
//	Prompt KillEclipseFilesLoc, "Kill Eclipse Files: yes (1) no (0) "
	String EclipseFilesToKillStrLoc  = StrVarOrDefault("EclipseFilesToKillStr", "")
	Prompt EclipseFilesToKillStrLoc, "List eclipse files to kill sep. by ; "
	
	
	DoPrompt CurrentDataFolder+ " General Parameters", velLoc, sigLoc, HVmonLoc, ZeroMonLoc, OffsetLoc, AtomLoc, LastEclipseLoc, LastFileLoc, EclipseFilesToKillStrLoc, MolFracLoc
	
	
	Variable/G velocity = velLoc
	Variable/G sigma = sigLoc
	Variable/G HVmon = HVmonLoc
	Variable/G ZeroMon = ZeroMonLoc
	Variable/G offset = OffsetLoc
	String/G Atom = AtomLoc
	Variable/G molFrac = molFracLoc
	
	Variable/G mass
	Variable/G molPol			// molPol from Tarnovsky et al (1993)
	Variable/G initPol
	If(stringmatch(Atom, "Sodium"))
		mass = 3.81754e-26	
		molPol = 40
		initPol = 24
	elseif(stringmatch(atom, "Potassium Avg Mass"))
		mass = 6.49242e-26	
		molPol = 77
		initPol = 43
	elseif(stringmatch(atom, "Potassium 39"))
		mass =  6.4761e-26
		molPol = 77
		initPol = 43	
	elseif(stringmatch(atom, "Potassium weighted"))
		Variable/G numIsotopes = 3
		Make/D/O massWave = {K39mass, K40mass, K41mass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {K39frac, K40frac, K41frac}
		molPol = 77
		initPol = 43	
	elseif(stringmatch(atom, "Potassium avg test"))
		Variable/G numIsotopes = 1
		Make/D/O massWave = {Kavgmass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {1}
		molPol = 77
		initPol = 43
	elseif(stringmatch(atom, "Rubidium Avg Mass"))
		mass =   1.41923e-25	
		molPol = 79
		initPol = 47
	elseif(stringmatch(atom, "Rubidium avg test"))
		Variable/G numIsotopes = 1
		Make/D/O massWave = {Rbavgmass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {1}
		molPol = 79
		initPol = 47
	elseif(stringmatch(atom, "Rubidium weighted"))
		Variable/G numIsotopes = 2
		Make/D/O massWave = {Rb85mass, Rb87mass}
		massWave*=amu2kg
		Make/D/O isotopeFracs = {Rb85frac, Rb87frac}
		molPol = 79
		initPol = 47
	EndIf
	
	Variable/G LastEclipse = LastEclipseLoc
	Variable/G LastFile = LastFileLoc
	//Variable/G KillEclipseFiles = KillEclipseFilesLoc
	
	String/G EclipseFilesToKillStr = EclipseFilesToKillStrLoc
	Variable i = 0
	Variable NumEclipseFilesToKill = ItemsInList(EclipseFilesToKillStr)
	if(NumEclipseFilesToKill != 0)
		Make/O/N=(NumEclipseFilesToKill) EclipseFilesToKill
		do
			EclipseFilesToKill[i] = Str2Num(StringFromList(i, EclipseFilesToKillStr))
			i+=1
		while(i < NumEclipseFilesToKill)
		Variable/G KillEclipseFiles = 1
	else
		If(WaveExists(EclipseFilesToKill))
			EclipseFilesToKill = NaN
		EndIf
		Variable/G KillEclipseFiles = 0
	EndIf
	
	
	
	//          HV phase unwrapping
	Variable PlusMinusLoc = NumVarOrDefault("plusminus", 1)
	Prompt plusminusloc, "Add 2pi (1) or subtract 2pi (-1): "
	Variable Unwrap1Loc = NumVarOrDefault("Unwrap1", 0)
//	Prompt Unwrap1Loc, "Turn on first unwrap? 1 or 0"
	Variable Unwrap1EndLoc = NumVarOrDefault("Unwrap1End", 0)
	Prompt Unwrap1EndLoc, "Last file for first unwrapping: "
	Variable Unwrap2Loc = NumVarOrDefault("Unwrap2", 0)
//	Prompt Unwrap2Loc, "Turn on second unwrap? 1 or 0"
	Variable Unwrap2EndLoc = NumVarOrDefault("Unwrap2End", 0)
	Prompt Unwrap2EndLoc, "Last file for second unwrapping: "
	Variable Unwrap3Loc = NumVarOrDefault("Unwrap3", 0)
//	Prompt Unwrap3Loc, "Turn on second unwrap? 1 or 0"
	Variable Unwrap3EndLoc = NumVarOrDefault("Unwrap3End", 0)
	Prompt Unwrap3EndLoc, "Last file for third unwrapping: "
	String HVFilesToKillStrLoc  = StrVarOrDefault("HVFilesToKillStr", "")
	Prompt HVFilesToKillStrLoc, "List HV files to kill separated by ; "
		
//	DoPrompt CurrentDataFolder, plusminusloc, Unwrap1Loc, Unwrap1EndLoc, Unwrap2loc, Unwrap2EndLoc, Unwrap3Loc, Unwrap3EndLoc
	DoPrompt CurrentDataFolder+" HV unwrapping", plusminusloc, Unwrap1EndLoc, Unwrap2EndLoc, Unwrap3EndLoc, HVFilesToKillStrLoc
	
	Variable/G plusminus = plusminusloc
	
	if(Unwrap1EndLoc != 0)
		Variable/G Unwrap1 = 1
		Variable/G Unwrap1End = Unwrap1EndLoc
	else
		Variable/G Unwrap1 = 0
		Variable/G Unwrap1End = 0
	endif
	
	if(Unwrap2EndLoc != 0)
		Variable/G Unwrap2 = 1
		Variable/G Unwrap2End = Unwrap2EndLoc
	else
		Variable/G Unwrap2 = 0
		Variable/G Unwrap2End = 0
	endif
	
	if(Unwrap3EndLoc != 0)
		Variable/G Unwrap3 = 1
		Variable/G Unwrap3End = Unwrap3EndLoc
	else
		Variable/G Unwrap3 = 0
		Variable/G Unwrap3End = 0
	endif
	
	String/G HVFilesToKillStr = HVFilesToKillStrLoc
	i = 0
	Variable NumHVFilesToKill = ItemsInList(HVFilesToKillStr)
	if(NumHVFilesToKill != 0)
		Make/O/N=(NumHVFilesToKill) HVFilesToKill
		do
			HVFilesToKill[i] = Str2Num(StringFromList(i, HVFilesToKillStr))
			i+=1
		while(i < NumHVFilesToKill)
		Variable/G KillHVFiles = 1
	else
		If(WaveExists(HVFilesToKill))
			HVFilesToKill = NaN
		EndIf
		Variable/G KillHVFiles = 0
	EndIf
	
	
	//    Reference phase unwrapping
	Variable PlusMinusRefLoc = NumVarOrDefault("plusminusref", 1)
	Prompt plusminusrefloc, "Add 2pi (1) or subtract 2pi (-1): "
	Variable UnwrapRef1Loc = NumVarOrDefault("UnwrapRef1", 0)
	Variable UnwrapRef1EndLoc = NumVarOrDefault("UnwrapRef1End", 0)
	Prompt UnwrapRef1EndLoc, "Last file for unwrapping: "
	String RefFilesToKillStrLoc  = StrVarOrDefault("RefFilesToKillStr", "")
	Prompt RefFilesToKillStrLoc, "List Ref files to kill separated by ; "
	Variable InvertPhaseLoc = NumVarOrDefault("InvertPhase",1)
	Prompt InvertPhaseLoc, "Multiply phase shift by -1? yes (1) no (0)"
	
	DoPrompt CurrentDataFolder+" Reference phase unwrapping", plusminusrefloc, UnwrapRef1EndLoc, RefFilesToKillStrLoc, InvertPhaseLoc
	
	Variable/G plusminusref = plusminusrefloc
	Variable/G InvertPhase = InvertPhaseLoc
	
	if(UnwrapRef1EndLoc != 0)
		Variable/G UnwrapRef1 = 1
		Variable/G UnwrapRef1End = UnwrapRef1EndLoc
	else
		Variable/G UnwrapRef1 = 0
		Variable/G UnwrapRef1End = 0
	endif
	
	String/G RefFilesToKillStr = RefFilesToKillStrLoc
	i = 0
	Variable NumRefFilesToKill = ItemsInList(RefFilesToKillStr)
	if(NumRefFilesToKill != 0)
		Make/O/N=(NumRefFilesToKill) RefFilesToKill
		do
			RefFilesToKill[i] = Str2Num(StringFromList(i, RefFilesToKillStr))
			i+=1
		while(i < NumRefFilesToKill)
		Variable/G KillRefFiles = 1
	else
		If(WaveExists(RefFilesToKill))
			RefFilesToKill = NaN
		EndIf
		Variable/G KillRefFiles = 0
	EndIf
	
End







// DOUBLE CHECK THAT YOUR DATA FOLDER IS SET CORRECTLY BEFORE RUNNING THESE PROCEDURES!!!
Function WriteGlobalVariables()
	Variable/G velocity = 1946.9
	Variable/G sigma = 118
	
	Variable/G HVmon = 1.49966
	Variable/G ZeroMon = 0.013754
	
	Variable/G offset = 3.0851E-5
	
	Variable/G mass = 6.49243e-26		//K natural abundance
End

// DOUBLE CHECK THAT YOUR DATA FOLDER IS SET CORRECTLY BEFORE RUNNING THESE PROCEDURES!!!
Function WriteUnwrapParams()
	Variable/G LastEclipse = 60
	Variable/G LastFile = 195

	Variable/G plusminus = -1		// plusminus = 1 if phase needs to be adjusted up, plusminus = -1 if phase needs to be adjusted down
	
	Variable/G UnwrapAll = 0
	
	Variable/G Unwrap1 = 1			// Turn on (1) or off (0) the first unwrapping
	Variable/G Unwrap1End = 115		// Last HV file index to unwrap
	
	Variable/G Unwrap2 = 0
	Variable/G Unwrap2End = 150		
	
	Variable/G Unwrap3 = 0 
	Variable/G Unwrap3End = 90
	
	
	Variable/G PlusMinusRef = 1
	
	Variable/G UnwrapRef1 = 1
	Variable/G UnwrapRef1End = 150	// Last Reference file to be unwrapped
	
	Variable/G InvertPhase = 1			// Multiple phase shift by minus 1 if this setting is true (1), otherwise set to false (0)
	
	Variable/G KillHVFiles = 0
	Make/O HVFilesToKill = {134,163}
	
	Variable/G KillRefFiles = 1
	Make/O RefFilesToKill = {180}
	
	Variable/G UnwrapRefSpecial = 0
	Make/O RefFilesToUnwrap = {34}
	
	Variable/G KillEclipseFiles = 1
	Make/O EclipseFilesToKill = {53,28,49,44,50}
End




Function PredictPhase(alpha, offset)
	Variable alpha
	Variable offset
	
	NVAR velocity, sigma
	
	Variable runtimer = 0
	If(runtimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Make/O/N=1 paramWave = {alpha}
	
	Make/O/N=20 position = .1 * x/1000+.0001, predictedPhase
	
	Variable/G num_pos = numpnts(position)
	Make/o/n=(num_pos,num_x, number_v) yy, xx   						// yy: along a column we increase in y position with each row. Along a row we step along x in each column.
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase, IFMprob, Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam

	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = offset + position					// account for y offset
	Calc_Efield(position_shifted)							// calculate E-field at each data point
	
	CalcPhiOverPol()										// calculate phase/polarizability
	
	Phi(paramWave, predictedPhase, position)
	
	//CleanupPolWaves()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End



Function PredictPhase2IFM(alpha, offset)
	Variable alpha
	Variable offset
	
	NVAR velocity, sigma
	
	Variable runtimer = 0
	If(runtimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Make/O/N=1 paramWave = {alpha}
	
	Make/O/N=20 position = .1 * x/1000+.0001, predictedPhase
	
	Variable/G num_pos = numpnts(position)
	Make/o/n=(num_pos,num_x, number_v) yy, xx   						// yy: along a column we increase in y position with each row. Along a row we step along x in each column.
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase, IFMprob, Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam

	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = offset + position					// account for y offset
	Calc_Efield2IFM(position_shifted)							// calculate E-field at each data point
	
	CalcPhiOverPol2IFM()										// calculate phase/polarizability
	
	Phi2IFM(paramWave, predictedPhase, position)
	
	//CleanupPolWaves()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End




Function PredictPhase2Elec(alpha, offset)
	Variable alpha
	Variable offset
	
	NVAR velocity, sigma
	
	Variable runtimer = 0
	If(runtimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Make/O/N=1 paramWave = {alpha}
	
	Make/O/N=20 position = .1 * x/1000, predictedPhase
	
	Variable/G num_pos = numpnts(position)
	Make/o/n=(num_pos,num_x, number_v) yy, xx   						// yy: along a column we increase in y position with each row. Along a row we step along x in each column.
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase, IFMprob, Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam

	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = offset + position					// account for y offset
	Calc_Efield(position_shifted)							// calculate E-field at each data point
	
	CalcPhiOverPol()										// calculate phase/polarizability
	
	Phi(paramWave, predictedPhase, position)
	
	CleanupPolWaves()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End



Function PredictPhase_Isotopes(alpha, offset)
	Variable alpha
	Variable offset
	
	NVAR velocity, sigma, mass
	
	Variable runtimer = 0
	If(runtimer)
		Variable timerRefNum = Startmstimer
	EndIf
	
	Make/O/N=1 paramWave = {alpha}
	
	Make/O/N=20 position = .1 * x/1000+.0001, predictedPhase
	
	Variable/G num_pos = numpnts(position)
	Make/o/n=(num_pos,num_x, number_v) yy, xx   						// yy: along a column we increase in y position with each row. Along a row we step along x in each column.
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase, IFMprob, Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam		
	
	CalculateConstants_Isotopes()
	
	NormalizeVelocities(velocity, sigma, n)
		
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = offset + position					// account for y offset
	Calc_Efield_Isotopes(position_shifted)							// calculate E-field at each data point
		
	CalcPhiOverPol_Isotopes()										// calculate phase/polarizability
		
	Phi_Isotopes(paramWave, predictedPhase, position)
	
	CleanupPolWaves()
	
	If(runtimer)
		Print "New: ", StopMSTimer(timerRefNum)*10^-6
	EndIf
End








//		PhaseFit(pos_polfit_pad_red, phase_polfit_red, phase_error_polfit_red)
//PhaseFit(pos_b_polfit_pad_red, phase_b_polfit_red, phase_error_b_polfit_red)
Function PhaseFitFunc(SeriesName)
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
	Make/o/n=(num_pos,num_x, number_v) yy, xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)						// The x distance squared
	Make/o/n=(num_pos) IFMphase, IFMprob, Contrast		// Waves for: total phase shift of monochrom. beam; phase shift of polychromatic beam
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = PositionWave + shift		// account for y offset
	Calc_Efield(position_shifted)							// calculate E-field at each data point

	CalcPhiOverPol()										// calculate phase/polarizability
	
	Variable LastFitPointLoc = NumVarOrDefault("LastFitPoint", 1000)
	
	Make/D/O W_coef={24} 				// Initial Guess
	Make/D/O Epsilonwave={1e-1}			// Amount to vary the fit parameter by in each iteration	
	FuncFit/N/M=2/TBOX=768 Phi W_coef PhaseWave[0,LastFitPointLoc]  /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave// /M=mask_g
	
	AppendFittedPolWave(PositionWave, PhaseWave)				

	CleanupPolWaves()

	if(runTimer)
		Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
	endIf
End



Function PhaseFit2IFM(SeriesName)
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
	Make/o/n=(num_pos,num_x, number_v) yy, xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)						// The x distance squared
	Make/o/n=(num_pos) IFMphase, IFMprob, Contrast		// Waves for: total phase shift of monochrom. beam; phase shift of polychromatic beam
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Variable shift = offset									// MatrixOP doesn't play well with constants, so convert to a local variable
	MatrixOP/O position_shifted = PositionWave + shift		// account for y offset
	Calc_Efield2IFM(position_shifted)							// calculate E-field at each data point

	CalcPhiOverPol2IFM()										// calculate phase/polarizability
	
	Variable LastFitPointLoc = NumVarOrDefault("LastFitPoint", 1000)
	
	Make/D/O W_coef={24} 				// Initial Guess
	Make/D/O Epsilonwave={1e-1}			// Amount to vary the fit parameter by in each iteration	
	FuncFit/N/M=2/TBOX=768 Phi2IFM W_coef PhaseWave[0,LastFitPointLoc]  /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave// /M=mask_g
	
	AppendFittedPolWave2IFM(PositionWave, PhaseWave)				

	CleanupPolWaves()

	if(runTimer)
		Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
	endIf
End





Function PhaseFitOffset(PositionWave, PhaseWave, ErrorWave)
	Wave PositionWave, PhaseWave, ErrorWave
	
	NVAR velocity, sigma
	
	Variable timerRefNum = Startmstimer
	
	Variable/G num_pos = numpnts(PositionWave)
	Make/o/n=(num_pos,num_x, number_v) yy, xx   
	xx[][][] = ((y * x_step) - (.5*num_x * x_step))
	MatrixOP/O xxsqr = magsqr(xx)	
	Make/o/n=(num_pos) IFMphase, IFMprob, Contrast						//Total phase shift of monochrom. beam; phase shift of polychromatic beam
	
	CalculateConstants()
	
	NormalizeVelocities(velocity, sigma, n)
	
	Make/D/O W_coef={24, 48.7e-6}	// initial guesses {alpha, shift}
	Make/D/O Epsilonwave={1e-1,1e-6}	// Amount to vary the fit parameters by in each iteration
	FuncFit/N/M=2/TBOX=768 PhiOffset W_coef PhaseWave /X=PositionWave /W=ErrorWave /I=1 /R /E=EpsilonWave 
	
	AppendFittedPolWave(PositionWave, PhaseWave)	
	
	CleanupPolWaves()
	
	//Print "New fit: ", StopMSTimer(timerRefNum)*10^-6
End




Function CalculateConstants()
//	Variable/G Vwire = (6035 * HVmon)	 +0.7							// Applied electrode voltage [V]
	
	NVAR HVmon, ZeroMon, mass
	
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

Function CalculateConstants_Isotopes()
//	Variable/G Vwire = (6035 * HVmon)	 +0.7							// Applied electrode voltage [V]
	Wave MassWave
	NVAR HVmon, ZeroMon
	
	Variable/G Vwire = GetHVKeithley(HVmon, ZeroMon)*1000
	Variable/G d = sqrt((R+gap)^2 - R^2)								// Distance from ground plane to the "infinite wire" [m]
	Variable/G lambda = Vwire*pi*eps0 / asinh(d/R)				// Charge density of the "infinite wire"
	Variable/G E0 = -lambda / (pi*eps0)		
	Duplicate/O MassWave SeparationWave
	SeparationWave = (L*hbar*2*pi)/ masswave/a				// beam separationS at center of grad e region
End




Function NormalizeVelocities(velocity, sigma, n)			//Find the normalization for the density profile. 
	Variable velocity, sigma, n

	NVAR num_pos
	Make/D/o/n=(number_v) vel, dens
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
		
	
	Make/D/o/n=(num_pos,num_x, number_v) v
	v[][][] =( (velocity) - (.5 * number_v * v_step ))+(v_step * z) // all rows and columns of a given layer have the same velocity
	MatrixOP/O density = powR(v,n)* exp (-magsqr(v-(velocity))/(2*magsqr(sigma))) / vel_normalization	//Probability density distribution. same layering as v
	MatrixOP/O vrecip = rec(v)
	
	Make/D/O/N=(num_pos, number_v) density2D = dens[q]
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
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqr0 = magsqr(Etot)
	
	MatrixOP/O DplusYY = d + yy - (s * vrecip)
	MatrixOP/O DminusYY =  s*vrecip - yy + d 
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrM = magsqr(Etot)

	s*=-1
	MatrixOP/O DplusYY = d + yy - (s * vrecip)
	MatrixOP/O DminusYY =  s*vrecip - yy + d 
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)
	MatrixOP/O EtotMagSqrP = magsqr(Etot)
End


Function Calc_EField_Isotopes(positionwave)			// Calculates the electric field of the gradient region and makes waves for phase, prob, and contrast
	Wave positionwave						// positions at which to calculate the E-field
	Nvar num_pos, Vwire, d, lambda, E0
	
	Wave yy, xx, xxsqr, vrecip, SeparationWave
	yy[][][] = positionwave[p]				// populate the yy wave with the positions. Each row contains a different position[p], but the position is the same across columns and layers.
	
	// First we calculate E^2 for the undiffracted beam.
	MatrixOP/O DplusYY = d+yy			
	MatrixOP/O DminusYY = -1*yy + d

	// The electric field magnitude of the primary and diffracted beams
	MatrixOP/O Etot= (E0) * sqrt(		magsqr(DplusYY/(magsqr(DplusYY)+xxsqr) + DminusYY/(magsqr(DminusYY)+xxsqr)) + magsqr(xx/(magsqr(DplusYY)+xxsqr) - xx/(magsqr(DminusYY)+xxsqr))   		)

	// Need the magnitude squared.
	MatrixOP/O EtotMagSqr = magsqr(Etot)
	
	String EsqrdWaveName
	Make/O/WAVE/N=(numpnts(SeparationWave)) EsqrdSepWaveRefs
	Variable i = 0, s
	do
		// For the diffracted beam now:
		s = SeparationWave[i]
		MatrixOP/O DplusYsep = d + yy - (s * vrecip)
		MatrixOP/O DminusYsep =  s*vrecip - yy + d 
	
		MatrixOP/O Etot_s= (E0)*sqrt(	magsqr(DplusYsep/(magsqr(DplusYsep)+xxsqr) + DminusYSep/(magsqr(DminusYSep)+xxsqr)) + magsqr(xx/(magsqr(DplusYsep)+xxsqr) - xx/(magsqr(DminusYsep)+xxsqr))		)

		MatrixOP/O Etot_sMagSqr = magsqr(Etot_s)
		
		EsqrdWaveName = "Etot_sMagSqrM" + num2str(i)
		Duplicate/O Etot_sMagSqr $EsqrdWaveName///WAVE=w
		EsqrdSepWaveRefs[i] = $EsqrdWaveName
		
		i+=1
	while(i<numpnts(SeparationWave))
	
	KillWaves DplusYY, DminusYY, Etot, DplusYsep, DminusYsep, Etot_s, Etot_sMagSqr
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
	Wave EtotMagSqr0, EtotMagSqrM, EtotMagSqrP, vrecip, IFMphase 
	// Find the phase shift for each velocity component of the main and diffracted beam
	MatrixOP/O dphi = EtotMagSqr0 * vrecip							//phase shift of first path:   alpha*E^2/v	// this is the +1 diffraction order
	MatrixOP/O dphi_sM = EtotMagSqrM * vrecip					//phase shift of second path a distance "sep" away		// this is the 0th order
	MatrixOP/O dphi_sP = EtotMagSqrP * vrecip				//phase shift of third path a distance "sep" away		// this is the -1 order
	
	Integrate/dim=1 dphi, dphi_sM, dphi_sP
	MatrixOP/O IFMphase_v_overpolM = dphi - dphi_sM				// phase shift divided by polarizability for minus ifm
	MatrixOP/O IFMphase_v_overpolP = dphi_sP - dphi			// phase shift divided by polarizability for plus ifm
End



Function CalcPhiOverPol_Isotopes()
	Wave EtotMagSqr, vrecip, IFMphase, SeparationWave
	// Find the phase shift for each velocity component of the main and diffracted beam
	MatrixOP/O dphi = EtotMagSqr * vrecip							//phase shift of first path:   alpha*E^2/v	
	Integrate/dim=1 dphi 
		
	Wave/Wave EsqrdSepWaveRefs
	//edit esqrdsepwaverefs[1]
	String IFMphase_v_overpolWaveName
	//Duplicate/O EsqrdSepWaveRefs IFMphase_v_overpol_Refs
	Make/O/WAVE/N=(numpnts(EsqrdSepWaveRefs)) IFMphase_v_overpol_Refs
	Variable i = 0
	//Wave LocalEsqrd=EsqrdSepWaveRefs[0]
	//edit localesqrd
	do
		// For the diffracted beam now:
		Wave LocalEsqrd =  EsqrdSepWaveRefs[i]
		MatrixOP/O dphi_s = LocalEsqrd* vrecip
		Integrate/dim=1 dphi_s
		MatrixOP/O IFMphase_v_overpol = dphi-dphi_s
		
		IFMphase_v_overpolWaveName = "IFMphase_v_overpolM" + num2str(i)
		Duplicate/O IFMphase_v_overpol $IFMphase_v_overpolWaveName
		IFMphase_v_overpol_Refs[i] = $IFMphase_v_overpolWaveName
		
		i+=1
	while(i<numpnts(SeparationWave))	
	
	KillWaves dphi, dphi_s
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





Function Phi_Isotopes(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
	
	Wave IFMprob, density, Contrast//, IFMphase_v_overpol39, IFMphase_v_overpol40, IFMphase_v_overpol41
	
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	Variable loc_v_step = v_step		
	
	Wave MassFracWave
	
	Wave/Wave IFMphase_v_overpol_Refs
	String TempWaveNameImag, TempWaveNamePrevImag, TempWaveNameReal, TempWaveNamePrevReal
	//Duplicate/O EsqrdSepWaveRefs IFMphase_v_overpol_Refs
	//Make/O/WAVE/N=(numpnts(EsqrdSepWaveRefs)) IFMphase_v_overpol_Refs

	//Wave LocalEsqrd=EsqrdSepWaveRefs[0]
	//edit localesqrd
	//KillWaves/z imagpart, realpart
	Wave LocalifmPhase =  IFMphase_v_overpol_Refs[0]
	MatrixOP/O imagpartM0 = MassFracWave[0] * sin(pol * LocalifmPhase)
	MatrixOP/O realpartM0 = MassFracWave[0] * cos(pol * LocalifmPhase)
	
	
		Variable i = 1
		do
			Wave LocalifmPhase =  IFMphase_v_overpol_Refs[i]
			
			TempWaveNameImag = "imagpartM"+num2str(i)
			TempWaveNamePrevImag = "imagpartM"+num2str(i-1)
			Wave TempPrevImag = $TempWaveNamePrevImag
			MatrixOP/O $TempWaveNameImag = MassFracWave[i] * sin(pol*LocalifmPhase)+TempPrevImag
			
			TempWaveNameReal = "realpartM"+num2str(i)
			TempWaveNamePrevReal = "realpartM"+num2str(i-1)
			Wave TempPrevReal = $TempWaveNamePrevReal
			MatrixOP/O $TempWaveNameReal = MassFracWave[i] * cos(pol*LocalifmPhase)+TempPrevReal
		
		
			//IFMphaseWaveName = "IFMphaseM" + num2str(i)
			//Duplicate/O IFMphase $IFMphaseWaveName
			//IFMphase_v_overpol_Refs[i] = $IFMphaseWaveName
		
			i+=1
		while(i<numpnts(IFMphase_v_overpol_Refs))	
	
	Wave TempImag = $TempWaveNameImag
	MatrixOP/O imagpart = TempImag*loc_v_step*density
	Wave TempReal = $TempWaveNameReal
	MatrixOP/O realpart = TempReal*loc_v_step*density
	Integrate/DIM=2 imagpart, realpart
	
	//MatrixOP/O IFMphase_v39 = pol * IFMphase_v_overpol39			// now we have the actual phase shift to take the real and imaginary parts of
	//MatrixOP/O IFMphase_v40 = pol * IFMphase_v_overpol40	
	//MatrixOP/O IFMphase_v41 = pol * IFMphase_v_overpol41	
	
	
	// We need to find the total velocity weighted phase shift
	// phi_total = ArcTan( Int( Prob(v) * exp(i * delta_phi) ) ) 
	//			 = ArcTan( Int(Prob(v)*sin(phi)) / Int(Prob(v)*cos(phi)) )
	// the code to implement this follows:
	
//	Variable/G p39 = K39frac
//	Variable/G p40 = K40frac
//	Variable/G p41 = K41frac
//	
//	Variable loc_v_step = v_step											
//	MatrixOP/O imagpart = loc_v_step * density * (p39*sin(IFMphase_v39)+p40*sin(IFMphase_v40)+p41*sin(IFMphase_v41))			// weight imaginary part by prob(v) 
//	MatrixOP/O realpart = loc_v_step * density * (p39*cos(IFMphase_v39)+p40*cos(IFMphase_v40)+p41*cos(IFMphase_v41))				// weight real part by prob(v)
//	Integrate/DIM=2 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift

	Contrast = sqrt( magsqr(imagpart[p][num_x-1][number_v-1]) + magsqr(realpart[p][num_x-1][number_v-1]) )
	IFMprob = atan2(imagpart[p][num_x-1][number_v-1],realpart[p][num_x-1][number_v-1])
	MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
	if(IFMprob[0]<0)
		IFMprob+=2*pi
	endif
	
	output = IFMprob
End



Function Phi2IFM(pw, output, y) :FitFunc
	Wave pw, output, y									// pw = fit parameters; output = predicted phase shift; y = grad E position relative to beam
	
	Wave IFMprob, density, Contrast, IFMphase_v_overpolM, IFMphase_v_overpolP
	
	Variable pol = pw[0]*4*pi*eps0*1E-30* x_step / (2 * hbar)	// pw[0] is alpha in cgs units but without the factor of 10^-24
	
	MatrixOP/O IFMphase_vM = pol * IFMphase_v_overpolM			// now we have the actual phase shift to take the real and imaginary parts of
	MatrixOP/O IFMphase_vP = pol * IFMphase_v_overpolP

	Variable loc_v_step = v_step											
	MatrixOP/O imagpart = 0.5*loc_v_step * density * (sin(IFMphase_vM) + sin(IFMphase_vP) )			// weight imaginary part by prob(v) 
	MatrixOP/O realpart = 0.5*loc_v_step * density * (cos(IFMphase_vM) + cos(IFMphase_vP) )				// weight real part by prob(v)
	Integrate/DIM=2 imagpart, realpart			// integrate real and imaginary parts along dv to find prob(v) weighted phase shift
	
	Contrast = sqrt( magsqr(imagpart[p][num_x-1][number_v-1]) + magsqr(realpart[p][num_x-1][number_v-1]) )
	IFMprob = atan2(imagpart[p][num_x-1][number_v-1],realpart[p][num_x-1][number_v-1])
	MatrixOP/O PhaseWrapped = IFMprob									// Duplicates IFMprob before unwrap for debugging
	
	Unwrap 2*Pi, IFMprob		//Correction from taking the arctan.  
	if(IFMprob[0]<0)
		IFMprob+=2*pi
	endif
	
	output = IFMprob
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





Function CleanupPolWaves()
	KillWaves/Z dens, density, DminusYsep, DminusYY, dphi, dphi_s, DplusYsep, DplusYY, Etot, EtotMagSqr, Etot_s, Etot_sMagSqr, v
	KillWaves/Z imagpart, realpart, xx, xxsqr, yy
	KillWaves/Z dphi_sM, dphi_sP, EtotMagSqr0, EtotMagSqrM, EtotMagSqrP, IFMphase, IFMphase_v, IFMphase_vM, IFMPhase_vP, IFMphase_v_overpol, IFMphase_v_overpolM, IFMphase_v_overpolP, vel, vrecip
	KillWaves/Z EtotMagSqrMmol, EtotMagSqrPmol, dphiMmol, dphiPmol, IFMphase_v_overpolMmol, IFMphase_v_overpolPmol, IFMphase_vMmol, IFMphase_vPmol
	KillWaves/z phiSag, phiSagAvg, phiSagLayered
	KillWaves/Z density2D, phiGrav, phiSagGrav, phiSagGravImag, phiSagGravLayered2D, phiSagGravReal, phiSagGravWeighted, phiSagImag, phiSagLayered2D, phiSagReal
	KillWaves/Z EtotMagSqrM0, EtotMagSqrMmol0, EtotMagSqrP0, EtotMagSqrPmol0, EtotMagSqrM1, EtotMagSqrMmol1, EtotMagSqrP1, EtotMagSqrPmol1,EtotMagSqrM2, EtotMagSqrMmol2, EtotMagSqrP2, EtotMagSqrPmol2
	KillWaves/Z  IFMphase_v_overpolM0, IFMphase_v_overpolP0,IFMphase_v_overpolMmol0, IFMphase_v_overpolPmol0, IFMphase_v_overpolM1, IFMphase_v_overpolP1,IFMphase_v_overpolMmol1, IFMphase_v_overpolPmol1, IFMphase_v_overpolM2, IFMphase_v_overpolP2,IFMphase_v_overpolMmol2, IFMphase_v_overpolPmol2
	KillWaves/Z EsqrdSepWaveRefs, IFMphase_v_overpol_Refs
	KillWaves/Z dphi_sM2, dphi_sP2, dphiM2mol, dphiP2mol, IFMphase_v_overpolM2, IFMphase_v_overpolP2, IFMphase_v_overpolM2mol, IFMphase_v_overpolP2mol
	KillWaves/Z EtotMagSqrM2, EtotMagSqrP2, EtotMagSqrM2mol, EtotMagSqrP2mol
	KillWaves/Z IFMphase_vM2, IFMphase_vP2, IFMphase_v_overpolM2, IFMphase_v_overpolP2
End




// Reduces the raw position, phase, and error waves to just the important parts
Function ReducePolWaves(position,phase,error)			
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
	SetAxis left 0,*
	SetAxis bottom 0,*
	if(1)	// set to true for finer ticks and measure distance in mm
		ModifyGraph nticks(bottom)=10,minor(bottom)=1,prescaleExp(bottom)=3
		Label bottom "Distance to Ground Plane (mm)"
		ModifyGraph minor(left)=1
		ModifyGraph sep(left)=4
	endif
	
	String DataFolderNameInit = GetDataFolder(0)
	String DataFolderName = ReplaceString("'", DataFolderNameInit, "")
	TextBox/C/N=text0/A=MC DataFolderName
	TextBox/C/N=text0/X=-16.00/Y=21.00
	WindowNamer(DataFolderName+" Phase vs Position")
//	String WindowName = "PP"+DataFolderName
//	DoWindow/C/T $WindowName, DataFolderName+": Phase vs Position"
	
	Print PositionReducedWaveName
	Print PhaseReducedWaveName
	Print ErrorReducedWaveName
End



Function ReducePhaseContrastWaves(position,phase,error, contrast, contrastError)			
	Wave position, phase, error, contrast, contrastError
	
	// Make new wave names
	String PositionReducedWaveName = NameOfWave(position) +"_red"
	String PhaseReducedWaveName= NameOfWave(phase) + "_red"
	String ErrorReducedWaveName= NameOfWave(error) + "_red"
	String contrastReducedWaveName= NameOfWave(contrast) + "_red"
	String ContrastErrorReducedWaveName= NameOfWave(contrasterror) + "_red"
	
	// The reduced waves are initially just copies of the raw waves
	Duplicate/O position $PositionReducedWaveName; Wave PositionReduced = $PositionReducedWaveName
	Duplicate/O phase $PhaseReducedWaveName; Wave PhaseReduced = $PhaseReducedWaveName
	Duplicate/O error $ErrorReducedWaveName; Wave ErrorReduced = $ErrorReducedWaveName
	Duplicate/O contrast $contrastReducedWaveName; Wave contrastReduced = $contrastReducedWaveName
	Duplicate/O Contrasterror $ContrastErrorReducedWaveName; Wave ContrastErrorReduced = $ContrastErrorReducedWaveName
	
	// If a phase, position, or error point is a NaN then we remove it and the corresponding points in the other waves
	RemoveNaNsXYZAB(PhaseReduced, PositionReduced, ErrorReduced, ContrastReduced, ContrastErrorReduced)
	
	// Sort the waves so that point 0 is closest to the ground plane (makes unwrapping easier)
	Sort PositionReduced,PositionReduced,PhaseReduced,ErrorReduced, ContrastReduced, ContrastErrorReduced
	
	Display/K=1/W=(40,100,940,700) PhaseReduced vs PositionReduced
	AppendToGraph/L=contrast ContrastReduced vs PositionReduced
	ModifyGraph axisEnab(left)={.40,1.00}
	ModifyGraph axisEnab(contrast)={0,.35}
	ModifyGraph freePos(contrast)=0
	ModifyGraph lblPosMode(left)=1,lblPosMode(contrast)=1,lblMargin(left)=10;DelayUpdate
	ModifyGraph lblMargin(contrast)=10;DelayUpdate
	Label contrast "Contrast"
	ModifyGraph mode=3,marker=8
	ErrorBars $PhaseReducedWaveName Y,wave=(ErrorReduced,ErrorReduced)
	ErrorBars $ContrastReducedWaveName Y,wave=(ContrastErrorReduced,ContrastErrorReduced)
	Label left "Phase"; Label bottom "Distance to Ground Plane (m)"
	ModifyGraph zero(left)=1
	//SetAxis left 0,*
	//SetAxis bottom 0,*
	ModifyGraph nticks(contrast)=10;DelayUpdate
	SetAxis contrast 0,*
	ModifyGraph gmSize=4
	if(1)	// set to true for finer ticks and measure distance in mm
		ModifyGraph nticks(bottom)=10,minor(bottom)=1,prescaleExp(bottom)=3
		Label bottom "Distance to Ground Plane (mm)"
		ModifyGraph minor(left)=1
		ModifyGraph sep(left)=4
	endif
	ModifyGraph gfSize=11
	
	String DataFolderNameInit = GetDataFolder(0)
	String DataFolderName = ReplaceString("'", DataFolderNameInit, "")
	TextBox/C/N=text0/A=MC DataFolderName
	TextBox/C/N=text0/X=-16.00/Y=21.00
	WindowNamer(DataFolderName+" Phase and Contrast vs Position")
//	String WindowName = "PP"+DataFolderName
//	DoWindow/C/T $WindowName, DataFolderName+": Phase vs Position"
	
	Print PositionReducedWaveName
	Print PhaseReducedWaveName
	Print ErrorReducedWaveName
	Print ContrastReducedWaveName
	Print ContrastErrorReducedWaveName
End








Function AppendFittedPolWave(PositionWave, PhaseWave)
	Wave PositionWave, PhaseWave
	Wave W_coef
	
	NVAR offset
	
	If(numpnts(W_coef)==1)						// Determine if the offset was fixed or a fit parameter
		PredictPhase(W_coef[0], offset)			// Run the predict phase routine with the fitted polarizability and fixed offset
	else
		PredictPhase(W_coef[0],W_coef[1])		// Run the predict phase routine with the fitted polarizaiblity and offset
	EndIf
		
	Wave position, predictedPhase
		
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
	Wave W_coef
	
	NVAR offset
	
	PredictPhase2IFM(W_coef[0], offset)			// Run the predict phase routine with the fitted polarizability and fixed offset
	
	Wave position, predictedPhase
		
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
	//CurveFitDialog/ f(x) = 0.5*amplitude*(1+Erf((x-x0)/sqrt(2*sigma^2)))+background
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = amplitude
	//CurveFitDialog/ w[1] = x0
	//CurveFitDialog/ w[2] = sigma
	//CurveFitDialog/ w[3] = background

	return 0.5*w[0]*(1+Erf((x-w[1])/sqrt(2*w[2]^2)))+w[3]
End



Function errfitneg(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = -0.5*amplitude*(Erf((x-x0)/sqrt(2*sigma^2))-1)+background
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = amplitude
	//CurveFitDialog/ w[1] = x0
	//CurveFitDialog/ w[2] = sigma
	//CurveFitDialog/ w[3] = background

	return -0.5*w[0]*(Erf((x-w[1])/sqrt(2*w[2]^2))-1)+w[3]
End




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



Function ExtractEclipseDataV3(SeriesName)
	String SeriesName		// ex: "e"     expects to find wave names with the convention in the next 5 lines
	
	
	String initposWaveName="pos_"+SeriesName
	String initcontrastWaveName = "contrast_"+SeriesName
	String initcontrast_errorWaveName = "contrast_error_"+SeriesName
	String initcountsWaveName = "avg_counts_"+SeriesName
	String initcounts_errorWaveName = "avg_counts_error_"+SeriesName
	
	Wave poswave = $initposWaveName
	Wave contrast = $initcontrastWaveName
	Wave contrast_error = $initcontrast_errorWaveName
	Wave counts = $initcountsWaveName
	Wave counts_error = $initcounts_errorWaveName

	NVAR LastEclipse
	Variable maxFile = LastEclipse		+15	// Last file with eclipse data plus 15 more to get through the first ref phase
	
	
	Variable lastEclipseIndex = maxFile/5
	Variable firstRefPhaseIndex = lastEclipseIndex+3
	String EclipsePosWaveName = NameOfWave(poswave) + "_eclipse"
	Make/O/N=(firstRefPhaseIndex) $EclipsePosWaveName
	Wave EclipsePos = $EclipsePosWaveName
	EclipsePos = poswave
	
	String EclipsePosPadWaveName = NameOfWave(poswave) + "_eclipse_pad"
	Make/O/N=(maxFile) $EclipsePosPadWaveName
	Wave EclipsePosPad = $EclipsePosPadWaveName
	Variable n=0; Variable i=0
	do
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
	String EclipseCCWaveName = "CountsContrast_"+SeriesName
	String EclipseCCErrorWaveName = "CountsContrast_error_"+SeriesName
	
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
	
	DeletePoints maxFile-15, 10, EclipseContrast, EclipseContrastError, EclipseCounts, EclipseCountsError, EclipseCountsContrast, EclipseCountsContrastError, EclipsePosPad
	
	NVAR KillEclipseFiles
	Wave EclipseFilesToKill
	
	Variable NumPointsToKill, Killindex
	
	If(KillEclipseFiles)
		NumPointsToKill = numpnts(EclipseFilesToKill)
		i = 0
		do
			Killindex = EclipseFilesToKill[i]-1
			if(Killindex>LastEclipse)
				EclipseContrast[Killindex-10] = NaN					// RemoveNaNs will kill the corresponding points in Counts....
			else
				EclipseContrast[Killindex] = NaN
			endif
			i+=1
		while(i<NumPointsToKill)
	EndIf
	
	
//	Keep4(EclipseContrast)
//	Keep4(EclipseContrastError)
//	Keep4(EclipseCounts)
//	Keep4(EclipseCountsError)
//	Keep4(EclipseCountsContrast)
//	Keep4(EclipseCountsContrastError)
//	Keep4(EclipsePosPad)
	
	Keep3(EclipseContrast)
	Keep3(EclipseContrastError)
	Keep3(EclipseCounts)
	Keep3(EclipseCountsError)
	Keep3(EclipseCountsContrast)
	Keep3(EclipseCountsContrastError)
	Keep3(EclipsePosPad)
	
	RemoveNaNsXYZABCD(EclipseContrast, EclipseContrastError,EclipseCounts,EclipseCountsError,EclipseCountsContrast,EclipseCountsContrastError,EclipsePosPad)
	
	if(1)
	//Normalize the Counts, Contrast, and C*C waves
	Variable pnts = numpnts(EclipseContrast)
	WaveStats/R=[pnts-1,pnts-4] EclipseContrast
	EclipseContrast /= V_avg; EclipseContrastError /= V_avg
	WaveStats/R=[pnts-1,pnts-4] EclipseCounts
	EclipseCounts /= V_avg; EclipseCountsError /= V_avg
	WaveStats/R=[pnts-1,pnts-4] EclipseCountsContrast
	EclipseCountsContrast /= V_avg; EclipseCountsContrastError /= V_avg
	endif
	
	//Make a fancy graph
	Display/K=1/W=(40,100,850,600) EclipseContrast vs EclipsePosPad
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
	
	String DataFolderNameInit = GetDataFolder(0)
	String DataFolderName = ReplaceString("'", DataFolderNameInit, "")
	TextBox/C/N=text0/A=MC DataFolderName
	String WindowName = DataFolderName+" Eclipse Data"
	WindowNamer(WindowName)
	
	//String WindowName = "Eclipse"+DataFolderName
	//DoWindow/C/T $WindowName, DataFolderName+": Phase vs Position"
	
//	Display/K=1/W=(40,350,440,575) EclipseCountsContrast vs EclipsePosPad
//	ModifyGraph mode=3,marker=8
//	ModifyGraph rgb($EclipseCCWaveName)=(1,16019,65535)
//	ErrorBars $EclipseCCWaveName Y,wave=(EclipseCountsContrastError,EclipseCountsContrastError)
//	Label left "Counts*Contrast"; Label bottom "Position"
	
	
	String tboxname = "CF_"+EclipseCCWaveName
	TextBox/C/N=$tboxname/X=(-5)/Y=50
	
	//Fit the counts*contrast data and report the 50% point
	Make/O/D/N=4 W_coef
	W_coef[0] = {2,40000,.0019200,.01}
	FuncFit/NTHR=1/TBOX=768 errfitneg W_coef  EclipseCountsContrast /X=EclipsePosPad /W=EclipseCountsContrastError /I=1 /D 
	
	Variable polfitCCoffset = gap-w_coef[2]
	String offsetboxtext = "Offset: " + num2str(polfitCCoffset) + " m"
	TextBox/C/N=CCoffsetbox/X=10/Y=70 offsetboxtext
	
	print "counts*contrast fit offset: " , polfitCCoffset
	
	
	//Fit the counts data only and report the 50% point
	Make/O/D/N=4 W_coef
	W_coef[0] = {2,40000,.0019200,.01}
	FuncFit/NTHR=1/TBOX=768 errfitneg W_coef  EclipseCounts /X=EclipsePosPad /W=EclipseCountsError /I=1 /D 
	
	Variable polfitCountsoffset = gap-w_coef[2]
	offsetboxtext = "Offset: " + num2str(polfitCountsoffset) + " m"
	TextBox/C/N=Countsoffsetbox/X=10/Y=30 offsetboxtext
	
	print "counts fit offset: " , polfitCountsoffset
	
	print "counts only - counts*contrast: ", (polfitCountsoffset- polfitCCoffset)*-1
End






//			ExtractEclipseDataV2(pos, contrast, contrast_error, avg_counts, avg_counts_error)
Function ExtractEclipseDataV2(poswave, contrast, contrast_error, counts, counts_error)
	Wave poswave, contrast, contrast_error, counts, counts_error

	NVAR LastEclipse
	Variable maxFile = LastEclipse			// Last file with eclipse data
	
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
	String EclipseCCWaveName = "CountsContrast_"
	String EclipseCCErrorWaveName = "CountsContrast_error_"
	
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
	
	NVAR KillEclipseFiles
	Wave EclipseFilesToKill
	
	Variable NumPointsToKill, Killindex
	
	If(KillEclipseFiles)
		NumPointsToKill = numpnts(EclipseFilesToKill)
		i = 0
		do
			Killindex = EclipseFilesToKill[i]-1
			EclipseContrast[Killindex] = NaN					// RemoveNaNs will kill the corresponding points in Counts....
			i+=1
		while(i<NumPointsToKill)
	EndIf
	
	
//	Keep4(EclipseContrast)
//	Keep4(EclipseContrastError)
//	Keep4(EclipseCounts)
//	Keep4(EclipseCountsError)
//	Keep4(EclipseCountsContrast)
//	Keep4(EclipseCountsContrastError)
//	Keep4(EclipsePosPad)
	
	Keep3(EclipseContrast)
	Keep3(EclipseContrastError)
	Keep3(EclipseCounts)
	Keep3(EclipseCountsError)
	Keep3(EclipseCountsContrast)
	Keep3(EclipseCountsContrastError)
	Keep3(EclipsePosPad)
	
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
	String DataFolderName = GetDataFolder(0)
	TextBox/C/N=text0/A=MC DataFolderName
	
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
// be careful about file vs index. 
// Use maxFile = 0 to automatically use the last file
// example: 		ExtractPhaseDataV2(pos, phase, phase_error)
Function ExtractPhaseDataV3(SeriesName)
	String SeriesName		// ex: "e"
	
	String initposWaveName="pos_"+SeriesName
	String initphaseWaveName = "phase_"+SeriesName
	String initphase_errorWaveName = "phase_error_"+SeriesName
	
	Wave poswave = $initposWaveName
	Wave phase = $initphaseWaveName
	Wave phase_error = $initphase_errorWaveName
	
	NVAR LastEclipse, LastFile
	Variable minFile = LastEclipse+1
	Variable maxFile = LastFile		
	
	Make/O/N=(maxFile-minFile+5+1) FileNumbers = LastEclipse-5+p+1
	
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
	String RefPhaseResidWN = "Res_"+NameOfWave(PolFitPhase)+"_ref"
	String FitTextBoxName = "CF_"+NameOfWave(PolFitPhase)+"_ref"
	
	Duplicate/O PolFitPhase $RefPhaseWaveName				//copy the phase data into a wave that will eventually only hold reference phase data
	Duplicate/O PolFitPhaseError $RefPhaseErrorWaveName
	
	Wave RefPhase = $RefPhaseWaveName
	Wave RefPhaseError = $RefPhaseErrorWaveName
	
		
		
	//	extract ref phase
	if(1)
		n=0; i=0
		do
			RefPhase[i]=NaN				//This is a HV phase, so erase it from the ref phase wave
			RefPhaseError[i]=NaN
			
			if(n==9)						// if n=9 we've gone through 10 HV phase measurements, so the next 5 will be reference phases
				n=0						
				i+=5						// so we'll skip the eraser over them in the next iteration
			else
				n+=1
			endif
			
			i+=1
		while(i<numpnts(RefPhase))
		
		i=0
		do										//remove reference phase measurements from phase data
			if(PolFitPhase[i]==RefPhase[i])	// returns true if the point is a reference phase measurement
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
	

	//Display/K=1/W=(40,100,550,600) 
	Display/K=1/W=(350,150,1250,650) PolFitPhase vs FileNumbers
	ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
	Label left "Phase"; Label bottom "File number"
	AppendToGraph/L='RefPhase' RefPhase vs FileNumbers
	ModifyGraph mode=3,marker=8
	ErrorBars $RefPhaseWaveName Y,wave=(RefPhaseError,RefPhaseError)
	Label RefPhase "RefPhase"
	ModifyGraph axisEnab(left)={.6,1}
	ModifyGraph axisEnab(RefPhase)={0,0.57}
	ModifyGraph freePos(RefPhase)={0,bottom}
	ModifyGraph lblPosMode(refphase)=1
	String DataFolderNameInit = GetDataFolder(0)
	String DataFolderName = ReplaceString("'", DataFolderNameInit, "")
	TextBox/C/N=text0/A=MC DataFolderName
	TextBox/C/N=text0/X=9.58/Y=45.65
	WindowNamer(DataFolderName+" Phase and Ref Phase")
	
	
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
	

	Variable maxUnwrapIndex, NumPointsToKill, Killindex, NumPointsToUnwrap, UnwrapIndex

	NVAR PlusMinusRef, UnwrapRef1, UnwrapRef1End, KillRefFiles, UnwrapRefSpecial
	Wave RefFilesToKill, RefFilesToUnwrap
	
	if(UnwrapRef1)
		//Unwrap
		i=0
		maxUnwrapIndex = UnwrapRef1End-minFile+5+1
		do
			RefPhase[i] += plusminusRef*2*pi
			i+=1
		while(i<maxUnwrapIndex)
	endif
	
	
	If(KillRefFiles)
		NumPointsToKill = numpnts(RefFilesToKill)
		i = 0
		do
			Killindex = RefFilesToKill[i]-minFile+5
			RefPhase[Killindex] = NaN
			i+=1
		while(i<NumPointsToKill)
	EndIf
	
	If(UnwrapRefSpecial)
		NumPointsToUnwrap = numpnts(RefFilesToUnwrap)
		i = 0
		do
			UnwrapIndex = RefFilesToUnwrap[i]-minFile+5
			RefPhase[UnwrapIndex] = NaN
			i+=1
		while(i<NumPointsToUnwrap)
	EndIf
	

	NVAR plusminus, UnwrapAll, Unwrap1, Unwrap1End, Unwrap2, Unwrap2End, Unwrap3, Unwrap3End

	If(UnwrapAll)
		PolFitPhase+=2*pi
	EndIf


	if(Unwrap1)
		//Unwrap
		i=0
		maxUnwrapIndex = Unwrap1End-minFile+5+1
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<maxUnwrapIndex)
	endif
	
	
	if(Unwrap2)
		//Unwrap again
		i=0
		Variable maxUnwrapIndex2 = Unwrap2End-minFile+5+1
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<maxUnwrapIndex2)
	endif
	
	NVAR KillHVFiles
	Wave HVFilesToKill
	
	if(KillHVFiles)
		NumPointsToKill = numpnts(HVFilesToKill)
		i = 0
		do
			Killindex = HVFilesToKill[i]-minFile+5
			PolFitPhase[Killindex] = NaN
			i+=1
		while(i<NumPointsToKill)
	endif

	
	if(0)
		//Unwrap special
		PolFitPhase[13] += 2*pi
	endif





	if(1)
		//Fit the reference phase to a quadratic polynomial and adjust the phase
		Make/O/N=3 w_coef
		CurveFit/NTHR=1/TBOX=768 poly 3,  RefPhase /W=RefPhaseError /X=FileNumbers /I=1 /D /R
		TextBox/C/N=$FitTextBoxName/X=10.7/Y=56
		PolFitPhase -= w_coef[0] + w_coef[1]*FileNumbers[p] + w_coef[2]*(FileNumbers[p])^2
		Wave RefPhaseResid = $RefPhaseResidWN
		RefPhaseResid/=RefPhase
		RefPhaseResid*=RefPhase
		ModifyGraph zero(Res_RefPhase)=1
	endif

	if(0)
		Display/K=1/W=(40,100,440,325) PolFitPhase vs PosPad
		ModifyGraph mode=3,marker=8
		ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
		Label left "Phase"; Label bottom "Distance to Ground Plane (m)"
	endif
	
	
	
	NVAR InvertPhase

	if(InvertPhase)							//Make the phase shift positive
		PolfitPhase*=-1
	endif


	ReducePolWaves(PosPad,PolFitPhase, PolFitPhaseError)
End




Function ExtractContrastDataV3(SeriesName)
	String SeriesName		// ex: "e"
	
	String initposWaveName="pos_"+SeriesName
	String initphaseWaveName = "contrast_"+SeriesName
	String initphase_errorWaveName = "contrast_error_"+SeriesName
	
	Wave poswave = $initposWaveName
	Wave phase = $initphaseWaveName
	Wave phase_error = $initphase_errorWaveName
	
	NVAR LastEclipse, LastFile
	Variable minFile = LastEclipse+1
	Variable maxFile = LastFile		
	
	Make/O/N=(maxFile-minFile+5+1) FileNumbers = LastEclipse-5+p+1
	
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
	String RefPhaseResidWN = "Res_"+NameOfWave(PolFitPhase)+"_ref"
	String FitTextBoxName = "CF_"+NameOfWave(PolFitPhase)+"_ref"
	
	Duplicate/O PolFitPhase $RefPhaseWaveName				//copy the phase data into a wave that will eventually only hold reference phase data
	Duplicate/O PolFitPhaseError $RefPhaseErrorWaveName
	
	Wave RefPhase = $RefPhaseWaveName
	Wave RefPhaseError = $RefPhaseErrorWaveName
	
		
		
	//	extract ref phase
	if(1)
		n=0; i=0
		do
			RefPhase[i]=NaN				//This is a HV phase, so erase it from the ref phase wave
			RefPhaseError[i]=NaN
			
			if(n==9)						// if n=9 we've gone through 10 HV phase measurements, so the next 5 will be reference phases
				n=0						
				i+=5						// so we'll skip the eraser over them in the next iteration
			else
				n+=1
			endif
			
			i+=1
		while(i<numpnts(RefPhase))
		
		i=0
		do										//remove reference phase measurements from phase data
			if(PolFitPhase[i]==RefPhase[i])	// returns true if the point is a reference phase measurement
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
	

	//Display/K=1/W=(40,100,550,600) 
	Display/K=1/W=(350,150,1250,650) PolFitPhase vs FileNumbers
	ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
	Label left "Contrast"; Label bottom "File number"
	AppendToGraph/L='RefPhase' RefPhase vs FileNumbers
	ModifyGraph mode=3,marker=8
	ErrorBars $RefPhaseWaveName Y,wave=(RefPhaseError,RefPhaseError)
	Label RefPhase "Ref Contrast"
	ModifyGraph axisEnab(left)={.6,1}
	ModifyGraph axisEnab(RefPhase)={0,0.57}
	ModifyGraph freePos(RefPhase)={0,bottom}
	ModifyGraph lblPosMode(refphase)=1
	String DataFolderNameInit = GetDataFolder(0)
	String DataFolderName = ReplaceString("'", DataFolderNameInit, "")
	TextBox/C/N=text0/A=MC DataFolderName
	TextBox/C/N=text0/X=9.58/Y=45.65
	WindowNamer(DataFolderName+" Contrast and Ref Contrast")
	
	
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
	

	Variable maxUnwrapIndex, NumPointsToKill, Killindex, NumPointsToUnwrap, UnwrapIndex

	NVAR PlusMinusRef, UnwrapRef1, UnwrapRef1End, KillRefFiles, UnwrapRefSpecial
	Wave RefFilesToKill, RefFilesToUnwrap
	
	
	If(KillRefFiles)
		NumPointsToKill = numpnts(RefFilesToKill)
		i = 0
		do
			Killindex = RefFilesToKill[i]-minFile+5
			RefPhase[Killindex] = NaN
			i+=1
		while(i<NumPointsToKill)
	EndIf
	


	NVAR plusminus, UnwrapAll, Unwrap1, Unwrap1End, Unwrap2, Unwrap2End, Unwrap3, Unwrap3End



	
	NVAR KillHVFiles
	Wave HVFilesToKill
	
	if(KillHVFiles)
		NumPointsToKill = numpnts(HVFilesToKill)
		i = 0
		do
			Killindex = HVFilesToKill[i]-minFile+5
			PolFitPhase[Killindex] = NaN
			i+=1
		while(i<NumPointsToKill)
	endif

	
	if(0)
		//Unwrap special
		PolFitPhase[13] += 2*pi
	endif





	if(0)
		//Fit the reference phase to a quadratic polynomial and adjust the phase
		Make/O/N=3 w_coef
		CurveFit/NTHR=1/TBOX=768 poly 3,  RefPhase /W=RefPhaseError /X=FileNumbers /I=1 /D /R
		TextBox/C/N=$FitTextBoxName/X=10.7/Y=56
		PolFitPhase -= w_coef[0] + w_coef[1]*FileNumbers[p] + w_coef[2]*(FileNumbers[p])^2
		Wave RefPhaseResid = $RefPhaseResidWN
		RefPhaseResid/=RefPhase
		RefPhaseResid*=RefPhase
		ModifyGraph zero(Res_RefPhase)=1
	endif

	if(0)
		Display/K=1/W=(40,100,440,325) PolFitPhase vs PosPad
		ModifyGraph mode=3,marker=8
		ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
		Label left "Contrast"; Label bottom "Distance to Ground Plane (m)"
	endif
	
	
	

	ReducePolWaves(PosPad,PolFitPhase, PolFitPhaseError)
End





Function ExtractPhaseContrast(SeriesName)
	String SeriesName		// ex: "e"
	
	String initposWaveName="pos_"+SeriesName
	String initphaseWaveName = "phase_"+SeriesName
	String initphase_errorWaveName = "phase_error_"+SeriesName
	String initcontrastWaveName = "contrast_"+SeriesName
	String initcontrast_errorWaveName = "contrast_error_"+SeriesName
	
	Wave poswave = $initposWaveName
	Wave phase = $initphaseWaveName
	Wave phase_error = $initphase_errorWaveName
	Wave contrast = $initcontrastWaveName
	Wave contrast_error = $initcontrast_errorWaveName
	
	NVAR LastEclipse, LastFile
	Variable minFile = LastEclipse+1
	Variable maxFile = LastFile		
	
	Make/O/N=(maxFile-minFile+5+1) FileNumbers = LastEclipse-5+p+1
	
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
	String contrastWaveName = NameOfWave(contrast) + "_polfit"
	String contrastErrorWaveName = NameOfWave(contrast_error) + "_polfit"

	Make/O/N=(maxFile-minFile+1) $PhaseWaveName
	Make/O/N=(maxFile-minFile+1) $PhaseErrorWaveName
	Make/O/N=(maxFile-minFile+1) $contrastWaveName
	Make/O/N=(maxFile-minFile+1) $contrastErrorWaveName
	
	Wave PolFitPhase = $PhaseWaveName
	Wave PolFitPhaseError = $PhaseErrorWaveName
	Wave PolFitcontrast = $contrastWaveName
	Wave PolFitcontrastError = $contrastErrorWaveName

	PolFitPhase = phase[x+minFile-1]
	PolFitPhaseError = phase_error[x+minFile-1]
	PolFitcontrast = contrast[x+minFile-1]
	PolFitcontrastError = contrast_error[x+minFile-1]
	
	String RefPhaseWaveName = NameOfWave(PolFitPhase)+"_ref"
	String RefPhaseErrorWaveName = NameOfWave(PolFitPhaseError)+"_ref"
	String RefPhaseResidWN = "Res_"+NameOfWave(PolFitPhase)+"_ref"
	String PhaseFitTextBoxName = "CF_"+NameOfWave(PolFitPhase)+"_ref"
	String RefcontrastWaveName = NameOfWave(PolFitcontrast)+"_ref"
	String RefcontrastErrorWaveName = NameOfWave(PolFitcontrastError)+"_ref"
	String RefcontrastResidWN = "Res_"+NameOfWave(PolFitcontrast)+"_ref"
	String ContrastFitTextBoxName = "CF_"+NameOfWave(PolFitcontrast)+"_ref"
	
	Duplicate/O PolFitPhase $RefPhaseWaveName				//copy the phase data into a wave that will eventually only hold reference phase data
	Duplicate/O PolFitPhaseError $RefPhaseErrorWaveName
	Duplicate/O PolFitcontrast $RefcontrastWaveName				//copy the contrast data into a wave that will eventually only hold reference contrast data
	Duplicate/O PolFitcontrastError $RefcontrastErrorWaveName
	
	Wave RefPhase = $RefPhaseWaveName
	Wave RefPhaseError = $RefPhaseErrorWaveName
	Wave Refcontrast = $RefcontrastWaveName
	Wave RefcontrastError = $RefcontrastErrorWaveName
	
		
		
	//	extract ref phase
	if(1)
		n=0; i=0
		do
			RefPhase[i]=NaN				//This is a HV phase, so erase it from the ref phase wave
			RefPhaseError[i]=NaN
			Refcontrast[i]=NaN				//This is a HV contrast, so erase it from the ref contrast wave
			RefcontrastError[i]=NaN
			
			if(n==9)						// if n=9 we've gone through 10 HV phase measurements, so the next 5 will be reference phases
				n=0						
				i+=5						// so we'll skip the eraser over them in the next iteration
			else
				n+=1
			endif
			
			i+=1
		while(i<numpnts(RefPhase))
		
		i=0
		do										//remove reference phase measurements from phase data
			if(PolFitPhase[i]==RefPhase[i])	// returns true if the point is a reference phase measurement
				PolFitPhase[i]=NaN
				PolFitPhaseError[i]=NaN
				PolFitcontrast[i]=NaN
				PolFitcontrastError[i]=NaN
			endif
			i+=1
		while(i<numpnts(PolFitPhase))
	endif



	// Insert last eclipse data points as a ref phase measurement, pad 
	InsertPoints 0, 5, $PhaseWaveName, $PhaseErrorWaveName, $RefPhaseWaveName, $RefPhaseErrorWaveName, $contrastWaveName, $contrastErrorWaveName, $RefcontrastWaveName, $RefcontrastErrorWaveName, $PosPadWaveName
	i=0
	do
		PolFitPhase[i]=NaN
		PolFitPhaseError[i]=NaN
		RefPhase[i]=phase[i+minFile-6]
		RefPhaseError[i]=phase_error[i+minFile-6]
		PolFitcontrast[i]=NaN
		PolFitcontrastError[i]=NaN
		Refcontrast[i]=contrast[i+minFile-6]
		RefcontrastError[i]=contrast_error[i+minFile-6]
		PosPad[i]=PosWave[FirstHVPosIndex-1]
		i+=1
	while(i<5)
	

	String DataFolderNameInit = GetDataFolder(0)
	String DataFolderName = ReplaceString("'", DataFolderNameInit, "")
	
	//Display/K=1/W=(40,100,550,600) 
	Display/K=1/W=(350,150,1250,650) PolFitPhase vs FileNumbers
	ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
	Label left "Phase"; Label bottom "File number"
	AppendToGraph/L='RefPhase' RefPhase vs FileNumbers
	ModifyGraph mode=3,marker=8
	ErrorBars $RefPhaseWaveName Y,wave=(RefPhaseError,RefPhaseError)
	Label RefPhase "RefPhase"
	ModifyGraph axisEnab(left)={.6,1}
	ModifyGraph axisEnab(RefPhase)={0,0.57}
	ModifyGraph freePos(RefPhase)={0,bottom}
	ModifyGraph lblPosMode(refphase)=1
	TextBox/C/N=text0/A=MC DataFolderName
	TextBox/C/N=text0/X=9.58/Y=45.65
	ModifyGraph gfSize=11
	WindowNamer(DataFolderName+" Phase and Ref Phase")
	
	
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
	Keep3(PolFitcontrast)
	Keep3(PolFitcontrastError)
	Keep3(Refcontrast)
	Keep3(RefcontrastError)
	

	Variable maxUnwrapIndex, NumPointsToKill, Killindex, NumPointsToUnwrap, UnwrapIndex

	NVAR PlusMinusRef, UnwrapRef1, UnwrapRef1End, KillRefFiles//, UnwrapRefSpecial
	Wave RefFilesToKill//, RefFilesToUnwrap
	
	if(UnwrapRef1)
		//Unwrap
		i=0
		maxUnwrapIndex = UnwrapRef1End-minFile+5+1
		do
			RefPhase[i] += plusminusRef*2*pi
			i+=1
		while(i<maxUnwrapIndex)
	endif
	
	
	If(KillRefFiles)
		NumPointsToKill = numpnts(RefFilesToKill)
		i = 0
		do
			Killindex = RefFilesToKill[i]-minFile+5
			RefPhase[Killindex] = NaN
			Refcontrast[Killindex] = NaN
			i+=1
		while(i<NumPointsToKill)
	EndIf
	
//	If(UnwrapRefSpecial)
//		NumPointsToUnwrap = numpnts(RefFilesToUnwrap)
//		i = 0
//		do
//			UnwrapIndex = RefFilesToUnwrap[i]-minFile+5
//			RefPhase[UnwrapIndex] = NaN
//			i+=1
//		while(i<NumPointsToUnwrap)
//	EndIf

	NVAR plusminus, Unwrap1, Unwrap1End, Unwrap2, Unwrap2End, Unwrap3, Unwrap3End

//	If(UnwrapAll)
//		PolFitPhase+=2*pi
//	EndIf


	if(Unwrap1)
		//Unwrap
		i=0
		maxUnwrapIndex = Unwrap1End-minFile+5+1
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<maxUnwrapIndex)
	endif
	
	
	if(Unwrap2)
		//Unwrap again
		i=0
		Variable maxUnwrapIndex2 = Unwrap2End-minFile+5+1
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<maxUnwrapIndex2)
	endif
	
	NVAR KillHVFiles
	Wave HVFilesToKill
	
	if(KillHVFiles)
		NumPointsToKill = numpnts(HVFilesToKill)
		i = 0
		do
			Killindex = HVFilesToKill[i]-minFile+5
			PolFitPhase[Killindex] = NaN
			PolFitcontrast[Killindex] = NaN
			i+=1
		while(i<NumPointsToKill)
	endif

	
	if(0)
		//Unwrap special
		PolFitPhase[13] += 2*pi
	endif





	if(1)
		//Fit the reference phase to a quadratic polynomial and adjust the phase
		Make/O/N=3 w_coef
		CurveFit/NTHR=1/TBOX=768 poly 3,  RefPhase /W=RefPhaseError /X=FileNumbers /I=1 /D /R
		TextBox/C/N=$PhaseFitTextBoxName/X=10.7/Y=56
		PolFitPhase -= w_coef[0] + w_coef[1]*FileNumbers[p] + w_coef[2]*(FileNumbers[p])^2
		Wave RefPhaseResid = $RefPhaseResidWN
		RefPhaseResid/=RefPhase
		RefPhaseResid*=RefPhase
		ModifyGraph zero(Res_RefPhase)=1
	endif
	
	
	Display/K=1/W=(350,150,1250,650) PolFitcontrast vs FileNumbers
	ErrorBars $contrastWaveName Y,wave=(PolFitcontrastError,PolFitcontrastError)
	Label left "Contrast"; Label bottom "File number"
	AppendToGraph/L='Refcontrast' Refcontrast vs FileNumbers
	ModifyGraph mode=3,marker=8
	ErrorBars $RefcontrastWaveName Y,wave=(RefcontrastError,RefcontrastError)
	Label Refcontrast "Ref contrast"
	ModifyGraph axisEnab(left)={.6,1}
	ModifyGraph axisEnab(Refcontrast)={0,0.57}
	ModifyGraph freePos(Refcontrast)={0,bottom}
	ModifyGraph lblPosMode(refcontrast)=1
	TextBox/C/N=text0/A=MC DataFolderName
	TextBox/C/N=text0/X=9.58/Y=45.65
	ModifyGraph gfSize=11
	WindowNamer(DataFolderName+" Contrast and Ref contrast")
	
//	if(1)
//		//Fit the reference contrast to a quadratic polynomial and adjust the contrast
//		Make/O/N=3 w_coef
//		CurveFit/NTHR=1/TBOX=768 poly 3,  Refcontrast /W=RefcontrastError /X=FileNumbers /I=1 /D /R
//		TextBox/C/N=$contrastFitTextBoxName/X=10.7/Y=56
//		if(0)
//			PolFitcontrast /= w_coef[0] + w_coef[1]*FileNumbers[p] + w_coef[2]*(FileNumbers[p])^2
//			PolFitcontrastError /= w_coef[0] + w_coef[1]*FileNumbers[p] + w_coef[2]*(FileNumbers[p])^2
//		else
//			//Variable AvgRefContrastLoc = numVarOrDefault("AvgRefContrast",0.2)
//			NVAR AvgRefContrast
//			PolFitContrast /= AvgRefContrast
//			PolFitContrastError /= AvgRefContrast
//		endif
//		Wave RefcontrastResid = $RefcontrastResidWN
//		RefcontrastResid/=Refcontrast
//		RefcontrastResid*=Refcontrast
//		ModifyGraph zero(Res_Refcontrast)=1
//	endif

	if(1)
		//Fit the reference contrast to a quadratic polynomial and adjust the contrast
		Make/O/N=2 w_coef
		CurveFit/NTHR=1/TBOX=768 line,  Refcontrast /W=RefcontrastError /X=FileNumbers /I=1 /D /R
		TextBox/C/N=$contrastFitTextBoxName/X=10.7/Y=56
		if(0)
			PolFitcontrast /= w_coef[0] + w_coef[1]*FileNumbers[p] 
			PolFitcontrastError /= w_coef[0] + w_coef[1]*FileNumbers[p] 
		else
			//Variable AvgRefContrastLoc = numVarOrDefault("AvgRefContrast",0.2)
			NVAR AvgRefContrast
			PolFitContrast /= AvgRefContrast
			PolFitContrastError /= AvgRefContrast
		endif
		Wave RefcontrastResid = $RefcontrastResidWN
		RefcontrastResid/=Refcontrast
		RefcontrastResid*=Refcontrast
		ModifyGraph zero(Res_Refcontrast)=1
	endif



	if(0)
		Display/K=1/W=(40,100,440,325) PolFitPhase vs PosPad
		ModifyGraph mode=3,marker=8
		ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
		Label left "Phase"; Label bottom "Distance to Ground Plane (m)"
	endif
	
	
	
	NVAR InvertPhase

	if(InvertPhase)							//Make the phase shift positive
		PolfitPhase*=-1
	endif

	ReducePhaseContrastWaves(PosPad,PolFitPhase, PolFitPhaseError, PolFitContrast, PolFitContrastError)
End







Function ExtractPhaseDataV2(poswave, phase, phase_error)
	Wave poswave, phase, phase_error
	
	NVAR LastEclipse, LastFile
	Variable minFile = LastEclipse+1
	Variable maxFile = LastFile		
	
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
	String RefPhaseResidWN = "Res_"+NameOfWave(PolFitPhase)+"_ref"
	
	Duplicate/O PolFitPhase $RefPhaseWaveName				//copy the phase data into a wave that will eventually only hold reference phase data
	Duplicate/O PolFitPhaseError $RefPhaseErrorWaveName
	
	Wave RefPhase = $RefPhaseWaveName
	Wave RefPhaseError = $RefPhaseErrorWaveName
	
		
		
	//	extract ref phase
	if(1)
		n=0; i=0
		do
			RefPhase[i]=NaN				//This is a HV phase, so erase it from the ref phase wave
			RefPhaseError[i]=NaN
			
			if(n==9)						// if n=9 we've gone through 10 HV phase measurements, so the next 5 will be reference phases
				n=0						
				i+=5						// so we'll skip the eraser over them in the next iteration
			else
				n+=1
			endif
			
			i+=1
		while(i<numpnts(RefPhase))
		
		i=0
		do										//remove reference phase measurements from phase data
			if(PolFitPhase[i]==RefPhase[i])	// returns true if the point is a reference phase measurement
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

	//Display/K=1/W=(40,100,550,600) 
	Display/K=1/W=(140,150,850,600) PolFitPhase
	ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
	Label left "Phase"; Label bottom "File number minus 1"
	AppendToGraph/L='RefPhase' RefPhase
	ModifyGraph mode=3,marker=8
	ErrorBars $RefPhaseWaveName Y,wave=(RefPhaseError,RefPhaseError)
	Label RefPhase "RefPhase"
	ModifyGraph axisEnab(left)={.7,1}
	ModifyGraph axisEnab(RefPhase)={0,0.67}
	ModifyGraph freePos(RefPhase)={0,bottom}
	ModifyGraph lblPosMode(refphase)=1
	String DataFolderNameInit = GetDataFolder(0)
	String DataFolderName = ReplaceString("'", DataFolderNameInit, "")
	TextBox/C/N=text0/A=MC DataFolderName
	TextBox/C/N=text0/X=9.58/Y=45.65
	
	
	
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
	

	Variable maxUnwrapIndex, NumPointsToKill, Killindex, NumPointsToUnwrap, UnwrapIndex

	NVAR PlusMinusRef, UnwrapRef1, UnwrapRef1End, KillRefFiles, UnwrapRefSpecial
	Wave RefFilesToKill, RefFilesToUnwrap
	
	if(UnwrapRef1)
		//Unwrap
		i=0
		maxUnwrapIndex = UnwrapRef1End-minFile+5+1
		do
			RefPhase[i] += plusminusRef*2*pi
			i+=1
		while(i<maxUnwrapIndex)
	endif
	
	
	If(KillRefFiles)
		NumPointsToKill = numpnts(RefFilesToKill)
		i = 0
		do
			Killindex = RefFilesToKill[i]-minFile+5
			RefPhase[Killindex] = NaN
			i+=1
		while(i<NumPointsToKill)
	EndIf
	
	If(UnwrapRefSpecial)
		NumPointsToUnwrap = numpnts(RefFilesToUnwrap)
		i = 0
		do
			UnwrapIndex = RefFilesToUnwrap[i]-minFile+5
			RefPhase[UnwrapIndex] = NaN
			i+=1
		while(i<NumPointsToUnwrap)
	EndIf
	

	NVAR plusminus, UnwrapAll, Unwrap1, Unwrap1End, Unwrap2, Unwrap2End, Unwrap3, Unwrap3End

	If(UnwrapAll)
		PolFitPhase+=2*pi
	EndIf


	if(Unwrap1)
		//Unwrap
		i=0
		maxUnwrapIndex = Unwrap1End-minFile+5+1
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<maxUnwrapIndex)
	endif
	
	
	if(Unwrap2)
		//Unwrap again
		i=0
		Variable maxUnwrapIndex2 = Unwrap2End-minFile+5+1
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<maxUnwrapIndex2)
	endif
	
	NVAR KillHVFiles
	Wave HVFilesToKill
	
	if(KillHVFiles)
		NumPointsToKill = numpnts(HVFilesToKill)
		i = 0
		do
			Killindex = HVFilesToKill[i]-minFile+5
			PolFitPhase[Killindex] = NaN
			i+=1
		while(i<NumPointsToKill)
	endif

	
	if(0)
		//Unwrap special
		PolFitPhase[13] += 2*pi
	endif




/// fix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!???????????	
	if(1)
		//Fit the reference phase to a quadratic polynomial and adjust the phase
		Make/O/N=3 w_coef
		CurveFit/NTHR=1/TBOX=768 poly 3,  RefPhase /W=RefPhaseError /I=1 /D /R
		TextBox/C/N=CF_phase_polfit_ref/X=7.95/Y=44.71
		PolFitPhase -= w_coef[0] + w_coef[1]*x + w_coef[2]*x^2
		Wave RefPhaseResid = $RefPhaseResidWN
		RefPhaseResid/=RefPhase
		RefPhaseResid*=RefPhase
//		Duplicate/O RefPhaseResid RefPhaseResidTemp
//		RefPhaseResid /= RefPhaseResid
//		RefPhaseResid *= RefPhaseResidTemp
		ModifyGraph zero(Res_RefPhase)=1
		TextBox/C/N=CF_phase_polfit_ref/X=5.40/Y=52.00
	endif

	if(0)
		Display/K=1/W=(40,100,440,325) PolFitPhase vs PosPad
		ModifyGraph mode=3,marker=8
		ErrorBars $PhaseWaveName Y,wave=(PolFitPhaseError,PolFitPhaseError)
		Label left "Phase"; Label bottom "Distance to Ground Plane (m)"
	endif
	
	
	
	NVAR InvertPhase

	if(InvertPhase)							//Make the phase shift positive
		PolfitPhase*=-1
	endif


	ReducePolWaves(PosPad,PolFitPhase, PolFitPhaseError)
End





// Extracts the phase data from a the standard data sets (10 HV on, 5 Ref...)
// minfile is the first hv file
// be careful about file vs index. 
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
	
	Duplicate/O PolFitPhase $RefPhaseWaveName				//copy the phase data into a wave that will eventually only hold reference phase data
	Duplicate/O PolFitPhaseError $RefPhaseErrorWaveName
	
	Wave RefPhase = $RefPhaseWaveName
	Wave RefPhaseError = $RefPhaseErrorWaveName
	
		
		
	//	extract ref phase
	if(1)
		n=0; i=0
		do
			RefPhase[i]=NaN				//This is a HV phase, so erase it from the ref phase wave
			RefPhaseError[i]=NaN
			
			if(n==9)						// if n=9 we've gone through 10 HV phase measurements, so the next 5 will be reference phases
				n=0						
				i+=5						// so we'll skip the eraser over them in the next iteration
			else
				n+=1
			endif
			
			i+=1
		while(i<numpnts(RefPhase))
		
		i=0
		do										//remove reference phase measurements from phase data
			if(PolFitPhase[i]==RefPhase[i])	// returns true if the point is a reference phase measurement
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
	
	
	
	//g	
//	if(1)
//		//Unwrap
//		i=0
//		Variable maxUnwrapIndex = 60
//		Variable plusminus = 1	// plusminus = 1 if phase needs to be adjusted up, plusminus = -1 if phase needs to be adjusted down
//		do
//			PolFitPhase[i] += plusminus*2*pi
//			i+=1
//		while(i<maxUnwrapIndex)
//	endif
//	
//	if(1)
//		//Unwrap again
//		i=0
//		Variable maxUnwrapIndex2 = 130
//		do
//			PolFitPhase[i] += plusminus*2*pi
//			i+=1
//		while(i<maxUnwrapIndex2)
//	endif
//	
//	if(1)
//		PolFitPhase[68]=nan
//	endif
//
//	
//	if(0)
//		//Unwrap special
//		PolFitPhase[118] -= 2*pi
//	endif



//h
//	if(1)
//		//Unwrap
//		i=0
//		Variable maxUnwrapIndex = 105-70+5
//		Variable plusminus = -1	// plusminus = 1 if phase needs to be adjusted up, plusminus = -1 if phase needs to be adjusted down
//		do
//			PolFitPhase[i] += plusminus*2*pi
//			i+=1
//		while(i<maxUnwrapIndex)
//	endif
//	
//	if(1)
//		//Unwrap again
//		i=0
//		Variable maxUnwrapIndex2 = 181-70+5
//		do
//			PolFitPhase[i] += plusminus*2*pi
//			i+=1
//		while(i<maxUnwrapIndex2)
//	endif
//	
//	if(1)
//		PolFitPhase[98,99]=nan
//	endif
//
//	
//	if(0)
//		//Unwrap special
//		PolFitPhase[118] -= 2*pi
//	endif


//i
//	i=0
//	Variable maxUnwrapIndex0 = 160-75
//	Variable plusminus0 = -1	// plusminus = 1 if phase needs to be adjusted up, plusminus = -1 if phase needs to be adjusted down
//	do
//		RefPhase[i] += plusminus0*2*pi
//		i+=1
//	while(i<maxUnwrapIndex0)
//	
//	RefPhase[92]-=2*pi
//	RefPhase[108] = nan
//
//
//	if(1)
//		//Unwrap
//		i=0
//		Variable maxUnwrapIndex = 81-minFile+5
//		Variable plusminus = 1	// plusminus = 1 if phase needs to be adjusted up, plusminus = -1 if phase needs to be adjusted down
//		do
//			PolFitPhase[i] += plusminus*2*pi
//			i+=1
//		while(i<maxUnwrapIndex)
//	endif
//	
//	if(1)
//		//Unwrap again
//		i=0
//		Variable maxUnwrapIndex2 = 156-minFile+5
//		do
//			PolFitPhase[i] += plusminus*2*pi
//			i+=1
//		while(i<maxUnwrapIndex2)
//	endif
//	
//	if(0)
//		PolFitPhase[98,99]=nan
//	endif
//
//	
//	if(1)
//		//Unwrap special
//		PolFitPhase[13] += 2*pi
//	endif




//j
	PolFitPhase+=2*pi

	if(1)
		//Unwrap
		i=0
		Variable maxUnwrapIndex = 121-minFile+5
		Variable plusminus = 1	// plusminus = 1 if phase needs to be adjusted up, plusminus = -1 if phase needs to be adjusted down
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<maxUnwrapIndex)
	endif
	
	if(0)
		//Unwrap again
		i=0
		Variable maxUnwrapIndex2 = 156-minFile+5
		do
			PolFitPhase[i] += plusminus*2*pi
			i+=1
		while(i<maxUnwrapIndex2)
	endif
	
	if(0)
		PolFitPhase[98,99]=nan
	endif

	
	if(0)
		//Unwrap special
		PolFitPhase[13] += 2*pi
	endif




/// fix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!???????????	
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
	

	if(stringmatch("h",series))							//Make the phase shift positive
		PolfitPhase*=-1
	endif
	
	if(stringmatch("g",series) && 0)
		PolFitPhase[25]=nan
	endif

	ReducePolWaves(PosPad,PolFitPhase, PolFitPhaseError)
End







///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






Static Function RemoveNaNsXYZAB(theXWave, theYWave, theZWave, theAWave, theBWave) //adapted from Remove Points.ipf
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





Static Function RemoveNaNsXYZABCD(theXWave, theYWave, theZWave, theAWave, theBWave, theCWave, theDWave) //adapted from Remove Points.ipf
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



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
	
	
	
	
	





// also possibly see local procedure for these parameters

// velocity parameters
//Constant velocity = 1452.5		//g
//Constant velocity = 1450.4		//h
//Constant velocity = 1449.95		//m
//Constant velocity = 1990.84		// b 4/7/09

//Constant sigma = 112.8 		 	//g
//Constant sigma = 112.5 		 	//h
//Constant sigma = 112.4 		 	//m
//Constant sigma = 154.4				//b 4/7/09

	
// y offset parameter
//Constant offset = 73.8e-6 	//g
//Constant offset = 44.9e-6	//g counts =50%
//Constant offset = 49.6e-6 	//h
//Constant offset =   54.6e-6 	//h 2IFM
//Constant offset =   43.3e-6 	//h 2IFM counts=70%
//Constant offset = 22.4e-6		//h counts=50%
//Constant offset = 28.8e-6		//m counts contrast = 50%
//Constant offset = 2.13e-5		//m counts=50%
//Constant offset =   5.35e-05		//m counts = 80%
//Constant offset = -35.5e-6	//q 10/16/08
//Constant offset = 4.39542e-05 // b 4/7/09 counts contrast= 50%
//Constant offset = 1.73e-5		// b 4/7/09 counts=50%

//Constant offset = 35e-6 

//Constant HVmon = 0.9165	//g
//Constant HVmon = 1.143		//h
//Constant HVmon = 0.961 		//m
//Constant HVmon = 1.010	//q 10/16/08
//Constant HVmon= 1.57		// b 4/7/09
//Constant HVmon = .75