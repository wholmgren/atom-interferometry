#pragma rtGlobals=1		// Use modern global access method.

#include ":WindowNamer"
#include ":Raw Beam Convolution"


Function ProcessRawBeam(fileName)
	String fileName
	
	Make/O/WAVE/N=6 WaveRefs
	// WaveRefs[0] = counts
	// WaveRefs[1] = position
	// WaveRefs[2] = counts error
	
	//String path = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:090630:diffscans:"
	String path = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:090727:diffscans:"
	
	LoadAndDisplayScan2(path, fileName, 0)
	Wave countsWave = WaveRefs[0]; Wave positionWave = WaveRefs[1]
	
	Variable maxPosition = 0.8e3
	RezeroScan2(positionWave, countsWave)
	
	CropScan(positionWave, countsWave, maxPosition)
	
	Wave errWave = GenErrWave(countsWave)
	WaveRefs[2] = errWave
	
		
		Variable resampleRate = .1
		resampleWaves(positionWave, countsWave, errWave, resampleRate)
		Wave resampledCounts = WaveRefs[3]; Wave resampledPosition = WaveRefs[4]
		Wave resampledErr = GenErrWave(resampledCounts)
		WaveRefs[5] = resampledErr
	
		DisplayScan(ResampledCounts, resampledPosition, resampledErr)
		
	//GoFitRawBeam(ResampledCounts, resampledPosition, resampledErr)
	GoFitRawBeamGauss(ResampledCounts, resampledPosition, resampledErr)
	//GoFitRawBeam3sGauss(ResampledCounts, resampledPosition, resampledErr)
		
		
		CleanUpGraph(resampledCounts, resampledErr)
	
	SetAxis left .1,*
End





Function LoadAndDisplayScan2(path, fileName, displayYesNo)
	String path, filename
	Variable displayYesNo
	//Variable motorCountsDivider

	String datay = fileName + "y"	// counts
	String datax = fileName + "x" 	// motor position
	
	Variable k=1
	do
		if(StringMatch(WaveList(datay,";",""),""))
			break
		else
			k+=1
			datay = fileName+"y"+num2str(k)
			datax = fileName+"x"+num2str(k)
		endif
	while(1)

	//Wave wave0; wave0=0; Wave wave1; wave1=0
	LoadWave/Q/G/D/W/N/O/A path + filename
	//Duplicate/O wave0 $datay; Duplicate/o wave1 $datax 		// copy counts and motor position into waves named $datay and $datax
	Duplicate/O counts $datay; Duplicate/o position $datax
	Wave positionWave = $datax; Wave countsWave = $datay;	// create references to the datay and datax waves
	
	if(WaveMax(positionWave)>(2^24-100))					// if we're counting down from 2^24 we'll run into trouble with the first points at position = 0
		positionWave-=2^24*(positionWave[p]>1000000)			// so add 2^24 to any positions that are less than a small number
	endif
	
	//positionWave *= .4		// 400 nm to 1 um
	
	If(displayYesNo)
		Display/k=1/W=(100,75,1000,400) countsWave vs positionWave	
		ModifyGraph mode=3,marker=8
		TextBox/C/N=text0/F=0/A=LT fileName
	EndIf
	
	if(waveexists(WaveRefs))
		Wave/Wave WaveRefs
		WaveRefs[0] = countsWave; WaveRefs[1] = positionWave;		
	endif
End




Function RezeroScan2(positionWave, countsWave)
	Wave positionWave, countsWave
	
	Variable/G K0 = 1,K1 = 200,K2 = 1600,K3 = 200;	//provide initial guesses so that the fit is less susceptible to noise spikes
	Make/O/T/N=1 RezeroConstraints = {"K0<10", "K1<1000", "K1>20", "K3>130", "K3<400"}
	
	Variable v_FitError=0
	Variable V_FitOptions=4 
	
	CurveFit/G/Q/NTHR=0 gauss  countsWave  /X=positionWave /C=RezeroConstraints
	Wave w_coef
//	print w_coef
	print GetRTError(1)
		
	if(v_fiterror!=0)
		K0 = 1;K1 = 500;K2 = 1600;K3 = 100;	//provide initial guesses so that the fit is less susceptible to noise spikes
		CurveFit/G/NTHR=0 gauss  countsWave  /X=positionWave /C=RezeroConstraints
		print GetRTError(1)
	endif
	if(v_fiterror!=0)	// 	if fit doens't work perfectly try a different initial guess
		K2*=-1
		v_FitError=0
		CurveFit/Q/G/NTHR=0 gauss  countsWave  /X=positionWave /C=RezeroConstraints
		print GetRTError(1)
	endif
	if(v_fiterror!=0)	// 	if fit doens't work perfectly try a different initial guess
		K2=-1400
		v_FitError=0
		CurveFit/G/Q/NTHR=0 gauss  countsWave  /X=positionWave /C=RezeroConstraints
		print GetRTError(1)
	endif
	if(v_fiterror!=0)	// 	if fit doens't work perfectly try a different initial guess
		K2*=-1
		v_FitError=0
		CurveFit/Q/G/NTHR=0 gauss  countsWave  /X=positionWave /C=RezeroConstraints
		print GetRTError(1)
	endif
	if(v_fiterror!=0)	// 	if fit doens't work perfectly try a different initial guess
		K2=-1000
		v_FitError=0
		CurveFit/Q/G/NTHR=0 gauss  countsWave  /X=positionWave /C=RezeroConstraints
		print GetRTError(1)
	endif
	if(v_fiterror!=0)	// 	if fit doens't work perfectly try a different initial guess
		K2=6709352
		V_FitError=0
		CurveFit/Q/G/NTHR=0 gauss  countsWave  /X=positionWave /C=RezeroConstraints
		print GetRTError(1)
		if(v_fitError!=0)
			print "could not rezero scan!"
		endif
	endif
	
	positionWave -= w_coef[2]
	
	//print w_coef
End




Function CropScan(pos, atomcounts, halfwidth)
	Wave pos, atomcounts
	Variable halfwidth
	
	Variable flip = 0
	if(pos[0] > pos[numpnts(pos)-1])
		pos*=-1
	endif
	
	Variable V_flag
	
	FindLevel/Q pos, -halfwidth
	Variable minpnt = floor(V_LevelX)
	FindLevel/Q pos, halfwidth
	Variable maxpnt = ceil(V_LevelX)
	
	if(V_flag==1)
		v_flag=0
		halfwidth*=0.8
		FindLevel/Q pos, -halfwidth
		 minpnt = floor(V_LevelX)
		FindLevel/Q pos, halfwidth
		 maxpnt = ceil(V_LevelX)
	endif
	
	
	Variable pnts = maxpnt-minpnt
	
	Make/D/O/N=(pnts) AtomCountsCrop, PosCrop
	
	AtomCountsCrop = atomcounts[p+minpnt]; Duplicate/O AtomCountsCrop $nameofwave(atomcounts)
	PosCrop = pos[p+minpnt]; Duplicate/O PosCrop $nameofwave(pos)
	
	If(flip)
		PosCrop*=-1
	Endif
	
	KillWaves AtomCountsCrop, PosCrop
End	




Function/WAVE GenErrWave(counts)
	Wave counts
	Wave/Wave WaveRefs
	
	Variable binLength = 2 //in milliseconds
	
	String ErrWaveName = nameofwave(counts)+"Err"
	Duplicate/O counts $ErrWaveName
	Wave CountsErr = $ErrWaveName
	
	CountsErr = sqrt(counts*binLength)/binLength
	
	//WaveRefs[2] = CountsErr
	
	Return CountsErr
End



Function ResampleWaves(positionWave, countsWave, errWave, rate)
	Wave positionWave, countsWave, errWave
	Variable rate
	
	String NewPosName = nameofwave(positionWave) + "_samp"
	String newCountsName = nameofwave(countsWave) + "_samp"
	String NewErrName = nameofwave(errWave) + "_samp"
	
	Duplicate/O positionWave $NewPosName; Wave posResampled = $NewPosName
	Duplicate/O countsWave $NewCountsName; Wave countsResampled = $NewCountsName
	Duplicate/O errWave $NewErrName; Wave errResampled = $NewErrName
	
	Resample/RATE=(rate) posResampled, countsResampled, errResampled
	
	Wave/Wave WaveRefs
	WaveRefs[3] = countsResampled
	WaveRefs[4] = posResampled
	WaveRefs[5] = errResampled
End



Function DisplayScan(counts, pos, err)
	Wave counts, pos, err
	
	Display/k=1/W=(100,75,1100,500) counts vs pos			// display a nice graph (will not ask for confirmation when closing)
	ModifyGraph mode=3,marker=8
	TextBox/C/N=text0/F=0/A=LT fileName
	ErrorBars $nameofwave(counts) Y, wave=(err, err)
	
	WindowNamer(nameofwave(counts))
End



Function CleanUpGraph(counts, error)
	Wave counts, error
	
	String CountsFitWaveName = "fit_"+nameofwave(counts)
	String ResWaveName = "Res_"+nameofwave(counts)
	
	ModifyGraph lsize($CountsFitWaveName)=2,rgb($CountsFitWaveName)=(0,0,65535)
	ModifyGraph zero(Res_Left)=1
	ModifyGraph zeroThick(Res_Left)=2
	ModifyGraph axisOnTop(Res_Left)=1
	ModifyGraph log(left)=1
	ModifyGraph msize($ResWaveName)=1
	SetAxis left 1,*
	Label bottom "Detector position (um)"
	Label left "Atom flux (kCounts/s)"
	ModifyGraph nticks(bottom)=10
	//ErrorBars $ResWaveName Y,wave=(error,error)
End



Function BatchProcessRawBeam(filesToProcess)
	Wave/T filesToProcess
	
	Variable NumFilesToProcess = numpnts(filestoProcess)
	Make/O/N=(5,NumFilesToProcess) results
	
	Wave w_coef
	
	Variable n=0
	do
		ProcessRawBeam(filesToProcess[n])
		results[][n] = w_coef[p]
		n+=1
	while(n<NumFilesToProcess)
	
End