#pragma rtGlobals=1		// Use modern global access method.

#include "RemoveNaNs"
#include "Pol Fit Routine"
#include "Fringe Fit"



Constant chopL = 1.27068   //1.269


Window ChopperFitControlTable() : Table
	PauseUpdate; Silent 1		// building window...
	Make/O/T ParameterNum = {"K0","K1","K2","K3","K4"}
	Make/O/T ParameterName = {"Flow vel", "v ratio", "c_0", "phi12avg / pi", "phi12diff / pi"}
	Make/O/D InitialGuesses = {2800,18,0.3,1,0}
	Make/O/D/N=(numpnts(InitialGuesses)) w_coef, w_sigma
	Make/O HoldWave = {0,0,0,0,1}
	Make/O/T Constraints1 = {"","","","",""}
	Make/O/T Constraints2 = {"","","","",""}
	Make/O/T/N=0 Constraints
	Edit/W=(5,44,850,251) ParameterNum,ParameterName,InitialGuesses,w_coef,W_sigma, HoldWave, Constraints1, Constraints2, Constraints
EndMacro








Function ChopperDCphase(series_name, OnekTwok)
	String series_name
	Variable OnekTwok
	
	String onektwokStr = ""
	if(OnekTwok==1)
		onektwokStr = "1k2k"
	elseif(OnekTwok == 2)
		onektwokStr = "2k2k"
	endif
	
	String inputphase = "phase" + onektwokStr + "_" + series_name
	String inputphase_error = "phase_error" + onektwokStr + "_" + series_name

	String phase_on_name = "phase" + onektwokStr + "_on_" + series_name
	String phase_error_on_name = "phase_error" + onektwokStr + "_on_" + series_name
	
	String phase_ref_name = "phase" + onektwokStr + "_ref_" + series_name
	String phase_error_ref_name = "phase_error" + onektwokStr + "_ref_" + series_name
	
	String fileNumbers_name = "file_index_" + series_name
	Wave fileNumbers = $fileNumbers_name
	
	Wave phase_in = $inputphase
	Wave phase_error_in = $inputphase_error
	
	Duplicate/O phase_in $phase_on_name, $phase_ref_name
	Duplicate/O phase_error_in $phase_error_on_name, $phase_error_ref_name
	
	Wave phase_on = $phase_on_name
	Wave phase_error_on = $phase_error_on_name
	Wave phase_ref = $phase_ref_name
	Wave phase_error_ref = $phase_error_ref_name
	
	Variable NumOnOff = 5
	
	phase_on *= (mod(Floor(p/NumOnOff),2)==1)	// keep every other set of (NumOnOff) points
	//phase_on *= (mod(p,NumOnOff)!=0)	// delete first points from each set
	phase_ref *= (mod(Floor(p/NumOnOff),2)==0)
	
	String killFilesName = "KillFiles_"+series_name
	Wave/Z KillFiles = $killFilesName
	
	phase_on[15]=NaN	// the voltage is coming up to full strength in this file
	
	If(waveexists(killfiles))
		Variable i
		Variable killpnt
		For(i=0; i<numpnts(killFiles); i+=1)
			killpnt = killFiles[i]-1
			phase_on[killpnt]=nan
			phase_ref[killpnt]=nan
		EndFor
	endif
	
	phase_on = phase_on==0 ? nan : phase_on
	phase_ref = phase_ref==0 ? nan : phase_ref
	
	phase_error_on = numtype(phase_on)==0 ? phase_error_on : nan
	phase_error_ref = numtype(phase_ref)==0 ? phase_error_ref : nan
	
	Display/K=1/W=(50,50,750,500) phase_on vs fileNumbers
	AppendToGraph phase_ref vs fileNumbers
	ModifyGraph zero(left)=1
	WindowNamer("Chop DC phase " + oneKTwoKStr + " " + series_name)
	ModifyGraph nticks(bottom)=10,minor(bottom)=1
	
	Variable refPolyOrder = 4
	Make/O w_coef
	CurveFit/NTHR=0 poly refPolyOrder,  phase_ref /W=phase_error_ref /I=1 /R
	
	//String fitwavename = "fit_"+phase_ref_name
	//Wave fitwave = $fitwavename
	String reswavename = "res_"+phase_ref_name
	Wave reswave = $reswavename
	
	Variable k
	for(k=0; k<refPolyOrder; k+=1)
		phase_on -= w_coef[k]*p^k
		phase_ref -= w_coef[k]*p^k
	endfor
//	phase_on -= w_coef[0]+w_coef[1]*p+w_coef[2]*p^2
//	phase_ref -= w_coef[0]+w_coef[1]*p+w_coef[2]*p^2
	//fitwave -= w_coef[0]+w_coef[1]*p+w_coef[2]*p^2
	
	ModifyGraph mode=3
	ModifyGraph marker=8
	Label left "phase shift " + onektwokStr
	Label bottom "file number"
	
	reswave/=phase_ref; reswave*=phase_ref
	
	Wavestats/Q/R=[NumOnOff,2*NumOnOff-1] phase_on
	Variable phiHalfV = V_avg
	Wavestats/Q/R=[3*NumOnOff, 4*NumOnOff-1] phase_on
	Variable phiFullV = V_avg
	
	if(phiHalfV > pi)
		phase_on -= 2*pi; phiHalfV -=2*pi; phiFullV -=2*pi
	elseif(phiHalfV < -pi)
		phase_on += 2*pi; phiHalfV +=2*pi; phiFullV +=2*pi
	elseif(phiHalfV < 0 && phiFullV > 0)
		phase_on -= 2*pi*(mod(Floor(p/NumOnOff),2)==1)*(p>=(numOnOff*3))
	elseif(phiHalfV > 0 && phiFullV < 0)
		phase_on += 2*pi*(mod(Floor(p/NumOnOff),2)==1)*(p>=(numOnOff*3))
	endif
	
	Extract phase_on, phase_fullon, abs(phase_on) > abs(phiHalfV*2)				// phihalfv*2 should still be 1/2 of full phase (since phi is proportional to v^2)
	Extract fileNumbers, fileNumbers_fullon, abs(phase_on) > abs(phiHalfV*2)
	
	Wavestats/Q/W/ALPH=.34 phase_fullon
	printf "phi: %g\r std dev: %g\r std err: %g\r", V_avg, V_sdev, V_sem
	printf "phi %%: %g\r std dev %%: %g\r std err %%: %g\r", V_avg/pi, V_sdev/pi, V_sem/pi
	printf "phase shift / pi = %.3f pm %.3f\r", V_avg/pi, V_sem/pi
End





Function GoFitChopperDCphase(series_name)
	String series_name

	String inputphase = "phase_on_" + series_name
	String inputphase_error = "phase_error_on_" + series_name
	String inputvoltage = "hvpad_" + series_name	
	
	Wave phaseWave = $inputphase
	Wave errWave = $inputphase_error
	Wave voltageWave = $inputvoltage
	
	Variable ttime= stopMSTimer(-2)
	
	Make/D/O w_coef={-1.75}
	
	FuncFit/M=2/TBOX=768 FitChopperDCphase W_coef phaseWave /X=voltageWave /R /I=1 /W=errWave 
	
	KillWaves/Z phi_volt_vel, phi_volt_vel_re, phi_volt_vel_im
End

Function FitChopperDCphase(pw,yw,xw) : FitFunc
	Wave pw, yw, xw
	
	Variable prefactor = pw[0]
	
	Variable/G vPnts = 50
	Variable nSigmas = 5
	Variable FlowV = 3000
	Variable vratio = 17
	Variable sigmaV = flowV/vratio
	Variable minVel = FlowV*(1-nSigmas/15)//; minVel = 0
	Variable maxVel = FlowV*(1+nSigmas/15)//; maxVel = 5000
	Variable velStep = (maxVel-minVel)/vPnts
	
	Make/D/O/N=(vPnts) velwave, velprobs, velprobsVel	
	SetScale/I x minVel,maxVel,"m/s", velwave, velprobs, velprobsvel
	velwave = x
	velprobs = x^3 * exp( - (x-flowV)^2/(2*sigmaV^2))
	Variable velprobsNorm = area(velprobs)
	velprobs/=velprobsNorm
	velprobsVel = velprobs*velwave
	Variable AvgV = area(velprobsVel)		// 	avgv = int[ v*P(v) dv, {0, inf}]

	Make/D/O/N=(numpnts(yw),vPnts) phi_volt_vel, phi_volt_vel_re, phi_volt_vel_im
	SetScale/I y minVel,maxVel,"m/s", phi_volt_vel, phi_volt_vel_re, phi_volt_vel_im
	
	phi_volt_vel = prefactor * xw[p]^2 / velwave[q]^2
	phi_volt_vel_re = velprobs[q]*cos(phi_volt_vel)
	phi_volt_vel_im = velprobs[q]*sin(phi_volt_vel)
	
	Integrate/DIM=1 phi_volt_vel_re, phi_volt_vel_im
	
	Duplicate/O yw phi_volt, phi_volt_flowV, c_volt
	
	phi_volt_flowv = prefactor*xw^2 / avgV^2
	
	phi_volt = atan2(phi_volt_vel_im[p][vPnts],phi_volt_vel_re[p][vPnts])
	c_volt = sqrt(phi_volt_vel_re[p][vpnts]^2 + phi_volt_vel_im[p][vpnts]^2)
	
	yw = phi_volt + floor((phi_volt_flowv - phi_volt)/pi)*pi
End



Function PhaseChoppers2IFM()
	Variable ttime= stopMSTimer(-2)
	
	Variable vPnts = 50
	Variable nSigmas = 5
	Variable FlowV =1860
	Variable vratio = 13
	Variable sigV = FlowV/vratio
	Variable VelNorm
	
	Variable vpi = flowv
	Variable AsymmetryFactor = 1.0
	
	Variable tpnts = 100
	Variable dutycycle = 1/2+.0001
	
	Variable posPnts =3
	Variable centerPos = .75e-3
	Variable beamWidth = 100e-6
	Variable minPos = centerPos-beamWidth/2
	Variable maxPos = centerPos+beamWidth/2

	Variable a = 1e-3	// ground plane - wire dist
	Variable radius = 0.775e-3
	Variable b = a*sqrt(1+2*radius/a)
	
	//Variable chopL = 1.4	//meter
	Variable chopL = 1.27069
	
	Variable fPnts = 200
	Variable minChopFreq = 2500
	Variable maxChopFreq = 2600
	//Variable maxOmega = 2*pi*maxChopFreq
	
	Variable/G minVel = FlowV*(1-nSigmas/vratio)//; minVel = 0
	Variable/G maxVel = FlowV*(1+nSigmas/vratio)//; maxVel = 5000
	Variable/G velStep = (maxVel-minVel)/vPnts
	
	Make/D/O/N=(vPnts) velwave = minVel+p*velStep
	Duplicate/O velwave velprobs 
	SetScale/I x minVel,maxVel,"m/s", velprobs
	Velprobs = exp(-(x-FlowV)^2/(2*SigV^2))*x^3
	VelNorm = sqrt(2*pi)*sigv*FlowV*(FlowV^2+3*sigv^2)
	//Variable velsum = sum(velprobs,-inf,inf)
	Velprobs/=velNorm
	
	Make/D/O/N=(fPnts, vPnts, tPnts) phi_fvt_re, phi_fvt_im, tTest, phi_fvt_diff
	Make/D/O/N=(fPnts, vPnts) phi_fv_re, phi_fv_im
	Make/D/O/N=(fPnts) phi_f_re, phi_f_im, phi_f, contrast_f
	SetScale/I x minChopFreq, maxChopFreq, "Hz", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_f_re, phi_f_im, phi_f, contrast_f, ttest, phi_fvt_diff
	SetScale/I y minVel, maxVel, "m/s", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, ttest, phi_fvt_diff
	SetScale/P z 0, 1/tpnts, "t*f", phi_fvt_re, phi_fvt_im, ttest, phi_fvt_diff
	
	Make/D/O/N=(posPnts) beamPos
	SetScale/i x minPos, maxPos, "m", beamPos
	beampos=x
	if(pospnts==1)
		beampos=centerpos
	endif
	Duplicate/O phi_f_re, phi_f_re_temp, phi_f_im_temp
	phi_f_re=0;phi_f_im=0;phi_f_re_temp=0; phi_f_im_temp=0
	
	Print "timer assign vars, make waves: ", (StopMSTimer(-2)-ttime)*10^-6	

	Variable phiint=0//100/pi

	Variable i = 0
	for(i=0; i<numpnts(beampos); i+=1)
		// foftsqrd( time(s), frequency(Hz) )
		//	MultiThread ttest = foftsqrd( z/x , x)
		//	MultiThread phi_fvt_diff = pi*(vpi/y)^2* ( foftsqrd( z/x , x) - foftsqrd( z/x + chopL/y, x) )	
		print "i= ", i, "; beampos[i] = ", beampos[i]
		MultiThread phi_fvt_diff =1.0*pi*(b^2-centerPos^2)/(b^2-beampos[i]^2)*(vpi/y)^2* ( AsymmetryFactor*foftsqrdexp( z/x , x) - foftsqrdexp( z/x + chopL/y, x) + phiint )
		//MultiThread phi_fvt_diff = pi*(b^2-centerPos^2)/(b^2-beampos[i]^2)*(vpi/y)^2* ( phiChopt2( z/x , x) - phiChopt2( z/x + chopL/y, x) + phiint )	
		//	MultiThread phi_fvt_diff = pi*( foftsqrd( z/x , x) - foftsqrd( z/x + chopL/y, x) )	
		//	MultiThread phi_fvt_diff = pi*( foftsqrdexp( z/x , x) - foftsqrdexp( z/x + chopL/y, x) )	
		Print "timer phi_fvt_diff: ", (StopMSTimer(-2)-ttime)*10^-6	
		//MultiThread phi_fvt_re =  cos(	phi_fvt_diff	)//*x
		//MultiThread phi_fvt_im = sin(	phi_fvt_diff	)//*x
		MatrixOP/O/S phi_fvt_re = cos(	phi_fvt_diff	)
		MatrixOP/O/S phi_fvt_im = sin(	phi_fvt_diff	)
		Integrate/DIM=2 phi_fvt_re, phi_fvt_im
	
		Print "timer cos&sin&Int phi_fvt_im: ", (StopMSTimer(-2)-ttime)*10^-6		

		phi_fv_re = phi_fvt_re[p][q][tPnts]*velprobs[q]
		phi_fv_im = phi_fvt_im[p][q][tPnts]*velprobs[q]
		Integrate/DIM=1 phi_fv_re, phi_fv_im
	
		phi_f_re_temp = phi_fv_re[p][vPnts]
		phi_f_im_temp = phi_fv_im[p][vPnts]
	
		phi_f_re+=phi_f_re_temp
		phi_f_im+=phi_f_im_temp
	endfor

	phi_f_re/=numpnts(beampos)
	phi_f_im/=numpnts(beampos)
	
	phi_f = atan2(phi_f_im, phi_f_re)
	//phi_f = atan(phi_f_im/phi_f_re)
	contrast_f = sqrt(phi_f_im^2 + phi_f_re^2)
	
	Variable thetaShift = .5
	duplicate/o phi_f_im phi_f_imShift, phi_f_reShift, phi_f_shift, phi_f_unshift
	phi_f_imShift = phi_f_im*cos(thetaShift)+phi_f_re*sin(thetaShift)
	phi_f_reShift = phi_f_re*cos(thetaShift)-phi_f_im*sin(thetaShift)
	phi_f_shift = atan2(phi_f_imShift, phi_f_reShift)
	phi_f_unshift = phi_f_shift-thetashift
	
	//phi_f_unshift += 7
	//unwrap 2*pi, phi_f_unshift
	//phi_f_unshift -=7
	
	phi_f = phi_f_unshift
	
	contrast_f = sqrt(phi_f_imshift^2+phi_f_reshift^2)
	
	//contrast_f*=cos(phi_f)

	Print "timer total: ", (StopMSTimer(-2)-ttime)*10^-6	
	
	KillWaves phi_fvt_re, phi_fvt_im, tTest, phi_fvt_diff, phi_fv_re, phi_fv_im, phi_f_re, phi_f_im
End




Function ExtractChopperData(series_name, onekTwok, startPntsToKill, stopPntsToKill)
	String series_name
	Variable onektwok, startPntsToKIll, stopPntsToKill
	
	String onektwokStr = ""
	if(OnekTwok==1)
		onektwokStr = "1k2k"
	elseif(OnekTwok == 2)
		onektwokStr = "2k2k"
	endif
	
	String newphase = "phase" + onektwokStr + "_tofit_" + series_name
	String newphase_error = "phase_error" + onektwokStr + "_tofit_" + series_name
	String newcontrast = "contrast" + onektwokStr + "_tofit_" + series_name
	String newcontrast_error = "contrast_error" + onektwokStr + "_tofit_" + series_name
	String newavg_counts = "avg_counts" + onektwokStr + "_tofit_" + series_name
	String newavg_counts_error = "avg_counts_error" + onektwokStr + "_tofit_" + series_name
	String newfreq = "freq_tofit_" + series_name
	String newfileNumbers = "file_index_tofit_" + series_name
	
	String inputphase = "phase" + onektwokStr + "_" + series_name
	String inputphase_error = "phase_error" + onektwokStr + "_" + series_name
	String inputcontrast = "contrast" + onektwokStr + "_" + series_name
	String inputcontrast_error = "contrast_error" + onektwokStr + "_" + series_name
	String inputavg_counts = "avg_counts" + onektwokStr + "_" + series_name
	String inputavg_counts_error = "avg_counts_error" + onektwokStr + "_" + series_name
	String inputfreq = "freq" + "_" + series_name	
	String inputFileNumbers = "file_index_"+series_name
	
	Wave inphase = $inputphase
	Wave inphase_error = $inputphase_error
	Wave incontrast = $inputcontrast
	Wave incontrast_error = $inputcontrast_error
	Wave inavg_counts = $inputavg_counts
	Wave inavg_counts_error = $inputavg_counts_error
	Wave infreq = $inputfreq
	Wave fileNumbers = $inputFileNumbers
	
	Duplicate/o inphase $newphase
	Duplicate/o inphase_error $newphase_error
	Duplicate/o incontrast $newcontrast
	Duplicate/o incontrast_error $newcontrast_error
	Duplicate/o inavg_counts $newavg_counts
	Duplicate/o inavg_counts_error $newavg_counts_error
	Duplicate/o infreq $newfreq
	Duplicate/o fileNumbers $newfileNumbers
	
	
	Wave phase_saved = $newphase
	Wave phase_error_saved = $newphase_error
	Wave contrast_saved = $newcontrast
	Wave contrast_error_saved = $newcontrast_error
	Wave avg_counts_saved = $newavg_counts
	Wave avg_counts_error_saved = $newavg_counts_error
	Wave freq_saved = $newfreq
	Wave files_saved = $newfileNumbers
	
	
	String killFilesName = "KillFiles_"+series_name
	Wave/Z KillFiles = $killFilesName
	
	If(waveexists(killfiles))
		Variable i
		Variable killpnt
		For(i=0; i<numpnts(killFiles); i+=1)
			killpnt = killFiles[i]-1
			contrast_saved[killpnt]=nan
		EndFor
	endif
	
//	if(ParamIsDefault(otherFilesToKill)==0)
//		Variable killindex
//		Variable NumPointsToKill = numpnts(otherFilesToKill)
//		Variable i
//		for(i=0; i<numPointsToKill; i+=1)
//			Killindex = otherFilesToKill[i]
//			contrast_saved[Killindex] = NaN
//		endfor
//	endif
	
	// Generate a manual tick wave for file numbers
	String freqTicksWaveName = NameOfWave(infreq) + "_chopfit_ticks"
	String fileNumbersTicksWaveName = NameOfWave(filenumbers) + "_chopfit_ticks"
	Extract/O FileNumbers, filenumbersTemp, mod(p,5)==0
	Extract/O freq_saved, freqticksTemp, mod(p,5)==0
	Duplicate/O freqTicksTemp $freqTicksWaveName
	Make/O/T/N=(numpnts(filenumbersTemp)) $fileNumbersTicksWaveName = num2str(filenumbersTemp)
	Wave FileNumbersTicks = $fileNumbersTicksWaveName
	Wave freqTicks = $freqTicksWaveName
	
	
	DeletePoints 0, startPntsToKill, phase_saved, phase_error_saved, contrast_saved, contrast_error_saved, avg_counts_saved, avg_counts_error_saved, freq_saved, files_saved
	
	Variable pnts = numpnts(phase_saved)
	
	DeletePoints pnts-stopPntsToKill, stopPntsToKill, phase_saved, phase_error_saved, contrast_saved, contrast_error_saved, avg_counts_saved, avg_counts_error_saved, freq_saved, files_saved
	
	RemoveNaNsXYZABCDE(phase_saved, phase_error_saved, contrast_saved, contrast_error_saved, avg_counts_saved, avg_counts_error_saved, freq_saved,files_saved)
	
	Display/K=1/W=(40,50,1000,600) contrast_saved vs freq_saved
	ErrorBars $newcontrast Y,wave=(contrast_error_saved,contrast_error_saved)
	Label bottom "frequency " + series_name
	Label left "contrast" + onektwokStr + " " + series_name
	AppendToGraph/l=phase/T phase_saved vs freq_saved
	ModifyGraph userticks(top)={FreqTicks, FileNumbersTicks}
	Label top "file number " + series_name
	ModifyGraph axisEnab(left)={.52,1},axisEnab(phase)={0,0.48},freePos(phase)=0
	Label phase "phase" + onektwokStr + " " + series_name
	SetAxis phase 0,6.29
	SetAxis left 0,.35
	ErrorBars $newphase Y,wave=(phase_error_saved,phase_error_saved)
	ErrorBars $newcontrast Y,wave=(contrast_error_saved,contrast_error_saved)
	ModifyGraph mode=3,marker=8
	ModifyGraph rgb($newcontrast)=(0,0,0),rgb($newphase)=(1,4,52428)
	ModifyGraph lblPosMode(phase)=1
	ModifyGraph lblPosMode(left)=1
	ModifyGraph gfSize=10
	ModifyGraph gmSize=4
	ModifyGraph nticks(bottom)=20,minor(bottom)=1
	ModifyGraph grid(bottom)=2
	ShowInfo
	WindowNamer("Chopper scan " + onektwokStr + " " + series_name)
End



Function FindKillables(sigmas)
	Variable sigmas
	
	String traces = TraceNameList("", ";", 1)							// Get all traces on the graph
	String resTraceName = StringFromList(0,GrepList(traces, "Res"))	// Find the name of the trace that starts with "Res"
	
	Wave resWave = $resTraceName
	
	wavestats/q resWave
	
	Duplicate/O resWave killables
	
	killables = abs(resWave) > sigmas*v_sdev
	
	String/G autoKill = ""
	String/G autoKillc = ""
	String/G autoKillfiles = ""
	String/G autoKillfilesSemi = ""
	
	Variable i
	for(i=0; i<numpnts(killables); i+=1)
		if(killables[i] == 1)
			//print i
			autoKill += num2str(i)+";"
			autoKillc += num2str(i)+","
			autoKillfiles += num2str(i+1)+"," 
			autoKillfilesSemi += num2str(i+1)+";" 
		endif
	endfor
	
	//print autokill
	print "points: ", autokillc
	print "files: ", autokillfiles
	print "files: ", autokillfilesSemi
End




Function FindKillablesChop(series_name, sigmas)
	String series_name
	Variable sigmas
	
	String traces = TraceNameList("", ";", 1)							// Get all traces on the graph
	String resTraceName = StringFromList(0,GrepList(traces, "Res"))	// Find the name of the trace that starts with "Res"
	
	Wave resWave = $resTraceName
	
	wavestats/q resWave
	
	Duplicate/O resWave killables
	
	killables = abs(resWave) > sigmas*v_sdev
	
	String/G autoKill = ""
	String/G autoKillc = ""
	String/G autoKillfiles = ""
	
	String filesWaveName = "file_index_tofit_"+series_name
	Wave files = $filesWaveName
	
	Variable i
	for(i=0; i<numpnts(killables); i+=1)
		if(killables[i] == 1)
			//print i
			autoKill += num2str(i)+";"
			autoKillc += num2str(i)+","
			autoKillfiles += num2str(files[i])+"," 
		endif
	endfor
	
	//print autokill
	print "points: ", autokillc
	print "files: ", autokillfiles
End


Function GoFitPhaseChoppers(series_name, onektwok)
	String series_name
	Variable onekTwok

	String onektwokStr = ""
	if(OnekTwok==1)
		onektwokStr = "1k2k"
	elseif(OnekTwok == 2)
		onektwokStr = "2k2k"
	endif

	String inputcontrast = "contrast"+onektwokStr+"_tofit_" + series_name
	String inputcontrast_error = "contrast_error"+onektwokStr+"_tofit_" + series_name
	String inputfreq = "freq_tofit_" + series_name	
	
	Wave contrastWave = $inputcontrast
	Wave errWave = $inputcontrast_error
	Wave frequencyWave = $inputfreq
	
	Variable ttime= stopMSTimer(-2)
	
	Duplicate/O InitialGuesses w_coef, eps
	//eps[0] = 1e-0; eps[1]=0.5
	eps = 1e-6
	eps[2]=.1
	
	Variable/G fPnts = numpnts(frequencyWave)
	MakeChopperFitHelperWavesV3(frequencyWave)
//	MakeChopperFitHelperWavesTrans(frequencyWave)

	String HoldString = 	ProcessHoldAndConstraintsWaves()
	
	//Wave mask_k
	
	FuncFit/M=2/H=HoldString PhaseChoppers2IFMfitv3 W_coef  contrastWave /X=frequencyWave /R /I=1 /W=errWave /C=Constraints ///M=mask_k// /E=eps
	//	FuncFit/N/M=2/L=400/ODR=0/H="0000" PhaseChoppers2IFMfitv2 W_coef  ywave /X=xwave /R /D /I=1 /W=errwave 
	//Add /N to turn off updates
	//Add /O to only graph initial guesses
	
	Wave M_Covar
	Duplicate/O M_Covar, CorMat	 
	CorMat = M_Covar[p][q]/sqrt(M_Covar[p][p]*M_Covar[q][q])
	
	Variable/G reducedChiSqrd = V_chisq/(V_npnts-numpnts(W_coef))
	printf "Reduced Chi-Squared: %g\r", reducedChiSqrd
	//TextBox/C/N=text1/Z=1/X=85.00/Y=30.00 "\\F'Symbol'c\\F'Geneva'\\S2\\M/dof = "+num2str(reducedChiSqrd)
	
	//ModifyGraph zero(Res_Left)=1
	
	Print "timer total v2: ", (StopMSTimer(-2)-ttime)*10^-6
	
	CalculateChopperBestFit(w_coef, frequencyWave, series_name+onektwokStr)
	
	String bestFitContrastStr = "contrastFit_"+series_name+onektwokStr
	String bestFitPhaseStr = "phaseFit_"+series_name+onektwokStr
	
	Wave bestFitContrastWave = $bestfitContrastStr
	Wave bestFitPhaseWave = $bestfitPhaseStr
	
	If(StringMatch(TraceInfo("", bestFitContrastStr,0),""))		// if the trace is not already on the graph then put it there
		AppendToGraph bestFitContrastWave
		//ModifyGraph rgb($fitPhaseWaveName)=(1,16019,65535)
		ModifyGraph zero(Res_Left)=1
		SetAxis/A Res_Left
		
		AppendToGraph/l=phase bestFitPhaseWave
	EndIf
	
	CleanupChopperFitWaves()
End



Function FindVpitester()
	make/o/n=100 vratiowave = p+5, rootwave, vpiwave,vpipercentwave

	variable i
	for(i=0;i<numpnts(vratiowave);i+=1)
		FindVpi(3000,vratiowave[i],pi)
		nvar theroot, vpi, vpipercent
		rootwave[i]=theroot
		vpiwave[i]=vpi
		vpipercentwave[i]=vpipercent
	endfor
End

Function FindVpi(FlowV, vRatio, phiMeasured)
	Variable FlowV, vRatio, phiMeasured
	
	Wave velwave, velprobs//, w_coef, InitialGuesses
	
	//Variable FlowV =InitialGuesses[0]
	//Variable vratio = InitialGuesses[1]
	Variable sigV = FlowV/vratio
	Variable velnorm
	
	
	Variable vCubed = 0
	if(vCubed)
		Velprobs = exp(-(x-FlowV)^2/(2*SigV^2))*x^3
		VelNorm = sqrt(2*pi)*sigv*FlowV*(FlowV^2+3*sigv^2)
	else
		Velprobs = exp(-(x-FlowV)^2/(2*SigV^2))
		VelNorm = sqrt(2*pi*SigV^2)
	endif
	//Variable velsum = sum(velprobs,-inf,inf)
	Velprobs/=velNorm
	
	
	Duplicate/O velprobs FofV, FofVc, phiIm, phiRe
	
	NVAR centerpos, b//, phiMeasured
	
	variable k = hbar*2*pi/gratingperiod/(kavgmass*amu2kg)*11*25.4e-3
	
	FofV = (2*centerpos*k/x-(k/x)^2)/((b^2-centerpos^2)*(b^2-(centerpos-k/x)^2))	/x	// x is the scaled index i.e. velocity
	Variable phiFlowVoverC = (2*centerpos*k/FlowV-(k/FlowV)^2)/((b^2-centerpos^2)*(b^2-(centerpos-k/FlowV)^2))	/ FlowV
	
	Make/d/o PhiRooterCoef = {phiMeasured,phiFlowVoverC}
	
	FindRoots/Q/L=(.0001)/H=(.001) PhiRooter, PhiRooterCoef
	
	Variable/g theRoot = V_root
	FofVc = FofV*theRoot
	FindLevel/Q FofVc, phiMeasured
	Variable/g vpi = v_levelx
	Variable/g vpipercent = (vpi-flowv)/flowv*100
	
	//print "the root: ", theRoot; print "vpi: ", vpi; print "% difference in vpi: ", vpipercent
End




Function PhiRooter(w,x)
	Wave w
	Variable x
	
	wave velprobs, FofV, phiIm, phiRe
	
	Variable phiMeasured = w[0]
	Variable phiFlowVoverC = w[1]
	
	Variable phiFlowV = phiFlowVoverC*x
	
	phiIm = velprobs*sin(x*FofV)
	phiRe = velprobs * cos(x*FofV)
	
	Variable imPart = area(phiIm)
	Variable rePart = area(phiRe)
	
	Variable phi = atan2(impart,repart)	// measured phase shift a beam with a given velocity distribution
	//print phi
	Variable phiDiff = phi-phiMeasured
	
	Variable phiUnwrapDelta = phiflowV-phi

	Variable phiUnwrapped = phi + 2*pi*round((phiFlowV-phi)/(2*pi))
	
	return phiUnWrapped-phiMeasured//+2*pi
	//return phiUnwrapped
End



Function MakeChopperFitHelperWavesV4(freq)
	Wave freq
	
	NVAR fPnts
	
	Wave InitialGuesses
	Variable FlowV = InitialGuesses[0]
	Variable SigV = FlowV/InitialGuesses[1]
	Variable phi12avg = InitialGuesses[3]
	Variable phi12diff = InitialGuesses[4]

	Variable/G vPnts = 100
	Variable nSigmas = 5
	
	Variable/G tpnts = 500
	//Variable dutycycle = 1/2
	
	Variable posPnts = 2
	Variable/G centerPos = .8e-3
	Variable beamWidth = 50e-6
	Variable minPos = centerPos-beamWidth/2
	Variable maxPos = centerPos+beamWidth/2

	Variable a = 1e-3	// ground plane - wire dist
	Variable radius = 0.785e-3
	Variable/G b = a*sqrt(1+2*radius/a)

	Variable minVel = FlowV*(1-nSigmas/15)//; minVel = 0
	Variable maxVel = FlowV*(1+nSigmas/15)//; maxVel = 5000
	Variable velStep = (maxVel-minVel)/vPnts
	
	Make/D/O/N=(vPnts) velwave, velprobs
	//	Duplicate/O velwave velprobs 
	SetScale/I x minVel,maxVel,"m/s", velwave, velprobs
	velwave = x
	
	Velprobs = exp(-(x-FlowV)^2/(2*SigV^2))*x^3
	Variable VelNorm = sqrt(2*pi)*sigv*FlowV*(FlowV^2+3*sigv^2)
	//Variable velsum = sum(velprobs,-inf,inf)
	Velprobs/=velNorm

	Make/D/O/N=(fPnts, vPnts, tPnts) phi_fvt_re, phi_fvt_im, phi_fvt_diff, foftsqrdexpChop1, foftsqrdexpChop2, sagPhaseWave//,ttest
	Make/D/O/N=(fPnts, vPnts) phi_fv_re, phi_fv_im
	Make/D/O/N=(fPnts) phi_f_re, phi_f_im, phi_f, con_f
	//SetScale/I x minChopFreq, maxChopFreq, "Hz", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_f_re, phi_f_im, phi_f, contrast_f, ttest, phi_fvt_diff
	SetScale/I y minVel, maxVel, "m/s", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_fvt_diff, foftsqrdexpChop1, foftsqrdexpChop2, sagPhaseWave//,ttest
	SetScale/P z 0, 1/tpnts, "t*f", phi_fvt_re, phi_fvt_im, phi_fvt_diff, foftsqrdexpChop1, foftsqrdexpChop2, sagPhaseWave
	
	Make/D/O/N=(posPnts) beamPos, beamChopProfile
	SetScale/i x minPos, maxPos, "m", beamPos
	beampos= numpnts(beampos) == 1 ? centerpos : x
	beamChopProfile = 1/numpnts(beampos)
	
	Duplicate/O phi_f_re, phi_f_re_temp, phi_f_im_temp
	phi_f_re=0;phi_f_im=0;phi_f_re_temp=0; phi_f_im_temp=0
	
	Variable phi1 =  .5 * (2*phi12avg + phi12diff)		// in units of pi
	Variable phi2 =  .5 * (2*phi12avg - phi12diff)	
	
	//MultiThread foftsqrdexpChop12 =( phi1*foftsqrdexp( z/freq[p] , freq[p])  - phi2* foftsqrdexp( z/freq[p] + chopL/y, freq[p])) * pi/y^2

	MultiThread foftsqrdexpChop1 = foftsqrdexp( z/freq[p] , freq[p]) / y^2
	MultiThread foftsqrdexpChop2 = foftsqrdexp( z/freq[p] + chopL/y, freq[p]) / y^2
	
	SagnacAndGravityPhaseChop()
End






Function MakeChopperFitHelperWavesV3(freq)
	Wave freq
	
	//NVAR fPnts
	Variable fPnts = numpnts(freq)
	
	Wave InitialGuesses
	Variable FlowV = InitialGuesses[0]
	Variable SigV = FlowV/InitialGuesses[1]
	Variable phi12avg = InitialGuesses[3]
	Variable phi12diff = InitialGuesses[4]

	Variable/G vPnts = 75	//was 100
	Variable nSigmas = 5
	
	Variable/G tpnts = 300	//was 500
	//Variable dutycycle = 1/2
	
	Variable posPnts = 1
	Variable/G centerPos = .8e-3					// distance from ground plane
	Variable beamWidth = 100e-6					// full width
	Variable minPos = centerPos-beamWidth/2
	Variable maxPos = centerPos+beamWidth/2

	Variable a = 1e-3								// ground plane to wire dist
	Variable radius = 0.785e-3						// wire radius
	Variable/G b = a*sqrt(1+2*radius/a)

	Variable minVel = FlowV*(1-nSigmas/15)//; minVel = 0
	Variable maxVel = FlowV*(1+nSigmas/15)//; maxVel = 5000
	Variable velStep = (maxVel-minVel)/vPnts
	
	Make/D/O/N=(vPnts) velwave, velprobs
	//	Duplicate/O velwave velprobs 
	SetScale/I x minVel,maxVel,"m/s", velwave, velprobs
	velwave = x
	
	Velprobs = exp(-(x-FlowV)^2/(2*SigV^2))*x^3
	Variable VelNorm = sqrt(2*pi)*sigv*FlowV*(FlowV^2+3*sigv^2)
	//Variable velsum = sum(velprobs,-inf,inf)
	Velprobs/=velNorm

	Make/D/O/N=(fPnts, vPnts, tPnts) phi_fvt_re, phi_fvt_im, phi_fvt_diff, foftsqrdexpChop1, foftsqrdexpChop2, sagPhaseWave//,ttest
	Make/D/O/N=(fPnts, vPnts) phi_fv_re, phi_fv_im
	Make/D/O/N=(fPnts) phi_f_re, phi_f_im, phi_f, con_f
	//SetScale/I x minChopFreq, maxChopFreq, "Hz", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_f_re, phi_f_im, phi_f, contrast_f, ttest, phi_fvt_diff
	SetScale/I y minVel, maxVel, "m/s", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_fvt_diff, foftsqrdexpChop1, foftsqrdexpChop2, sagPhaseWave//,ttest
	SetScale/P z 0, 1/tpnts, "t*f", phi_fvt_re, phi_fvt_im, phi_fvt_diff, foftsqrdexpChop1, foftsqrdexpChop2, sagPhaseWave
	
	Make/D/O/N=(posPnts) beamPos, beamChopProfile, deltaPos
	SetScale/i x minPos, maxPos, "m", beamPos
	beampos= numpnts(beampos) == 1 ? centerpos : x
	beamChopProfile = 1/numpnts(beampos)
	SetScale/i x -beamWidth/2, beamWidth/2, deltaPos
	deltaPos = x
	
	Duplicate/O phi_f_re, phi_f_re_temp, phi_f_im_temp
	phi_f_re=0;phi_f_im=0;phi_f_re_temp=0; phi_f_im_temp=0
	
	Variable phi1 =  .5 * (2*phi12avg + phi12diff)		// in units of pi
	Variable phi2 =  .5 * (2*phi12avg - phi12diff)	
	
	//MultiThread foftsqrdexpChop12 =( phi1*foftsqrdexp( z/freq[p] , freq[p])  - phi2* foftsqrdexp( z/freq[p] + chopL/y, freq[p])) * pi/y^2
	
	MultiThread foftsqrdexpChop1 = foftsqrdexp( z/freq[p] , freq[p]) / y^2
	MultiThread foftsqrdexpChop2 = foftsqrdexp( z/freq[p] + chopL/y, freq[p]) / y^2
	
	SagnacAndGravityPhaseChop()
End






Function vAvgTester()
	variable vpnts=200
	variable nsigmas = 5
	
	Variable FlowV = 3000
	Variable Vratio = 20
	Variable Vavg = 3000
	
	Variable minVel = FlowV*(1-nSigmas/15)//; minVel = 0
	Variable maxVel = FlowV*(1+nSigmas/15)//; maxVel = 5000
	Variable velStep = (maxVel-minVel)/vPnts
	
	Make/D/O/N=(vPnts) velwave, velprobs
	//	Duplicate/O velwave velprobs 
	SetScale/I x minVel,maxVel,"m/s", velwave, velprobs
	velwave = x
	
	duplicate/o velprobs velprobsVavg, velint
	
	variable velnorm

	
	velprobs = 	exp(-vratio^2/2*(x/FlowV-1)^2)*x^3
	velnorm = area(velprobs)
	velprobs/=velnorm
	//VelNorm = sqrt(2*pi)*sigv*FlowV*(FlowV^2+3*sigv^2)
	velprobsVavg = exp(-vratio^2/2* ( x/vavg * (1+1/vratio^2+2/(3+vratio^2)) -1 )^2 )
	velnorm = area(velprobsVavg)
	velprobsVavg/=velnorm
	
	velint = velprobs*x
	Integrate velint
	print "velprobs avg: ", velint[vpnts]
	
	velint = velprobsvavg*x
	integrate velint
	print "velprobsvavg: ", velint[vpnts]
End
	


Function MakeChopperFitHelperWavesTrans(freq)
	Wave freq
	
	//NVAR fPnts
	variable fPnts = numpnts(freq)
	Variable/G vPnts = 100	//was 100
	Variable/G tpnts = 500	//was 500
	Variable/G zPnts = 100
	
	variable numPerFreqChunk = 20
	variable numFreqChunks = floor(numpnts(freq)/numPerFreqChunk)
	variable extraFreqs = mod(numpnts(freq),numPerFreqChunk)
	variable startp, endp
	variable i
	
	
	Wave InitialGuesses
	Variable FlowV = InitialGuesses[0]
	Variable SigV = FlowV/InitialGuesses[1]
	Variable phi12avg = InitialGuesses[3]
	Variable phi12diff = InitialGuesses[4]


	Variable nSigmas = 5
	
	//Variable dutycycle = 1/2
	
	Variable posPnts = 1
	Variable/G centerPos = .8e-3					// distance from ground plane
	Variable beamWidth = 200e-6					// full width
	Variable minPos = centerPos-beamWidth/2
	Variable maxPos = centerPos+beamWidth/2

	Variable a = 1e-3								// ground plane to wire dist
	Variable radius = 0.785e-3						// wire radius
	Variable/G b = a*sqrt(1+2*radius/a)

	Variable minVel = FlowV*(1-nSigmas/15)//; minVel = 0
	Variable maxVel = FlowV*(1+nSigmas/15)//; maxVel = 5000
	Variable velStep = (maxVel-minVel)/vPnts
	
	Make/D/O/N=(vPnts) velwave, velprobs
	//	Duplicate/O velwave velprobs 
	SetScale/I x minVel,maxVel,"m/s", velwave, velprobs
	velwave = x
	
	Velprobs = exp(-(x-FlowV)^2/(2*SigV^2))*x^3
	Variable VelNorm = sqrt(2*pi)*sigv*FlowV*(FlowV^2+3*sigv^2)
	//Variable velsum = sum(velprobs,-inf,inf)
	Velprobs/=velNorm
	

	Make/D/O/N=(numPerFreqChunk, vPnts, tPnts, zPnts) voltageOn
	Make/D/O/N=(fPnts, vPnts, tPnts) phi_fvt_re, phi_fvt_im, phi_fvt_diff, foftsqrdexpChop1, foftsqrdexpChop2, sagPhaseWave//,ttest
	Make/D/O/N=(fPnts, vPnts) phi_fv_re, phi_fv_im
	Make/D/O/N=(fPnts) phi_f_re, phi_f_im, phi_f, con_f
	//SetScale/I x minChopFreq, maxChopFreq, "Hz", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_f_re, phi_f_im, phi_f, contrast_f, ttest, phi_fvt_diff
	SetScale/I y minVel, maxVel, "m/s", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_fvt_diff, foftsqrdexpChop1, foftsqrdexpChop2, sagPhaseWave, voltageOn//,ttest
	SetScale/P z 0, 1/tpnts, "t*f", phi_fvt_re, phi_fvt_im, phi_fvt_diff, foftsqrdexpChop1, foftsqrdexpChop2, sagPhaseWave, voltageOn
	SetScale/i t -lc/2, lc/2, "m" voltageOn
	
	Make/D/O/N=(posPnts) beamPos, beamChopProfile, deltaPos
	SetScale/i x minPos, maxPos, "m", beamPos
	beampos= numpnts(beampos) == 1 ? centerpos : x
	beamChopProfile = 1/numpnts(beampos)
	SetScale/i x -beamWidth/2, beamWidth/2, deltaPos
	deltaPos = x
	
	Duplicate/O phi_f_re, phi_f_re_temp, phi_f_im_temp
	phi_f_re=0;phi_f_im=0;phi_f_re_temp=0; phi_f_im_temp=0
	
	Variable phi1 =  .5 * (2*phi12avg + phi12diff)		// in units of pi
	Variable phi2 =  .5 * (2*phi12avg - phi12diff)	
	

	Variable x0 = 0.75e-3	//approximate beam location
	Variable sep = 0.01e-3	//approximate beam separation
	
	Make/O/D/N=(zPnts) esqrdofZ
	SetScale/i x -lc/2, lc/2, "m" esqrdofZ
	
	esqrdofZ = ( ((x0+b)^2 + x^2)* ((x0-b)^2+x^2) )^-1 - ( ((x0-sep+b)^2 + x^2)* ((x0-sep-b)^2+x^2) )^-1
	
	Duplicate/O esqrdOfZ esqrdofZint
	Integrate esqrdofZint
	Variable/G intnorm = esqrdofZint[zPnts]
	esqrdofZint/=intnorm
	esqrdofZ/=intnorm
	
	

	for(i=0; i<=numFreqChunks; i+=1)
		if(i<numFreqChunks)
			startp = i*numPerFreqChunk
			endp = (i+1)*numPerFreqChunk-1
			Duplicate/O/R=[startp, endp] freq freqChunk
		elseif(i==numFreqChunks && extraFreqs!=0)
			startp = i*numPerFreqChunk
			endp = startp+extraFreqs-1
			Duplicate/O/R=[startp, endp] freq freqChunk
			Redimension/N=(extraFreqs, vPnts, tPnts, zPnts) voltageOn
		else
			break
		endif
		
		Multithread voltageOn = mod((t+lc/2+z*y/freqChunk[p]), y/freqChunk[p]) < y/freqChunk[p]/2 ? 1 : 0	// t starts at -lc/2
		Multithread voltageOn*=esqrdofZ[s]
		Integrate/DIM=3 voltageOn
		foftsqrdexpChop1[startp, endp][][] = voltageOn[p-startp][q][r][zpnts] / y^2
	
		//	MultiThread foftsqrdexpChop2 = foftsqrdexp( z/freq[p] + chopL/y, freq[p]) / y^2
		Multithread voltageOn = mod((t+lc/2+z/freqChunk[p]*y+chopL), y/freqChunk[p]) < y/freqChunk[p]/2 ? 1 : 0	// t starts at -lc/2
		Multithread voltageOn*=esqrdofZ[s]
		Integrate/DIM=3 voltageOn	
		foftsqrdexpChop2[startp, endp][][] = voltageOn[p-startp][q][r][zpnts] / y^2	
	endfor	

	SagnacAndGravityPhaseChop()
End


function testp()
	make/o/n=20 testw
	make/o/n=(40,20) test2d = q*p
	make/o/n=(20,20,20) test3d = q*p*r
	testw[0,9]=test2d[p][20]
	testw[10,19]=test2d[p][20]
	test2d[0,19][]=test3d[p][q][20]
	test2d[20,39][]=test3d[p-20][q][20]
end


//		Multithread voltageOn = mod((t+lc/2+z*y/freq[p]), y/freq[p]) < y/freq[p]/2 ? 1 : 0	// t starts at -lc/2
//		Multithread voltageOn*=esqrdofZ[s]
//		Integrate/DIM=3 voltageOn
//		foftsqrdexpChop1 = voltageOn[p][q][r][zpnts] / y^2
//	
//		//	MultiThread foftsqrdexpChop2 = foftsqrdexp( z/freq[p] + chopL/y, freq[p]) / y^2
//		Multithread voltageOn = mod((t+lc/2+z/freq[p]*y+chopL), y/freq[p]) < y/freq[p]/2 ? 1 : 0	// t starts at -lc/2
//		Multithread voltageOn*=esqrdofZ[s]
//		Integrate/DIM=3 voltageOn	
//		foftsqrdexpChop2 = voltageOn[p][q][r][zpnts] / y^2	



Function test()
	variable pnts1=10
	variable pnts2=10
	variable pnts3=10
	variable pnts4=10
	
	make/o/n=(pnts1,pnts2,pnts3,pnts4) fourd = s*p*q*r
	make/o/n=(pnts1,pnts2,pnts3) threed
	
	duplicate/o fourd fourdint
	integrate/dim=3 fourdint
	
	threed=fourdint[p][q][r][pnts4]
end




Function MakeChopperFitHelperWaves(freq)
	Wave freq
	
	NVAR fPnts
	
	Wave InitialGuesses
	Variable FlowV = InitialGuesses[0]
	Variable phi12avg = InitialGuesses[3]
	Variable phi12diff = InitialGuesses[4]

	Variable/G vPnts = 50
	Variable nSigmas = 5
	
	Variable/G tpnts = 500
	//Variable dutycycle = 1/2
	
	Variable posPnts = 1
	Variable/G centerPos = .8e-3
	Variable beamWidth = 50e-6
	Variable minPos = centerPos-beamWidth/2
	Variable maxPos = centerPos+beamWidth/2

	Variable a = 1e-3	// ground plane - wire dist
	Variable radius = 0.785e-3
	Variable/G b = a*sqrt(1+2*radius/a)

	Variable minVel = FlowV*(1-nSigmas/15)//; minVel = 0
	Variable maxVel = FlowV*(1+nSigmas/15)//; maxVel = 5000
	Variable velStep = (maxVel-minVel)/vPnts
	
	Make/D/O/N=(vPnts) velwave = minVel+p*velStep
	Duplicate/O velwave velprobs 
	SetScale/I x minVel,maxVel,"m/s", velprobs
	


		
	Make/D/O/N=(fPnts, vPnts, tPnts) phi_fvt_re, phi_fvt_im, phi_fvt_diff, foftsqrdexpChop12, sagPhaseWave//,ttest
	Make/D/O/N=(fPnts, vPnts) phi_fv_re, phi_fv_im
	Make/D/O/N=(fPnts) phi_f_re, phi_f_im, phi_f, con_f
	//SetScale/I x minChopFreq, maxChopFreq, "Hz", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_f_re, phi_f_im, phi_f, contrast_f, ttest, phi_fvt_diff
	SetScale/I y minVel, maxVel, "m/s", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_fvt_diff, foftsqrdexpChop12, sagPhaseWave//,ttest
	SetScale/P z 0, 1/tpnts, "t*f", phi_fvt_re, phi_fvt_im, phi_fvt_diff, foftsqrdexpChop12, sagPhaseWave
	
	Make/D/O/N=(posPnts) beamPos
	SetScale/i x minPos, maxPos, "m", beamPos
	beampos=x
	beampos=centerpos
	Duplicate/O phi_f_re, phi_f_re_temp, phi_f_im_temp
	phi_f_re=0;phi_f_im=0;phi_f_re_temp=0; phi_f_im_temp=0
	
	Variable phi1 =  .5 * (2*phi12avg + phi12diff)		// in units of pi
	Variable phi2 =  .5 * (2*phi12avg - phi12diff)	
	
	MultiThread foftsqrdexpChop12 =( phi1*foftsqrdexp( z/freq[p] , freq[p])  - phi2* foftsqrdexp( z/freq[p] + chopL/y, freq[p])) * pi/y^2
End






Function PhaseChoppers2IFMfitv3(pw, yw, freq) : FitFunc
	Wave pw, yw, freq
	
	// pw0 = flow velocity
	// pw1 = vratio
	// pw2 = phi12err
	// pw3 = contrast normalization constant
	// pw4 = delta phase
	
	// xw is the frequency
	// yw will be the contrast
	
	Variable ttime= stopMSTimer(-2)
	
	
	Variable FlowV =pw[0]
	Variable vratio = pw[1]
	Variable c0 = pw[2]
	Variable phi12avg = pw[3]
	Variable phi12diff = abs(pw[4])
	
	Variable sigV = FlowV/vratio
	Variable VelNorm
	
	//Wave HoldWave
	
	Wave Velprobs, phi_fvt_diff, beampos, beamChopProfile, phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_f_re, phi_f_im, phi_f, con_f, phi_f_re_temp, phi_f_im_temp, foftsqrdexpChop1, foftsqrdexpChop2, phiSagGravChop, deltaPos
	NVAR b, centerpos, tPnts, vPnts, fPnts, AvgSagGravPhi, AvgSagGravCon
	
	//Velprobs = exp(-(x-FlowV)^2/(2*SigV^2))*x^3
	Variable vCubed = 0
	if(vCubed)
		Velprobs = exp(-vratio^2/2*(x/FlowV-1)^2)*x^3
		VelNorm = sqrt(2*pi)*sigv*FlowV*(FlowV^2+3*sigv^2)
	else
		Velprobs = exp(-(x-FlowV)^2/(2*SigV^2))
		VelNorm = sqrt(2*pi*SigV^2)
	endif
	//Velprobs = 1/sqrt(2*pi*sigv)*exp(-(x-FlowV)^2/(2*SigV^2))
	//Variable velsum = sum(velprobs,-inf,inf)
	Velprobs/=velNorm
	
	SagnacAndGravityPhaseChop()		// Since the average sagnac phase shift and contrast changes a little bit with different velprobs, this should technically be run each time. in practice, it makes little to no difference for reasonable vratios and initial guesses
	
	phi_fvt_diff=0
	phi_fvt_re=0
	phi_fvt_im=0
	phi_f_re = 0
	phi_f_im = 0	
	
	//Print "timer assign vars, make waves: ", (StopMSTimer(-2)-ttime)*10^-6	
	
	Variable phi1VpiSqrd
	Variable phi2VpiSqrd
	
	// turn on for vpi correction
	if(1)
		Variable phi1 = (phi12avg + phi12diff/2)	 * pi
		Variable phi2 = (phi12avg - phi12diff/2)	 * pi
		FindVpi(FlowV, vRatio, phi1)
		NVAR vpi
		phi1VpiSqrd =  phi1*vpi^2
		FindVpi(FlowV, vRatio, phi2)
		phi2VpiSqrd =  phi2*vpi^2
	else
		phi1VpiSqrd =  (phi12avg + phi12diff/2)	 * pi *flowv^2
		phi2VpiSqrd =  (phi12avg - phi12diff/2)	 * pi *flowv^2
	endif
	
	//	Variable flowvPhi12err = flowv^2*phi12avg
	
	Variable sagSign = 1
	if(numpnts(beampos)==0 || numpnts(beampos)==1)
		MatrixOP/O/S phi_fvt_diff = phi1VpiSqrd*foftsqrdexpChop1 - phi2VpiSqrd*foftsqrdexpChop2 + sagSign*phiSagGravChop		// Still 5 times faster than a multithreaded wave assignment
		MatrixOP/O/S phi_fvt_re = cos(	phi_fvt_diff	)
		MatrixOP/O/S phi_fvt_im = sin(	phi_fvt_diff	)
	else
		//ttime= stopMSTimer(-2)
		phi_fvt_diff = 0
		Variable i, j
		Variable positionPhase1, positionPhase2
		Variable positionFactor
		Variable diffPos = centerPos+10e-6
		Variable c = (b^2-centerPos^2)*(b^2-diffPos^2)/(centerPos^2-diffPos^2)
		for(i=0; i<numpnts(beampos); i+=1)
			for(j=0; j<numpnts(beampos); j+=1)
				//		print "i= ", i, "; beampos[i] = ", beampos[i]
				//positionPhase = (b^2-centerPos^2)/(b^2-beampos[i]^2)
				positionPhase1 = c * ( (centerPos + deltapos[i] )^2 - (diffPos + deltapos[i])^2 ) / (b^2 - (centerPos+deltaPos[i])^2) / (b^2 - (diffPos+deltaPos[i])^2)
				positionPhase2 = c * ( (centerPos + deltapos[j] )^2 - (diffPos + deltapos[j])^2 ) / (b^2 - (centerPos+deltaPos[j])^2) / (b^2 - (diffPos+deltaPos[j])^2)
				positionFactor = beamChopProfile[i]*beamChopProfile[j] 
				//printf "positionPhase: %g		positionFactor: %g\r" positionPhase, positionFactor
				//MatrixOP/O/S phi_fvt_diff = positionFactor * ( phi1FlowV*foftsqrdexpChop1 - phi2flowV*foftsqrdexpChop2) +phi_fvt_diff 	// doesn't work because MatrixOP doesn't allow 3D waves to be both a source and a destination wave.
				MatrixOP/O/S phi_fvt_diff = ( positionPhase1 * phi1VpiSqrd*foftsqrdexpChop1 - positionPhase2 * phi2VpiSqrd*foftsqrdexpChop2)	+ sagSign*phiSagGravChop  // Multithread is slower than matrixOP by a factor of 5
				Multithread phi_fvt_re += positionFactor*cos(	phi_fvt_diff	)
				Multithread phi_fvt_im += positionFactor*sin(	phi_fvt_diff	)
				//Print "timer phi_fvt_diff: ", (StopMSTimer(-2)-ttime)*10^-6	
			endfor
		endfor
	endif
	

	Integrate/DIM=2 phi_fvt_re, phi_fvt_im
	
	Variable multiIFMyesno=0
	Variable secondIFMweight = 0.05
	if(multiIFMyesno==1)
		MatrixOP/O/S phi_fvt_diff = phi1VpiSqrd*foftsqrdexpChop1 - phi2VpiSqrd*foftsqrdexpChop2 - phiSagGravChop	
		MatrixOP/O/S phi_fvt_diffM = -1*phi1VpiSqrd*foftsqrdexpChop1 + phi2VpiSqrd*foftsqrdexpChop2 - phiSagGravChop	
		MatrixOP/O/S phi_fvt_diff2 = 2*phi1VpiSqrd*foftsqrdexpChop1 - 2*phi2VpiSqrd*foftsqrdexpChop2- phiSagGravChop	
		MatrixOP/O/S phi_fvt_re = (1-secondIFMweight)* cos(	phi_fvt_diff	)+ secondIFMweight*cos(phi_fvt_diff2)
		MatrixOP/O/S phi_fvt_im = (1-secondIFMweight)*sin(	phi_fvt_diff	)+ secondIFMweight*sin(phi_fvt_diff2)
		Integrate/DIM=2 phi_fvt_re, phi_fvt_im
	endif
	
	Variable molFrac = 0
	if(molFrac!=0)
		phi1VpiSqrd *= 5/6
		phi2VpiSqrd *= 5/6
		MatrixOP/O/S phi_fvt_diff_mol = phi1VpiSqrd*foftsqrdexpChop1 - phi2VpiSqrd*foftsqrdexpChop2 -phiSagGravChop
		MatrixOP/O/S phi_fvt_re =(1-molFrac)* cos(	phi_fvt_diff	)+ molFrac*cos(phi_fvt_diff_mol)
		MatrixOP/O/S phi_fvt_im =(1-molFrac)*  sin(	phi_fvt_diff	)+ molFrac*sin(phi_fvt_diff_mol)
		Integrate/DIM=2 phi_fvt_re, phi_fvt_im
	endif
	
	//Print "timer cos&sin&Int phi_fvt_im: ", (StopMSTimer(-2)-ttime)*10^-6		

	phi_fv_re = phi_fvt_re[p][q][tPnts]*velprobs[q]
	phi_fv_im = phi_fvt_im[p][q][tPnts]*velprobs[q]
	Integrate/DIM=1 phi_fv_re, phi_fv_im
	
	phi_f_re = phi_fv_re[p][vPnts]
	phi_f_im = phi_fv_im[p][vPnts]
	
	//phi_f_re+=phi_f_re_temp
	//phi_f_im+=phi_f_im_temp
	//endfor

	//phi_f_re/=numpnts(beampos)
	//phi_f_im/=numpnts(beampos)
	
	variable unwrapShift = pi+1
	phi_f = atan2(phi_f_im, phi_f_re)
	phi_f +=unwrapShift 
//	phi_f = mod(phi_f, 2*pi)	
	phi_f-=unwrapShift
	phi_f-=sagSign*AvgSagGravPhi
	
	con_f = sqrt(phi_f_im^2 + phi_f_re^2)
	
	yw = con_f*c0/AvgSagGravCon
	
	//Print "timer total v2: ", (StopMSTimer(-2)-ttime)*10^-6
End


Function PhaseChoppers2IFMfitv3mol(pw, yw, freq) : FitFunc
	Wave pw, yw, freq
	
	// pw0 = flow velocity
	// pw1 = vratio
	// pw2 = phi12err
	// pw3 = contrast normalization constant
	// pw4 = delta phase
	
	// xw is the frequency
	// yw will be the contrast
	
	//Variable ttime= stopMSTimer(-2)
	
	
	Variable FlowV =pw[0]
	Variable vratio = pw[1]
	Variable c0 = pw[2]
	Variable phi12avg = pw[3]
	Variable phi12diff = abs(pw[4])
	Variable mfrac = 0
	
	Variable sigV = FlowV/vratio
	Variable VelNorm
	
	Wave HoldWave
	
	Wave Velprobs, phi_fvt_diff, beampos, phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_f_re, phi_f_im, phi_f, con_f, phi_f_re_temp, phi_f_im_temp, foftsqrdexpChop1, foftsqrdexpChop2, phiSagGravChop
	NVAR b, centerpos, tPnts, vPnts, fPnts
	
	Velprobs = exp(-(x-FlowV)^2/(2*SigV^2))*x^3
	VelNorm = sqrt(2*pi)*sigv*FlowV*(FlowV^2+3*sigv^2)
	//Variable velsum = sum(velprobs,-inf,inf)
	Velprobs/=velNorm
	
	phi_f_re = 0
	phi_f_im = 0	
	
	//Print "timer assign vars, make waves: ", (StopMSTimer(-2)-ttime)*10^-6	
	
	Variable phi1FlowV =  .5 * (2*phi12avg + phi12diff)	 * pi *flowv^2
	Variable phi2FlowV =  .5 * (2*phi12avg - phi12diff)	 * pi *flowv^2

	//	Variable flowvPhi12err = flowv^2*phi12avg
	
	Wave phiSagGravChop
	MatrixOP/O/S phi_fvt_diff = phi1FlowV*foftsqrdexpChop1 - phi2flowV*foftsqrdexpChop2-phiSagGravChop
	
	phi1FlowV *= 5/6
	phi2FlowV *= 5/6
	MatrixOP/O/S phi_fvt_diff_mol = phi1FlowV*foftsqrdexpChop1 - phi2flowV*foftsqrdexpChop2 -phiSagGravChop

	//Variable i = 0
	//for(i=0; i<numpnts(beampos); i+=1)
	// foftsqrd( time(s), frequency(Hz) )
	//	MultiThread ttest = foftsqrd( z/x , x)
	//	MultiThread phi_fvt_diff = pi*(vpi/y)^2* ( foftsqrd( z/x , x) - foftsqrd( z/x + chopL/y, x) )	
	//print "i= ", i, "; beampos[i] = ", beampos[i]
	//MultiThread phi_fvt_diff = pi*(b^2-centerPos^2)/(b^2-beampos[i]^2)*(vpi/y)^2* ( foftsqrdexp( z/freq[p] , freq[p]) - foftsqrdexp( z/freq[p] + chopL/y, freq[p]) )	
	//MultiThread phi_fvt_diff = pi*(b^2-centerPos^2)/(b^2-beampos[i]^2)*(vpi/y)^2* ( foftsqrd( z/freq[p] , freq[p]) - deltaphi*foftsqrd( z/freq[p] + chopL/y, freq[p]) )	
	//	MultiThread phi_fvt_diff = pi*( foftsqrd( z/x , x) - foftsqrd( z/x + chopL/y, x) )	
	//	MultiThread phi_fvt_diff = pi*( foftsqrdexp( z/x , x) - foftsqrdexp( z/x + chopL/y, x) )	
	//	Print "timer phi_fvt_diff: ", (StopMSTimer(-2)-ttime)*10^-6	
	//MultiThread phi_fvt_re =  cos(	phi_fvt_diff	)
	//MultiThread phi_fvt_im = sin(	phi_fvt_diff	)
	MatrixOP/O/S phi_fvt_re = (1-mfrac) * cos(	phi_fvt_diff	) + mfrac * cos( phi_fvt_diff_mol )
	MatrixOP/O/S phi_fvt_im = (1-mfrac) * sin(	phi_fvt_diff	) + mfrac * sin(phi_fvt_diff_mol)
	Integrate/DIM=2 phi_fvt_re, phi_fvt_im
	
	//Print "timer cos&sin&Int phi_fvt_im: ", (StopMSTimer(-2)-ttime)*10^-6		

	phi_fv_re = phi_fvt_re[p][q][tPnts]*velprobs[q]
	phi_fv_im = phi_fvt_im[p][q][tPnts]*velprobs[q]
	Integrate/DIM=1 phi_fv_re, phi_fv_im
	
	phi_f_re = phi_fv_re[p][vPnts]
	phi_f_im = phi_fv_im[p][vPnts]
	
	//phi_f_re+=phi_f_re_temp
	//phi_f_im+=phi_f_im_temp
	//endfor

	//phi_f_re/=numpnts(beampos)
	//phi_f_im/=numpnts(beampos)
	
	phi_f = atan2(phi_f_im, phi_f_re)
	phi_f += pi+1
	phi_f = mod(phi_f, 2*pi)	
	phi_f-=pi+1
	
	con_f = sqrt(phi_f_im^2 + phi_f_re^2)
	
	yw = con_f*c0
	
	//Print "timer total v2: ", (StopMSTimer(-2)-ttime)*10^-6
End





Function PhaseChoppers2IFMfitv2(pw, yw, freq) : FitFunc
	Wave pw, yw, freq
	
	// pw0 = flow velocity
	// pw1 = vratio
	// pw2 = phi12err
	// pw3 = contrast normalization constant
	// pw4 = delta phase
	
	// xw is the frequency
	// yw will be the contrast
	
	//Variable ttime= stopMSTimer(-2)
	
	
	Variable FlowV =pw[0]
	Variable vratio = pw[1]
	Variable sigV = FlowV/vratio
	Variable VelNorm
	
	Variable phi12err = pw[2]
	Variable deltaphi = pw[4]
	
	Wave HoldWave
	
	Wave Velprobs, phi_fvt_diff, beampos, phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_f_re, phi_f_im, phi_f, con_f, phi_f_re_temp, phi_f_im_temp, foftsqrdexpChop12
	NVAR b, centerpos, tPnts, vPnts, fPnts
	
	Velprobs = exp(-(x-FlowV)^2/(2*SigV^2))*x^3
	VelNorm = sqrt(2*pi)*sigv*FlowV*(FlowV^2+3*sigv^2)
	//Variable velsum = sum(velprobs,-inf,inf)
	Velprobs/=velNorm
	
	phi_f_re = 0
	phi_f_im = 0	
	
	//Print "timer assign vars, make waves: ", (StopMSTimer(-2)-ttime)*10^-6	

	Variable flowvPhi12err = flowv^2*phi12err
	MatrixOP/O/S phi_fvt_diff = flowvPhi12err*foftsqrdexpChop12

	//Variable i = 0
	//for(i=0; i<numpnts(beampos); i+=1)
	// foftsqrd( time(s), frequency(Hz) )
	//	MultiThread ttest = foftsqrd( z/x , x)
	//	MultiThread phi_fvt_diff = pi*(vpi/y)^2* ( foftsqrd( z/x , x) - foftsqrd( z/x + chopL/y, x) )	
	//print "i= ", i, "; beampos[i] = ", beampos[i]
	//MultiThread phi_fvt_diff = pi*(b^2-centerPos^2)/(b^2-beampos[i]^2)*(vpi/y)^2* ( foftsqrdexp( z/freq[p] , freq[p]) - foftsqrdexp( z/freq[p] + chopL/y, freq[p]) )	
	//MultiThread phi_fvt_diff = pi*(b^2-centerPos^2)/(b^2-beampos[i]^2)*(vpi/y)^2* ( foftsqrd( z/freq[p] , freq[p]) - deltaphi*foftsqrd( z/freq[p] + chopL/y, freq[p]) )	
	//	MultiThread phi_fvt_diff = pi*( foftsqrd( z/x , x) - foftsqrd( z/x + chopL/y, x) )	
	//	MultiThread phi_fvt_diff = pi*( foftsqrdexp( z/x , x) - foftsqrdexp( z/x + chopL/y, x) )	
	//	Print "timer phi_fvt_diff: ", (StopMSTimer(-2)-ttime)*10^-6	
	//MultiThread phi_fvt_re =  cos(	phi_fvt_diff	)
	//MultiThread phi_fvt_im = sin(	phi_fvt_diff	)
	MatrixOP/O/S phi_fvt_re = cos(	phi_fvt_diff	)
	MatrixOP/O/S phi_fvt_im = sin(	phi_fvt_diff	)
	Integrate/DIM=2 phi_fvt_re, phi_fvt_im
	
	//Print "timer cos&sin&Int phi_fvt_im: ", (StopMSTimer(-2)-ttime)*10^-6		

	phi_fv_re = phi_fvt_re[p][q][tPnts]*velprobs[q]
	phi_fv_im = phi_fvt_im[p][q][tPnts]*velprobs[q]
	Integrate/DIM=1 phi_fv_re, phi_fv_im
	
	phi_f_re = phi_fv_re[p][vPnts]
	phi_f_im = phi_fv_im[p][vPnts]
	
	//phi_f_re+=phi_f_re_temp
	//phi_f_im+=phi_f_im_temp
	//endfor

	//phi_f_re/=numpnts(beampos)
	//phi_f_im/=numpnts(beampos)
	
	phi_f = atan2(phi_f_im, phi_f_re)
	phi_f += pi+1
	phi_f = mod(phi_f, 2*pi)	
	phi_f-=pi+1
	
	con_f = sqrt(phi_f_im^2 + phi_f_re^2)
	
	yw = con_f*pw[3]
	
	//Print "timer total v2: ", (StopMSTimer(-2)-ttime)*10^-6
End




Function CleanupChopperFitWaves()
	KillWaves/Z phi_fvt_diff, phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_f_re, phi_f_im, phi_f_re_temp, phi_f_im_temp, foftsqrdexpChop12, foftsqrdexpChop1, foftsqrdexpChop2, phi_fvt_diff_mol, sagPhaseWave, phiSagGravChop
	KillWaves/z voltageon, phi_fvt_diffM, phi_fvt_diff2
End



Function SagnacAndGravityPhaseChop()
	Wave velwave, velprobs
	
	Variable verbose=0	// set = 1 if you want phase shift and contrast printed to history
	
	Duplicate/O velWave phiGrav, phiSag, phiSagGrav, phiSagGravWeighted, phiSagGravImag, phiSagGravReal
	
	Variable GratingTilt = 0e-3	// in radians         // 0.5*pi/180
	Variable L1g2g = 0.94		// meters
	Variable g = 9.8				// meters/s^2
	Variable gFactor = g*sin(GratingTilt)*L1g2g^2/a*2*pi
	phiGrav = gFactor * velwave^-2

	Variable Latitude = 32.2	//degrees
	Variable OmegaEarth = 2*pi/(24*3600)*cos((90-Latitude)*pi/180)    //; print OmegaEarth	// earth rotation rate for a given latitude
	Variable SagFactor = OmegaEarth*L1g2g^2*4*pi/1e-7
	//print sagfactor/3000
	phiSag = SagFactor / velwave
	
	phiSagGrav = phiGrav+phiSag
	phiSagGravWeighted = velprobs*phiSagGrav
	Variable phiSagGravAvg = area(phiSagGravWeighted); 
	
	// Calculate how much we'll need to unwrap the phase when we take the atan2 of the integrated e^iphi
	Variable/G k = -1
	do
		k+=1
	while(phiSagGravAvg-k*2*pi > pi)
	//print k
	
	phiSagGravImag = velProbs * sin(phiSagGrav)
	phiSagGravReal = velProbs * cos(phiSagGrav)
	Variable phiSagGravImagInt = area(phiSagGravImag)
	Variable phiSagGravRealInt = area(phiSagGravReal)
	
	Variable/G AvgSagGravPhi = atan2(phiSagGravImagInt, phiSagGravRealInt) + 2*pi*k
	Variable/G AvgSagGravCon = sqrt(phiSagGravImagInt^2+phiSagGravRealInt^2)
	
	if(verbose)
		printf "Improper velocity averaged Sagnac and gravity phase shift: %2.7f\r", phiSagGravAvg
		printf "Velocity averaged Sagnac and gravity phase shift: %2.7f\r", avgSagGravphi
		printf "Velocity averaged Sagnac and gravity contrast: %1.3f\r", avgSagGravCon
	endif
	
	Duplicate/O foftsqrdexpChop1 phiSagGravChop
	phiSagGravChop = phiSagGrav[q]
End




Function CalculateChopperBestFit(w_coef, frequencyWave, series_name)
	Wave w_coef, frequencyWave
	String series_name
	
//	Variable bestFitPnts = numpnts(frequencyWave)//200
	Variable bestFitPnts = 200
	
	Make/O/N=(bestFitPnts) contrastFit, freqFit, phaseFit
	
	Variable minf = wavemin(frequencyWave)
	Variable maxf = wavemax(frequencyWave)
	SetScale/I x minf, maxf, "Hz" contrastfit, freqfit, phasefit
	freqFit = x
	
//	NVAR fpnts
//	fpnts = bestFitPnts	
	MakeChopperFitHelperWavesV3(freqFit)
//	MakeChopperFitHelperWavesTrans(freqFit)
	
	PhaseChoppers2IFMfitv3(w_coef, contrastFit, freqFit)
	
	Wave phi_f
	phaseFit = phi_f
	phaseFit+=pi
	
	String name = "contrastFit_"+series_name
	Duplicate/O contrastFit $name
	
	name = "phaseFit_"+series_name
	Duplicate/O phaseFit $name
	
	CleanupChopperFitWaves()
End






Function AppendBestFitToGraph2(w_coef, freqFit)
	Wave w_coef, freqFit
	
	Duplicate/O freqFit contrastFit2
	
	NVAR fpnts
	fpnts = numpnts(freqFit)	
	MakeChopperFitHelperWaves(freqFit)
	
	PhaseChoppers2IFMfitv2(w_coef, contrastFit2, freqFit)
End






Function AppendResults(params)
	Wave/T params
	Wave w_coef, w_sigma
	
	String text = "", textTotal = "\\Zr100"
	Variable n = 0
	do
		sprintf text, "  %s = %g ± %g  ", params[n], w_coef[n], w_sigma[n]
		print text
		textTotal += text
		if(n<numpnts(w_coef)-1)
			textTotal+="\r"
		endif
		n+=1
	while(n<numpnts(w_coef))
	//print textTotal
	TextBox/C/N=ResultsTBox /A=RT textTotal
End






Function/S ProcessHoldAndConstraintsWaves()
	Variable n
	Wave HoldWave = HoldWave
	Wave/T Constraints = Constraints, Constraints1= Constraints1, Constraints2 = Constraints2
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
	
	Concatenate/O/NP/T {Constraints1Dup, Constraints2Dup}, $nameofwave(Constraints)
	Variable conLength = numpnts(Constraints)
	For(n=0;n<conLength;n+=1)
		If (stringmatch(Constraints[n],"")==1 && numpnts(constraints)>0)
			DeletePoints n,1, Constraints
			n -= 1
		Endif
	EndFor
	
	KillWaves Constraints1Dup, Constraints2Dup
	
	Return hold
end







Function PhaseChoppers2IFMfit(pw, yw, xw) : FitFunc
	Wave pw, yw, xw
	
	// pw0 = flow velocity
	// pw1 = vratio
	// pw2 = vpi/vflow
	// pw3 = normalization constant
	
	Variable ttime= stopMSTimer(-2)
	
	Variable vPnts = 100
	Variable nSigmas = 5
	Variable FlowV =pw[0]
	Variable vratio = pw[1]
	Variable sigV = FlowV/vratio
	Variable VelNorm
	
	Variable vpi = flowv*pw[2]
	
	Variable tpnts = 100
	Variable dutycycle = 1/2
	
	Variable posPnts = 1
	Variable centerPos = .8e-3
	Variable beamWidth = 50e-6
	Variable minPos = centerPos-beamWidth/2
	Variable maxPos = centerPos+beamWidth/2

	Variable a = 1e-3	// ground plane - wire dist
	Variable radius = 0.785e-3
	Variable b = a*sqrt(1+2*radius/a)
	
	//Variable chopL = 1.4	//meter
	Variable chopL = 1.269
	
	Variable fPnts = 500
	Variable minChopFreq = 1
	Variable maxChopFreq = 20e3
	//Variable maxOmega = 2*pi*maxChopFreq
	
	Variable/G minVel = FlowV*(1-nSigmas/15)//; minVel = 0
	Variable/G maxVel = FlowV*(1+nSigmas/15)//; maxVel = 5000
	Variable/G velStep = (maxVel-minVel)/vPnts
	
	Make/D/O/N=(vPnts) velwave = minVel+p*velStep
	Duplicate/O velwave velprobs 
	SetScale/I x minVel,maxVel,"m/s", velprobs
	Velprobs = exp(-(x-FlowV)^2/(2*SigV^2))*x^3
	VelNorm = sqrt(2*pi)*sigv*FlowV*(FlowV^2+3*sigv^2)
	//Variable velsum = sum(velprobs,-inf,inf)
	Velprobs/=velNorm
	
	Make/D/O/N=(fPnts, vPnts, tPnts) phi_fvt_re, phi_fvt_im, tTest, phi_fvt_diff
	Make/D/O/N=(fPnts, vPnts) phi_fv_re, phi_fv_im
	Make/D/O/N=(fPnts) phi_f_re, phi_f_im, phi_f, contrast_f
	SetScale/I x minChopFreq, maxChopFreq, "Hz", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, phi_f_re, phi_f_im, phi_f, contrast_f, ttest, phi_fvt_diff
	SetScale/I y minVel, maxVel, "m/s", phi_fvt_re, phi_fvt_im, phi_fv_re, phi_fv_im, ttest, phi_fvt_diff
	SetScale/P z 0, 1/tpnts, "t*f", phi_fvt_re, phi_fvt_im, ttest, phi_fvt_diff
	
	Make/D/O/N=(posPnts) beamPos
	SetScale/i x minPos, maxPos, "m", beamPos
	beampos=x
	beampos=centerpos
	Duplicate/O phi_f_re, phi_f_re_temp, phi_f_im_temp
	phi_f_re=0;phi_f_im=0;phi_f_re_temp=0; phi_f_im_temp=0
	
	//Print "timer assign vars, make waves: ", (StopMSTimer(-2)-ttime)*10^-6	

	Variable i = 0
	for(i=0; i<numpnts(beampos); i+=1)
		// foftsqrd( time(s), frequency(Hz) )
		//	MultiThread ttest = foftsqrd( z/x , x)
		//	MultiThread phi_fvt_diff = pi*(vpi/y)^2* ( foftsqrd( z/x , x) - foftsqrd( z/x + chopL/y, x) )	
		//print "i= ", i, "; beampos[i] = ", beampos[i]
		MultiThread phi_fvt_diff = pi*(b^2-centerPos^2)/(b^2-beampos[i]^2)*(vpi/y)^2* ( foftsqrdexp( z/x , x) - foftsqrdexp( z/x + chopL/y, x) )	
		//	MultiThread phi_fvt_diff = pi*( foftsqrd( z/x , x) - foftsqrd( z/x + chopL/y, x) )	
		//	MultiThread phi_fvt_diff = pi*( foftsqrdexp( z/x , x) - foftsqrdexp( z/x + chopL/y, x) )	
		//	Print "timer phi_fvt_diff: ", (StopMSTimer(-2)-ttime)*10^-6	
		//MultiThread phi_fvt_re =  cos(	phi_fvt_diff	)
		//MultiThread phi_fvt_im = sin(	phi_fvt_diff	)
		MatrixOP/O/S phi_fvt_re = cos(	phi_fvt_diff	)
		MatrixOP/O/S phi_fvt_im = sin(	phi_fvt_diff	)
		Integrate/DIM=2 phi_fvt_re, phi_fvt_im
	
		//Print "timer cos&sin&Int phi_fvt_im: ", (StopMSTimer(-2)-ttime)*10^-6		

		phi_fv_re = phi_fvt_re[p][q][tPnts]*velprobs[q]
		phi_fv_im = phi_fvt_im[p][q][tPnts]*velprobs[q]
		Integrate/DIM=1 phi_fv_re, phi_fv_im
	
		phi_f_re_temp = phi_fv_re[p][vPnts]
		phi_f_im_temp = phi_fv_im[p][vPnts]
	
		phi_f_re+=phi_f_re_temp
		phi_f_im+=phi_f_im_temp
	endfor

	phi_f_re/=numpnts(beampos)
	phi_f_im/=numpnts(beampos)
	
	phi_f = atan2(phi_f_im, phi_f_re)
	contrast_f = sqrt(phi_f_im^2 + phi_f_re^2)
	
	yw = contrast_f(xw[p])*pw[4]
	
	Print "timer total v1: ", (StopMSTimer(-2)-ttime)*10^-6
End




ThreadSafe Function foftsqrd(t, f)
	Variable t, f
	
	Variable dutycycle = 0.5
	
	Variable n = floor(t*f)
	
	if(t-n/f < dutycycle/f)
		return 1
	else
		return 0
	endif
End






ThreadSafe Function foftsqrdexp(t, f)
	Variable t, f
	// t is in units of seconds
	// f is in units of Hz
	
	Variable dutycycle = 0.5
	Variable frc = 500e3
	Variable wrc = 2*pi*frc	// 2*pi*frc
	
	Variable n = floor(t*f)
	//print n
	
	Variable tshifted = t-n/f
	
	if(tshifted < dutycycle/f)
		return (1-exp(-tshifted*wrc))^2
	else
		return (exp( -( (tshifted - (dutycycle/f)) * wrc)) )^2
		//return (1-exp(-(dutycycle/f)*wrc/f))*exp(-((tshifted - (dutycycle/f))*wrc/f))
	endif
End




 Function foftsqrdtrans(t, chopfreq, velocity, intnorm)
	Variable t, chopfreq, velocity, intnorm
	// t is in units of seconds
	// f is in units of Hz
	
	Variable dutycycle = 0.5
	Variable frc = 500e3
	Variable wrc = 2*pi*frc	// 2*pi*frc
	
	Variable tshifted = mod(t,1/chopfreq)		// map t back into the first chopper period
	
	Wave voltageOn,esqrdofZ,esqrdofZintGated
	//NVAR lc,intnorm,zPnts
	Variable toff = 1/(2*chopfreq)
	Variable zoff = velocity * toff
	voltageOn = mod((x+lc/2+t*chopfreq*2*zoff), 2*zoff) < zoff ? 1 : 0
	esqrdofZintGated=esqrdofZ*voltageOn
	Integrate esqrdofZintGated
	esqrdofZintGated/=intnorm
	return esqrdofZintGated[zPnts-1]
End



Function Testfoftsqrdtrans()
	variable intnorm = CalcEsqrdOfZChop()
	
	print foftsqrdtrans(.01, 1e5, 1000,intnorm)
End


Function TestHelperTrans()
	make/o/n=50 freqtest = (p+1)*200
	
	MakeChopperFitHelperWavesTrans(freqtest)
End

Function TestHelperV3()
	make/o/n=10 freqtest = (p+1)*500
	
	MakeChopperFitHelperWavesV3(freqtest)
End


Constant zPnts = 1000
	
Constant lc = 0.02	//distance to integrate esqrd over


Function CalcEsqrdOfZChop()
	Variable phase
	
	Variable a = 1e-3		//wire to ground plane
	Variable D = 1.57e-3	//wire diameter
	Variable b = a*sqrt(1+D/a)
	
	Variable x0 = 0.75e-3	//approximate beam location
	Variable s = 0.01e-3	//approximate beam separation
	
	Make/O/D/N=(zPnts) esqrdofZ
	SetScale/i x -lc/2, lc/2, esqrdofZ
	
	esqrdofZ = ( ((x0+b)^2 + x^2)* ((x0-b)^2+x^2) )^-1 - ( ((x0-s+b)^2 + x^2)* ((x0-s-b)^2+x^2) )^-1
	
	Duplicate/O esqrdOfZ esqrdofZint
	
	Integrate esqrdofZint
	Variable/G intnorm = esqrdofZint[zPnts]
	esqrdofZint/=intnorm
	
	Duplicate/O esqrdofZint voltageOn
	
	Variable velocity = 1000	
	Variable chopfreq = 30e3	// chopper frequency
	Variable toff = 1/(2*chopfreq)
	Variable zoff = velocity * toff
	//Variable phase = pi/2
	voltageOn = mod((x+lc/2+phase/(2*pi)*2*zoff), 2*zoff) < zoff ? 1 : 0
	
	Duplicate/O esqrdofZ esqrdofZintGated
	
	esqrdofZintGated*=voltageOn
	Integrate esqrdofZintGated
	esqrdofZintGated/=intnorm
	//print esqrdofZintGated[zPnts]
	
	return intnorm
	
	//return esqrdofZintGated[zPnts-1]
End





Function CalcEsqrdOfZ(phase)
	Variable phase
	
	Variable a = 1e-3		//wire to ground plane
	Variable D = 1.57e-3	//wire diameter
	Variable b = a*sqrt(1+D/a)
	
	Variable x0 = 0.75e-3	//approximate beam location
	Variable s = 0.01e-3	//approximate beam separation
	
	Variable zPnts = 1000
	
	Variable lc = 0.02	//distance to integrate esqrd over
	
	Make/O/D/N=(zPnts) esqrdofZ
	SetScale/i x -lc/2, lc/2, esqrdofZ
	
	esqrdofZ = ( ((x0+b)^2 + x^2)* ((x0-b)^2+x^2) )^-1 - ( ((x0-s+b)^2 + x^2)* ((x0-s-b)^2+x^2) )^-1
	
	Duplicate/O esqrdOfZ esqrdofZint
	
	Integrate esqrdofZint
	
	Variable intnorm = esqrdofZint[zPnts]
	esqrdofZint/=intnorm
	
	Duplicate/O esqrdofZint voltageOn
	
	Variable velocity = 1000	
	Variable chopfreq = 30e3	// chopper frequency
	Variable toff = 1/(2*chopfreq)
	Variable zoff = velocity * toff
	//Variable phase = pi/2
	voltageOn = mod((x+lc/2+phase/(2*pi)*2*zoff), 2*zoff) < zoff ? 1 : 0
	
	Duplicate/O esqrdofZ esqrdofZintGated
	
	esqrdofZintGated*=voltageOn
	Integrate esqrdofZintGated
	esqrdofZintGated/=intnorm
	print esqrdofZintGated[zPnts]
	return esqrdofZintGated[zPnts]
End



Function gogater()
	Variable i
	Variable imax=500
	
	Make/O/N=(imax) phaseout, chopphase, sqr
	phaseout=nan
	chopphase=p/imax*2*pi
	//setscale/i x 0,2*pi
	sqr = mod(p/imax*2*pi+.6, 2*pi) < pi ? 1 : 0
	
	for(i=0; i<imax; i+=1)
		phaseout[i] =CalcEsqrdofZ(chopphase[i])
		DoUpdate
		//Sleep/T 1
	endfor
End






ThreadSafe Function phiChopt2(t, f)
	Variable t, f
	// t is in units of seconds
	// f is in units of Hz
	
	Variable dutycycle = 0.5
	Variable frc = 60e3
	Variable wrc = 2*pi*frc	// 2*pi*frc
	
	Variable n = floor(t*f)
	//print n
	
	Variable tshifted = t-n/f
	
	return pi*(tshifted*f)^2
End





