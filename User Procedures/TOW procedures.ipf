#pragma rtGlobals=1		// Use modern global access method.

#include <Remove Points>
#include ":RemoveNaNs"
#include ":Fringe Fit"
#include ":PhysicalConstants"
#include ":LabConstants"
#include ":WindowNamer"

/////////////////////////////////////////////////////////////////////////////////////
//
//
//	This procedure file contains all of the functions specific to TOW/MZW data processing.
//	It was originally written by Will Holmgren in 2012.
//	
//	extractTOW(series_name) is the main function of this procedure. It takes the results of Fringe Fit
//		and processes them into phase and contrast vs wavelength data. Very long, lots of options.
//
//	Important helper functions for extractTOW() are:
//		avgPhases(seriesName)  
//		Binner(phase_on_Active, phaseBinned_Active, phaseBinnedStdErr_Active, phaseBinnedStdDev_Active, wavelengths, wavelengthsBinned, wavelengthBinWidth, filterType, phaseSetStdDevFilter, filesPerWavelengthThreshold, trimmedMeanPercent, refStdDev, useRefStdDev)
//	
//	getWavelengthData() parses wavelength data files.
//
//	goFitTOWsaf(series_name,all) and goFitTOWsafR(series_name,all) fits the phase vs wavelength data using the full atomic theory. 
//		The parameter 'all' controls which type of processed phase data to use (spline, poly, binned). see the code for details.
//		(Saf stands for Safronova, just because I originally took the equation from one of her papers.)
//
//	FitTOWsaf(pw,yw,xw) and FitTOWsafR(pw,yw,xw) generate the data needed for the the goFitTOW functions.
//
//
//	Other functions of note:
//		getEtalonData() parses etalon data files.
//		FilterTOWfit(initialize, StdDevFilter) filters the data according to its distance from the best fit curve.
//		CalcR(TOW) calculate R given a TOW.
//
/////////////////////////////////////////////////////////////////////////////////////



// Take the results of fringe fit for a specified series and process them into phase vs wavelength.
// The variable wattsPerVolt sets the scale of the laser power calibration.
// You can optionally make a wave named "KillFiles_"+series_name to manually kill certain files.
// There are a lot of lines of code, but a large part is just to deal with the wave names and displaying plots
// There's also extra code to compare the results when using a spline or a polynomial to fit the reference phase.
// Eventually, extractTOW calls the function Binner(). Binner() does the hard work of turning the delta phase
// measurements into binned and filtered measurements of phase vs wavelength. 
// Binner has two parameters that you may want to change: 	
//		Variable minCutoffWavelength = 768.970-0.100
//		Variable maxCutoffWavelength = 768.970+0.100

Function extractTOW(series_name)
	String series_name
	
	Variable wattsPerVolt =  .05 //1.4
	//Variable wattsPerVolt =  1.4/10 
	
	// decide to use 1k2k or 2k2k fringe data
	Variable OnekTwok =1
	String onektwokStr = ""
	if(OnekTwok==1)
		onektwokStr = "1k2k"
	elseif(OnekTwok == 2)
		onektwokStr = "2k2k"
	endif
	
	// set up wave names
	String inputphase = "phase" + onektwokStr + "_" + series_name
	String inputphase_error = "phase_error" + onektwokStr + "_" + series_name

	String phase_on_name = "phase" + onektwokStr + "_on_" + series_name
	String phase_error_on_name = "phase_error" + onektwokStr + "_on_" + series_name
	
	String phase_ref_name = "phase" + onektwokStr + "_ref_" + series_name
	String phase_error_ref_name = "phase_error" + onektwokStr + "_ref_" + series_name
	
	String inputcontrast = "contrast" + onektwokStr + "_" + series_name
	String inputcontrast_error = "contrast_error" + onektwokStr + "_" + series_name

	String contrast_on_name = "contrast" + onektwokStr + "_on_" + series_name
	String contrast_error_on_name = "contrast_error" + onektwokStr + "_on_" + series_name
	
	String contrast_ref_name = "contrast" + onektwokStr + "_ref_" + series_name
	String contrast_error_ref_name = "contrast_error" + onektwokStr + "_ref_" + series_name
	
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

	Wave contrast_in = $inputcontrast
	Wave contrast_error_in = $inputcontrast_error	
	Duplicate/O contrast_in $contrast_on_name, $contrast_ref_name
	Duplicate/O contrast_error_in $contrast_error_on_name, $contrast_error_ref_name
	
	Wave contrast_on = $contrast_on_name
	Wave contrast_error_on = $contrast_error_on_name
	Wave contrast_ref = $contrast_ref_name
	Wave contrast_error_ref = $contrast_error_ref_name
	
		
	// provide switches for the number of on/off files and whether or not the data starts with the reference phase.
	Variable NumOnOff = 1
	
	Variable refPhaseStartsSeries = 1	// 1 for yes, 0 for no
	Variable phaseOnBool, phaseOffBool
	if(refPhaseStartsSeries)
		phaseOnBool = 1
		phaseOffBool = 0
	elseif(refPhaseStartsSeries==0)
		phaseOnBool = 0
		phaseOffBool = 1
	else
		Abort "Which files are the reference files?"
	endif
	

	
	// assign the on and reference conditions. NaN the wrong points.
	
	phase_on *= (mod(Floor(p/NumOnOff),2)==phaseOnBool)	// keep every other set of (NumOnOff) points
	//phase_on *= (mod(p,NumOnOff)!=0)	// delete first points from each set
	phase_ref *= (mod(Floor(p/NumOnOff),2)==phaseOffBool)
	
	contrast_on *= (mod(Floor(p/NumOnOff),2)==phaseOnBool)	// keep every other set of (NumOnOff) points
	contrast_ref *= (mod(Floor(p/NumOnOff),2)==phaseOffBool)	
	
	phase_on = phase_on==0 ? nan : phase_on
	phase_ref = phase_ref==0 ? nan : phase_ref
	
	phase_error_on = numtype(phase_on)==0 ? phase_error_on : nan
	phase_error_ref = numtype(phase_ref)==0 ? phase_error_ref : nan

	contrast_on = contrast_on==0 ? nan : contrast_on
	contrast_ref = contrast_ref==0 ? nan : contrast_ref
	
	contrast_error_on = numtype(contrast_on)==0 ? contrast_error_on : nan
	contrast_error_ref = numtype(contrast_ref)==0 ? contrast_error_ref : nan
	
	
	
	// Manual elimination of files. 
	
	String killFilesName = "KillFiles_"+series_name
	Wave/Z KillFiles = $killFilesName
	
	If(waveexists(killfiles))
		Variable i
		Variable killpnt
		For(i=0; i<numpnts(killFiles); i+=1)
			killpnt = killFiles[i]-1
			phase_on[killpnt]=nan
			phase_ref[killpnt]=nan
			contrast_on[killpnt]=nan
			contrast_ref[killpnt]=nan
		EndFor
	endif


	//Displays a figure with phase and reference phase. Fits Ref phase to a polynomial and subtracts the fit
	//from both the phase and ref phase. Therefore, the new ref phase is the fit residuals.
	Display/K=1/W=(50,50,750,500) phase_on vs fileNumbers
	AppendToGraph phase_ref vs fileNumbers
	ModifyGraph zero(left)=1
	WindowNamer("Phase on ref " + oneKTwoKStr + " " + series_name)
	ModifyGraph nticks(bottom)=10,minor(bottom)=1
	ModifyGraph rgb($phase_ref_name)=(0,0,0)
	ErrorBars $phase_on_name Y,wave=($phase_error_on_name,$phase_error_on_name)
	ErrorBars $phase_ref_name Y,wave=($phase_error_ref_name,$phase_error_ref_name)

	Variable refPolyOrder = 6
	Make/O w_coef
	CurveFit/NTHR=0 poly refPolyOrder,  phase_ref /W=phase_error_ref /I=1 /R /A=0
	String reswavename = "res_"+phase_ref_name
	Wave reswave = $reswavename
	//ErrorBars $reswavename Y,wave=($phase_error_ref_name,$phase_error_ref_name)
	reswave/=phase_ref; reswave*=phase_ref	// NaNs points in the residual wave if a point is added to killFiles
	
	ModifyGraph mode=3
	ModifyGraph marker=8
	Label left "phase shift (rad)"
	Label bottom "file number"
	
	// subtract out the reference phase
	Variable k
	for(k=0; k<refPolyOrder; k+=1)
		phase_on -= w_coef[k]*p^k
		phase_ref -= w_coef[k]*p^k
	endfor

	
	
	//Displays a figure with phase and reference phase. Fits Ref phase with a smoothing spline and subtracts the fit
	//from both the phase and ref phase. Therefore, the new ref phase is the fit residuals.
	Variable displayRefPhaseFigures = 1
	if(displayRefPhaseFigures)
		String phase_on_spline_name = "phase" + onektwokStr + "_on_spline_" + series_name
		String phase_ref_spline_name  = "phase" + onektwokStr + "_ref_spline_" + series_name
		String phase_error_ref_spline_name = "phase_error" + onektwokStr + "_ref_spline_" + series_name
		Duplicate/O phase_on $phase_on_spline_name
		Duplicate/o phase_ref $phase_ref_spline_name
		Duplicate/O phase_error_ref $phase_error_ref_spline_name
		Wave phaseOnSpline = $phase_on_spline_name
		Wave phaseRefSpline = $phase_ref_spline_name
		Wave phaseErrorRefSpline = $phase_error_ref_spline_name
		
		String fileNumbers_ref_name = "file_index_ref_" + series_name
		Duplicate/o fileNumbers $fileNumbers_ref_name
		Wave fileNumbersRef = $fileNumbers_ref_name
		
		RemoveNaNsXYZ(phaseRefSpline, phaseErrorRefSpline, fileNumbersRef)
		
		String phaseRefSplineSSName = "phase_ref_spline_SS_"+series_name
		Variable smoothing = 2
		Interpolate2/T=3/N=2000/F=(smoothing)/SWAV=phaseErrorRefSpline/Y=$phaseRefSplineSSName fileNumbersRef, phaseRefSpline
		Wave phaseRefSplineSS = $phaseRefSplineSSName
		
		phaseOnSpline -= phaseRefSplineSS(fileNumbers[p])
		phaseRefSpline -= phaseRefSplineSS(fileNumbersRef[p])
		
		Display/K=1/W=(650,50,1350,500) phaseRefSpline vs fileNumbersRef
		AppendToGraph phaseOnSpline vs fileNumbers
		ModifyGraph mode=3
		ModifyGraph marker=8
		ModifyGraph zero(left)=1
		Label left "phase shift (rad)"
		Label bottom "file number"
		WindowNamer("Phase on ref spline " + series_name)
		ModifyGraph nticks(bottom)=10,minor(bottom)=1
		ModifyGraph rgb($phase_ref_spline_name)=(0,0,0)
		ErrorBars $phase_on_spline_name Y,wave=($phase_error_on_name,$phase_error_on_name)
		ErrorBars $phase_ref_spline_name Y,wave=($phase_error_ref_spline_name,$phase_error_ref_spline_name)
	endif
	
	// display reference contrast if desired.
	Variable displayRefContrastFigure = 1
	if(displayRefContrastFigure)
		Display/K=1/W=(150,250,950,700) contrast_on vs fileNumbers
		AppendToGraph contrast_ref vs fileNumbers
		ModifyGraph zero(left)=1
		WindowNamer("contrast on ref " + oneKTwoKStr + " " + series_name)
		ModifyGraph nticks(bottom)=10,minor(bottom)=1
		ModifyGraph rgb($contrast_ref_name)=(0,0,0)
		ErrorBars $contrast_on_name Y,wave=($contrast_error_on_name,$contrast_error_on_name)
		ErrorBars $contrast_ref_name Y,wave=($contrast_error_ref_name,$contrast_error_ref_name)
		ModifyGraph mode=3
		ModifyGraph marker=8
		Label left "contrast"
		Label bottom "file number"
	endif
	
	
	
	/////	The next section deals with the power normalization	//////

	//	call up the avg power wave for this series
	String power_name = "avg_power_" + series_name			// avg photodiode voltage in each file (comes from fringe fit procedure)
	Wave avg_power_loc = $power_name
	String powerStdDevName = "avg_powerStdDev_"+series_name
	Wave powerStdDev = $powerStdDevName
		
	String powerW_name = "powerW_" + series_name			// power in W series name
	String powerStdDevW_name = "powerStdDevW_" + series_name			// power in W series name
	String ref_powerV_name = "ref_powerV_"+series_name	// reference power voltage
	String powerVeto_name = "powerVeto_"+series_name
	Duplicate/O avg_power_loc $powerW_name, $ref_powerV_name
	Duplicate/O powerStdDev $powerStdDevW_name, $powerVeto_name
	Wave powerW = $powerW_name
	Wave ref_powerV = $ref_powerV_name
	Wave powerStdDevW = $powerStdDevW_name
	Wave powerVeto = $powerVeto_name
	
	//Subtract off the power measurement without the light on
	ref_powerV = numtype(phase_ref)==0 ? ref_powerV : nan	// NaN the non reference condition points
	Wavestats/Q ref_powerV
	Variable avg_refV = V_avg
	powerW-= avg_refV
	
	// Convert voltage to Watts
	powerW *= wattsPerVolt		
	powerStdDevW *= wattsPerVolt								

	Variable powerVetoCut = 0.010	//.005
	powerVeto = (powerStdDevW/powerW > powerVetoCut && mod(p,2)!=0) ? powerW : NaN
	powerW = powerVeto==powerW ? NaN : powerW

	Display/K=1/W=(150,250,850,700) powerW vs fileNumbers
	ErrorBars $powerW_name Y,wave=($powerStdDevW_name,$powerStdDevW_name)
	Label left "Power (W)"
	Label bottom "file number"
	AppendToGraph powerVeto vs fileNumbers
	ErrorBars $powerVeto_name Y,wave=($powerStdDevW_name,$powerStdDevW_name)
	ModifyGraph mode=3
	ModifyGraph marker=8
	ModifyGraph marker($powerVeto_name)=1,rgb($powerVeto_name)=(0,0,0)
	WindowNamer("Power " + series_name)
	
	// make new waves to hold the normalized data
	String phase_on_Pnorm_name = "phase_on_Pnorm_" + series_name
	String phase_on_SS_Pnorm_name = "phase_on_SS_Pnorm_" + series_name
	String phase_error_on_Pnorm_name = "phase_error_on_Pnorm_" + series_name
	Duplicate/O phase_on $phase_on_Pnorm_name, $phase_error_on_Pnorm_name
	Duplicate/O phaseOnSpline $phase_on_SS_Pnorm_name
	Wave phase_on_Pnorm = $phase_on_Pnorm_name
	Wave phase_on_SS_Pnorm = $phase_on_SS_Pnorm_name
	Wave phase_error_on_Pnorm = $phase_error_on_Pnorm_name
	
	// do the normalization if desired.
	Variable normalizePowerBool = 1
	if(normalizePowerBool)
		phase_on_Pnorm /= powerW	
		phase_on_SS_Pnorm /= powerW
		phase_error_on_Pnorm /=powerW
	else
		//do nothing to the Pnorm waves
	endif
	

	
	
	/// separate the good data (taken at a fixed laser frequency) from the bad data (taken while tuning the laser).
	//   These functions also plot the processed data. Lots of repeated code. It's not ideal, but it was fast at the time.
	//   3 options:
	//   0: Call the AvgPhases function to do smart binning/averaging. Recommended.
	//	1: extractMessy2() relies on 2 waves, extract_start_* and extract_end_* to pull out the the good data
	//	2: this option assumes that it takes x files out of a n*x long set to tune the laser frequency. This ends up throwing away a lot of good points.
	Variable extractSetting=0
	if(extractSetting==0)
		AvgPhases(series_Name)
	elseif(extractSetting==1)
		extractMessy2(series_name)
	elseif(extractSetting==2)
		Variable numsets = floor(numpnts(phase_on)/numOnOff/2)
		
		String phase_on_avg_name = "phase" +onektwokStr+ "_on_" + "avg_" + series_name
		String phase_error_on_avg_name = "phase_error"+onektwokStr + "_on_" + "avg_" + series_name
		String phase_error_on_avg_sdev_name = "phase_error"+onektwokStr + "_on_" + "avg_sdev_" + series_name
		
		Make/O/D/N=(numsets) $phase_on_avg_name, $phase_error_on_avg_name, $phase_error_on_avg_sdev_name
		
		Wave phase_on_avg = $phase_on_avg_name
		Wave phase_error_on_avg = $phase_error_on_avg_name
		Wave phase_error_on_avg_sdev = $phase_error_on_avg_sdev_name
		
		Duplicate/O phase_on_avg phase_on_noPnorm_avg, phase_error_on_noPnorm_avg
		
		Variable n
		For(n=0; n<numsets; n+=1)
			Wavestats/Q/R=[n*numOnOff*2,(n+1)*numOnOff*2-1] phase_on
			phase_on_avg[n] = V_avg
			phase_error_on_avg[n] = v_sem
			phase_error_on_avg_sdev[n] = v_sdev
			
			Wavestats/Q/R=[n*numOnOff*2,(n+1)*numOnOff*2-1] phase_on_noPnorm
			phase_on_noPnorm_avg[n] = V_avg
			phase_error_on_noPnorm_avg[n] = v_sem
		EndFor

		String wavelength_name = "wavelength_"+series_name
		Wave/Z wavelength = $wavelength_name
	
		if(waveexists(wavelength))
			String wavelength_pad_name = "wavelength_pad_"+series_name
			Make/O/N=(numpnts(phase_on)) $wavelength_pad_name
			Wave wavelength_pad = $wavelength_pad_name
			
			wavelength_pad = wavelength[floor(p/(NumOnOff*2))]
			
			Display/K=1/W=(550,50,1250,500) phase_on vs wavelength_pad
			ModifyGraph mode=3
			ModifyGraph marker=8
			if(normalizePowerBool)
				Label left "phase shift (rad/W)"
			else
				Label left "phase shift (rad)"
			endif
			Label bottom "wavelength (nm)"
			ModifyGraph zero(left)=1
			ErrorBars $phase_on_name Y,wave=($phase_error_on_name,$phase_error_on_name)
			ModifyGraph nticks(bottom)=10,minor(bottom)=1
			WindowNamer("Phase vs wavelength " + series_name)
			
			
			Display/K=1/W=(550,450,1250,900) phase_on_avg vs wavelength
			ModifyGraph mode=3
			ModifyGraph marker=8
			if(normalizePowerBool)
				Label left "phase shift (rad/W)"
			else
				Label left "phase shift (rad)"
			endif
			Label bottom "wavelength (nm)"
			ModifyGraph zero(left)=1
			ErrorBars $phase_on_avg_name Y,wave=($phase_error_on_avg_name,$phase_error_on_avg_name)
			ModifyGraph nticks(bottom)=10,minor(bottom)=1
			WindowNamer("Phase avg vs wavelength " + series_name)
		endif
	endif
	
End




// uses waves named "extract_start_"+series_name and "extract_end_"+series_name to separate out the good data from the data taken while tuning the laser.
// It also makes plots of the processed data.
Function extractMessy2(series_name)
	String series_name

	String extract_start_name = "extract_start_"+series_name
	String extract_end_name = "extract_end_"+series_name
	Wave extract_start = $extract_start_name
	Wave extract_end = $extract_end_name
	Variable numsets=numpnts(extract_start)

	String wavelengths_name = "wavelengths_"+series_name
	Wave wavelengths = $wavelengths_name	
	
	String wavelength_name = "wavelength_"+series_name
	Make/O/N=(numsets) $wavelength_name
	Wave wavelength = $wavelength_name	
	
	String phase_on_avg_name = "phase1k2k_on_" + "avg_" + series_name
	String phase_error_on_avg_name = "phase_error1k2k_on_" + "avg_" + series_name
	String phase_error_on_avg_sdev_name = "phase_error1k2k_on_" + "avg_sdev_" + series_name
	String phase_on_avg_spline_name = "phase_on_" + "avg_spline_" + series_name
	String phase_error_on_avg_spline_name = "phase_error_on_" + "avg_spline_" + series_name
	String phase_error_on_avg_SS_sdev_name = "phase_error_on_" + "avg_SS_sdev_" + series_name
	String phase_on_all_spline_name = "phase1k2k_on_all_spline_"+ series_name
	
	String phase_on_all_name = "phase1k2k_on_all_"+ series_name
	String phase_error_on_all_name = "phase_error1k2k_on_all_"+series_name
	
	String contrast_on_all_name = "contrast1k2k_on_all_"+ series_name
	String contrast_error_on_all_name = "contrast_error1k2k_on_all_"+series_name
	
	String wavelength_pad_name = "wavelength_pad_"+series_name
	
	Make/O/D/N=(numsets) $contrast_on_all_name, $contrast_error_on_all_name, $phase_on_avg_name, $phase_error_on_avg_name, $phase_error_on_avg_sdev_name, $phase_on_avg_spline_name, $phase_error_on_avg_spline_name, $phase_error_on_avg_SS_sdev_name
	Make/O/D/N=0 $phase_on_all_name, $phase_error_on_all_name, $wavelength_pad_name, $phase_on_all_spline_name
	
	Wave phase_on_avg = $phase_on_avg_name
	Wave phase_error_on_avg = $phase_error_on_avg_name
	Wave phase_error_on_avg_sdev = $phase_error_on_avg_sdev_name
	Wave phase_on_all = $phase_on_all_name
	Wave phase_error_on_all = $phase_error_on_all_name
	Wave contrast_on_all = $contrast_on_all_name
	Wave contrast_error_on_all = $contrast_error_on_all_name
	Wave wavelength_pad = $wavelength_pad_name
	Wave phase_on_spline_avg = $phase_on_avg_spline_name
	Wave phase_error_on_spline_avg = $phase_error_on_avg_spline_name
	Wave phase_error_on_spline_avg_sdev = $phase_error_on_avg_SS_sdev_name
	Wave phase_on_all_spline = $phase_on_all_spline_name
	
	Duplicate/O phase_on_avg phase_on_noPnorm_avg, phase_error_on_noPnorm_avg
	
	String phase_on_name = "phase1k2k_on_"+series_name
	String phase_error_on_name = "phase_error1k2k_on_"+series_name
	Wave phase_on =$phase_on_name
	Wave phase_error_on = $phase_error_on_name
	String phase_on_spline_name = "phase1k2k_on_spline_"+series_name
	Wave phase_on_spline = $phase_on_spline_name

	String contrast_on_name = "contrast1k2k_on_"+series_name
	String contrast_error_on_name = "contrast_error1k2k_on_"+series_name
	Wave contrast_on =$contrast_on_name
	Wave contrast_error_on = $contrast_error_on_name
	
	// do the extraction
	Variable n
	For(n=0; n<numsets; n+=1)
		Make/O/N=0 thisSet, thisErrSet, thisWaveSet, thisSplineSet, thisConSet, thisConErrSet
		
		Extract phase_on, thisSet, ( p>(extract_start[n]-1) && p<(extract_end[n]-1) )
		Extract phase_error_on, thisErrSet, ( p>(extract_start[n]-1) && p<(extract_end[n]-1) )
		
		Concatenate/NP {thisSet}, phase_on_all
		Concatenate/NP {thisErrSet}, phase_error_on_all
		
		Wavestats/Q thisSet
		phase_on_avg[n] = V_avg
		phase_error_on_avg[n] = v_sem
		phase_error_on_avg_sdev[n] = v_sdev
		
		Extract wavelengths, thisWaveSet, ( p>(extract_start[n]-1) && p<(extract_end[n]-1) )
		Wavestats/Q thisWaveSet
		wavelength[n] = v_avg
		Concatenate/NP {thisWaveSet}, wavelength_pad

		Extract phase_on_spline, thisSplineSet, ( p>(extract_start[n]-1) && p<(extract_end[n]-1) )
		Wavestats/Q thisSplineSet
		phase_on_spline_avg[n] = V_avg
		phase_error_on_spline_avg[n] = v_sem
		phase_error_on_spline_avg_sdev[n] = v_sdev
		Concatenate/NP {thisSplineSet}, phase_on_all_spline
		
		
		Extract contrast_on, thisConSet, ( p>(extract_start[n]-1) && p<(extract_end[n]-1) )
		Extract contrast_error_on, thisConErrSet, ( p>(extract_start[n]-1) && p<(extract_end[n]-1) )
		
		Concatenate/NP {thisConSet}, contrast_on_all
		Concatenate/NP {thisConErrSet}, contrast_error_on_all
		
//		Wavestats/Q thisSet
//		phase_on_noPnorm_avg_br[n] = V_avg
//		phase_error_on_noPnorm_avg_br[n] = v_sem
	EndFor
	
	
	String mask_all_name = "mask_all_"+series_name
	String mask_all_spline_name = "mask_all_spline_"+series_name
	String mask_name = "mask_"+series_name
	Duplicate/o phase_on_all $mask_all_name, $mask_all_spline_name
	Duplicate/O phase_on_avg $mask_name
	Wave mask_all = $mask_all_name
	Wave mask_all_spline = $mask_all_spline_name
	Wave mask = $mask_name
	mask_all=1
	mask_all_spline=1
	mask=1
	
	
	// It's all about displaying the results from here on.
	
	String phase_on_avg_nameTrace2 = phase_on_avg_name + "#1"
	
	Display/K=1/W=(550,450,1250,900) phase_on_avg vs wavelength
	ModifyGraph mode=3
	ModifyGraph marker=8
	//	if(normalizePowerBool)
	//		Label left "phase shift (rad/W)"
	//	else
	//		Label left "phase shift (rad)"
	//	endif
	Label left "norm. phase shift (rad/W)"
	Label bottom "wavelength (nm)"
	ModifyGraph zero(left)=1
	ErrorBars $phase_on_avg_name Y,wave=($phase_error_on_avg_name,$phase_error_on_avg_name)
	ModifyGraph nticks(bottom)=10,minor(bottom)=1
	AppendToGraph phase_on_avg vs wavelength
	ModifyGraph mode($phase_on_avg_nameTrace2) =2, rgb($phase_on_avg_nameTrace2)=(1,12815,52428)
	ErrorBars $phase_on_avg_nameTrace2 Y,wave=($phase_error_on_avg_sdev_name,$phase_error_on_avg_sdev_name)
	WindowNamer("Phase avg vs wavelength " + series_name)
	
	
	Display/K=1/W=(250,150,950,600) phase_on_all vs wavelength_pad
	ModifyGraph mode=3,marker=8
	ErrorBars $phase_on_all_name Y,wave=($phase_error_on_all_name,$phase_error_on_all_name)
	Label left "norm. phase shift (rad/W)"
	Label bottom "wavelength (nm)"
	ModifyGraph zero(left)=1
	ModifyGraph nticks(bottom)=10,minor(bottom)=1
	ModifyGraph zColor($phase_on_all_name)={$mask_all_name,0,1,RedWhiteBlue,1}
	WindowNamer("Phase all vs wavelength " + series_name)
	
	
	String phase_on_avg_nameTrace2SS = phase_on_avg_spline_name + "#1"
	
	Display/K=1/W=(550,450,1250,900) phase_on_spline_avg vs wavelength
	ModifyGraph mode=3
	ModifyGraph marker=8
	//	if(normalizePowerBool)
	//		Label left "phase shift (rad/W)"
	//	else
	//		Label left "phase shift (rad)"
	//	endif
	Label left "norm. phase shift (rad/W)"
	Label bottom "wavelength (nm)"
	ModifyGraph zero(left)=1
	ErrorBars $phase_on_avg_spline_name Y,wave=($phase_error_on_avg_spline_name,$phase_error_on_avg_spline_name)
	ModifyGraph nticks(bottom)=10,minor(bottom)=1
	AppendToGraph phase_on_spline_avg vs wavelength
	ModifyGraph mode($phase_on_avg_nameTrace2SS) =2, rgb($phase_on_avg_nameTrace2SS)=(1,12815,52428)
	ErrorBars $phase_on_avg_nameTrace2SS Y,wave=($phase_error_on_avg_SS_sdev_name,$phase_error_on_avg_SS_sdev_name)
	WindowNamer("Phase avg SS vs wavelength " + series_name)
	
	
	Display/K=1/W=(250,150,950,600) phase_on_all_spline vs wavelength_pad
	ModifyGraph mode=3,marker=8
	ErrorBars $phase_on_all_spline_name Y,wave=($phase_error_on_all_name,$phase_error_on_all_name)
	Label left "norm. phase shift (rad/W)"
	Label bottom "wavelength (nm)"
	ModifyGraph zero(left)=1
	ModifyGraph nticks(bottom)=10,minor(bottom)=1
	ModifyGraph zColor($phase_on_all_spline_name)={$mask_all_spline_name,0,1,RedWhiteBlue,1}
	WindowNamer("Phase all SS vs wavelength " + series_name)
	
	
	
	Display/K=1/W=(250,150,950,600) contrast_on_all vs wavelength_pad
	ModifyGraph mode=3,marker=8
	ErrorBars $contrast_on_all_name Y,wave=($contrast_error_on_all_name,$contrast_error_on_all_name)
	Label left "contrast"
	Label bottom "wavelength (nm)"
	ModifyGraph zero(left)=1
	ModifyGraph nticks(bottom)=10,minor(bottom)=1
	ModifyGraph zColor($contrast_on_all_name)={$mask_all_name,0,1,RedWhiteBlue,1}
	WindowNamer("contrast vs wavelength " + series_name)
End



// Will's most advanced method to process the MZW data. Recommended over the others.
Function avgPhases(seriesName)
	String seriesName
	
	String wavelengthsName = "wavelengths_"+seriesName
	Wave/Z wavelengths = $wavelengthsName
	
	if(waveExists(wavelengths)==0)
		getWavelengthData()
		Wave wavelengths = $wavelengthsName
	endif
	
	Variable minWavelength = wavemin(wavelengths)
	Variable maxWavelength = wavemax(wavelengths)
	//print maxwavelength-minwavelength
	
	Variable wavelengthBinWidth = 1		//units of pm
	wavelengthBinWidth*=1e-3
	
	Variable wavelengthBinnedPnts = ceil((maxWavelength-minWavelength)/wavelengthBinWidth)
	//print wavelengthbinnedallpnts
	
	String wavelengthsBinnedName = "wavelengthsBinned_"+seriesName
	String phaseBinnedName = "phaseBinned_"+seriesName
	String phaseBinnedStdErrName = "phaseBinnedStdErr_"+seriesName
	String phaseBinnedStdDevName = "phaseBinnedStdDev_"+seriesName
	Make/D/O/N=(wavelengthBinnedPnts) $wavelengthsBinnedName, $phaseBinnedName, $phaseBinnedStdErrName, $phaseBinnedStdDevName
	Wave wavelengthsBinned = $wavelengthsBinnedName
	Wave phaseBinned = $phaseBinnedName
	Wave phaseBinnedStdErr = $phaseBinnedStdErrName
	Wave phaseBinnedStdDev = $phaseBinnedStdDevName
	
	String phaseBinnedSSName = "phaseBinnedSS_"+seriesName
	String phaseBinnedSSStdErrName = "phaseBinnedSSStdErr_"+seriesName
	String phaseBinnedSSStdDevName = "phaseBinnedSSStdDev_"+seriesName
	Make/D/O/N=(wavelengthBinnedPnts) $phaseBinnedSSName, $phaseBinnedSSStdErrName, $phaseBinnedSSStdDevName
	Wave phaseBinnedSS = $phaseBinnedSSName
	Wave phaseBinnedSSStdErr = $phaseBinnedSSStdErrName
	Wave phaseBinnedSSStdDev = $phaseBinnedSSStdDevName
	
	
	// initialize the waves.
	wavelengthsBinned = minWavelength + p*wavelengthBinWidth
	phaseBinned = nan
	phaseBinnedStdErr = nan
	phaseBinnedStdDev = nan
	
	phaseBinnedSS = nan
	phaseBinnedSSStdErr = nan
	phaseBinnedSSStdDev = nan
	
	// wavelengths with fewer files than this won't be included in the final data.
	Variable filesPerWavelengthThreshold = 5	//3 for 10 s x 20 files/lamba //5 for 5s x 40 files/lambda
	
	Variable filterType = 2		// 0 for no filtering, 1 for std dev filtering, 2 for trimmed mean
	Variable phaseSetStdDevFilter = 2
	Variable trimmedMeanPercent = 20
	
	
//	String phase_on_name = "phase1k2k_on_"+ seriesName
	String phase_on_name = "phase_on_Pnorm_" + seriesName
	Wave phase_on = $phase_on_name
	
//	String phase_on_ss_name = "phase1k2k_on_spline_"+ seriesName
	String phase_on_ss_name = "phase_on_SS_Pnorm_" + seriesName
	Wave phase_on_ss = $phase_on_ss_name
	
//	String phase_ref_name = "phase1k2k_ref_"+ seriesName
	String phase_ref_name = "phase1k2k_ref_"+ seriesName
	Wave phase_ref = $phase_ref_name
	Wavestats/Q phase_ref
	Variable phase_ref_StdDev = v_sdev
	
//	String phase_ref_ss_name = "phase1k2k_ref_spline_"+ seriesName
	String phase_ref_ss_name = "phase1k2k_ref_spline_"+ seriesName
	Wave phase_ref_ss = $phase_ref_ss_name
	Wavestats/Q phase_ref_ss
	Variable phase_ref_ss_StdDev = v_sdev
	
	Variable useRefStdDev=0	// don't use this feature with power normalization. I can't remember why it's here a year later...
	
	// Call the binning function. Do it for each kind of reference phase fit.
	Binner(phase_on, phaseBinned, phaseBinnedStdErr, phaseBinnedStdDev, wavelengths, wavelengthsBinned, wavelengthBinWidth, filterType, phaseSetStdDevFilter, filesPerWavelengthThreshold, trimmedMeanPercent, phase_ref_StdDev, useRefStdDev)
	Binner(phase_on_ss, phaseBinnedSS, phaseBinnedSSStdErr, phaseBinnedSSStdDev, wavelengths, wavelengthsBinned, wavelengthBinWidth, filterType, phaseSetStdDevFilter, filesPerWavelengthThreshold, trimmedMeanPercent, phase_ref_ss_StdDev, useRefStdDev)
	
	
	// Make plots of processed phase vs wavelength
	
	String phaseBinnedNameTrace2 = phaseBinnedName + "#1"
	
	Display/K=1/W=(550,450,1250,900) phaseBinned vs wavelengthsBinned
	ModifyGraph mode=3
	ModifyGraph marker=8
	Label left "norm. phase shift (rad/W)"
	Label bottom "wavelength (nm)"
	ModifyGraph zero(left)=1
	ErrorBars $phaseBinnedName Y,wave=($phaseBinnedStdErrName,$phaseBinnedStdErrName)
	ModifyGraph nticks(bottom)=10,minor(bottom)=1
	AppendToGraph phaseBinned vs wavelengthsBinned
	ModifyGraph mode($phaseBinnedNameTrace2) =2, rgb($phaseBinnedNameTrace2)=(1,12815,52428)
	ErrorBars $phaseBinnedNameTrace2 Y,wave=($phaseBinnedStdDevName,$phaseBinnedStdDevName)
	WindowNamer("Phase binned " + seriesName)
	
	
	String phaseBinnedSSNameTrace2 = phaseBinnedSSName + "#1"
	
	Display/K=1/W=(550,450,1250,900) phaseBinnedSS vs wavelengthsBinned
	ModifyGraph mode=3
	ModifyGraph marker=8
	Label left "norm. phase shift (rad/W)"
	Label bottom "wavelength (nm)"
	ModifyGraph zero(left)=1
	ErrorBars $phaseBinnedSSName Y,wave=($phaseBinnedSSStdErrName,$phaseBinnedSSStdErrName)
	ModifyGraph nticks(bottom)=10,minor(bottom)=1
	AppendToGraph phaseBinnedSS vs wavelengthsBinned
	ModifyGraph mode($phaseBinnedSSNameTrace2) =2, rgb($phaseBinnedSSNameTrace2)=(1,12815,52428)
	ErrorBars $phaseBinnedSSNameTrace2 Y,wave=($phaseBinnedSSStdDevName,$phaseBinnedSSStdDevName)
	WindowNamer("Phase binned SS " + seriesName)
End



// This does the hard work of binning the data. 
Function Binner(phase_on_Active, phaseBinned_Active, phaseBinnedStdErr_Active, phaseBinnedStdDev_Active, wavelengths, wavelengthsBinned, wavelengthBinWidth, filterType, phaseSetStdDevFilter, filesPerWavelengthThreshold, trimmedMeanPercent, refStdDev, useRefStdDev)
	Wave phase_on_Active, phaseBinned_Active, phaseBinnedStdErr_Active, phaseBinnedStdDev_Active, wavelengths, wavelengthsBinned
	Variable wavelengthBinWidth, filterType, phaseSetStdDevFilter, filesPerWavelengthThreshold, trimmedMeanPercent, refStdDev, useRefStdDev
	
	// data outside this range will be discarded.
	Variable minCutoffWavelength = 768.970-0.100
	Variable maxCutoffWavelength = 768.970+0.100
	
	// For each wavelength bin, find the relavant phase data, filter it, and then assign the result to the output wave
	Variable i
	For(i=0; i<numpnts(wavelengthsBinned); i+=1)
		Make/O/D/N=0 thisPhaseSet
		Extract/O phase_on_Active, thisPhaseSet, ( wavelengths[p] > (wavelengthsBinned[i] -wavelengthBinWidth/2) && wavelengths[p] <= (wavelengthsBinned[i]+wavelengthBinWidth/2) && wavelengths[p] > minCutoffWavelength && wavelengths[p] < maxCutoffWavelength)
		
		if(numpnts(thisPhaseSet)!=0)
			
			WaveStats/Z/Q thisPhaseSet
			
			// do the filtering.
			if(filterType==0)
				//do nothing -- no filtering
			elseif(filterType==1)	// std. dev. filtering
				thisPhaseSet = (thisPhaseSet[p] > v_avg-phaseSetStdDevFilter*v_sdev && thisPhaseset[p] < v_avg+phaseSetStdDevFilter*v_sdev) ? thisPhaseSet[p] : NaN
				Wavestats/Z/Q thisPhaseSet
			elseif(filterType==2)	// trimed mean computation
				Sort thisPhaseSet, thisPhaseSet		// sorts from lowest to highest, puts nans at the end
				Variable thisPhaseSetPnts = v_npnts	// number of pnts that are not nan's or inf's
				Variable pntsCut = round(thisPhaseSetPnts*trimmedMeanPercent/100)	// number of pnts to cut from each end of the set. round up.
				if(pntsCut==0)	// always throw out at least one point high and low.
					pntsCut=1
				endif
				//print "trimmed points from each end of the set: ", pntsCut
				thisPhaseSet[0,pntsCut-1] = NaN				//NaNs the lowest n=pntsCut pnts
				thisPhaseSet[v_npnts-pntsCut,v_npnts] = NaN	//NaNs the largest n=pntsCut pnts. the wave may already have NaNed points at the end of it due to the previous sort operation, so we need to be careful to start at the right place
				Wavestats/Z/Q thisPhaseSet
			endif
				
			// do the final wave assignments	
			if(v_npnts>= filesPerWavelengthThreshold)	
				phaseBinned_Active[i] = v_avg
				if(useRefStdDev)
					phaseBinnedStdErr_Active[i] = refStdDev/sqrt(v_npnts)
					phaseBinnedStdDev_Active[i] = refStdDev
				else
					phaseBinnedStdErr_Active[i] = v_sem//0.00889148/sqrt(v_npnts)//v_sem
					phaseBinnedStdDev_Active[i] = v_sdev//0.00889148//v_sdev
				endif
			endif

		endif
	EndFor
End





Function extractMessy(series_name)
	String series_name

	String wavelength_name = "wavelength_"+series_name
	Wave wavelength = $wavelength_name
	Variable numsets=numpnts(wavelength)
	
	String phase_on_avg_name = "phase1k2k_on_" + "avg_" + series_name
	String phase_error_on_avg_name = "phase_error1k2k_on_" + "avg_" + series_name
	String phase_error_on_avg_sdev_name = "phase_error1k2k_on_" + "avg_sdev_" + series_name
	String phase_on_all_name = "phase1k2k_on_all_"+ series_name
	String phase_error_on_all_name = "phase_error1k2k_on_all_"+series_name
	String wavelength_pad_name = "wavelength_pad_"+series_name
	
	Make/O/D/N=(numsets) $phase_on_avg_name, $phase_error_on_avg_name, $phase_error_on_avg_sdev_name 
	Make/O/D/N=0 $phase_on_all_name, $phase_error_on_all_name, $wavelength_pad_name
	
	Wave phase_on_avg = $phase_on_avg_name
	Wave phase_error_on_avg = $phase_error_on_avg_name
	Wave phase_error_on_avg_sdev = $phase_error_on_avg_sdev_name
	Wave phase_on_all = $phase_on_all_name
	Wave phase_error_on_all = $phase_error_on_all_name
	Wave wavelength_pad = $wavelength_pad_name
	
	
	Duplicate/O phase_on_avg phase_on_noPnorm_avg, phase_error_on_noPnorm_avg
	
	String phase_on_name = "phase1k2k_on_"+series_name
	String phase_error_on_name = "phase_error1k2k_on_"+series_name
	Wave phase_on =$phase_on_name
	Wave phase_error_on = $phase_error_on_name
	
	String extract_start_name = "extract_start_"+series_name
	String extract_end_name = "extract_end_"+series_name
	Wave extract_start = $extract_start_name
	Wave extract_end = $extract_end_name

	
	Variable n
	For(n=0; n<numsets; n+=1)
		Make/O/N=0 thisSet, thisErrSet	
		
		Extract phase_on, thisSet, ( p>(extract_start[n]-1) && p<(extract_end[n]-1) )
		Extract phase_error_on, thisErrSet, ( p>(extract_start[n]-1) && p<(extract_end[n]-1) )
		
		Wavestats/Q thisSet
		phase_on_avg[n] = V_avg
		phase_error_on_avg[n] = v_sem
		phase_error_on_avg_sdev[n] = v_sdev
		
		Concatenate/NP {thisSet}, phase_on_all
		Concatenate/NP {thisErrSet}, phase_error_on_all
		
		
		Make/O/N=(numpnts(thisSet)) thisWavelength = wavelength[n]
		Concatenate/NP {thisWavelength}, wavelength_pad
		
//		Wavestats/Q thisSet
//		phase_on_noPnorm_avg_br[n] = V_avg
//		phase_error_on_noPnorm_avg_br[n] = v_sem
	EndFor
	
	String phase_on_avg_nameTrace2 = phase_on_avg_name + "#1"
	
	Display/K=1/W=(550,450,1250,900) phase_on_avg vs wavelength
	ModifyGraph mode=3
	ModifyGraph marker=8
	//	if(normalizePowerBool)
	//		Label left "phase shift (rad/W)"
	//	else
	//		Label left "phase shift (rad)"
	//	endif
	Label left "norm. phase shift (rad/W)"
	Label bottom "wavelength (nm)"
	ModifyGraph zero(left)=1
	ErrorBars $phase_on_avg_name Y,wave=($phase_error_on_avg_name,$phase_error_on_avg_name)
	ModifyGraph nticks(bottom)=10,minor(bottom)=1
	AppendToGraph phase_on_avg vs wavelength
	ModifyGraph mode($phase_on_avg_nameTrace2) =2, rgb($phase_on_avg_nameTrace2)=(1,12815,52428)
	ErrorBars $phase_on_avg_nameTrace2 Y,wave=($phase_error_on_avg_sdev_name,$phase_error_on_avg_sdev_name)
	WindowNamer("Phase avg vs wavelength " + series_name)
	
	
	Display/K=1/W=(250,150,950,600) phase_on_all vs wavelength_pad
	ModifyGraph mode=3,marker=8
	ErrorBars $phase_on_all_name Y,wave=($phase_error_on_all_name,$phase_error_on_all_name)
	Label left "norm. phase shift (rad/W)"
	Label bottom "wavelength (nm)"
	ModifyGraph zero(left)=1
	ModifyGraph nticks(bottom)=10,minor(bottom)=1
	WindowNamer("Phase all vs wavelength " + series_name)
End




Function getWavelengthData()
	String/G directory_name// = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:120507:fringe data:"

	string/G series_name// = "e"
	
	NewPath/O/Q myPath directory_name						// create a path for the desired directory
	String fileList = IndexedFile(myPath, -1, "????")			// create a list containing all .txt files in the specified symbolic path.
	fileList = GrepList(fileList, "^"+series_name+"\\d*_wavelengths.txt", 0, ";")	// search the list for all files of the form "start of string + series_name + digit"
	//print fileList
	Variable numFiles = ItemsInList(fileList, ";")				// count the number of items in that list

	String wavelengthsName = "wavelengths_"+series_name
	Make/D/O/N=(numFiles) $wavelengthsName
	Wave wavelengths = $wavelengthsName
	
	String wavelengthsStdDevName = "wavelengthsStdDev_"+series_name
	Duplicate/O wavelengths $wavelengthsStdDevName
	Wave wavelengthsStdDev = $wavelengthsStdDevName
	
	String fileName
	Variable i
	For(i=1; i<=numFiles; i+=1)
		fileName = series_name+num2str(i)+"_wavelengths.txt"
		LoadWave/Q/G/D/W/N/O/P=myPath fileName
		Wave wavelength
		
		Wavestats/Q wavelength
		wavelengths[i-1] = v_avg
		wavelengthsStdDev[i-1] = v_sdev
	EndFor
	
	Variable stdDevFilter = 500e-6
	
	wavelengths = (wavelengthsStdDev > stdDevFilter )? NaN : wavelengths
	
	String fileIndexName = "file_index_"+series_name
	Wave fileIndex = $fileIndexName
	
	Display/K=1 wavelengths vs fileIndex
	WindowNamer("Wavelength data " + series_name)
End




Function getEtalonData()
	String/G directory_name = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:120118:fringe data:"

	string/G series_name = "c"
	
	NewPath/O/Q myPath directory_name						// create a path for the desired directory
	String fileList = IndexedFile(myPath, -1, "????")			// create a list containing all .txt files in the specified symbolic path.
	fileList = GrepList(fileList, "^"+series_name+"_wavelength_"+"\\d", 0, ";")	// search the list for all files of the form "start of string + series_name + digit"
	Variable numFiles = ItemsInList(fileList, ";")				// count the number of items in that list
	
	String this10setGrepStr
	String this10set
	String lastInThis10set
	Variable i
	
	
	String thisFile, thisFileNumberStr
	Variable thisFileNumber
	Make/O/N=(numFiles) etalonFileNums
	String regExp = "^[[:alpha:]]+_wavelength_(\\d+)$"
	
	For(i=0; i<numFiles; i+=1)
		thisFile = StringFromList(i,fileList,";")
		SplitString/E=(regExp) thisFile, thisFileNumberStr
		thisFileNumber=str2num(thisFileNumberStr)
		etalonFileNums[i]=thisFileNumber
	EndFor
	
	Sort etalonFileNums,etalonFileNums
	
	Variable lastEtalonFileNum = etalonFileNums[numFiles]
	
	Variable numFilesPerLayout = 12
	Variable layoutNum = 2	//used if additional layouts are needed
	String thisEtalonWindowName = "WNEtalonData"+series_name+"1_1"
	String etalonSummaryLayoutName = "Etalon Summary " + series_name
	NewLayout/K=1 as etalonSummaryLayoutName
	
	loadAndDisplayEtalonData(directory_Name, series_name, 1)
	
	AppendLayoutObject graph $thisEtalonWindowName
	
	For(i=1; i<numFiles; i+=1)
//		if(i==0)
//			this10setGrepStr = "_\\d$"
//		else
//			this10setGrepStr = "_"+num2str(i*20/10)+"\\d$"
//		endif
//			
//		this10set = GrepList(fileList, this10setGrepStr,0,";")
//		
//		lastInThis10set = StringFromList(ItemsInList(this10set)-1,this10set,";")
		
		thisFileNumber = etalonFileNums[i]
		loadAndDisplayEtalonData(directory_Name, series_name, thisFileNumber)
		
		thisEtalonWindowName = "WNEtalonData"+series_name+num2str(etalonFileNums[i])+"_1"
		
		if(mod(i,numFilesPerLayout)==0)
			Execute "Tile"	//Tiles the previous layout
			etalonSummaryLayoutName = "Etalon Summary " + series_name + " cont. (" + num2str(layoutNum) + ")"
			NewLayout/K=1 as etalonSummaryLayoutName	//starts a new layout
			layoutNum+=1
		endif	
		
		AppendLayoutObject graph $thisEtalonWindowName
		
		//print this10set
		
//		print StringFromList(ItemsInList(this10set)-1,this10set,";")

	//	For(j=0; j<ItemsInList(this10set,";"); j+=1)
	//		print StringFromList(j,this10set,";")
	//	EndFor
	
	EndFor
	
	Execute "Tile"

End





Function loadAndDisplayEtalonData(directoryName, seriesName, fileNumber)
	String directoryName, seriesName
	Variable fileNumber
	
//	SplitString/E=(regExp) fileName, seriesName, fileNumberStr
	
	String fileNumberStr = num2str(fileNumber)
	
	NewPath/O/Q myPath directoryName
	String fileName = seriesName+"_wavelength_"+fileNumberStr
	LoadWave/Q/J/D/W/N/O/K=0/L={0,5,0,1,0} /P=myPath fileName
	
	String gs_name = "grating_spec_"+seriesName+fileNumberStr
	String sync_name = "sync_"+seriesName+fileNumberStr
	String abs_pb_name = "abs_pb_"+seriesName+fileNumberStr
	String et100_name = "et100_"+seriesName+fileNumberStr
	String et15_name = "et15_"+seriesName+fileNumberStr
	String abs_nopb_name = "abs_nopb_"+seriesName+fileNumberStr
	String fp_name = "fp_"+seriesName+fileNumberStr
	
	Duplicate/O X_0_ $gs_name
	Duplicate/O X_1_ $abs_nopb_name
	Duplicate/O X_2_ $et15_name
	Duplicate/O X_3_ $et100_name
	Duplicate/O X_4_ $fp_name
	Duplicate/O X_5_ $sync_name
	
	Wave thisGS_wave = $gs_name
	Wave thisSync_wave =$sync_name
	//Wave thisAbs_pb_wave= $abs_pb_name
	Wave thisEt100_wave =$et100_name
	Wave thisEt15_wave =$et15_name
	Wave thisAbs_nopb_wave= $abs_nopb_name
	Wave thisFP_wave = $fp_name
	
	thisAbs_nopb_wave*=4
	thisFP_wave*=10
	
	Variable etPnts = numpnts(thisSync_wave)
	//SetScale/P x 0,0.1,"ms", thisGS_wave, thisSync_wave, thisAbs_pb_wave, thisEt100_wave, thisEt15_wave, thisAbs_nopb_wave
	
	FindLevel/Q/EDGE=2 thisSync_wave, 2
	Variable TTLlow_start = ceil(v_levelx)
	FindLevel/Q/EDGE=1/R=[TTLlow_start, etPnts] thisSync_wave, 2
	Variable TTLlow_end = floor(v_levelx)
	
	Variable numPntsToDelete = etPnts-TTLlow_end
	DeletePoints TTLlow_end, numPntsToDelete, thisGS_wave, thisSync_wave, thisEt100_wave, thisEt15_wave, thisAbs_nopb_wave, thisFP_wave
	DeletePoints 0, TTLlow_start, thisGS_wave, thisSync_wave, thisEt100_wave, thisEt15_wave, thisAbs_nopb_wave, thisFP_wave
	
	Display/K=1 thisGS_wave, thisSync_wave, thisFP_wave, thisEt100_wave, thisEt15_wave, thisAbs_nopb_wave
	//ModifyGraph lsize=2
	SetAxis left 0,7
	ModifyGraph rgb($gs_name)=(0,0,0),rgb($sync_name)=(1,4,52428),rgb($fp_name)=(0,43690,65535),rgb($et100_name)=(2,39321,1),rgb($et15_name)=(43690,43690,43690)
	TextBox/C/N=text0/F=0/A=MC seriesName+fileNumberStr
	TextBox/C/N=text0/A=RT/X=13.50/Y=7.01
	WindowNamer("Etalon Data "+seriesName+fileNumberStr)
End






Function loadAndDisplayEtalonDataOld(directoryName, seriesName, fileNumber)
	String directoryName, seriesName
	Variable fileNumber
	
//	SplitString/E=(regExp) fileName, seriesName, fileNumberStr
	
	String fileNumberStr = num2str(fileNumber)
	
	NewPath/O/Q myPath directoryName
	String fileName = seriesName+"_wavelength_"+fileNumberStr
	LoadWave/Q/J/D/W/N/O/K=0/L={0,5,0,1,0} /P=myPath fileName
	
	String gs_name = "grating_spec_"+fileNumberStr
	String sync_name = "sync_"+fileNumberStr
	String abs_pb_name = "abs_pb_"+fileNumberStr
	String et100_name = "et100_"+fileNumberStr
	String et15_name = "et15_"+fileNumberStr
	String abs_nopb_name = "abs_nopb_"+fileNumberStr
	String fp_name = "fp_"+fileNumberStr
	

	Duplicate/O X_0_ $gs_name
	Duplicate/O X_1_ $sync_name
	Duplicate/O X_2_ $abs_pb_name
	Duplicate/O X_3_ $et100_name
	Duplicate/O X_4_ $et15_name
	Duplicate/O X_5_ $abs_nopb_name
	
	Wave thisGS_wave = $gs_name
	Wave thisSync_wave =$sync_name
	Wave thisAbs_pb_wave= $abs_pb_name
	Wave thisEt100_wave =$et100_name
	Wave thisEt15_wave =$et15_name
	Wave thisAbs_nopb_wave= $abs_nopb_name
	
	thisAbs_nopb_wave*=4
	
	Variable etPnts = numpnts(thisSync_wave)
	//SetScale/P x 0,0.1,"ms", thisGS_wave, thisSync_wave, thisAbs_pb_wave, thisEt100_wave, thisEt15_wave, thisAbs_nopb_wave
	
	FindLevel/Q/EDGE=2 thisSync_wave, 2
	Variable TTLlow_start = ceil(v_levelx)
	FindLevel/Q/EDGE=1/R=[TTLlow_start, etPnts] thisSync_wave, 2
	Variable TTLlow_end = floor(v_levelx)
	
	Variable numPntsToDelete = etPnts-TTLlow_end
	DeletePoints TTLlow_end, numPntsToDelete, thisGS_wave, thisSync_wave, thisAbs_pb_wave, thisEt100_wave, thisEt15_wave, thisAbs_nopb_wave
	DeletePoints 0, TTLlow_start, thisGS_wave, thisSync_wave, thisAbs_pb_wave, thisEt100_wave, thisEt15_wave, thisAbs_nopb_wave
	
	Display/K=1 thisGS_wave, thisSync_wave, thisAbs_pb_wave, thisEt100_wave, thisEt15_wave, thisAbs_nopb_wave
	//ModifyGraph lsize=2
	SetAxis left 0,11
	ModifyGraph rgb($gs_name)=(0,0,0),rgb($sync_name)=(1,4,52428),rgb($abs_pb_name)=(0,43690,65535),rgb($et100_name)=(2,39321,1),rgb($et15_name)=(43690,43690,43690)
	TextBox/C/N=text0/F=0/A=MC fileNumberStr
	TextBox/C/N=text0/A=RT/X=13.50/Y=7.01
	WindowNamer("Etalon Data "+seriesName+fileNumberStr)
End











Function etalonTheory()
	Variable pntsPerCycle = 100
	
	Variable f0D1 = c/KD1wavelength	//d1 freq in hz
	Variable f0D2 = c/KD2wavelength	//d2 freq
	
	Variable d1d2diff = (f0D2-f0D1)/1e9	//d1 d2 freq diff in ghz
	//print/D d1d2diff
	
	Variable minf = -50
	Variable maxf = d1d2diff+50
	Variable minwave = KD2wavelength*.999*1e9
	Variable maxwave = KD1wavelength*1.001*1e9
	
	Variable/G et100f = 97.8
	Variable/G et15f = 15.855//15.855 (109.125 fringes) works. also 16.00 (108.125) and 15.71 (110.125)
	Variable/G et50f = 36.81
	Variable/G et1571f = 15.71
	Variable/G et16f = 16.00
	Variable/G fp15f = 1.500
	
	Variable/G et100phase = pi/2*1.05
	Variable/G et15phase = -pi*3/8
	Variable/G et50phase = pi*3/4*1.1
	Variable/G et1571phase = et15phase
	Variable/G et16phase = et15phase
	Variable/G fp15phase = pi/2
	
	Variable et100cycles = (d1d2diff+100)/et100f
	Variable et50cycles =  (d1d2diff+100)/et50f
	Variable et15cycles =  (d1d2diff+100)/et15f
	Variable et1571cycles =  (d1d2diff+100)/et1571f
	Variable et16cycles = (d1d2diff+100)/et16f
	Variable fp15cycles = (d1d2diff+100)/fp15f
	
	Variable et100pnts = round(pntsPerCycle*et100cycles)
	Variable et50pnts = round(pntsPerCycle*et50cycles)
	Variable et15pnts = round(pntsPerCycle*et15cycles)
	Variable et1571pnts = round(pntsPerCycle*et1571cycles)
	Variable et16pnts = round(pntsPerCycle*et16cycles)
	Variable fp15pnts = round(pntsPerCycle*fp15cycles)
	
	Make/O/N=(et100pnts) et100
	Make/O/N=(et50pnts) et50
	Make/O/N=(et15pnts) et15
	Make/O/N=(et1571pnts) et1571
	Make/O/N=(et16pnts) et16
	Make/O/N=(fp15pnts) fp15

	SetScale/i x minf, maxf, et100, et50, et15, et1571, fp15, et16
	
	et100 = sin(2*pi*x/et100f+et100phase)
	et50 = sin(2*pi*x/et50f+et50phase)
	et15 = sin(2*pi*x/et15f+et15phase)
	et1571 = sin(2*pi*x/et1571f+et1571phase)
	et16 = sin(2*pi*x/et16f+et16phase)
	fp15 = sin(2*pi*x/fp15f+fp15phase)
	
	fp15 = fp15>=0 ? fp15 : 0
	fp15*=fp15
	fp15*=fp15
	fp15*=fp15
	
	if(0)
	fp15=0
	variable n
	for(n=1; n<=40;n+=1)
		fp15 += sin(n*2*pi*x/fp15f+fp15phase)
		fp15 += cos(n*2*pi*x/fp15f+fp15phase)
	endfor
	fp15/=(n-1)
	endif
	
	et50+=4
	et15-=3
	et1571 -=6
	et16-=2
	
	fp15-=6
	
	Make/O/N=2 fd1, fd2
	SetScale/p x 0,0.000000000001, fd1
	fd1={-6,5}
	SetScale/p x d1d2diff,0.000000000001, fd2
	//SetScale/i x d1d2diff,d1d2diff, fd2
	fd2={-6,5}
	
	
	Variable makeGraph = 0
	if(makeGraph)
		Display/W=(53,116,1371,731) et100,et50,et15,fp15,fd1,fd2
		ModifyGraph gfSize=14
		ModifyGraph lSize=2
		ModifyGraph rgb(fd1)=(0,0,0),rgb(fd2)=(0,0,0)
		ModifyGraph rgb(et15)=(43690,43690,43690)
		ModifyGraph rgb(et100)=(2,39321,1)
		ModifyGraph rgb(fp15)=(0,43690,65535)
		ModifyGraph zero(bottom)=1
		ModifyGraph nticks(bottom)=20
		ModifyGraph minor(bottom)=1
		ModifyGraph zeroThick(bottom)=2
		Label bottom "Detuning from D1 line (GHz)"
		ShowInfo
		Legend/C/N=text0/J/F=0/X=17.60/Y=18.01 "\\s(et50) 40-50 GHz\r\\s(et100) 97.8 GHz\r\\s(et15) 15.855 GHz\r\\s(fp15) 1.5 GHz FP"
		TextBox/C/N=text1/F=0/X=4.19/Y=22.06 "d2-d1 = 1730.1 GHz"
		TextBox/C/N=text2/F=0/X=92.68/Y=22.79 "D1 line"
		TextBox/C/N=text3/F=0/X=-1.64/Y=23.35 "D2 line"
	endif
End






Function safPol()
	Variable red_p12_D_ns = 4.102	//reduced matrix elements. unc (5)
	Variable red_p32_D_ns = 5.800	//(8)
	Variable deltaE_np12 = 0.059165	//
	Variable deltaE_np32 = 0.059428
	Variable A_ns = 6.26
	Variable deltaA_ns = 0.33

	Variable omega0D1 = 2*pi*c/KD1wavelength*(hbar/hartree)
	Variable omega0D2 = 2*pi*c/KD2wavelength*(hbar/hartree)
	Variable naturalLinewidth = 6.1e6*(hbar/hartree)
	
	Variable omegaPnts = 10000
	Variable minOmega = omega0D1*.999
	Variable maxOmega = omega0D2*1.001
	
	make/o/n=(omegaPnts) alpha_omega
	SetScale/i x minOmega, maxOmega, alpha_omega
	
	alpha_omega = 1/3 * deltaE_np12 *red_p12_D_ns^2/(deltaE_np12^2-x^2) + 1/3*deltaE_np32*red_p32_D_ns^2/(deltaE_np32^2-x^2) + A_ns		// in atomic units (a0^3)
	
	Duplicate/o alpha_omega phiSaf
	
	phiSaf = alpha_omega*AUtoSI
End



Function safPolf()
	Variable red_p12_D_ns = 4.102*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//reduced matrix elements. unc (5)
	Variable red_p32_D_ns = 5.800*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//(8)
	Variable deltaE_np12 = 0.059165*hartree	//
	Variable deltaE_np32 = 0.059428*hartree
	Variable A_ns = 6.26*AUtoSI
	Variable deltaA_ns = 0.33*AUtoSI
	
	Variable phiNorm = -3.5e34

	Variable f0D1 = c/KD1wavelength
	Variable f0D2 = c/KD2wavelength
	
	Variable fPnts = 50000
	Variable minf = f0D1*.999
	Variable maxf = f0D2*1.001
	Variable minwave = KD2wavelength*.999*1e9
	Variable maxwave = KD1wavelength*1.001*1e9
	
	make/o/n=(fPnts) alpha_f, alpha_wave
	SetScale/i x minf, maxf, alpha_f
	
	alpha_f = 1/3 * deltaE_np12 *red_p12_D_ns^2/(deltaE_np12^2-(x*h)^2) + 1/3*deltaE_np32*red_p32_D_ns^2/(deltaE_np32^2-(x*h)^2) + A_ns-10*deltaA_ns		// in atomic units (a0^3)
	
	Duplicate/o alpha_f phiSaf_f
	
	SetScale/i x minwave, maxwave, alpha_wave
	alpha_wave = 1/3 * deltaE_np12 *red_p12_D_ns^2/(deltaE_np12^2-((c/x*1e9)*h)^2) + 1/3*deltaE_np32*red_p32_D_ns^2/(deltaE_np32^2-((c/x*1e9)*h)^2) + A_ns-0*deltaA_ns
	
	alpha_f = abs(alpha_f)>2e-35 ? NaN : alpha_f
	alpha_wave = abs(alpha_wave)>2e-35 ? NaN : alpha_wave
	
	Duplicate/o alpha_f phiSaf_f
	Duplicate/o alpha_wave phiSaf_wave
	
	phiSaf_f = alpha_f*phiNorm
	phiSaf_wave = alpha_wave*phiNorm
End


////////////////////////

Function goFitTOWsafR(series_name,all)
	String series_name
	Variable all
	
	Variable oneKtwoK=1
	
	String onektwokStr = ""
	if(OnekTwok==1)
		onektwokStr = "1k2k"
	elseif(OnekTwok == 2)
		onektwokStr = "2k2k"
	endif
	
	String phase_on_avg_name = "phase" + onektwokStr + "_on_avg_" + series_name
	String phase_error_on_avg_name = "phase_error" + onektwokStr + "_on_avg_" + series_name
	String phase_error_on_avg_sdev_name = "phase_error" + onektwokStr + "_on_avg_sdev_" + series_name
	String phase_error_tot_name = "phase_error_tot_" + series_name
	String phase_on_name = "phase" + onektwokStr + "_on_all_" + series_name
	String phase_on_spline_name = "phase" + onektwokStr + "_on_all_spline_" + series_name
	String phase_error_on_name = "phase_error" + onektwokStr + "_on_all_" + series_name
	String wavelength_pad_name = "wavelength_pad_"+series_name
	String wavelength_name = "wavelength_"+series_name
	
	String maskWaveAvgName = "mask_avg_"+series_name
	String maskWaveAllName = "mask_all_"+series_name
	String maskWaveAllSplineName = "mask_all_spline_"+series_name
	
	String phase_bin_name = "PhaseBinned_"+series_name
	String phase_error_bin_name = "PhaseBinnedStdErr_"+series_name
	String wavelength_bin_name = "wavelengthsBinned_"+series_name
	String maskWaveBinnedName = "mask_binned_"+series_name

	String phase_binSS_name = "PhaseBinnedSS_"+series_name
	String phase_error_binSS_name = "PhaseBinnedSSStdErr_"+series_name
	String wavelength_binSS_name = "wavelengthsBinned_"+series_name
	String maskWaveBinnedSSName = "mask_binned_SS_"+series_name
	

	
	if(0)
		phase_error_on_avg_name = phase_error_tot_name
	endif
	
	
	Variable firstPntToFit
	Variable lastPntToFit
	
	// choose to fit the entire data set or the averaged data set
	//Variable all = 0
	if(all==0)
		Wave phase_on_fit = $phase_on_avg_name
		Wave phase_error_on_fit = $phase_error_on_avg_name
		if(0)
			Wave phase_error_on_fit = $phase_error_on_avg_sdev_name
		endif
		Wave wavelength_fit = $wavelength_name
		Wave/Z maskWave = $maskWaveAvgName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==1)
		Wave phase_on_fit = $phase_on_name
		Wave phase_error_on_fit = $phase_error_on_name
		Wave wavelength_fit = $wavelength_pad_name
		Wave/Z maskWave = $maskWaveAllName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==2)
		Wave phase_on_fit = $phase_on_spline_name
		Wave phase_error_on_fit = $phase_error_on_name
		Wave wavelength_fit = $wavelength_pad_name
		Wave/Z maskWave = $maskWaveAllSplineName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==3)
		Wave phase_on_fit = $phase_bin_name
		Wave phase_error_on_fit = $phase_error_bin_name
		Wave wavelength_fit = $wavelength_bin_name
		Wave/Z maskWave = $maskWaveBinnedName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==4)
		Wave phase_on_fit = $phase_binSS_name
		Wave phase_error_on_fit = $phase_error_binSS_name
		Wave wavelength_fit = $wavelength_binSS_name
		Wave/Z maskWave = $maskWaveBinnedSSName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0//315
		lastPntToFit=1000//534
	else
		print "all, binned, ss?"
		abort
	endif
	
	Make/O/D w_coef = {2,1e36}
	
	Variable suppressUpdateYesNo=0
	
	Duplicate/O phase_error_on_fit phase_error_on_fit2
	phase_error_on_fit2*=2
	
	Make/O epsWave = {1e-4,1e30}

	
	//FuncFit/N=(suppressUpdateYesNo) FitTOWsaf w_coef phase_on_fit[firstPntToFit,lastPntToFit] /X=wavelength_fit /W=phase_error_on_fit /I=1 /D /R 
	FuncFit/N=(suppressUpdateYesNo) FitTOWsafR w_coef phase_on_fit[firstPntToFit,lastPntToFit] /X=wavelength_fit /W=phase_error_on_fit /I=1 /D /R  /M=maskWave /E=epsWave ///TBOX=(256+512)
	//FuncFit/N=(suppressUpdateYesNo)/ODR=2 FitTOWsaf w_coef phase_on_fit[firstPntToFit,lastPntToFit] /X=wavelength_fit /W=phase_error_on_fit /XW=wavelength_err_q /I=1 /D /R
	
	String reswavename = "res_"+nameofwave(phase_on_fit)
	Wave reswave = $reswavename
	ErrorBars $reswavename Y,wave=($nameofwave(phase_error_on_fit),$nameofwave(phase_error_on_fit))
//	if(firstPntToFit!=0)
//		reswave[0,firstPntToFit-1]=NaN
//	endif
//	if(lastPntToFit<numpnts(resWave)-1)
//		reswave[lastPntToFit+1,numpnts(reswave)]=NaN
//	endif
	reswave = (numtype(phase_on_fit)!= 0 || maskWave==0) ? NaN : resWave
	
	String chiStr
	sprintf chiStr "chi sqrd / dof = %.2f", V_chisq/(V_npnts-1)
	print chiStr
	TextBox/C/N=chibox chiStr
	
	
	String fitAlphaName = "afit_"+nameOfWave(phase_on_fit)	//a is for alpha
	
	Wave w_sigma
	Make/O/D/N=200 fitTOWsafphi, fitTOWsafwave, fitTOWsafphifull, fitTOWsafwavefull, $fitAlphaName
	Wave fitAlpha = $fitAlphaName
	
	Wavestats/Q/R=[firstPntToFit,lastPntToFit] wavelength_fit
	SetScale/I x, v_min, v_max, fitTOWsafwave, fitTOWsafphi, fitAlpha
	fitTOWsafwave = x
	
	Duplicate/O w_coef w_coefNoNorm
	w_coefNoNorm[1]= w_coef[1]/abs(w_coef[1])
	print w_coefNoNorm
	FitTOWsafR(w_coefNoNorm, fitAlpha, fitTOWsafwave)
	
	FitTOWsafR(w_coef, fitTOWsafphi, fitTOWsafwave)
	FindLevel/Q fitTOWsafphi, 0
	Variable TOWcen = V_LevelX
	
	w_coef[0]+=w_sigma[0]
	FitTOWsafR(w_coef, fitTOWsafphi, fitTOWsafwave)
	FindLevel/Q fitTOWsafphi, 0	
	Variable TOWpluserr = V_LevelX
	
	w_coef[0]-=2*w_sigma[0]
	FitTOWsafR(w_coef, fitTOWsafphi, fitTOWsafwave)
	FindLevel/Q fitTOWsafphi, 0	
	Variable TOWminuserr = V_LevelX
	
	printf "TOW = %.4f (+%.4f, %.4f) nm\r", TOWcen, TOWpluserr-TOWcen, TOWminuserr-TOWcen
	Variable avgTOWfiterr = (abs(TOWpluserr-TOWcen)+abs(TOWminuserr-TOWcen))/2
	printf "TOW = %.4f ± %.4f nm\r", TOWcen, avgTOWfiterr
	
	
	Wavestats/Q wavelength_fit
	SetScale/I x, v_min, v_max, fitTOWsafwavefull, fitTOWsafphifull
	fitTOWsafwavefull=x
	w_coef[0]+=w_sigma[0]
	FitTOWsafR(w_coef, fitTOWsafphifull, fitTOWsafwavefull)
	
	fitTOWsafphifull = abs(fitTOWsafphifull) > wavemax(phase_on_fit)*2 ? NaN : fitTOWsafphifull
	
	String fitCoefsStr
	sprintf fitCoefsStr "Coefficient values ± one standard deviation\r\tR\t\t=%.4f ± %.4f\r\tphiNorm\t=%.2g ± %.2g", w_coef[0], w_sigma[0], w_coef[1], w_sigma[1]
	TextBox/C/N=fitCoefs fitCoefsStr
	
	String TOWerrstr
	sprintf TOWerrstr	"TOW = %.4f ± %.4f nm", TOWcen, avgTOWfiterr
	TextBox/C/N=TOWbox TOWerrstr

	TOWcen-=0.00056
	Variable Rdopp = CalcR(TOWcen)
	sprintf TOWerrstr	"TOW Dopp. = %.4f ± %.4f nm\rR Dopp. = %.4f ± %.4f", TOWcen, avgTOWfiterr, Rdopp, w_sigma[0]
	TextBox/C/N=TOWboxDop TOWerrstr
	printf "TOW Doppler corrected = %.4f ± %.4f nm\r", TOWcen, avgTOWfiterr
	
	//FitTOWresultsAlgebraic()
End	




Function FitTOWsafR(pw,yw,xw) : FitFunc
	Wave pw, yw, xw	
	
	//assumes xw is in nm
	
	Variable ratio=pw[0]
	Variable phaseNorm = pw[1]
	
	//Parameters from the Saf06 theory for K
	Variable red_p12_D_ns = 4.112*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//reduced matrix elements. unc (5)
	Variable red_p32_D_ns = 5.814*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//(8)
	Variable deltaE_np12 = 0.059165*hartree	//
	Variable deltaE_np32 = 0.059428*hartree
	Variable A_ns = 6.26*AUtoSI
//	Variable A_ns = 6.394*AUtoSI	//from Saf K-Magic personal communications
	//Variable deltaA_ns = 0.33*AUtoSI
	
	
	deltaE_np12 =  c/kd1wavelength*h
	deltaE_np32 =  c/kd2wavelength*h
	
	
	yw=0
	
	yw = A_ns + 1e0*(1/3 * deltaE_np12 *red_p12_D_ns^2/(deltaE_np12^2-((c/xw*1e9)*h)^2) + 1/3*deltaE_np32*ratio*red_p12_D_ns^2/(deltaE_np32^2-((c/xw*1e9)*h)^2))
	
	yw*=phaseNorm
End






Function CalcR(TOW)
	Variable TOW
	
	Variable ratio
	
	//Parameters from the Saf06 theory for K
	Variable red_p12_D_ns = 4.112*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//reduced matrix elements. unc (5)
	Variable red_p32_D_ns = 5.814*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//(8)
	Variable deltaE_np12 = 0.059165*hartree	//
	Variable deltaE_np32 = 0.059428*hartree
	Variable A_ns = 6.26*AUtoSI
//	Variable A_ns = 6.394*AUtoSI	//from Saf K-Magic personal communications
	//Variable deltaA_ns = 0.33*AUtoSI
	
	print/d deltaE_np12
	print/d c/kd1wavelength*h
	
	print/d deltaE_np32
	print/d c/kd2wavelength*h
	
	deltaE_np12 =  c/kd1wavelength*h
	deltaE_np32 =  c/kd2wavelength*h
	
	Variable fTOW = c/TOW*1e9
	
	ratio = -deltaE_np12/deltaE_np32 * (deltaE_np32^2-(h*fTOW)^2)/(deltaE_np12^2-(h*fTOW)^2) - 3*A_ns*(deltaE_np32^2-(h*fTOW)^2)/(deltaE_np32*red_p12_D_ns^2)
	
	print/d "solved from fit eqn first: ", ratio

	Variable fd1 =  c/kd1wavelength
	Variable fd2 = c/kd2wavelength
	Variable ratio2 = -fd1/fd2*(fd2^2-fTOW^2)/(fd1^2-fTOW^2)
	print/d "no A contribution: ", ratio2
	
	
	Variable ratio3 = -(fd2^2-fTOW^2)/fd2 * (fd1/(fd1^2-fTOW^2)+ 3*h*A_ns/red_p12_D_ns^2)
	print/d "paper eqn: ", ratio3
	
	
//	yw = A_ns + 1/3*(deltaE_np12 * red_p12_D_ns^2 / (deltaE_np12^2-((c/xw*1e9)*h)^2) + deltaE_np32*ratio*red_p12_D_ns^2/(deltaE_np32^2-((c/xw*1e9)*h)^2))
	Variable ratio4 = -fd1/fd2 * (fd2^2-(fTOW)^2)/(fd1^2-(fTOW)^2) - 3*A_ns*h*(fd2^2-(fTOW)^2)/(fd2*red_p12_D_ns^2)
	print/d "solved from fit eqn: ", ratio4
	
	return ratio4
//	Variable d1cont = 1/3*deltaE_np12 * red_p12_D_ns^2 / (deltaE_np12^2-((c/TOW*1e9)*h)^2
//	print "d1 contribution = ", d1cont
//	print "A_ns = ", A_ns
//	print "A_ns/d1cont = ", A_ns/d1cont
//	print "Saf au ratio = ", 6.26/32171
End



Function compL()
	print "Saf wavenumber"
	print "Falke fs comb"
	print "NIST ASD"
	print " "
	print/d 1/12985.2*1e7
	print/d c/389286058.716*1e3
	print/d kd1wavelength*1e9
	print " "
	print/d 1/13042.9*1e7
	print/d c/391016170.03*1e3
	print/d kd2wavelength*1e9
End






Function goFitTOWsaf(series_name,all)
	String series_name
	Variable all
	
	Variable oneKtwoK=1
	
	String onektwokStr = ""
	if(OnekTwok==1)
		onektwokStr = "1k2k"
	elseif(OnekTwok == 2)
		onektwokStr = "2k2k"
	endif
	
	String phase_on_avg_name = "phase" + onektwokStr + "_on_avg_" + series_name
	String phase_error_on_avg_name = "phase_error" + onektwokStr + "_on_avg_" + series_name
	String phase_error_on_avg_sdev_name = "phase_error" + onektwokStr + "_on_avg_sdev_" + series_name
	String phase_error_tot_name = "phase_error_tot_" + series_name
	String phase_on_name = "phase" + onektwokStr + "_on_all_" + series_name
	String phase_on_spline_name = "phase" + onektwokStr + "_on_all_spline_" + series_name
	String phase_error_on_name = "phase_error" + onektwokStr + "_on_all_" + series_name
	String wavelength_pad_name = "wavelength_pad_"+series_name
	String wavelength_name = "wavelength_"+series_name
	
	String maskWaveAvgName = "mask_avg_"+series_name
	String maskWaveAllName = "mask_all_"+series_name
	String maskWaveAllSplineName = "mask_all_spline_"+series_name
	
	String phase_bin_name = "PhaseBinned_"+series_name
	String phase_error_bin_name = "PhaseBinnedStdErr_"+series_name
	String wavelength_bin_name = "wavelengthsBinned_"+series_name
	String maskWaveBinnedName = "mask_binned_"+series_name

	String phase_binSS_name = "PhaseBinnedSS_"+series_name
	String phase_error_binSS_name = "PhaseBinnedSSStdErr_"+series_name
	String wavelength_binSS_name = "wavelengthsBinned_"+series_name
	String maskWaveBinnedSSName = "mask_binned_SS_"+series_name
	

	
	if(0)
		phase_error_on_avg_name = phase_error_tot_name
	endif
	
	
	Variable firstPntToFit
	Variable lastPntToFit
	
	// choose to fit the entire data set or the averaged data set
	//Variable all = 0
	if(all==0)
		Wave phase_on_fit = $phase_on_avg_name
		Wave phase_error_on_fit = $phase_error_on_avg_name
		if(0)
			Wave phase_error_on_fit = $phase_error_on_avg_sdev_name
		endif
		Wave wavelength_fit = $wavelength_name
		Wave/Z maskWave = $maskWaveAvgName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==1)
		Wave phase_on_fit = $phase_on_name
		Wave phase_error_on_fit = $phase_error_on_name
		Wave wavelength_fit = $wavelength_pad_name
		Wave/Z maskWave = $maskWaveAllName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==2)
		Wave phase_on_fit = $phase_on_spline_name
		Wave phase_error_on_fit = $phase_error_on_name
		Wave wavelength_fit = $wavelength_pad_name
		Wave/Z maskWave = $maskWaveAllSplineName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==3)
		Wave phase_on_fit = $phase_bin_name
		Wave phase_error_on_fit = $phase_error_bin_name
		Wave wavelength_fit = $wavelength_bin_name
		Wave/Z maskWave = $maskWaveBinnedName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==4)
		Wave phase_on_fit = $phase_binSS_name
		Wave phase_error_on_fit = $phase_error_binSS_name
		Wave wavelength_fit = $wavelength_binSS_name
		Wave/Z maskWave = $maskWaveBinnedSSName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0//315
		lastPntToFit=1000//534
	else
		print "all, binned, ss?"
		abort
	endif
	
	Make/O/D w_coef = {1e-38,1e35}
	
	Variable suppressUpdateYesNo=0
	
	Duplicate/O phase_error_on_fit phase_error_on_fit2
	phase_error_on_fit2*=2
	
	Make/O epsWave = {1e-40,1e30}

	
	//FuncFit/N=(suppressUpdateYesNo) FitTOWsaf w_coef phase_on_fit[firstPntToFit,lastPntToFit] /X=wavelength_fit /W=phase_error_on_fit /I=1 /D /R 
	FuncFit/N=(suppressUpdateYesNo) FitTOWsaf w_coef phase_on_fit[firstPntToFit,lastPntToFit] /X=wavelength_fit /W=phase_error_on_fit /I=1 /D /R  /M=maskWave /E=epsWave ///TBOX=(256+512)
	//FuncFit/N=(suppressUpdateYesNo)/ODR=2 FitTOWsaf w_coef phase_on_fit[firstPntToFit,lastPntToFit] /X=wavelength_fit /W=phase_error_on_fit /XW=wavelength_err_q /I=1 /D /R
	
	String reswavename = "res_"+nameofwave(phase_on_fit)
	Wave reswave = $reswavename
	ErrorBars $reswavename Y,wave=($nameofwave(phase_error_on_fit),$nameofwave(phase_error_on_fit))
//	if(firstPntToFit!=0)
//		reswave[0,firstPntToFit-1]=NaN
//	endif
//	if(lastPntToFit<numpnts(resWave)-1)
//		reswave[lastPntToFit+1,numpnts(reswave)]=NaN
//	endif
	reswave = (numtype(phase_on_fit)!= 0 || maskWave==0) ? NaN : resWave
	
	String chiStr
	sprintf chiStr "chi sqrd / dof = %.2f", V_chisq/(V_npnts-1)
	print chiStr
	TextBox/C/N=chibox chiStr
	
	
	String fitAlphaName = "afit_"+nameOfWave(phase_on_fit)	//a is for alpha
	
	Wave w_sigma
	Make/O/D/N=200 fitTOWsafphi, fitTOWsafwave, fitTOWsafphifull, fitTOWsafwavefull, $fitAlphaName
	Wave fitAlpha = $fitAlphaName
	
	Wavestats/Q/R=[firstPntToFit,lastPntToFit] wavelength_fit
	SetScale/I x, v_min, v_max, fitTOWsafwave, fitTOWsafphi, fitAlpha
	fitTOWsafwave = x
	
	Duplicate/O w_coef w_coefNoNorm
	w_coefNoNorm[1]= w_coef[1]/abs(w_coef[1])
	print w_coefNoNorm
	FitTOWsaf(w_coefNoNorm, fitAlpha, fitTOWsafwave)
	
	FitTOWsaf(w_coef, fitTOWsafphi, fitTOWsafwave)
	FindLevel/Q fitTOWsafphi, 0
	Variable TOWcen = V_LevelX
	
	w_coef[0]+=w_sigma[0]
	FitTOWsaf(w_coef, fitTOWsafphi, fitTOWsafwave)
	FindLevel/Q fitTOWsafphi, 0	
	Variable TOWpluserr = V_LevelX
	
	w_coef[0]-=2*w_sigma[0]
	FitTOWsaf(w_coef, fitTOWsafphi, fitTOWsafwave)
	FindLevel/Q fitTOWsafphi, 0	
	Variable TOWminuserr = V_LevelX
	
	printf "TOW = %.3f (+%.3f, %.3f) nm\r", TOWcen, TOWpluserr-TOWcen, TOWminuserr-TOWcen
	Variable avgTOWfiterr = (abs(TOWpluserr-TOWcen)+abs(TOWminuserr-TOWcen))/2
	printf "TOW = %.3f ± %.3f nm\r", TOWcen, avgTOWfiterr
	
	
	Wavestats wavelength_fit
	SetScale/I x, v_min, v_max, fitTOWsafwavefull, fitTOWsafphifull
	fitTOWsafwavefull=x
	w_coef[0]+=w_sigma[0]
	FitTOWsaf(w_coef, fitTOWsafphifull, fitTOWsafwavefull)
	
	fitTOWsafphifull = abs(fitTOWsafphifull) > wavemax(phase_on_fit)*2 ? NaN : fitTOWsafphifull
	
	String fitCoefsStr
	sprintf fitCoefsStr "Coefficient values ± one standard deviation\r\tA_ns\t\t=%.2g ± %.2g\r\tphiNorm\t=%.2g ± %.2g", w_coef[0], w_sigma[0], w_coef[1], w_sigma[1]
	TextBox/C/N=fitCoefs fitCoefsStr
	
	String TOWerrstr
	sprintf TOWerrstr	"TOW = %.3f ± %.3f nm", TOWcen, avgTOWfiterr
	TextBox/C/N=TOWbox TOWerrstr
	
	FitTOWresultsAlgebraic()
End	





Function FitTOWresultsAlgebraic()
	Wave w_coef, w_sigma
	Variable A_ns = w_coef[0]
	Variable delta_A_ns = w_sigma[0]
	Variable red_p12_D_ns = 4.102*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//reduced matrix elements. unc (5)
	Variable red_p32_D_ns = 5.800*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//(8)
	Variable deltaE_np12 = 0.059165*hartree	//
	Variable deltaE_np32 = 0.059428*hartree
	
	Variable aquad = 3*A_ns
	Variable bquad = -3*A_ns*deltaE_np12^2 - 3*A_ns*deltaE_np32^2 - deltaE_np12*red_p12_D_ns^2 - deltaE_np32*red_p32_D_ns^2
	Variable cquad = 3*A_ns*deltaE_np12^2*deltaE_np32^2 + deltaE_np12*red_p12_D_ns^2*deltaE_np32^2 + deltaE_np32*red_p32_D_ns^2*deltaE_np12^2
	
	Variable fhsqrdquad = (-bquad-sqrt(bquad^2-4*aquad*cquad))/(2*aquad) 
	Variable TOWquad = h*c/sqrt(fhsqrdquad)*1e9
	print TOWquad
	
	A_ns= w_coef[0]+w_sigma[0]
	aquad = 3*A_ns
	bquad = -3*A_ns*deltaE_np12^2 - 3*A_ns*deltaE_np32^2 - deltaE_np12*red_p12_D_ns^2 - deltaE_np32*red_p32_D_ns^2
	cquad = 3*A_ns*deltaE_np12^2*deltaE_np32^2 + deltaE_np12*red_p12_D_ns^2*deltaE_np32^2 + deltaE_np32*red_p32_D_ns^2*deltaE_np12^2
	
	fhsqrdquad = (-bquad-sqrt(bquad^2-4*aquad*cquad))/(2*aquad) 
	Variable TOWquadpluserr = h*c/sqrt(fhsqrdquad)*1e9 - TOWquad
	print TOWquadpluserr
	
	A_ns=w_coef[0]-w_sigma[0]
	aquad = 3*A_ns
	bquad = -3*A_ns*deltaE_np12^2 - 3*A_ns*deltaE_np32^2 - deltaE_np12*red_p12_D_ns^2 - deltaE_np32*red_p32_D_ns^2
	cquad = 3*A_ns*deltaE_np12^2*deltaE_np32^2 + deltaE_np12*red_p12_D_ns^2*deltaE_np32^2 + deltaE_np32*red_p32_D_ns^2*deltaE_np12^2
	
	fhsqrdquad = (-bquad-sqrt(bquad^2-4*aquad*cquad))/(2*aquad) 
	Variable TOWquadminuserr = h*c/sqrt(fhsqrdquad)*1e9 - TOWquad
	print TOWquadminuserr
	
	printf "TOW = %.4f (+%.4f, %.4f) nm\r", TOWquad, TOWquadpluserr, TOWquadminuserr
	Variable avgTOWquadfiterr = (abs(TOWquadpluserr)+abs(TOWquadminuserr))/2
	printf "TOW = %.4f ± %.4f nm\r", TOWquad, avgTOWquadfiterr
End


Function FitTOWsaf(pw,yw,xw) : FitFunc
	Wave pw, yw, xw	
	
	//assumes xw is in nm
	
	Variable A_ns=pw[0]
	Variable phaseNorm = pw[1]
	
	//Parameters from the Saf06 theory for K
	Variable red_p12_D_ns = 4.102*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//reduced matrix elements. unc (5)
	Variable red_p32_D_ns = 5.800*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//(8)
	Variable deltaE_np12 = 0.059165*hartree	//
	Variable deltaE_np32 = 0.059428*hartree
	//Variable A_ns = 6.26*AUtoSI
	//Variable deltaA_ns = 0.33*AUtoSI
	
	yw=0
	
	yw = A_ns + 1e0*(1/3 * deltaE_np12 *red_p12_D_ns^2/(deltaE_np12^2-((c/xw*1e9)*h)^2) + 1/3*deltaE_np32*red_p32_D_ns^2/(deltaE_np32^2-((c/xw*1e9)*h)^2))
	
	yw*=phaseNorm
End



Function alpha0()
	//Wave pw, yw, xw	
	
	//assumes xw is in nm
	
	Variable A_ns=-5.51e-39
	//Variable phaseNorm = pw[1]
	
	//Parameters from the Saf06 theory for K
	Variable red_p12_D_ns = 4.102*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//reduced matrix elements. unc (5)
	Variable red_p32_D_ns = 5.800*sqrt(hartree*AUtoSI)//(a0tom)^(3/2)//*hbar	//(8)
	Variable deltaE_np12 = 0.059165*hartree	//
	Variable deltaE_np32 = 0.059428*hartree
	//Variable A_ns = 6.26*AUtoSI
	//Variable deltaA_ns = 0.33*AUtoSI
	
	//yw=0
	
	Variable alpha0 = A_ns + 1/3 *red_p12_D_ns^2/deltaE_np12 + 1/3*red_p32_D_ns^2/deltaE_np32
	
	print alpha0
	print alpha0/AUtoSI
	//yw*=phaseNorm
End



Function alpha02()
	//Wave pw, yw, xw	
	
	//assumes xw is in nm
	
	Variable A_ns=-5.51e-39/AUtoSI
	print a_ns
	//Variable phaseNorm = pw[1]
	
	//Parameters from the Saf06 theory for K
	Variable red_p12_D_ns = 4.102//(a0tom)^(3/2)//*hbar	//reduced matrix elements. unc (5)
	Variable red_p32_D_ns = 5.800//(a0tom)^(3/2)//*hbar	//(8)
	Variable deltaE_np12 = 0.059165	//
	Variable deltaE_np32 = 0.059428
	//Variable A_ns = 6.26*AUtoSI
	//Variable deltaA_ns = 0.33*AUtoSI
	
	//yw=0
	
	Variable alpha0 = 1/3 *red_p12_D_ns^2/deltaE_np12 + 1/3*red_p32_D_ns^2/deltaE_np32
	
	print alpha0
	//print alpha0/AUtoSI
	//yw*=phaseNorm
End





Function TOWwaveErr(series_name)
	String series_name
	
	String phase_error_wave_name = "phase_error1k2k_on_avg_"+series_name
	
	wave phase_error_wave = $phase_error_wave_name
	wave w_coef
	
	String wavelength_name = "wavelength_"+series_name
	Wave wavelength_fit = $wavelength_name
	
	Make/O/D/N=1000 fitTOWsafphifull, fitTOWsafwavefull
	Wavestats/Q wavelength_fit
	SetScale/I x, v_min-.1, v_max+.1, fitTOWsafwaveFull, fitTOWsafphiFull
	fitTOWsafwaveFull = x
		
	FitTOWsaf(w_coef, fitTOWsafphifull, fitTOWsafwavefull)
	
	Duplicate/O fitTOWsafPhiFull phiFitDer, phiWaveErrFull
	Differentiate phiFitDer
	
	phiWaveErrFull = phiFitDer*(fitTOWsafwaveFull)^2/(c*1e9)*1.6e9
	
	String phase_error_wavelength_name = "phase_error_wavelength_"+series_name
	String phase_error_tot_name = "phase_error_tot_"+series_name
	Duplicate/O wavelength_fit $phase_error_wavelength_name, $phase_error_tot_name
	Wave phase_error_wavelength = $phase_error_wavelength_name
	Wave phase_error_tot = $phase_error_tot_name
	
	
	phase_error_wavelength = phiWaveErrFull(wavelength_fit[p])
	phase_error_tot = sqrt(phase_error_wave^2 + phase_error_wavelength^2)
End




Function FilterTOWfit(initialize, StdDevFilter)
	Variable initialize, StdDevFilter
	
	String traces = TraceNameList("", ";", 1)							// Get all traces on the graph
	String resTraceName = StringFromList(0,GrepList(traces, "Res"))	// Find the name of the trace that starts with "Res"
	
	Wave resWave = $resTraceName
	
	String series_name
	SplitString /E="\\S*_(\\S*)$" resTraceName, series_name

	String maskWaveName
	
	//Figure out which mask wave to reference
	if(GrepString(resTraceName,"all"))
		if(GrepString(resTraceName,"spline"))
			maskWaveName = "mask_all_spline_"+series_name
		else
			maskWaveName = "mask_all_"+series_name
		endif
	else
		maskWaveName = "mask_"+series_name
	endif
	
	Wave maskWave = $maskWaveName
	
	//Reinitialize mask wave or perform filter. Return NaN or number of filtered points.
	if(initialize==1)
		maskWave = 1
		Return NaN
	else
		Wavestats/Q resWave
		maskWave = (abs(resWave) > StdDevFilter*v_sdev || maskWave==0) ? 0 : 1
		Variable totalFilteredPoints = numpnts(maskWave)-sum(maskWave)
		//print "Total filtered points = ", totalFilteredPoints
		Return totalFilteredPoints
	endif
End




//////////////////////


Function goFitTOWpow(series_name)
	String series_name
	Variable oneKtwoK=1
	
	String onektwokStr = ""
	if(OnekTwok==1)
		onektwokStr = "1k2k"
	elseif(OnekTwok == 2)
		onektwokStr = "2k2k"
	endif
	
	String phase_on_name = "phase" + onektwokStr + "_on_" + series_name
	String phase_error_on_name = "phase_error" + onektwokStr + "_on_" + series_name
	String wavelength_pad_name = "wavelength_pad_"+series_name
	Wave phase_on = $phase_on_name
	Wave phase_error_on = $phase_error_on_name
	Wave wavelength_pad = $wavelength_pad_name
		
	Variable numOrders = 3
	Make/O/D/N=(numOrders+2) w_coef = 0.01
	w_coef[0]=769
	w_coef[1]=0.1
	
	Variable suppressUpdateYesNo=0
	
	Variable firstPntToFit=60
	Variable lastPntToFit=289
	
	FuncFit/N=(suppressUpdateYesNo)/TBOX=(256+512) FitTOWpow w_coef phase_on[firstPntToFit,lastPntToFit] /X=wavelength_pad /W=phase_error_on /I=1 /D /R 
	
	
	String reswavename = "res_"+phase_on_name	
	Wave reswave = $reswavename
	ErrorBars $reswavename Y,wave=($phase_error_on_name,$phase_error_on_name)
	reswave[0,firstPntToFit-1]=NaN
	reswave[lastPntToFit+1,numpnts(reswave)]=NaN
End




Function FitTOWpow(pw,yw,xw) : FitFunc
	Wave pw, yw, xw
	
	Variable TOW=pw[0]
	Variable phaseNorm = pw[1]
	
	yw=0
	
	Variable numOrders = numpnts(pw)-2
	Variable n
	For(n=1; n<=numOrders;n+=1)
		yw+=pw[1+n]*(xw-TOW)^n
	EndFor
	
	yw*=phaseNorm
End







Function goFitTOWlin(series_name,all)
	String series_name
	Variable all
	
	Variable oneKtwoK=1
	
	String onektwokStr = ""
	if(OnekTwok==1)
		onektwokStr = "1k2k"
	elseif(OnekTwok == 2)
		onektwokStr = "2k2k"
	endif
	
	String phase_on_avg_name = "phase" + onektwokStr + "_on_avg_" + series_name
	String phase_error_on_avg_name = "phase_error" + onektwokStr + "_on_avg_" + series_name
	String phase_error_on_avg_sdev_name = "phase_error" + onektwokStr + "_on_avg_sdev_" + series_name
	String phase_error_tot_name = "phase_error_tot_" + series_name
	String phase_on_name = "phase" + onektwokStr + "_on_all_" + series_name
	String phase_on_spline_name = "phase" + onektwokStr + "_on_all_spline_" + series_name
	String phase_error_on_name = "phase_error" + onektwokStr + "_on_all_" + series_name
	String wavelength_pad_name = "wavelength_pad_"+series_name
	String wavelength_name = "wavelength_"+series_name
	
	String maskWaveAvgName = "mask_avg_"+series_name
	String maskWaveAllName = "mask_all_"+series_name
	String maskWaveAllSplineName = "mask_all_spline_"+series_name
	
	String phase_bin_name = "PhaseBinned_"+series_name
	String phase_error_bin_name = "PhaseBinnedStdErr_"+series_name
	String wavelength_bin_name = "wavelengthsBinned_"+series_name
	String maskWaveBinnedName = "mask_binned_"+series_name

	String phase_binSS_name = "PhaseBinnedSS_"+series_name
	String phase_error_binSS_name = "PhaseBinnedSSStdErr_"+series_name
	String wavelength_binSS_name = "wavelengthsBinned_"+series_name
	String maskWaveBinnedSSName = "mask_binned_SS_"+series_name
	

	
	if(0)
		phase_error_on_avg_name = phase_error_tot_name
	endif
	
	
	Variable firstPntToFit
	Variable lastPntToFit
	
	// choose to fit the entire data set or the averaged data set
	//Variable all = 0
	if(all==0)
		Wave phase_on_fit = $phase_on_avg_name
		Wave phase_error_on_fit = $phase_error_on_avg_name
		if(0)
			Wave phase_error_on_fit = $phase_error_on_avg_sdev_name
		endif
		Wave wavelength_fit = $wavelength_name
		Wave/Z maskWave = $maskWaveAvgName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==1)
		Wave phase_on_fit = $phase_on_name
		Wave phase_error_on_fit = $phase_error_on_name
		Wave wavelength_fit = $wavelength_pad_name
		Wave/Z maskWave = $maskWaveAllName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==2)
		Wave phase_on_fit = $phase_on_spline_name
		Wave phase_error_on_fit = $phase_error_on_name
		Wave wavelength_fit = $wavelength_pad_name
		Wave/Z maskWave = $maskWaveAllSplineName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==3)
		Wave phase_on_fit = $phase_bin_name
		Wave phase_error_on_fit = $phase_error_bin_name
		Wave wavelength_fit = $wavelength_bin_name
		Wave/Z maskWave = $maskWaveBinnedName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0
		lastPntToFit=1000
	elseif(all==4)
		Wave phase_on_fit = $phase_binSS_name
		Wave phase_error_on_fit = $phase_error_binSS_name
		Wave wavelength_fit = $wavelength_binSS_name
		Wave/Z maskWave = $maskWaveBinnedSSName
		if(waveexists(maskWave))
		else
			Duplicate/O phase_on_fit maskWave
			maskWave = 1
		endif
		firstPntToFit=0//315
		lastPntToFit=1000//534
	else
		print "all, binned, ss?"
		abort
	endif
	
	Make/O/D w_coef = {768.97,1}
	
	Variable suppressUpdateYesNo=0
	
	Duplicate/O phase_error_on_fit phase_error_on_fit2
	phase_error_on_fit2*=2
	
//	Make/O epsWave = {1e-4,1e30}

	
	//FuncFit/N=(suppressUpdateYesNo) FitTOWsaf w_coef phase_on_fit[firstPntToFit,lastPntToFit] /X=wavelength_fit /W=phase_error_on_fit /I=1 /D /R 
	FuncFit/N=(suppressUpdateYesNo) FitTOWlin w_coef phase_on_fit[firstPntToFit,lastPntToFit] /X=wavelength_fit /W=phase_error_on_fit /I=1 /D /R  /M=maskWave ///E=epsWave ///TBOX=(256+512)
	//FuncFit/N=(suppressUpdateYesNo)/ODR=2 FitTOWsaf w_coef phase_on_fit[firstPntToFit,lastPntToFit] /X=wavelength_fit /W=phase_error_on_fit /XW=wavelength_err_q /I=1 /D /R
	
	String reswavename = "res_"+nameofwave(phase_on_fit)
	Wave reswave = $reswavename
	ErrorBars $reswavename Y,wave=($nameofwave(phase_error_on_fit),$nameofwave(phase_error_on_fit))
//	if(firstPntToFit!=0)
//		reswave[0,firstPntToFit-1]=NaN
//	endif
//	if(lastPntToFit<numpnts(resWave)-1)
//		reswave[lastPntToFit+1,numpnts(reswave)]=NaN
//	endif
	reswave = (numtype(phase_on_fit)!= 0 || maskWave==0) ? NaN : resWave
	
	String chiStr
	sprintf chiStr "chi sqrd / dof = %.2f", V_chisq/(V_npnts-1)
	print chiStr
	TextBox/C/N=chibox chiStr
	
	Wave w_sigma
	
	String fitCoefsStr
	sprintf fitCoefsStr "Coefficient values ± one standard deviation\r\tTOW\t\t\t=%.4f ± %.4f\r\trad/Watt/NM\t=%.2f ± %.2f", w_coef[0], w_sigma[0], w_coef[1], w_sigma[1]
	TextBox/C/N=fitCoefs fitCoefsStr
	
	String TOWlinErrstr
	sprintf TOWlinErrstr	"TOWlin = %.4f ± %.4f nm", w_coef[0], w_sigma[0]
	TextBox/C/N=TOWboxLin TOWlinErrstr
	
	printf "TOW = %.4f ± %.4f nm\r", w_coef[0], w_sigma[0]
End	





Function FitTOWlin(w,wavelength) : FitFunc
	Wave w
	Variable wavelength

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(wavelength) = radPerWattPerNM*(wavelength-TOW)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ wavelength
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = TOW
	//CurveFitDialog/ w[1] = radPerWattPerNM

	return w[1]*(wavelength-w[0])
End





Function FitTOWcube(w,wavelength) : FitFunc
	Wave w
	Variable wavelength

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(wavelength) = phaselin*(wavelength-TOW)+phasecube*(wavelength-TOW)^3
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ wavelength
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = TOW
	//CurveFitDialog/ w[1] = phaselin
	//CurveFitDialog/ w[2] = phasecube

	return w[1]*(wavelength-w[0])+w[2]*(wavelength-w[0])^3
End

Function FitTOWd1d2(w,wavelength) : FitFunc
	Wave w
	Variable wavelength

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(wavelength) = c*phinorm*(3/8*(1/wavelength-1/770.1083)+5/8*(1/wavelength-1/770.1083))+offset
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ wavelength
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = phinorm
	//CurveFitDialog/ w[1] = offset

	return c*w[0]*(3/8*(1/wavelength-1/770.1083)+5/8*(1/wavelength-1/770.1083))+w[1]
End

Function FitTOW0linquadcube(w,wavelength) : FitFunc
	Wave w
	Variable wavelength

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(wavelength) = phase0+phaselin*(wavelength-TOW)+phasequad*(wavelength-TOW)^2+phasecube*(wavelength-TOW)^3
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ wavelength
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = TOW
	//CurveFitDialog/ w[1] = phase0
	//CurveFitDialog/ w[2] = phaselin
	//CurveFitDialog/ w[3] = phasequad
	//CurveFitDialog/ w[4] = phasecube

	return w[1]+w[2]*(wavelength-w[0])+w[3]*(wavelength-w[0])^2+w[4]*(wavelength-w[0])^3
End
