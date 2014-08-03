#pragma rtGlobals=1		// Use modern global access method.

#include ":WindowNamer"
#include <Remove Points>
#include ":RemoveNaNs"

// 		INSTRUCTIONS
// Set directory with fringe data
// set series name, laser calib name, start index and stop index
// if you want to analyze every file in a series you can set the stop index to 0
// set the time of flight parameter
// run FringeFit()
// if you want a copy of your data, run SaveFringeData(appended_name, yesno).  yesno should be set to 0 unless you want to enable overwrites, in which case it should be 1




// FringeFitTrim() is essentially a more advanced version of the original FringeFit()
// It bins the data, then finds the trimmed mean and std err of each bin
// filter?
// center the laser calib and voltage on 0?


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function FringeFitTrim()

	// Choose your directory (uncomment the one you want, recomment the previously used one)

	//string/G directory_name = "C:Users:Alex:Desktop:DAQ:Data:2010:100929:fringe data:"
	String/G directory_name = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:120502:fringe data:"
	//String/G directory_name = "Macintosh HD:Users:holmgren:Documents:School:UA:Cronin labs:Dropbox:111010:fringe data:"

	// Choose your data series and laser calibration file
	string/G series_name = "g"
	string laser_calib_series = "f"
	string laser_calib_number = "1"
	string laser_calib_name = "laser_calib_"+laser_calib_series+laser_calib_number+".txt"	
	string appended_name = ""
	//sprintf laser_calib_name, ("laser_calib_"+series_name+"1.txt")
	
	// Choose the data file range within the specified series
	variable/G start_index=1
	variable/G stop_index=1		//Specify stop_index = 0 if you want all the data for this series analyzed
	
	Variable overwriteYesNo = 1
	Variable display1k1kYesNo = 0
	Variable display2k2kYesNo
	
	Variable suppressUpdateYesNo =1
	
	Variable/G fringePlotMin = 0
	Variable/G fringePlotMax = 400
	
	if(stop_index==0)
		stop_index = GetTotalFileNumber(directory_name, series_name)
	endif
	
	Variable DisplaySmallGraphs = 0

	string fname,ind							// file name and index strings
	variable/G Lmax, Lmin, Lavg, Lamp 			// for laser calibration
	
	variable tof = 0								// tof = time of flight in integer ms. compensates for delay between when an atom is incident on the grating and when it is counted. 
	
	Variable advancedFilterYesNo = 0
	Variable filterStdDevYesNo = 0
	
	Variable resolution = 1		// errors are not reported correctly if resample != 1
	
	Variable numStdDev = 3 				// number of standard deviations a point can be from the best fit line before it is marked for removal
	
	// prepare laser calibration 
	load(directory_name+laser_calib_name) 					// load laser calibration file
	
	duplicate/o voltage, laser_calib_voltage
	wavestats/q laser_calib_voltage  							// get wavestats and then calculate max, min, avg, amp.
	Lmax = V_max ;  Lmin = V_min ;  Lavg = (Lmax+Lmin)/2  ;   Lamp = (Lmax-Lmin)/2//; print "Lavg = ", Lavg; print "Lamp  = ", Lamp

	// make waves in which to store the contrast, phase, etc of each data file
	make/o/d/n=(stop_index-start_index+1) contrast, contrast_error, phase, phase_error,  avg_counts, avg_counts_error,  file_index,  simple_avg, chisqrdDOF, avg_power, avg_powerStdDev
	Duplicate/o contrast contrast1k2k, contrast_error1k2k, phase1k2k, phase_error1k2k, avg_counts1k2k, avg_counts_error1k2k, phase2k2k, phase_error2k2k, contrast2k2k, contrast_error2k2k, chisqrdDOF2k
	
	Make/D/N=4/O w_coef, w_sigma				// w[0]= avg_counts, w[1]=phase, w[2]=contrast, w[3]=grating period
	W_coef[0] = {20, 0, .2, 100}		// initial guesses for the fit parameters
	
	Make/D/N=6/O w_coef2k//, w_sigma2k
	w_coef2k = {20,0,.25,100,0,.05}
	
	// analyze each data file
	variable index
	String indexString
	
	variable displayYesNo = 1
	variable displayn=1					// display fringe plots? 1 = yes, 0 = no
	
	For(index=start_index; index<=stop_index; index+=1)		// repeat loop through the rest of the specified data files
	
		indexString = num2str(index)
		sprintf fname,(directory_name+series_name+"%g.txt"),index
									
		load(fname)										// load a single data file
		
		ReindexToTOF(tof)										// shifts voltage and counts into alignment
		if(advancedFilterYesNo)
			AdvancedFilter()   											// kill noise spikes from atom counts (and delete laser corresponding point too)
		Endif
		//FilterRawWithStdDevCut(counts,2)
		fringeTrim(voltagetof,counts,Lavg,Lamp, resolution)	// transform counts vs time into counts vs postion
		
		if(displayn == 1 && displayYesNo)								
			display_fringes2(indexString)								// display the fringe plot
			displayn=0											// if displayn=0 no additional fringe plots will be displayed. The first will be updated though.
		endif
		
		FitCountsVsGratingPosition(index,start_index,stop_index, 1, suppressUpdateYesNo)				// fit the atom fringes (avg counts vs grating position) with a sine function
		
		if(filterStdDevYesNo)											// filter data by comparing to best-fit ?
			FilterCountsWithBestFit(numStdDev)							// filter
			FitCountsVsGratingPosition(index,start_index,stop_index, 1, suppressUpdateYesNo)	// refit
		endif

		wavestats/q counts
		simple_avg[index-start_index] = v_avg
		
		Wave/Z power
		if(waveexists(power))
			wavestats/q/r=[500,4800] power
			avg_power[index-start_index] = v_avg
			avg_powerStdDev[index-start_index] = v_sdev
		endif
		

		if(displayYesNo)
			TextBox/C/N=text0/X=0/Y=0 indexString						// update textbox
		endif
	endfor								

	file_index = x+start_index									// populate the file_index wave
	
	if(DisplaySmallGraphs)
		display_avg_counts_results()								// display results
		display_phase_results()	
		display_contrast_results()	
	endif
	
	//SaveFringeData(series_name+appended_name, overwriteYesNo,display1k1kYesNo)
	//SaveFringeData2k(series_name+appended_name, overwriteYesNo,display2k2kYesNo)
	
	//GetHeader(directory_name, series_name, start_index, stop_index)
end
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function FringeFitNew()

	// Choose your directory (uncomment the one you want, recomment the previously used one)

	//string/G directory_name = "C:Users:Alex:Desktop:DAQ:Data:2010:100929:fringe data:"
	String/G directory_name = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:110930:fringe data:"
	//String/G directory_name = "Macintosh HD:Users:holmgren:Documents:School:UA:Cronin labs:Dropbox:111010:fringe data:"

	// Choose your data series and laser calibration file
	string/G series_name = "c"
	string laser_calib_series = "a"
	string laser_calib_number = "1"
	string laser_calib_name = "laser_calib_"+laser_calib_series+laser_calib_number+".txt"	
	string appended_name = ""
	//sprintf laser_calib_name, ("laser_calib_"+series_name+"1.txt")
	
	// Choose the data file range within the specified series
	variable/G start_index=14
	variable/G stop_index=14		//Specify stop_index = 0 if you want all the data for this series analyzed
	
	Variable overwriteYesNo = 1
	Variable display1k1kYesNo = 0
	Variable display2k2kYesNo
	
	Variable suppressUpdateYesNo =1
	
	Variable/G fringePlotMin = 0
	Variable/G fringePlotMax = 400
	
	if(stop_index==0)
		stop_index = GetTotalFileNumber(directory_name, series_name)
	endif
	
	Variable DisplaySmallGraphs = 0

	string fname,ind							// file name and index strings
	variable Lmax, Lmin, Lavg, Lamp 			// for laser calibration
	
	variable tof = 0								// tof = time of flight in integer ms. compensates for delay between when an atom is incident on the grating and when it is counted. 
	
	Variable filterYesNo = 1
	
	Variable resampleRate = 1		// errors are not reported correctly if resample != 1
	
	Variable numStdDev = 3 				// number of standard deviations a point can be from the best fit line before it is marked for removal
	
	// prepare laser calibration 
	load(directory_name+laser_calib_name) 					// load laser calibration file
	
	duplicate/o voltage, laser_calib_voltage
	wavestats/q laser_calib_voltage  							// get wavestats and then calculate max, min, avg, amp.
	Lmax = V_max ;  Lmin = V_min ;  Lavg = (Lmax+Lmin)/2  ;   Lamp = (Lmax-Lmin)/2//; print "Lavg = ", Lavg; print "Lamp  = ", Lamp

	// make waves in which to store the contrast, phase, etc of each data file
	make/o/d/n=(stop_index-start_index+1) contrast, contrast_error,   phase, phase_error,   avg_counts, avg_counts_error,  file_index,  simple_avg, chisqrdDOF, avg_power, avg_powerStdDev
	Duplicate/o contrast contrast1k2k, contrast_error1k2k, phase1k2k, phase_error1k2k, avg_counts1k2k, avg_counts_error1k2k, phase2k2k, phase_error2k2k, contrast2k2k, contrast_error2k2k, chisqrdDOF2k
	
	Make/D/N=4/O w_coef, w_sigma				// w[0]= avg_counts, w[1]=phase, w[2]=contrast, w[3]=grating period
	W_coef[0] = {20, 0, .2, 100}		// initial guesses for the fit parameters
	
	Make/D/N=6/O w_coef2k//, w_sigma2k
	w_coef2k = {20,0,.25,100,0,.05}
	
	// analyze each data file
	variable index
	String indexString
	
	variable displayYesNo =1
	variable displayn=1					// display fringe plots? 1 = yes, 0 = no
	
	For(index=start_index; index<=stop_index; index+=1)		// repeat loop through the rest of the specified data files
	
		indexString = num2str(index)
		sprintf fname,(directory_name+series_name+"%g.txt"),index
									
		load(fname)										// load a single data file
		
		ReindexToTOF(tof)										// shifts voltage and counts into alignment
		//filter()   											// kill noise spikes from atom counts (and delete laser corresponding point too)
		FilterRawWithStdDevCut(counts,2)
		fringe2(voltagetof,counts,Lmax,Lmin,Lavg,Lamp, resampleRate)	// transform counts vs time into counts vs postion
		
		if(displayn == 1 && displayYesNo)								
			display_fringes2(indexString)								// display the fringe plot
			displayn=0											// if displayn=0 no additional fringe plots will be displayed. The first will be updated though.
		endif
		
		FitCountsVsGratingPosition(index,start_index,stop_index, 0, suppressUpdateYesNo)				// fit the atom fringes (avg counts vs grating position) with a sine function
		
		if(filterYesNo)											// filter data by comparing to best-fit ?
			FilterCountsWithBestFit(numStdDev)							// filter
			FitCountsVsGratingPosition(index,start_index,stop_index, 1, suppressUpdateYesNo)	// refit
		endif

		wavestats/q counts
		simple_avg[index-start_index] = v_avg
		
		Wave/Z power
		if(waveexists(power))
			wavestats/q/r=[500,4800] power
			avg_power[index-start_index] = v_avg
			avg_powerStdDev[index-start_index] = v_sdev
		endif
		

		if(displayYesNo)
			TextBox/C/N=text0/X=0/Y=0 indexString						// update textbox
		endif
	endfor								

	file_index = x+start_index									// populate the file_index wave
	
	if(DisplaySmallGraphs)
		display_avg_counts_results()								// display results
		display_phase_results()	
		display_contrast_results()	
	endif
	
	SaveFringeData(series_name+appended_name, overwriteYesNo,display1k1kYesNo)
	SaveFringeData2k(series_name+appended_name, overwriteYesNo,display2k2kYesNo)
	
	GetHeader(directory_name, series_name, start_index, stop_index)
end
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








Function FringeFit()

	// Choose your directory (uncomment the one you want, recomment the previously used one)

	//string/G directory_name = "C:Users:Alex:Desktop:DAQ:Data:2010:100929:fringe data:"
	String/G directory_name = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:110930:fringe data:"

	// Choose your data series and laser calibration file
	string/G series_name = "c"
	string laser_calib_name = "laser_calib_a1.txt"	
	//sprintf laser_calib_name, ("laser_calib_"+series_name+"1.txt")
	
	// Choose the data file range within the specified series
	variable/G start_index=14
	variable/G stop_index=14		//Specify stop_index = 0 if you want all the data for this series analyzed
	//doesn't yet work for series "l" or for "aa" etc.
	
	if(stop_index==0)
		stop_index = GetTotalFileNumber(directory_name, series_name)
	endif

	variable displayn=1						// display fringe plots? 1 = yes, 0 = no
	string fname,ind							// file name and index strings
	variable Lmax, Lmin, Lavg, Lamp 			// for laser calibration
	variable bins = 50   						// the number of grating position bins the counts will be placed in
	
	variable tof = 0								// tof = time of flight in integer ms. compensates for delay between when an atom is incident on the grating and when it is counted. 
	
	
	// prepare laser calibration 
	load(directory_name+laser_calib_name) 					// load laser calibration file
	
	duplicate/o voltage, laser_calib_voltage
	wavestats/q laser_calib_voltage  							// get wavestats and then calculate max, min, avg, amp.
	Lmax = V_max ;  Lmin = V_min ;  Lavg = (Lmax+Lmin)/2  ;   Lamp = (Lmax-Lmin)/2

	// make waves in which to store the contrast, phase, etc of each data file
	make/o/d/n=(stop_index-start_index+1) contrast, contrast_error,   phase, phase_error,   avg_counts, avg_counts_error,  file_index,  simple_avg
	
	// analyze each data file
	variable index=start_index
	do
		sprintf ind,"%g",index
		sprintf fname,(directory_name+series_name+"%g.txt"),index
									
		load(fname)										// load a single data file
		
		ReindexToTOF(tof)										// shifts voltage and counts into alignment
		//filter()   											// kill noise spikes from atom counts (and delete laser corresponding point too)
		fringe(voltagetof,counts,Lmax,Lmin,Lavg,Lamp,bins)	// transform counts vs time into counts vs postion
		sinfit(index,start_index,stop_index,bins)				// fit the atom fringes (avg counts vs grating position) with a sine function

		wavestats/q counts
		simple_avg[index] = v_avg
		
		if(displayn)								
			display_fringes(ind)								// display the fringe plot
			displayn=0											// if displayn=0 no additional fringe plots will be displayed. The first will be updated though.
		endif

		TextBox/C/N=text0/X=0/Y=0 ind
		index += 1
	while(index<=stop_index)									// repeat loop through the rest of the specified data files

	file_index = x+start_index									// populate the file_index wave
	
	display_avg_counts_results()								// display results
	display_phase_results()	
	display_contrast_results()	
	
	SaveFringeData(series_name,1,0)
end





Function modIt(wavetomod, phaseShift)
	Wave wavetoMod
	Variable phaseShift
	
	waveToMod+=phaseShift
	waveToMod=mod(waveToMod,2*pi)
End



Function FilterRawWithStdDevCut(wavein, numStdDev)
	Wave wavein
	Variable numStdDev
	
	Wavestats/Q wavein
	
	Variable avg = V_avg
	Variable stdDev = V_sdev
	
	wavein = wavein > (v_avg+numStdDev*stdDev) ? NaN : wavein
End



Function FilterCountsWithBestFit(numStdDev)
	Variable numStdDev
	
	Wave counts, counts_error, gratingPosition, w_coef
	
	Duplicate/O counts bestFitCounts, countsResiduals, maskHi, maskLo
	
	bestFitCounts = abs(w_coef[0])*(1 + abs(w_coef[2])*sin(2*pi/abs(w_coef[3])*gratingPosition + w_coef[1]))
	
	countsResiduals = bestFitCounts - counts
	
	// maskHi(Lo) = 1 if point is ok, = 0 if point should be removed
	maskHi = counts < bestFitCounts + counts_error*numStdDev /// 1.5
	maskLo = counts > bestFitCounts - counts_error*numStdDev
	
	counts = counts / maskHi / maskLo
	
	RemoveNaNsInfsXYZ(counts, counts_error, gratingPosition)
	
	Sort gratingPosition counts, counts_error, gratingPosition
	
	Variable minGratingPosition = wavemin(gratingPosition)
	Variable maxGratingPosition = wavemax(gratingPosition)
	
	Variable cutoffDistance = 0 //nm
	
	FindLevel/Q/P gratingPosition maxGratingPosition-cutoffDistance
	Variable startPoint = round(V_LevelX)
	DeletePoints startPoint, numpnts(gratingPosition)-startPoint, gratingPosition, counts, counts_error
	
	FindLevel/Q/P gratingPosition minGratingPosition+cutoffDistance
	DeletePoints 0, round(V_LevelX), gratingPosition, counts, counts_error	
	
	
	//countsResiduals = bestFitCounts - counts
End





Function GetTotalFileNumber(directory_name, series_name)
	String directory_name, series_name

	NewPath/O/Q myPath directory_name						// create a path for the desired directory
	String fileList = IndexedFile(myPath, -1, ".txt")			// create a list containing all .txt files in the specified symbolic path.
	fileList = GrepList(fileList, "^"+series_name+"\\d*.txt", 0, ";")	// search the list for all files of the form "start of string + series_name + digit"
	Variable numFiles = ItemsInList(fileList, ";")				// count the number of items in that list

	return numFiles		
End





Function GetHeader(directory_name, series_name, start_index, stop_index)
	String directory_name, series_name
	Variable start_index, stop_index
	
	//make some wavenames, waves, and wave assignments to save the header info in
	String timeWaveName = "time_"+series_name
	String timeDiffWaveName = "timeDiff_" + series_name
	String freqWaveName = "freq_"+series_name
	String HVpWaveName = "HVp_"+series_name
	String HVmWaveName = "HVm_"+series_name
	String PosWaveName = "Pos_"+series_name

	Make/O/D/N=(stop_index-start_index+1) $timeWaveName, $timeDiffWaveName, $freqWaveName, $HVpWaveName, $HVmWaveName, $PosWaveName
	
	Wave timeWave = $timeWaveName
	Wave timeDiffWave = $timeDiffWaveName
	Wave freqWave = $freqWaveName
	Wave HVpWave = $HVpWaveName 
	Wave HVmWave = $HVmWaveName
	Wave PosWave = $PosWaveName

	// make variables for header info from each file
	Variable hours, minutes, seconds, timetotal
	Variable frequency
	Variable hvp, hvm
	Variable pos

	NewPath/O/Q myPath directory_name
	
	String greped = ""
	String header = ""
	String thisFile
	
	// Search each file for the header information
	Variable i
	for(i = start_index; i<=stop_index; i+=1)
		thisFile = series_name+num2str(i)+".txt"
	
		Grep/LIST/Q/P=myPath/E="#" thisFile		// pulls out all lines starting with a # and puts each line into a ; separated string
		header = S_value

		//find time
		greped = GrepList(header, "#\\d*:\\d*:\\d"); //print greped		// find the part of the list item with the time stamp
		sscanf greped, "#%f%*[:]%f%*[:]%f", hours, minutes, seconds		// from that list item, find three numbers
		timetotal = hours*3600+minutes*60+seconds
		timeWave[i-start_index]=timetotal				// seconds since midnight
		
		//find chopper frequency
		greped = GrepList(header, "# Chopper")		// find the list item with the chopper frequency
		sscanf greped, "# Chopper frequency %f", frequency
		freqWave[i-start_index]= frequency
		
		//find postive pillar voltage
		greped = GrepList(header, "# HV\+")
		sscanf greped, "# HV+ %f", hvp
		HVpWave[i-start_index]=hvp
        
        	//find negative pillar voltage
        	greped = GrepList(header, "# HV-")
		sscanf greped, "# HV- %f", hvm
		HVmWave[i-start_index]=hvm
		
		//find motor position
		greped = GrepList(header, "# Motor position")
		sscanf greped, "# Motor position %f", pos
		PosWave[i-start_index]=pos
		
		if(0)	// print found variables for debugging
			print header
			printf "%f:%f:%f	%f\r", hours, minutes, seconds, timetotal
			print "freq: ", frequency
			print "hvp: ", hvp
			print "hvm: ", hvm
			print "pos: ", pos		
		endif
	endfor
	
	// calculate a wave with seconds since first file was written
	Variable t0 = timeWave[0]			// seconds since midnight for first file
	timeDiffWave =  timeWave - t0
	
//	Duplicate/O timeWave times			// generic waves without the series name appended
//	Duplicate/O freqWave freq
	
	print nameofwave(timeWave)
	print nameofwave(timeDiffWave)
	print nameofwave(freqWave)
	print nameofwave(hvpwave)
	print nameofwave(hvmwave)
	print nameofwave(poswave)
End







Function GetFrequencies(directory_name, series_name, start_index, stop_index)
	String directory_name, series_name
	Variable start_index, stop_index
	
	String freqWaveName = "freq_"+series_name
	Make/O/N=(stop_index-start_index+1) $freqWaveName
	Wave freqWave = $freqWaveName

	NewPath/O/Q myPath directory_name

	String greped = ""
	Variable f
	String thisFile
	
	Variable i
	for(i = start_index; i<=stop_index; i+=1)
		thisFile = series_name+num2str(i)+".txt"
	
		Grep/LIST/Q/P=myPath/E="# Chopper frequency" thisFile
		greped = S_value
		//print greped
		sscanf greped, "# Chopper frequency %f", f
		//print f
		freqWave[i-start_index]=f
	endfor
	
	Duplicate/O freqWave freq
	
	print nameofwave(freqWave)
End




Function GetTimes(directory_name, series_name, start_index, stop_index)
	String directory_name, series_name
	Variable start_index, stop_index
	
	String timeWaveName = "time_"+series_name
	Make/O/N=(stop_index-start_index+1) $timeWaveName
	Wave timeWave = $timeWaveName

	NewPath/O/Q myPath directory_name

	String greped = ""
	Variable hours, minutes, seconds, timetotal
	String thisFile
	
	Variable i
	for(i = start_index; i<=stop_index; i+=1)
		thisFile = series_name+num2str(i)+".txt"
	
		Grep/LIST/Q/P=myPath/E="#\\d\:" thisFile
		greped = S_value
		sscanf greped, "#%f:%f:%f", hours, minutes, seconds
		timetotal = hours*3600+minutes*60+seconds
		//printf "%f:%f:%f	%f\r", hours, minutes, seconds, timetotal
		timeWave[i-start_index]=timetotal
	endfor
	
	//Differentiate timeWave
	
	Duplicate/O timeWave times
	
	print nameofwave(timeWave)
End



Function GetVoltagesOld(directory_name, series_name, start_index, stop_index)
    String directory_name, series_name
    Variable start_index, stop_index
    
    String HVWaveName = "HV_"+series_name
    Make/O/N=(stop_index-start_index+1) $HVWaveName
    Wave HVWave = $HVWaveName

    NewPath/O/Q myPath directory_name

    String greped = ""
    Variable f
    String thisFile
    
    Variable i
    Variable j = 0
    for(i = start_index; i<=stop_index; i+=1)
        thisFile = series_name+num2str(i)+".txt"
    
        Grep/LIST/Q/P=myPath/E="# HV" thisFile
        greped = S_value
        sscanf greped, "# HV %f", f
        //print f
        HVWave[j]=f
        j+=1
    endfor
    
    print nameofwave(hvwave)
End



Function GetVoltages(directory_name, series_name, start_index, stop_index)
    String directory_name, series_name
    Variable start_index, stop_index
    
    String HVpWaveName = "HVp_"+series_name
    String HVmWaveName = "HVm_"+series_name
    Make/O/N=(stop_index-start_index+1) $HVpWaveName, $HVmWaveName
    Wave HVpWave = $HVpWaveName; Wave HVmWave = $HVmWaveName

    NewPath/O/Q myPath directory_name

    String greped = ""
    Variable f
    String thisFile
    
    Variable i
    Variable j = 0
    for(i = start_index; i<=stop_index; i+=1)
        thisFile = series_name+num2str(i)+".txt"
    
        Grep/LIST/Q/P=myPath/E="# HV+" thisFile
        greped = S_value
        sscanf greped, "# HV+ %f", f
        //print f
        HVpWave[j]=f
      //  j+=1
        
        Grep/LIST/Q/P=myPath/E="# HV-" thisFile
        greped = S_value
        sscanf greped, "# HV- %f", f
        //print f
        HVmWave[j]=f
        j+=1
    endfor
    
    print nameofwave(hvpwave)
    print nameofwave(hvmwave)
End



Function GetPositions(directory_name, series_name, start_index, stop_index)
    String directory_name, series_name
    Variable start_index, stop_index
    
    String PosWaveName = "Pos_"+series_name
    Make/O/N=(stop_index-start_index+1) $PosWaveName
    Wave PosWave = $PosWaveName

    NewPath/O/Q myPath directory_name

    String greped = ""
    Variable f
    String thisFile
    
    Variable i
    Variable j = 0
    for(i = start_index; i<=stop_index; i+=1)
        thisFile = series_name+num2str(i)+".txt"
    
        Grep/LIST/Q/P=myPath/E="# Motor position" thisFile
        greped = S_value
        sscanf greped, "# Motor position %f", f
        //print f
        PosWave[j]=f
        j+=1
    endfor
    
    print nameofwave(poswave)
End



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function load(dirstr)
	string dirstr
	LoadWave/q/a/G/D/w/o dirstr 			// loads a single data file and assigns each column to a wave. Counts -> counts. Laser IFM voltage -> voltage. All other named columns also loaded by name
	If(V_flag==0)
		If(stringmatch(dirstr, "*laser*"))
			Abort "Could not locate laser_calib file. Check path and file names."
		Else
			Abort "Could not locate file. Check path and file names.\r" + dirstr
		EndIf
	EndIf 
end
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



Function ReindexToTOF(tof)
	variable tof
	Wave voltage
	
	duplicate/o voltage voltageTOF
	voltageTOF= voltage[x-tof]
	
	variable n = numpnts(voltageTOF)
	DeletePoints (n-5),5,voltagetof,counts		// Delete the first and last 5 points of each wave.
	DeletePoints 0,5,voltagetof,counts	
End



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function AdvancedFilter()   
	Wave counts
	
	Variable lowPassStart = 0.1
	Variable rejectStart = 0.1
	Variable coefficients = 5
	
	Duplicate/o counts counts_noFilter
	
	Make/O/D/N=0 coefs
	FilterFIR/DIM=0/LO={0.1,0.1,5}/WINF=KaiserBessel20/COEF coefs, counts
End






function filter()
	variable imax, discard_spike=0		// discard_spike keeps a tally of many spikes have been deleted for this data file
	wave counts
	variable c,i,thresh= 1.2  			// set this threshold higher to kill more spikes
	for(c=0; c<3; c=c+1)				// run kill spike loop 3 times
		variable curvature=0
		for(i=0; i< numpnts(counts); i=i+1)							// for each element of counts...
			curvature = (counts[i+1] - 2*counts[i] + counts[i-1])		// calculate the curvature at point i...
			if(curvature < -thresh * min(min(counts[i-1],counts[i]), counts[i+1]) ) // find the minimum of counts[i] and its nearest neighbors, compare to curvature
				DeletePoints i, 1, counts								// delete the offending point in counts
				DeletePoints i, 1, voltagetof							// ... should be i-tof for the laser voltage signal
				discard_spike = discard_spike+1 						// increament spike kill tally
			endif
		endfor
	endfor
	print "killspikes deleted ", discard_spike, " points"
end



//function filter2()
//	wave atom
//	wavestats/q atom
//	variable thresh=2.5,i=0
//	for(i=0;i<numpnts(atom);i=i+1)
//		if(atom[i]>(thresh*V_avg))
//			atom[i] = V_avg
//		endif
//	endfor
//end
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function fringe(laser, atom, Lmax, Lmin, Lavg, Lamp, bins)
	wave laser, atom
	variable Lmax, Lmin, Lavg, Lamp, bins

	//determine grating position
	duplicate/o laser gratingPosition	// creates a new wave with same dimensions as laser
	variable k_o =   0.00186539  		// optical grating calibration done with HeNe //variable d,lambda=632.8,k_o //d=lambda/(sin(atan2(7.125,37.25))) //k_o = 2*pi/d
	gratingPosition = 1/k_o*asin((laser - Lavg)/Lamp)	// horizontal position of grating (units of nm). Index is time in ms.
	

	wavestats/q gratingPosition
	variable gratingMaxPos = V_max ;  variable gratingMinPos = V_min 				// max and min grating positions in nm

	// we currently have atom intensity vs time
	// we will now transform our data to atom intensity vs grating position
	make/o/d/n=(bins) binnedposition, atombin, atombin_error 
	binnedposition = (x/bins)* (gratingMaxPos - gratingMinPos) + gratingMinPos	// binnedposition ranges from gratingMinPos to gratingMaxPosition and has "bins" number of points
	variable j,n, i
	atombin = 0																		// initialize atombin, the average number of counts/ms at each binned grating position
	for(j=0;j<bins;j+=1)															// for each bin...
		n=0																				// create a tally for number of times counts are dumped into this jth bin
		for(i=0;i<numpnts(gratingPosition);i+=1)										// for each ms of data...
			if(gratingPosition[i]>=binnedposition[j] && gratingPosition[i]<binnedposition[j+1])// if the grating position falls within this jth bin
				atombin[j] += atom[i]															// place the counts for this ms into the counts for this position
				n+=1																			// increment the tally for number of times counts were put into this bin
			endif																			// do nothing otherwise
		endfor																			// continue with the next ms of data until done
		atombin_error[j] = sqrt(atombin[j])/n										// calculate the statistical error in number of counts/ms at this grating position
		atombin[j]=atombin[j]/n														// calculate the average number of counts/ms at this grating position
	endfor																			// continue with the next grating position until done
end





function fringe2(laser, atom, Lmax, Lmin, Lavg, Lamp, resampleRate)
	wave laser, atom
	variable Lmax, Lmin, Lavg, Lamp, resampleRate

	//determine grating position
	duplicate/o laser gratingPosition	// creates a new wave with same dimensions as laser
	variable k_o =   0.00186539  		// optical grating calibration done with HeNe //variable d,lambda=632.8,k_o //d=lambda/(sin(atan2(7.125,37.25))) //k_o = 2*pi/d
	gratingPosition = 1/k_o*asin((laser - Lavg)/Lamp)	// horizontal position of grating (units of nm). Index is time in ms.
	
	duplicate/O atom counts_error
	counts_error = sqrt(atom)
			
	Resample/RATE=(resampleRate) gratingposition, atom, counts_error
	
	//counts_error *=2
	//counts_error = sqrt(atom)
end




	// we currently have atom flux vs time
	// fringeTrim transforms our data to atom flux vs grating position
function fringeTrim(laser, atom, Lavg, Lamp, resolution)
	wave laser, atom
	variable Lavg, Lamp, resolution

	//determine grating position
	duplicate/o laser gratingPosition	// creates a new wave with same dimensions as laser
	variable k_o =   0.00186539  		// optical grating calibration done with HeNe //variable d,lambda=632.8,k_o //d=lambda/(sin(atan2(7.125,37.25))) //k_o = 2*pi/d
	gratingPosition = 1/k_o*asin((laser - Lavg)/Lamp)	// horizontal position of grating (units of nm). Index is time in ms (if sample rate is 1 kHz as usual).
			
	wavestats/q gratingPosition
	variable gratingMaxPos = round(V_max) ;  variable gratingMinPos = round(V_min) 				// max and min grating positions in nm
	
	Variable bins = round((gratingMaxPos-gratingMinPos) / resolution)
	
	Variable filterType=2
	Variable countSetStdDevFilter = 2
	Variable trimmedMeanPercent = 0
	Variable samplesPerPositionThreshold = 5

	make/o/d/n=(bins) binnedposition, atombin, atombin_error 
	binnedposition = gratingMinPos + resolution*p	// binnedposition ranges from gratingMinPos to gratingMaxPosition and has "bins" number of points
	atombin = NaN
	atombin_error = NaN
	
	variable i																	// initialize atombin, the average number of counts/ms at each binned grating position
	for(i=0;i<bins;i+=1)															// for each bin...
		Make/O/D/N=0 thisCountSet
		Extract/O atom, thisCountSet, ( gratingPosition[p] > (binnedposition[i] -resolution/2) && gratingPosition[p] <= (binnedposition[i]+resolution/2) ) 
		
		if(numpnts(thisCountSet)!=0)
			WaveStats/Z/Q thisCountSet
			
			if(filterType==0)
				//do nothing -- no filtering
			elseif(filterType==1)	// std. dev. filtering
				thisCountSet = (thisCountSet[p] > v_avg-countSetStdDevFilter*v_sdev && thisCountSet[p] < v_avg+countSetStdDevFilter*v_sdev) ? thisCountSet[p] : NaN
				Wavestats/Z/Q thisCountSet
			elseif(filterType==2)	// trimed mean computation
				Sort thisCountSet, thisCountSet		// sorts from lowest to highest, puts nans at the end
				Variable thisCountSetPnts =v_npnts	// number of pnts that are not nan's or inf's
				Variable pntsCut = round(thisCountSetPnts*trimmedMeanPercent/100)	// number of pnts to cut from each end of the set. round up.
				if(pntsCut==0)	// always throw out at least one point
					pntsCut=1
				endif
				//print "trimmed points from each end of the set: ", pntsCut
				thisCountSet[0,pntsCut-1] = NaN
				thisCountSet[v_npnts-pntsCut,v_npnts] = NaN
				Wavestats/Z/Q thisCountSet
			endif
				
				
			if(v_npnts>= samplesPerPositionThreshold)	
				atombin[i] = v_avg
				atombin_error[i] = v_sem
			endif

		endif
	endfor																			// continue with the next grating position until done
end
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








Function randomInitGuessesFringeFit()
	Wave w_coef, w_coef2k
	w_coef[0]=enoise(50)+75
	w_coef[1]=enoise(pi/2)+pi	
	w_coef[2]=enoise(0.1)+0.15
	
	w_coef2k[0]=enoise(50)+75
	w_coef2k[1]=enoise(pi/2)+pi	
	w_coef2k[2]=enoise(0.1)+0.15
	w_coef2k[4]=enoise(pi/2)+pi	
	w_coef2k[5]=enoise(0.02)+0.03
end




Function FitCountsVsGratingPosition(ind, start, stop, twoKyesno, suppressUpdateYesNo)
	variable ind, start, stop, twoKyesno, suppressUpdateYesNo
	wave contrast, contrast_error, counts, gratingposition, counts_error//, w_sigma
	wave phase, phase_error, avg_counts, avg_counts_error, chisqrdDOF
	
	Variable countsPnts = numpnts(counts)
	
	//Redimension/N=4 w_coef, w_sigma
	wave w_coef, w_sigma

	Make/O/T Constraints1k2k = {"K1>-0.01", "K1<6.3", "K2>0.0001", "K2<0.5"}
	
	Make/O eps1k = {0.1,1,0.001,0.01}
	
	Make/O eps1k2k = {0.1,1,0.001,0.01,0.01,0.001}
	
	//Make/O/T Constraints = {"K2>.001"}
	
	//randomInitGuessesFringeFit() // randomize initial guesses. causes fit to be much slower, but perhaps more accurate?
		
	// call fitting function "interference". do not display results in history. let k0,k1,k2 vary; hold k3 constant. w_coef is the fitting coefficients wave
	
	//FuncFit/Q/H="0001" interferenceAAO w_coef Counts /X=gratingPosition /I=1 /D ///C=Constraints
	FuncFit/Q/H="0001"/N=(suppressUpdateYesNo) interferenceAAO w_coef Counts /X=gratingPosition /W=counts_error /I=1 /D /R /E=eps1k ///C=Constraints
	// "inteferenceAAO" fits to a sine wave of the form counts = avg_counts*(1 + contrast * sin ( 2pi/period * x + phase )
	// AAO stands for all-at-once fit. this is more efficient than the standard point-by-point fits when dealing with larger numbers of points because it eliminates most of the function call overhead.
	
	Variable bestFitPhase = w_coef[1]
	bestFitPhase = mod(bestFitPhase,2*pi)
	bestFitPhase += 2*pi
	bestFitPhase = mod(bestFitPhase,2*pi)
	
	// store the fit coefficients and their uncertainties 
	contrast[ind-start] = abs(w_coef[2])
	contrast_error[ind-start] = w_sigma[2]
	phase[ind-start] = bestFitPhase
	phase_error[ind-start] = w_sigma[1]
	avg_counts[ind-start] = w_coef[0]
	avg_counts_error[ind-start] = w_sigma[0] 
	
	chisqrdDOF[ind-start]= v_chisq/(numpnts(counts)-1-3)
	
	if(twokyesno)
		wave contrast1k2k, contrast_error1k2k, phase1k2k, phase_error1k2k, avg_counts1k2k, avg_counts_error1k2k, phase2k2k, phase_error2k2k, contrast2k2k, contrast_error2k2k, w_coef2k, chisqrdDOF2k//, w_sigma2k
	
		//randomInitGuessesFringeFit() 
	
		FuncFit/Q/H="000100"/N=(suppressUpdateYesNo) interferenceAAO2 w_coef2k Counts /X=gratingPosition /W=counts_error /I=1 /D /R /E=eps1k2k ///C=Constraints1k2k
	
		if(w_sigma[1] > 0.1)	//	if phase error is greater than 100 mrad then refit data with new initial guesses
			//randomInitGuessesFringeFit() 
			FuncFit/Q/H="000100"/N=(suppressUpdateYesNo) interferenceAAO2 w_coef2k Counts /X=gratingPosition /W=counts_error /I=1 /D /R  /E=eps1k2k// /C=Constraints1k2k
			//print ind
		endif
	
		bestFitPhase = w_coef2k[1]
		bestFitPhase = mod(bestFitPhase,2*pi)
		bestFitPhase += 2*pi
		bestFitPhase = mod(bestFitPhase,2*pi)
	
		contrast1k2k[ind-start] = abs(w_coef2k[2])
		contrast_error1k2k[ind-start] = w_sigma[2]
		phase1k2k[ind-start] = bestFitPhase
		phase_error1k2k[ind-start] = w_sigma[1]
		avg_counts1k2k[ind-start] = w_coef2k[0]
		avg_counts_error1k2k[ind-start] = w_sigma[0] 
		
		bestFitPhase = w_coef2k[4]
		bestFitPhase = mod(bestFitPhase,2*pi)
		bestFitPhase += 2*pi
		bestFitPhase = mod(bestFitPhase,2*pi)
		phase2k2k[ind-start] = bestfitphase
		phase_error2k2k[ind-start] = w_sigma[4]

		contrast2k2k[ind-start] = abs(w_coef2k[5])
		contrast_error2k2k[ind-start] = w_sigma[5]
		
		chisqrdDOF2k[ind-start] = v_chisq/(numpnts(counts)-1-5)
	endif
	
	Variable minScanAmplitude = 175	//nm
	Wavestats/Q gratingPosition
	Variable scanAmplitude = v_max-v_min
	//print scanamplitude
	if(scanAmplitude < minScanAmplitude)
		if(twokyesno)
			contrast1k2k[ind-start] = NaN
			contrast_error1k2k[ind-start] = NaN
			phase1k2k[ind-start] = NaN
			phase_error1k2k[ind-start] = NaN
			avg_counts1k2k[ind-start] = NaN
			avg_counts_error1k2k[ind-start] = NaN
	
			phase2k2k[ind-start] = NaN
			phase_error2k2k[ind-start] = NaN
			contrast2k2k[ind-start] = NaN
			contrast_error2k2k[ind-start] = NaN
		
			contrast[ind-start] = NaN
			contrast_error[ind-start] = NaN
			phase[ind-start] = NaN
			phase_error[ind-start] = NaN	
			avg_counts[ind-start] = NaN
			avg_counts_error[ind-start] = NaN 
		else	
			contrast[ind-start] = NaN
			contrast_error[ind-start] = NaN
			phase[ind-start] = NaN
			phase_error[ind-start] = NaN
			avg_counts[ind-start] = NaN
			avg_counts_error[ind-start] = NaN 
		endif
	endif
	
	//print v_chisq/(countsPnts-1)
	
	//printf "index = %g     phase = %g      contrast = %g\r" ind, bestFitPhase, abs(w_coef[2])
end




Function interferenceAAO(pw, yw, xw) : FitFunc
	Wave pw, yw, xw
	
	Variable avgCounts = abs(pw[0])
	Variable phase = pw[1]
	Variable contrast = abs(pw[2])
	Variable period = pw[3]
	
	Multithread yw = avgCounts*(1 + contrast*sin(2*pi/period * xw + phase))
End


Function interferenceAAO2(pw, yw, xw) : FitFunc
	Wave pw, yw, xw
	
	Variable avgCounts = abs(pw[0])
	Variable phase = pw[1]
	Variable contrast = abs(pw[2])
	Variable period = pw[3]
	Variable phase2 = pw[4]
	Variable contrast2 = pw[5]
	
	Multithread yw = avgCounts*(1 + contrast*sin(2*pi/period * xw + phase) + contrast2*sin(2*pi/period*xw*2 + phase2) )
End



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function sinfit(ind, start, stop, bins)
	variable ind, start, stop, bins
	wave contrast, contrast_error, atombin, binnedposition, atombin_error//, w_sigma
	wave phase, phase_error, avg_counts, avg_counts_error
	
	// make a wave with four elements, one each fit parameter
	Make/D/N=4/O w_coef, w_sigma				// w[0]= avg_counts, w[1]=phase, w[2]=contrast, w[3]=grating period
	W_coef[0] = {20, 0, .2, 100}		// initial guesses for the fit parameters
	
	Make/O/T/N=4 Constraints
	Constraints = {"K0>.01", "K2>.001"}
	
	// call fitting function "interference". do not display results in history. let k0,k1,k2 vary; hold k3 constant. w_coef is the fitting coefficients wave
	// fit atombin vs binnedposition. do not fit the first or last two bins. weight points by atombin_error. atombin_error contains std devs. autoname the fit wave
	variable pntsToExclude = 2	//number of points to exclude from the sides of the binned fringe fit
	
	FuncFit/Q/H="0001" interference w_coef atombin[pntsToExclude, bins-2-pntsToExclude] /X=binnedposition /W=atombin_error /I=1 /D ///C=Constraints
	// "inteference" fits to a sine wave of the form intensity = avg_counts*(1 + contrast * sin ( 2pi/period * x + phase )
	
	Variable bestFitPhase = w_coef[1]
	bestFitPhase = mod(bestFitPhase,2*pi)
	bestFitPhase += 2*pi
	bestFitPhase = mod(bestFitPhase,2*pi)	

	// store the fit coefficients and their uncertainties 
	contrast[ind-start] = abs(w_coef[2])
	contrast_error[ind-start] = w_sigma[2]
	phase[ind-start] = bestFitPhase
	phase_error[ind-start] = w_sigma[1]
	avg_counts[ind-start] = w_coef[0]
	avg_counts_error[ind-start] = w_sigma[0] 
	//print ind
end


Function interference(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = i_o*(1 + abs(visibility)*sin(2*pi/tau*x + phi))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = i_o
	//CurveFitDialog/ w[1] = phi
	//CurveFitDialog/ w[2] = visibility
	//CurveFitDialog/ w[3] = tau
	
	return abs(w[0])*(1 + abs(w[2])*sin(2*pi/abs(w[3])*x + w[1]))
End
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function display_fringes2(ind)
	string ind
	
	NVAR fringePlotMax, fringePlotMin
	
	display/k=1/W=(40,50,800,500) counts vs gratingposition
	ModifyGraph mode=3,marker=8
	ErrorBars counts Y,wave=(counts_error, counts_error)
	SetAxis left fringePlotMin,fringePlotMax  //  set range for fringes
	ModifyGraph grid(bottom)=1
	SetAxis bottom -500,500 
	ModifyGraph nticks(bottom)=8
	TextBox/C/N=text0/X=0/Y=0 ind	
	Label bottom "grating position [nm]"
	Label left "atom kcounts/sec"
	make/o/n=0 fit_counts; appendToGraph fit_counts  // when fit_counts is populated by the best-fit function it will already be on the graph
	ModifyGraph rgb(fit_counts)=(0,0,0)
	ModifyGraph mode(counts)=3
end




Function display_fringes(ind)
	string ind
	display/k=1 atombin vs binnedposition
	ModifyGraph mode=3,marker=8
	ErrorBars atombin Y,wave=(atombin_error,atombin_error)
	SetAxis left 0,120  //  set range for fringes
	ModifyGraph grid(bottom)=1
	SetAxis bottom -500,500 
	ModifyGraph nticks(bottom)=8
	TextBox/C/N=text0/X=100.00/Y=100.00 ind
	Label bottom "grating position [nm]"
	Label left "atom kcounts/sec"
	duplicate/o atombin fit_atombin; appendToGraph fit_atombin  // for use later
	ModifyGraph rgb(fit_atombin)=(0,0,0)
	ModifyGraph mode(atombin)=3
end


Function display_avg_counts_results()
	display/k=1 avg_counts vs file_index
	ModifyGraph mode=3,marker=8;DelayUpdate
	ErrorBars avg_counts Y,wave=(avg_counts_error,avg_counts_error)
	Label bottom "file index"
	Label left "avg_counts"	
end


Function display_phase_results()
	display/k=1 phase vs file_index
	SetAxis left 0,6.29
	ModifyGraph mode=3,marker=8;DelayUpdate
	ErrorBars phase Y,wave=(phase_error,phase_error)
	Label bottom "file index"
	Label left "phase"
end


Function display_contrast_results()
	display/k=1 contrast vs file_index
	SetAxis left 0,.25
	ModifyGraph mode=3,marker=8;DelayUpdate
	ErrorBars contrast Y,wave=(contrast_error,contrast_error)
	Label bottom "file index"
	Label left "contrast"
end
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




Function SaveFringeData(appended_name, overwriteYesno, displayYesNo)
	String appended_name
	Variable overwriteyesno, displayYesNo
	
	String newphase = "phase_" + appended_name
	String newphase_error = "phase_error_" + appended_name
	String newcontrast = "contrast_" + appended_name
	String newcontrast_error = "contrast_error_" + appended_name
	String newavg_counts = "avg_counts_" + appended_name
	String newavg_counts_error = "avg_counts_error_" + appended_name
	String newfile_index = "file_index_" + appended_name
	
	if(overwriteyesno==0)
		Duplicate phase $newphase
		Duplicate phase_error $newphase_error
		Duplicate contrast $newcontrast
		Duplicate contrast_error $newcontrast_error
		Duplicate avg_counts $newavg_counts
		Duplicate avg_counts_error $newavg_counts_error
		Duplicate file_index $newfile_index
	endif
	if(overwriteyesno==1)
		Duplicate/o phase $newphase
		Duplicate/o phase_error $newphase_error
		Duplicate/o contrast $newcontrast
		Duplicate/o contrast_error $newcontrast_error
		Duplicate/o avg_counts $newavg_counts
		Duplicate/o avg_counts_error $newavg_counts_error
		Duplicate/o file_index $newfile_index
	endif
	
	Wave phase_saved = $newphase
	Wave phase_error_saved = $newphase_error
	Wave contrast_saved = $newcontrast
	Wave contrast_error_saved = $newcontrast_error
	Wave avg_counts_saved = $newavg_counts
	Wave avg_counts_error_saved = $newavg_counts_error
	Wave file_index_saved = $newfile_index
	
	if(0)
		display/k=1 avg_counts_saved vs file_index_saved
		ModifyGraph mode=3,marker=8;DelayUpdate
		ErrorBars $newavg_counts Y,wave=(avg_counts_error_saved,avg_counts_error_saved)
		Label bottom "file index " + appended_name
		Label left "avg_counts "	+ appended_name
		WindowNamer("Counts " + appended_name)
	
		display/k=1 phase_saved vs file_index_saved
		SetAxis left 0,6.29
		ModifyGraph mode=3,marker=8;DelayUpdate
		ErrorBars $newphase Y,wave=(phase_error_saved,phase_error_saved)
		Label bottom "file index " + appended_name
		Label left "phase " + appended_name
		WindowNamer("Phase " + appended_name)
	
		display/k=1 contrast_saved vs file_index_saved
		SetAxis left 0,.25
		ModifyGraph mode=3,marker=8;DelayUpdate
		ErrorBars $newcontrast Y,wave=(contrast_error_saved,contrast_error_saved)
		Label bottom "file index " + appended_name
		Label left "contrast " + appended_name
		WindowNamer("Contrast " + appended_name)
	endif
	
	if(displayYesNo)
	Display/K=1/W=(40,50,1200,650) avg_counts_saved vs file_index_saved
	ErrorBars $newavg_counts Y,wave=(avg_counts_error_saved,avg_counts_error_saved)
	Label bottom "file index " + appended_name
	Label left "avg_counts "	+ appended_name
	AppendToGraph/l=phase phase_saved vs file_index_saved
	AppendToGraph/l=contrast contrast_saved vs file_index_saved
	ModifyGraph axisEnab(left)={0,0.35},axisEnab(phase)={0.35,0.65},freePos(phase)=0,axisEnab(contrast)={0.70,1},freePos(contrast)=0
	Label phase "phase " + appended_name
	Label contrast "contrast " + appended_name
	SetAxis phase 0,6.29
	SetAxis contrast 0,.3
	ErrorBars $newphase Y,wave=(phase_error_saved,phase_error_saved)
	ErrorBars $newcontrast Y,wave=(contrast_error_saved,contrast_error_saved)
	ModifyGraph mode=3,marker=8
	ModifyGraph rgb($newavg_counts)=(0,0,0),rgb($newphase)=(1,4,52428)
	ModifyGraph lblPosMode(phase)=1
	ModifyGraph lblPosMode(contrast)=1
	ModifyGraph gfSize=10
	ModifyGraph gmSize=4
	ModifyGraph nticks(bottom)=20,minor(bottom)=1
	ModifyGraph grid(bottom)=2
	ShowInfo
	WindowNamer("Fringe summary " + appended_name)
	endif
End



Function SaveFringeData2k(appended_name, overwriteyesno, display2kYesNo)
	String appended_name
	Variable overwriteyesno, display2kYesNo
	
	// 1k2k = 1k component of a 2k fit
	// 2k2k = 2k component of a 2k fit
	String newphase = "phase1k2k_" + appended_name
	String newphase_error = "phase_error1k2k_" + appended_name
	String newcontrast = "contrast1k2k_" + appended_name
	String newcontrast_error = "contrast_error1k2k_" + appended_name
	String newavg_counts = "avg_counts1k2k_" + appended_name
	String newavg_counts_error = "avg_counts_error1k2k_" + appended_name
	String newchisqrdDOF = "chisqrdDOF1k2k_" + appended_name
	String newavg_power = "avg_power_" + appended_name
	String newavg_powerStdDev = "avg_powerStdDev_" + appended_name
	String newphase2k = "phase2k2k_" + appended_name
	String newphase_error2k = "phase_error2k2k_" + appended_name
	String newcontrast2k = "contrast2k2k_" + appended_name
	String newcontrast_error2k = "contrast_error2k2k_" + appended_name
	String newfile_index = "file_index_" + appended_name
	
	if(overwriteyesno==0)
		Duplicate phase1k2k $newphase
		Duplicate phase_error1k2k $newphase_error
		Duplicate contrast1k2k $newcontrast
		Duplicate contrast_error1k2k $newcontrast_error
		Duplicate phase2k2k $newphase2k
		Duplicate phase_error2k2k $newphase_error2k
		Duplicate contrast2k2k $newcontrast2k
		Duplicate contrast_error2k2k $newcontrast_error2k
		Duplicate avg_counts1k2k $newavg_counts
		Duplicate avg_counts_error1k2k $newavg_counts_error
		Duplicate file_index $newfile_index
		Duplicate chisqrdDOF $newchisqrdDOF
		Duplicate avg_power $newavg_power
		Duplicate avg_powerStdDev $newavg_powerStdDev
	endif
	if(overwriteyesno==1)
		Duplicate/O phase1k2k $newphase
		Duplicate/O phase_error1k2k $newphase_error
		Duplicate/O contrast1k2k $newcontrast
		Duplicate/O contrast_error1k2k $newcontrast_error
		Duplicate/O phase2k2k $newphase2k
		Duplicate/O phase_error2k2k $newphase_error2k
		Duplicate/O contrast2k2k $newcontrast2k
		Duplicate/O contrast_error2k2k $newcontrast_error2k
		Duplicate/O avg_counts1k2k $newavg_counts
		Duplicate/O avg_counts_error1k2k $newavg_counts_error
		Duplicate/O file_index $newfile_index
		Duplicate/O chisqrdDOF $newchisqrdDOF
		Duplicate/O avg_power $newavg_power
		Duplicate/O avg_powerStdDev $newavg_powerStdDev
	endif
	
	Wave phase_saved = $newphase
	Wave phase_error_saved = $newphase_error
	Wave contrast_saved = $newcontrast
	Wave contrast_error_saved = $newcontrast_error
	Wave phase_saved2k = $newphase2k
	Wave phase_error_saved2k = $newphase_error2k
	Wave contrast_saved2k = $newcontrast2k
	Wave contrast_error_saved2k = $newcontrast_error2k
	Wave avg_counts_saved = $newavg_counts
	Wave avg_counts_error_saved = $newavg_counts_error
	Wave file_index_saved = $newfile_index
	Wave chisqrdDOF_saved = $newchisqrdDOF
	Wave avg_power_saved = $newavg_power
	Wave avg_powerStdDev_saved = $newavg_powerStdDev
	
	if(0)
		display/k=1 avg_counts_saved vs file_index_saved
		ModifyGraph mode=3,marker=8;DelayUpdate
		ErrorBars $newavg_counts Y,wave=(avg_counts_error_saved,avg_counts_error_saved)
		Label bottom "file index " + appended_name
		Label left "avg_counts "	+ appended_name
		WindowNamer("Counts " + appended_name)
	
		display/k=1 phase_saved vs file_index_saved
		SetAxis left 0,6.29
		ModifyGraph mode=3,marker=8;DelayUpdate
		ErrorBars $newphase Y,wave=(phase_error_saved,phase_error_saved)
		Label bottom "file index " + appended_name
		Label left "phase " + appended_name
		WindowNamer("Phase " + appended_name)
	
		display/k=1 contrast_saved vs file_index_saved
		SetAxis left 0,.25
		ModifyGraph mode=3,marker=8;DelayUpdate
		ErrorBars $newcontrast Y,wave=(contrast_error_saved,contrast_error_saved)
		Label bottom "file index " + appended_name
		Label left "contrast " + appended_name
		WindowNamer("Contrast " + appended_name)
	endif
	
	
	// for 1k component of the 2k included fit
	if(1)	// 1 for including chi sqrd / dof
	// display the traces
	Display/K=1/W=(40,50,1200,650) avg_counts_saved vs file_index_saved
	AppendToGraph/l=phase phase_saved vs file_index_saved
	AppendToGraph/l=contrast contrast_saved vs file_index_saved
	AppendToGraph/l=chi chisqrddof_saved vs file_index_saved
	
	// set the axes sizes
	ModifyGraph axisEnab(left)={0,0.33},axisEnab(phase)={0.37,0.65},freePos(phase)=0,axisEnab(contrast)={0.70,1}
	ModifyGraph standoff(chi)=0,axisEnab(left)={0,0.18},axisEnab(phase)={0.2,0.52}
	ModifyGraph axisEnab(contrast)={0.54,0.78},axisEnab(chi)={0.8,1}
	ModifyGraph freePos(chi)={0,bottom},freePos(contrast)=0
	
	Label bottom "file index " + appended_name
	Label left "avg counts "	+ appended_name
	Label phase "phase 1k2k " + appended_name
	Label contrast "contrast 1k2k " + appended_name
	Label chi "Chi Sqrd / DOF " + appended_name
	ModifyGraph lblPosMode(phase)=1
	ModifyGraph lblPosMode(contrast)=1
	ModifyGraph lblPosMode(chi)=1
	
	SetAxis phase 0,6.29
	SetAxis contrast 0,.3
	
	ErrorBars $newphase Y,wave=(phase_error_saved,phase_error_saved)
	ErrorBars $newcontrast Y,wave=(contrast_error_saved,contrast_error_saved)
	ErrorBars $newavg_counts Y,wave=(avg_counts_error_saved,avg_counts_error_saved)
	
	ModifyGraph mode=3,marker=8
	ModifyGraph rgb($newavg_counts)=(0,0,0),rgb($newphase)=(1,4,52428)
	ModifyGraph rgb($newchisqrdDOF)=(2,39321,1)
	
	ModifyGraph gfSize=10
	ModifyGraph gmSize=4
	
	ModifyGraph nticks(bottom)=20,minor(bottom)=1
	ModifyGraph grid(bottom)=2
	ShowInfo
	WindowNamer("Fringe summary 1k2k " + appended_name)
	else
	Display/K=1/W=(40,50,1200,650) avg_counts_saved vs file_index_saved
	ErrorBars $newavg_counts Y,wave=(avg_counts_error_saved,avg_counts_error_saved)
	Label bottom "file index " + appended_name
	Label left "avg_counts "	+ appended_name
	AppendToGraph/l=phase phase_saved vs file_index_saved
	AppendToGraph/l=contrast contrast_saved vs file_index_saved
	ModifyGraph axisEnab(left)={0,0.33},axisEnab(phase)={0.37,0.65},freePos(phase)=0,axisEnab(contrast)={0.70,1},freePos(contrast)=0
	Label phase "phase 1k2k " + appended_name
	Label contrast "contrast 1k2k " + appended_name
	SetAxis phase 0,6.29
	SetAxis contrast 0,.3
	ErrorBars $newphase Y,wave=(phase_error_saved,phase_error_saved)
	ErrorBars $newcontrast Y,wave=(contrast_error_saved,contrast_error_saved)
	ModifyGraph mode=3,marker=8
	ModifyGraph rgb($newavg_counts)=(0,0,0),rgb($newphase)=(1,4,52428)
	ModifyGraph lblPosMode(phase)=1
	ModifyGraph lblPosMode(contrast)=1
	ModifyGraph gfSize=10
	ModifyGraph gmSize=4
	ModifyGraph nticks(bottom)=20,minor(bottom)=1
	ModifyGraph grid(bottom)=2
	ShowInfo
	WindowNamer("Fringe summary 1k2k " + appended_name)
	endif
	
	if(display2kYesNo)
	//for 2k component of the 2k included fit
	Display/K=1/W=(40,50,1200,650) avg_counts_saved vs file_index_saved
	ErrorBars $newavg_counts Y,wave=(avg_counts_error_saved,avg_counts_error_saved)
	Label bottom "file index " + appended_name
	Label left "avg counts "	+ appended_name
	AppendToGraph/l=phase phase_saved2k vs file_index_saved
	AppendToGraph/l=contrast contrast_saved2k vs file_index_saved
	ModifyGraph axisEnab(left)={0,0.35},axisEnab(phase)={0.35,0.65},freePos(phase)=0,axisEnab(contrast)={0.70,1},freePos(contrast)=0
	Label phase "phase 2k2k " + appended_name
	Label contrast "contrast 2k2k " + appended_name
	SetAxis phase 0,6.29
	SetAxis contrast 0,.1
	ErrorBars $newphase2k Y,wave=(phase_error_saved2k,phase_error_saved2k)
	ErrorBars $newcontrast2k Y,wave=(contrast_error_saved2k,contrast_error_saved2k)
	ModifyGraph mode=3,marker=8
	ModifyGraph rgb($newavg_counts)=(0,0,0),rgb($newphase2k)=(1,4,52428)
	ModifyGraph lblPosMode(phase)=1
	ModifyGraph lblPosMode(contrast)=1
	ModifyGraph gfSize=10
	ModifyGraph gmSize=4
	ModifyGraph nticks(bottom)=20,minor(bottom)=1
	ModifyGraph grid(bottom)=2
	ShowInfo
	WindowNamer("Fringe summary 2k2k " + appended_name)
	endif
End





Function DisplayFringeData(appended_name)
	String appended_name
	Variable yesno
	
	String newphase = "phase_" + appended_name
	String newphase_error = "phase_error_" + appended_name
	String newcontrast = "contrast_" + appended_name
	String newcontrast_error = "contrast_error_" + appended_name
	String newavg_counts = "avg_counts_" + appended_name
	String newavg_counts_error = "avg_counts_error_" + appended_name
	String newfile_index = "file_index_" + appended_name
	
	Wave phase_saved = $newphase
	Wave phase_error_saved = $newphase_error
	Wave contrast_saved = $newcontrast
	Wave contrast_error_saved = $newcontrast_error
	Wave avg_counts_saved = $newavg_counts
	Wave avg_counts_error_saved = $newavg_counts_error
	Wave file_index_saved = $newfile_index
	
	display avg_counts_saved vs file_index_saved
	ModifyGraph mode=3,marker=8;DelayUpdate
	ErrorBars $newavg_counts Y,wave=(avg_counts_error_saved,avg_counts_error_saved)
	Label bottom "file index " + appended_name
	Label left "avg_counts "	+ appended_name
	
	display phase_saved vs file_index_saved
	SetAxis left 0,6.29
	ModifyGraph mode=3,marker=8;DelayUpdate
	ErrorBars $newphase Y,wave=(phase_error_saved,phase_error_saved)
	Label bottom "file index " + appended_name
	Label left "phase " + appended_name
	
	display contrast_saved vs file_index_saved
	SetAxis left 0,.25
	ModifyGraph mode=3,marker=8;DelayUpdate
	ErrorBars $newcontrast Y,wave=(contrast_error_saved,contrast_error_saved)
	Label bottom "file index " + appended_name
	Label left "contrast " + appended_name
End