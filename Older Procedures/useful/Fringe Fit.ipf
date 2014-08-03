#pragma rtGlobals=1		// Use modern global access method.

#include ":WindowNamer"

// 		INSTRUCTIONS
// Set directory with fringe data
// set series name, laser calib name, start index and stop index
// if you want to analyze every file in a series you can set the stop index to 1 (this isn't a good idea when analyzing polarizability data since we don't want the extra files at the end confusing us)
// set the time of flight parameter
// run FringeFit()
// if you want a copy of your data, run SaveFringeData(appended_name, yesno).  yesno should be set to 0 unless you want to enable overwrites, in which case it should be 1





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function FringeFit()
	//pauseupdate ; silent 1

	// Choose your directory (uncomment the one you want, recomment the previously used one)
	//string directory_name = laptop directory goes here
	//string directory_name = "C:Documents and Settings:pamandalex:Desktop:mar7:df:"
	string directory_name = "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:100915 desktop:fringe data:"
	//string directory_name = vincents directory

	// Choose your data series and laser calibration file
	string series_name = "c"
	string laser_calib_name = "laser_calib_c1.txt"	
	//sprintf laser_calib_name, ("laser_calib_"+series_name+"1.txt")
	
	// Choose the data file range within the specified series
	variable start_index=1
	variable stop_index=160		//Specify stop_index = 0 if you want all the data for this series analyzed
									//doesn't yet work for series "l" or for "aa" etc.
	
	if(stop_index==0)
		stop_index = GetTotalFileNumber(directory_name, series_name)
	endif

	variable displayn=1						// display fringe plots? 1 = yes, 0 = no
	string fname,ind							// file name and index strings
	variable Lmax, Lmin, Lavg, Lamp 			// for laser calibration
	variable bins = 50   						// the number of grating position bins the counts will be placed in
	
	variable tof = 3								// tof = time of flight in integer ms. compensates for delay between when an atom is incident on the grating and when it is counted. 
	
	
	// prepare laser calibration 
	load(directory_name+laser_calib_name) 					// load laser calibration file
	
	duplicate/o voltage, laser_calib_voltage
	wavestats/q laser_calib_voltage  							// get wavestats and then calculate max, min, avg, amp.
	Lmax = V_max ;  Lmin = V_min ;  Lavg = (Lmax+Lmin)/2  ;   Lamp = (Lmax-Lmin)/2

	// make waves in which to store the contrast, phase, etc of each data file
	make/o/d/n=(stop_index-start_index+1) contrast contrast_error   phase phase_error   avg_counts avg_counts_error  file_index  simple_avg
	
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

	file_index = x+start_index									// set the next file to be analyzed
	display_avg_counts_results()								// display results
	display_phase_results()	
	display_contrast_results()	
	
	SaveFringeData(series_name,0)
end
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




Function GetTotalFileNumber(directory_name, series_name)
	String directory_name, series_name

	NewPath/O/Q myPath directory_name

	// Get a list containing all text files in the specified symbolic path.
	String fileList = IndexedFile(myPath, -1, ".txt")
 
	// Filter the list of files to include
	// only files that start with series_name
	fileList = ListMatch(fileList, series_name+"*", ";")

	Variable numFiles = ItemsInList(fileList, ";")

	if(stringmatch("l", series_name))
		return NaN
	else
		return numFiles	
	endif		
End






/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function load(dirstr)
	string dirstr
	LoadWave/q/a/G/D/w/o dirstr 			// loads a single data file and assigns each column to a wave. Counts -> Counts. Laser IFM voltage -> voltage.
	If(V_flag==0)
		Abort "Could not locate file. Check path and file names."
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
																// Does this work if tof != 0 and spikes have been removed???

	wavestats/q gratingPosition
	variable gratingMaxPos = V_max ;  variable gratingMinPos = V_min 				// max and min grating positions in nm

	// we currently have atom intensity vs time
	// we will now transform our data to atom intensity vs grating position
	make/o/d/n=(bins) binnedposition atombin atombin_error 
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function sinfit(ind, start, stop, bins)
	variable ind, start, stop, bins
	wave contrast, contrast_error, atombin, binnedposition, atombin_error//, w_sigma
	wave phase, phase_error, avg_counts, avg_counts_error
	
	// make a wave with four elements, one each fit parameter
	Make/D/N=4/O w_coef, w_sigma				// w[0]= avg_counts, w[1]=phase, w[2]=contrast, w[3]=grating period
	W_coef[0] = {12,0,.1,100}		// initial guesses for the fit parameters
	
	Make/O/T/N=4 Constraints
	Constraints = {"K0>.01", "K2>.001"}
	
	// call fitting function "interference". do not display results in history. let k0,k1,k2 vary; hold k3 constant. w_coef is the fitting coefficients wave
	// fit atombin vs binnedposition. do not fit the first or last two bins. weight points by atombin_error. atombin_error contains std devs. autoname the fit wave
	variable pntsToExclude = 2	//number of points to exclude from the sides of the binned fringe fit
	
	FuncFit/Q/H="0001" interference w_coef atombin[pntsToExclude, bins-2-pntsToExclude] /X=binnedposition /W=atombin_error /I=1 /D ///C=Constraints
	// "inteference" fits to a sine wave of the form intensity = avg_counts*(1 + contrast * sin ( 2pi/period * x + phase )
	
	// store the fit coefficients and their uncertainties 
	contrast[ind-start] = abs(w_coef[2])
	contrast_error[ind-start] = w_sigma[2]
	phase[ind-start] = w_coef[1]
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
Function display_fringes(ind)
			string ind
			display/k=1 atombin vs binnedposition
			ModifyGraph mode=3,marker=8
			ErrorBars atombin Y,wave=(atombin_error,atombin_error)
			SetAxis left 0,100  //  set range for fringes
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
	Wave phase
	phase = mod(phase,2*pi)
	phase += 2*pi
	phase = mod(phase,2*pi)
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




Function SaveFringeData(appended_name, yesno)
	String appended_name
	Variable yesno
	
	String newphase = "phase_" + appended_name
	String newphase_error = "phase_error_" + appended_name
	String newcontrast = "contrast_" + appended_name
	String newcontrast_error = "contrast_error_" + appended_name
	String newavg_counts = "avg_counts_" + appended_name
	String newavg_counts_error = "avg_counts_error_" + appended_name
	String newfile_index = "file_index_" + appended_name
	
	if(yesno==0)
	Duplicate phase $newphase
	Duplicate phase_error $newphase_error
	Duplicate contrast $newcontrast
	Duplicate contrast_error $newcontrast_error
	Duplicate avg_counts $newavg_counts
	Duplicate avg_counts_error $newavg_counts_error
	Duplicate file_index $newfile_index
	endif
	if(yesno==1)
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
	
	Display/K=1/W=(40,50,1200,700) avg_counts_saved vs file_index_saved
	ErrorBars $newavg_counts Y,wave=(avg_counts_error_saved,avg_counts_error_saved)
	Label bottom "file index " + appended_name
	Label left "avg_counts "	+ appended_name
	AppendToGraph/l=phase phase_saved vs file_index_saved
	AppendToGraph/l=contrast contrast_saved vs file_index_saved
	ModifyGraph axisEnab(left)={0,0.3},axisEnab(phase)={0.35,0.65},freePos(phase)=0,axisEnab(contrast)={0.70,1},freePos(contrast)=0
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
	WindowNamer("Contrast Phase and Avg Counts " + appended_name)
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