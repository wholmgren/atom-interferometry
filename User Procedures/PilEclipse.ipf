#pragma rtGlobals=1		// Use modern global access method.

#include "Pol2Electrodes"
#include "Phase Choppers"
#include "Fringe Fit"


function PilEclipse()
	
	Wave/T seriesNames 
	//make/o/t seriesNames = {"d"}
//	Make/O/T fileNames = {"fLGstart","fLGend","gLGstart","gLGend","kLGstart","kLGend","lLGstart","lLGend","pLGstart","pLGend","qLGstart","qLGend"}
	
	Make/O/N=(numpnts(seriesNames)) offsetWaveStart, offsetWaveEnd
	
	variable i, j
	for(j=0; j<2; j+=1)
	for(i=0; i<numpnts(seriesNames); i+=1)
	
		String fileName = SelectString(j, seriesNames[i]+"LGstart", seriesNames[i]+"LGend")
	
		LoadWave/G/D/W/N/O/Q "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:110927:"+fileName
	
		String countsName = "counts"+fileName
		String countsErrName = "countsErr"+fileName
		String posName = "pos"+fileName
	
		Duplicate/O counts $countsName, $countsErrName
		Wave countsD = $countsName			//d for dummy, just a pointer to refer to the real wave
		Wave countsErr = $countsErrName
	
		Duplicate/O position $posName
		Wave pos = $posName
	
		countsErr = sqrt(countsD)
	
		display/k=1 countsD vs pos
	
	
		if(j==0)
			Make/O/D W_coef = {50,.02,-200,3}
			FuncFit/NTHR=0/TBOX=768 errfitneg W_coef  countsD /X=pos /W=countsErr /I=1 /D 
			offsetWaveStart[i] = w_coef[2]
		else
			Make/O/D W_coef = {50,.02,-250,3}
			FuncFit/NTHR=0/TBOX=768 errfit W_coef  countsD /X=pos /W=countsErr /I=1 /D 
			offsetWaveEnd[i] = w_coef[2]
		endif
	
	
		String fitName = "fit_"+countsName
		ModifyGraph lsize($fitName)=2,rgb($fitName)=(0,0,65535)
	
	endfor
	endfor	
end	






function PilEclipsePol()
	
	// Assumes that the current data folder is named eg '110927d'
	String CurrentDataFolder = GetDataFolder(0)
	String dateStr
	String seriesStr
	SplitString/E="'([0-9]+)([a-z]+)'" CurrentDataFolder, dateStr, seriesStr
	
	Variable/G offsetStartum, offsetEndum	// we will save the 50% flux best fit parameters in these variables
	
	variable j
	for(j=0; j<2; j+=1)
	
		String fileName = SelectString(j, seriesStr+"LGstart", seriesStr+"LGend")
	
		LoadWave/G/D/W/N/O/Q "Macintosh HD:Users:holmgren:Desktop:alcro2 downloads:" + dateStr+ ":" + fileName
	
		String countsName = "counts"+fileName
		String countsErrName = "countsErr"+fileName
		String posName = "pos"+fileName					// position wave comes in in microns
	
		Duplicate/O counts $countsName, $countsErrName
		Wave countsD = $countsName			//d for dummy, just a pointer to refer to the real wave
		Wave countsErr = $countsErrName
	
		Duplicate/O position $posName
		Wave pos = $posName
	
		countsErr = sqrt(countsD)
	
		display/k=1 countsD vs pos
		
		wavestats/Q pos
		
		// figure out which way the motor is moving and whether the counts are increasing or decreasing
		Variable motorDirection = pos[V_npnts]/abs(pos[V_npnts])		// +1 for positive direction, -1 for negative direction
		Variable countsDirection = (countsD[V_npnts]-countsD[0])/abs((countsD[V_npnts]-countsD[0]))
		
		// use the appropriate fit function with the appropriate initial guesses
		if(countsDirection*motorDirection==1)
			Make/O/D W_coef = {50,.02,200*motorDirection,3}
			FuncFit/NTHR=0/TBOX=768 errfit W_coef  countsD /X=pos /W=countsErr /I=1 /D 
		elseif(countsDirection*motorDirection==-1)
			Make/O/D W_coef = {50,.02,200*motorDirection,3}
			FuncFit/NTHR=0/TBOX=768 errfitneg W_coef  countsD /X=pos /W=countsErr /I=1 /D 	
		endif
		
		Variable/G reducedChiSqrd = V_chisq/(V_npnts-(numpnts(W_coef)))
		printf "Reduced Chi-Squared: %g\r\r", reducedChiSqrd
	
	
		if(j==0)
			offsetStartum = w_coef[2]
		elseif(j==1)
			offsetEndum = w_coef[2]
		endif	
	
		String fitName = "fit_"+countsName
		ModifyGraph lsize($fitName)=2,rgb($fitName)=(0,0,65535)

	endfor	
	
	
	//Calculate the eclipse and phase measurement positions
	NVAR firstPosLG, lastPosLG
	
	Variable/G firstPosum, lastPosum, startPosum, endPosum, startPosAbs, endPosAbs
	
	firstPosum = firstPosLG*0.4				// first phase measurement position relative to arbitrary LG 0 position in microns
	//offsetStartum = offsetStartLG*0.4			// 50% flux position relative to same arbitrary LG 0 position in microns
	startPosum = firstPosum-offsetStartum		// first phase measurement position relative to 50% flux position
	startPosAbs = 3810/2*-1*motorDirection+startPosum			// first phase measurement position relative to center of gradE region

	print "first phase measurement position relative to arbitrary LG 0 position in microns: ", firstPosum
	print "start 50% flux position relative to same arbitrary LG 0 position in microns: ", offsetStartum
	print "first phase measurement position relative to 50% flux position: ", startPosum
	print "first phase measurement position relative to center of gradE region: ", startPosAbs

	printf "\r"

	lastPosum = lastPosLG*0.4				// first phase measurement position relative to arbitrary LG 0 position in microns
	//offsetEndum = offsetEndLG*0.4			// 50% flux position relative to same arbitrary LG 0 position in microns
	endPosum = lastPosum-offsetEndum		// first phase measurement position relative to 50% flux position
	endPosAbs = 3810/2*motorDirection-offsetEndum			// first phase measurement position relative to center of gradE region

	print "last phase measurement position relative to arbitrary LG 0 position in microns: ", lastPosum
	//print "end 50% flux position relative to last phase measurement position: ", offsetEndum  //	n/a
	print "last phase measurement position relative to 50% flux position: ", offsetEndum
	print "last phase measurement position relative to center of gradE region: ", endPosAbs


end	