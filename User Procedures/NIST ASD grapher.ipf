#pragma rtGlobals=1		// Use modern global access method.

//	by Will Holmgren, 10/26/2011
//	Cronin group, Dept. of Physics, University of Arizona
//	holmgren@email.arizona.edu
//
//	Loads the output from the NIST, processes it, and displays a graph of the energy levels. Put cursors on two levels to calculate the transition wavelength
//
//	Instructions for generating the NIST ASD data:
//	1. go to: 	http://physics.nist.gov/PhysRefData/ASD/levels_form.html
// 	2. enter atom and ionization state (e.g. k i for neutral potassium)
//	3. set "Level Units" to eV
//	4. set "Format ouput" to ASCII
//	5. deselect Lande-g and Leading percentages check boxes under "Level information"
//	6. select, copy, and paste the output into a text file. Start with the column names line (	| Configuration | Term | J | Level |	). 
//		Do not include the line of "---" at the very top. End at the end of the next ionization state line (	| k ii | Limit | -- | 4.34 |	 ).
//	7. save the plain text file
//	8. set the fileToLoad path in the LoadASD() function
//	9. run the following command: goASD()
//	10. put cursors on the levels to calculate the associated transition wavelength.


Static Constant h =  6.62607e-34
Static Constant eCharge = 1.60217646e-19
Static Constant c = 299792458



Function LoadASD()
	String fileToLoad = "Macintosh HD:Users:holmgren:Documents:School:UA:Cronin labs:NIST Levels:balevels.txt"
	String columnFormatting = "C=3,F=-2; C=1,F=0;"
	LoadWave/A/O/J/D/W/E=1/K=0/V={",|"," $",0,0}/B=columnFormatting fileToLoad
	
	Variable i
	Wave/T Configuration, Term, JW
	Configuration = TrimWS(Configuration[p])
	Term = TrimWS(Term[p])
	JW = TrimWS(JW[p])
End


Function goASD()
	LoadASD()
	CleanupASD()
	MakeConfigTerm()
	LevelLines()
End


Function CleanupASD()
	Wave/T Configuration, Term, JW
	Wave Level//, Lande
	
	//Duplicate/O Level Level2
	
	Variable i,j
	Variable numNaNs = 0
	For(j=0;j<2;j+=1)
	For(i=0; i<numpnts(Level); i+=1)
		if(numtype(Level[i])==2)
			DeletePoints i, 1, Configuration, Term, JW, Level//, Lande, Level2
			numNans+=1
		endif
	EndFor
	EndFor
End


Function MakeConfigTerm()
	Wave/T Configuration, JW, Term
	
	String c0, cf
	String s1, s2
	
	Duplicate/O/T Configuration SimpleConfig, ConfigTerm, SimpleTerm, SimpleConfigTerm
	
	Variable i
	for(i=0; i<numpnts(configuration); i+=1)
		if(stringmatch(configuration[i],""))
			configuration[i] = configuration[i-1]
		endif
		SplitString /E="([[:alpha:][:digit:]]+)\.([[:alpha:][:digit:]]+)" configuration[i] , s1, s2
		if(stringmatch(s1,""))
			SplitString /E="([[:digit:]]+[[:alpha:]])([[:digit:]]+)" configuration[i] , s1, s2
			SimpleConfig[i] = configuration[i]
			ConfigTerm[i] =configuration[i] + " " + JW[i]
		else
			if(0)
				SimpleConfig[i] = s1+s2
				ConfigTerm[i] = s1+s2 + " " + JW[i]
			else
				SimpleConfig[i] = s1+"."+s2
				ConfigTerm[i] = s1+"."+s2 + " " + JW[i]
			endif
		endif
	endfor
	
	for(i=0; i<numpnts(Term); i+=1)
		if(stringmatch(Term[i],""))
			Term[i] = Term[i-1]
		endif
		if(stringmatch(Term[i],"Limit"))
			SplitString /E="([[:alpha:]]+) ([[:alpha:]]+)" Configuration[i] , s1, s2
			SimpleTerm[i]=s1+" Limit "+s2
			SimpleConfigTerm[i]=s1+" Limit "+s2
		else
			SplitString /E="([[:digit:]][[:alpha:]])" Term[i] , s1
			SimpleTerm[i] = s1+JW[i]
			SimpleConfigTerm[i] = SimpleConfig[i]+" "+ SimpleTerm[i]
		endif
	endfor
End




Function MakeConfigTermSingleValence()
	Wave/T Configuration, JW
	
	String c0, cf
	String s1, s2
	
	Duplicate/O/T Configuration SimpleConfig, ConfigTerm, SimpleConfigTerm
	
	Variable i
	for(i=0; i<numpnts(configuration); i+=1)
		if(stringmatch(configuration[i],""))
			configuration[i] = configuration[i-1]
		endif
		SplitString /E="([[:alpha:][:digit:]]+)\.([[:alpha:][:digit:]]+)" configuration[i] , s1, s2
		if(stringmatch(s1,""))
			SplitString /E="([[:digit:]]+[[:alpha:]])([[:digit:]]+)" configuration[i] , s1, s2
			SimpleConfig[i] = s1
			ConfigTerm[i] = s1 + " " + JW[i]
		else
			SimpleConfig[i] = s2
			ConfigTerm[i] = s2 + " " + JW[i]
		endif
	endfor
End



Function LevelLines()
	Wave Level
	Wave/T ConfigTerm, SimpleTerm, SimpleConfig, SimpleConfigTerm
	
	Variable pnts = numpnts(Level)
	String nameString, nLevel, LLevel, spinConfig, JLevel
	
	Make/O/N=(pnts)/WAVE waveOfConfigTermWaves
	
	Display/K=1/W=(30,70,800,600)
	
	Variable xpos = .5
	
	Variable spinpos, spinTripletExists=0
	
	Variable i
	For(i=0; i<pnts; i+=1)
		Make/O/N=2 $SimpleConfigTerm[i]
		Wave configtermWave = $SimpleConfigTerm[i]
		configtermWave = Level[i]
		waveOfConfigTermWaves[i] = configtermWave
		
		SplitString /E="([[:digit:]]+)([[:alpha:]])([[:digit:]]+)" SimpleTerm[i] , spinConfig, LLevel, JLevel

		AppendToGraph configtermWave
		
		if(str2num(spinConfig)==1)
			spinpos = 0
		elseif(str2num(spinConfig) == 3)
			spinpos = 4
			spinTripletExists =1
		endif
		
		if(stringmatch(LLevel,"s"))
			SetScale/p x 0+spinpos,1, configtermWave
			SetDrawEnv xcoord= bottom,ycoord= left; DrawText 0+spinpos+xpos,Level[i], SimpleConfigTerm[i]
		elseif(stringmatch(LLevel,"p"))
			SetScale/p x 1+spinpos,1, configtermWave
			SetDrawEnv xcoord= bottom,ycoord= left; DrawText 1+spinpos+xpos,Level[i], SimpleConfigTerm[i]
		elseif(stringmatch(LLevel,"d"))
			SetScale/p x 2+spinpos,1, configtermWave	
			SetDrawEnv xcoord= bottom,ycoord= left; DrawText 2+spinpos+xpos,Level[i], SimpleConfigTerm[i]
		elseif(stringmatch(LLevel,"f"))
			SetScale/p x 3+spinpos,1, configtermWave	
			SetDrawEnv xcoord= bottom,ycoord= left; DrawText 3+spinpos+xpos,Level[i], SimpleConfigTerm[i]
		elseif(stringmatch(SimpleTerm[i],"*Limit*"))
			SetScale/p x 0,8, configtermWave	
			SetDrawEnv xcoord= bottom,ycoord= left; DrawText xpos,Level[i], SimpleConfigTerm[i]
			ModifyGraph lsize($nameofwave(configtermWave))=2,rgb($nameofwave(configtermWave))=(0,0,0)
		endif	
		
	EndFor
	
	Label left "Level energy (eV)"
	Make/O/N=8 tickxPos =0.5+p
	Make/O/T tickLabels = {"S", "P", "D","F","S", "P", "D","F"}
	ModifyGraph userticks(bottom)={tickxPos,tickLabels}
	if(SpinTripletExists)
		SetAxis bottom 0,*
	else
		SetAxis bottom 0,4
	endif
	SetAxis left -0.1, *
	
	ShowInfo
	
	DoWindow/C/T kwTopWin, "energy levels"	
	
	setLambdaWindowHook()
End







Function LevelLinesSingleValence()
	Wave Level
	Wave/T ConfigTerm, SimpleTerm
	
	Variable pnts = numpnts(Level)
	String nameString, nLevel, LLevel
	
	Make/O/N=(pnts)/WAVE waveOfConfigTermWaves
	
	Display/K=1/W=(30,70,600,400)
	
	Variable xpos = .5
	
	Variable i
	For(i=0; i<pnts; i+=1)
		Make/O/N=2 $configterm[i]
		Wave configtermWave = $configterm[i]
		configtermWave = Level[i]
		waveOfConfigTermWaves[i] = configtermWave
		
		SplitString /E="([[:digit:]]+)([[:alpha:]])" ConfigTerm[i] , nLevel, LLevel
		
		AppendToGraph configtermWave
		
		if(stringmatch(LLevel,"s"))
			SetDrawEnv xcoord= bottom,ycoord= left; DrawText 0+xpos,Level[i], configterm[i]
		elseif(stringmatch(LLevel,"p"))
			SetScale/p x 1,1, configtermWave
			SetDrawEnv xcoord= bottom,ycoord= left; DrawText 1+xpos,Level[i], configterm[i]
		elseif(stringmatch(LLevel,"d"))
			SetScale/p x 2,1, configtermWave	
			SetDrawEnv xcoord= bottom,ycoord= left; DrawText 2+xpos,Level[i], configterm[i]
		elseif(stringmatch(LLevel,"f"))
			SetScale/p x 3,1, configtermWave	
			SetDrawEnv xcoord= bottom,ycoord= left; DrawText 3+xpos,Level[i], configterm[i]
		endif	
		
	EndFor
	
	Label left "Level energy (eV)"
	Make/O tickxPos = {0.5,1.5,2.5,3.5}
	Make/O/T tickLabels = {"S", "P", "D","F"}
	ModifyGraph userticks(bottom)={tickxPos,tickLabels}
	SetAxis bottom 0,4.25
	SetAxis left -0.1, 4
End





Function setLambdaWindowHook()
	GetWindow kwTopwin activeSW
	String activeSubwindow = S_value
	SetWindow $activeSubwindow, hook(CursorLambdaHook) = CursorLambdaHookFunc
	
	String/G levelA, levelB
	Variable/G energyA, energyB, energyDifference, Lambda
End


Function CursorLambdaHookFunc(s)
	Struct WMWinHookStruct &s
	
	SVAR levelA, levelB
	NVAR energyA, energyB, energyDifference, Lambda
	
	switch(s.eventCode)
		case 7:
		//print s.winName
		String movedCursor = s.cursorName
		String waveNameWithCursor
		SplitString /E="'([[:alpha:][:digit:]\s/\.]+)'" s.traceName , waveNameWithCursor	// traceName might have ' ' around it - need to get rid of them
		Wave newCursorWave = $waveNameWithCursor
		
		if(stringmatch(movedCursor, "A"))
			energyA = newCursorWave
			levelA = waveNameWithCursor
		elseif(stringmatch(movedCursor,"B"))
			energyB = newCursorWave
			levelB = waveNameWithCursor
		endif
		
		Variable/G energyDifference = abs(energyA-energyB)
		Variable/G Lambda = h*c/(energyDifference*echarge)*1e9
		
		String tboxStr
		if(energyB > energyA)
			sprintf tboxStr "%s - %s\rTransition wavelength = %.2f nm", levelA, levelB, Lambda
			printf "%s energy = %.3g eV\r", levelA, energyA
			printf "%s energy = %.3g eV\r", levelB, energyB
		else
			sprintf tboxStr "%s - %s\rTransition wavelength = %.2f nm", levelB, levelA, Lambda
			printf "%s energy = %.3g eV\r", levelB, energyB
			printf "%s energy = %.3g eV\r", levelA, energyA
		endif
		TextBox/C/N=lambdaBox/F=0/A=MC tboxStr
		TextBox/C/N=lambdaBox/X=18.00/Y=-30.00
		TextBox/C/N=lambdaBox/F=2
		

		printf "Energy difference = %.3f eV\r", energyDifference
		printf "Transition wavelength = %.2f nm\r\r", Lambda
	endswitch
End






// ==================================================================

Function/T TrimWS(str)
	// TrimWhiteSpace (code from Jon Tischler)
	String str
	return TrimWSL(TrimWSR(str))
End

// ==================================================================

Function/T TrimWSL(str)
	// TrimWhiteSpaceLeft (code from Jon Tischler)
	String str
	Variable i, N=strlen(str)
	for (i=0;char2num(str[i])<=32 && i<N;i+=1) // find first non-white space
	endfor
	return str[i,Inf]
End

// ==================================================================

Function/T TrimWSR(str)
	// TrimWhiteSpaceRight (code from Jon Tischler)
	String str
	Variable i
	for (i=strlen(str)-1; char2num(str[i])<=32 && i>=0; i-=1) // find last non-white space
	endfor
	return str[0,i]
End
