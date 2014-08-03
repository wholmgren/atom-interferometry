#pragma rtGlobals=1		// Use modern global access method.

#include ":NIST ASD Grapher"


//	by Will Holmgren, 2/27/2012
//	Cronin group, Dept. of Physics, University of Arizona
//	holmgren@email.arizona.edu
//
//	Loads the output from the NIST, processes it, and displays a graph of the energy levels. Put cursors on two levels to calculate the transition wavelength
//
//	Instructions for generating the NIST ASD Lines data:
//	1. go to: 	http://physics.nist.gov/PhysRefData/ASD/lines_form.html
// 	2. enter atom and ionization state (e.g. k i for neutral potassium)
//	3. enter lower/upper wavelength limits
//	4. set "Units" to nm
//	4. set "Format ouput" to ASCII
//	3. set "Energy Level Units" to eV
//	4. set "Maximum lower energy" to 0 for transitions from ground state only
//	5. deselect Lande-g and Leading percentages check boxes under "Level information"
//	5. deselect Bibliographic information
//	5. select include fik
//	5. select "Wavelengths in: Vacuum (all wavelengths)"
//	5. deselect "relative intensity"
//	6. select, copy, and paste the output into a text file. Start with the column names line (	| Configuration | Term | J | Level |	). 
//		Do not include the line of "---" at the very top. End at the end of the next ionization state line (	| k ii | Limit | -- | 4.34 |	 ).
//	7. save the plain text file
//	8. set the fileToLoad path in the LoadASD() function
//	9. run the following command: goASD()
//	10. put cursors on the levels to calculate the associated transition wavelength.

// Observed, Ritz, Aki, fik, Acc, EiEk, Configurations, Terms, JiJk, gigk, Type


Static Constant h =  6.62607e-34
Static Constant hbar = 1.054571726e-34
Static Constant eCharge = 1.60217646e-19
Static Constant c = 299792458
Static Constant eMass = 9.10938188e-31
Static Constant hartree = 4.35974434e-18
Static Constant AUtoSI = 1.64878e-41

Function goASDLines()
	LoadASDLines()
	CleanupASDLines()
	MakeConfigTermLines()
	DisplayLines()
	CalcAlphaASD()
	
	AppendToGraph/R alphaASD
End




Function LoadASDLines()
	String fileToLoad = "Macintosh HD:Users:holmgren:Documents:School:UA:Cronin labs:NIST Levels:rblines.txt"
	String columnFormatting = "C=5,F=0; C=6,F=-2;"
	LoadWave/A/O/J/D/W/K=0/V={",|"," $",0,0}/B=columnFormatting/L={1,6,0,0,0}  fileToLoad
	
	Wave/T Ei___________Ek, Configurations, Terms, Ji___Jk, gi___gk, Type
	Duplicate/O Acc_ Acc
	Duplicate/O/T Ei___________Ek, EiEk
	Duplicate/O/T Ji___Jk, JiJk
	Duplicate/O/T gi___gk, gigk
	
	EiEk = TrimWS(EiEk[p])
	Configurations = TrimWS(Configurations[p])
	Terms = TrimWS(Terms[p])
	JiJk = TrimWS(JiJk[p])
	gigk = TrimWS(gigk[p])
	Type = TrimWS(Type[p])
	
	Edit/W=(50,100,1100,400)/K=1 Observed, Ritz, Aki, fik, Acc, EiEk, Configurations, Terms, JiJk, gigk, Type
End




Function CleanupASDLines()
	Wave/T EiEk, Configurations, Terms, JiJk, gigk, Type
	Wave Observed, Ritz, Aki, fik, Acc
	
	Variable i,j
	Variable numNaNs = 0
	For(j=0;j<2;j+=1)
	For(i=0; i<numpnts(Observed); i+=1)
		if(numtype(Observed[i])==2)
			DeletePoints i, 1, Observed, Ritz, Aki, fik, Acc, EiEk, Configurations, Terms, JiJk, gigk, Type
			numNans+=1
		endif
	EndFor
	EndFor
End




Function MakeConfigTermLines()
	Wave/T Configurations, JiJk, Terms, gigk
	
	String c0, cf
	String s1, s2, s3, s4
	
	Duplicate/O/T Configurations ConfigGround,ConfigExcited, TermGround, TermExcited, ConfigTermTrans
	Make/O/D/N=(numpnts(configurations)) gGround, gExcited
	
	Variable i
	for(i=0; i<numpnts(configurations); i+=1)
		if(stringmatch(configurations[i],""))
			configurations[i] = configurations[i-1]
		endif
		SplitString /E="(\\S+) - (\\S+)" configurations[i] , s1, s2
		ConfigGround[i] = s1
		ConfigExcited[i] = s2
	endfor
	
	for(i=0; i<numpnts(Terms); i+=1)
		if(stringmatch(Terms[i],""))
			Terms[i] = Terms[i-1]
		endif
		SplitString /E="([[:alpha:][:digit:]]+)\\s+\-\\s+([[:alpha:][:digit:]]+)" Terms[i] , s1, s2
		SplitString /E="(\\S+)\\s+\-\\s+(\\S+)" JiJk[i], s3, s4
		TermGround[i] = s1+s3
		TermExcited [i] = s2+s4
	endfor
	
	
	for(i=0; i<numpnts(Terms); i+=1)
		SplitString /E="([[:alpha:][:digit:]]+)\\s+\-\\s+([[:alpha:][:digit:]]+)" gigk[i] , s1, s2
		gGround[i] = str2num(s1)
		gExcited [i] = str2num(s2)
	endfor
	
	
	for(i=0; i<numpnts(Terms); i+=1)
		ConfigTermTrans[i] = ConfigGround[i]+" "+ TermGround[i] + " - " + ConfigExcited[i]+" "+ TermExcited[i] 
	endfor
	
	AppendToTable ConfigGround,ConfigExcited, TermGround, TermExcited, ConfigTermTrans, gGround, gExcited
End




Function DisplayLines()
	Wave/T EiEk, Configurations, Terms, JiJk, gigk, Type, ConfigTermTrans
	Wave Observed, Ritz, Aki, fik, Acc
	
	Variable pnts = numpnts(observed)
	Variable i
	
	Duplicate/O Observed BestWavelength
	BestWavelength = numtype(Observed[p])==0 ? Observed : Ritz

if(0)
	Display/K=1/W=(100,100,800,500) Aki vs BestWavelength
	ModifyGraph log(left)=1
	SetAxis left 1,*
	ModifyGraph mode=1
	ModifyGraph zColor(Aki)={BestWavelength,355,830,SpectrumBlack,0}
	Label bottom "Wavelength (nm)"
	Label left "Aik (Hz)"
	ModifyGraph lsize(Aki)=3
	
	For(i=0; i<pnts; i+=1)
		SetDrawEnv xcoord= bottom,ycoord= left; DrawText Observed[i],Aki[i], ConfigTermTrans[i]
	EndFor
else
	//Display/K=1/W=(100,100,900,500) fik vs BestWavelength
	AppendToGraph fik vs BestWavelength
	ModifyGraph log(left)=1
	//SetAxis left 1,*
	ModifyGraph mode(fik)=1
	ModifyGraph zColor(fik)={BestWavelength,355,830,SpectrumBlack,0}
	Label bottom "Wavelength (nm)"
	Label left "fik"
	ModifyGraph lsize(fik)=3

	For(i=0; i<pnts; i+=1)
		SetDrawEnv xcoord= bottom,ycoord= left; DrawText Observed[i],fik[i], ConfigTermTrans[i]
	EndFor
endif
	

End



Function calcAlphaASD()
	Wave bestWavelength, Observed, Aki, fik, gExcited, gGround
	Wave/Z fikNew
	
	Make/O/D/N=500000 alphaASD
	
	Variable minWave = 300//wavemin(bestWavelength)*.98
	Variable maxWave = 900//wavemax(bestWavelength)*1.02
	
	SetScale/i x minWave, maxWave, alphaASD
	
	Variable AkiFiltered, fikFiltered, thisPrefactor, wavelengthCalc, alphaASDstatic
	alphaASD=0
	Variable i 
	For(i=0; i<numpnts(bestWavelength); i+=1)
		AkiFiltered = Aki[i]
		fikFiltered = fik[i]
		//fikFiltered = fikNew[i]
		wavelengthCalc = bestWavelength[i]
		if(numtype(AkiFiltered)!=0 || numtype(fikFiltered)!=0)
			AkiFiltered=0
			fikFiltered=0
		endif
		thisPrefactor=  fikFiltered
		alphaASD += thisPrefactor / ( (h*c/wavelengthCalc/hartree)^2 - (2*pi*c/x/(hartree/hbar))^2 )*1e-18	// 1e-18 is for converting nm^2 to m^2
		alphaASDstatic += thisPrefactor / ( (h*c/wavelengthCalc/hartree)^2 )*1e-18
//		thisPrefactor = gExcited[i]/gGround[i] * eps0*eMass*c/(2*pi*eCharge^2)* (wavelengthCalc)^2 * AkiFiltered//*4*pi*eps0
//		alphaASD += thisPrefactor / ( (1/wavelengthCalc)^2 - (1/x)^2 ) * 1/(2*pi*c)^2 //* (1/bestWavelength[i])
//		alphaASDstatic += thisPrefactor / ( (1/wavelengthCalc)^2) * 1/(2*pi*c)^2 
	EndFor
	
	Variable alphaCore = 5.8
	alphaASD+=alphaCore
	alphaASDstatic+=alphaCore
	
	print "alphaASDstatic a.u = ", alphaASDstatic
	
	Variable staticFilterMultiplier = 1000
	alphaASD = abs(alphaASD) > staticFilterMultiplier*alphaASDstatic ? NaN : alphaASD
	
	Duplicate/O alphaASD alphaASDder
	Differentiate alphaASDder
	
	//print "alphaASDstatic a.u.= ", alphaASDstatic/AUtoSI
	
	//Variable overallPrefactor = echarge^2/eMass*1/3/(h*c)
	//Variable overallPrefactor = eps0*eMass*c/(2*pi*eCharge^2)
	
	//alphaASD*=overallPrefactor
	
//	Wavestats/Q
	
//	alphaASD = alphaASD > V_max ? NaN : alphaASD
End





//Function calcAlphaASDomega()
//	Wave bestWavelength, Observed, Aki, fik, gExcited, gGround, bestOmega
//	
//	Make/O/D/N=1000000 alphaASDomega
//	
//	Variable minWave = wavemin(bestWavelength)*.98
//	Variable maxWave = wavemax(bestWavelength)*1.02
//
//	Variable minOmega = 2*pi*c/maxWave*1e9
//	Variable maxOmega = 2*pi*c/minWave*1e9
//	
//	SetScale/i x minOmega, maxOmega, alphaASDomega
//	
//	Variable AkiFiltered
//	alphaASDomega=0
//	Variable i 
//	For(i=0; i<numpnts(bestOmega); i+=1)
//		AkiFiltered = Aki[i]
//		if(numtype(AkiFiltered)!=0)
//			AkiFiltered=0
//		endif
//		alphaASDomega += 1/3*gExcited[i]/gGround[i]*AkiFiltered  *  (bestOmega[i])/( bestOmega[i]^2 - x^2 )// * (1/bestWavelength[i])
//	EndFor
//	
//	
//	
//	//alphaASD = abs(alphaASD) > 1e27 ? NaN : alphaASD
//End