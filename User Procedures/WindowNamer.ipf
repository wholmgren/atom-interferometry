#pragma rtGlobals=1		// Use modern global access method.

Function WindowNamer(baseName)
	String baseName
	
	String WindowName = baseName
	
	sprintf WindowName, "%.29s", WindowName		// ensures that the name is less than 31 characters
	
	WindowName = "WN"+WindowName	// ensures that the name is valid
	WindowName = ReplaceString(" ", WindowName, "") // ensures that the name is valid
	String WindowName29 = WindowName
	WindowName+="_1"
	//print WindowName
	
	//print WindowName
	// WindowName is something Igor uses internally and we don't usually see
	// baseName+" "+num2str(k) is title of the graph that we see

	Variable k=1
	do
		if(StringMatch(WinList(WindowName,";",""),""))
			break
		else
			k+=1
			WindowName = WindowName29+"_"+num2str(k)
		endif
	while(1)
	
	if(k==1)
		DoWindow/C/T $WindowName, baseName	
	else
		DoWindow/C/T $WindowName, baseName + " " + num2str(k)
	endif
	
End