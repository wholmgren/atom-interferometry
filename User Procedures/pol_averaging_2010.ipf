#pragma rtGlobals=1		// Use modern global access method.

constant a0tocm = 0.148184
constant NaP = 24.11
constant NaPstat = 0.06
constant NaPsys = 0.06
constant NaPSysStatFracErr = 0.0035194

constant beamWidthCorrection = 0.04 //%
constant kIsoCorrection = 0.04//%
constant rbIsoCorrection = 0.02//%

////////////////////////////////////////////////////////////////////////////////////////////
//
//	A few useful functions for averaging polarizability data sets. Change the names of the waves to suit your purposes.
//
//	PolStat() - calculate the statistical error
//	PolRatios() - calculate the polarizability ratios and their errors
//	PolAvg() - calculate average polarizabilities, with and without reference to Ekstrom/Pritchard
//	PolAvgTrim() - same as PolAvg(), but with waves used for a trimmed mean.
//
////////////////////////////////////////////////////////////////////////////////////////////

Function PolStat()
	Wave NaCmol, KCmol, RbCmol
	Wave NaCmolStat, KCmolStat, RbCmolStat
	
	NVAR NaAvg, KAvg, RbAvg, NaStdErr, KStdErr, RbStdErr//, NaPercentError, KPercentError, RbPercentError
	
	Variable/G RbPercentError = RbStdErr/RbAvg*100
	Variable/G KPercentError = KStdErr/KAvg*100
	Variable/G NaPercentError = NaStdErr/NaAvg*100
	print "Rb= ", RbCmol[1000], " +- ", RbCmolStat[1000], " (", RbCmolStat[1000]/RbCmol[1000]*100, "%)"
	print "K= ", KCmol[1000], " +- ", KCmolStat[1000], " (", KCmolStat[1000]/KCmol[1000]*100, "%)"
	print "Na= ", NaCmol[1000], " +- ", NaCmolStat[1000], " (", NaCmolStat[1000]/NaCmol[1000]*100, "%)"
	
	String RbStr, KStr, NaStr, LegendString
	
	sprintf RbStr, "Rb = %4.2f ± %1.2f (%1.2f%)\r", RbAvg, RbStdErr, RbPercentError
	sprintf KStr, "K = %4.2f ± %1.2f (%1.2f%)\r", KAvg, KStdErr, KPercentError
	sprintf NaStr, "Na = %4.2f ± %1.2f (%1.2f%)", NaAvg, NaStdErr, NaPercentError
			
	LegendString = "\\s(RbCMol) "+RbStr+"\\s(KcMol) "+KStr+"\\s(NaCmol) "+NaStr
//	print LegendString
	
	Legend/C/N=PolAbsTBox/J "\\s(RbCMol) "+RbStr+"\\s(KcMol) "+KStr+"\\s(NaCmol) "+NaStr
End


Function PolRatios()
	Wave NaCmol, KCmol, RbCmol
	Wave NaCmolStat, KCmolStat, RbCmolStat
	
	Variable ratio, ratioErr, fracError
	
	ratio = RbCmol[1000]/NaCmol[1000]; fracError = sqrt((RbCmolStat[1000]/RbCmol)^2+(NaCmolStat[1000]/NaCmol)^2); ratioErr = fracError*ratio
	print "Rb/Na = ", ratio, " +- ", ratioErr, " (", fracError*100, "%)"
	
	ratio = KCmol[1000]/NaCmol[1000]; fracError = sqrt((KCmolStat[1000]/KCmol)^2+(NaCmolStat[1000]/NaCmol)^2); ratioErr = fracError*ratio
	print "K/Na = ", ratio, " +- ", ratioErr, " (", fracError*100, "%)"
	
	ratio = RbCmol[1000]/KCmol[1000]; fracError = sqrt((RbCmolStat[1000]/RbCmol)^2+(KCmolStat[1000]/KCmol)^2); ratioErr = fracError*ratio
	print "Rb/K = ", ratio, " +- ", ratioErr, " (", fracError*100, "%)"
	
//	print "Rb/Na = ", RbCmol[1000]/NaCmol[1000], " +- ", RbCmol[1000]/NaCmol[1000]*sqrt((RbCmolStat[1000]/RbCmol)^2+(NaCmolStat[1000]/NaCmol)^2)
//	print "K/Na = ", KCmol[1000]/NaCmol[1000], " +- ", KCmol[1000]/NaCmol[1000]*sqrt((KCmolStat[1000]/KCmol)^2+(NaCmolStat[1000]/NaCmol)^2)
//	print "Rb/K = ", RbCmol[1000]/KCmol[1000], " +- ", RbCmol[1000]/KCmol[1000]*sqrt((RbCmolStat[1000]/RbCmol)^2+(KCmolStat[1000]/KCmol)^2)
End



Function PolAvg()
	Wave NaCmol, KCmol, RbCmol
	Wave NaCmolStat, KCmolStat, RbCmolStat
	Wave NaCmolAvg, KCmolAvg, RbCmolAvg
	
	Variable points, stderr
	
	points = numpnts(NaCmol)
	Wavestats/Q/R=[0,points-2] NaCmol
	Variable/G NaAvg = V_avg
	NaCmolAvg = V_avg
	NaCmol[points] = V_avg
	Variable/G NaStdErr = V_sdev/sqrt(V_npnts)
	NaCmolStat[points] = NaStdErr
	
	points = numpnts(KCmol)
	Wavestats/Q/R=[0,points-2] KCmol
	Variable/G KAvg = V_avg
	KCmolAvg = V_avg
	KCmol[points] = V_avg
	Variable/G KStdErr = V_sdev/sqrt(V_npnts)
	KCmolStat[points] = KStdErr
	
	points = numpnts(RbCmol)
	Wavestats/Q/R=[0,points-2] RbCmol
	Variable/G RbAvg = V_avg
	RbCmolAvg = V_avg
	RbCmol[points] = V_avg
	Variable/G RbStdErr = V_sdev/sqrt(V_npnts)
	RbCmolStat[points] = RbStdErr
	
	NaAvg*=(1+beamWidthCorrection/100)
	KAvg*=(1+beamWidthCorrection/100)*(1+kIsoCorrection/100)
	RbAvg*=(1+beamWidthCorrection/100)*(1+rbIsoCorrection/100)
	
	Variable/G RbPercentError = RbStdErr/RbAvg*100
	Variable/G KPercentError = KStdErr/KAvg*100
	Variable/G NaPercentError = NaStdErr/NaAvg*100
	
	printf "Rb = %2.4f ± %1.4f (%1.3f%)\r", RbAvg, RbStdErr, RbPercentError
	printf "K = %2.4f ± %1.4f (%1.3f%)\r", KAvg, KStdErr, KPercentError
	printf "Na = %2.4f ± %1.4f (%1.3f%)\r\r", NaAvg, NaStdErr, NaPercentError
	
	String RbStr, KStr, NaStr, LegendString
	
	sprintf RbStr, "Rb = %4.2f ± %1.2f (%1.2f%)\r", RbAvg, RbStdErr, RbPercentError
	sprintf KStr, "K = %4.2f ± %1.2f (%1.2f%)\r", KAvg, KStdErr, KPercentError
	sprintf NaStr, "Na = %4.2f ± %1.2f (%1.2f%)", NaAvg, NaStdErr, NaPercentError
			
	LegendString = "\\s(RbCMol) "+RbStr+"\\s(KcMol) "+KStr+"\\s(NaCmol) "+NaStr
	
	Legend/C/N=PolAbsTBox/J "\\s(RbCMol) "+RbStr+"\\s(KcMol) "+KStr+"\\s(NaCmol) "+NaStr
	
	
	Variable ratio, ratioErr, fracError
	Variable RbP, RbPErr, RbPFracErr, KP, KPErr, KPFracErr, NaRatio = NaAvg/NaP
	
	ratio = RbAvg/NaAvg; fracError = sqrt((RbStdErr/RbAvg)^2+(NaStdErr/NaAvg)^2); ratioErr = fracError*ratio; fracError*=100
	sprintf RbStr, "Rb/Na = %2.3f ± %2.3f (%1.2f%)\r", ratio, ratioErr, fracError
	printf "Rb/Na = %2.4f ± %2.4f (%1.3f%)\r", ratio, ratioErr, fracError
	RbP = ratio*NaP; RbPFracErr = sqrt((fracError/100)^2+(NaPStat/NaP)^2+(NaPSys/NaP)^2); RbPErr = RbP*RbPFracErr; RbPFracErr*=100
	ratio = KAvg/NaAvg; fracError = sqrt((KStdErr/KAvg)^2+(NaStdErr/NaAvg)^2); ratioErr = fracError*ratio; fracError*=100
	sprintf KStr, "K/Na = %2.3f ± %2.3f (%1.2f%)\r", ratio, ratioErr, fracError
	printf "K/Na = %2.4f ± %2.4f (%1.3f%)\r", ratio, ratioErr, fracError
	KP = ratio*NaP; KPFracErr = sqrt((fracError/100)^2+(NaPStat/NaP)^2+(NaPSys/NaP)^2); KPErr = KP*KPFracErr; KPFracErr*=100
	ratio = RbAvg/KAvg; fracError = sqrt((RbStdErr/RbAvg)^2+(KStdErr/KAvg)^2); ratioErr = fracError*ratio; fracError*=100
	sprintf NaStr, "Rb/K = %2.3f ± %2.3f (%1.2f%)", ratio, ratioErr, fracError
	printf "Rb/K = %2.4f ± %2.4f (%1.3f%)\r\r", ratio, ratioErr, fracError
	
	TextBox/C/N=PolRatioTBox RbStr+KStr+NaStr
	
	printf "Rb_Pritchard = %2.4f ± %2.4f (%1.3f%)\r", RbP, RbPErr, RbPFracErr
	printf "K_Pritchard = %2.4f ± %2.4f (%1.3f%)\r", KP, KPErr, KPFracErr
	printf "NaP/NaH = %2.4f\r", NaRatio
	
	String RbPStr, KPStr, NaRatioStr
	sprintf RbPStr, "Rb_Pritchard = %2.2f ± %2.2f (%1.2f%)\r", RbP, RbPErr, RbPFracErr
	sprintf KPStr, "K_Pritchard = %2.2f ± %2.2f (%1.2f%)\r", KP, KPErr, KPFracErr
	sprintf NaRatioStr, "NaPrit./NaMeas = %2.4f", NaRatio
	TextBox/C/N=PolPritchardTBox RbPStr+KPStr+NaRatioStr
End


Function ScaleAxis(scaleFactor)
	Variable scaleFactor
	
	NVAR NaAvg, KAvg, RbAvg
	
	SetAxis rb RbAvg*(1-scaleFactor), RbAvg*(1+scaleFactor)
	SetAxis K KAvg*(1-scaleFactor), KAvg*(1+scaleFactor)
	SetAxis left NaAvg*(1-scaleFactor), NaAvg*(1+scaleFactor)
End
	
	
Function PolAvgTrim()	
	Wave naTrimMean, naTrimSem, kTrimMean, kTrimSem, rbTrimMean, rbTrimSem

	Variable naAvg = naTrimMean[0]
	Variable naStdErr= naTrimSem[0]
	Variable kAvg = kTrimMean[0]
	Variable kStdErr= kTrimSem[0]
	Variable rbAvg = rbTrimMean[0]
	Variable rbStdErr= rbTrimSem[0]

	NaAvg*=(1+beamWidthCorrection/100)
	KAvg*=(1+beamWidthCorrection/100)*(1+kIsoCorrection/100)
	RbAvg*=(1+beamWidthCorrection/100)*(1+rbIsoCorrection/100)
	
	Variable/G RbPercentError = RbStdErr/RbAvg*100
	Variable/G KPercentError = KStdErr/KAvg*100
	Variable/G NaPercentError = NaStdErr/NaAvg*100
	
	printf "Rb = %2.4f ± %1.4f (%1.3f%)\r", RbAvg, RbStdErr, RbPercentError
	printf "K = %2.4f ± %1.4f (%1.3f%)\r", KAvg, KStdErr, KPercentError
	printf "Na = %2.4f ± %1.4f (%1.3f%)\r\r", NaAvg, NaStdErr, NaPercentError
	
	String RbStr, KStr, NaStr, LegendString
	
	sprintf RbStr, "Rb = %4.2f ± %1.2f (%1.2f%)\r", RbAvg, RbStdErr, RbPercentError
	sprintf KStr, "K = %4.2f ± %1.2f (%1.2f%)\r", KAvg, KStdErr, KPercentError
	sprintf NaStr, "Na = %4.2f ± %1.2f (%1.2f%)", NaAvg, NaStdErr, NaPercentError
			
	LegendString = "\\s(RbCMol) "+RbStr+"\\s(KcMol) "+KStr+"\\s(NaCmol) "+NaStr
	
	Legend/C/N=PolAbsTBox/J "\\s(RbCMol) "+RbStr+"\\s(KcMol) "+KStr+"\\s(NaCmol) "+NaStr

	Variable ratio, ratioErr, fracError
	Variable RbP, RbPErr, RbPFracErr, KP, KPErr, KPFracErr, NaRatio = NaAvg/NaP
	
	ratio = RbAvg/NaAvg; fracError = sqrt((RbStdErr/RbAvg)^2+(NaStdErr/NaAvg)^2); ratioErr = fracError*ratio; fracError*=100
	sprintf RbStr, "Rb/Na = %2.3f ± %2.3f (%1.2f%)\r", ratio, ratioErr, fracError
	printf "Rb/Na = %2.4f ± %2.4f (%1.3f%)\r", ratio, ratioErr, fracError
	RbP = ratio*NaP; RbPFracErr = sqrt((fracError/100)^2+(NaPStat/NaP)^2+(NaPSys/NaP)^2); RbPErr = RbP*RbPFracErr; RbPFracErr*=100
	ratio = KAvg/NaAvg; fracError = sqrt((KStdErr/KAvg)^2+(NaStdErr/NaAvg)^2); ratioErr = fracError*ratio; fracError*=100
	sprintf KStr, "K/Na = %2.3f ± %2.3f (%1.2f%)\r", ratio, ratioErr, fracError
	printf "K/Na = %2.4f ± %2.4f (%1.3f%)\r", ratio, ratioErr, fracError
	KP = ratio*NaP; KPFracErr = sqrt((fracError/100)^2+(NaPStat/NaP)^2+(NaPSys/NaP)^2); KPErr = KP*KPFracErr; KPFracErr*=100
	ratio = RbAvg/KAvg; fracError = sqrt((RbStdErr/RbAvg)^2+(KStdErr/KAvg)^2); ratioErr = fracError*ratio; fracError*=100
	sprintf NaStr, "Rb/K = %2.3f ± %2.3f (%1.2f%)", ratio, ratioErr, fracError
	printf "Rb/K = %2.4f ± %2.4f (%1.3f%)\r\r", ratio, ratioErr, fracError
	
	TextBox/C/N=PolRatioTBox RbStr+KStr+NaStr
	
	printf "Rb_Pritchard = %2.4f ± %2.4f (%1.3f%)\r", RbP, RbPErr, RbPFracErr
	printf "K_Pritchard = %2.4f ± %2.4f (%1.3f%)\r", KP, KPErr, KPFracErr
	printf "NaP/NaH = %2.4f\r", NaRatio
	
	String RbPStr, KPStr, NaRatioStr
	sprintf RbPStr, "Rb_Pritchard = %2.2f ± %2.2f (%1.2f%)\r", RbP, RbPErr, RbPFracErr
	sprintf KPStr, "K_Pritchard = %2.2f ± %2.2f (%1.2f%)\r", KP, KPErr, KPFracErr
	sprintf NaRatioStr, "NaPrit./NaMeas = %2.4f", NaRatio
	TextBox/C/N=PolPritchardTBox RbPStr+KPStr+NaRatioStr
End