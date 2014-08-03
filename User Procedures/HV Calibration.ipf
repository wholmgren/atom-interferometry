#pragma rtGlobals=1		// Use modern global access method.
#pragma version = 1.1

Function GetHVKeithley(Vmon, VmonZero)
	Variable Vmon, VmonZero		//in Volts
	
	// Fit to y = mx+b
	// from calibration data "m" on 6/9/09 with VmonZeroBase=.013696
	Variable m = 6.01085		// kV/V
	Variable b = 0.00128738	// kV
	
	Variable HV = m*Vmon+b	// in kV
	
	// Fit HVatVmon2 vs VmonZero-VmonZeroBase to a quadratic polynomial
	// fit coefficients:
	Variable c0 = 12.0233
	Variable c1 = 34.5716
	Variable c2 = 57683.9
	
	Variable VmonZeroBase=.013696
	Variable DiffZero = VmonZero-VmonZeroBase
	Variable HVatVmon2Correct = c0+c1*DiffZero+c2*DiffZero^2
	Variable HVatVmon2wrong = m*2+b
	
	Variable HVatVmon2ratio = HVatVmon2Correct/HVatVmon2wrong
	
	Variable TrueHV = HV*HVatVmon2ratio
	
	//print HVatVmon2Correct

	//print HV
	
	//print "HV = ", TrueHV, "kV"
	
	return TrueHV
End


Function GetHVWavetek(Vmon)
	Variable Vmon		//in Volts
	
	// Fit keithleycalib vs wavetekcalib to y = mx+b
	// from calibration data on 6/12/09 with VmonZeroBase=.013687
	Variable m =     1.00225		// V/V
	Variable b =       0.000224631	// V

	if(0)
	wave w_coef
	 m = w_coef[1]	
	 b =  w_coef[0]	
	endif
	
	Variable VmonKeithley = m*Vmon+b
	Print VmonKeithley
	
	Variable TrueHV = GetHVKeithley(VmonKeithley, .013683)
	
	return TrueHV
End