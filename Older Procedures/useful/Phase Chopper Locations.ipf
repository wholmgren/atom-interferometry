#pragma rtGlobals=1		// Use modern global access method.

#include ":Constants"



Function PhaseVsPosition()
	Variable vflow = 2800
	Variable vratio = 17
	Variable sigv = vflow/vratio
	Variable nsigs = 4
	Variable vmax = vflow*(1+nsigs/vratio)
	Variable vmin = vflow*(1-nsigs/vratio)
	
	Variable voltage = 4.6e3
	Variable position = .75 //[mm]
	position /=1000	// convert to m
	
	Variable Lc1 = 11*25.4/1000

	CalculateConstantsC1C2(voltage)
	NVAR b, PhasePrefactor

//	Variable posPnts = 40
//	Make/O/N=(posPnts) phaseshift, contrast
	
	Variable velPnts = 100
	Make/O/N=(velPnts) probvel, phase_v, realpart, imagpart, phase_v, contrast_v
	SetScale/i x vmin, vmax, probvel, phase_v, realpart, imagpart, phase_v, contrast_v
	probvel = x^3*exp(-0.5*vratio^2*(x/vflow-1)^2)
	Variable probnorm = area(probvel)
	probvel/=probnorm
	
	Variable polarizability = AlphaNaSI
	Variable mass = Na23mass*amu2kg
	
	Variable SeperationTimesV = Lc1*hbar*2*pi/mass/gratingPeriod
	
	phase_v = PhasePrefactor*polarizability/x  * b* (-1/(b^2-position^2) + 1/(b^2-(position+SeperationTimesV/x)^2))
	
	realpart = probvel*cos(phase_v)
	imagpart = probvel*sin(phase_v)
	
	Variable realpartsum = area(realpart)
	Variable imagpartsum = area(imagpart)
	
	Variable phaseshift = atan2(imagpartsum, realpartsum)
	Variable contrast = sqrt(imagpartsum^2+realpartsum^2)
	
	if(phaseshift<0)
		phaseshift+=2*pi
	endif
	
	print "phase shift = ", phaseshift
	print "contrast = ", contrast
End




Function FindVoltage()
	Variable vflow = 3000
	Variable vratio = 15
	Variable sigv = vflow/vratio
	Variable nsigs = 4
	Variable vmax = vflow*(1+nsigs/vratio)
	Variable vmin = vflow*(1-nsigs/vratio)
	
	Variable voltage = 5e3
	Variable position = .75 //[mm]
	position /=1000	// convert to m
	
	Variable Lc1 = 11*25.4/1000

	CalculateConstantsC1C2(voltage)
	NVAR b, PhasePrefactor

//	Variable posPnts = 40
//	Make/O/N=(posPnts) phaseshift, contrast
	
	Variable velPnts = 100
	Make/O/N=(velPnts) probvel, phase_v, realpart, imagpart, phase_v, contrast_v
	SetScale/i x vmin, vmax, probvel, phase_v, realpart, imagpart, phase_v, contrast_v
	probvel = x^3*exp(-vratio^2*(x/vflow-1)^2)
	Variable probnorm = area(probvel)
	probvel/=probnorm
	
	Variable polarizability = 24.1e-30*4*pi*eps0
	
	phase_v = PhasePrefactor*polarizability/x  * b* (1/(b^2-position^2) - 1/(b^2-(position+Lc1*hbar*2*pi/(Na23mass*amu2kg)/gratingPeriod/x)^2))
	
	realpart = probvel*cos(phase_v)
	imagpart = probvel*sin(phase_v)
	
	Variable realpartsum = area(realpart)
	Variable imagpartsum = area(imagpart)
	
	Variable phaseshift = atan2(imagpartsum, realpartsum)
	Variable contrast = sqrt(imagpartsum^2+realpartsum^2)
	
	print "phase shift = ", phaseshift
	print "contrast = ", contrast
End




Function CalculateConstantsC1C2(voltage)
	Variable voltage
	
	Variable a = .001 		// distance from wire to ground plane [m]
	Variable D = .00155	// diameter of wire [m]
	Variable R = D/2
	
	Variable/G b = a*sqrt(1+2*R/a)	
	Variable/G lambda = 2*pi*eps0*voltage*(ln((a+R+b)/(a+R-b)))^-1
	Variable/G PhasePrefactor = lambda^2/(pi*eps0^2*hbar)
End