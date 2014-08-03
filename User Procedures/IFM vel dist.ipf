#pragma rtGlobals=1		// Use modern global access method.

#include ":PhysicalConstants"


function go()
	variable/g vint, vint2
	make/o/n=6 entire_vavg, detected_vavg, xvratio, vavg_shift
	xvratio = {5, 8, 10,15,20,30}
	variable i=0
	do
		vdist(3000, xvratio[i]); entire_vavg[i] = vint2; detected_vavg[i] = vint; vavg_shift[i] = (vint - vint2) / vint * 100
		i+=1
	while(i<=6)
end


function vdist(vavg, vratio)
	variable vratio
	variable vavg  // m/s
	//variable vratio = 8

	variable sigmav = vavg/vratio
	make/o/n=100 vel, pv, pvd
	vel = x* 2* vavg / numpnts(vel)
	pv = vel^3 * exp( - (vel - vavg)^2 /  (2*sigmav^2)  )
	wavestats/q pv ; pv /= v_max

	make/o/n=2000 vflux, dflux
	vflux= 0
	variable xo = 0
	variable detwidth = 60
	variable i=1
	do
		beam(vel[i])
		Wave iflux
		vflux += iflux //* pv[i]
		dflux = iflux*(x>1000+xo-detwidth/2)*(x<1000+xo+detwidth/2)
		pvd[i] = sum(dflux)*pv[i]
		i += 1
	while(i<numpnts(vel))

	wavestats/q pvd ; pvd /= v_max

	duplicate/o pvd pvdx
	pvdx = pvd * vel
	wavestats/q pvdx
	variable/g vint = v_sum
	wavestats/q pvd
	vint /= v_sum
	print "avg vel equals ",vint

	duplicate/o pv pvx
	pvx = pv * vel
	wavestats/q pvx
	variable/g vint2 = v_sum
	wavestats/q pv
	vint2 /= v_sum
	print "avg vel of entire beam equals ",vint2
	variable FractionalShift = (vint2 - vint) / vint2 * 100
	print "fractional vshift is " , (vint2 - vint) / vint2 * 100 ,  " %"
	MatrixOP/O pvdiff = pv - pvd
	
	return FractionalShift
end




function beam(vel)
	variable vel
	//silent -1; pauseupdate
	variable s1w = 36 // microns
	variable s2w = 47
	variable s1s2 = 0.889e6
	variable g1g2 = .93927e6
	variable g3det = 0.492e6
	variable g1det = g1g2*2+g3det
	Variable OneGTwoSdist = 0.127E6
	variable s2det = g1det+OneGTwoSdist
	variable open_fraction = .37
	
	// calculate diffraction amplitudes
	variable e0, e1, e2
	e0 = open_fraction* sinc(0*open_fraction*pi)
	e1 = open_fraction* sinc(1*open_fraction*pi)
	e2 = open_fraction* sinc(2*open_fraction*pi)
	
	// calculate diffraction intensities
	//print e0^2 / e0^2, e1^2/ e0^2, e2^2/ e0^2 

	// calculate IFM intensities
	variable ifmw1, ifmw2 ,ifmw3, ifmw4, ifmw5, ifmw6, ifmw7, ifmw8
	ifmw1 = e1*e1*e2  *  e0*e1*e1
	ifmw2 = e1*e1*e1  *  e0*e1*e0
	ifmw3 = e1*e1*e0  *  e0*e1*e1
	ifmw4 = e1*e1*e1  *  e0*e1*e2
	ifmw5 = ifmw4
	ifmw6 = ifmw3
	ifmw7 = ifmw2
	ifmw8 = ifmw1

	variable ifmw21, ifmw22, ifmw23, ifmw24, ifmw25, ifmw26, ifmw27, ifmw28
	ifmw21 = e2*e1*e2  *  e1*e1*e1   // second order from 1g weights
	ifmw22 = e2*e1*e1  *  e1*e1*e0
	ifmw23 = e2*e1*e0  *  e1*e1*e1
	ifmw24 = e2*e1*e1  *  e1*e1*e2
	ifmw25 = ifmw24
	ifmw26 = ifmw23
	ifmw27 = ifmw22
	ifmw28 = ifmw21
	
	//print "IFM weights 1-4: ", ifmw1, ifmw2 ,ifmw3, ifmw4
	//print "IFM2 weights 1-4: ",  ifmw21, ifmw22, ifmw23, ifmw24
	//print "IFM ratios 1-4: ",  ifmw21/ifmw1, ifmw22/ifmw2, ifmw23/ifmw3, ifmw24/ifmw4

	if(0)		// set to 1 to turn off 2nd order ifm
		ifmw21 = 0
		ifmw22 = 0
		ifmw23 = 0
		ifmw24 = 0
		ifmw25 = 0
		ifmw26 = 0
		ifmw27 = 0
		ifmw28 = 0
	endif


	variable x2 = 163  *  1000 / vel	// location of first order IFM at 3g
	variable x3 = g3det / g1g2 * x2		// displacement of first order diffraction from 3g at detector plane

	// calculate location of first order IFM outputs at detector plane
	variable ifmx1 = -x2 -  2*x3
	variable ifmx2 = -x2 -  1*x3
	variable ifmx3 = -x2  +0*x3   // used to be an error here:  +1*x3
	variable ifmx4 = -x2 + 1*x3   // used to be an error here:  +2*x3
	variable ifmx5 = x2 -  1*x3
	variable ifmx6 = x2 -  0*x3
	variable ifmx7 = x2 + 1*x3
	variable ifmx8 = x2 + 2*x3

	// calculate location of second order IFM outputs at detector plane
	variable ifmx21 = -3*x2 -  3*x3  // 2nd order from 1g positon
	variable ifmx22 = -3*x2 -  2*x3
	variable ifmx23 = -3*x2 - 1*x3
	variable ifmx24 = -3*x2 - 0*x3
	variable ifmx25 = 3*x2 +  0*x3
	variable ifmx26 = 3*x2 +  1*x3
	variable ifmx27 = 3*x2 + 2*x3
	variable ifmx28 = 3*x2 + 3*x3

	// points defining the trapazoidal beam
	//variable g1int = .80257e6; s2det=G1int		// decomment to find t1,t2 at interaction region
	variable t1 = 0.5*(s1w+(s1w+s2w)*s2det/s1s2)//s2w + (s2w +s1w) * g1det / s1s2
	variable t2 = 0.5*abs((s2w+(s2w-s1w)*s2det/s1s2))//s2w + (s2w -s1w) * g1det / s1s2
	variable t3 = -t2
	variable t4 = -t1

	make/o/n=2000 flux, xdet, iflux
	xdet = x - 1000
	
	// make trapazoidal beam shape
	flux = min((xdet - t1) / (t2 - t1), 1)
	flux = min(flux, (xdet - t4) / (t3 - t4) ) 
	flux = max(0,flux)

	// calculate mulitple IFM flux at det plane
	iflux = 0
	iflux += flux(x-ifmx1)  * ifmw1
	iflux += flux(x-ifmx2)  * ifmw2
	iflux += flux(x-ifmx3)  * ifmw3
	iflux += flux(x-ifmx4)  * ifmw4
	iflux += flux(x-ifmx5)  * ifmw5
	iflux += flux(x-ifmx6)  * ifmw6
	iflux += flux(x-ifmx7)  * ifmw7
	iflux += flux(x-ifmx8)  * ifmw8

	iflux += flux(x-ifmx21)  * ifmw21
	iflux += flux(x-ifmx22)  * ifmw22
	iflux += flux(x-ifmx23)  * ifmw23
	iflux += flux(x-ifmx24)  * ifmw24
	iflux += flux(x-ifmx25)  * ifmw25
	iflux += flux(x-ifmx26)  * ifmw26
	iflux += flux(x-ifmx27)  * ifmw27
	iflux += flux(x-ifmx28)  * ifmw28

	variable normalize = sum(iflux)
	iflux /= normalize

end


   
Window Graph6() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(343.5,58.25,738,372.5) vavg_shift vs xvratio
	AppendToGraph vas50i4 vs xvratio
	AppendToGraph vas50i4ofp5 vs xvratio
	AppendToGraph vavg_shift_wrong vs xvratio
	AppendToGraph vas50i2 vs xvratio
	AppendToGraph vas50i4off100 vs xvratio
	AppendToGraph vas50i4off200 vs xvratio
	AppendToGraph vas50i4off300 vs xvratio
	ModifyGraph lSize(vavg_shift)=2,lSize(vas50i4)=2,lSize(vas50i4ofp5)=2,lSize(vavg_shift_wrong)=2
	ModifyGraph lStyle(vavg_shift_wrong)=3
	ModifyGraph rgb(vas50i4)=(0,15872,65280),rgb(vas50i4ofp5)=(0,65280,65280),rgb(vavg_shift_wrong)=(0,0,0)
	ModifyGraph zero(left)=1
	ModifyGraph fSize=14
	Label left "velocity correction (%)"
	Label bottom "velocity ratio"
	SetAxis/E=1 left -0.5,4.2017617
	Legend/C/N=text0/J/A=MC/X=-7.85/Y=21.31 "\\s(vavg_shift) vavg_shift\r\\s(vas50i4) vas50i4\r\\s(vas50i4ofp5) vas50i4ofp5\r\\s(vavg_shift_wrong) vavg_shift_wrong"
	AppendText "\\s(vas50i2) vas50i2\r\\s(vas50i4off100) vas50i4off100\r\\s(vas50i4off200) vas50i4off200\r\\s(vas50i4off300) vas50i4off300"
EndMacro
