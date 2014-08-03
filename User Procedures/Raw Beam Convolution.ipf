#pragma rtGlobals=1		// Use modern global access method.





Function TestRawBeam()
	Variable testpnts = 401
	Variable halfScanWidth = 400
	Variable resolution = 2*halfScanWidth/testpnts
	Make/D/O/N=(testpnts) testY, testX
	SetScale/P x -(testpnts-1)/2, 1, "um", testX, testY
	testX = x
	
	Make/D/O/N=1 w_coef={615, 50,40,0, 0}
	
	FitRawBeam(w_coef, testY, testX)
	
	Duplicate/O testY testYErr
	testYErr = sqrt(testY)
End



Function FitRawBeam(pw, yw, xw) : FitFunc
	Wave pw, yw, xw
	
	//pw[0] = amplitude
	//pw[1] = 1s width
	//pw[2] = 2s width
	//pw[3] = offset
	//pw[4] = background

	Variable detWidth = 70
	
	Variable s1 = pw[1]
	Variable s2 = pw[2]
	
	Variable L12 = .889E6   
	Variable L1d = 2.39736E6
	Variable L2d = L1d - L12
	Variable a = L2d/L1d
	
	Variable pp = (1/2)*abs(s2+(s2-s1)*a)
	Variable dd = (1/2)*(s2+(s2+s1)*a)
	//printf "amp = %d, s1= %d, p = %d, d = %d\r" pw[0], s1, pp, dd
	
	Variable I0 = pw[0]
	
	Variable dX = 1	//resolution
	
	Variable x0 = pw[3]//; xw -= x0
	
	Variable BeamPnts = 401
	Make/D/O/N=(BeamPnts) RawBeam
	SetScale/P x  -dX*((BeamPnts-1)/2),dX,"um", RawBeam
	
	RawBeam = I0*( (dd+x)/(dd-pp)*(-dd<=x && x<-pp) + 1*(-pp<=x && x<pp) + (dd-x)/(dd-pp)*(pp<=x && x<=dd) )
	//RawBeam = I0*( (dd+x-x0)/(dd-pp)*(-dd<=x-x0 && x-x0<-pp) + 1*(-pp<=x-x0 && x-x0<pp) + (dd-x-x0)/(dd-pp)*(pp<=x-x0 && x-x0<=dd) )
	
	Variable detpnts = 101
	Make/D/O/N=(detpnts) Detector
	SetScale/P x -dX*((DetPnts-1)/2),dX, "um", Detector
	
	Detector = (x >= -detWidth/2 && x<= detWidth/2)
	//Detector = (x-x0 >= -detWidth/2 && x-x0<= detWidth/2)
	Variable sumdet = sum(Detector, -inf, inf)
	Detector /= sumdet
	
	Duplicate/O RawBeam RawBeamConvolved
	
	Convolve/A Detector, RawBeamConvolved
	
	Variable background = pw[4]
	yw = RawBeamConvolved(xw[p]-x0)+background
End



Function GoFitRawBeam(ywave, xwave, errwave)
	Wave ywave, xwave, errwave
	
	Make/O/T ParamsRawBeam={"Amplitude", "1s width", "2s width", "x offset", "background"}
	Make/D/O/N=1 W_coef={550, 50,2,0, 50}
	Make/T/O/N=1 Constraints = {"K0>1", "K0<1000", "K1>5", "K1<100"}
	
	FuncFit/L=600/H="00000" FitRawBeam W_coef  ywave /X=xwave /R /D /W=errwave /I=1 ///C=Constraints 
	
	AppendResults(ParamsRawBeam)
End



///////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function TestRawBeamGauss()
	Variable testpnts = 401
	Variable halfScanWidth = 400
	Variable resolution = 2*halfScanWidth/testpnts
	Make/D/O/N=(testpnts) testY, testX
	SetScale/P x -(testpnts-1)/2, 1, "um", testX, testY
	testX = x
	
	Make/D/O/N=1 initGuesses={615, 50, 40, 0, 0, 10, 50}
	Wave w_coef; initGuesses = w_coef
	initGuesses[1]=37; initguesses[2]=42
	
	FitRawBeamGauss(initGuesses, testY, testX)
	
	Duplicate/O testY testYErr
	testYErr = sqrt(testY)
End



Function GoFitRawBeamGauss(ywave, xwave, errwave)
	Wave ywave, xwave, errwave
	
	Make/O/T ParamsRawBeamGauss={"Amplitude", "1s width", "2s width", "x offset", "background", "gauss amp ratio", "gauss width"}
	Make/D/O/N=1 W_coef={550, 30, 45, 0, 1, .002, 200}
	//Make/T/O/N=1 Constraints = {"K5>0.0001", "K6>0.0001"}
	Make/T/O/N=1 Constraints = {"K4>0.0001", "K5>0.0001", "K6>0.0001"}
	Duplicate/O  w_coef, eps
	eps = 1e-6; eps[1]=1e0; eps[2]=1e0
	
	FuncFit/M=2/L=600/H="0000000" FitRawBeamGauss W_coef  ywave /X=xwave /R /D /W=errwave /I=1 /E=eps /C=Constraints
	
	Wave M_Covar
	Duplicate/O M_Covar, CorMat	 
	CorMat = M_Covar[p][q]/sqrt(M_Covar[p][p]*M_Covar[q][q])
	
	//NVAR numHolds
	Variable reducedChiSqrd = V_chisq/(V_npnts-(numpnts(W_coef)))
	printf "Reduced Chi-Squared: %g\r", reducedChiSqrd
	TextBox/C/N=text1/Z=1/X=85.00/Y=30.00 "\\F'Symbol'c\\F'Geneva'\\S2\\M/dof = "+num2str(reducedChiSqrd)	
	
	AppendResults(ParamsRawBeamGauss)
End



Function FitRawBeamGauss(pw, yw, xw) : FitFunc
	Wave pw, yw, xw
	
	//pw[0] = amplitude
	//pw[1] = 1s width
	//pw[2] = 2s width
	//pw[3] = offset
	//pw[4] = background
	//pw[5] = gaussian amplitude
	//pw[6] = gaussian width
	
	Variable I0 = pw[0] 
	Variable s1 = pw[1]
	Variable s2 = pw[2]
	Variable x0 = pw[3]
	Variable background = pw[4]
	Variable gaussAmpRatio = pw[5]
	Variable gaussWidth = pw[6]

	Variable detWidth = 60

	Variable L12 = .889E6   
	Variable L1gDet = 2.39736E6
	Variable OneGTwoSdist = 0.127E6
	Variable L2d = L1gDet+OneGTwoSdist
	Variable a = L2d/L12
	
	Variable pp = (1/2)*abs(s2+(s2-s1)*a)
	Variable dd = (1/2)*(s2+(s2+s1)*a)
	//printf "amp = %d, s1= %d, p = %d, d = %d\r" pw[0], s1, pp, dd
	
	Variable dX = 1
	
	Variable BeamPnts = 401// max(numpnts(xw)*2+1, 201)
	Make/D/O/N=(BeamPnts) RawBeam
	SetScale/P x  -dX*((BeamPnts-1)/2),dX,"um", RawBeam
	
	RawBeam = I0*( (dd+x)/(dd-pp)*(-dd<=x && x<-pp) + 1*(-pp<=x && x<pp) + (dd-x)/(dd-pp)*(pp<=x && x<=dd) )
	//RawBeam = I0*( (dd+x-x0)/(dd-pp)*(-dd<=x-x0 && x-x0<-pp) + 1*(-pp<=x-x0 && x-x0<pp) + (dd-x-x0)/(dd-pp)*(pp<=x-x0 && x-x0<=dd) )
	
	Variable detpnts = 101
	Make/D/O/N=(detpnts) Detector
	SetScale/P x -dX*((DetPnts-1)/2),dX, "um", Detector
	
	Detector = (x >= -detWidth/2 && x<= detWidth/2)
	//Detector = (x-x0 >= -detWidth/2 && x-x0<= detWidth/2)
	Variable sumdet = sum(Detector, -inf, inf)
	Detector /= sumdet
	
	Duplicate/O RawBeam RawBeamConvolved
	
	Convolve/A Detector, RawBeamConvolved
	
	yw = RawBeamConvolved(xw[p]-x0) + background + I0*gaussAmpRatio*exp(-(xw[p]-x0)^2/(2*gaussWidth^2)) ///(sqrt(2*pi)*GaussWidth)
End




///////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function TestRawBeam3sGauss()

	Variable halfScanWidth = 400		//microns
	Variable dX = 1						// microns/point (ie distance between points)
	
	Variable testpnts = halfScanWidth*2/dX+1

	Make/D/O/N=(testpnts) testY, testX
	SetScale/P x -halfScanWidth, dX, "um", testX, testY
	testX = x
	
	Make/D/O/N=1 testparams={40000, 40, 30, -5, 1, .0001, 75,50}
	
	FitRawBeam3sGauss(testparams, testY, testX)
	
	Duplicate/O testY testYErr
	testYErr = sqrt(testY)
	
	
End



Function GoFitRawBeam3sGauss(ywave, xwave, errwave)
	Wave ywave, xwave, errwave
	
	Make/O/T ParamsRawBeamGauss={"Amplitude", "1s width", "2s width", "x offset", "background", "gauss amp ratio", "gauss width", "s3 width"}
	Make/D/O/N=1 w_coef={32000, 10, 10, -5, 1.2, .01, 75, 50}
	Make/T/O/N=1 Constraints = {"K1>1", "K1<50", "K2>1", "K4>0.0001", "K5>0.000001", "K6>10", "K6<1000","K2<K1*3", "K2>K1*3"}		// "K2<K1*3", "K2>K1*3",
	Duplicate/O w_coef eps; eps = 1e-6; eps[5]=w_coef[5]/10; eps[7]=10; eps[1]=1; eps[2]=1
	
	FuncFit/M=2/L=600/H="00000001" FitRawBeam3sGauss W_coef  ywave /X=xwave /R /D /W=errwave /I=1 /E=eps /C=Constraints 
	
	Wave M_Covar
	Duplicate/O M_Covar, CorMat	 
	CorMat = M_Covar[p][q]/sqrt(M_Covar[p][p]*M_Covar[q][q])
	
	NVAR numHolds
	Variable reducedChiSqrd = V_chisq/(V_npnts-(numpnts(W_coef)-numHolds))
	printf "Reduced Chi-Squared: %g\r", reducedChiSqrd
	TextBox/C/N=text1/Z=1/X=85.00/Y=30.00 "\\F'Symbol'c\\F'Geneva'\\S2\\M/dof = "+num2str(reducedChiSqrd)	

	AppendResults(ParamsRawBeamGauss)
End



Function FitRawBeam3sGauss(pw, yw, xw) : FitFunc
	Wave pw, yw, xw
	
	//pw[0] = amplitude
	//pw[1] = 1s width
	//pw[2] = 2s width
	//pw[3] = offset
	//pw[4] = background
	//pw[5] = gaussian amplitude
	//pw[6] = gaussian width
	//pw[7] = 3s width
	
	Variable I0 = pw[0] 
	Variable s1 = pw[1]
	Variable s2 = pw[2]
	Variable x0 = pw[3]
	Variable background = pw[4]
	Variable gaussAmp = pw[5]
	Variable gaussWidth = pw[6]
	Variable s3 = pw[7]


	Variable detWidth = 70

	Variable L12 = 780e3
	Variable L1d = 2400e3
	Variable L2d = L1d - L12
	Variable a = L2d/L1d
	
	Variable pp = (1/2)*abs(s2+(s2-s1)*a)
	Variable dd = (1/2)*(s2+(s2+s1)*a)
	//printf "amp = %d, s1= %d, p = %d, d = %d\r" pw[0], s1, pp, dd
	
	Variable dX = 1	// microns/point 
	
	Variable BeamHalfwidth = 1000
	Variable BeamPnts = beamHalfwidth*2/dX+1
	Make/D/O/N=(BeamPnts) RawBeam
	SetScale/P x  -BeamHalfwidth, dX, "um", RawBeam
	
	RawBeam = ( (dd+x)/(dd-pp)*(-dd<=x && x<-pp) + 1*(-pp<=x && x<pp) + (dd-x)/(dd-pp)*(pp<=x && x<=dd) )
	Variable rawBeamNorm = area(RawBeam)
	RawBeam /= rawBeamNorm
	
	//RawBeam = I0*( (dd+x-x0)/(dd-pp)*(-dd<=x-x0 && x-x0<-pp) + 1*(-pp<=x-x0 && x-x0<pp) + (dd-x-x0)/(dd-pp)*(pp<=x-x0 && x-x0<=dd) )
	
	Variable DetHalfwidth = 50
	Variable detpnts = DetHalfwidth*2/dX+1
	Make/D/O/N=(detpnts) Detector
	SetScale/P x -DetHalfwidth, dX, "um", Detector
	
	Detector = (x >= -detWidth/2 && x<= detWidth/2)
	Variable sumdet = area(Detector)
	Detector /= sumdet
	
	//3s
	Variable s3Halfwidth = 500
	Variable s3pnts = s3Halfwidth*2/dX+1
	Make/D/O/N=(s3pnts) s3wave
	SetScale/P x -s3Halfwidth, dX, "um", s3wave
	
	s3wave = (x >= -s3/2 && x<= s3/2)
	Variable sums3 = area(s3wave)
	s3wave /= sums3
	
	Duplicate/O RawBeam RawBeamConvolved
	
	Convolve/A s3wave, RawBeamConvolved
	Convolve/A Detector, RawBeamConvolved
	
	yw = I0 * ( RawBeamConvolved(xw[p]-x0) + gaussAmp*exp(-(xw[p]-x0)^2/(2*gaussWidth^2)) )           + background
End




/////////////////////////////////////////////////////////////////////////////////

Function TwoSlits()
	Variable Slitpnts = 101
	Variable dX=1
	Make/D/O/N=(Slitpnts) OneS, TwoS
	SetScale/P x -dX*((SlitPnts-1)/2),dX, "um", OneS, TwoS
	
	Variable OneSWidth = 50
	OneS = (x >= -OneSWidth/2 && x<= OneSWidth/2)

	Variable TwoSWidth = 50
	TwoS = (x >= -TwoSWidth/2 && x<= TwoSWidth/2)
	
	Convolve/A OneS, TwoS
End




////////////////////////////////////////////////////////////////////////////////////



Function GoFitDoubleGauss(positionWave, countsWave, errorWave)
	Wave positionWave, countsWave, errorWave
	
	Make/D/O W_coef = {.4,366,0,.008,.04,.8}	//initial guesses
	Make/T/O dgconstraints = {"K1 > 0","K3 > 0","K4 > 0","K5 > 0"}
	
	FuncFit/TBOX=256/H="100000" DoubleGauss W_coef countsWave /X=positionWave /I=1 /D /W=errorWave /C=dgconstraints 
	
	Wavestats/Q countsWave; print "Chi-Squared per DOF: ", V_chisq/(V_npnts-13-1)
	
	//PrintProfileResults()
End



Function DoubleGauss(params, output, pos) : FitFunc
	Wave params, output, pos
	
	Variable bkg = params[0] 		// background counts
	Variable scale = params[1]		// scaling factor
	Variable x0 = params[2]			// mean
	Variable peakRatio = params[3]	// ratio of pedestal to peak
	Variable peakSig = params[4]	// std dev of peak gaussian
	Variable pedSig = params[5]		// std dev of pedestal
	
	output = 0																							// set output = 0
	output += scale * Exp (-(pos - x0)^2/(2 * peakSig^2))			// add central peak
	output += scale * peakRatio * Exp (-(pos - x0)^2/(2 * pedSig^2))	// add pedestal
	output += bkg																						// add background
End





Function DoubleGaussTest()
	Make/O/N=600 testwave = x/100 - 3
	Duplicate/O testwave testwaveout
	Make/O/N=5 testParams = {1,366,0,.008,.025,.8}
	
	//Display/K=1 testwaveout vs testwave
	//ModifyGraph mode=4,marker=8,msize=2
	DoubleGauss(testParams,testwaveout,testwave)
End

