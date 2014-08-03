#pragma rtGlobals=1		// Use modern global access method.




Function WeightedAverage()
	Variable num1, num2, sig1, sig2, avg, stddev
	
	num1= 3037.9
	sig1 = 3.9

	num2 = 3044.8
	sig2 = 2.2
	
	avg = (num1/sig1^2+num2/sig2^2)/(1/sig1^2+1/sig2^2)
	stddev = sqrt(1/(1/sig1^2+1/sig2^2))
	print avg, " ± ", stddev
End


Function WeightedAverage2(num1, sig1, num2, sig2)
	Variable num1, sig1, num2, sig2
	
	Variable avg, stddev
	
	avg = (num1/sig1^2+num2/sig2^2)/(1/sig1^2+1/sig2^2)
	stddev = sqrt(1/(1/sig1^2+1/sig2^2))
	//print/d avg, " ± ", stddev
	printf "%3.4f ± %.4f\r" avg, stddev
End



Function WeightedAverageWaves(numbers, sigmas)
	Wave numbers, sigmas
	
	Variable avg, stddev, sigsum
	
	Variable i
	for(i=0; i<numpnts(numbers); i+=1)
		avg += numbers[i]/sigmas[i]^2
		sigsum += 1/sigmas[i]^2
	endfor
	avg/=sigsum
	stddev = sqrt(1/sigsum)
	
	print avg, " ± ", stddev
End
