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
	print avg, " ± ", stddev
End

