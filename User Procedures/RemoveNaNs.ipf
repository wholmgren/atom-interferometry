#pragma rtGlobals=1		// Use modern global access method.


Function RemoveNaNsInfsXY(theXWave, theYWave)
	Wave theXWave
	Wave theYWave

	Variable p, numPoints, numNaNs
	Variable xval, yval
	
	numNaNs = 0
	p = 0											// the loop index
	numPoints = numpnts(theXWave)			// number of times to loop

	do
		xval = theXWave[p]
		yval = theYWave[p]
		if ((numtype(xval)!=0) %| (numtype(yval)!=0))		// either is NaN?
			numNaNs += 1
		else										// if not an outlier
			theXWave[p - numNaNs] = xval		// copy to input wave
			theYWave[p - numNaNs] = yval		// copy to input wave
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numNaNs, numNaNs, theXWave, theYWave
	
	return(numNaNs)
End


Function RemoveNaNsXYZ(theXWave, theYWave, theZWave) //adapted from Remove Points.ipf
	Wave theXWave
	Wave theYWave
	Wave theZWave

	Variable p, numPoints, numNaNs
	Variable xval, yval, zval
	
	numNaNs = 0
	p = 0											// the loop index
	numPoints = numpnts(theXWave)			// number of times to loop

	do
		xval = theXWave[p]
		yval = theYWave[p]
		zval = theZWave[p]
		if ((numtype(xval)==2) %| (numtype(yval)==2) %| (numtype(zval)==2))		// either is NaN?
			numNaNs += 1
		else										// if not an outlier
			theXWave[p - numNaNs] = xval		// copy to input wave
			theYWave[p - numNaNs] = yval		// copy to input wave
			theZWave[p - numNaNs] = zval
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numNaNs, numNaNs, theXWave, theYWave, theZWave
End




Function RemoveNaNsInfsXYZ(theXWave, theYWave, theZWave)
	Wave theXWave
	Wave theYWave
	Wave theZWave

	Variable p, numPoints, numNaNs
	Variable xval, yval, zval
	
	numNaNs = 0
	p = 0											// the loop index
	numPoints = numpnts(theXWave)			// number of times to loop

	do
		xval = theXWave[p]
		yval = theYWave[p]
		zval = theZWave[p]
		if ((numtype(xval)!=0) %| (numtype(yval)!=0) %| (numtype(zval)!=0))		// either is NaN?
			numNaNs += 1
		else										// if not an outlier
			theXWave[p - numNaNs] = xval		// copy to input wave
			theYWave[p - numNaNs] = yval		// copy to input wave
			theZWave[p - numNaNs] = zval
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numNaNs, numNaNs, theXWave, theYWave, theZWave
	
	return(numNaNs)
End




Function RemoveNaNsManyWaves(waveOfWaves)
	Wave/WAVE waveOfWaves
	
	Variable numWaves = numpnts(waveOfWaves)
	
	Print numWaves
End





Function RemoveNaNsXYZAB(theXWave, theYWave, theZWave, theAWave, theBWave) //adapted from Remove Points.ipf
	Wave theXWave
	Wave theYWave
	Wave theZWave
	Wave theAWave
	Wave theBWave

	Variable p, numPoints, numNaNs
	Variable xval, yval, zval, aval, bval
	
	numNaNs = 0
	p = 0											// the loop index
	numPoints = numpnts(theXWave)			// number of times to loop

	do
		xval = theXWave[p]
		yval = theYWave[p]
		zval = theZWave[p]
		aval = theAWave[p]
		bval = theBWave[p]
		if ((numtype(xval)==2) %| (numtype(yval)==2) %| (numtype(zval)==2) %| (numtype(aval)==2) %| (numtype(bval)==2) )		// either is NaN?
			numNaNs += 1
		else										// if not an outlier
			theXWave[p - numNaNs] = xval		// copy to input wave
			theYWave[p - numNaNs] = yval		// copy to input wave
			theZWave[p - numNaNs] = zval
			theAWave[p - numNaNs] = aval
			theBWave[p - numNaNs] = bval
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numNaNs, numNaNs, theXWave, theYWave, theZWave, theAWave, theBWave
End

Function RemoveNaNsXYZABC(theXWave, theYWave, theZWave, theAWave, theBWave, theCWave) //adapted from Remove Points.ipf
	Wave theXWave
	Wave theYWave
	Wave theZWave
	Wave theAWave
	Wave theBWave
	Wave theCWave

	Variable p, numPoints, numNaNs
	Variable xval, yval, zval, aval, bval, cval
	
	numNaNs = 0
	p = 0											// the loop index
	numPoints = numpnts(theXWave)			// number of times to loop

	do
		xval = theXWave[p]
		yval = theYWave[p]
		zval = theZWave[p]
		aval = theAWave[p]
		bval = theBWave[p]
		cval = theCWave[p]

		if ((numtype(xval)==2) %| (numtype(yval)==2) %| (numtype(zval)==2) %| (numtype(aval)==2) %| (numtype(bval)==2) %| (numtype(cval)==2))		// either is NaN?
			numNaNs += 1
		else										// if not an outlier
			theXWave[p - numNaNs] = xval		// copy to input wave
			theYWave[p - numNaNs] = yval		// copy to input wave
			theZWave[p - numNaNs] = zval
			theAWave[p - numNaNs] = aval
			theBWave[p - numNaNs] = bval
			theCWave[p - numNaNs] = cval
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numNaNs, numNaNs, theXWave, theYWave, theZWave, theAWave, theBWave, theCWave
End




Function RemoveNaNsXYZABCD(theXWave, theYWave, theZWave, theAWave, theBWave, theCWave, theDWave) //adapted from Remove Points.ipf
	Wave theXWave
	Wave theYWave
	Wave theZWave
	Wave theAWave
	Wave theBWave
	Wave theCWave
	Wave theDWave

	Variable p, numPoints, numNaNs
	Variable xval, yval, zval, aval, bval, cval, dval
	
	numNaNs = 0
	p = 0											// the loop index
	numPoints = numpnts(theXWave)			// number of times to loop

	do
		xval = theXWave[p]
		yval = theYWave[p]
		zval = theZWave[p]
		aval = theAWave[p]
		bval = theBWave[p]
		cval = theCWave[p]
		dval = theDWave[p]
		if ((numtype(xval)==2) %| (numtype(yval)==2) %| (numtype(zval)==2) %| (numtype(aval)==2) %| (numtype(bval)==2) %| (numtype(cval)==2) %| (numtype(dval)==2))		// either is NaN?
			numNaNs += 1
		else										// if not an outlier
			theXWave[p - numNaNs] = xval		// copy to input wave
			theYWave[p - numNaNs] = yval		// copy to input wave
			theZWave[p - numNaNs] = zval
			theAWave[p - numNaNs] = aval
			theBWave[p - numNaNs] = bval
			theCWave[p - numNaNs] = cval
			theDWave[p - numNaNs] = dval
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numNaNs, numNaNs, theXWave, theYWave, theZWave, theAWave, theBWave, theCWave, theDWave
End



Function RemoveNaNsXYZABCDE(theXWave, theYWave, theZWave, theAWave, theBWave, theCWave, theDWave, theEWave) //adapted from Remove Points.ipf
	Wave theXWave
	Wave theYWave
	Wave theZWave
	Wave theAWave
	Wave theBWave
	Wave theCWave
	Wave theDWave
	Wave theEWave
	
	Variable p, numPoints, numNaNs
	Variable xval, yval, zval, aval, bval, cval, dval, eval
	
	numNaNs = 0
	p = 0											// the loop index
	numPoints = numpnts(theXWave)			// number of times to loop

	do
		xval = theXWave[p]
		yval = theYWave[p]
		zval = theZWave[p]
		aval = theAWave[p]
		bval = theBWave[p]
		cval = theCWave[p]
		dval = theDWave[p]
		eval = theEWave[p]
		if ((numtype(xval)==2) %| (numtype(yval)==2) %| (numtype(zval)==2) %| (numtype(aval)==2) %| (numtype(bval)==2) %| (numtype(cval)==2) %| (numtype(dval)==2) %| (numtype(eval)==2))		// either is NaN?
			numNaNs += 1
		else										// if not an outlier
			theXWave[p - numNaNs] = xval		// copy to input wave
			theYWave[p - numNaNs] = yval		// copy to input wave
			theZWave[p - numNaNs] = zval
			theAWave[p - numNaNs] = aval
			theBWave[p - numNaNs] = bval
			theCWave[p - numNaNs] = cval
			theDWave[p - numNaNs] = dval
			theEWave[p - numNaNs] = eval
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numNaNs, numNaNs, theXWave, theYWave, theZWave, theAWave, theBWave, theCWave, theDWave
End
