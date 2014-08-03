#pragma rtGlobals=1		// Use modern global access method.



Function SetScaleTest()
	Variable xPnts = 10
	Variable yPnts = 10
	
	Variable xDimMin = -0.5
	Variable xDimMax = 0.5
	
	Variable yDimMin = 100
	Variable yDimMax = 200
	
	Make/O/N=(xPnts) oneDwave
	Make/O/N=(xPnts,yPnts) twoDwave1, twoDwave2
	
	SetScale/I x xDimMin, xDimMax, oneDwave, twoDwave1, twoDwave2		// works as you'd expect
	SetScale/I y yDimMin, yDimMax, oneDwave, twoDwave1					// twoDwave1 still gets its proper scaling despite the odd presence of oneDwave in the command
	SetScale/I y yDimMin, yDimMax, twoDwave2, oneDwave					// twoDwave2 does not get scaled due to the presence of oneDwave in the command
	
	oneDwave = x
	twoDwave1 = y
	twoDwave2 = y
End