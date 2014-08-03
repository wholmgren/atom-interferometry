#pragma rtGlobals=3		// Use modern global access method and strict wave access.


function doForecast(timewave,csreference,windx,windy,dp)
	wave timewave,windx,windy,csreference
	variable dp//number of points ahead
	variable n = numpnts(timewave)
	wave clouds,latitude,longitude
	make/o /n=(n) w_measurement=0,w_forecast=0,w_forecast_temp=0,w_clouds1=0,w_clouds2=0,w_persistence=0
	duplicate/o timewave w_forecast_time
	make/o/n=2 location,windvector
	//windvector={-.06,.04}

	location={0,5}
	variable i
	variable t0,t1,dt
	t0=timewave[0]
	t1=timewave[1]
	dt=(t1-t0)*dp //t1-t0 is in units of seconds
	for(i=0;i<n;i+=1)
		windvector={windx[i],windy[i]}
		if(mod(i,1000)==0 && i!=0)
			print windvector
		endif
		windvector/=1000 //m/s -> km per second
		getAllSystemsByTime(timewave[i])


		w_forecast_temp[i]=inerpLatLong(location[0]-windvector[0]*dt,location[1]-windvector[1]*dt,longitude,latitude,clouds)//forecast
		w_forecast_time[i]=timewave[i]+dt

		w_clouds1[i-dp]=inerpLatLong(location[0],location[1],longitude,latitude,clouds)//measurement

		w_clouds2[i]=inerpLatLong(location[0],location[1],longitude,latitude,clouds)//persistance
	endfor
	for(i=1;i<=dp;i+=1)//make sure measurement is available for all times for which forecast is evaluated
		getAllSystemsByTime(timewave[n-1]+i*dt)
	
		w_clouds1[n-1+i-dp]=inerpLatLong(location[0],location[1],longitude,latitude,clouds)
	endfor


	w_measurement=interp(w_forecast_time[p],timewave,csreference)*(1-w_clouds1)//*(w_measurement<daycsreference)

	w_forecast=interp(w_forecast_time[p],timewave,csreference)*(1-w_forecast_temp)//*(w_forecast<daycsreference)

	w_persistence=interp(w_forecast_time[p],timewave,csreference)*(1-w_clouds2)


end



function doForecastByDate(startdate,enddate,dt,velmethod)
	variable startdate,enddate,dt,velmethod
	getReferenceByDate(startdate,enddate,cloudyonly=1)
	wave dayreferencetime,daycsreference,bestWindx,bestWindy,bestWinddate
	make/o/n=(numpnts(dayreferencetime)) windx,windy

	if(velmethod==0)
		windx=getWindNCDC(dayreferencetime,0)
		windy=getWindNCDC(dayreferencetime,1)
	elseif(velmethod==1)
		windx=getWindLeuth(dayreferencetime+dt*60/2,0)
		windy=getWindLeuth(dayreferencetime+dt*60/2,1)
		windx*=-1.5
		windy*=-1.5
	elseif(velmethod==2)
		windx=getWindMeas(dayreferencetime,0)
		windy=getWindMeas(dayreferencetime,1)
	elseif(velmethod==3)

		windx=interp(floor(dayreferencetime/3600)*3600,bestWinddate,bestwindx)

		windy=interp(floor(dayreferencetime/3600)*3600,bestWinddate,bestwindy)
	elseif(velmethod==4)

		windx=interp(floor(dayreferencetime/3600-2)*3600,bestWinddate,bestwindx)

		windy=interp(floor(dayreferencetime/3600-2)*3600,bestWinddate,bestwindy)
	elseif(velmethod==5)

		windx=interp(floor(dayreferencetime/3600)*3600,bestWinddate,kalman_wind_u)

		windy=interp(floor(dayreferencetime/3600)*3600,bestWinddate,kalman_wind_v)
	elseif(velmethod==6)
		windx=interp(dayreferencetime,bestWinddate,kalman_wind_u)
		windy=interp(dayreferencetime,bestWinddate,kalman_wind_v)
	elseif(velmethod==7)
		windx=gnoise(3)
		windy=gnoise(3)

	endif

	doForecast(dayreferencetime,daycsreference,windx,windy,dt/15)


	wave w_measurement,w_forecast,w_clouds,w_persistence
	print "clearness Index:"
	print sum(w_measurement)/sum(daycsreference)
	print "forecast:"
	ErrorStats(w_measurement,w_forecast)
	print "persistence"
	ErrorStats(w_measurement,w_persistence)
end



function getAllSystemsByTime(seconds)
	variable seconds
	String connectionStr = "DSN=tfsdata;UID=root;PWD="
	string latlonselector=" && latitude>'32.0' && latitude<'32.7' && longitude<'-110.7' && longitude > '-111.2'"
	string idexclude=" && pm.id<>31 && pm.id<>7 && pm.id<>97 && pm.id<>21"
	string statement = "select pm.id as ids,random,finalyield,if(csreference>0.01,cloudderating,0) as clouds,latitude,longitude Êfrom processedmeasurements as pm join (locations as l, systemparameters as s) on (l.id=pm.id && s.id=pm.id) where longitude<>0 && s.errorcode=0 Ê&& pm.outage=0 && pm.shadederating=0 && time='"+Secs2Date(seconds,-2)+" "+Secs2Time(seconds,3)+"'"+latlonselector
	//print statement
	SQLHighLevelOp /CSTR={connectionStr,SQL_DRIVER_COMPLETE} /O statement
	wave latitude,longitude //{-110.984,32.3384}
	latitude=(latitude-32.2352)*2*pi*6353/360
	longitude=(longitude+110.944)*2*pi*6353*cos(32*pi/180)/360
	// ÊÊÊÊÊTextBox/C/N=text0 Secs2Date(seconds,-2)+" "+Secs2Time(seconds,2)
end



function getReferenceByDate(startDate,endDate,[cloudyonly])
	variable startDate,endDate,cloudyonly
	if( ParamIsDefault(cloudyonly))
		cloudyonly=0
	endif
	string startdatestring =Secs2Date(startDate,-2)+" "+Secs2Time(startDate,3)
	string enddatestring = Secs2Date(endDate,-2)+" "+Secs2Time(endDate,3)
	String connectionStr = "DSN=tfsdata;UID=root;PWD="
	string timeselector=""//" && hour(time)>=5 && hour(time)<=19"
	string cloudyselector
	if(cloudyonly==0)
		cloudyselector=""
	else
		cloudyselector=" Ê&& iscloudy=1 "
	endif

	string statement = "select time as dayreferencetime,csreference as daycsreference from reference join iscloudy on (iscloudy.date=date(reference.time)) where Êtime>='"+startdatestring+"' && time<'"+enddatestring+"'"+cloudyselector+timeselector
	SQLHighLevelOp /CSTR={connectionStr,SQL_DRIVER_COMPLETE} /O statement
end




function Kalman(measurement,mX,mH,mA,mQ,mP,mR)
	// this function contains the inner workings of the Kalman algorithm
	// for more details see "An Introduction to the Kalman Filter" by ÊGreg Welch and Gary Bishop
	// the inputs are:
	// measurement: the most recent measurement
	// mX the most recent state estimate
	// mH the operator that converts the state into a measurement: meas = mH x mX
	// mA is the propagator: ÊmX1 = mA x mX0
	// mQ is the covariance matrix of the process noise
	// mR is the covariance matrix of the measurement noise
	// mP is the covariance matrix of the state estimate errors

	wave measurement,mH,mA,mQ,mP,mX,mR
	wave Z = measurement

	variable nParam=numpnts(mX)
	matrixop/o mPminus= mA x mP x mA^t +mQ // apriori estimate of state estimate errors

	//kalaman update
	// mK is the weighting factor that determines the relative importance of the measurement and the extrapolation from previous measurements
	// this matrix is determined based on the apriori estimate of state estimate errors and the measurement error
	matrixop/o mK = mPminus x mH^t x inv( mH x mPminus x mH^t +mR)

	// make estimate of present state based on previous state
	propagate(mX,mA,1)

	// present state estimate gets updated with information from measurement
	matrixop/o mX=mX+ mK x (- mH x mX+Z)
	//estimate of state estimate errors gets updated
	matrixOp/o mP=(identity(nParam)-mK x mH) x mPminus

end

function Propagate(mX,mA,n)
	// propagate state vector
	// note that the wave gets altered!!!

	wave mX,mA
	variable n

	variable i=0
	for (i=0;i<n;i+=1)
		matrixop/o mX= mA x mX
	endfor
end

function inerpLatLong(long,lat,longwave,latwave,fywave)
	variable lat,long
	wave latwave,longwave,fywave

	duplicate/o latwave interpdistances_temp,out,ones
	ones=1
	interpdistances_temp =sqrt( (latwave-lat)^2 + (longwave-long)^2)
	//variable n =numpnts(fywave)
	// ÊÊÊÊÊmatrixop/o out = sum(fywave*exp(-distances/.005))/sum(exp(-distances/.005))
	//matrixop/o out = Êsum( Êfywave Ê/ Ê(distances+.005) Ê) //1/ Êsum(ones Ê/ Ê(distances+.005) Ê)
	//
	//return out
	variable n=4,i
	make/o/n=4 values
	for(i=0;i<n;i+=1)
		wavestats/q interpdistances_temp
		values[i] = fywave[v_minloc]//*(distances[v_minloc]<.1)
		interpdistances_temp[v_minloc]=100
	endfor
	sort/q values values
	return values(2.5)

end