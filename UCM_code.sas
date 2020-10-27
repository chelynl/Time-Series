/* UCM Modeling in SAS  */

*Create library;
libname time "/home/u37604878/Time Series";

*Append training/validation;
proc append base=time.training data=time.valid force;
run;

/* VISUALIZE OZONE CONCENTRATION */
%macro ts(s=maxOzoneC);
	proc sgplot data=time.training;
		series x=date y=&s;
	run;
%mend;
quit;
%ts;


/* FITTING FOR OZONE CONCENTRATION */
*Start with level model;
proc ucm data=time.training plots=all;
	level;
	irregular;
	model maxOzonec;
	forecast back=28 lead=28 plot=forecasts;
run;

*Check potential trig modeling for seasonality;
proc ucm data=time.training plots=all;
	level;
	irregular;
	model maxOzonec;
	season length=365 type=trig keeph=1 to 7;
	estimate plot=wn;
	forecast back=28 lead=28 plot=forecasts;
run;

proc ucm data=time.training plots=all;
	level;
	irregular;
	model maxOzonec;
	estimate plot=wn;
	splineseason length=365 knots=2, 3;
	forecast back=28 lead=28 plot=forecasts;
run;

*Evaluate trend component;
proc ucm data=time.training plots=residuals(loess);
	level;
	irregular;
	model maxOzonec;
	estimate plot=wn;
	slope;
	forecast back=28 lead=28 plot=forecasts;
run;

*See plots of potential predictors;
%ts;
%ts(s=maxNOconc);
%ts(s=maxSO2conc);
%ts(s=maxCOconc);
%ts(s=TMAX);
%ts(s=PRCP);
%ts(s=WSF5);
%ts(s=AWND);

*Regression model for predictors of ozone. Use regressors from ARIMA section;
%let vars = maxNOconc maxSO2conc maxCOconc TMAX TAVG PRCP WSF2 WSF5 AWND;
proc glmselect data=time.training plots=all;
	Backward: model maxOzonec=&vars/selection=backward select=SL slstay=0.01;
run;

*Model Xs from model;
proc ucm data=time.training plots=residuals(loess);
	level;
	irregular;
	model maxOzonec = maxNOconc maxCOconc TMAX PRCP WSF5;
	estimate plot=wn;
	forecast back=28 lead=28 plot=forecasts;
run;

*Test various AR, MA, seasonal terms based on autocorrelation plots for white noise;
proc ucm data=time.training plots=residuals(acf pacf loess);
	level;
	irregular p=12; * Tested multiples of AR, MA, and seasonal terms: p = x q = x sp = x sq = x;
	model maxOzonec=maxNOconc maxCOconc TMAX PRCP WSF5;
	estimate plot=wn;
	forecast back=28 lead=28 plot=forecasts;
run;


*Append training/test;
proc append base=time.training data=time.test force;
run;

*Forecast and export validation/test;
ods excel file="/home/u37604878/Time Series/UCM_VT.xlsx";
proc ucm data=time.training plots=none;
	level;
	irregular p=12;
	model maxOzonec=maxNOconc maxCOconc TMAX PRCP WSF5;
	estimate plot=wn;
	forecast back=42 lead=42 plot=forecasts;
run;
ods excel close;

ods excel file="/home/u37604878/Time Series/UCM_T.xlsx";
proc ucm data=time.training plots=none;
	level;
	irregular p=12;
	model maxOzonec=maxNOconc maxCOconc TMAX PRCP WSF5;
	estimate plot=wn;
	forecast back=14 lead=14 plot=forecasts;
run;
ods excel close;
