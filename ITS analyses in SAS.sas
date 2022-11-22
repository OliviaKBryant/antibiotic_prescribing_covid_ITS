/******************************************************************************
   Project: ITS analysis for impact of COVID-19 on antibiotic prescribing
   Author: Wallis Lau
   Date created: 31/10/2022
   Date updated: 22/11/2022
   Notes: ITS analysis of COVID-19 on antibiotic prescriptions
******************************************************************************/


/** Import CSV dataset **/
PROC IMPORT OUT= ITS 
            DATAFILE= "H:\ITS_processed_dataset.csv" 
            DBMS=CSV REPLACE;
            GETNAMES=YES;
            DATAROW=2; 
RUN;

/** 1. Analyses for Prescribing rate **/
data ITS2;
set ITS;
ln_appointments=log(Appointments); 
run;

proc genmod data=ITS2;
	class studymonth; /** set month as a categorical variable using the "class" statement **/
    model Item_all_antibiotics = time /** baseline trend **/
                                 intervention1 intervention2 /** step change **/
                                 time_after_inter1 time_after_inter2 /** slope change **/
                                 studymonth /** seasonal effect **/
                                / dist=negbin /** negative binomial **/
                                  link=log /** loglink function **/ 
                                  offset=ln_appointments   
  ;
  ods output ParameterEstimates=Result_All_AB_per_appt; 
  output out=residual_All_AB_per_appt resraw = Resraw  RESDEV= RESDEV ;  
  run;

  /** Check residual **/
  proc autoreg data=residual_All_AB_per_appt all plots;
  model RESDEV=time / method=ml;
  title "Residual for [All Antibiotics]"; 
  run;

  /** Adjust for 1st order autocorrelation **/
  proc sort data=residual_All_AB_per_appt;
    by time;
  run;
  data residual_All_AB_per_appt;
  set residual_All_AB_per_appt;
    lag1_RESDEV=lag(RESDEV); 
  run;

  proc genmod data=residual_All_AB_per_appt;
	class studymonth;
    model Item_all_antibiotics = time /** baseline trend **/
                                 intervention1 intervention2 /** step change **/
                                 time_after_inter1 time_after_inter2 /** slope change **/
                                 studymonth /** seasonal effect **/
					             lag1_RESDEV /** lag_residual to adjust for 1st order autocorrelation **/
                                 / dist=negbin link=log offset=ln_appointments   
  ;
  ods output ParameterEstimates=Result_All_AB_per_appt2;
  output out=residual_All_AB_per_appt2 resraw = Resraw2  RESDEV= RESDEV2;
  run;

  /** Check residual **/
  proc autoreg data=residual_All_AB_per_appt2 all plots;
    model  RESDEV2=time / method=ml ;
    title "Residual for [All Antibiotics] lag1 ";
  run; 

  /** Output results **/
  data ITS_results;
  set Result_All_AB_per_appt2; 
    IRR=round(exp(estimate),.0001);
    IRR_lci=round(exp(LowerWaldCL),.0001);
    IRR_uci=round(exp(UpperWaldCL),.0001);
    keep parameter ProbChiSq IRR IRR_lci IRR_uci;
  run;

/** 2. Absolute number of prescriptions **/
proc genmod data=ITS;
	class studymonth; /** set month as a categorical variable using the "class" statement **/
    model Item_all_antibiotics = time /** baseline trend **/
                                 intervention1 intervention2 /** step change **/
                                 time_after_inter1 time_after_inter2 /** slope change **/
                                 studymonth /** seasonal effect **/
                                / dist=negbin /** negative binomial **/
                                  link=log; /** loglink function **/ 
                                  
  ods output ParameterEstimates=Result_All_AB_abs; 
  output out=residual_All_AB_abs resraw = Resraw  RESDEV= RESDEV ;  
  run;

  /** Check residual **/
  proc autoreg data=residual_All_AB_abs all plots;
  model RESDEV=time / method=ml;
  title "Residual for [All Antibiotics] prescriptions"; 
  run;

  /** Adjust for 1st order autocorrelation **/
  proc sort data=residual_All_AB_abs;
    by time;
  run;
  data residual_All_AB_abs;
  set residual_All_AB_abs;
    lag1_RESDEV=lag(RESDEV); 
  run;

  proc genmod data=residual_All_AB_abs;
	class studymonth;
    model Item_all_antibiotics = time /** baseline trend **/
                                 intervention1 intervention2 /** step change **/
                                 time_after_inter1 time_after_inter2 /** slope change **/
                                 studymonth /** seasonal effect **/
				       	         lag1_RESDEV /** lag_residual to adjust for 1st order autocorrelation **/
                                / dist=negbin link=log   
  ;
    ods output ParameterEstimates=Result_All_AB_abs2;
    output out=residual_All_AB_abs2 resraw = Resraw2  RESDEV= RESDEV2;
  run;

  /** Check residual again **/
  proc autoreg data=residual_All_AB_abs2 all plots;
    model RESDEV2=time / method=ml ;
    title "Residual for [All Antibiotics] prescriptions lag1 ";
  run;

  /** Output results **/
  data ITS_results_prescription;
  set Result_All_AB_abs2; 
    IRR=round(exp(estimate),.0001);
    IRR_lci=round(exp(LowerWaldCL),.0001);
    IRR_uci=round(exp(UpperWaldCL),.0001);
    keep parameter ProbChiSq IRR IRR_lci IRR_uci;
  run;

