/***************************************************************************************************
****************************************************************************************************;
Macro Name: Plots_CIF
Created Date/Author/Contact: Jan 2019/  Chao Zhang, Yaqi Jia, and Yuan Liu
Current Version: V2
Working Environment: SAS 9.4 English version TS1M3   

Contact Information: chao.zhang2@emory.edu

Purpose:  :  To create a cumulative incidence plot. The SAS macro %Plots_CIF implements appropriate
Nonparametric methods for estimating cumulative incidence functions. The macro also implements     
Gray’s method (Gray 1988) for testing differences between these functions in multiple groups. And 
also, the number at risk of specified time-point with accumulated number of event of interest, 
accumulated number of competing by the specified time-point, and censored information was shown below
the CIF plot.    

Note: All calculation is based on PROC Lifetest in SAS 9.4 with EVENTCODE option. See SAS help manual 
for more technique details. 

References: 1 SAS Institute Inc, 2017 SAS/STAT® 14.3 User’s Guide. Cary, NC: SAS Institute Inc. User’s Guide the Lifetest Procedure.  
            2 Sanjay Matange, Annotate your SGPLOT Graphs. Paper CC01-2014

Parameters: 

DSN            The name of the data set to be analyzed(required).   
 
GRPLIST        The variable list that defines the groups for comparison (optional) separated by space.

TIME_EVENT     Name of time to event outcome variable (required).                                                             

CENSOR         Name of censoring indicator variable.  Values of 0 indicate censored observations(required).

EVENTCODE      The value in CENSOR that indicate event of interest, and this value will appear 
               EVENTCODE= option. The default value is 1.

YAXISVALUE     Specify the ticket value of y axis. 

XAXISVALUE	   Usage xaxisvalue = 0 12 24 36 ...; specify the ticket value of X axis. Number atRisk, 
               accumulated #events of interest, accumulated #competing shown correspond to the values of X axis.

NatRisk        List the number at risk of specified time-point, accumulated number of events of interest, and accumulated number of competing.
               The default value is F.

TIMELIST       List of time points separated by spaces to report survival estimates and 95% CI 
               (optional).  

UNITS          Units of the time variable, i.e. days, months, etc. The default value is none.
OUTPATH        Path for output table to be stored.
FNAME          File name for output table.
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted in
               debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

******************************************************************************************************/

%macro Plots_CIF(dsn=, grplist=, time_event=, censor=, eventcode=1, xaxisvalue=, yaxisvalue=, NatRisk=F, timelist=,units=,outpath=,filename=,debug=F);

    /* Capitalize */
    
    %let NatRisk = %UPCASE(&NatRisk);
   /* Count number of time estimates requested */
   %let time_cnt = %sysfunc(countw(&timelist,' '));

   /*count number of &grplist*/
    %let var_cnt = %sysfunc(countw(&grplist,' '));

   ods rtf file="&OUTPATH.&FILENAME &SYSDATE..DOC.doc" startpage=NO image_dpi=800 ;
   
   %if &var_cnt = 0 %then %do;

		ODS SELECT NONE;
           proc lifetest data= &dsn plots=cif(test) outcif= _est ; 
              time &time_event*&censor(0)/eventcode= &eventcode;
        run;
   
		ODS SELECT ALL;

	data  _est;
	  set _est;
	  label CIF = "Cumulative Incidence";
	run;

    ods select all;
	  
	 title height=12pt font=Arial  "CIF Plot for &time_event" ;
      proc sgplot data= _est;
	      %if &yaxisvalue ~= %str() %then %do; yaxis values=(&yaxisvalue);  %end;
          %if &xaxisvalue ~= %str() %then %do; xaxis values= (&xaxisvalue); %end;
		step y=CIF x=&time_event ;
     run;
   title ;

     ODS SELECT none ;
       proc lifetest data= &dsn plots=cif(test) timelist=&timelist outcif= _est_list REDUCEOUT; 
              time &time_event*&censor(0)/eventcode= &eventcode;
       run;

      data _est_list2;
	    set _est_list;
        /* Convert to character */

		cif_c = PUT(cif,best5.3 );  lcl_c = PUT(CIF_LCL, best5.3); ucl_c = PUT(CIF_UCL, best5.3);

		if cif = . then cif_c= "NA"; if cif_LCL = . then lcl_c= "NA"; if cif_Ucl =. then ucl_c= "NA"; 
       
       rate = TRIM(cif_c) || " (" || TRIM(LEFT(LCL_c)) || ", " || TRIM(LEFT(UCL_c)) || ")";

     run;

    ods select all;
     proc report nowd data=_est_list2 style(report)={rules=GROUPS}  
            style(header)={BACKGROUNDCOLOR=none} LS=256 STYLE(COLUMN) = {JUST = C};
         col TIMELIST  rate;
        * define &grp/order style={just=l};
		 %if &units ~= %str() %then %do;
         LABEL TIMELIST = "TIME (&units)"
		 %end;
		 %else %do;
		 LABEL TIMELIST="Time"
		 %end;
            rate = "CIF Estimate (95% CI)"; 
      run;
    %end; /*End of %if &var_cnt = 0 %then %do;*/




  %else %if &var_cnt > 0 %then %do;
   %let n=1;   
       %DO %UNTIL (%scan(&grplist,&n)= );
          %let grp=%scan(&grplist,&n);
       
	  
     ODS SELECT none ;
		 ods output GrayTest =GrayTest; 
         proc lifetest data= &dsn plots=cif(test) outcif= _est ; 
              time &time_event*&censor(0)/eventcode= &eventcode;
              strata &grp / order=internal;
         run;

       ods output Survivalplot=Survivalplot;
          proc lifetest data=&dsn plots=survival( %if &xaxisvalue ~= %str() %then atrisk= &xaxisvalue ; %else atrisk;);
           time &time_event*&censor(0);
            strata &grp / order=internal  ;
         run;


          Proc sql noprint;
             select ProbChiSq format=pvalue6.4 into:pval1
             from GrayTest;
          quit;

		 data _est1;
		    set _est;
			by stratum;
			if first.stratum then do ; total_event=0; total_AllEventTypes =0; end;
			 total_event + event; total_AllEventTypes + AllEventTypes ;
			 total_competing = total_AllEventTypes - total_event;
             total_event_ = lag(total_event); total_competing_ =lag( total_competing);
            if last.stratum then do; total_event_ = total_event; total_competing_ = total_competing; end;
			if first.stratum then delete;
		 run;

          data _est1;
		    set _est1;
			if censored = 1 then CIF_Censored = CIF;
			rename stratum =stratumNum;
			label CIF_Censored ='Censored';
			risk = TRIM(AtRisk) || " (" || TRIM(LEFT(total_event_)) || ", " || TRIM(LEFT(total_competing_)) || ")";
            
		  run;

		  
 PROC SQL;
    CREATE TABLE _EST_plot AS
      SELECT A.*, B.tAtRisk,B.Stratum
      from _EST1 A LEFT JOIN  Survivalplot B
       ON _EST1.StratumNum = Survivalplot.StratumNum and _EST1.AtRisk = Survivalplot.AtRisk ;
QUIT;


  data _est_plot;
    set _est_plot;
    label CIF = "Cumulative Incidence";
  run; 
    

     data Anno_AtRisk;
          
        length function $15 x1space $12 label $30 ; 
        set _EST_plot ;
        by stratum;
        retain id 0 function x1space  y1space  textcolor textsize textweight anchor width widthunit y1;

        if first.stratum then do;
            id+1; textcolor="GraphData" || put(id, 1.0) || ":contrastcolor"; 
            function="text"; x1space="datavalue"; y1space="graphpercent"; widthunit="percent"; 
            textsize=9; textweight="normal"; anchor="center"; width=30; y1=(id-1)*4+3;  end;
        
        if tatrisk ne . then do; x1=tatrisk; label=risk; output; end; 
      run;

	 
     ods select all;
	
      /* title height=12pt font=Arial  "CIF Plot With Number At Risk" ;*/
      proc sgplot data= _est_plot  %if &NatRisk = T %then %do; pad=(bottom=15% right=8%)  sganno = Anno_AtRisk;%END;;
	      %if &yaxisvalue ~= %str() %then %do; yaxis values=(&yaxisvalue);%end;
		  %else %do;
          yaxis values=(0.0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0);
		  %end;
          %if &xaxisvalue ~= %str() %then %do; xaxis values= (&xaxisvalue);%end;
          step y=CIF x=&time_event /group=&grp name ="CIF";
		  scatter x=&time_event y=CIF_Censored / markerattrs=(symbol=plus) name='censored';
          scatter x=&time_event y=CIF_Censored / group=&grp markerattrs=(symbol=plus);
          inset ("Gray's Test p="="%cmpres(&PVAL1)")/position=topleft border LABELALIGN=LEFT;
		  keylegend "CIF";
		  keylegend "censored" / location=inside position=topright;
	  run;
	  
	

	 ODS SELECT none ;
       proc lifetest data= &dsn plots=cif(test) timelist=&timelist outcif= _est_list REDUCEOUT; 
              time &time_event*&censor(0)/eventcode= &eventcode;
              strata &grp / order=internal;
         run;

      data _est_list2;
	    set _est_list;
        /* Convert to character */

		/*cif_c = PUT(cif,percent8.2);  lcl_c = PUT(CIF_LCL,percent8.2); ucl_c = PUT(CIF_UCL,percent8.2);*/
		cif_c = PUT(cif,best5.3 );  lcl_c = PUT(CIF_LCL, best5.3); ucl_c = PUT(CIF_UCL, best5.3);

		if cif = . then cif_c= "NA"; if cif_LCL = . then lcl_c= "NA"; if cif_Ucl =. then ucl_c= "NA"; 
       
       rate = TRIM(cif_c) || " (" || TRIM(LEFT(LCL_c)) || ", " || TRIM(LEFT(UCL_c)) || ")";

     run;

    ods select all;
     proc report nowd data=_est_list2 style(report)={rules=GROUPS}  
            style(header)={BACKGROUNDCOLOR=none} LS=256 STYLE(COLUMN) = {JUST = C};
         col &grp TIMELIST  rate;
         define &grp/order style={just=l};
		 %if &units ~= %str() %then %do;
         LABEL TIMELIST = "TIME (&units)"
		 %end;
		 %else %do;
		 LABEL TIMELIST="Time"
		 %end;
            rate = "CIF Estimate (95% CI)"; 
      run;


%let n=%eval(&n+1);
%end;
%end;
ods rtf close;
/* If not in debug mode */
   %if &debug = F %then %do;
      /* Delete intermediate datasets */
      PROC DATASETS lib=work noprint;
        delete anno_atrisk Graytest  _est  _est1  _est_list  _est_list2 _est_plot;
	  QUIT;
   %end;
%mend Plots_CIF ;


***example ***;

proc format;
value DiseaseGroup 1='ALL' 2='AML-Low Risk' 3='AML-High Risk';
value sex 0= 'F' 1='M';
run;
data bmt;
  input Group T Status WaitTime @@;
  logWaittime=log(WaitTime);
  label T = ' disease-free survival time (days)' ;
datalines;
1 2081 0 98 1 1602 0 1720 1 1496 0 127  1 1462 0 168
1 1433 0 93 1 1377 0 2187 1 1330 0 1006 1 996  0 1319
1 226  0 208 1 1199 0 174  1 1111 0 236  1 530  0  151 
1 1182 0 203 1 1167 0 191 1 418 2 110 1 383 1 824
1 276  2 146 1 104 1 85   1 609 1 187  1 172 2 129
1 487 2 128  1 662 1 84   1 194 2 329 1 230 1 147 
1 526 2 943 1 122 2 2616  1 129 1 937 1 74 1 303
1 122 1 170 1 86 2 239    1 466 2 508 1 192 1 74
1 109 1 393 1 55 1 331    1   1 2 196 1 107 2 178
1 110 1 361 1 332 2 834   2 2569 0 270 2 2506 0 60
2 2409 0 120 2 2218 0 60 2 1857 0 90 2 1829 0 210
2 1562 0 90 2 1470 0 240 2 1363 0 90 2 1030 0 210
2 860 0 180 2 1258 0 180 2 2246 0 105 2 1870 0 225
2 1799 0 120 2 1709 0 90 2 1674 0 60 2 1568 0 90
2 1527 0 450 2 1324 0 75 2 957 0 90 2 932 0 60
2 847 0 75 2 848 0 180 2 1850 0 180 2 1843 0 270
2 1535 0 180 2 1447 0 150 2 1384 0 120 2 414 2 120
2 2204 2 60 2 1063 2 270 2 481 2 90 2 105 2 120
2 641 2 90 2 390 2 120 2 288 2 90 2 421 1 90
2 79 2 90 2 748 1 60 2 486 1 120 2 48 2 150
2 272 1 120 2 1074 2 150 2 381 1 120 2 10 2 240
2 53 2 180 2 80 2 150 2 35 2 150 2 248 1 30
2 704 2 105 2 211 1 90 2 219 1 120 2 606 1 210
3 2640 0 750 3 2430 0 24 3 2252 0 120 3 2140 0 210
3 2133 0 240 3 1238 0 240 3 1631 0 690 3 2024 0 105
3 1345 0 120 3 1136 0 900 3 845 0 210 3 422 1 210
3 162 2 300 3 84 1 105 3 100 1 210 3 2 2 75
3 47 1 90 3 242 1 180 3 456 1 630 3 268 1 180
3 318 2 300 3 32 1 90 3 467 1 120 3 47 1 135
3 390 1 210 3 183 2 120 3 105 2 150 3 115 1 270
3 164 2 285 3 93 1 240 3 120 1 510 3 80 2 780
3 677 2 150 3 64 1 180 3 168 2 150 3 74 2 750
3 16 2 180 3 157 1 180 3 625 1 150 3 48 1 210
3 273 1 240 3 63 2 360 3 76 1 330 3 113 1 240
3 363 2 180
;
run;

%let dir = H:\Chao_022216\try_fine_gray\output;

data bmt;
  set bmt;
  call streaminit(123);      
  u = rand("Uniform");      
  output;
run;

data bmt;
  set bmt;
   if u>=0.7 then sex=1; else if u<0.5 then sex=0;
   if u>=0.4 then race='white'; else race='AA';
   Dftime= t / 365;
   
   label sex    = 'gender'
         group  = 'patient group'
         Dftime = 'Survival(Years)'
	;
 format Group DiseaseGroup. SEX SEX.;
run;

*title 'Fig 1 Plot of CIF for different patient groups ';
%Plots_CIF(dsn= bmt,
          grplist= group race ,
          time_event=Dftime,
          censor=status,
          eventcode=1,
          xaxisvalue=0 1 2 3 4 5 6  ,
          yaxisvalue= , 
          timelist=  2 4 6   ,
          NatRisk = T,
          units= Years,
		  filename= CIF PLOTS,
		  outpath= &dir.\ ,
          debug=T)
