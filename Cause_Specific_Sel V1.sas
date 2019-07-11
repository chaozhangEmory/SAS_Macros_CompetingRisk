/***************************************************************************************************
****************************************************************************************************;
Macro Name: CauseSpecific_Sel
Created Date/Author/Contact: Jan 2019/  Chao Zhang,  Yuan Liu, and Yaqi Jia
Current Version: V1
Working Environment: SAS 9.4 English version  

Contact Information: chao.zhang2@emory.edu

Purpose:  :  To conduct multivariate backward selection Cause-Specific Hazard model on competing risk
survival data. For a particular cause of interest, we treat all the competing events as censored 
observations in the analysis.  The risk set at each time includes only those subjects who have not 
failed from competing events or are truly censored.    


References:
    1 SAS Institute Inc, 2017 SAS/STAT® 14.3 User’s Guide. Cary, NC: SAS Institute Inc. User’s Guide the Phreg Procedure.
    2 Changbin Guo and Ying So (2018), Cause-Specific Analysis of Competing Risks Using the Phreg Procedure. SAS2159-2018 SAS Institute Inc.
    3 PHREG_SEL V23 created by Dana Nickleach

Parameters: 

DSN            The name of the data set to be analyzed(required).   
 
EVENT          Name of time to event outcome variable (required).                                                             

CENSOR         Name of censoring indicator variable.  Values of 0 indicate censored observations(required).

CENSORED_VALUE The value in CENSOR that were specified the Values of Censored(required). For example, 0 = Censored, 1= Relapse, 2= Death,
               if relapse is the event of interest, type 0 and type 2 were specified censored in Cause-Specified Hazard Model 
               of competing-Risk survival data.

VAR            The list of variables on interest in the initial model that would be eliminated during the backward selection 
               procedure separated by spaces. The order of variables in this list will be preserved in the final report(required).

CVAR           The list of categorical variables that are in VAR(required). If need to change the reference group, you can 
               follow each variable name by (DESC) where needed and separate terms by *.  See code example. 

INC            Number of variables to include in the model. The first n variables in the var parameter will be included 
               in every model(optional).  The default value is 0. 

SLSTAY         The significance level for removing variables from the model (optional).  The default value is 0.2.

ID             Variable to be used in the ID statement in PHREG (optional).  Refer to SAS Help and Documentation for proper use of this 
               option. ID and COVSAGG will be used together in the model to compute the robust sandwich covariance matrix estimate for 
               cluster data. The ID statement identifies the variable that represents the clusters.  If observations have more than 
               one record in the data file, the number of observations reported in the table footer will be the number of records, 
               not unique observations. 

ORIENTATION    Orientation of the output Word table. Default is portrait, can be changed to landscape.

OUTPATH        Path for output table to be stored(required).

FNAME          File name for output table(required).

DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted in
               debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

******************************************************************************************************/

Options macrogen symbolgen mlogic mprint mfile;
%macro CauseSpecific_sel(dsn=, event=, censor=, censored_value=, var=,cvar=,inc=0,slstay=.2,id=,covsagg=f, ORIENTATION = portrait,outpath=,filename=,debug=F);
    /* Save current options */
    PROC OPTSAVE out=_options;   RUN;

    /* Upper case T/F */
    %let debug   = %UPCASE(&debug);
    %let var     = %UPCASE(&var);
    %let cvar    = %UPCASE(&cvar);
	%let covsagg = %UPCASE(&covsagg);
	/* Initialize error variable */
    %LET __Macro_Err = 0;

   /* Time Should not be missing value*/
   
    %if &censor = %STR()  %then %do;
       %put ERROR: censor should be specified.; 
       %let __Macro_Err=1;
    %end;
  
	 %if %INDEX(&cvar,*) > 0 OR %INDEX(&cvar,(DESC)) > 0 %then %let cvar_cnt = %sysfunc(countw(&cvar,*));
     %else %let cvar_cnt = %sysfunc(countw(&cvar,' ')); 

	 
      %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
     %end;

     ODS SELECT NONE;
      ODS OUTPUT   ModelBuildingSummary(nowarn)=removed NObs = numobs
          ParameterEstimates = MLEC ModelInfo=modelinf %if &cvar_cnt ~=0 %then %do; ModelANOVA = type3 ClassLevelFreq=clfreq %end;;

          proc phreg data=&dsn   %if &covsagg=T %then %do; COVS(AGGREGATE) %end;simple;
		   
             class %if %sysevalf(%superq(CVAR)~=,boolean) %then %do; 
               %sysfunc(TRANSLATE(&cvar," ","*")) %end;/ order=internal param=glm;

               model &event*&censor(&censored_value) = &var/rl selection=backward include=&inc 
               slstay=&slstay hierarchy=single;
        

         %if &id ~= %STR() %then %do;
            id &id;
         %end;

      run;
	  ODS SELECT ALL;

  
       /* Save order */
       DATA MLEC2;
          set MLEC;
          order = _n_;
       RUN;
       
	   /* Check to see if data set containing a removed variable exists */
      /* If no variables were removed then the data set will not be created and the selection */
      /* process is done */
      %let continue = %sysfunc(exist(removed));

      %if &continue = 1 %then %do;
         /* Get name of variable removed */
         /* Need SCAN because variable label is included after variable name in the same column */
         PROC SQL noprint;
            select UPCASE(SCAN(EffectRemoved,1,' '))
            into :remove
            from removed;
         quit;

         %let remlist =  &remove;

       /* If any variables were removed from the model create a list */
      %if &remlist ~= %STR() %then %do;
         /* Get variable labels for footnote */
         PROC CONTENTS DATA = &dsn (keep=&remlist) out=cont noprint;
         RUN;

         DATA cont;
            set cont end=last;
             if label = ' ' then label = name;
             if _N_ ~= 1 and last then label = 'and ' || label;
         RUN;

         PROC SQL noprint;
            select label into :remLab separated by ', '
             from cont;
         QUIT;

           %let foottext = The following variables were removed from the model: &remLab..;

      %end;
      %else %let foottext = No variables were removed from the model.;
      %end;

	  
/* Get dataset name and outcome name and use them in the output*/
   PROC SQL noprint;
      select cvalue into :event
      from modelinf
      where description = 'Dependent Variable';
   QUIT;

   /* Get outcome label */
   data _NULL_; 
      set &dsn (obs=1); 
      call symput('eventvv', VLABEL(&event));
   run;

	  ** get number for continuous vars **;
      proc sql noprint;
        select NObsRead into :numread from numobs; 
        select NObsUsed into :numused from numobs; 
      quit;

    %if &cvar_cnt ~=0 %then %do;
	  **merge type3 and mlec datasets **;
	  data type3;
	     set type3;
		rename ProbChiSq = type3_p ;
		
	  run;

       proc sql;
	     create table MergeData AS 
	     select A.* , B.effect, B.type3_p  
		 from MLEC2 A, type3 B
		 where MLEC2.Parameter = type3.effect;
       QUIT;
      

     **add frequency to mergedata **;
        %let continue2 = %sysfunc(exist(Clfreq));

        %if &continue2 = 1 %then %do;

         data Clfreq(rename=(value=level));
            length parameter $100;
            set Clfreq;
            retain;
            if class ~="" then parameter = upcase(class);
            else if class ="" then do;
            retain class;
            end;
            drop class control_var;
			format freq 8.0 ;
         run; 

	   data mergedata;
	      length parameter $100;
		  set mergedata;
		  parameter  = UPCASE(parameter);
		  rename ClassVal0 = level;
	   run;

       proc sql;
	     create table mergedata_ as
		    select A.*, B.FREQ
			from mergedata A LEFT JOIN Clfreq B
			  on A.parameter = B.parameter AND A.level =B.level;
			   
		QUIT;


    proc sort data=mergedata_; by effect order; run;
    data mergedata_ ;
	  set mergedata_;
	  if first.effect  then p_type3 = type3_p ;
	  by effect order;
    run;
 
	 data mergedata_ ;
	    set  mergedata_ ;
        N = freq;
        if freq =. then N = &numused;
	  run;

	proc sort data=mergedata_; by order; run;

   %end;
 %end;
  %else %do;
	    data mergedata_;
		  set mlec2;
		  rename Parameter = effect;
		 run;
  %end;
 
DATA mergedata_;
      set MergeData_;
      /* HR and 95% CI */
      if HazardRatio ~= . then HR = TRIM(LEFT(PUT(HazardRatio,8.2))) || " (" || 
         TRIM(LEFT(PUT(HRLowerCL,8.2))) || "-" || 
         TRIM(LEFT(PUT(HRUpperCL,8.2))) || ")";
      else HR = '-';
      /* Number of observations read and used - save in data set */
      numused = &numused;
      numread = &numread;
	  parameter = upcase(parameter);
RUN;

	**catch type of interest **;
	
   proc freq data=&DSN;
     table status/out = freqCount;
   run;

   data freqCount;
     set freqCount;
	 where &censor ~ in (&censored_value);
   run;

    Proc sql noprint;
         select &censor  into:interest
         from freqCount;
	quit;
  %PUT interest: &interest;

	 /* part 3 create the word format table  */;

*---- table template -----;  
   ODS PATH WORK.TEMPLAT(UPDATE) SASUSR.TEMPLAT(UPDATE) SASHELP.TMPLMST(READ);

   PROC TEMPLATE;
      DEFINE STYLE STYLES.TABLES;
      NOTES "MY TABLE STYLE"; 
      PARENT=STYLES.MINIMAL;

       STYLE SYSTEMTITLE /FONT_SIZE = 12pt   FONT_FACE = "TIMES NEW ROMAN";

       STYLE HEADER /
           FONT_FACE = "TIMES NEW ROMAN"
            CELLPADDING=8
            JUST=C
            VJUST=C
            FONT_SIZE = 10pt
           FONT_WEIGHT = BOLD; 

     STYLE TABLE /
            FRAME=HSIDES           
           RULES=GROUP              
           CELLPADDING=6            
            CELLSPACING=6           
            JUST=C
            FONT_SIZE = 10pt
           BORDERWIDTH = 0.5pt;  

     STYLE DATAEMPHASIS /
           FONT_FACE = "TIMES NEW ROMAN"
           FONT_SIZE = 10pt
           FONT_WEIGHT = BOLD;

     STYLE DATA /
           FONT_FACE = "TIMES NEW ROMAN" 
           FONT_SIZE = 10pt;

     STYLE SYSTEMFOOTER /FONT_SIZE = 9pt FONT_FACE = "TIMES NEW ROMAN" JUST=C;
END;

RUN; 



   **** PRINT THE TABLE ************;
   OPTIONS ORIENTATION=&ORIENTATION MISSING = "-" NODATE;
   ODS rtf STYLE=TABLES file= "&OUTPATH.&FILENAME &SYSDATE..DOC"; 

   PROC REPORT DATA=mergedata_ HEADLINE HEADSKIP CENTER STYLE(REPORT)={JUST=CENTER} SPLIT='~' 
     nowd SPANROWS lS=256; 
      COLUMN effect %if &cvar_cnt ~=0 %then %do;  Level N %end;("&eventvv" '------------------------------------------'
          (HR  ProbChiSq  %if &cvar_cnt ~=0 %then %do; p_type3 %end;)); 
      DEFINE effect/order order=data "Covariate"  STYLE(COLUMN) = {JUST = C CellWidth=25%}; 
      %if &cvar_cnt ~=0 %then %do;
      DEFINE Level/DISPLAY "Level" STYLE(COLUMN) = {JUST =L CellWidth=20%}; 
      DEFINE N/DISPLAY "N" STYLE(COLUMN) = {JUST =C CellWidth=8%}FORMAT=8.0; 
      %end;
       %if &cvar_cnt ~=0 %then %do;
         DEFINE p_type3/display   MISSING "Global P-value" STYLE(COLUMN) = {JUST=C CellWidth=8%} 
               FORMAT=PVALUE8.3; 
       %end;
      DEFINE HR/DISPLAY "Hazard Ratio" STYLE(COLUMN) = {JUST = C CellWidth=15%}; 
      DEFINE ProbChiSq/DISPLAY "HR P-value" STYLE(COLUMN) = {JUST = C CellWidth=8%} FORMAT=PVALUE8.3; 

      COMPUTE ProbChiSq; 
         IF . < ProbChiSq <=0.05 THEN CALL DEFINE("ProbChiSq", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
       
      ENDCOMP; 

       
       %if &cvar_cnt ~=0 %then %do;
         COMPUTE p_type3; 
            IF . < p_type3 <0.05 THEN CALL DEFINE("p_type3", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
         ENDCOMP;   
       %end;

       compute after effect; line '';   endcomp;  

       compute after _page_;
         length text1 - text3 $200.;
         text1 = "*Number of observations in the original data set = %TRIM(&numread). Number of observations used = %TRIM(&numused).";
		 line @0 text1 $200.;
          
         text2 = "*Backward selection with an alpha level of removal of &slstay was used.  &foottext" ;
		 
          line @0 text2 $200.;

		 text3 = "*Cause-Specific Hazard Model:  Type of Censored= &censored_value;  Type of Interest= %TRIM(&interest)." ;
		 
         line @0 text3 $200.;
       ENDCOMP;
 
     RUN;
   ods rtf close;

   
   %if &debug = F %then %do;
      *--- DELETE ALL TEMPORARY DATASETS that were created; 
      proc datasets lib=work memtype=data noprint;  
           delete %if cvar ~= %STR() %then %do; CLfreq  %end; mergedata mlec modelinf numobs type3 REMOVED Cont Freqcount Mergedata_ Mlec2  ;
       quit;
    %end;
   
%mend CauseSpecific_sel;
/*
**example **;

data bmt;
  input Group T Status WaitTime @@;
  logWaittime=log(WaitTime);
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

%let dir = H:\Chao_paper_2019124\Informative_censoring\output;

proc format;
value DiseaseGroup 1='ALL' 2='AML-Low Risk' 3='AML-High Risk';
value sex 0= 'F' 1='M';
run;

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
        
   label sex = 'gender'
         group = 'patient group'
		 label T = ' disease-free survival time (days)' ;
	;
 format Group DiseaseGroup. SEX SEX.;
run;


title 'Table 1 Multivariable Cause-Specific Hazard Regression for Relapse' ;
%CauseSpecific_sel(dsn=bmt, 
                 event=t, censor=status, 
                 censored_value=0 2, 
                 var=Group logWaitTime sex,
                 cvar=Group   sex,
                 inc=3,
                 slstay=.3,
                 ORIENTATION = portrait, 
                 outpath= &dir.\,
                 filename=aaa,
                 debug=t);
title '';

title 'Table 2 Multivariable Cause-Specific Hazard Regression for Death' ;
%CauseSpecific_sel(dsn=bmt, 
                 event=t, censor=status, 
                 censored_value=0 1, 
                 var=Group logWaitTime sex,
                 cvar=Group  sex,
                 inc=0,
                 slstay=.3,
                 ORIENTATION = portrait, 
                 outpath= &dir.\,
                 filename=bbb,
                 debug=t);
title '';


