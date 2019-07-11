/**********************************************************************************************
Macro: %FineGray_sel
Created Date/Author/Contact: Sep 21, 2016/Chao Zhang, Yuan Liu, and Yaqi Jia;
Last Update Date/Person/Contact: Oct 2016/Yuan Liu
Current Version: V4
Working Environment: SAS 9.4 English version

Contact: Dr. Yuan Liu/email: yliu31@emory.edu

Purpose:  To conduct multivariable analysis for competing risk model. The proportional 
subdistribution hazards model as proposed by Fine and Gray (1999) was used followed by 
backward elimination. 

Notes: The model runs using PROC PHREG.  The final list of variables selected is saved in two
macro variables: &_finalcvar and &_finalvar, and is also written into the log.The macro
“MUTLIPLE_PHREG V21.sas” or later is also required.

Parameters: 

DSN            The name of the data set to be analyzed.

EVENT          Name of time to event outcome variable.   

CENSOR         Name of censoring indicator variable.  Values of 0 indicate censored. 

EVENT_CODE     The value in CENSOR that indicate event of interest, and this value will 
			   appear EVENTCODE= option.

VAR            The list of variables on interest in the initial model that would be 
			   eliminated during the backward selection procedure separated by spaces. 
			   The order of variables in this list will be preserved in the final report.

CVAR           The list of categorical variables that are in VAR. If need to 
			   change the reference group, you can follow each variable name by (DESC) or 
			   by (ref = “Ref level in formatted value”) where needed and separate terms by *.
			   See code example. 

INC            Number of variables to include in the model (optional).  The first n variables 
               in the VAR parameter will be included in every model.  The default value is 0.   

ALPHA          the criterion for retaining a variable in the backward elimination procedure.

TYPE3          Set to F to suppress type III p-values from being reported in the table 
               (optional).  The default value is T.  This only has an effect if REPORT = T.

ID             Variable to be used in the ID statement in PHREG (optional).  Refer to SAS
			   Help and Documentation for proper use of this option.


CLNUM          Set to T if you want to see the number of observations for each level of covariates. The default is T.

ORIENTATION	   orientation of the output Word table. Default is portrait, can be changed to landscape.

REPORT 		   Set it to T if a results summary table is desired. Otherwise check Log for variable 
			   selected by the backward elimination.

FILENAME       File name for output table.  This is necessary if report=T.

OUTPATH        File path for output table to be stored.  This is necessary if report=T.

DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
* For more details, please see the related documentation
**********************************************************************************************/

%macro FINEGRAY_SEL(dsn=,  event=,  censor=,  event_code=1, var=,  cVar=, inc= 0, alpha=.2,
         type3=T, id=,  clnum=T, ORIENTATION = portrait, report =T, outpath=, filename=, debug=F) ;
     
   %local cnt cnt_ pos pos_ L1 L2 clist_ _cnt_ cvar_cnt i cnt_bef  cvar_ref cnt2 len2 cntw pos2 pos3 bef cnt_bef cvar_ nn;

   /* Capitalize */
    
    %let debug      = %UPCASE(&debug);
	%let type3      = %UPCASE(&type3);
	%let report     = %UPCASE(&report);
    
   /* Macros for final variable lists */
   %global _finalvar _finalcvar;

  /* Save current options */
   PROC OPTSAVE out=_options;   RUN;

  /* Count number of variables in var */
    %let var_cnt = %sysfunc(countw(&var,' '));



  /* It could be longer if there are interaction terms */
   %let length = 64; 


   /* Time Should not be missing value*/
   
   %if &event = %STR()  %then %do;
      %put ERROR: EVENT should be specified.; 
      %goto exit;
   %end;



/*Check CVAR whether in the proper format, such as separate by * if reference level is specified*/
  /* Count number of class variables */

  %if %superq(cvar) = %str( ) %then %let cvar_cnt = 0;
  %if %superq(cvar) ~= %str( ) %then %do;
  		
        %let cnt = %qsysfunc(countw(%superq(cvar),'*'));
		%let cnt_= %qsysfunc(countw(%superq(cvar),' '));

        %let pos=%qsysfunc(countc(%superq(cvar), '(' ));
        %let pos_ = %qsysfunc(countc(%superq(cvar), '*' ));

		%if %superq(pos) = 1 %then %do;
		%let L1=%qsysfunc(findc(%superq(cvar),  '(' ));
	    %let L2=%qsysfunc(findc(%superq(cvar),  ')' ));
		%let clist_=%qsysfunc(SUBSTR(%superq(cvar), %superq(L1),%eval(%superq(L2) - %superq(L1) +1)));
	    %let _cnt_ = %qsysfunc(countw(%superq(clist_),' '));
        %end;

        %if %superq(cnt) = 1 and  %superq(cnt_) = 1 %then %let cvar_cnt = 1;
		%else %if %superq(pos_) > 0 %then %let cvar_cnt = %superq(cnt);
		%else %if %superq(pos_) = 0 and %superq(pos) = 0 %then %let cvar_cnt = %superq(cnt_);
		%else %if (%superq(pos) = 1 and %superq(pos_) = 0) and %superq(cnt_) = %superq(_cnt_) %then  %let cvar_cnt = 1;
		%else %do; %put ERROR: The categorical variables in CVAR should be separated by * if you specify the reference level.; 
          		   %goto exit; %end;

  %end;

/*Build up categorical variable list without reference level specified*/
  %let cvarlist =; 

  %if &cvar_cnt = 1 and &pos = 1 %then %let cvarlist = %sysfunc(SUBSTR(%superq(cvar), 1,%eval(&L1-1)));
  %else %if &cvar_cnt = 1 and &pos = 0 %then %let cvarlist = %superq(cvar);
  %if &cvar_cnt > 1  and %superq(pos_) = 0 %then  %let cvarlist = %superq(cvar) ; 

  %if &cvar_cnt > 1  and %superq(pos_) > 0 %then %do; 
 			
      %do i = 1 %to &cvar_cnt; 

				%let cvar_ref = %SCAN(&cvar, &i, '*'); 

				%let cnt2=%qsysfunc(countc(&cvar_ref, '(' ));
				
				%let len2=%LENGTH(&cvar_ref);
				%let cntw=%qsysfunc(countw(&cvar_ref,' '));

				/*%put cvar_ref = &cvar_ref cnt2= &cnt2  len2=&len2 cntw =&cntw ;*/

				%if %superq(cnt2) = 1 %then %do;

	                %let pos2=%sysfunc(findc(&cvar_ref, '(' ));
					%let pos3=%sysfunc(findc(&cvar_ref, ')' ));

					%let bef = %sysfunc(SUBSTR(&cvar_ref, 1,%eval(&pos2-1)));
			        %let cnt_bef = %qsysfunc(countw(&bef,' '));

					%if %superq(cnt_bef) = 1 and  %superq(pos3) = %superq(len2) %then %let cvar_=&bef;
					%else %do; %put ERROR: The categorical variables in CVAR should be separated by * if you specify the reference level.; 
          		                %goto exit;%end;
					%end;
					
                %if %superq(cnt2) = 0 and %superq(cntw) = 1 %then %LET cvar_ = &cvar_ref;

				%if %superq(cnt2) > 1 or ( %superq(cntw) > 1 and %superq(cnt2) = 0) %then %do; %put ERROR: The categorical variables in CVAR should be separated by * if you specify the reference level.; 
          		%goto exit;  %end;	 

			%let cvarlist = %superq(cvarlist) %superq(cvar_); 
 %end;
		 
 %end;

 /* Make sure that CVAR are either in VAR or FORCEINVAR*/
 %put &var   ;
 %put &cvarlist;

   %do i = 1 %to &cvar_cnt;
 		%let _check_ = 1;
     %do j = 1 %to &var_cnt;
         %if %upcase(%SCAN(&var, &j)) = %upcase(%SCAN(&cvarlist, &i)) %then %let _check_ = 0;;
     %end;
   %end;

   %if &_check_ = 1 %then %do; 
   %put ERROR: The varaible appears in CVAR should also be in VAR.; 
   %goto exit;  %end;


  /* Check for variable names that are longer than the default and increase length as needed */
   %DO i = 1 %to &cvar_cnt; 
      %IF %LENGTH(%SCAN(&var, &i,' ')) > &length %then %let length = %LENGTH(%SCAN(&var, &i,' ');
   %END;


   /*Fine&Gray model with backwards elimination */ 

   %let remlist =;

  %do %while (&var_cnt > 0);

	%if &cvar_cnt > 0 %then %let  cvar_ref = %sysfunc(tranwrd(%superq(cvar),*,%str()));
    
   ods select none;
   ODS OUTPUT  NObs = numobs ClassLevelFreq =clfreq ParameterEstimates = MLEC ModelInfo=modelinf ModelANOVA = type3  ;
   proc phreg data=&dsn simple namelen=&length %if &id ~= %STR() %then %do; COVS(AGGREGATE) %end;;
        %if &cvar_cnt > 0 %then %do; 
        class &cvar_ref /order=internal param=glm;%end;
        model &event*&censor(0) = &var / eventcode=&event_code rl;
		%if &id ~= %STR() %then %do;
            id &id;
         %end;
    run;
	ods select all;

	/*Identifiy next candidate to be removed*/

     proc sql noprint; select count(*) into :nn from type3; quit;

	 %if &nn = &inc %then %goto exit2;
	 %else %do;

	* in Type3 data remove variables that should not be fixed in model;
    	Data typ3_n; set type3; if _N_ > &inc;run;
	* for the rest variables, see which one has the largest p value and compare it to alpha; 

	 proc sql noprint; 
	    select Effect into: var_ separated by ' ' from typ3_n ;quit; * keep the original order;

     proc sort data=typ3_n; by desending ProbChiSq;run;
     proc sql noprint;
	   select ProbChiSq into: rplist separated by ' ' from typ3_n ;
       select Effect into: rvlist separated by ' ' from typ3_n ;
     quit;


	%let premove_p = %scan(&rplist,1, ' '); 
	%let premove = %scan(&rvlist,1, ' ');

	%put &premove_p &premove;

	* if the current largest pvalue is less than alpha, then stop, else update &var and &cvar;
	%if &premove_p < &alpha %then %goto exit2;

	%else %do;
	     %let remlist = &remlist &premove;
		 /* update &var list and keep the same input order;*/
	     %let var = %sysfunc(tranwrd(%upcase(&var),%upcase(&premove),%str())); 

		 /* update &cvar;*/
		 %if &cvar_cnt >  0 %then %do;* update &cvar;

		    %let pos_=%qsysfunc(findc(%superq(cvar), '*' ));
			

		    %let newcvar = ;
		 	%do i = 1 %to &cvar_cnt;
/*				%put pos_ = &pos_ cvar= &cvar;*/
  				%if &pos_ > 0 or &cvar_cnt = 1 %then %let scanvar = %SCAN(&cvar, &i, '*');
				%else %let scanvar = %SCAN(&cvar, &i, ' ');

				%let pos = %qsysfunc(findc(%superq(scanvar), '(' ));

				%if &pos > 0 %then %let trimvar=%qsysfunc(SUBSTR(%superq(scanvar), 1,%eval(&pos-1)));
				%else %let trimvar = &scanvar;

/*				%put trimvar = &trimvar scanvar=&scanvar premove=&premove;*/

				%if %upcase(&trimvar) = %upcase(&premove) %then %let newcvar = &newcvar; 
				%else  %let newcvar = &newcvar*&scanvar; 

			%end;
			%let cvar = &newcvar;
		 %end;
	%end;

    %let var_cnt = %sysfunc(countw(&var,' ')); * update number of variable in &var; 
    %let cvar_cnt = %sysfunc(countw(&cvar,'*'));  * update number of variable in &cvar; 

%end;/*%do %if (&var_cnt > 0);*/

   /* Save in macro variables */
   %let _finalvar = &var;
   %let _finalcvar = &cvar;

%end; /*end of 	 %if &nn = &inc %then %goto exit2;  %else %do;*/

%exit2:

/* Produce report of final model */
   %if &report = T %then %do;
      %let forceInVar = .;
	  %if &inc >0 %then %do;

	    Data typ3_n; set type3; if _N_ <= &inc;run;
		proc sql noprint;
		select Effect into :inforcevar separated by ' '
		from Typ3_n;

        PROC CONTENTS DATA = &dsn (keep=&inforcevar) out=cont noprint; RUN;

		 DATA cont;
            set cont end=last;
             if label = ' ' then label = name;
             if _N_ ~= 1 and last then label = 'and ' || label;
         RUN;

         PROC SQL noprint;
            select label into :inLab separated by ', '
             from cont;
         QUIT;
		 %end;
		 %else %let inLab = None;

	     
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

           %let foottext = The following variables were forced in the model: &inLab . The following variables were removed from the model: &remLab..;

      %end;
      %else %let foottext = The following variables were forced in the model: &inLab . No variables were removed from the model.;

      %MULTIPLE_PHREG(        
          OUTPATH = &outpath, 
          FNAME = &filename,
		  ORIENTATION = &ORIENTATION,
		  clnum=&clnum,
          TYPE3 = &type3,
          debug= &debug,
          FOOTNOTE="** Backward selection with an alpha level of removal of &alpha was used.  &foottext"
        );

  %end;


   %if &debug = F %then %do;
      *--- DELETE ALL TEMPORARY DATASETS that were created; 
      proc datasets lib=work memtype=data noprint;  
           delete %if cvar ~= %STR() %then %do; freq  %end; mergedata mlec modelinf numobs type3 typ3_n keptlist Cont
                  result1 result2 result2_ result3 Ori_list Removed _options _cont ;
       quit;
    %end;
   
/* Final variables selected (get rid of double spaces */
   %put Categorical variables selected: &_finalcvar;
   %put All variables selected: &_finalvar;

%exit:
%mend finegray_sel;

**example **;
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
   Dftime= t / 365.25;
   
   label sex = 'gender'
         group = 'patient group'
	;
 format Group DiseaseGroup. SEX SEX.;
run;

%include "H:\Macros\MULTIPLE_PHREG V21.sas";
%let dir = H:\Chao_022216\try_fine_gray\output;

Title "Table 3 Multivariable Survival Analysis for Fine and Gray’s Model";
%finegray_sel (dsn=bmt,  
       event=t,    
       censor=Status,  
       var= race Group  logWaitTime sex  , 
       cVar= Group(desc)*  race *sex,
       event_code=1,
	   inc = 1,
	   alpha=0.2,
	   Type3=t,
	   debug=F,
	   outpath =  &dir.\ ,
	   filename = Multivariable Fine and Gray) ;
Title;

Title "Table 3 Multivariable Survival Analysis for Fine and Gray’s Model";
%finegray_sel (dsn=bmt,  
       event=t,    
       censor=Status,  
       var= race Group  logWaitTime sex  , 
       cVar= Group(REF="AML-High Risk")*  race *sex,
       event_code=1,
	   inc = 1,
	   alpha=0.2,
	   Type3=t,
	   debug=F,
	   outpath =  &dir.\ ,
	   filename = Multivariable Fine and Gray) ;
Title;
