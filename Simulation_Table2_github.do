****************************
* time-fixed exposure      *
*    long-term effect      *
* Estimation of ATE and ATT*
* using alternative methods*
*                          *
* as in Table 2           *
*                          *
* BDS  26/01/26            *
****************************

*************************************************************************************
/* SELECT THE FOLDER WHERE DATA and Do files ARE HELD */
cd 		"C:\Users\sejjbld\OneDrive - University College London\_GRANTS\Ruth Gilbert\HOPE\WP3&4\papers\IJPDS Notes\data_analyses"

/* SAVE RESULTS */
cap log close
log using log\Simulation_Table2_github.log, replace
*************************************************************************************
/* READ and DESCRIBE THE DATA */
use 		"data\Simulation_timevar_SEND_Study",clear
describe

/* RESHAPE THE DATA FROM LONG TO WIDE FORMAT TO FACILITATE SOME CALCULATIONS */
reshape 	wide Hosp SEN Y  Hosp0_ Hosp1_ SEN0_ SEN1_ den Y0_ Y1_   , i(id) j(t)

*************************************************************************************
/* 
   NEW VARIABLES

   RENAMING INTERVENTION, OUTCOME, TIME-VARYING CONFOUNDER, AND DENOMINATOR WITH GENERIC NAMES HELPS UNDERSTANDING THE CODE 
*/
 
*THE  INTERVENTION  AT TIME 1 (BASELINE)
	rename SEN1 A

*THE  TIME-VARYING CONFOUNDER AT TIME 1 (BASELINE)
	rename Hosp1 H
	
*THE OUTCOME: FOR TABLE 2 THE OUTCOME IS THE CUMULATIVE RATE UP TO END OF YEAR 4 
*observed outcome
	gen double 	Y=Y2+Y3+Y4
*observed rate
	gen double 	D=den2+den3+den4
	gen double 	rate=Y/D

*mean denominators (to derive some of the estimates)
	su D 
	scalar mean=r(mean)
	su D if A==1
	scalar mean_att=r(mean)
	


*************************************************************************************
/*ALTERNATIVE ESTIMATION METHODS */
/* (1) G-COMPUTATION */
*************************************************************************************


/*INCORRECT outcome MODEL 1 : 
USING "margins" */

* the (incorrect) baseline confounders 
	vl create base_conf=(male white idacibin eyfsp_bin H)

*ATE
qui poisson 		Y $base_conf i.A , exp(D) nolog robust
margins 			A, post
*ATE_ratio, count and rate
nlcom 				_b[1.A]/ _b[0bn.A]
*ATE rate, dif
nlcom 				(_b[1.A]/mean)- (_b[0bn.A]/mean)
*ATT
qui poisson 		Y $base_conf i.A , exp(D) nolog robust 
margins 			A, post vce(unconditional)  subpop(A)
*ATT_ratio, count and rate
nlcom 				_b[1.A]/ _b[0bn.A]
*ATT rate, dif
nlcom 				(_b[1.A]/mean_att)- (_b[0bn.A]/mean_att)

/*INCORRECT outcome MODEL 2:  
USING "margins" */

* the (incorrect) baseline confounders 
	vl modify base_conf=(male white idacibin H)

*ATE
qui poisson 		Y $base_conf i.eyfsp_bin##i.A , exp(D) nolog robust
margins 			A, post
*ATE_ratio, count and rate
nlcom 				_b[1.A]/ _b[0bn.A]
*ATE rate, dif
nlcom 				(_b[1.A]/mean)- (_b[0bn.A]/mean)
*ATT
qui poisson 		Y $base_conf i.eyfsp_bin##i.A , exp(D) nolog robust
margins 			A , post vce(unconditional)  subpop(A)
*ATT_ratio, count and rate
nlcom 				_b[1.A]/ _b[0bn.A]
*ATT dif, rate
nlcom 				(_b[1.A]/mean_att)- (_b[0bn.A]/mean_att)


/*CORRECT MODEL 3:  
USING "margins" */

*all interactions included
*ATE
#delimit ;
qui poisson 		Y i.(male white idacibin eyfsp_bin H)##i.A 
							i.(male white idacibin eyfsp_bin)##i.H##i.A 
							i.(male white idacibin)##i.eyfsp_bin##i.A  
							i.(male white)##i.idacibin##i.A 
							i.male##i.white##i.A, 
							exp(D) nolog robust
							;
#delimit cr
margins 			A, post
*ATE_ratio, count and rate
nlcom 				_b[1.A]/ _b[0bn.A]
*ATE rate, dif
nlcom 				(_b[1.A]/mean)- (_b[0bn.A]/mean)
*ATT
#delimit ;
qui poisson 		Y i.(male white idacibin eyfsp_bin H)##i.A 
							i.(male white idacibin eyfsp_bin)##i.H##i.A 
							i.(male white idacibin)##i.eyfsp_bin##i.A  
							i.(male white)##i.idacibin##i.A 
							i.male##i.white##i.A, 
							exp(D) nolog robust
							;
#delimit cr
margins 			A , post vce(unconditional)  subpop(A)
*ATT_ratio, count and rate
nlcom 				_b[1.A]/ _b[0bn.A]
*ATT dif, rate
nlcom 				(_b[1.A]/mean_att)- (_b[0bn.A]/mean_att)


/*CORRECT MODEL 3:  
USING own ado and  bootstraps  where the baseline confounders are: male white idacibin eyfsp_bin H and all their interactions 

The ADO parses the names of: outcome, intervention, denominator and time-varying confounder (Y A D H)

The time varying confounder is parsed because the same ado can be used to produce the estimates for short term effects (see "Simulations_Table3_github.do")

The estimated effects are the ATE and ATT on the ratio and difference scale
*/

bootstrap ATE_rate_ratio=r(ATE_rate_ratio)  ATT_rate_ratio=r(ATT_rate_ratio) ATE_rate_dif=r(ATE_rate_dif) ATT_rate_dif=r(ATT_rate_dif)  , ///
         reps(1000) nodrop seed(1212): SEND_gcomp Y A D H



*************************************************************************************
/*ALTERNATIVE ESTIMATION METHODS */
/* (2) IPW */
**********************************************************************************

/*USING own ado and  bootstraps (teffects does not support Poisson distributions for the rate)

here the baseline confounders for the PS are: male white idacibin eyfsp_bin H (ie no interaction)

The ADO parses the names of: outcome, intervention, denominator and time-varying confounder (Y A D H)

The time varying confounder is parsed because the same ado can be used to produce the estimates for short term effects (see "Simulations_Table3_github.do")

The estimated effects are the ATE and ATT on the ratio and difference scale

*/

/*INCORRECT MODEL 1:  */
* the (incorrect) baseline confounders 
	vl modify base_conf=(male white idacibin eyfsp_bin)

*CHECK THE OVERLAP
logit A $base_conf H, nolog
cap drop den
predict den
twoway (kdensity den if A==1, lcol(red)) (kdensity den if A==0,lcol(blue)), name(PS_incorrect, replace) xtitle(PS) legend(order(1 "Exposed" 2 "unexoposed"))
*NOTE: THE OVERLAP IS PARTICULARLY POOR BECAUSE ALL CONFOUNDERS ARE BINARY
*/

bootstrap  ATE_rr=r(ATE_rr)  ATE_rd=r(ATE_rd)  ATT_rr=r(ATT_rr)  ATT_rd=r(ATT_rd)  , reps(1000) nodrop seed(1212) : SEND_IPW Y A D H



/*CORRECT MODEL 2:  USING own ado and  bootstraps*/
* the (correct) baseline confounders 
gen 		inter=(idacibin ==1)*(eyfsp_bin==0)
vl modify 	base_conf=(male white idacibin eyfsp_bin inter)

*CHECK THE OVERLAP
logit A $base_conf H, nolog
cap drop den
predict den
twoway (kdensity den if A==1, lcol(red)) (kdensity den if A==0,lcol(blue)), name(PS_correct, replace) xtitle(PS) legend(order(1 "Exposed" 2 "unexoposed"))
*NOTE: THE OVERLAP IS PARTICULARLY POOR BECAUSE ALL CONFOUNDERS ARE BINARY
*/

*with  bootstraps
bootstrap  ATE_rr=r(ATE_rr)  ATE_rd=r(ATE_rd)  ATT_rr=r(ATT_rr)  ATT_rd=r(ATT_rd)  , reps(1000) nodrop seed(1212) : SEND_IPW Y A D H



**********************************************************************************
******* AIPW 
**********************************************************************************

*********model 1* incorrect PS and Y
vl modify 	base_conf=(male white idacibin eyfsp_bin  H)

*ATE
teffects aipw (rate $base_conf,poisson)   (A $base_conf)
*ATE ratio
nlcom ( _b[POmean:0.A]+_b[ATE:r1vs0.A])/_b[POmean:0.A]
*ATT
teffects aipw (rate $base_conf,poisson)   (A  $base_conf), atet
*ATT ratio
nlcom ( _b[POmean:0.A]+_b[ATET:r1vs0.A])/_b[POmean:0.A]



*model 2, correct PS, incorrect Y
vl modify 	base_conf=(male white idacibin eyfsp_bin H inter)
*ATE
teffects aipw (rate $base_conf, poisson)   (A  $base_conf)
*ATE ratio
nlcom ( _b[POmean:0.A]+_b[ATE:r1vs0.A])/_b[POmean:0.A]
*ATT
teffects aipw (rate $base_conf, poisson)   (A  $base_conf), atet
*ATT ratio
nlcom ( _b[POmean:0.A]+_b[ATET:r1vs0.A])/_b[POmean:0.A]
*/


*model 3, incorrect PS, right Y 
*ATE
teffects aipw (rate i.(male white idacibin eyfsp_bin H) i.(male white idacibin i.eyfsp_bin)##i.H i.(male white idacibin)##i.eyfsp_bin  i.(male white)##i.idacibin  i.male##i.white, poisson)   (A  $base_conf)
*ATE ratio
nlcom ( _b[POmean:0.A]+_b[ATE:r1vs0.A])/_b[POmean:0.A]

*ATT
teffects aipw (rate i.(male white idacibin eyfsp_bin H) i.(male white idacibin i.eyfsp_bin)##i.H i.(male white idacibin)##i.eyfsp_bin  i.(male white)##i.idacibin  i.male##i.white,poisson )   (A  $base_conf), atet
*ATT ratio
nlcom ( _b[POmean:0.A]+_b[ATET:r1vs0.A])/_b[POmean:0.A]
*/

*model 4, correct PS, and correct Y
*ATE
teffects aipw (rate i.(male white idacibin eyfsp_bin H) i.(male white idacibin i.eyfsp_bin)##i.H i.(male white idacibin)##i.eyfsp_bin  i.(male white)##i.idacibin  i.male##i.white,poisson)   (A  $base_conf)
*ATE ratio
nlcom ( _b[POmean:0.A]+_b[ATE:r1vs0.A])/_b[POmean:0.A]

*ATT
teffects aipw (rate i.(male white idacibin eyfsp_bin H) i.(male white idacibin i.eyfsp_bin)##i.H i.(male white idacibin)##i.eyfsp_bin  i.(male white)##i.idacibin  i.male##i.white , poisson)   (A  $base_conf), atet
*ATT ratio
nlcom ( _b[POmean:0.A]+_b[ATET:r1vs0.A])/_b[POmean:0.A]

**********************************************************************************
*******Wald-type IV for rates 
**********************************************************************************
*
*IV_bin
ivregress 2sls rate  (A= region )
ivregress 2sls rate $base_conf (A= region )
ivregress gmm  rate $base_conf  (A= region )


**********************************************************************************
*SIMPLE analyses
**********************************************************************************
*wrong as it ignores interactions
poisson 	Y A male white idacibin eyfsp_bin H  , exp(D) nolog irr
*ATE,ratio, rate
nlcom 		exp(_b[A])
*ATE,dif,  rate
nlcom 		exp(_b[_cons]+_b[A]) -exp(_b[_cons])
*/




exit

*short-cut IPW

/*INCORRECT MODEL 1:  : USING "teffects" */
/*
vl modify 	base_conf=(male white idacibin eyfsp_bin)
*ATE
teffects 	ipw (rate)   (A  $base_conf H )
*ATE ratio
nlcom 		( _b[POmean:0.A]+_b[ATE:r1vs0.A])/_b[POmean:0.A]
*ATT
teffects 	ipw (rate)   (A $base_conf H), atet
*ATT ratio
nlcom 		( _b[POmean:0.A]+_b[ATET:r1vs0.A])/_b[POmean:0.A]
*/

/*CORRECT MODEL 2:  : USING "teffects" */
/*
vl modify 	base_conf=(male white idacibin eyfsp_bin inter)
*ATE
teffects ipw (rate)   (A  $base_conf H)
*ATE ratio
nlcom ( _b[POmean:0.A]+_b[ATE:r1vs0.A])/_b[POmean:0.A]
*ATT
teffects ipw (rate)   (A  $base_conf H), atet
*ATT ratio
nlcom ( _b[POmean:0.A]+_b[ATET:r1vs0.A])/_b[POmean:0.A]
*/

*
