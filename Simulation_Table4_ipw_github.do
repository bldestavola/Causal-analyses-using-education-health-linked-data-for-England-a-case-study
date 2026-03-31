****************************
* time-varying exposure    *
*                          *
* Estimation of ATE and ATT*
* using g-computation      *
*                          *
* as in Table 4, gcomp     *
*                          *
* BDS  27/01/26            *
****************************

**********************Q1********************************************
* what is the effect of SEN at t=1 and t=2 on Cum Y to end of year 4?
********************************************************************
/* 
exposure SEN at t=1 and t=2 and t=3
outcome: cumulative Y up to end of year 3 (Y4)
confounders: male, white, idacibin, eyfsp Hosp at t=1
time-var conf: Hosp at t=2 and t=3
*/

*************************************************************************************
/* SELECT THE FOLDER WHERE DATA and Do files ARE HELD */
cd 		"C:\Users\sejjbld\OneDrive - University College London\_GRANTS\Ruth Gilbert\HOPE\WP3&4\papers\IJPDS case study\analyses"

/* SAVE RESULTS */
cap log close
log using log\Simulation_Table4_ipw_github_310326.log, replace
*************************************************************************************


/* READ and DESCRIBE THE DATA */
use 		"data\Simulation_timevar_SEND_Study",clear
describe

*drop POs to use same names for estimation
drop Y0_ Y1_ Y2_* Y3_* Y4_* SEN0_ SEN1_
	


**************************************************************************************
*create new and lagged values
xtset id t
gen SEN_lag1=L.SEN
gen Hosp_lag1=L.Hosp
gen Y_lag1=L.Y

gen SEN_lag2=L2.SEN
gen Hosp_lag2=L2.Hosp

gen SEN_lag3=L3.SEN
gen Hosp_lag3=L3.Hosp

sort id t
qui by id: gen Hosp_base=Hosp[1] 

replace den=. if t==1
replace SEN_lag2=0 if t<=2

gen inter_sen=(eyfsp_bin==0)*SEN_lag1
gen inter=(idacibin ==1)*(eyfsp_bin==0)

*needed for the bootstrap
gen newid=id

*create cumulative outcome and cumulative denominator
gen cumY=0
sort id t
qui by id: replace cumY=cumY[_n-1]+Y if t>1
qui by id: gen cumden=sum(den)
qui by id: gen cumden4=cumden[3]

/*
***********************************************************************************
*standard regression
***********************************************************************************
poisson Y  SEN_lag1 SEN_lag2 SEN_lag3  if  t==4, e(den) nolog
di "Risk difference"
nlcom exp(_b[_cons]+_b[SEN_lag1]+_b[SEN_lag2]+_b[SEN_lag3])-exp(_b[_cons])
di "Risk ratio"
nlcom exp(_b[SEN_lag1]+_b[SEN_lag2]+_b[SEN_lag3])

poisson Y  male  idacibin eyfsp_bin Hosp_lag1 Hosp_lag2 Hosp_lag3 SEN_lag1 SEN_lag2 SEN_lag3  inter if  t==4, e(den) nolog
di "Risk difference"
nlcom exp(_b[_cons]+_b[SEN_lag1]+_b[SEN_lag2]+_b[SEN_lag3])-exp(_b[_cons])
di "Risk ratio"
nlcom exp(_b[SEN_lag1]+_b[SEN_lag2]+_b[SEN_lag3])

*/


********************************************************************************
*Using IPW to estimate the sustained effct of SEN up to year 4
********************************************************************************

reshape long
replace den=. if t==1
replace SEN_lag2=0 if t<=2



*************************************************
* confounders 
  vl create base_conf=(male white Hosp_base idacibin eyfsp_bin inter)
  vl create tv_conf=(Hosp Hosp_lag1)
*************************************************


ta t	  


su SEN Y den
su id newid 
xtset newid t

isid newid t
duplicates report newid t
xtset, clear

Sustained_IPW_stab SEN Y den newid
bootstrap rate_000=r(rate_000) rate_111=r(rate_111) ATE_RR=r(ATE_RR) ATE_RD=r(ATE_RD), ///
       reps(1000)  seed(3009) cluster(id) idcluster(newid):                              ///
	   Sustained_IPW_stab SEN Y den newid
*nodrop

exit



















**IPW
*PS for receiving SEN in year 1 (as above)
cap drop wt*
cap drop ps*
logit SEN male white idacibin eyfsp Hosp if t==1, nolog
predict ps1 if t==1

*PS for receiving SEN in year 2
logit SEN male white idacibin eyfsp Hosp_lag1 Hosp SEN_lag1 if t==2, nolog
predict ps2 if t==2


*PS for receiving SEN in year 3
logit SEN male white idacibin eyfsp Hosp_lag1 Hosp SEN_lag1 if t==3, nolog
predict ps3 if t==3

gen ps=ps1 if t==1
replace ps=ps2 if t==2
replace ps=ps3 if t==3
gen wt=cond(SEN, 1/ps, 1/(1-ps))

tabstat ps* wt*, by(SEN) s(min max) 
br id t SEN ps* wt


twoway (kdensity ps if SEN==1 & t==1, lcol(red)) (kdensity ps if SEN==0 & t==1,lcol(blue)) (kdensity ps if SEN==1 & t==2, lcol(red) lpat(dash)) (kdensity ps if SEN==0 & t==2,lcol(blue) lpat(dash)) (kdensity ps if SEN==1 & t==3, lcol(red) lpat(dot)) (kdensity ps if SEN==0 & t==3,lcol(blue) lpat(dot)), name(byhand, replace) xtitle(PS)
* maybe SEN1 should not be included in the PS2 model? No no we need it to get the rigth results!


keep id t Hosp SEN Y cum* wt rate3_obs ps cumY3_obs cumden3  cumY3_sust_000 cumY3_sust_111 rate3_sust_000 rate3_sust_111 PO_dif_sust_rate3 PO_dif_sust_tot3  
br id t Hosp SEN Y rate* ps wt cumY3* cumden3 

reshape wide Hosp SEN Y   ps wt , i(id) j(t)

cap drop cum_wt
gen cum_wt=wt1*wt2*wt3
l id Hosp1 Hosp2 SEN1 SEN2 SEN3 cumY3_obs cumden3 wt1 wt2 wt3 cum_wt in 1/2, noobs sepby(id)


/*Gaussian- count
regress cumY3_obs SEN1 SEN2 SEN3 [pw=cum_wt]
di "base_rate"
nlcom (_b[_cons])
di "exposed rate"
nlcom (_b[_cons]+_b[SEN1]+_b[SEN2]+_b[SEN3])
di "ace_ipw (Risk difference)"
nlcom (_b[SEN1]+_b[SEN2]+_b[SEN3])
*/

*Poisson- rates
poisson cumY3_obs SEN1 SEN2 SEN3 [pw=cum_wt], exp(cumden3) nolog
di "base_rate"
nlcom exp(_b[_cons])
di "exposed rate"
nlcom exp(_b[_cons]+_b[SEN1]+_b[SEN2]+_b[SEN3])
di "ace_ipw (Risk difference)"
nlcom exp(_b[_cons]+_b[SEN1]+_b[SEN2]+_b[SEN3])-exp(_b[_cons])
di "ace_ipw (Risk ratio)"
nlcom exp(_b[SEN1]+_b[SEN2]+_b[SEN3])

ex
*Gaussian- rate
regress rate3_obs SEN1 SEN2 SEN3 [pw=cum_wt]
di "base_rate"
nlcom (_b[_cons])
di "exposed rate"
nlcom (_b[_cons]+_b[SEN1]+_b[SEN2]+_b[SEN3])
di "ace_ipw (Risk difference)"
nlcom (_b[SEN1]+_b[SEN2]+_b[SEN3])
di " rate ratio"
nlcom (_b[_cons]+_b[SEN1]+_b[SEN2]+_b[SEN3])/(_b[_cons])



exit













*Using G-comp****************************************

use "data\Simulation_timevar_SEN_v2.dta",clear
des
br

xtset id t
gen SEN_lag=L.SEN
gen Hosp_lag=L.Hosp
gen Y_lag=L.Y

gen SEN_lag2=L2.SEN
gen Hosp_lag2=L2.Hosp
gen cumden3=den*2


preserve
expand 10

sort id
by id:gen original=_n==1


**************************************************************
*Potential Outcomes at time 2
**************************************************************
*unexposed
poisson Y  male white idacibin eyfsp Hosp_lag SEN_lag if SEN_lag==0 & t==2 & original 
predict Y2_0 if t==2
*exposed
poisson Y  male white idacibin eyfsp Hosp_lag SEN_lag if SEN_lag==1 & t==2 & original
predict Y2_1 if t==2

sort id t
qui by id: replace Y2_0=Y2_0[2]
qui by id: replace Y2_1=Y2_1[2]

di "The observed and potential mean number of events at time 2 are: "
tabstat Y Y2_0 Y2_1,by(t) s(count mean sd min max) c(s)


**************************************************************
*potential confounder at time 2
**************************************************************
logit Hosp male white idacibin eyfsp Hosp_lag SEN_lag Y if t==2  & original

cap drop p2_*
*unexposed
gen p2_00= _b[_cons] + _b[male] * male+ _b[white] * white + _b[idacibin] * (idacibin==1) + _b[eyfsp]*eyfsp+ _b[Hosp_lag] * Hosp_lag +_b[SEN_lag] * 0  + _b[Y] * Y2_0
gen H2_00= runiform() < 1/(1+ exp(-p2_00))

*exposed
gen p2_11= _b[_cons] + _b[male] * male+ _b[white] * white + _b[idacibin] * (idacibin==1) + _b[eyfsp]*eyfsp+ _b[Hosp_lag] * Hosp_lag +_b[SEN_lag] * 1  + _b[Y] * Y2_1
gen H2_11= runiform() < 1/(1+ exp(-p2_11))

di "The observed and potential hospitalizaiton time 2 are: "
su Hosp H2_00 H2_11 if t==2 
drop p2_*


**************************************************************
*potential outcomes at time 3
**************************************************************
poisson Y  male white idacibin eyfsp Hosp_lag Hosp_lag2 SEN_lag SEN_lag2 Y_lag if  t==3  & original

cap drop m3_*
*unexposed
gen m3_000= exp(_b[_cons] + _b[male] * male + _b[white] * white + _b[idacibin] * (idacibin==1) + _b[eyfsp]*eyfsp +  _b[Hosp_lag] * H2_00  + _b[Hosp_lag2] * Hosp_lag2 +_b[SEN_lag] * 0 +_b[SEN_lag2] * 0 + _b[Y_lag] * Y2_0 ) 

*exposed
gen m3_111= exp(_b[_cons] + _b[male] * male + _b[white] * white + _b[idacibin] * (idacibin==1) + _b[eyfsp]*eyfsp +  _b[Hosp_lag] * H2_11  + _b[Hosp_lag2] * Hosp_lag2 +_b[SEN_lag] * 1 +_b[SEN_lag2] * 1 + _b[Y_lag] * Y2_1 ) 

su m3_* if t==3


gen Y3_000=rpoisson(m3_000)
gen Y3_111=rpoisson(m3_111)
di "The observed and potential mean number of events at time 3 are: "
su Y Y3_000 Y3_111 if  t==3

*Estimands
*difference in mean rates: 
cap drop rate*
gen rate_111=Y3_111/cumden3
gen rate_000=Y3_000/cumden3
gen rate_dif_3=Y3_111/cumden3-Y3_000/cumden3
gen rate_ratio_3=Y3_111/Y3_000
su rate* if t==3

restore

**************************************************************
	
