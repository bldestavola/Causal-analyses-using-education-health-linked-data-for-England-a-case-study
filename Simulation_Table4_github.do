****************************
* time-varying exposure    *
*                          *
* Estimation of ATE and ATT*
* using g-computation      *
*                          *
* as in Table 4            *
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
cd 		"C:\Users\sejjbld\OneDrive - University College London\_GRANTS\Ruth Gilbert\HOPE\WP3&4\papers\IJPDS Notes\data_analyses"

/* SAVE RESULTS */
cap log close
log using log\Simulation_Table4_github_ipw.log, replace
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

gen inter=(eyfsp_bin==0)*SEN_lag1
gen inter_ps=(idacibin ==1)*(eyfsp_bin==0)

*needed for the bootstrap
gen newid=id

*create cumulative outcome and cumulative denominator
gen cumY=0
sort id t
qui by id: replace cumY=cumY[_n-1]+Y if t>1
qui by id: gen cumden=sum(den)
qui by id: gen cumden4=cumden[3]


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


***********************************************************************************
*g-computation
***********************************************************************************
*simple specification (incorrect)
*bootstrap rate_000=r(rate_000)  rate_111=r(rate_111) ATE_RD=r(ate_RD) ATE_RR=r(ate_RR)  , reps(1000) nodrop seed(1012): Sustained_SEND_gcomp_simple
*full specification (correct)
*bootstrap rate_000=r(rate_000)  rate_111=r(rate_111) ATE_RD=r(ate_RD) ATE_RR=r(ate_RR)  , reps(1000) nodrop seed(1012): Sustained_SEND_gcomp



********************************************************************************
*IPW 
********************************************************************************

* confounders 
	vl create base_conf=(male white Hosp_base idacibin eyfsp_bin  inter_ps)
	vl create tv_conf=(Hosp Hosp_lag1)
	

drop Hosp0_ Hosp1_ Hosp3_00 Hosp3_11 
cap drop cumden 

reshape wide Y inter  Hosp SEN den SEN_lag1 Hosp_lag1 Y_lag1 SEN_lag2 Hosp_lag2 SEN_lag3 Hosp_lag3 cumden cumY,i(id) j(t)

Sustained_IPW_stab  SEN1 SEN2 SEN3 SEN4 Y1 Y2 Y3 Y4 den1 den2 den3 den4
ex
bootstrap rate_000=r(rate_000) rate_111=r(rate_111) ATE_RR=r(ATE_RR) ATE_RD=r(ATE_RD), reps(10) nodrop seed(3009) cluster(id): Sustained_IPW_stab SEN Y den 

bootstrap rate_000=r(rate_000) rate_111=r(rate_111) ATE_RR=r(ATE_RR) ATE_RD=r(ATE_RD), reps(1000) nodrop seed(3009) cluster(id) idcluster(newid): Sustained_IPW_stab SEN Y den 


ex




