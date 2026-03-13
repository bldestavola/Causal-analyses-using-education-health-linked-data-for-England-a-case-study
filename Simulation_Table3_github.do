****************************
* time-fixed exposure      *
*        short-term effect *
* Estimation of ATE and ATT*
* using alternative methods*
*                          *
* as in Table 3           *
*                          *
* BDS  15/01/26            *
****************************

*************************************************************************************
/* SELECT THE FOLDER WHERE DATA and Do files ARE HELD */
cd 		"C:\Users\sejjbld\OneDrive - University College London\_GRANTS\Ruth Gilbert\HOPE\WP3&4\papers\IJPDS Notes\data_analyses"

/* SAVE RESULTS */
cap log close
log using log\Simulation_Table3_github.log, replace
*************************************************************************************
/* READ and DESCRIBE THE DATA */
use 		"data\Simulation_timevar_SEND_Study",clear
describe


*************************************************************************************
/* LOOP FROM EXPOSURE IN YEAR 1 TO OUTCOME IN YEAR 2, THEN MOVE ON TO EXPOSURE IN YEAR 2 TO OUTCOME IN YEAR 3, ETC*/

	
foreach j of numlist 1/3{                                /* j is the index for exposure time    */
	local i=`j'+1                                        /* i is the index for the outcome time */
 	di "__________________________________________________________"
	di "Short term effect of SEN at time `j' on Y at time `i'"
	di
	
	preserve
	
	
	keep 		if t==`j' | t==`i'
	sort 		id t
	qui by id: 	gen A=SEN[1]
	qui by id: 	gen H=Hosp[1]
	qui by id: 	gen Y_2=Y[2]
	qui by id: 	gen D_2=den[2]
*    qui 		gen inter=A*(eyfsp_bin==0)
	qui  		gen 		inter=(idacibin ==1)*(eyfsp_bin==0)

	gen 		rate=Y_2/D_2
	
	reshape wide

	qui su D_2 
	scalar mean=r(mean)
	qui su D_2 if A==1
	scalar mean_att=r(mean)




di in red"________________________________________"
di in red "Gcomp for SEND at time `j' on Y at time `i'___________________________"
di in blue "ATE"
*ATE
qui poisson Y_2 i.(male white idacibin eyfsp_bin H)##i.A i.(male white idacibin i.eyfsp_bin)##i.H##i.A i.(male white idacibin)##i.eyfsp_bin##i.A  i.(male white)##i.idacibin##i.A , exp(D_2) nolog robust
margins A, post
*ATE_ratio, count and rate
nlcom _b[1.A]/ _b[0bn.A]
*ATE rate, dif
nlcom (_b[1.A]/mean)- (_b[0bn.A]/mean)


di
di in blue "ATT"
*ATT
qui poisson Y_2 i.(male white idacibin eyfsp_bin H)##i.A i.(male white idacibin i.eyfsp_bin)##i.H##i.A i.(male white idacibin)##i.eyfsp_bin##i.A  i.(male white)##i.idacibin##i.A, exp(D_2) nolog robust
margins A , post vce(unconditional)  subpop(A)
*ATT_ratio, count and rate
nlcom _b[1.A]/ _b[0bn.A]
*ATT dif, rate
nlcom (_b[1.A]/mean_att)- (_b[0bn.A]/mean_att)

di "IPW for SEND at time `j' on Y at time `i'___________________________________________________"

*IPW
vl create 	base_conf=(male white idacibin eyfsp_bin inter)

bootstrap  ATE_rr=r(ATE_rr)  ATE_rd=r(ATE_rd)  ATT_rr=r(ATT_rr)  ATT_rd=r(ATT_rd)  , reps(1000) nodrop seed(1212) : SEND_IPW Y_2 A D_2 H


di "AIPW for SEND at time `j' on Y at time `i'___________________________________________________"
*AIPW

*model 4, correct PS, and correct Y
vl modify 	base_conf=(male white idacibin eyfsp_bin H inter)
*ATE
*ATE
teffects aipw (rate i.(male white idacibin eyfsp_bin H) i.(male white idacibin i.eyfsp_bin)##i.H i.(male white idacibin)##i.eyfsp_bin  i.(male white)##i.idacibin  i.male##i.white,poisson)   (A  $base_conf)
*ATE ratio
nlcom ( _b[POmean:0.A]+_b[ATE:r1vs0.A])/_b[POmean:0.A]

*ATT
teffects aipw (rate i.(male white idacibin eyfsp_bin H) i.(male white idacibin i.eyfsp_bin)##i.H i.(male white idacibin)##i.eyfsp_bin  i.(male white)##i.idacibin  i.male##i.white , poisson)   (A  $base_conf), atet
*ATT ratio
nlcom ( _b[POmean:0.A]+_b[ATET:r1vs0.A])/_b[POmean:0.A]


di "___________________________________________________"
di "___________________________________________________"


restore

}


di 
di
di "______________________________________"
di "assuming constant short-term effect"



reshape long
 
*THE TIME-VARYING INTERVENTION 
	rename SEN A

*THE  TIME-VARYING CONFOUNDER 
	rename Hosp H
	

sort id t
qui by id: 	gen A_lag=A[_n-1]
qui by id: 	gen H_lag=H[_n-1]
gen 		D=den

******g-comp WITH MARGINS
*ATE
qui poisson 	Y i.A_lag##i.(male white idacibin eyfsp_bin H_lag t) if t>1, exp(D) nolog irr cluster(id)
margins 	A_lag, post
*ATE_ratio, count and rate
nlcom _b[1.A_lag]/ _b[0bn.A_lag]
*ATE_dif
nlcom _b[1.A_lag]/mean- _b[0bn.A_lag]/mean

*ATT
qui poisson Y i.A_lag##i.(male white idacibin eyfsp_bin H_lag t) if t>1, exp(D) nolog irr cluster(id) robust
margins A_lag , post vce(unconditional)  subpop(A)
*ATT_ratio, count and rate
nlcom _b[1.A_lag]/ _b[0bn.A_lag]

*ATT_dif
nlcom _b[1.A_lag]/mean_att- _b[0bn.A_lag]/mean_att


ex

