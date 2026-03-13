cap program drop SEND_gcomp
program SEND_gcomp,rclass
args Y A D H

#delimit;
poisson 	`Y' i.(male white idacibin eyfsp_bin)##i.`H' 
				i.(male white idacibin)##i.eyfsp_bin
				i.(male white)##i.idacibin
				i.male##i.white
				if `A'==0 , exp(`D') nolog
				;
#delimit cr	
predict 		y_0

#delimit;
poisson 	`Y' i.(male white idacibin eyfsp_bin)##i.`H' 
				i.(male white idacibin)##i.eyfsp_bin
				i.(male white)##i.idacibin
				i.male##i.white
				if `A'==1 , exp(`D') nolog
				;
#delimit cr	
predict 		y_1

qui su 			`D'
local 			meanfup=r(mean)


*for the ATE
su 				y_0
local 			Y0=r(mean)
su 				y_1
local 			Y1=r(mean)
local 			ATE_rate_ratio=`Y1'/`Y0'
local 			ATE_rate_dif=`Y1'/`meanfup'-`Y0'/`meanfup'

*for the ATT
su 				y_0 if `A'==1
local 			Y0_ATT=r(mean)
su 				y_1 if `A'==1
local 			Y1_ATT=r(mean)
local 			ATT_rate_ratio=`Y1_ATT'/`Y0_ATT'
qui su 			`D' if `A'==1
local 			meanfup=r(mean)
local 			ATT_rate_dif=`Y1_ATT'/`meanfup'-`Y0_ATT'/`meanfup'

 scalar Y0=`Y0'
 scalar Y1=`Y1'
 scalar Y0_ATT=`Y0_ATT'
 scalar Y1_ATT=`Y1_ATT'
 
/* scalar ATE_rate_ratio=`ATE_rate_ratio'
 scalar ATT_rate_ratio=`ATT_rate_ratio'
 scalar ATE_rate_dif=`ATE_rate_dif'
 scalar ATT_rate_dif=`ATT_rate_dif'
scalar list
*/

*post results
return scalar ATE_rate_ratio=`ATE_rate_ratio'
return scalar ATT_rate_ratio=`ATT_rate_ratio'
return scalar ATE_rate_dif=`ATE_rate_dif'
return scalar ATT_rate_dif=`ATT_rate_dif'
cap drop y_0 
cap drop y_1
end

