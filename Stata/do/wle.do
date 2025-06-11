program drop _all


program define wle ,rclass
args  thres resp iterate
tempvar sumscore 
tokenize `resp' 
egen `sumscore' =rowtotal(`resp')
qui gen wle=.

local N=_N

local maxi= colsof(`thres')
*Récupérer le nombre d'item

local J= rowsof(`thres')

forvalues n=1/`N'{
*Récupérer le nombre de modalités
local centralterm=0
forvalues j=1/`J'{
	local sumthres=``j''[`n']+1
	local centralterm=`centralterm'+DM[`j',`sumthres']
}

local sumscore_n = `sumscore'[`n']
local theta=0
wle_indiv `theta' `sumscore_n' `centralterm' `thres'
forvalues i=1/`iterate'{
	local der1=r(der1)

	local der2=r(der2)
	local theta=`theta'-(`der1'/`der2')
	wle_indiv `theta' `sumscore_n' `centralterm' `thres'
}
qui replace wle=`theta' if _n==`n'
}



end











******************************************************************
*Valeur individuel de lnv

program define wle_indiv,rclass
args theta sumscore_n centralterm thres

tempname testinfo


local J= rowsof(`thres')

local maxi= colsof(`thres')
*Récupérer le nombre d'item


* DM : design matrix "sum of each threshold for each item"
mat DM=J(`J',1,0)

forvalues i=1/`maxi'{
mat DM=DM,DM[1...,`i']+`thres'[1...,`i']
}



local lastterm=0
forvalues j=1/`J'{
	local denomin=0
forvalues l= 0/`maxi'{
	local l1=`l'+1
	local denomin =`denomin'+exp(`l'*`theta'-DM[ `j',`l1'])
}
local lastterm=`lastterm'+ln(`denomin')
}

testinfo `theta' `thres'
local testinfo=r(testinfo)

local lnv = 0.5*ln(`testinfo')+`theta'*`sumscore_n'-`centralterm'-`lastterm'
 
local h=0.0001
local theta_h1= `theta'+`h'
local theta_h2=`theta'+2*`h'

testinfo `theta_h1' `thres'
local testinfo_h1=r(testinfo)

testinfo `theta_h2' `thres'
local testinfo_h2=r(testinfo)

local lastterm_h1=0
forvalues j=1/`J'{
	local denomin=0
forvalues l= 0/`maxi'{
	local l1=`l'+1
	local denomin =`denomin'+exp(`l'*`theta_h1'-DM[ `j',`l1'])
}
local lastterm_h1=`lastterm_h1'+ln(`denomin')
}

local lastterm_h2=0
forvalues j=1/`J'{
	local denomin=0
forvalues l= 0/`maxi'{
	local l1=`l'+1
	local denomin =`denomin'+exp(`l'*`theta_h2'-DM[ `j',`l1'])
}
local lastterm_h2=`lastterm_h2'+ln(`denomin')
}



local derivative=  ((0.5*ln(`testinfo_h1')+`theta_h1'*`sumscore_n'-`centralterm'-`lastterm_h1')-(`lnv'))/`h'

local secDerivative=(`lnv'-2*(0.5*ln(`testinfo_h1')+`theta_h1'*`sumscore_n'-`centralterm'-`lastterm_h1')+(0.5*ln(`testinfo_h2')+`theta_h2'*`sumscore_n'-`centralterm'-`lastterm_h2'))/`h'^2
 
return scalar lnv = `lnv'
return scalar der1=`derivative'
return scalar der2=`secDerivative'
end 







******************************************************************
****

program define testinfo,rclass	
args theta thres
local modamax= colsof(`thres')+1
local maxi= colsof(`thres')
local J= rowsof(`thres')

local testinfo=0

forvalues j=1/`J'{
	


local part1=0
forvalues x=0/`maxi'{
	pcmformula `theta' `thres' `x' `j'
	local proba= r(proba)
	local part1=`part1'+`proba'*`x'^2
}

local part2=0
forvalues x=0/`maxi'{
	pcmformula `theta' `thres' `x' `j'
	local proba= r(proba)
	local part2=`part2'+`proba'*`x'
}

local testinfo=`testinfo'+(`part1'-`part2'^2)
}

return scalar testinfo = `testinfo'

end

****************************************************************************************************************

program define pcmformula,rclass
args theta thres x j

local modamax= colsof(`thres')+1
local maxi= colsof(`thres')
local J= rowsof(`thres')
* DM : design matrix "sum of each threshold for each item"
mat DM=J(`J',1,0)

forvalues i=1/`maxi'{
mat DM=DM,DM[1...,`i']+`thres'[1...,`i']
}
	local denomin=0
	forvalues l= 0/`maxi'{
		local l1=`l'+1
		local denomin =`denomin'+exp(`l'*`theta'-DM[ `j',`l1'])
	}
	local x1=`x'+1
	 local  proba=exp(`x'*`theta'-DM[ `j',`x1'])/`denomin'
	return scalar proba= `proba'

end