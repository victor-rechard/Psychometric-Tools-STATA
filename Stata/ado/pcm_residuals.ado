*Version 1.0*
**************************************************************************************************************
*Ce module calcule les résidus d'un PCM pour des items et un trait latent 																																
*INPUT																																																													
*   			item : chaine de caractère contenant le nom des items séparé par des espaces 
*				seuils: matrice qui contient les paramètres de seuil des items 
*              prediction: prediction individuelle pour une variable latent								
*OUTPUT
*				NONE 
* 				génère des variables contenant dans chacune d'elle les résidus standardisé  pour les items d'entrés et un trait latent estimé																																																        
**************************************************************************************************************

program define pcm_residuals, rclass
version 18
args  item seuils prediction 
tokenize `item'
local  nbitem : word count `item'  

forvalues p=1/`nbitem' {
quietly{
gen residus_``p''=.
local modamax= colsof(`seuils')
local N= _N 
 forvalues k=1/`N'{
    local denomin=1
	local var=0
    forvalues j = 1/`modamax' {
        local somme=0
        forvalues i=1/`j' {
            local somme=`somme'+`seuils'[`p',`i']
        }
        local denomin=`denomin'+exp(`j'*(`prediction'[`k'])-`somme')
    }
    local esperance=0
    forvalues j = 1/`modamax' {
        local somme=0
        forvalues i=1/`j' {
            local somme=`somme'+`seuils'[`p',`i']
        }
        local esperance=`esperance'+`j'*exp(`j'*(`prediction'[`k'])-`somme')/`denomin'
    }
	forvalues j = 0/`modamax' {
        local somme=0
        forvalues i=1/`j' {
            local somme=`somme'+`seuils'[`p',`i']
        }
        local var=`var'+((`j'-`esperance')^2)*exp(`j'*(`prediction'[`k'])-`somme')/`denomin'
    }
	di `esperance'
    replace residus_``p'' = (``p''[`k'] - `esperance')/sqrt(`var') if _n == `k'
}
}
di `esperance'
di "________________________________"
di "Résidus de l'Item ``p''   " 
*sum residus_``p'' 
di ""
di ""
}


display "Les résidus ont été calculés avec  succès !"

			 
	end