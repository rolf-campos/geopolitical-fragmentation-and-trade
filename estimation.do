********************************************************************************
*                                                                              *
*            R. Campos, J. Estefania-Flores, D. Furceri, J. Timini             *
*                    "Geopolitical fragmentation and trade"                    *
*                                                                              *
*                       Journal of Comparative Economics                       *
*              Volume 51, Issue 4, December 2023, Pages 1289-1315              *
*                                                                              *
*                              REPLICATION CODE                                *
********************************************************************************

cd "<insert path here>"
clear

use "geofragtrade.dta"

* Generate fixed effects
egen exp = group(iso_o)
egen imp = group(iso_d)
egen exp_year = group(iso_o year)
egen imp_year = group(iso_d year)
egen dir_pair = group(iso_o iso_d)

* Table 3. MATR trade effects – Main estimates
ppmlhdfe trade matr, abs(exp_year imp_year dir_pair) cluster(exp imp year)
ppmlhdfe trade matr ta, abs(exp_year imp_year dir_pair) cluster(exp imp year)
ppmlhdfe trade matr ta gattwto, abs(exp_year imp_year dir_pair) cluster(exp imp year)

* Table 4. MATR trade effects – Robustness tests.
ppmlhdfe trade matr ta gatt wto, abs(exp_year imp_year dir_pair) cluster(exp imp year)
ppmlhdfe trade matr ta lntar, abs(exp_year imp_year dir_pair) cluster(exp imp year)
ppmlhdfe trade matr ta_noeu gattwto eu, abs(exp_year imp_year dir_pair) cluster(exp imp year)
ppmlhdfe trade matr_nt ta gattwto, abs(exp_year imp_year dir_pair) cluster(exp imp year)

* Generate time-varying dummy INTL_BRDR
gen INTL_BRDR = 0
replace INTL_BRDR = 1 if iso_o != iso_d

forval year = 1949/2019 {
	gen INTL_BRDR_`year' = INTL_BRDR * (year == `year')
	label var INTL_BRDR_`year' "INTL_BRDR $\times$ `year'"
}

drop INTL_BRDR_1949

* Generate time-varying distance
forval year = 1949/2019 {
	gen lndist_`year' = lndist * (year == `year')
	label var lndist_`year' "lndist $\times$ `year'"
}

drop lndist_1949

ppmlhdfe trade matr ta gattwto INTL_BRDR_*, abs(exp_year imp_year dir_pair) cluster(exp imp year)
ppmlhdfe trade matr ta gattwto lndist_*, abs(exp_year imp_year dir_pair) cluster(exp imp year)
ppmlhdfe trade matr ta gattwto INTL_BRDR_* lndist_*, abs(exp_year imp_year dir_pair) cluster(exp imp year)

