/* =================================================================== */

// Prepare Inputs (table_0, table_t, table_j, table_tj, table_i, table_ti, table_ji, table_tji)

// Check consistency 

// Copyright (C) 2015 by Yun Ling

/* =================================================================== */

clear all
set mem 10g

gl tb "D:\Dropbox\Dissertation\Data\CrunchBase\Raw Table"
gl tbc "D:\Dropbox\Dissertation\Data\CrunchBase\Table"
gl ext "D:\Dropbox\Dissertation\Data\CrunchBase\external"
gl py "D:\Dropbox\Dissertation\Data\CrunchBase\py"
gl s1 "D:\Dropbox\Dissertation\Data\CrunchBase\Inputs\1 screening results"
gl s2 "D:\Dropbox\Dissertation\Data\CrunchBase\Inputs\2 input tables"
gl cpp "D:\Dropbox\Dissertation\Data\CrunchBase\Inputs\2 input tables\cpp"


// Time Horizon
gl begy = 1998
gl endy = 2014

// Frequency is monthly data
set more off
set seed 1

/* =============================== table_0 ============================ */

use "$s1\\s1(b)", clear
codebook s
codebook k
codebook ym
duplicates r k ym
sort k ym
// t
g t = ym - ym(${begy},1)
sum t
// j
egen j = group(k)
replace j = j-1
sum j
preserve
keep j k 
duplicates drop
save "$s2\\j", replace
restore
// match
g match = s
sum match
// status_code
// r_obs
g r_obs = ret
order t j status_code match r_obs
keep t j status_code match r_obs
xtset j t
tsfill
replace match = l.match if match==. & l.match<.
replace match = -1 if match==.
replace r_obs = -999 if r_obs==.
replace status_code = 0 if status_code==.
sort t j
codebook j
codebook match
save "$s2\\table_0", replace

/* =============================== table_t ============================ */

/* table_t: 
   t, mktrf, smb, hml, rf */
   
use "$ext\\ff_us", clear
foreach v in mktrf smb hml rf {
	ren `v' `v'_us
}
merge m:m t using "$ext\\ff_global"
foreach v in mktrf smb hml wml rf {
	ren `v' `v'_gl
}
drop _merge
replace t = t - ym(${begy},1)
format t %8.0g
keep if t>=0 & t<= ym(${endy},12)-ym(${begy},1)
drop rf_gl
ren rf_us rf
order t mktrf* smb* hml* wml* rf
drop wml_gl 
save "$s2\\table_t", replace

/* =============================== table_j and table_tj ============================ */

/* table_j: 
   j, 
   location, 
   category, 
   product 
*/

set more off

// (1) locations
use "$s2\\j", clear
joinby k using "$tb\\organization_offices"
g l = city_l
replace l = region_l if l==.
replace l = country_l if l==.
drop if l==.
keep k l
forvalues i = 4(-1)1 {
	merge m:m l using "$tb\\locations_`i'"
	drop if _merge==2
	g l`i' = l if _merge==3
	g name`i' = name if _merge==3
	replace l = parent_l if _merge==3
	keep k l l* l name*
	drop name
}
drop l
duplicates drop
g l = l4
replace l = l3 if l==.
replace l = l2 if l==.
replace l = l1 if l==.
g d_ca = (name3=="California")
g d_ny = (name3=="New York")
g d_tx = (name3=="Texas")
g d_ous = (name2=="United States") * (d_ca==0 & d_ny==0 & d_tx==0)
g d_ona = (name1=="North America") * (d_ca==0 & d_ny==0 & d_tx==0 & d_ous==0)
g d_as = (name1=="Asia")
g d_eu = (name1=="Europe")
egen n_loc = count(l), by(k)
preserve
egen N = count(k), by(name4)
keep name4 N
duplicates drop
drop if name4==""
gsort -N
g n = _n
drop if _n>20
tab n, generate(loc)
ren n loc
save "$s2\\loc_dummies_j", replace
restore
merge m:m name4 using "$s2\\loc_dummies_j"
keep k n_loc loc* d_*
drop loc
forvalues i = 1/20 {
	replace loc`i' = 0 if loc`i'==.
	egen temp = max(loc`i'), by(k)
	drop loc`i'
	ren temp loc`i'
}
foreach v of varlist d_* {
	egen temp = max(`v'), by(k)
	drop `v'
	ren temp `v'
}
drop location_type
duplicates drop
egen n_loc_major = rowtotal(loc*)
order k n* d* loc*
replace n_loc = 1 if n_loc==.
merge m:m k using "$s2\\j"
replace n_loc = 1 if n_loc==.
replace n_loc_major = 0 if n_loc_major==.
foreach v of varlist d_* loc* {
	replace `v' = 0 if `v'==.
}
sort j 
drop k
order j
drop _merge
count
codebook j
tempfile 1
save "`1'", replace

// (2) categories
use "$s2\\j", clear
joinby k using "$tb\\organization_categories"
merge m:m c using "$tb\\categories"
drop if _merge==2
drop _merge
keep k c name
egen n_cat = count(c), by(k)
preserve
egen N = count(k), by(name)
keep name N
duplicates drop
drop if name==""
gsort -N
g n = _n
drop if _n>20
tab n, generate(cat)
ren n cat
save "$s2\\cat_dummies_j", replace
restore
merge m:m name using "$s2\\cat_dummies_j"
keep k n_cat cat*
drop cat
forvalues i = 1/20 {
	replace cat`i' = 0 if cat`i'==.
	egen temp = max(cat`i'), by(k)
	drop cat`i'
	ren temp cat`i'
}
duplicates drop
merge m:m k using "$s2\\j"
replace n_cat = 1 if n_cat==. | n_cat==0
foreach v of varlist cat* {
	replace `v' = 0 if `v' ==.
}
drop _merge
sort j
order j
drop k
codebook j
tempfile 2
save "`2'", replace

// (3) product
use "$s2\\j", clear
joinby k using "$tb\\organization_products"
egen n_prod = count(d), by(k)
drop d
duplicates drop
merge m:m k using "$s2\\j"
replace n_prod = 0 if n_prod==.
drop _merge
sort j
order j
drop k
codebook j
tempfile 3
save "`3'", replace

// (4) combine
use "`1'", clear
merge m:m j using "`2'"
drop _merge
merge m:m j using "`3'"
drop _merge
order j n*
count
codebook j
foreach v of varlist _all {
	ren `v' `v'_st
}
ren j_st j
save "$s2\\table_j", replace


/* table_tj:
   t, j, 
   age, 
   funding factors, 
   human factors 
*/

// (A-1) people degree
use "$tb\\person_degrees", clear
ren organization_name school
ren degree_subject subject
ren degree_type_name type
keep p school subject type 
replace subject = lower(subject)
g d_engr = max(regexm(subject,"engineer"),regexm(subject,"computer"),regexm(subject,"cs"), ///
	regexm(subject,"ee"),regexm(subject,"electric"),regexm(subject,"mechani"), ///
	regexm(subject,"indus"))
g d_bus = max(regexm(subject,"econ"),regexm(subject,"fin"),regexm(subject,"bus"), ///
	regexm(subject,"acc"),regexm(subject,"market"),regexm(subject,"manage"),regexm(subject,"mgmt"), ///
	regexm(subject,"strategy"),regexm(subject,"entrep"))
g d_lp = max(regexm(subject,"law"),regexm(subject,"politic"),regexm(type,"jd"))
drop subject
replace type = subinstr(lower(type),".","",.)
g d_phd = max(regexm(type,"phd"),regexm(type,"doctor"))
g d_ms = max(regexm(type,"ma"),regexm(type,"ms"),regexm(type,"master"))
g d_mba = regexm(type,"mba")
drop type
replace school = subinstr(lower(school),",","",.)
g d_top20 = max(regexm(school,"harvard"),regexm(school,"stanford"),regexm(school,"mit"), ///
	regexm(school,"massachusetts institute of technology"),regexm(school,"sloan"), /// 
	regexm(school,"ucb"),regexm(school,"hass"),regexm(school,"berkeley"), ///
	regexm(school,"oxford"),regexm(school,"cambridge"), ///
	regexm(school,"princeton"),regexm(school,"university of chicago"), ///
	regexm(school,"columbia"),regexm(school,"duke"),regexm(school,"zurich"), ///
	regexm(school,"university of pennsylvania"),regexm(school,"wharton"), ///
	regexm(school,"caltech"),regexm(school,"california institute of technology"), ///
	regexm(school,"ucla"),regexm(school,"california")*regexm(school,"los angeles"), ///
	regexm(school,"yale"),regexm(school,"johns hopkins"),regexm(school,"cornell"), ///
	regexm(school,"northwestern"),regexm(school,"kellogg"),regexm(school,"tel aviv"), ///
	regexm(school,"new york university"),regexm(school,"nyu"),regexm(school,"stern"), ///
	regexm(school,"cmu"),regexm(school,"carnegie mellon university"))
g d_stanford = regexm(school,"stanford")
g d_harvard = regexm(school,"harvard")
g d_mit = max(regexm(school,"mit"),regexm(school,"massachusetts institute of technology"),regexm(school,"sloan"))
g d_ucb = max(regexm(school,"ucb"),regexm(school,"berkeley"),regexm(school,"hass"))
g d_oxford = regexm(school,"oxford")
g d_cambridge = regexm(school,"cambridge")
g d_princeton = regexm(school,"princeton")
g d_chicago = regexm(school,"university of chicago")
g d_columbia = regexm(school,"columbia")
g d_duke = regexm(school,"duke")
g d_zurich = regexm(school,"zurich")
g d_upenn = max(regexm(school,"university of pennsylvania"),regexm(school,"wharton"))
g d_caltech = max(regexm(school,"ucla"),regexm(school,"california institute of technology"))
g d_ucla = max(regexm(school,"university of california los angeles"),regexm(school,"ucla"))
g d_yale = regexm(school,"yale")
g d_jh = regexm(school,"johns hopkins")
g d_cornell = regexm(school,"cornell university")
g d_northwestern = max(regexm(school,"northwestern"),regexm(school,"kellogg"))
g d_telaviv = regexm(school,"tel aviv")
g d_nyu = max(regexm(school,"new york university"),regexm(school,"nyu"),regexm(school,"stern"))
g d_cmu = max(regexm(school,"cmu"),regexm(school,"carnegie mellon university"))
drop school
foreach v of varlist d_* {
	egen temp = max(`v'), by(p)
	drop `v'
	ren temp `v'
}
duplicates drop	
save "$s2\\degrees", replace

// (A-2) people employment
use "$tb\\organization_current_team", clear
g started_t = ym(year(started_on), month(started_on)) - ym(${begy},1)
g current = 1
tempfile 2
save "`2'", replace
use "$tb\\organization_past_team", clear
g started_t = ym(year(started_on), month(started_on)) - ym(${begy},1)
g ended_t = ym(year(ended_on), month(ended_on)) - ym(${begy},1)
g current = 0
append using "`2'"
keep k p *_t current
sort p started_t
duplicates drop
preserve
use "$s2\\table_0", clear
egen birth_t = min(t), by(j)
egen death_t = max(t), by(j)
keep j birth_t death_t 
duplicates drop
merge m:m j using "$s2\\j"
drop _merge
tempfile 3
save "`3'", replace
restore
joinby k using "`3'"
sort p started_t
replace ended_t = death_t if ended_t==. & current==1
replace ended_t = death_t if ended_t==. | ended_t>death_t
replace started_t = birth_t if started_t==. | started_t<birth_t
drop if started_t >= ended_t
keep k p started_t ended_t j 
duplicates drop
sort j p 
order j p
codebook j
save "$s2\\j-people", replace

// (B) start
set more off

// (1) age and funding
use "$s2\\table_0", clear
egen birth_t = min(t), by(j)
g age = (t - birth_t)/12
drop birth_t
g fund = 1 if status==0 & match>-1 
sort j t
bysort j: g n_rounds = sum(fund)
replace fund = 1 if status==-2
sort j t
bysort j n_rounds: g t_fromlast = _n
xtset j t 
g temp = l.t_fromlast
replace temp = 0 if temp==.
drop t_fromlast
ren temp t_fromlast
replace t_fromlast = t_fromlast/12 // in years
g t2_fromlast = t_fromlast^2
keep t j age n_rounds t_fromlast t2_fromlast
sort t j
count
tempfile 1
save "`1'", replace

// (2) people founded companies
use "$tb\\person_founded_companies", clear
merge m:m k using "$tb\\organization"
drop if _merge==2
keep p k founded_on
g founded_t = ym(year(founded_on),month(founded_on)) - ym(${begy},1)
drop if founded_on==.
drop founded_on
sort p founded_t
joinby p using "$s2\\j-people"
keep if founded_t < started_t // founded previously
drop founded_t j
egen t1 = min(started_t), by(k p)
egen t2 = max(ended_t), by(k p)
drop started_t ended_t
duplicates drop // count people, not # of founded companies 
reshape long t, i(k p) j(index)
egen id = group(k p)
xtset id t
tsfill
replace k = l.k if k==.
keep k t
egen n_founded = count(t), by(k t)
keep k t n_founded
duplicates drop
format t %tm
joinby k using "$s2\\j"
merge m:m j t using "$s2\\table_0"
drop if _merge==1
keep j t n_founded
replace n_founded = 0 if n_founded==.
sort t j 
order t j
codebook j
count
tempfile 2
save "`2'", replace

// (3) people degrees
use "$s2\\degrees", clear
joinby p using "$s2\\j-people"
egen t1 = min(started_t), by(p j)
egen t2 = max(ended_t), by(p j)
drop started_t ended_t
duplicates drop
reshape long t, i(p j) j(index)
order p j t
egen id = group(p j)
xtset id t
tsfill
replace p = l.p if p==.
replace j = l.j if j==.
foreach v of varlist d* {
	replace `v' = l.`v' if `v'==.
}
drop p 
collapse (max) d*, by(j t)
tempfile 3
save "`3'", replace

// (4) combine
use "`1'", clear
merge m:m t j using "`2'"
drop if _merge==2
drop _merge
merge m:m t j using "`3'"
drop if _merge==2
drop _merge
keep t - d_top20
foreach v of varlist _all {
	ren `v' `v'_st
}
foreach v of varlist d* {
	replace `v' = 0 if `v'==.
}
ren t_st t
ren j_st j
count
save "$s2\\table_tj", replace

/* =============================== table_i and table_ti ============================ */

use "$s2\\table_0", clear
keep match
keep if match~=-1
duplicates drop
ren match s
merge m:m s using "$s1\\s1(c)"
keep if _merge==3
drop _merge deal_no s
duplicates drop
sort k
codebook k
count
di r(N)*(${endy}-${begy}+1)*12
save "$s2\\i", replace

/* table_i: 
   i, 
   location, 
*/

set more off

use "$s2\\i", clear
joinby k using "$tb\\organization_offices"
g l = city_l
replace l = region_l if l==.
replace l = country_l if l==.
duplicates drop
keep k l
forvalues i = 4(-1)1 {
	merge m:m l using "$tb\\locations_`i'"
	drop if _merge==2
	g l`i' = l if _merge==3
	g name`i' = name if _merge==3
	replace l = parent_l if _merge==3
	keep k l l* l name*
	drop name
}
drop l
duplicates drop
g l = l4
replace l = l3 if l==.
replace l = l2 if l==.
replace l = l1 if l==.
g d_ca = (name3=="California")
g d_ny = (name3=="New York")
g d_tx = (name3=="Texas")
g d_ous = (name2=="United States") * (d_ca==0 & d_ny==0 & d_tx==0)
g d_ona = (name1=="North America") * (d_ca==0 & d_ny==0 & d_tx==0 & d_ous==0)
g d_as = (name1=="Asia")
g d_eu = (name1=="Europe")
egen n_loc = count(l), by(k)
preserve
egen N = count(k), by(name4)
keep name4 N
duplicates drop
drop if name4==""
gsort -N
g n = _n
drop if _n>20
tab n, generate(loc)
ren n loc
save "$s2\\loc_dummies_i", replace
restore
merge m:m name4 using "$s2\\loc_dummies_i"
keep k n_loc loc* d_*
drop loc
forvalues i = 1/20 {
	replace loc`i' = 0 if loc`i'==.
	egen temp = max(loc`i'), by(k)
	drop loc`i'
	ren temp loc`i'
}
foreach v of varlist d_* {
	egen temp = max(`v'), by(k)
	drop `v'
	ren temp `v'
}
drop location_type
replace n_loc = 1 if n_loc==. | n_loc==0
duplicates drop
egen n_loc_major = rowtotal(loc*)
merge m:m k using "$s2\\i"
drop if _merge==1
drop _merge
replace n_loc = 1 if n_loc==.
replace n_loc_major = 0 if n_loc_major==.
foreach v of varlist d_* loc* {
	replace `v' = 0 if `v'==.
}
sort k
ren k i
sort i
count
foreach v of varlist _all {
	ren `v' `v'_vc
}
ren i_vc i
sort i
save "$s2\\table_i", replace

/* table_ti:
   t, i,
   funding factors, 
   human factors 
*/

// (A) employment
use "$tb\\organization_current_team", clear
g started_t = ym(year(started_on), month(started_on)) - ym(${begy},1)
g current = 1
tempfile 2
save "`2'", replace
use "$tb\\organization_past_team", clear
g started_t = ym(year(started_on), month(started_on)) - ym(${begy},1)
g ended_t = ym(year(ended_on), month(ended_on)) - ym(${begy},1)
g current = 0
append using "`2'"
keep k p *_t current
sort p started_t
preserve
use "$s2\\i", clear
merge m:m k using "$tb\\organization"
keep if _merge==3
keep k founded_on closed_on
g birth_t = ym(year(founded_on),month(founded_on)) - ym(${begy},1)
g death_t = ym(year(closed_on), month(closed_on)) - ym(${begy},1)
keep k birth_t death_t 
duplicates drop
tempfile 3
save "`3'", replace
restore
merge m:m k using "`3'"
keep if _merge==3
drop _merge
codebook k
sort p started_t
replace ended_t = death_t if ended_t==. & current==1
replace ended_t = death_t if ended_t==. | ended_t>death_t
replace started_t = birth_t if started_t==. | started_t<birth_t
keep k p started_t ended_t
replace started_t = max(started_t, 0)
replace ended_t = min(ended_t, ym(${endy},12) - ym(${begy},1))
keep if started_t < ended_t
duplicates drop
ren k i
sort i p
save "$s2\\i-people", replace

// (B) tables

set more off

// (1) funding rounds
use "$s2\\table_0", clear
keep if match~=-1
ren match s
merge m:m s using "$s1\\s1(c)"
drop _merge
keep t k j
order k t 
g fund = 1
collapse (count) fund, by(k t)
xtset k t
count
loc last = r(N)+2
set obs `last'
loc last = r(N)+1
replace t = 0 in `last'
loc last = `last'+1
replace t = ym(${endy},12) - ym(${begy},1) in `last'
tsfill, full
bysort k: g n_fund = sum(fund)
drop if k==.
keep k t n_fund
ren k i
duplicates drop
codebook i
count
tempfile 1
save "`1'", replace

// (b) funded categories
use "$s2\\table_0", clear
keep if match~=-1
ren match s 
merge m:m s using "$s1\\s1(c)"
drop _merge
keep t k j 
ren k i
ren j k
joinby k using "$tb\\organization_categories"
replace c = 0 if c==. // add the dummy category to missing ones
drop k 
order i c t
sort i c t
egen min_t = min(t), by(i c) // find the earliest funding time for a category
format min_t %tm
keep i c min_t
ren min_t t
duplicates drop
duplicates r i c
egen id = group(i c)
xtset id t
count
loc last = r(N)+2
set obs `last'
loc last = r(N)+1
replace t = 0 in `last'
loc last = `last'+1
replace t = ym(${endy},12) - ym(${begy},1) in `last'
tsfill, full
g fund = 1 if i~=. & c~=.
foreach v in i c {
	egen temp = mean(`v'), by(id)
	drop `v'
	ren temp `v'
}
order i c
bysort id: g cat_fund = sum(fund)
collapse (sum) cat_fund, by(i t)
drop if i==.
duplicates drop
codebook i
count
tempfile 2
save "`2'", replace

// (3) people factors
use "$s2\\degrees", clear
joinby p using "$s2\\i-people"
egen t1 = min(started_t), by(p i)
egen t2 = max(ended_t), by(p i)
drop started_t ended_t
duplicates drop
reshape long t, i(p i) j(index)
order p i t
egen id = group(p i)
xtset id t
tsfill
replace p = l.p if p==.
replace i = l.i if i==.
foreach v of varlist d* {
	replace `v' = l.`v' if `v'==.
}
drop p
collapse (max) d*, by(i t)
duplicates drop
codebook i
count
tempfile 3
save "`3'", replace

// combine
use "`1'", clear
merge m:m t i using "`2'"
drop _merge
merge m:m t i using "`3'"
drop _merge
foreach v in n_fund cat_fund {
	replace `v' = 0 if `v'==.
}
foreach v of varlist d_* {
	replace `v' = 0 if `v'==.
}
foreach v of varlist _all {
	ren `v' `v'_vc
}
ren t_vc t
ren i_vc i
order t i
sort t i
codebook i
count
save "$s2\\table_ti", replace

/* =============================== table_ji and table_tji ============================ */

/* table_ji:
   j, i,
   location factors
*/

set more off

// i-j pairs
use "$s2\\i", clear
ren k i
g x = 1
tempfile 1
save "`1'", replace
use "$s2\\j", clear
g x = 1
tempfile 2
save "`2'", replace
joinby x using "`1'"
drop x
codebook i
codebook j
count
save "$s2\\i-j", replace

// (1) i locations
use "$s2\\i", clear
joinby k using "$tb\\organization_offices"
g l = city_l
replace l = region_l if l==.
replace l = country_l if l==.
duplicates drop
keep k l
forvalues i = 4(-1)1 {
	merge m:m l using "$tb\\locations_`i'"
	drop if _merge==2
	g l`i' = l if _merge==3
	g name`i' = name if _merge==3
	replace l = parent_l if _merge==3
	keep k l l* l name*
	drop name
}
drop l
duplicates drop
g l = l4
replace l = l3 if l==.
replace l = l2 if l==.
replace l = l1 if l==.
ren k i
forvalues i = 1/4 {
	replace location_type = `i' if l`i'<.
}
foreach v of varlist _all {
	ren `v' `v'_vc
}
ren i_vc i
drop if l_vc==.
duplicates drop
foreach v in vc {
	g name_`v' = name1_`v'
	replace name_`v' = name2_`v'+", "+name_`v' if l2_`v'<.
	replace name_`v' = name3_`v'+", "+name_`v' if l3_`v'<.
	replace name_`v' = name4_`v'+", "+name_`v' if l4_`v'<.
}
save "$s2\\i-locs", replace

// (2) j locations
use "$s2\\j", clear
joinby k using "$tb\\organization_offices"
g l = city_l
replace l = region_l if l==.
replace l = country_l if l==.
keep k j l
forvalues i = 4(-1)1 {
	merge m:m l using "$tb\\locations_`i'"
	drop if _merge==2
	g l`i' = l if _merge==3
	g name`i' = name if _merge==3
	replace l = parent_l if _merge==3
	keep k l l* l name* j
	drop name
}
drop l
duplicates drop
g l = l4
replace l = l3 if l==.
replace l = l2 if l==.
replace l = l1 if l==.
forvalues i = 1/4 {
	replace location_type = `i' if l`i'<.
}
foreach v of varlist _all {
	ren `v' `v'_st
}
ren j_st j
ren k_st k
drop if l_st==.
duplicates drop
foreach v in st {
	g name_`v' = name1_`v'
	replace name_`v' = name2_`v'+", "+name_`v' if l2_`v'<.
	replace name_`v' = name3_`v'+", "+name_`v' if l3_`v'<.
	replace name_`v' = name4_`v'+", "+name_`v' if l4_`v'<.
}
save "$s2\\j-locs", replace

// (3) combine
/*
use "$s2\\i-locs", clear
keep name_vc l_vc
ren name_vc name_st
ren l_vc l_st
append using "$s2\\j-locs"
keep name_st l_st
ren name_st name
ren l_st l
duplicates drop
count
outsheet using "$s2\\py\\locations.csv", comma quotes replace
*/
use "$s2\\i-locs", clear
keep i name_vc
ren name_vc name
merge m:m name using "$s2\\py\\locations_latlong"
keep if _merge==3
drop _merge name
ren latitude lat_vc
ren longtitude long_vc
append using "$s2\\py\\i-locs-add"
sort i
codebook i
outsheet i lat_vc long_vc using "$s2\\py\\i.csv", comma nonames replace
use "$s2\\j-locs", clear
keep j name_st
ren name_st name
merge m:m name using "$s2\\py\\locations_latlong"
keep if _merge==3
drop _merge name
ren latitude lat_st
ren longtitude long_st
append using "$s2\\py\\j-locs-add"
sort j 
codebook j
outsheet j lat_st long_st using "$s2\\py\\j.csv", comma nonames replace
// calculate_distance
insheet using "$s2\\py\\distance.csv", comma names clear
merge m:m i j using "$s2\\i-j"
drop _merge k
g d = 0 if dist<=100 // 50 mph for driving
replace d = 1 if dist<=1000 & d==. // 500 mph for flight
replace d = 2 if d==. // over 1-day round trip
sort j i
order j i
save "$s2\\table_ji", replace


/* table_tji:
   t, j, i,
   funding factors,
   category factors,
   human factors
*/

// (A) table index
use "$s1\\s1(c)", clear
keep s k 
ren k i
duplicates drop
ren s match
joinby match using "$s2\\table_0"
ren match s
keep t i j
order t i j
save "$s2\\t-i-j", replace

set more off

// (1) funding
/*
use "$s2\\t-i-j", clear
collapse (min) t, by(i j)
replace t = t+1
drop if t > ym(${endy},12) - ym(${begy},1)
qui count
loc last = r(N) + 1
set obs `last'
replace t = ym(${endy},12) - ym(${begy},1) in `last'
egen id = group(i j)
xtset id t
g funded = 1
tsfill, full
drop if id==.
foreach v in i j {
	egen temp = mean(`v'), by(id)
	drop `v'
	ren temp `v'
}
sort id t
bysort id: g temp = sum(funded)
drop funded
ren temp funded
drop id
drop if funded==0
duplicates drop
save "$s2\\table_tji_sparse(1)", replace
*/

// (2) category
use "$s2\\t-i-j", clear
merge m:m j using "$s2\\j"
keep if _merge==3
drop _merge
joinby k using "$tb\\organization_categories"
drop if c==.
drop k
egen t_min = min(t)
g t_mask = t - t_min
foreach v in i j c {
	egen `v'_mask = group(`v')
	replace `v'_mask = `v'_mask - 1
}
tempfile 2
save "`2'", replace
keep i_mask j_mask t_mask
duplicates drop
sort i_mask j_mask t_mask
outsheet i_mask j_mask t_mask using "$cpp\\investments.txt", nonames replace
use "`2'", clear
keep j_mask c_mask
duplicates drop
sort j_mask c_mask
outsheet j_mask c_mask using "$cpp\\j-category.txt", nonames replace
cd "$cpp"
//shell "$cpp\\category.exe"
insheet using "$cpp\\i-j category-match.txt", clear
ren v1 t_mask
ren v2 i_mask
ren v3 j_mask
ren v4 n_same_cat
merge m:m t_mask i_mask j_mask using "`2'"
drop _merge
egen temp = mean(t_min)
replace t = t_mask + temp
drop temp t_min t_mask
foreach v in i j c {
	egen temp = mean(`v'), by(`v'_mask)
	drop `v'
	ren temp `v'
	drop `v'_mask
}
keep t i j n_same_cat
duplicates drop
merge 1:1 i j t using "$s2\\t-i-j"
drop _merge
replace n_same_cat = 1 if n_same_cat==.
sort t j i 
order t j i n_same_cat
replace t = t + 1 
drop if t > ym(${endy},12)
save "$s2\\table_tji_sparse(2)", replace

// (3) common people
use "$s2\\i-people", clear
append using "$s2\\j-people"
egen p_mask = group(p)
replace p_mask = p_mask - 1
egen i_mask = group(i) if i<.
replace i_mask = i_mask - 1
egen j_mask = group(j) if j<.
replace j_mask = j_mask - 1
egen min_t = min(started_t)
g started_mask = started_t - min_t
g ended_mask = ended_t - min_t
outsheet i_mask p_mask started_mask ended_mask if i<. using "$cpp\\i-people.txt", nonames replace
outsheet j_mask p_mask started_mask ended_mask if i==. using "$cpp\\j-people.txt", nonames replace
keep i* j* min_t
duplicates drop
preserve
keep if i<.
keep i i_mask
duplicates drop
tempfile i
save "`i'"
restore
preserve
keep if i==.
keep j j_mask
duplicates drop
tempfile j
save "`j'"
restore
sum min_t
loc min_t = r(min)
cd "$cpp"
shell "people.exe" "i-people.txt" "j-people.txt" "i-j people-match.txt"
insheet using "$cpp\\i-j people-match.txt", clear
ren v1 i_mask
ren v2 j_mask
ren v3 t_mask
ren v4 n_people
merge m:m i_mask using "`i'"
keep if _merge==3
drop _merge
merge m:m j_mask using "`j'"
keep if _merge==3
drop _merge
g min_t = `min_t'
g t = t_mask + min_t
keep i j t n_people
sort t j i 
order t j i 
save "$s2\\table_tji_sparse(3)", replace

// (4) common school
use "$s2\\degrees", clear
cap drop n_degrees
drop d_engr - d_top20
g school = ""
order d*, last
loc j = 1
loc S1 = `j'
foreach v of varlist d_stanford - d_cmu {
	ren `v' d_school`j'
	loc j = `j'+1
}
loc S2 = `j' - 1
forvalues s = `S1'/`S2' {
	replace school = school+"" if d_school`s'==0
	replace school = school+",`s'" if d_school`s'==1
}
replace school = substr(school,2,length(school))
drop d_*
split school, p(",")
drop school 
foreach v in school* {
	destring `v', replace
}
egen n_school = rownonmiss(school*)
drop if n_school==0
expand n_school
bysort p: g n = _n
tostring n, replace
foreach v of varlist school* {
	replace `v' = . if regexm("`v'",n)==0
}
egen sch = rowmin(school*)
drop school*
keep p sch
drop if sch==.
duplicates drop
tempfile 1
save "`1'", replace
use "$s2\\i-people", clear
append using "$s2\\j-people"
sort i j p started_t ended_t
joinby p using "`1'"
egen t_sch1 = min(started_t), by(sch j)
egen t_sch2 = max(ended_t), by(sch j)
sum started_t
loc min_t = r(min)
foreach v of varlist t_* {
	replace `v' = `v' - `min_t'
}
egen sch_mask = group(sch)
replace sch_mask = sch_mask-1
egen i_mask = group(i) if i<.
replace i_mask = i_mask-1
egen j_mask = group(j) if i==.
replace j_mask = j_mask-1
duplicates drop
preserve
keep if i<.
keep i_mask sch_mask t_sch1 t_sch2 
order i_mask sch_mask t_sch1 t_sch2 
duplicates drop
sort i_mask sch_mask 
outsheet using "$cpp\\i-school.txt", nonames replace
restore
preserve
keep if i==.
keep j_mask sch_mask t_sch1 t_sch2
order j_mask sch_mask t_sch1 t_sch2
duplicates drop
sort j_mask sch_mask
outsheet using "$cpp\\j-school.txt", nonames replace
restore
// i:i_mask
preserve
keep if i<.
keep i i_mask
duplicates drop
tempfile i
save "`i'", replace
restore
// j:j_mask
preserve
keep if i==.
keep j j_mask
duplicates drop
tempfile j
save "`j'", replace
restore
// school-match
shell "people.exe" "i-school.txt" "j-school.txt" "i-j school-match.txt"
insheet using "$cpp\\i-j school-match.txt", clear
ren v1 i_mask
ren v2 j_mask
ren v3 t_mask
ren v4 n_school
merge m:m i_mask using "`i'"
keep if _merge==3
drop _merge
merge m:m j_mask using "`j'"
keep if _merge==3
drop _merge
g min_t = `min_t'
g t = t_mask + min_t
keep i j t n_school
sort t j i
order t j i
save "$s2\\table_tji_sparse(4)", replace

/* =============================== aggregate table_i and table_ti ============================ */

/* table_i:
   i,
   size,
   location (table_i aggregate)
*/

set more off

// (1) size
use "$s1\\s1(c)", clear
drop deal_no
duplicates drop
egen size_vc = count(k), by(s)
keep s size_vc k
preserve
keep s size_vc
duplicates drop
tempfile 1
save "`1'", replace
restore

// (2) location
joinby k using "$tb\\organization_offices"
g l = city_l
replace l = region_l if l==.
replace l = country_l if l==.
duplicates drop
keep s k l
forvalues i = 4(-1)1 {
	merge m:m l using "$tb\\locations_`i'"
	drop if _merge==2
	g l`i' = l if _merge==3
	g name`i' = name if _merge==3
	replace l = parent_l if _merge==3
	keep s k l l* l name*
	drop name
}
drop l
duplicates drop
g l = l4
replace l = l3 if l==.
replace l = l2 if l==.
replace l = l1 if l==.
duplicates drop
egen m_loc_vc = count(l), by(s)
drop name* l*
preserve
tempfile 2
save "`2'",replace
restore
merge m:m s using "`2'"
drop _merge
ren k i
merge m:m i using "$s2\\table_i"
drop i n_loc_major _merge
forvalues i = 1/20 {
	egen temp = max(loc`i'_vc), by(s)
	drop loc`i'_vc
	ren temp loc`i'_vc
}
foreach v of varlist d_* {
	egen temp = max(`v'), by(s)
	drop `v'
	ren temp `v'
}
egen n_loc_major_vc = rowtotal(loc*)
drop if s==.
drop n_loc_vc
ren m_loc_vc n_loc_vc
merge m:m s using "`1'"
drop _merge
duplicates drop
replace n_loc_vc = 1 if n_loc_vc==.
replace n_loc_major_vc = 0 if n_loc_major_vc==.
forvalues i = 1/20 {
	replace loc`i'_vc = 0 if loc`i'_vc==.
}
foreach v of varlist d_* {
	replace `v' = 0 if `v'==.
}
order s size
order loc*, last
codebook s
sort s
count
save "$s2\\table_i[a]", replace


/* table_ts:
   t, s,
   cooperated,
   funding, human (table_ti aggregate)
*/

set more off

// (1) cooperated
use "$s1\\s1(c)", clear
ren s match
merge m:m match using "$s2\\table_0"
drop if _merge==2
ren match s
keep deal_no k t j
sort deal_no t k j
duplicates drop
egen mask_no = group(deal_no)
replace mask_no = mask_no-1
egen k_mask = group(k)
replace k_mask = k_mask-1 
egen min_t = min(t)
g t_mask = t - min_t
sort mask_no t k_mask
outsheet mask_no k_mask t_mask using "$cpp\\deals.txt", nonames replace
cd "$cpp"
shell "$cpp\\deals.exe"
preserve
insheet using "$cpp\\cooperated.txt", names clear
ren deal_no mask_no
tempfile 1
save "`1'", replace
restore
merge m:m mask_no using "`1'"
drop _merge
replace earliest_cooperated_time = . if earliest_cooperated_time==-1
replace earliest_cooperated_time = earliest_cooperated_time + min_t
replace earliest_cooperated_time = t if earliest_cooperated_time ==.
keep deal_no earliest_cooperated_time
merge m:m deal_no using "$s1\\s1(c)"
egen temp = min(earliest_cooperated_time), by(s)
drop earliest_cooperated_time
ren temp begin_t
keep s begin_t
duplicates drop
sort s
replace begin_t = begin_t+1
ren begin_t t
g cooperated = 1
tempfile 2
save "`2'", replace
keep s 
duplicates drop
g x = 1
tempfile s
save "`s'", replace
clear
loc last = (${endy}-${begy}+1)*12
set obs `last'
g t = _n-1
g x = 1
joinby x using "`s'"
merge 1:1 s t using "`2'"
drop _merge x
sort s t
bysort s: g temp = sum(cooperated)
drop cooperated
ren temp cooperated
save "$s2\\table_ts(1)", replace

// (2) funding (aggregate table_ti)
// # of funded rounds
use "$s2\\table_0", clear
keep if match~=-1
ren match s 
sort j s t
bysort j s: g n = _n
keep if n==1
drop n
merge m:m s using "$s1\\s1(c)"
drop _merge
keep t k j
order k t 
sort j k t
bysort j k: g n = _n
keep if n==1
drop n
collapse (count) j, by(k t)
xtset k t
count
loc last = r(N)+2
set obs `last'
loc last = r(N)+1
replace t = 0 in `last'
loc last = `last'+1
replace t = ym(${endy},12) - ym(${begy},1) in `last'
tsfill, full
drop if k==.
bysort k: g n_fund = sum(j)
drop j
joinby k using "$s1\\s1(c)"
collapse (median) n_fund, by(s t)
count
tempfile 1
save "`1'", replace
// # of funded categories
use "$s2\\table_0", clear
keep if match~=-1
ren match s 
merge m:m s using "$s1\\s1(c)"
drop _merge
ren k i
merge m:m j using "$s2\\j"
drop if _merge==2
drop _merge
merge m:m k using "$tb\\organization_categories"
drop if _merge==2
drop _merge
replace c = 0 if c==. // add the dummy category to missing ones
drop k 
order i c t
sort i c t
egen min_t = min(t), by(i c) // find the earliest funding time for a category
keep i c min_t
ren min_t t
duplicates drop
ren i k 
joinby k using "$s1\\s1(c)"
keep c s t
egen min_t = min(t), by(s c) // find the earliest funding time for a category for s
format min_t %tm
keep s c min_t
ren min_t t
duplicates drop
duplicates r s c
egen id = group(s c)
xtset id t
count
loc last = r(N)+2
set obs `last'
loc last = r(N)+1
replace t = 0 in `last'
loc last = `last'+1
replace t = ym(${endy},12) - ym(${begy},1) in `last'
tsfill, full
g fund = 1 if s~=. & c~=.
foreach v in s c {
	egen temp = mean(`v'), by(id)
	drop `v'
	ren temp `v'
}
order s c
bysort id: g cat_fund = sum(fund)
collapse (median) cat_fund, by(s t)
drop if s==.
merge m:m s t using "`1'"
keep if _merge==3
drop _merge
count
order t s n_fund cat_fund 
sort t s
duplicates drop
save "$s2\\table_ts(2)", replace

// (3) people (aggregate table_ti)
// # degrees
// simple maximum for people factors
use "$s2\\table_ti", clear
drop n_fund cat_fund
ren i k
joinby k using "$s1\\s1(c)"
drop deal_no k 
order s t
collapse (max) d*, by(t s)
order t s
sort t s
merge 1:1 t s using "$s2\\table_ts(1)"
keep if _merge==3
drop _merge
merge 1:1 t s using "$s2\\table_ts(2)"
keep if _merge==3
drop _merge
count
order t s cooperated n_fund cat_fund
foreach v in cooperated n_fund cat_fund {
	ren `v' `v'_vc
}
codebook s
count
save "$s2\\table_ti[a]", replace
rm "$s2\\table_ts(1).dta"
rm "$s2\\table_ts(2).dta"

/* =============================== aggregate table_ji and table_tji ============================ */

/* table_js:
   j, s,
   location (table_ji) aggregate
*/

set more off

use "$s1\\s1(c)", clear
keep s
duplicates drop
g x = 1
tempfile 1
save "`1'", replace
use "$s2\\table_0", clear
keep j
duplicates drop
g x = 1
joinby x using "`1'"
drop x
count
save "$s2\\s-j", replace

use "$s2\\i-locs", clear
keep i name_vc
ren name_vc name
merge m:m name using "$s2\\py\\locations_latlong"
keep if _merge==3
drop _merge name
ren latitude lat_vc
ren longtitude long_vc
append using "$s2\\py\\i-locs-add"
sort i
codebook i
ren i k
joinby k using "$s1\\s1(c)"
codebook k
codebook s
keep s lat_vc long_vc
sort s
order s
duplicates drop
outsheet s lat_vc long_vc using "$s2\\py\\i.csv", comma nonames replace
use "$s2\\j-locs", clear
keep j name_st
ren name_st name
merge m:m name using "$s2\\py\\locations_latlong"
keep if _merge==3
drop _merge name
ren latitude lat_st
ren longtitude long_st
append using "$s2\\py\\j-locs-add"
sort j 
codebook j
outsheet j lat_st long_st using "$s2\\py\\j.csv", comma nonames replace
// calculate_distance
insheet using "$s2\\py\\distance.csv", comma names clear
ren i s
merge m:m s j using "$s2\\s-j"
drop _merge
g travel = 0 if dist<=100 // 50 mph for driving
replace travel = 1 if dist<=1000 & travel==. // 500 mph for flight
replace travel = 2 if travel==. // over 1-day round trip
count
save "$s2\\table_ji[a]", replace

/* table_tjs:
   t, j, s
   funding, category, people (table_tji) aggregate
*/

set more off

// (1) funding

// all-funded
use "$s1\\s1(c)", clear
keep s k 
duplicates drop
ren s match
joinby match using "$s2\\table_0"
ren match s
keep s t j
duplicates drop
collapse (min) t, by (s j)
replace t = t+1
drop if t > ym(${endy},12) - ym(${begy},1)
qui count
loc last = r(N) + 1
set obs `last'
replace t = ym(${endy},12) - ym(${begy},1) in `last'
egen id = group(s j)
xtset id t
g all_funded = 1
tsfill, full
drop if id==.
foreach v in s j {
	egen temp = mean(`v'), by(id)
	drop `v'
	ren temp `v'
}
sort id t
bysort id: g temp = sum(all_funded)
drop all_funded
ren temp all_funded
drop id
drop if all_funded==0
duplicates drop
save "$s2\\table_tji_sparse(1a)[a]", replace

// max_funded
use "$s1\\s1(c)", clear
keep s k 
ren k i
duplicates drop
ren s match
joinby match using "$s2\\table_0"
ren match s
keep t i j
order t i j
collapse (min) t, by(i j)
tempfile 2
save "`2'", replace
use "$s1\\s1(c)", clear
keep s k 
duplicates drop
ren k i 
joinby i using "`2'"
collapse (min) t, by(s j)
replace t = t+1
drop if t > ym(${endy},12) - ym(${begy},1)
qui count
loc last = r(N) + 1
set obs `last'
replace t = ym(${endy},12) - ym(${begy},1) in `last'
egen id = group(s j)
xtset id t
g max_funded = 1
tsfill, full
drop if id==.
foreach v in s j {
	egen temp = mean(`v'), by(id)
	drop `v'
	ren temp `v'
}
sort id t
bysort id: g temp = sum(max_funded)
drop max_funded
ren temp max_funded
drop id
drop if max_funded==0
duplicates drop
save "$s2\\table_tji_sparse(1b)[a]", replace

// (2) category - abort for large dataset
use "$s1\\s1(c)", clear
keep s k 
duplicates drop
ren k i 
joinby i using "$s2\\table_tji_sparse(2)"
drop i 
collapse (sum) sum_same_cat = n_same_cat (max) max_same_cat = n_same_cat, by(t j s)
order t j s
duplicates drop
save "$s2\\table_tjs_sparse(2)[a]", replace

// (3) common people
use "$s1\\s1(c)", clear
drop deal_no
duplicates drop
ren k i
joinby i using "$s2\\table_tji_sparse(3)"
collapse (max) n_people, by(t j s)
order t j s
save "$s2\\table_tji_sparse(3)[a]", replace

// (4) common school - very large
use "$s2\\degrees", clear
cap drop n_degrees
drop d_engr - d_top20
g school = ""
order d*, last
loc j = 1
loc S1 = `j'
foreach v of varlist d_stanford - d_cmu {
	ren `v' d_school`j'
	loc j = `j'+1
}
loc S2 = `j' - 1
forvalues s = `S1'/`S2' {
	replace school = school+"" if d_school`s'==0
	replace school = school+",`s'" if d_school`s'==1
}
replace school = substr(school,2,length(school))
drop d_*
split school, p(",")
drop school 
foreach v in school* {
	destring `v', replace
}
egen n_school = rownonmiss(school*)
drop if n_school==0
expand n_school
bysort p: g n = _n
tostring n, replace
foreach v of varlist school* {
	replace `v' = . if regexm("`v'",n)==0
}
egen sch = rowmin(school*)
drop school*
keep p sch
drop if sch==.
duplicates drop
tempfile 1
save "`1'", replace
use "$s2\\i-people", clear
ren i k
joinby k using "$s1\\s1(c)"
keep s started_t ended_t p
duplicates drop
order s p
append using "$s2\\j-people"
sort s j p started_t ended_t
joinby p using "`1'"
egen t_sch1 = min(started_t), by(sch j)
egen t_sch2 = max(ended_t), by(sch j)
sum started_t
loc min_t = r(min)
foreach v of varlist t_* {
	replace `v' = `v' - `min_t'
}
egen sch_mask = group(sch)
replace sch_mask = sch_mask-1
egen s_mask = group(s) if s<.
replace s_mask = s_mask-1
egen j_mask = group(j) if s==.
replace j_mask = j_mask-1
duplicates drop
preserve
keep if s<.
keep s_mask sch_mask t_sch1 t_sch2 
order s_mask sch_mask t_sch1 t_sch2 
duplicates drop
sort s_mask sch_mask 
outsheet using "$cpp\\s-school.txt", nonames replace
restore
preserve
keep if s==.
keep j_mask sch_mask t_sch1 t_sch2
order j_mask sch_mask t_sch1 t_sch2
duplicates drop
sort j_mask sch_mask
outsheet using "$cpp\\j-school.txt", nonames replace
restore
// s:s_mask
preserve
keep if s<.
keep s s_mask
duplicates drop
tempfile s
save "`s'", replace
restore
// j:j_mask
preserve
keep if s==.
keep j j_mask
duplicates drop
tempfile j
save "`j'", replace
restore
// school-match
//shell "people.exe" "s-school.txt" "j-school.txt" "s-j school-match.txt"
insheet using "$cpp\\s-j school-match.txt", clear
ren v1 s_mask
ren v2 j_mask
ren v3 t_mask
ren v4 n_school
merge m:m s_mask using "`s'"
keep if _merge==3
drop _merge
merge m:m j_mask using "`j'"
keep if _merge==3
drop _merge
g min_t = `min_t'
replace t_mask = t_mask + min_t
keep s j t_mask n_school
ren t_mask t
order t s j n_school
save "$s2\\table_tji_sparse(4)[a]", replace
