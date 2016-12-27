/* =================================================================== */

// Filter Organized Database (startup, vc, vc-syndicate, ipo, ma, funding, last_event)

// Copyright (C) 2015 by Yun Ling

/* =================================================================== */

clear all
set mem 2g
gl tb "D:\Dropbox\Dissertation\Data\CrunchBase\Raw Table"
gl tbc "D:\Dropbox\Dissertation\Data\CrunchBase\Table"
gl ext "D:\Dropbox\Dissertation\Data\CrunchBase\external"
gl py "D:\Dropbox\Dissertation\Data\CrunchBase\py"
gl s1 "D:\Dropbox\Dissertation\Data\CrunchBase\Inputs\1 screening results"

// Time Horizon
gl begy = 1998
gl endy = 2014

// Frequency is monthly data
set more off
set seed 1

/* =============================== startup ============================ */

// (I) whole pool
use "$tb\\organization", clear
// (i) keep if primary_role == company
keep if primary_role == "company"
drop primary_role
// (ii) founded and closed ym in horizon
drop if founded_on == . 
drop if year(founded_on) < ${begy} | year(founded_on) > ${endy}
g founded_ym = ym(year(founded_on), month(founded_on))
format founded_ym %tm
drop if year(closed_on) < ${begy}
drop if year(closed_on) > ${endy} & closed_on < .
g closed_ym = ym(year(closed_on), month(closed_on))
format closed_ym %tm
// (iii) current_team, websites, categories, founders, headquarters
foreach v in current_team categories founders websites headquarters {
	keep if `v' > 0 
}
// (v) select variables
keep k path founded_ym closed_ym
count
save "$s1\\valid_startup", replace


// (II) ipo
use "$tb\\ipo", clear
// (i) ipo ym 
drop if went_public_on == .
g ipo_ym = ym(year(went_public_on), month(went_public_on))
format ipo_ym %tm
// (ii) market value
g ipo_mv = opening_valuation
replace ipo_mv = money_raised if ipo_mv==.
replace ipo_mv = opening_share_price * shares_sold if ipo_mv == .
// (iii) ticker for external validation
ren stock_symbol ipo_ticker
drop if ipo_ticker == "" & ipo_mv == .
// (v) select variables
keep k ipo_no ipo_ym ipo_ticker ipo_mv
count
// (vi) internal validation wtih "`valid_startup'"
merge m:m k using "$s1\\valid_startup"
keep if _merge==3
drop if ipo_ym <= founded_ym
drop if ipo_ym >= closed_ym
keep k ipo_no ipo_ym ipo_ticker ipo_mv
count
// (v) external information form SDC for ipo_mv
g x = (ipo_mv==.)
preserve
insheet using "$s1\\ipo\\IPO.csv", names clear
keep k ipo_mv
ren ipo_mv ipo_mv_2
tempfile 1
save "`1'"
restore
merge m:m k using "`1'"
replace ipo_mv = . if ipo_mv < 10000000
replace ipo_mv = ipo_mv_2 if ipo_mv==.
sort ipo_mv
drop ipo_mv_2 _merge x
save "$s1\\ipo", replace


// (III) ma
use "$tb\\acquisition", clear
// (i) ma ym
g ma_ym = ym(year(announced_on), month(announced_on))
format ma_ym %tm
// (ii) ma type
drop if acquisition_type == "LBO"
// (iii) ma status
drop if acquisition_status == "Cancelled"
// (iv) select variables
keep ma_no acquirer_k acquiree_k ma_ym price
// (v) get the names of acquirer and acquiree
foreach v in acquirer acquiree {
	ren `v'_k k
	merge m:m k using "$tb\\organization"
	keep if _merge==3
	drop _merge
	ren k `v'_k
	merge m:m uuid using "$tb\\organization_name"
	keep if _merge==3
	drop _merge uuid
	ren name `v'_name
}
duplicates t *_k, generate(x)
sort x acquiree_name
drop if _n==_N
keep ma_no - ma_ym *_name
// (vi) internal validation wtih "`valid_startup'"
ren acquiree_k k
merge m:m k using "$s1\\valid_startup"
keep if _merge==3
drop if ma_ym <= founded_ym
drop if ma_ym >= closed_ym
keep ma_no - ma_ym *_name
sort ma_ym
duplicates r k
count
// (v)  external information form SDC for ma_mv (hand-match acquiree name)
preserve
use "$ext\\sdc_ma", clear
order tgt_name ma_ym
outsheet using "$ext\\ma_right.csv", comma replace
restore
g x = (price<1000000) | (price==.)
foreach v in acquirer acquiree {
	replace `v'_name = trim(subinstr(lower(`v'_name),"!","",.))
}
preserve
drop if x == 0
order acquiree_name ma_ym
outsheet using "$ext\\ma_left.csv", comma replace
//shell python.exe "$py\\fuzmatch.py" "ma_left.csv" "ma_right.csv" "merge.csv" // python
insheet using "$ext\\merge.csv", comma names clear
foreach v in ma_ym ma_ym_m2 {
	cap drop y m
	g y = substr(`v',1,4)
	g m = substr(`v',6,2)
	destring y, replace
	destring m, replace
	drop `v'
	g `v' = ym(y,m)
	format `v' %tm
	cap drop y m
}
tempfile ma_2
save "`ma_2'"
restore
keep if x == 0 
append using "`ma_2'"
replace x = 0.5 if _merge2==3
replace x = 1 if _score1<0.000001 & _score1~=.
replace x = 1 if _score2<0.05 & _score2~=.
sort x // hand-screening and replace x = 1
order acquiree_name tgt_name_m1 tgt_name_m2 x
replace price = tgt_mv_m2 if x==0.5
replace x = 0 if x==0.5
drop *_m1 *_m2 _merge1 _merge2 _score1 _score2
replace x = (price==.)
order x 
sort x
order acquiree_name, after(acquirer_name)
drop x
ren price ma_mv
replace ma_mv = . if ma_mv==0
replace ma_mv = ma_mv * 1000000 if ma_mv < 1000
g x = (ma_mv == .)
preserve
insheet using "$s1\\ma\\MA.csv", comma names clear
keep k ma_mv_2
tempfile 1
save "`1'", replace
restore
merge m:m k using "`1'"
drop if _merge==2
replace ma_mv = ma_mv_2 if ma_mv == .
drop _merge ma_mv_2 x
sort k
replace ma_mv = . if ma_mv < 100000
count if ma_mv < .
save "$s1\\ma", replace

/* =============================== funding ============================ */

// (I) funded startup
use "$tb\\funding-round", clear
// (i) fund ym
drop if announced_on==.
g fund_ym = ym(year(announced_on), month(announced_on))
format fund_ym %tm
// (ii) fund type
keep if funding_type=="seed" | funding_type=="venture" | funding_type=="undisclosed" ///
	| funding_type=="convertible_note" | funding_type=="private_equity"
// (iii) select variables
keep deal_no uuid funded_k fund_ym series funding_type money_raised post_money_valuation
replace series = lower(series)
order deal_no fund_ym 
order uuid funding_type, last
ren money_raised invest
ren post_money_valuation post_mv
ren funding_type fund_type
// (iv) internal validation wtih "`valid_startup'"
ren funded_k k
merge m:m k using "$s1\\valid_startup"
keep if _merge==3
drop if fund_ym <= founded_ym
drop if fund_ym >= closed_ym
keep deal_no - fund_type 
drop if post_mv <= invest // drop if pre_mv<=0
count
// (v) get the names of startup
drop uuid
merge m:m k using "$tb\\organization"
keep if _merge==3
drop _merge
keep deal_no - fund_type uuid
merge m:m uuid using "$tb\\organization_name"
keep if _merge==3
keep deal_no - fund_type name
duplicates r deal_no
duplicates r fund_ym k
duplicates drop
// (vi) collapse (k, fund_ym)
duplicates t fund_ym k, generate(x)
gsort x k fund_ym -deal_no
bysort x k fund_ym: g n = _n
keep if n==1
duplicates r fund_ym k
sort name fund_ym
save "$s1\\funding", replace

/* =============================== vc ============================ */

// (I) vc
use "$tb\\funding-round_investors", clear
// (i) only vc
keep if organization==1
drop organization p
ren k investor_k
drop money_invested_currency_code
// (ii) have valid funding histories
merge m:m deal_no using "$s1\\funding"
keep if _merge==3
drop _merge
// (iii) select variables
keep deal_no investor_k fund_ym k
ren k funded_k
ren investor_k k
// (iv) get the name of vc
merge m:m k using "$tb\\organization"
keep if _merge==3
drop _merge
//keep if primary_role == "investor" // can be relaxed
keep deal_no - funded_k uuid founded_on
merge m:m uuid using "$tb\\organization_name"
keep if _merge==3
keep deal_no - funded_k name founded_on
g founded_ym = ym(year(founded_on), month(founded_on))
format founded_ym %tm
drop founded_on
// (v) number of investments vc make
egen N_investments = count(deal_no), by(k)
tab N_investments
keep if N_investments >= 10
save "$s1\\funding_vc", replace

/* =============================== combine startup, ipo, ma, funding, funding_vc first ============================ */

// build a delete subset of valid startup
use "$s1\\ipo", clear
keep if ipo_mv ==. // ipo_mv missing
keep k
g delete_from = "ipo missing"
preserve
use "$s1\\ma", clear
keep if ma_mv ==. // ma_mv missing
keep k
g delete_from = "ma missing"
tempfile 1
save "`1'", replace
restore
append using "`1'"
preserve
use "$s1\\funding", clear
ren k funded_k
merge m:m funded_k using "$s1\\funding_vc"
keep if _merge==1 | invest == . // missing investor or investment
keep funded_k
ren funded_k k
duplicates drop 
g delete_from = "fund missing"
tempfile 2
save "`2'", replace
restore
append using "`2'"
tab delete_from
qui count
loc N = r(N)+1
set obs `N'
replace k =  266628 if _n==_N
save "$s1\\delete_k", replace

// find valid startups
use "$s1\\delete_k", clear
merge m:m k using "$s1\\valid_startup"
keep if _merge==2
drop delete_from _merge
merge m:m k using "$s1\\ipo"
drop if _merge==2
drop _merge
merge m:m k using "$s1\\ma"
drop if _merge==2
drop _merge
merge m:m k using "$s1\\funding"
drop if _merge==2
egen fund_n = count(fund_ym) if _merge==3, by(k)
egen post_mv_n = count(post_mv) if _merge==3, by(k)
egen fund_ym_min = min(fund_ym) if _merge==3, by(k)
egen fund_ym_max = max(fund_ym) if _merge==3, by(k)
format fund_ym_min %tm
format fund_ym_max %tm
drop deal_no fund_ym series invest post_mv fund_type
duplicates drop
g funded = (_merge==3)
drop _merge
drop x
duplicates drop
// delete if founded_ym >= min(ipo_ym, ma_ym, closed_ym, fund_ym_min)
drop if founded_ym >= min(ipo_ym, ma_ym, closed_ym, fund_ym_min)
// delete if fund_ym_max >= min(ipo_ym, ma_ym, closed_ym)
drop if fund_ym_max >= min(ipo_ym, ma_ym, closed_ym) & fund_ym_max<.
// status unknown 
count if ipo_ym < .
count if ma_ym < .
count if closed_ym < .
g final_status = "ipo" if ipo_ym < . & ipo_ym < ma_ym 
replace final_status = "ma" if ma_ym < . & final_status==""
replace final_status = "closed" if closed_ym < . & final_status==""
codebook k
count

// add news to detect any walking dead
merge m:m k using "$tb\\organization_news"
drop if _merge==2
g posted_ym = ym(min(2014,year(posted_on)),month(posted_on))
egen posted_ym_max = max(posted_ym), by(k)
format posted_ym_max %tm
drop url _merge posted_ym posted_on
duplicates drop
codebook k
count

// add websites to detect final status
merge m:m k using "$tb\\organization_websites"
drop if _merge==2
drop _merge
g y=(title=="homepage")
egen ymax = max(y), by(k)
// no homepage then closed
replace final_status = "closed" if ymax==0 & final_status==""
replace url = "" if ymax==0
replace title = "" if ymax==0
keep if title=="homepage" | ymax==0
drop y ymax
duplicates drop
codebook k
count
preserve
keep url k funded
duplicates drop
order url k
sort k
drop if url==""
outsheet using "$ext\\check_websites.csv", comma names replace
//shell python.exe "$py\\check_websites.py" "check_websites.csv" "valid_websites.csv"
insheet using "$ext\\valid_websites.csv", comma names clear
egen valid_max = max(valid), by(k)
keep k valid_max
duplicates drop
ren valid_max valid
tempfile 1
save "`1'", replace
restore
merge m:m k using "`1'"
drop if _merge==2
drop _merge
drop url
duplicates drop
codebook k
count 
// replace to death if valid in [0,4]
replace valid = valid/1000
gl dead = 4
replace final_status = "closed" if valid<${dead} & final_status==""

// assign closed time
g last_ym = max(founded_ym, posted_ym_max, fund_ym_max) if final_status=="closed" & closed_ym==.
format last_ym %tm
// dead in 6-30 months since last_ym
replace closed_ym = min(ym(${endy},12),last_ym+floor(runiform()*25)+6) if final_status=="closed" & closed_ym==.
g closed_mv = 0.1 if final_status=="closed"
// remaining alive
replace final_status="alive" if final_status==""

// generate final_ym and final_mv
g final_ym = ipo_ym if final_status=="ipo"
replace final_ym = ma_ym if final_status=="ma"
replace final_ym = closed_ym if final_status=="closed"
replace final_ym = ym(${endy},12) if final_status=="alive"
g final_mv = ipo_mv if final_status=="ipo"
replace final_mv = ma_mv if final_status=="ma"
replace final_mv = closed_mv if final_status=="closed"
g no = ma_no if final_status=="ma"
replace no = ipo_no if final_status=="ipo"

// add description to detect growth
merge m:m k using "$tb\\organization"
drop if _merge==2
drop _merge images-primary_image funding_rounds-past_team headquarters-news board_* ipo-closed_on
codebook k
count

// summarize and take subsample
tab final_status
tab final_status if funded==1
count if funded==1
count if funded==0
// (a) delete some samples that's funded but not ipo/ma/dead
// delete from those delete~=0 & post_mv_n==0 (funded==1), delete ratio = 0.75
g delete = 0 if final_status~="alive" & funded==1
gl ratio = 0.75
qui sum fund_n
forvalues i = 1/8 {
	replace delete = (runiform() < ${ratio}) if fund_n==`i' & post_mv_n==0 & delete~=0 & funded==1
}
replace delete = 0 if delete~=1 & funded==1
tab post_mv_n fund_n
tab post_mv_n fund_n if delete==0
tab final_status if funded==1 & delete==0
// (b) delete some samples that're not funded but not ipo/ma/dead
// delete from those delete~=0 & funded==0, delete ratio = 0.6
replace delete = 0 if final_status~="alive" & funded==0
count if delete~=0 & funded==0
gl ratio = 0.6
replace delete = (runiform() < ${ratio}) if funded==0
replace delete = 0 if delete~=1 & funded==0
count if delete==0 & funded==0
tab delete funded
// summarize again
drop if delete==1
drop delete
tab final_status
tab final_status if funded==1
count if funded==1
count if funded==0

// double check
order k name path founded_ym final_status final_ym final_mv no last_ym ///
	funded valid current_team categories offices investments acquisitions description
keep k name path founded_ym final_status final_ym final_mv no last_ym ///
	funded valid current_team categories offices investments acquisitions description
tab final_status
tab final_status funded
codebook k
count
drop if k==92615
save "$s1\\s1(a)", replace

/* ---------------------- Funding, VC, and final_mv ------------------- */

// assign final_mv according to some descriptive characteristics
use "$s1\\s1(a)", clear
sum final_mv // biggest one is Alibaba's IPO
count if final_mv==.
count if final_status=="alive"
tab funded if final_status=="alive"
tab funded final_status
// funded
preserve
keep if funded==1
count
joinby k using "$s1\\funding"
codebook k
count
drop if invest==. | invest==0
replace invest = 400000 if deal_no==8223
replace invest = 230300 if deal_no==23247
drop x n
count if post_mv<.
sort k fund_ym deal_no
// drop duplicated rounds
drop if invest==invest[_n-1] & post_mv==post_mv[_n-1] & k == k[_n-1]
g status = "funding"
ren fund_ym ym
replace no = deal_no
drop deal_no 
order status ym no series invest post_mv fund_type, last
tempfile funding
save "`funding'", replace
restore
// append birth
preserve 
g status = "birth"
g ym = founded_ym
g mv = 1
order status ym mv, last
tempfile birth
save "`birth'", replace
restore
// append final
g status = final_status
g ym = final_ym
g mv = final_mv
order status ym mv, last
tempfile final
save "`final'", replace
// combine
use "`birth'", clear
append using "`funding'"
append using "`final'"
sort k ym
codebook k
duplicates tag k ym, gen(tag)
egen tag_max = max(tag), by(k)
drop if tag_max>0
drop tag tag_max
duplicates r k ym
drop name
count if mv==.
tab status if mv==.
tab status funded if mv==.
format ym %tm
count
g status_code = -2 if status=="birth"
replace status_code = 0 if status=="alive" | status=="funding"
replace status_code = 1 if status=="ipo" | status=="ma"
replace status_code = -1 if status=="closed"
count if post_mv<=invest & invest <.

// impute mv's
tab status if mv<.
tab status if mv==.
tab final_status if final_mv==. & funded ==1
// some post_mv <.
g exist = (post_mv<.)
egen temp = max(exist), by(k)
drop exist
ren temp exist
codebook k if final_status~="alive"
codebook k if exist==1
codebook k if exist==0 & final_status=="alive" & funded==1 // B
codebook k if final_status=="alive" & funded==0

// A. impute for funded ones with either (a) final_mv <. (final_status~="alive"), or (b) some post_mv <.
// - need to impute both status=="funding" and status=="alive"
// impute total_invest until final_ym
sort k ym
bysort k: g round = _n if status=="funding"
bysort k: g total_invest = sum(invest)
replace round = round-1
egen round_max = max(round), by(k)
g invest_speed = total_invest/(ym-founded_ym) if round==round_max
egen temp = mean(invest_speed), by(k)
drop invest_speed
ren temp invest_speed
replace total_invest = invest_speed*(ym-founded_ym) if ym==final_ym
drop invest_speed last_ym
// impute mv
tab status if mv==.
g imputed_mv = mv
replace imputed_mv = post_mv if mv==.
g imputed_mv_per_invest = imputed_mv/total_invest
egen temp = mean(imputed_mv_per_invest), by(k)
drop imputed_mv_per_invest
ren temp imputed_mv_per_invest
replace imputed_mv = imputed_mv_per_invest*total_invest if imputed_mv==.
count if final_status~="alive" & imputed_mv==.
count if exist==1 & imputed_mv==.
codebook k if (final_status~="alive" | exist==1) & funded==1 // A
codebook k if funded==1 & imputed_mv==. // B
codebook k if funded==1

// B. impute for funded ones with final_mv==. (final_status=="alive") and all post_mv ==.
// - need to impute both status=="funding" and status=="alive"
// use valid, current_team, offices, investments, acquisitions to gauge the well-being
codebook k if exist==0 & final_status=="alive" & funded==1 // B
foreach v in valid current_team offices investments acquisitions {
	sum `v', d
	sum `v' if status=="ipo" | status=="ma", d
}
g state=""
codebook k if funded==1 & final_status=="alive"
// (1) very well 
tab final_status if current_team>13 & funded==1
tab final_status if offices>5 & funded==1
tab final_status if investments>1 & funded==1
tab final_status if acquisitions>2  & funded==1
tab final_status funded if current_team>13 | offices>5 | investments>1 | acquisitions>2
replace state = "excellent" if state=="" & (current_team>13 | offices>5 | investments>1 | acquisitions>2) & funded==1
codebook k if state=="excellent" & funded==1 & final_status=="alive"
// (2) well 
tab final_status if current_team>9 & funded==1
tab final_status if offices>2 & funded==1
tab final_status if investments>0 & funded==1
tab final_status if acquisitions>0  & funded==1
tab final_status funded if current_team>9 | offices>2 | investments>0 | acquisitions>0
replace state = "well" if state=="" & (current_team>9 | offices>2 | investments>0 | acquisitions>0) & funded==1
codebook k if state=="well" & funded==1 & final_status=="alive"
// (3) survival
replace state = "survival" if state=="" & funded==1
codebook k if state=="survival" & funded==1 & final_status=="alive"
// impute values
g log_imputed_mv_per_invest = log(imputed_mv_per_invest) if final_status~="alive"
loc v log_imputed_mv_per_invest
foreach s in survival well excellent {
	sum `v' if state=="`s'" & status=="birth" & final_status~="closed"
	hist `v' if state=="`s'" & status=="birth" & final_status~="closed"
}
g temp = .
foreach s in survival well excellent {
	qui sum `v' if state=="`s'" & status=="birth" & final_status~="closed"
	loc mean = r(mean)
	loc sd = r(sd)/2
	loc min = r(min)
	loc max = r(max)
	replace temp = rnormal(`mean',`sd') if (exist==0 & final_status=="alive" & funded==1) & state=="`s'"
	replace temp = max(`min',temp) if (exist==0 & final_status=="alive" & funded==1) & state=="`s'"
	replace temp = min(`max',temp) if (exist==0 & final_status=="alive" & funded==1) & state=="`s'"
	replace imputed_mv_per_invest = exp(temp) if imputed_mv_per_invest==. & (exist==0 & final_status=="alive" & funded==1) & state=="`s'"
}
drop temp
egen temp = mean(imputed_mv_per_invest), by(k)
drop imputed_mv_per_invest
ren temp imputed_mv_per_invest
replace imputed_mv = imputed_mv_per_invest*total_invest if imputed_mv==. 
count if funded==1 & imputed_mv==.

// C. impute for unfunded ones with final_mv==. (final_status=="alive")
// - need to impute only status=="alive"
codebook k if status=="alive" & funded==0
// (1) excellent
tab final_status if current_team>23 & funded==0
tab final_status if offices>6 & funded==0
tab final_status if investments>3 & funded==0
tab final_status if acquisitions>5 & funded==0
tab final_status funded if current_team>23 | offices>6 | investments>3 | acquisitions>5
replace state = "excellent" if (current_team>23 | offices>6 | investments>3 | acquisitions>5) & funded==0
codebook k if state=="excellent" & status=="alive" & funded==0
// (2) well
tab final_status if current_team>18 & funded==0
tab final_status if offices>4 & funded==0
tab final_status if investments>2 & funded==0
tab final_status if acquisitions>1 & funded==0
tab final_status if (current_team>18 | offices>4 | investments>2 | acquisitions>1) & funded==0
replace state = "well" if state=="" & (current_team>18 | offices>4 | investments>2 | acquisitions>1) & funded==0
codebook k if state=="well" & status=="alive" & funded==0
// (3) average
replace state = "survival" if state=="" & funded==0
codebook k if state=="survival" & status=="alive" & funded==0
// impute values
tab final_status state if funded==0
g log_final_mv = log(final_mv) if final_status~="alive" & funded==0
codebook k if (final_status=="ipo" | final_status=="ma") & funded==0
codebook k if final_status=="alive" & funded==0
codebook k if final_status=="closed" & funded==0
loc v log_final_mv
g temp = .
foreach s in survival well excellent {
	sum `v' if state=="`s'" & status=="birth" & funded==0
	loc mean = r(mean)
	loc sd = r(sd)/2
	loc min = r(min)
	loc max = r(max)
	replace temp = rnormal(`mean',`sd') if (status=="alive" & funded==0) & state=="`s'"
	replace temp = max(`min',temp) if (status=="alive" & funded==0) & state=="`s'"
	replace temp = min(`max',temp) if (status=="alive" & funded==0) & state=="`s'"
	replace final_mv = exp(temp) if final_mv==. & (status=="alive" & funded==0) & state=="`s'"
}
drop temp
replace imputed_mv = final_mv if status=="alive" & funded==0
count if imputed_mv==.
drop state imputed_mv_per_invest log_final_mv
sum imputed_mv if status=="alive" & funded==0

// generate cumulative return
g ret = log(imputed_mv)
//replace ret = -(rbeta(3,3)*2+2) if status=="closed"
tab status if ret==.
sum ret
replace ret = -999 if ret==.
order k ym status_code ret 
codebook k
codebook k if funded==1
codebook k if funded==0
duplicates drop
save "$s1\\s1(b)", replace

/* =============================== vc syndicate ============================ */

// (I) vc syndicate
use "$s1\\s1(b)", clear
keep if status=="funding"
ren no deal_no
keep deal_no 
codebook deal_no
count 
joinby deal_no using "$s1\\funding_vc"
duplicates drop
egen size = count(k), by(deal_no)
sort deal_no k
keep deal_no k
bysort deal_no: g n = _n
keep deal_no k n
reshape wide k, i(deal_no) j(n)
egen s = group(k*), missing
replace s = s-1
reshape long k, i(deal_no s) j(n)
drop if k==.
drop n
sort s deal_no k
codebook s
codebook deal_no
save "$s1\\s1(c)", replace

use "$s1\\s1(b)", clear
g deal_no = no if status=="funding"
merge m:m deal_no using "$s1\\s1(c)"
count if status=="funding" & _merge==1
codebook k if status=="funding"
codebook k if status=="funding" & _merge==3
g checked = 1 if status=="funding" & _merge==3
egen temp = mean(checked), by(k)
drop checked
ren temp checked
codebook k if checked==1
codebook k if checked==. & funded==1
drop if checked==. & funded==1
codebook k if funded==1
codebook k if funded==0
drop if status=="funding" & _merge==1
drop _merge checked
count if s==. & status=="funding"
sort k ym
order k ym status_code ret s
duplicates drop
save "$s1\\s1(b)", replace
