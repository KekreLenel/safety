* computeMoments.do
* Computes all moments in sections 4 through 6 of paper and organizes data for 
*  plotData.m and estimateIRFs.m;
#delimit ;

clear all;
set maxvar 10000;

timer on 1;

local direc "..\Empirical";

scalar startyear = 1995;
scalar startyear_nfa = 2006; 
scalar startquarter = 1;
scalar endyear = 2019;
scalar endquarter = 4;
scalar endmonth = 12;

matrix Table1 = J(5,1,0);
matrix colnames Table1 = "Data";
matrix rownames Table1 = 
 "rho(p)"
 "sigma(rho)/E(rho)"
 "skew(omega)"
 "rho(omega)" 
 "rho(p,omega)";
matrix Table2 = J(15,1,0);
matrix colnames Table2 = "Data";
matrix rownames Table2 =
 "y*/(sy)"
 "sigma(dlog c)"
 "sigma(dlog y*)"
 "rho(y*/y,L4.y*/y)"
 "sigma(dlog x)"
 "4r"
 "nfa/(4y)"
 "rho-1(re,rF-r)"
 "4(re-r)"
 "b(dnfa/y,re-r)"
 "bHs*/(4y)"
 "l"
 "l*"
 "log P/P-1"
 "log P*/P*-1";
matrix Table3 = J(6,1,.);
matrix colnames Table3 = "Data";
matrix rownames Table3 =
 "b(rF-r,log y/y-4)"
 ""
 "b(rF-r,re)"
 ""
 "b(dnfa/y,rF-r)"
 "";
matrix Table4 = J(6,1,0);
matrix colnames Table4 = "Data";
matrix rownames Table4 =
 "sigma(4r)"
 "sigma(4(re-r))"
 "sigma(4(rF-r))"
 "sigma(dlog q)"
 "sigma(dlog E)"
 "rho(dlog q,dlogc*/c)";
matrix Table5 = J(2,1,0);
matrix colnames Table5 = "Data";
matrix rownames Table5 =
 "sigma(dlog y)"
 "sigma(dlog y*)";
matrix Table6 = J(6,1,0);
matrix colnames Table6 = "Data";
matrix rownames Table6 =
 "sigma(dnfa/y)"
 "sigma(nx/y)"
 "sigma((dnfa-nx)/y)" 
 "dnfa/y"
 "nx/y"
 "(dnfa-nx)/y";
matrix Table11 = J(6,2,.);
matrix colnames Table11 = "re" "rF";
matrix rownames Table11 =
 "L3.dp"
 ""
 "L3.i* L3.i" 
 ""
 "log L3.y/L15.y"
 ""; 

********** Annualized spreads vs U.S. Treasuries  **********;

use "`direc'\data\Monthly_Data.dta", replace; 
qui gen yearmonth = year + (month-1)/12;
qui egen cip_govt_3m_2020 = rowmean(cip_govt_AUD_3m_2020 cip_govt_CAD_3m_2020 cip_govt_CHF_3m_2020 
 cip_govt_DKK_3m_2020 cip_govt_EUR_3m_2020 cip_govt_GBP_3m_2020 cip_govt_JPY_3m_2020 cip_govt_NOK_3m_2020 cip_govt_NZD_3m_2020 cip_govt_SEK_3m_2020);
qui gen dcpnf3m_y = 100*(dcpnf3m - y_USD_3m_2020);
export delimited yearmonth year month cip_govt_3m_2020 dcpnf3m_y using "`direc'\output\fig_1.csv" if year >= 1997 & yearmonth <= 2020.5, replace;

********** Second moments of disaster risk and safety shocks **********;

use "`direc'\data\Monthly_Data.dta", replace;
qui egen cip_govt_3m_2020 = rowmean(cip_govt_AUD_3m_2020 cip_govt_CAD_3m_2020 cip_govt_CHF_3m_2020 
 cip_govt_DKK_3m_2020 cip_govt_EUR_3m_2020 cip_govt_GBP_3m_2020 cip_govt_JPY_3m_2020 cip_govt_NOK_3m_2020 cip_govt_NZD_3m_2020 cip_govt_SEK_3m_2020);

qui correl p_all_bl L3.p_all_bl if ym >= ym(startyear,startquarter*3) & ym <= ym(endyear,endquarter*3); 
matrix Table1[1,1] = el(r(C),2,1);

qui tabstat p_all_bl if (month == 3 | month == 6 | month == 9 | month == 12) & ym >= ym(startyear,startquarter*3) & ym <= ym(endyear,endquarter*3), statistics(mean sd) save; 
matrix Table1[2,1] = el(r(StatTotal),2,1)/el(r(StatTotal),1,1);

qui tabstat cip_govt_3m_2020 if (month == 3 | month == 6 | month == 9 | month == 12) & ym >= ym(startyear,startquarter*3) & ym <= ym(endyear,endquarter*3), statistics(skewness) save; 
matrix Table1[3,1] = el(r(StatTotal),1,1);

qui correl cip_govt_3m_2020 L3.cip_govt_3m_2020 if ym >= ym(startyear,startquarter*3) & ym <= ym(endyear,endquarter*3); 
matrix Table1[4,1] = el(r(C),2,1);

qui correl cip_govt_3m_2020 p_all_bl if ym >= ym(startyear,startquarter*3) & ym <= ym(endyear,endquarter*3); 
matrix Table1[5,1] = el(r(C),2,1);

********** Conditional and unconditional relationship between excess Foreign bond return and equity **********;

use "`direc'\data\Monthly_Data.dta", replace;
qui gen AUD_eom = 1/inv_AUD_eom;
qui gen EUR_eom = 1/inv_EUR_eom;
qui gen GBP_eom = 1/inv_GBP_eom;
qui gen NZD_eom = 1/inv_NZD_eom;

foreach currency in AUD CAD CHF DKK EUR GBP JPY NOK NZD SEK {;
	qui gen dlog_`currency'_eom = log(`currency'_eom)-log(L3.`currency'_eom);
	qui replace diffy_`currency'_3m_2020 = diffy_`currency'_3m_2020/4;
	qui gen uip_`currency'_3m = L3.diffy_`currency'_3m_2020 - 100*dlog_`currency'_eom;
};
qui egen diffy_3m_2020 = rowmean(diffy_AUD_3m_2020 diffy_CAD_3m_2020 diffy_CHF_3m_2020
	diffy_DKK_3m_2020 diffy_EUR_3m_2020 diffy_GBP_3m_2020 diffy_JPY_3m_2020 diffy_NOK_3m_2020 diffy_NZD_3m_2020 diffy_SEK_3m_2020);
qui egen uip_3m = rowmean(uip_AUD_3m uip_CAD_3m uip_CHF_3m 
	uip_DKK_3m uip_EUR_3m uip_GBP_3m uip_JPY_3m uip_NOK_3m uip_NZD_3m uip_SEK_3m);
	
qui gen realequityreturn3mo = 100*(log(msci_tri_gross) - log(L3.msci_tri_gross) - (log(cpiurban_s_upd) - log(L3.cpiurban_s_upd))); 
	
qui gen d12indpro = indpro - L12.indpro;

// Table 11;
qui ivreg2 realequityreturn3mo L3.msci_dyield_gross if ym >= ym(startyear,3*startquarter) & ym <= ym(endyear,3*endquarter), kernel(tru) bw(4);
matrix Table11[1,1] = el(e(b),1,1);
matrix Table11[2,1] = sqrt(el(e(V),1,1));
qui predict resid_e_base, residuals;
	
qui ivreg2 uip_3m L3.diffy_3m_2020 L3.d12indpro if ym >= ym(startyear,3*startquarter) & ym <= ym(endyear,3*endquarter), kernel(tru) bw(4);
matrix Table11[3,2] = el(e(b),1,1);
matrix Table11[4,2] = sqrt(el(e(V),1,1));
matrix Table11[5,2] = el(e(b),1,2);
matrix Table11[6,2] = sqrt(el(e(V),2,2));
qui predict resid_F_base, residuals;

// Table 2;
qui correl resid_e_base resid_F_base;
matrix Table2[8,1] = el(r(C),2,1);

// Table 3;
qui ivreg2 uip_3m realequityreturn3mo if ym >= ym(startyear,3*startquarter) & ym <= ym(endyear,3*endquarter), kernel(tru) bw(4);
matrix Table3[3,1] = el(e(b),1,1);
matrix Table3[4,1] = sqrt(el(e(V),1,1));
	
// Figure 8;
qui gen eF_base = resid_e_base*resid_F_base;
qui gen yearmonth = year + (month-1)/12;
export delimited yearmonth year month eF_base using "`direc'\output\fig_8.csv" if ym >= ym(startyear,3*startquarter) & ym <= ym(endyear,3*endquarter), replace;

********** Relative size **********;

use "`direc'\data\Annual_Data.dta", replace;
qui gen gdp_usG10 = gdp_usAT + gdp_usAU + gdp_usBE + gdp_usCA + gdp_usCH + 
 gdp_usDE + gdp_usDK + gdp_usES + gdp_usFI + gdp_usFR + gdp_usGB + 
 gdp_usGR + gdp_usIE + gdp_usIT + gdp_usJP + gdp_usLU + gdp_usNL + 
 gdp_usNO + gdp_usNZ + gdp_usPT + gdp_usSE;
qui gen rel_gdp_g10 = gdp_usG10/gdp_usUS;
qui tabstat rel_gdp_g10 if year >= startyear & year <= endyear, statistics(mean) save;
matrix Table2[1,1] = el(r(StatTotal),1,1);

********** Second moments of consumption, investment, and output **********;

// Table 2, U.S. consumption and investment, and Table 5, U.S. output;
use "`direc'\data\Quarterly_Data.dta", replace;	
qui gen gdpdef = us_gdp/us_rgdp;
qui gen y_pc = us_rgdp/us_pop;
qui gen dlog_y_pc = log(y_pc) - log(L.y_pc);
qui gen c_pc = (us_cn + us_cs)/(us_pop*gdpdef);
qui gen dlog_c_pc = log(c_pc) - log(L.c_pc);
qui gen x_pc = (us_fi + us_cd)/(us_pop*gdpdef);
qui gen dlog_x_pc = log(x_pc) - log(L.x_pc);
qui tabstat dlog_c_pc dlog_x_pc dlog_y_pc if yq >= yq(startyear,startquarter) & yq <= yq(endyear,endquarter), statistics(sd) save;
matrix Table2[2,1] = el(r(StatTotal),1,1);
matrix Table2[5,1] = el(r(StatTotal),1,2);	
matrix Table5[1,1] = el(r(StatTotal),1,3);

// Table 2, G10 output and G10/U.S. output;
foreach currency in USD AUD CAD CHF DKK EUR GBP JPY NOK NZD SEK{;
	qui gen gdpdef_`currency' = gdp_`currency'/rgdp_`currency';
	qui gen y_pc_`currency' = rgdp_`currency'/pop_`currency';
	qui gen log_y_pc_`currency' = log(y_pc_`currency');
	qui gen dlog_y_pc_`currency' = log_y_pc_`currency' - L.log_y_pc_`currency';
	qui gen rel_y_pc_`currency' = y_pc_`currency'/y_pc_USD;
	qui regress rel_y_pc_`currency' yq;
	qui _predict rel_y_pc_`currency'_dt, residuals;	
};	
qui tabstat dlog_y_pc_AUD dlog_y_pc_CAD dlog_y_pc_CHF dlog_y_pc_DKK dlog_y_pc_EUR 
 dlog_y_pc_GBP dlog_y_pc_JPY dlog_y_pc_NOK dlog_y_pc_NZD dlog_y_pc_SEK 
 if yq >= yq(startyear,startquarter) & yq <= yq(endyear,endquarter), statistics(sd) save;
matrix Table2[3,1] = (el(r(StatTotal),1,1) + el(r(StatTotal),1,2) + el(r(StatTotal),1,3)
 + el(r(StatTotal),1,4) + el(r(StatTotal),1,5) + el(r(StatTotal),1,6) + el(r(StatTotal),1,7)
 + el(r(StatTotal),1,8) + el(r(StatTotal),1,9) + el(r(StatTotal),1,10))/10;
matrix Table5[2,1] = el(Table2,3,1);
 
foreach currency in AUD CAD CHF DKK EUR GBP JPY NOK NZD SEK{;
	qui correl rel_y_pc_`currency'_dt L4.rel_y_pc_`currency'_dt if yq >= yq(startyear,startquarter) & yq <= yq(endyear,endquarter);
	matrix Table2[4,1] = el(Table2,4,1)+el(r(C),2,1)/10;
};	

********** First and second moments of asset prices **********; 
 
// Table 2 and Table 4, asset returns;
use "`direc'\data\Monthly_Data.dta", replace;
qui gen realyield3mo = L3.yield3mo_upd_eom/100 - 4*(log(cpiurban_s_upd) - log(L3.cpiurban_s_upd));
qui gen equityreturn3mo = 4*(log(msci_tri_gross) - log(L3.msci_tri_gross));
qui gen excessequity3mo = equityreturn3mo - L3.yield3mo_upd_eom/100;
qui gen AUD_eom = 1/inv_AUD_eom;
qui gen EUR_eom = 1/inv_EUR_eom;
qui gen GBP_eom = 1/inv_GBP_eom;
qui gen NZD_eom = 1/inv_NZD_eom;	
foreach currency in AUD CAD CHF DKK EUR GBP JPY NOK NZD SEK {;
	qui gen excessForeign3mo_`currency' = L3.diffy_`currency'_3m_2020/100 - 4*(log(`currency'_eom) - log(L3.`currency'_eom));
};
qui egen excessForeign3mo = rowmean(excessForeign3mo_AUD excessForeign3mo_CAD excessForeign3mo_CHF
	excessForeign3mo_DKK excessForeign3mo_EUR excessForeign3mo_GBP excessForeign3mo_JPY
	excessForeign3mo_NOK excessForeign3mo_NZD excessForeign3mo_SEK);

qui tabstat excessequity3mo if ym >= ym(startyear,startquarter*3) & ym <= ym(endyear,endquarter*3), statistics(mean) save;
matrix Table2[9,1] = el(r(StatTotal),1,1);

qui tabstat realyield3mo excessequity3mo excessForeign3mo if (month == 3 | month == 6 | month == 9 | month == 12) & ym >= ym(startyear,startquarter*3) & ym <= ym(endyear,endquarter*3), statistics(sd) save;
matrix Table4[1,1] = el(r(StatTotal),1,1);
matrix Table4[2,1] = el(r(StatTotal),1,2);
matrix Table4[3,1] = el(r(StatTotal),1,3);
 
foreach currency in AUD CAD CHF DKK EUR GBP JPY NOK NZD SEK {;
	qui gen dlog_`currency'_eom = log(`currency'_eom)-log(L3.`currency'_eom);
};
qui egen dlog_e_eom = rowmean(dlog_AUD_eom dlog_CAD_eom dlog_CHF_eom 
 dlog_DKK_eom dlog_EUR_eom dlog_GBP_eom dlog_JPY_eom dlog_NOK_eom dlog_NZD_eom dlog_SEK_eom);
qui tabstat dlog_e_eom if (month == 3 | month == 6 | month == 9 | month == 12) & ym >= ym(startyear,startquarter*3) & ym <= ym(endyear,endquarter*3), statistics(sd) save;
matrix Table4[4,1] = el(r(StatTotal),1,1);

qui gen dlog_cpi_oecd_USD = log(cpi_oecd_USD)-log(L3.cpi_oecd_USD);
foreach currency in CAD CHF DKK EUR GBP JPY NOK SEK {;
	qui gen dlog_cpi_oecd_`currency' = log(cpi_oecd_`currency')-log(L3.cpi_oecd_`currency');
	qui gen dlog_q_`currency'_eom = dlog_`currency'_eom + (dlog_cpi_oecd_USD - dlog_cpi_oecd_`currency');
};	
qui keep if month == 3 | month == 6 | month == 9 | month == 12;
qui gen quarter = month/3;
qui gen yq = yq(year,quarter);
format yq %tq;
qui drop month;
qui keep yq year quarter dlog_e_eom dlog_q_CAD_eom dlog_q_CHF_eom dlog_q_DKK_eom dlog_q_EUR_eom 
 dlog_q_GBP_eom dlog_q_JPY_eom dlog_q_NOK_eom dlog_q_SEK_eom dlog_AUD_eom dlog_NZD_eom dlog_cpi_oecd_USD;
qui merge 1:1 yq using "`direc'\data\Quarterly_Data.dta";
tsset yq;
foreach currency in AUD NZD {;
	qui gen dlog_cpi_oecd_`currency' = log(cpi_oecd_`currency')-log(L1.cpi_oecd_`currency');
	qui gen dlog_q_`currency'_eom = dlog_`currency'_eom + (dlog_cpi_oecd_USD - dlog_cpi_oecd_`currency');
};
qui egen dlog_q_eom = rowmean(dlog_q_AUD_eom dlog_q_CAD_eom dlog_q_CHF_eom 
 dlog_q_DKK_eom dlog_q_EUR_eom dlog_q_GBP_eom dlog_q_JPY_eom dlog_q_NOK_eom dlog_q_NZD_eom dlog_q_SEK_eom);	

qui tabstat dlog_q_eom if yq >= yq(startyear,startquarter) & yq <= yq(endyear,endquarter), statistics(sd) save;
matrix Table4[5,1] = el(r(StatTotal),1,1);	

qui gen gdpdef = us_gdp/us_rgdp;
qui gen c_pc = (us_cn + us_cs)/(us_pop*gdpdef);
qui gen dlog_c_pc = log(c_pc) - log(L.c_pc);

foreach currency in CAD DKK EUR GBP JPY NOK NZD SEK{;
	qui gen gdpdef_`currency' = gdp_`currency'/rgdp_`currency';
	qui gen c_pc_`currency' = (cn_`currency' + cs_`currency')/(pop_`currency'*gdpdef_`currency');
	qui gen dlog_c_pc_`currency' = log(c_pc_`currency') - log(L.c_pc_`currency');
	
	qui gen ddlog_c_pc_`currency' = dlog_c_pc_`currency' - dlog_c_pc;	
	qui correl dlog_q_`currency'_eom ddlog_c_pc_`currency' if yq >= yq(startyear,startquarter) & yq <= yq(endyear,endquarter);
	matrix Table4[6,1] = Table4[6,1] + el(r(C),2,1)/8;
};	

********** First and second moments of NFA **********;

// Table 2, average NFA/GDP;
use "`direc'\data\Quarterly_Data.dta", replace;
qui keep if quarter == 4;
qui keep year us_gdp;
qui merge 1:1 year using "`direc'\data\Annual_Data.dta";
tsset year;
qui gen us_iip_ann = us_iip;
qui gen us_iip_gdp_ann = us_iip/(1000*us_gdp);
qui tabstat us_iip_gdp_ann if year >= startyear & year <= endyear, statistics(mean) save; 
matrix Table2[7,1] = el(r(StatTotal),1,1);

qui gen quarter = 4;
qui gen yq = yq(year,quarter);
tsset yq;
qui keep yq us_iip_ann us_iip_gdp_ann;
tempfile iip_gdp_ann;
save `iip_gdp_ann';

// Tables 2 and 3, projecting change in NFA/GDP on asset returns;
use "`direc'\data\Monthly_Data.dta", replace;
qui gen AUD_eom = 1/inv_AUD_eom;
qui gen EUR_eom = 1/inv_EUR_eom;
qui gen GBP_eom = 1/inv_GBP_eom;
qui gen NZD_eom = 1/inv_NZD_eom;
qui keep if month == 3 | month == 6 | month == 9 | month == 12;
qui gen yq = qofd(dofm(ym));
format yq %tq;
qui keep yq diffy*_2020 AUD_eom CAD_eom CHF_eom DKK_eom EUR_eom GBP_eom JPY_eom NOK_eom NZD_eom SEK_eom
 msci_tri_gross yield3mo_upd_eom cpi_oecd_USD;
order yq;
qui merge 1:1 yq using "`direc'\data\Quarterly_Data.dta";
drop _merge;
tsset yq;
qui egen diffy_3m_2020 = rowmean(diffy_AUD_3m_2020 diffy_CAD_3m_2020 diffy_CHF_3m_2020
 diffy_DKK_3m_2020 diffy_EUR_3m_2020 diffy_GBP_3m_2020 diffy_JPY_3m_2020 diffy_NOK_3m_2020 diffy_NZD_3m_2020 diffy_SEK_3m_2020);
qui replace diffy_3m_2020 = diffy_3m_2020/100;
qui gen dlog_cpi_USD = log(cpi_oecd_USD) - log(L.cpi_oecd_USD);
foreach currency in AUD CAD CHF DKK EUR GBP JPY NOK NZD SEK {;	
	qui gen dlog_`currency'_eom = log(`currency'_eom)-log(L1.`currency'_eom);
	qui gen excessForeign_`currency' = L1.diffy_`currency'_3m_2020/400 - dlog_`currency'_eom; 
};
qui egen excessForeign = rowmean(excessForeign_AUD excessForeign_CAD excessForeign_CHF excessForeign_DKK
	 excessForeign_EUR excessForeign_GBP excessForeign_JPY excessForeign_NOK excessForeign_NZD excessForeign_SEK);
qui gen excessequity = (log(msci_tri_gross) - log(L.msci_tri_gross)) - L.yield3mo_upd_eom/400 - dlog_cpi_USD;

qui gen us_iip_gdp = us_iip/(1000*us_gdp/4);
qui gen dus_iip_gdp = (us_iip-L.us_iip)/(1000*us_gdp/4);
qui gen us_tb_gdp = (us_ex_trans - us_im_trans)/(1000*us_gdp/4);
qui gen resid_gdp = dus_iip_gdp - us_tb_gdp;

qui regress dus_iip_gdp excessequity excessForeign if yq >= yq(startyear_nfa,startquarter) & yq <= yq(endyear,endquarter);

matrix Table2[10,1] = el(e(b),1,1);
matrix Table3[5,1] = el(e(b),1,2);
matrix Table3[6,1] = sqrt(el(e(V),2,2));

// Table 6;
qui tabstat dus_iip_gdp us_tb_gdp resid_gdp if yq >= yq(startyear_nfa,startquarter) & yq <= yq(endyear,endquarter), statistics(sd) save;
matrix Table6[1,1] = el(r(StatTotal),1,1); 
matrix Table6[2,1] = el(r(StatTotal),1,2); 
matrix Table6[3,1] = el(r(StatTotal),1,3);

qui merge 1:1 yq using "`iip_gdp_ann'", nogen;
qui gen dus_iip_gdp_ann = (us_iip_ann - L4.us_iip_ann)/(1000*us_gdp);
qui gen us_tb_gdp_ann = (us_ex_trans + L1.us_ex_trans + L2.us_ex_trans + L3.us_ex_trans - us_im_trans - L1.us_im_trans - L2.us_im_trans - L3.us_im_trans)/(1000*us_gdp) if quarter(dofq(yq)) == 4;
qui gen resid_gdp_ann = dus_iip_gdp_ann - us_tb_gdp_ann;

qui tabstat dus_iip_gdp_ann us_tb_gdp_ann resid_gdp_ann if yq >= yq(startyear,startquarter) & yq <= yq(endyear,endquarter), statistics(mean) save;
matrix Table6[4,1] = el(r(StatTotal),1,1); 
matrix Table6[5,1] = el(r(StatTotal),1,2); 
matrix Table6[6,1] = el(r(StatTotal),1,3);

********** Tbills held by foreigners plus central bank liquidity swaps **********;

use "`direc'\data\Monthly_Data.dta", replace;
qui keep if month == 12;
qui gen quarter = month/3;
qui gen yq = yq(year,quarter);
format yq %tq;
qui keep yq Tbill_foreign swpt;
qui merge 1:1 yq using "`direc'\data\Quarterly_Data.dta", nogen;
* qui keep if quarter == 4;
qui gen Tbill_foreign_swpt_gdp = (Tbill_foreign+swpt)/(1000*us_gdp);
qui tabstat Tbill_foreign_swpt_gdp if yq >= yq(startyear,startquarter) & yq <= yq(endyear,endquarter), statistics(mean) save;
matrix Table2[11,1] = el(r(StatTotal),1,1);

********** Predictability of dollar carry trade **********;

use "`direc'\data\Monthly_Data.dta", replace;
qui gen AUD_eom = 1/inv_AUD_eom;
qui gen EUR_eom = 1/inv_EUR_eom;
qui gen GBP_eom = 1/inv_GBP_eom;
qui gen NZD_eom = 1/inv_NZD_eom;
foreach currency in AUD CAD CHF DKK EUR GBP JPY NOK NZD SEK {;
	qui gen dF3log_`currency'_eom = log(F3.`currency'_eom)-log(`currency'_eom);
	qui gen uip_`currency'_3m = diffy_`currency'_3m_2020/4 - 100*dF3log_`currency'_eom; 
};
qui egen uip_3m = rowmean(uip_AUD_3m uip_CAD_3m uip_CHF_3m 
 uip_DKK_3m uip_EUR_3m uip_GBP_3m uip_JPY_3m uip_NOK_3m uip_NZD_3m uip_SEK_3m);

qui gen d12indpro = indpro - L12.indpro;
qui ivreg2 uip_3m d12indpro if ym >= ym(startyear,startquarter*3) & ym <= ym(endyear,endquarter*3), kernel(tru) bw(4);

matrix Table3[1,1] = el(e(b),1,1);
matrix Table3[2,1] = sqrt(el(e(V),1,1));

********** Add additional exogenous targets / normalizations to Table 2 **********;

matrix Table2[6,1] = 0.02;
matrix Table2[12,1] = 1;
matrix Table2[13,1] = 1;
matrix Table2[14,1] = 0;
matrix Table2[15,1] = 0;

********** Save all Tables to Excel **********;

putexcel set "`direc'\output\table_1.xlsx", replace;
putexcel A1 = matrix(Table1), rownames;

putexcel set "`direc'\output\table_2.xlsx", replace;
putexcel A1 = matrix(Table2), rownames;

putexcel set "`direc'\output\table_3.xlsx", replace;
putexcel A1 = matrix(Table3), rownames;

putexcel set "`direc'\output\table_4.xlsx", replace;
putexcel A1 = matrix(Table4), rownames;

putexcel set "`direc'\output\table_5.xlsx", replace;
putexcel A1 = matrix(Table5), rownames;

putexcel set "`direc'\output\table_6.xlsx", replace;
putexcel A1 = matrix(Table6), rownames;

putexcel set "`direc'\output\table_11.xlsx", replace;
putexcel A1 = matrix(Table11), rownames;

********** Save data for recursive VAR **********;

use "`direc'\data\Monthly_Data.dta", replace;
qui gen equityreturn = 100*(log(msci_tri_gross) - log(L1.msci_tri_gross));
qui gen excessequity = equityreturn - L1.yield1mo_upd_eom/12;
qui gen equityreturn3mo = 100*(log(msci_tri_gross) - log(L3.msci_tri_gross));
qui gen excessequity3mo = equityreturn3mo - L3.yield3mo_upd_eom/4;

* now create series for countries with monthly IP data available;
qui gen EUR_eom = 1/inv_EUR_eom;
qui gen GBP_eom = 1/inv_GBP_eom;
qui gen NZD_eom = 1/inv_NZD_eom;	
qui gen dlog_cpi_oecd_USD = 100*(log(cpi_oecd_USD)-log(L1.cpi_oecd_USD));
foreach currency in CAD DKK EUR GBP JPY NOK SEK {;
	qui gen dlog_`currency'_eom = 100*(log(`currency'_eom) - log(L1.`currency'_eom));
	qui gen dlog_cpi_oecd_`currency' = 100*(log(cpi_oecd_`currency')-log(L1.cpi_oecd_`currency'));
	qui gen dlog_q_`currency'_eom = dlog_`currency'_eom + (dlog_cpi_oecd_USD - dlog_cpi_oecd_`currency');	
};
qui egen dlog_e_eom = rowmean(dlog_CAD_eom dlog_DKK_eom dlog_EUR_eom dlog_GBP_eom dlog_JPY_eom dlog_NOK_eom dlog_SEK_eom);
qui gen log_e_eom = sum(dlog_e_eom);
qui egen dlog_q_eom = rowmean(dlog_q_CAD_eom dlog_q_DKK_eom dlog_q_EUR_eom dlog_q_GBP_eom dlog_q_JPY_eom dlog_q_NOK_eom dlog_q_SEK_eom);
qui gen log_q_eom = sum(dlog_q_eom);
qui egen cip_govt_3m_2020 = rowmean(cip_govt_CAD_3m_2020 cip_govt_DKK_3m_2020 cip_govt_EUR_3m_2020 cip_govt_GBP_3m_2020 cip_govt_JPY_3m_2020 cip_govt_NOK_3m_2020 cip_govt_SEK_3m_2020);
qui replace cip_govt_3m_2020 = cip_govt_3m_2020/100;
qui egen diffy_3m_2020 = rowmean(diffy_CAD_3m_2020 diffy_DKK_3m_2020 diffy_EUR_3m_2020 diffy_GBP_3m_2020 diffy_JPY_3m_2020 diffy_NOK_3m_2020 diffy_SEK_3m_2020);
foreach currency in CAD DKK EUR GBP JPY NOK SEK USD {;
	qui gen dlog_ip_oecd_`currency' = 100*(log(ip_oecd_`currency')-log(L1.ip_oecd_`currency'));
};	
qui egen dlog_ip_oecd = rowmean(dlog_ip_oecd_CAD dlog_ip_oecd_DKK dlog_ip_oecd_EUR dlog_ip_oecd_GBP dlog_ip_oecd_JPY dlog_ip_oecd_NOK dlog_ip_oecd_SEK);
qui gen log_ip_oecd = sum(dlog_ip_oecd);
qui gen log_ip_oecd_USD = sum(dlog_ip_oecd_USD);
qui gen diff_log_ip = log_ip_oecd - log_ip_oecd_USD;

qui keep year month excessequity excessequity3mo log_e_eom log_q_eom cip_govt_3m_2020 diffy_3m_2020 log_ip_oecd_USD diff_log_ip; 
export delimited "`direc'\output\monthly_var.csv", replace;

qui keep if month == 3 | month == 6 | month == 9 | month == 12;
export delimited "`direc'\output\quarterly_var.csv", replace;

********** Save data for Great Recession model vs. data comparison **********;

use "`direc'\data\Monthly_Data.dta", replace; 
qui egen cip_govt_3m_2020 = rowmean(cip_govt_AUD_3m_2020 cip_govt_CAD_3m_2020 cip_govt_CHF_3m_2020 
 cip_govt_DKK_3m_2020 cip_govt_EUR_3m_2020 cip_govt_GBP_3m_2020 cip_govt_JPY_3m_2020 cip_govt_NOK_3m_2020 cip_govt_NZD_3m_2020 cip_govt_SEK_3m_2020);
qui gen AUD_eom = 1/inv_AUD_eom;
qui gen EUR_eom = 1/inv_EUR_eom;
qui gen GBP_eom = 1/inv_GBP_eom;
qui gen NZD_eom = 1/inv_NZD_eom;
keep if month == 3 | month == 6 | month == 9 | month == 12;
qui gen quarter = month/3;
qui gen yq = yq(year,quarter);
format yq %tq;
drop month;
keep yq diffy* AUD_eom CAD_eom CHF_eom DKK_eom EUR_eom GBP_eom JPY_eom NOK_eom NZD_eom SEK_eom
 cpi_oecd_CAD cpi_oecd_CHF cpi_oecd_DKK cpi_oecd_EUR cpi_oecd_GBP cpi_oecd_JPY cpi_oecd_NOK cpi_oecd_SEK cpi_oecd_USD
 yield3mo_upd_eom p_all_bl cip_govt_3m_2020;
order yq;
qui merge 1:1 yq using "`direc'\data\Quarterly_Data.dta";
drop _merge;
tsset yq;
qui egen diffy_3m_2020 = rowmean(diffy_AUD_3m_2020 diffy_CAD_3m_2020 diffy_CHF_3m_2020
 diffy_DKK_3m_2020 diffy_EUR_3m_2020 diffy_GBP_3m_2020 diffy_JPY_3m_2020 diffy_NOK_3m_2020 diffy_NZD_3m_2020 diffy_SEK_3m_2020);
qui gen dlog_cpi_USD = log(cpi_oecd_USD) - log(L.cpi_oecd_USD);
foreach currency in AUD CAD CHF DKK EUR GBP JPY NOK NZD SEK {;
	qui gen dlog_`currency'_eom = log(`currency'_eom)-log(L.`currency'_eom);
	qui gen dlog_cpi_`currency' = log(cpi_oecd_`currency')-log(L.cpi_oecd_`currency');
	qui gen dlog_q_`currency'_eom = dlog_`currency'_eom + dlog_cpi_USD - dlog_cpi_`currency';
};
qui egen dlog_q_eom = rowmean(dlog_q_AUD_eom dlog_q_CAD_eom dlog_q_CHF_eom 
 dlog_q_DKK_eom dlog_q_EUR_eom dlog_q_GBP_eom dlog_q_JPY_eom dlog_q_NOK_eom dlog_q_NZD_eom dlog_q_SEK_eom);

qui gen gdpdef = us_gdp/us_rgdp;
qui gen y_pc = us_rgdp/us_pop;
qui gen dlog_y_pc = log(y_pc) - log(L.y_pc);
qui gen c_pc = (us_cn + us_cs)/(us_pop*gdpdef);
qui gen dlog_c_pc = log(c_pc) - log(L.c_pc);	
foreach currency in AUD CAD CHF DKK EUR GBP JPY NOK NZD SEK USD{;
	qui gen gdpdef_`currency' = gdp_`currency'/rgdp_`currency';
	qui gen y_pc_`currency' = rgdp_`currency'/pop_`currency';	
	qui gen dlog_y_pc_`currency' = log(y_pc_`currency') - log(L.y_pc_`currency');	
};
foreach currency in CAD DKK EUR GBP JPY NOK NZD SEK USD{;
	qui gen c_pc_`currency' = (cn_`currency' + cs_`currency')/(pop_`currency'*gdpdef_`currency');
	qui gen dlog_c_pc_`currency' = log(c_pc_`currency') - log(L.c_pc_`currency');
};
egen dlog_y_pc_avg = rowmean(dlog_y_pc_AUD dlog_y_pc_CAD dlog_y_pc_CHF 
 dlog_y_pc_DKK dlog_y_pc_EUR dlog_y_pc_GBP dlog_y_pc_JPY dlog_y_pc_NOK dlog_y_pc_NZD dlog_y_pc_SEK);	
qui egen dlog_c_pc_avg = rowmean(dlog_c_pc_CAD dlog_c_pc_DKK dlog_c_pc_EUR 
 dlog_c_pc_GBP dlog_c_pc_JPY dlog_c_pc_NOK dlog_c_pc_NZD dlog_c_pc_SEK);
gen yearquarter = year + (quarter-1)/4;
qui gen nfa_y = us_iip/(1000*us_gdp);
qui gen dnfa_y = (us_iip - L.us_iip)/(1000*us_gdp) if year >= 2006;
qui gen nx_y = (us_ex_trans - us_im_trans)/(1000*us_gdp);
tsset yq;
keep if year >= 1995 & year < .;
export delimited yearquarter year quarter p_all_bl cip_govt_3m_2020
 dlog_q_eom dlog_c_pc dlog_c_pc_avg dlog_y_pc dlog_y_pc_avg 
 nfa_y dnfa_y nx_y yield3mo_upd_eom diffy_3m_2020 using "`direc'\output\modelvsdata.csv", replace;

timer off 1;
timer list 1;