% -------------------------------------------------------------------------
% extract_irfs.m: prepares impulse response series
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/safety
% -------------------------------------------------------------------------
function [irf_data, irf_idxs, irf_titles] = extract_irfs(temp_series, trg_prms)


global current_lmbd
n_irf = size(temp_series,1);
extract_series;

% E_series needs to be adjusted because only changes not levels are stored
% relative to no shock impulse respone
E_series = 1.0 + E_series - exp(E_change_series(2)*(0:n_irf-1)');

irf_data = zeros(198,14);
irf_data(:,1)   = 10000*(p_series(3:200,1)-p_series(1,1));                      irf_idxs.p          = 1 ;  irf_titles{1 }  = "$p$";
irf_data(:,2)   = 10000*(exc_retA_series(3:200,1)-exc_retA_series(1,1));        irf_idxs.excA       = 2 ;  irf_titles{2 }  = "$$r^e - r$$";          
irf_data(:,3)   = 10000*(exc_rf_series(3:200,1)-exc_rf_series(1,1));            irf_idxs.excrf      = 3 ;  irf_titles{3 }  = "$r^\ast - \Delta \log q - r$";
irf_data(:,4)   = 10000*(exc_ret_series(3:200,1)-exc_ret_series(1,1));          irf_idxs.excrk      = 4 ;  irf_titles{4 }  = "$r^k - r$";
irf_data(:,5)   = 10000*(E_rfh_series(3:200,1)-E_rfh_series(1,1));              irf_idxs.Erh        = 5 ;  irf_titles{5 }  = "$E[r]$";   
irf_data(:,6)   = 10000*(E_rff_f_series(3:200,1)-E_rff_f_series(1,1));          irf_idxs.Erf        = 6 ;  irf_titles{6 }  = "$E[r^\ast]$";   
irf_data(:,7)   = 10000*(log(qx_series(3:200,1))-log(qx_series(1,1)));          irf_idxs.log_qx     = 7 ;  irf_titles{7 }  = "$\log q$";
irf_data(:,8)   = 10000*(log(E_series(3:200,1))-log(E_series(1,1)));            irf_idxs.log_E      = 8 ;  irf_titles{8 }  = "$\log E $"; 
irf_data(:,9)   = 10000*(tht_h_series(3:200,1)-tht_h_series(1,1));              irf_idxs.thth       = 9 ;  irf_titles{9 }  = "$\theta$";  
irf_data(:,10)  = 10000*(log(yh_series(3:200,1))-log(yh_series(1,1)));          irf_idxs.log_yh     = 10;  irf_titles{10}  = "$\log y $";
irf_data(:,11)  = 10000*(log(yf_series(3:200,1))-log(yf_series(1,1)));          irf_idxs.log_yf     = 11;  irf_titles{11}  = "$\log y^\ast $";
irf_data(:,12)  = 10000*(log(c_h_series(3:200,1))-log(c_h_series(1,1)));        irf_idxs.log_ch     = 12;  irf_titles{12}  = "$\log c $";
irf_data(:,13)  = 10000*(log(c_f_series(3:200,1))-log(c_f_series(1,1)));        irf_idxs.log_cf     = 13;  irf_titles{13}  = "$\log c^\ast $";
irf_data(:,14)  = 10000*(log(lh_series(3:200,1))-log(lh_series(1,1)));          irf_idxs.log_lh     = 14;  irf_titles{14}  = "$\log \ell $";
irf_data(:,15)  = 10000*(log(lf_series(3:200,1))-log(lf_series(1,1)));          irf_idxs.log_lf     = 15;  irf_titles{15}  = "$\log \ell^\ast $";
irf_data(:,16)  = 10000*(log(wh_series(3:200,1))-log(wh_series(1,1)));          irf_idxs.log_wh     = 16;  irf_titles{16}  = "$\log w $";
irf_data(:,17)  = 10000*(log(wf_series(3:200,1))-log(wf_series(1,1)));          irf_idxs.log_wf     = 17;  irf_titles{17}  = "$\log w^\ast $";
irf_data(:,18)  = 10000*(log(infl_h_series(3:200,1))-log(infl_h_series(1,1)));  irf_idxs.log_infl_h = 18;  irf_titles{18}  = "$\log P/P_{-1} $";
irf_data(:,19)  = 10000*(log(infl_f_series(3:200,1))-log(infl_f_series(1,1)));  irf_idxs.log_infl_f = 19;  irf_titles{19}  = "$\log P^\ast/P^\ast_{-1} $";
irf_data(:,20)  = 10000*(log(inv_series(3:200,1))-log(inv_series(1,1)));        irf_idxs.log_x      = 20;  irf_titles{20}  = "$\log x $";
irf_data(:,21)  = 10000*(omg_series(3:200,1)-omg_series(1,1));                  irf_idxs.omg        = 21;  irf_titles{21}  = "$\omega $"; 
irf_data(:,22)  = 10000*(log(q_series(3:200,1))-log(q_series(1,1)));            irf_idxs.log_qk     = 22;  irf_titles{22}  = "$\log q^k $"; 
irf_data(:,23)  = 10000*(log(h_kap_series(3:200,1))-log(h_kap_series(1,1)));    irf_idxs.log_kh     = 23;  irf_titles{23}  = "$\log \kappa $"; 
irf_data(:,24)  = 10000*(nfa_series(3:200,1)-nfa_series(1,1))./yh_series(1);    irf_idxs.nfa_ybar   = 24;  irf_titles{24}  = "$nfa/y_0$"; 
irf_data(:,25)  = 10000*(valb_series(3:200,1)-valb_series(1,1))./yh_series(1);  irf_idxs.val_ybar   = 25;  irf_titles{25}  = "$val/y_0$"; 
irf_data(:,26)  = 10000*(nx_series(3:200,1)-nx_series(1,1))./yh_series(1);      irf_idxs.nx_ybar    = 26;  irf_titles{26}  = "$nx/y_0$"; 
irf_data(:,27)  = 10000*(log(z_series(3:200,1))-log(z_series(1,1)));            irf_idxs.log_z      = 27;  irf_titles{27}  = "$\log z $"; 
irf_data(:,28)  = 10000*((zf_series(3:200,1))-(zf_series(1,1)));                irf_idxs.log_zf     = 28;  irf_titles{28}  = "$\log z_f $";
irf_data(:,29)  = 10000*(nom_ih_series(3:200,1)-nom_ih_series(1,1));            irf_idxs.ih         = 29;  irf_titles{29}  = "$i $";
irf_data(:,30)  = 10000*(nom_if_series(3:200,1)-nom_if_series(1,1));            irf_idxs.if         = 30;  irf_titles{30}  = "$i^\ast $";
irf_data(:,31)  = 10000*(bg_rel_series(3:200,1)-bg_rel_series(1,1));            irf_idxs.bg         = 31;  irf_titles{31}  = "$-b_{H,s}^g/(c + \zeta^*q^{-1}c^*)$";
irf_data(:,32)  = 10000*[0.1*trg_prms.omg_rho.^[0:max(size(p_series(3:end))-1)]];   irf_idxs.matteo   = 32;  irf_titles{32}  = "$\varsigma$";
irf_data(:,33)  = 10000*[0.1*trg_prms.omg_rho.^[0:max(size(p_series(3:end))-1)]];   irf_idxs.pog      = 33;  irf_titles{33}  = "$\log\gamma^\ast$";
irf_data(:,34)  = 10000*(E_excA_series(3:200,1)-E_excA_series(1,1));            irf_idxs.EexcA      = 34 ;     irf_titles{34}  = "$E[r^e-r]$";   
irf_data(:,35)  = 10000*(rA_series(3:200,1)-rA_series(1,1));                    irf_idxs.retA       = 35 ;  irf_titles{35}  = "$$r^e$$"; 
irf_data(:,36)  = 10000*(nfa_rel_series(3:200,1)-nfa_rel_series(1,1));          irf_idxs.nfa_y      = 36;  irf_titles{36}  = "$nfa/y$"; 
irf_data(:,37)  = 10000*(nx_rely_series(3:200,1)-nx_rely_series(1,1));          irf_idxs.nx_y       = 37;  irf_titles{37}  = "$nx/y$"; 
end
