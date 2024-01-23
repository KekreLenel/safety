% prepare series 

get_var_indices; % script that collects variable indices 

% shock series
p_shocks          = temp_series(:,jx_p_shock);
z_shocks          = temp_series(:,jx_z_shock);
zf_shocks        = temp_series(:,jx_zf_shock);
omg_shocks        = temp_series(:,jx_omg_shock);

% cumulative growth (z in paper)
z_series          = exp(cumsum(z_shocks));

% only disaster series
dis_series        = zeros(size(z_shocks));
dis_series(z_shocks >= -trg_prms.disast_shock - sqrt(eps) & z_shocks <= -trg_prms.disast_shock + sqrt(eps) )  ...
    = z_shocks(z_shocks >= -trg_prms.disast_shock - sqrt(eps) & z_shocks <= -trg_prms.disast_shock + sqrt(eps) );
z_shocks_no_dis   = z_shocks-dis_series;
z_shocks_no_dis   = z_shocks-dis_series;

% capital stock
k_series          = temp_series(:,jx_k).*z_series;
k_series_unadj    = temp_series(:,jx_k);
k_chosen_series   = k_series./exp(dis_series); % capital chosen at t
k_chosen_series   = [k_chosen_series(2:end); k_chosen_series(end)];
E_kpa_series      = temp_series(:,jx_Ekpa).*z_series;
E_kpa_series_unadj= temp_series(:,jx_Ekpa);


% wealth shares
tht_h_series      = temp_series(:,jx_tht_h);

% real wage chosen last period
wh_state_series   = temp_series(:,jx_whs);
wf_state_series   = temp_series(:,jx_wfs);

% other state series
p_series          = exp(temp_series(:,jx_p));
zf_series        = temp_series(:,jx_zf);
rho_ih_series     = temp_series(:,jx_rho_ih);
rho_if_series     = temp_series(:,jx_rho_if);

omg_series        = exp(temp_series(:,jx_omg)) + trg_prms.omg_shift; 

% goods price: P/P_h and P*/P_h*
P_Phh_series      = temp_series(:,jx_P_P_hh);
P_Phf_series      = temp_series(:,jx_P_P_hf);

% labor series
lh_series         = temp_series(:,jx_lh);
lf_series         = temp_series(:,jx_lf);

% tot
s_series          = temp_series(:,jx_s);

% output in goods units
yh_series         = temp_series(:,jx_yh).*z_series;
yf_series         = temp_series(:,jx_yf).*z_series;

% expected return series
E_rfh_series        = temp_series(:,jx_Erfh);
E_rfh_f_series      = temp_series(:,jx_Erfh_f);
E_rff_series        = temp_series(:,jx_Erff);
E_rff_f_series      = temp_series(:,jx_Erff_f);
E_rk_f_series       = temp_series(:,jx_Erk_f);
E_rk_series         = temp_series(:,jx_Erk);
E_exc_series        = E_rk_series - E_rfh_series;
E_excrf_series      = E_rff_series - E_rfh_series;
E_rfh10_series      = temp_series(:,jx_Eq20hw);
E_rff10_series      = temp_series(:,jx_Eq20f);
E_rfh10_omg_series  = temp_series(:,jx_Eq20h);
E_rfh10_f_series    = temp_series(:,jx_fEq20hw);
E_rff10_f_series    = temp_series(:,jx_fEq20f);
E_rfh10_omg_f_series  = temp_series(:,jx_fEq20h);
E_rf_spread_series = temp_series(:,jx_Erf_spread);

corr_rf_spread_m1_series  = temp_series(:, jx_corr_rf_spread_mh);
corr_rf_spread_m2_series  = temp_series(:, jx_corr_rf_spread_mf);
corr_rf_spread_m1B_series  = temp_series(:, jx_corr_rf_spread_mhB);
corr_rf_spread_m2B_series  = temp_series(:, jx_corr_rf_spread_mfB);
corr_p_spread_m1B_series  = temp_series(:, jx_corr_p_mhB);
corr_p_spread_m2B_series  = temp_series(:, jx_corr_p_mfB);

% price of capital
q_series          = temp_series(:,jx_q);
q_f_series        = temp_series(:,jx_q).*P_Phh_series./P_Phf_series;

% inflation and  exchange rates
infl_h_series     = temp_series(:,jx_infl_h);
P_h_series        = cumprod(infl_h_series);
infl_f_series     = temp_series(:,jx_infl_f);
P_f_series        = cumprod(infl_f_series);
qx_series         = temp_series(:,jx_qx);
E_change_series   = -log(qx_series(1:end-1)./qx_series(2:end) .* infl_h_series(2:end)./infl_f_series(2:end));
E_change_series   = [E_change_series(1); E_change_series];
qx_change_series   = -log(qx_series(1:end-1)./qx_series(2:end));
qx_change_series   = [qx_change_series(1); qx_change_series];
E_series          = [1; exp(cumsum(E_change_series(2:end)))];
P_h_series        = [1; (cumprod(infl_h_series(2:end)))];
P_f_series        = [1; (cumprod(infl_f_series(2:end)))];

% wages (expressed in countries basket)
wh_series         = temp_series(:,jx_wh).*z_series;
wf_series         = temp_series(:,jx_wf).*z_series;

% investment (expressed in capital units)
inv_series        = temp_series(:,jx_inv).*z_series;
% spending on investment (expressed in each country's basket)
inv_h_series      = temp_series(:,jx_invh).*z_series;
inv_f_series      = temp_series(:,jx_invf).*z_series;

% aggregate consumption
c_h_series       = temp_series(:,jx_ch).*z_series;
c_f_series       = temp_series(:,jx_cf).*z_series;

% home consumption of home good, foreign good consumption in home good units
chh_series        = c_h_series./ ( (1.0-trg_prms.varsigma_h)/trg_prms.varsigma_h * s_series.^(trg_prms.sigma-1.0) + 1.0); 
chf_series        = c_h_series - chh_series; 

% foreign consumption of home good and foreign consumption in home good units
cfh_series        = c_f_series.*P_Phf_series./P_Phh_series./( (1.0-trg_prms.varsigma_f)/trg_prms.varsigma_f * s_series.^(trg_prms.sigma-1.0) + 1.0);
cff_series        = c_f_series.*P_Phf_series./P_Phh_series - cfh_series;

% portfolio shares
sh_h_series      = temp_series(:,jx_sh_h);
sh_f_series      = temp_series(:,jx_sh_f);
ssh_h_series      = temp_series(:,jx_ssh_h);
ssh_f_series      = temp_series(:,jx_ssh_f);


aggr_wealth_series = temp_series(:,jx_aggrw).*z_series;
bg_series          = trg_prms.bg_yss*(c_h_series + c_f_series./qx_series) +temp_series(:,jx_bg).*z_series;
bg_xtra_series     = temp_series(:,jx_bg).*z_series;
bg_rel_series      = bg_series./(c_h_series + c_f_series./qx_series);
bg_xtra_rel_series = bg_xtra_series./(c_h_series + c_f_series./qx_series);
bg_nom_series      = bg_series.*P_h_series;

% adjust wedge for govrnment debt issuance 
omg_series = omg_series - trg_prms.b_lmbd*bg_xtra_rel_series;


% inflation and nominal rates
nom_ih_series     = temp_series(:,jx_nom_i_h);
nom_if_series     = temp_series(:,jx_nom_i_f);
pi_series         = temp_series(:,jx_pi);

% savings series home in home basket units
h_sav_series      = temp_series(:,jx_h_sav).*z_series;
f_sav_series      = temp_series(:,jx_f_sav).*z_series;
h_ksav_series     = temp_series(:,jx_h_ksav).*z_series;
h_b_sav_series    = h_sav_series - h_ksav_series; 
h_bh_sav_series   = temp_series(:,jx_h_bhsav).*z_series;
h_bf_sav_series   = h_b_sav_series - h_bh_sav_series; 

% deployed capital series
h_kap_series      = temp_series(:,jx_h_kap).*z_series;
f_kap_series      = temp_series(:,jx_f_kap).*z_series;

% nominal bonds
q1_h_series       = temp_series(:,jx_q1h);
q2_h_series       = temp_series(:,jx_q2h);
q3_h_series       = temp_series(:,jx_q3h);
q4_h_series       = temp_series(:,jx_q4h);
q19_h_series      = temp_series(:,jx_q19h);
q20_h_series      = temp_series(:,jx_q20h);
q1_hw_series      = temp_series(:,jx_q1hw);
q2_hw_series      = temp_series(:,jx_q2hw);
q3_hw_series      = temp_series(:,jx_q3hw);
q4_hw_series      = temp_series(:,jx_q4hw);
q19_hw_series     = temp_series(:,jx_q19hw);
q20_hw_series     = temp_series(:,jx_q20hw);
q1_f_series       = temp_series(:,jx_q1f);
q2_f_series       = temp_series(:,jx_q2f);
q3_f_series       = temp_series(:,jx_q3f);
q4_f_series       = temp_series(:,jx_q4f);
q19_f_series      = temp_series(:,jx_q19f);
q20_f_series      = temp_series(:,jx_q20f);

rq1_h_series        = log(1./infl_h_series(2:end)./q1_h_series(1:end-1));
rq2_h_series        = log(q1_h_series(2:end)./infl_h_series(2:end)./q2_h_series(1:end-1));
rq3_h_series        = log(q2_h_series(2:end)./infl_h_series(2:end)./q3_h_series(1:end-1));
rq4_h_series        = log(q3_h_series(2:end)./infl_h_series(2:end)./q4_h_series(1:end-1));
rq20_h_series       = log(q19_h_series(2:end)./infl_h_series(2:end)./q20_h_series(1:end-1));
rq1_hw_series       = log(1./infl_h_series(2:end)./q1_hw_series(1:end-1));
rq2_hw_series       = log(q1_hw_series(2:end)./infl_h_series(2:end)./q2_hw_series(1:end-1));
rq3_hw_series       = log(q2_hw_series(2:end)./infl_h_series(2:end)./q3_hw_series(1:end-1));
rq4_hw_series       = log(q3_hw_series(2:end)./infl_h_series(2:end)./q4_hw_series(1:end-1));
rq20_hw_series      = log(q19_hw_series(2:end)./infl_h_series(2:end)./q20_hw_series(1:end-1));
rq1_f_series        = log(1./infl_f_series(2:end)./q1_f_series(1:end-1));
rq2_f_series        = log(q1_f_series(2:end)./infl_f_series(2:end)./q2_f_series(1:end-1));
rq3_f_series        = log(q2_f_series(2:end)./infl_f_series(2:end)./q3_f_series(1:end-1));
rq4_f_series        = log(q3_f_series(2:end)./infl_f_series(2:end)./q4_f_series(1:end-1));      
rq20_f_series       = log(q19_f_series(2:end)./infl_f_series(2:end)./q20_f_series(1:end-1));      
rq1_fh_series       = rq1_f_series  + log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1)); 
rq2_fh_series       = rq2_f_series  + log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq3_fh_series       = rq3_f_series  + log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq4_fh_series       = rq4_f_series  + log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq20_fh_series      = rq20_f_series + log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));

rq1_hf_series   = rq1_h_series   - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq2_hf_series   = rq2_h_series   - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq3_hf_series   = rq3_h_series   - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq4_hf_series   = rq4_h_series   - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq20_hf_series  = rq20_h_series  - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq1_hwf_series  = rq1_hw_series  - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq2_hwf_series  = rq2_hw_series  - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq3_hwf_series  = rq3_hw_series  - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq4_hwf_series  = rq4_hw_series  - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rq20_hwf_series = rq20_hw_series - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));

rq1_h_series  = [rq1_h_series(1); rq1_h_series ];
rq2_h_series  = [rq2_h_series(1); rq2_h_series ];
rq3_h_series  = [rq3_h_series(1); rq3_h_series ];
rq4_h_series  = [rq4_h_series(1); rq4_h_series ];
rq20_h_series  = [rq20_h_series(1); rq20_h_series ];
rq1_hw_series = [rq1_hw_series(1); rq1_hw_series];
rq2_hw_series = [rq2_hw_series(1); rq2_hw_series];
rq3_hw_series = [rq3_hw_series(1); rq3_hw_series];
rq4_hw_series = [rq4_hw_series(1); rq4_hw_series];
rq20_hw_series = [rq20_hw_series(1); rq20_hw_series];
rq1_hf_series  = [rq1_hf_series(1); rq1_hf_series ];
rq2_hf_series  = [rq2_hf_series(1); rq2_hf_series ];
rq3_hf_series  = [rq3_hf_series(1); rq3_hf_series ];
rq4_hf_series  = [rq4_hf_series(1); rq4_hf_series ];
rq20_hf_series  = [rq20_hf_series(1); rq20_hf_series ];
rq1_hwf_series = [rq1_hwf_series(1); rq1_hwf_series];
rq2_hwf_series = [rq2_hwf_series(1); rq2_hwf_series];
rq3_hwf_series = [rq3_hwf_series(1); rq3_hwf_series];
rq4_hwf_series = [rq4_hwf_series(1); rq4_hwf_series];
rq20_hwf_series = [rq20_hwf_series(1); rq20_hwf_series];
rq1_f_series  = [rq1_f_series(1); rq1_f_series ];
rq2_f_series  = [rq2_f_series(1); rq2_f_series ];
rq3_f_series  = [rq3_f_series(1); rq3_f_series ];
rq4_f_series  = [rq4_f_series(1); rq4_f_series ];
rq20_f_series  = [rq20_f_series(1); rq20_f_series ];
rq1_fh_series = [rq1_fh_series(1); rq1_fh_series];
rq2_fh_series = [rq2_fh_series(1); rq2_fh_series];
rq3_fh_series = [rq3_fh_series(1); rq3_fh_series];
rq4_fh_series = [rq4_fh_series(1); rq4_fh_series];
rq20_fh_series = [rq20_fh_series(1); rq20_fh_series];

% nominal bonds, who prices them
who_q1_h_series       = temp_series(:,jx_Pq1h);
who_q2_h_series       = temp_series(:,jx_Pq2h);
who_q3_h_series       = temp_series(:,jx_Pq3h);
who_q4_h_series       = temp_series(:,jx_Pq4h);
who_q1_hw_series       = temp_series(:,jx_Pq1h);
who_q2_hw_series       = temp_series(:,jx_Pq2h);
who_q3_hw_series       = temp_series(:,jx_Pq3h);
who_q4_hw_series       = temp_series(:,jx_Pq4h);
who_q1_f_series       = temp_series(:,jx_Pq1f);
who_q2_f_series       = temp_series(:,jx_Pq2f);
who_q3_f_series       = temp_series(:,jx_Pq3f);
who_q4_f_series       = temp_series(:,jx_Pq4f);

v1_series       = temp_series(:,jx_v1).*z_series;
v2_series       = temp_series(:,jx_v2).*z_series;

ce1_series       = temp_series(:,jx_ce1);
ce2_series       = temp_series(:,jx_ce2);

% derived return series
rfh_series        = log(nom_ih_series(1:end-1)./infl_h_series(2:end)) ;
rfh_series        = [rfh_series(1); rfh_series];
rfh_f_series      = log(nom_ih_series(1:end-1)./infl_h_series(2:end)) - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rfh_f_series        = [rfh_f_series(1); rfh_f_series];
rfh_omg_series    = log(nom_ih_series(1:end-1)./infl_h_series(2:end)./(1 - omg_series(1:end-1)));
rfh_omg_series    = [rfh_omg_series(1); rfh_omg_series];
rfh_f_omg_series   = log(nom_ih_series(1:end-1)./infl_h_series(2:end)./(1 - omg_series(1:end-1))) - log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rfh_f_omg_series    = [rfh_f_omg_series(1); rfh_f_omg_series];
rff_series        = log((nom_if_series(1:end-1))./infl_f_series(2:end)) ;
rff_series        = [rff_series(1); rff_series];
rff_h_series      = rff_series(2:end)+log(P_Phf_series(2:end)./P_Phf_series(1:end-1)./P_Phh_series(2:end).*P_Phh_series(1:end-1));
rff_h_series      = [rff_h_series(1); rff_h_series];
rk_series         = log(exp(dis_series(2:end)).*(pi_series(2:end) + (1.0-trg_prms.delta)*q_series(2:end))./q_series(1:end-1));
rk_series         = [rk_series(1); rk_series];
rk_f_series       = rk_series(2:end)+log(P_Phh_series(2:end)./P_Phh_series(1:end-1)./P_Phf_series(2:end).*P_Phf_series(1:end-1));
rk_f_series       = [rk_f_series(1); rk_f_series];

rfh_lev_series  = log(1/(1+trg_prms.zeta)*exp(rq20_h_series)  + trg_prms.zeta/(1+trg_prms.zeta)*exp(rq20_fh_series));
rff_lev_series  = (1/(1+trg_prms.zeta)*rq20_hf_series + trg_prms.zeta/(1+trg_prms.zeta)*rq20_f_series);
Erfh_lev_series = (1/(1+trg_prms.zeta)*E_rfh10_omg_series   + trg_prms.zeta/(1+trg_prms.zeta)*E_rff10_series);
Erff_lev_series = (1/(1+trg_prms.zeta)*E_rfh10_omg_f_series + trg_prms.zeta/(1+trg_prms.zeta)*E_rff10_f_series);

rA_series         = log((1 + trg_prms.debt_to_equity )*exp(rk_series) - trg_prms.debt_to_equity*exp(rfh_lev_series));
rA_f_series       = (1 + trg_prms.debt_to_equity )*rk_f_series - trg_prms.debt_to_equity*rff_lev_series;


exc_ret_series    = rk_series - rfh_series;
exc_ret_omg_series  = rk_series - rfh_omg_series;

exc_retA_series      = rA_series - rfh_series;
exc_retA_omg_series  = rA_series - rfh_omg_series;

exc_ret_f_series    = (rk_f_series) - (rff_series);
exc_retA_f_series   = rA_f_series - (rff_series);
exc_ret_fh_series   = rk_series - (rff_h_series);
exc_retA_fh_series  = rA_series - (rff_h_series);

exc_rf_series = rff_h_series - rfh_series;
exc_rf_omg_series = rff_h_series - rfh_omg_series;

h_wealth_series   = aggr_wealth_series.*tht_h_series;

% investor's returns

ret_h_series     = (1.0-sh_h_series(1:end-1)).*exp(rk_series(2:end)) + (ssh_h_series(2:end).*exp(rfh_series(2:end))./(1-omg_series(1:end-1)) + (1.0-ssh_h_series(2:end) -(1- sh_h_series(2:end))).*exp(rff_h_series(2:end)));
ret_f_series     = (1.0-sh_f_series(1:end-1)).*exp(rk_series(2:end)) + (ssh_f_series(2:end).*exp(rfh_series(2:end))./(1-omg_series(1:end-1)) + (1.0-ssh_f_series(2:end) -(1- sh_f_series(2:end))).*exp(rff_h_series(2:end)));

ret_h_series     = [ret_h_series(1); ret_h_series];
ret_f_series     = [ret_f_series(1); ret_f_series];

% nfa series
nfa_series       = h_sav_series(1:end-1) - h_kap_series(2:end).*exp(-dis_series(2:end)).*q_series(1:end-1);
nfa_series       = [nfa_series; nfa_series(end)];
nfa_k_series      = h_ksav_series(1:end-1)- h_kap_series(2:end).*exp(-dis_series(2:end)).*q_series(1:end-1);
nfa_k_series      = [nfa_k_series; nfa_k_series(end)];

nfa_rel_series        = nfa_series./yh_series;
nfa_rel_growth_series = (nfa_series(2:end)-nfa_series(1:end-1))./yh_series(2:end);
nfa_rel_growth_series = [nfa_rel_growth_series(1); nfa_rel_growth_series];

nfa_relk_series    = nfa_series./h_sav_series;
nfa_k_rel_series  = nfa_k_series./h_sav_series;
nfa_b_rel_series  = nfa_series./h_sav_series - nfa_k_series./h_sav_series;
nfa_bh_rel_series  = nfa_series./h_sav_series - nfa_k_series./h_sav_series - h_bf_sav_series./h_sav_series;
nfa_bf_rel_series  = nfa_series./h_sav_series - nfa_k_series./h_sav_series - h_bh_sav_series./h_sav_series;
nfa_bh_rely_series  = nfa_series./yh_series - nfa_k_series./yh_series - h_bf_sav_series./yh_series;
nfa_bf_rely_series  = nfa_series./yh_series - nfa_k_series./yh_series - h_bh_sav_series./yh_series;
nfa_b_rely_series = (nfa_series - nfa_k_series)./yh_series;
nfa_k_rely_series = (nfa_k_series)./yh_series;


nfa_return_series = (mean(nfa_k_rel_series(1:end-1)).*rk_series(2:end) + mean(nfa_bh_rel_series(1:end-1)).*rfh_series(2:end)./(1-omg_series(1:end-1)) + mean(nfa_bf_rel_series(1:end-1)).*rff_h_series(2:end))./abs(mean(nfa_relk_series(1:end-1))); 

nfa_return_series = [nfa_return_series(1); nfa_return_series];
exc_retnfa_series = nfa_return_series - rfh_series;


E_rA_series        = (1 + trg_prms.debt_to_equity )*E_rk_series - trg_prms.debt_to_equity*Erfh_lev_series;
E_excA_series      = (log(E_rA_series) - log(E_rfh_series));
% E_excA_omg_series  = (E_rA_series - E_rfh_omg_series);
E_excA_f_series    = (1 + trg_prms.debt_to_equity)*E_rk_f_series - trg_prms.debt_to_equity*Erff_lev_series - E_rff_f_series; 

% unlevered claim on equity
divk_series_1       = pi_series(1:end)  - trg_prms.delta*q_series(1:end) ;
divk_price_series_1 = divk_series_1./q_series; 
divk_series_1f      = divk_series_1.*P_Phh_series./P_Phf_series ;

% unlevered claim on equity
divk_series_2  = k_series(1:end-1).*(pi_series(1:end-1) + (1 - trg_prms.delta)*q_series(1:end-1)) ... 
               - q_series(1:end-1).*k_series(2:end);
divk_series_2       = [divk_series_2; divk_series_2(end)]; 
divk_price_series_2 = divk_series_2./(q_series.*k_chosen_series); 
divk_series_2f      = divk_series_2.*P_Phh_series./P_Phf_series ;

% dividend series: profits of one unit of capital held every period, financed with fixed leverage

div_series_1 = pi_series(2:end) + (1 - trg_prms.delta)*q_series(2:end) - ...
                trg_prms.debt_to_equity/(1.0+trg_prms.debt_to_equity).*q_series(1:end-1).*(exp(rfh_lev_series(2:end)))   ...
               - 1/(1 + trg_prms.debt_to_equity).*q_series(2:end);
div_series_1f = pi_series(2:end).*P_Phh_series(2:end)./P_Phf_series(2:end)+ (1 - trg_prms.delta)*q_series(2:end).*P_Phh_series(2:end)./P_Phf_series(2:end) - ...
                trg_prms.debt_to_equity/(1.0+trg_prms.debt_to_equity).*q_series(1:end-1).*P_Phh_series(1:end-1)./P_Phf_series(1:end-1).*(exp(rff_lev_series(2:end))) ...
               - 1/(1 + trg_prms.debt_to_equity).*q_series(2:end).*P_Phh_series(2:end)./P_Phf_series(2:end);            

div_series_1 = [div_series_1(1); div_series_1];

div_series_1f = [div_series_1f(1); div_series_1f];

qE_series_1   = 1/(1 + trg_prms.debt_to_equity).*q_series(1:end);

% dividend series : aggregate profits, financed with fixed leverage, aggregate investment
div_series_2= k_series(2:end-1).*(pi_series(2:end-1) + (1 - trg_prms.delta)*q_series(2:end-1) ... 
            - trg_prms.debt_to_equity/(1.0+trg_prms.debt_to_equity).*q_series(1:end-2).*(exp(rfh_lev_series(2:end-1))) ) ...
               - 1/(1 + trg_prms.debt_to_equity).*q_series(2:end-1).*k_series(3:end);
div_series_2f = k_series(2:end-1).*(pi_series(2:end-1).*P_Phh_series(2:end-1)./P_Phf_series(2:end-1) + (1 - trg_prms.delta)*q_series(2:end-1).*P_Phh_series(2:end-1)./P_Phf_series(2:end-1) ... 
               - trg_prms.debt_to_equity/(1.0+trg_prms.debt_to_equity).*q_series(1:end-2).*P_Phh_series(1:end-2)./P_Phf_series(1:end-2).*(exp(rff_lev_series(2:end-1))) ) ...
               - 1/(1 + trg_prms.debt_to_equity).*q_series(2:end-1).*P_Phh_series(2:end-1)./P_Phf_series(2:end-1).*k_series(3:end);
div_series_2 = [div_series_2(1); div_series_2; div_series_2(end)];
div_series_2f = [div_series_2f(1); div_series_2f; div_series_2f(end)];


qE_series_2 = 1/(1 + trg_prms.debt_to_equity).*q_series.*k_chosen_series;

div_price_series_1 = div_series_1./qE_series_1; 
div_price_series_2 = div_series_2./qE_series_2; 


div_price_smoothed_series_1 = movsum(div_series_1,[3,0])./qE_series_1;
div_price_smoothed_series_1(1:3) = div_price_smoothed_series_1(1)*4;
div_price_smoothed_series_2 = movsum(div_series_2,[3,0])./qE_series_2;
div_price_smoothed_series_2(1:3) = div_price_smoothed_series_2(1)*4;

yield1_h_series = -log(q1_h_series);
yield2_h_series = -log(q2_h_series)/2;
yield3_h_series = -log(q3_h_series)/3;
yield4_h_series = -log(q4_h_series)/4;
yield1_hw_series = -log(q1_hw_series);
yield2_hw_series = -log(q2_hw_series)/2;
yield3_hw_series = -log(q3_hw_series)/3;
yield4_hw_series = -log(q4_hw_series)/4;
yield1_f_series = -log(q1_f_series);
yield2_f_series = -log(q2_f_series)/2;
yield3_f_series = -log(q3_f_series)/3;
yield4_f_series = -log(q4_f_series)/4;

kh_series = h_ksav_series./q_series; % kh saved in t

imports =  chf_series./s_series./P_Phh_series;
exports =  cfh_series(1:end-1)./P_Phh_series(1:end-1) + inv_h_series(1:end-1)./P_Phh_series(1:end-1)  + q_series(1:end-1).*(1.0-trg_prms.delta).*h_kap_series(1:end-1) - q_series(1:end-1).*exp(-dis_series(2:end)).*h_kap_series(2:end);
exports = [exports; exports(end)];
nx_series = exports - imports; 

nx_rely_series = nx_series./yh_series; 


collected_data = [ h_bh_sav_series(2:end-1), h_bf_sav_series(2:end-1), h_ksav_series(2:end-1), h_kap_series(2:end-1), ...
                    exp(rfh_series(2:end-1)),                exp(rff_series(2:end-1)),         qx_series(2:end-1),     exp(rk_series(2:end-1)), ...
                     q_series(2:end-1),                      z_series(2:end-1),                z_shocks_no_dis(2:end-1),        dis_series(2:end-1),    exports(2:end-1), ...
                     imports(2:end-1),                       tht_h_series(2:end-1),            k_series(3:end)./exp(dis_series(3:end)), omg_series(2:end-1)./(1-omg_series(2:end-1))];


val_series =  h_bf_sav_series(1:end-1).*(exp(rff_h_series(2:end))-exp(rfh_omg_series(2:end))) ... 
              +  ( h_ksav_series(1:end-1) - q_series(1:end-1).*h_kap_series(2:end) ).*(exp(rk_series(2:end))-exp(rfh_omg_series(2:end)));
val_series = [val_series(1); val_series]; % h_asset_income_series - return on kappa assets

valE_series =  h_bf_sav_series(1:end-1).*(exp(rff_h_series(2:end))-exp(rfh_omg_series(2:end))) ... 
              +  ( h_ksav_series(1:end-1) - q_series(1:end-1).*E_kpa_series(1:end-1) ).*(exp(rk_series(2:end))-exp(rfh_omg_series(2:end)));
valE_series = [valE_series(1); valE_series]; % h_asset_income_series - return on kappa assets

val_rff_series = h_bf_sav_series(1:end-1).*(exp(rff_h_series(2:end))-exp(rfh_omg_series(2:end)));
val_rff_series = [val_rff_series(1); val_rff_series];

val_rk_series = ( h_ksav_series(1:end-1) - q_series(1:end-1).*h_kap_series(2:end) ).*(exp(rk_series(2:end))-exp(rfh_omg_series(2:end)));
val_rk_series = [val_rk_series(1); val_rk_series];

valE_rk_series = ( h_ksav_series(1:end-1) - q_series(1:end-1).*E_kpa_series(1:end-1) ).*(exp(rk_series(2:end))-exp(rfh_omg_series(2:end)));
valE_rk_series = [valE_rk_series(1); valE_rk_series];

valb_series =    (h_bh_sav_series(1:end-1)+ h_bf_sav_series(1:end-1)).*(exp(rfh_omg_series(2:end))-exp(rk_series(2:end))) ... 
              +  h_bf_sav_series(1:end-1).*(exp(rff_h_series(2:end))-exp(rfh_omg_series(2:end)));
valb_series = [valb_series(1); valb_series]; % h_asset_income_series - return on kappa assets

valb_rf_series = (h_bh_sav_series(1:end-1)+ h_bf_sav_series(1:end-1)).*(exp(rfh_omg_series(2:end))-exp(rk_series(2:end)));
valb_rf_series = [valb_rf_series(1); valb_rf_series];

valb_rff_series =  h_bf_sav_series(1:end-1).*(exp(rff_h_series(2:end))-exp(rfh_omg_series(2:end)));
valb_rff_series = [valb_rff_series(1); valb_rff_series];

c_grth_h_series = c_h_series(2:end)./c_h_series(1:end-1);
c_grth_f_series = c_f_series(2:end)./c_f_series(1:end-1);
c_grth_h_series = [c_grth_h_series(1); c_grth_h_series];
c_grth_f_series = [c_grth_f_series(1); c_grth_f_series];
        
c_h_change_series = log(c_h_series(2:end)) - log(c_h_series(1:end-1));
c_h_change_series = [0; c_h_change_series];
c_f_change_series = log(c_f_series(2:end)) - log(c_f_series(1:end-1));
c_f_change_series = [0; c_f_change_series];
inv_change_series = log(inv_series(2:end)) - log(inv_series(1:end-1));
inv_change_series = [inv_change_series(1); inv_change_series];
y_h_change_series = log(yh_series(2:end)) - log(yh_series(1:end-1));
y_h_change_series = [0; y_h_change_series];
y_f_change_series = log(yf_series(2:end)) - log(yf_series(1:end-1));
y_f_change_series = [0; y_f_change_series];

% regressions using private bonds
y_growth_series = 100*(log(yh_series(5:end)) - log(yh_series(1:end-4))); 
y_growth_series = [zeros(4,1);y_growth_series];
uip_pvt_series  = 100*(yield1_f_series(1:end-1) - yield1_hw_series(1:end-1) - E_change_series(2:end));
uip_pvt_series  = [0.0; uip_pvt_series];

mc_h_series = (( 1 + (1/trg_prms.ies_h - 1)*trg_prms.chi0_h * trg_prms.chi/(1+trg_prms.chi)*(lh_series.^( (1.0 + trg_prms.chi)/trg_prms.chi )))./c_h_series).^(1/trg_prms.ies_h);
mc_f_series = (( 1 + (1/trg_prms.ies_f - 1)*trg_prms.chi0_f * trg_prms.chi/(1+trg_prms.chi)*(lf_series.^( (1.0 + trg_prms.chi)/trg_prms.chi )))./c_f_series).^(1/trg_prms.ies_f);
sdf_h_series = trg_prms.bbeta_h*mc_h_series(2:end)./mc_h_series(1:end-1).*(v1_series(2:end)./ce1_series(1:end-1)).^(1/trg_prms.ies_h - trg_prms.gma_h);
sdf_f_series = trg_prms.bbeta_f*mc_f_series(2:end)./mc_f_series(1:end-1).*(v2_series(2:end)./ce2_series(1:end-1)).^(1/trg_prms.ies_f - trg_prms.gma_f);
sdf_h_series =[sdf_h_series(1);sdf_h_series];
sdf_f_series =[sdf_f_series(1);sdf_f_series];
