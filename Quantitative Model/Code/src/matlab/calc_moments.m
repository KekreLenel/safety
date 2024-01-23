% -------------------------------------------------------------------------
% calc_moments.m: calculate moments for tables 
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/safety
% -------------------------------------------------------------------------

%% Table 1 

    % see collect_moments

%% Table 2: targeted moments 
    if ~dis
        
        % regression: regress excess returns on price dividend ratio
        y_vec   = 4*100*exc_retA_series(5:end); % start at t=4 bc of dividend smoothing
        x_vec   = div_price_smoothed_series_1(4:end-1);
        [~,~,resid_dp,~,~] = regress(y_vec,[ ones(size(x_vec)), x_vec]);
        
        % regression: excess foreign bond returns on yield difference and past output growth
        fama_yield_pvt_series  = 100*4*(yield1_f_series - yield1_h_series);
        x_vec  = uip_pvt_series(5:end);        
        y_vec  = fama_yield_pvt_series(4:end-1); 
        y_vec2 = y_growth_series(4:end-1);       
        [~,~,resid_carry,~,~] = regress(x_vec,[ ones(size(y_vec)), y_vec, y_vec2]);

        % regression: nfa growth on excess equity returns and excess foreign bond excess returns
        x_vec  = (exc_retA_series(3:end-1));
        x2_vec = (exc_rf_series(3:end-1));
        y_vec = (( nfa_rel_growth_series(3:end-1) ));
        reg_coeffs1 = regress(y_vec,[ ones(size(x_vec)), x_vec, x2_vec]);

        table_2_moms_tmp(1,sss)  = mean(yf_series./(s_series.*yh_series));
        table_2_moms_tmp(2,sss)  = 100*std(diff(log(c_h_series)));
        table_2_moms_tmp(3,sss)  = 100*std(diff(log(yf_series)));   
        table_2_moms_tmp(4,sss)  = corr(yf_series(1:end-1)./yh_series(1:end-1),yf_series(2:end)./yh_series(2:end));
        table_2_moms_tmp(5,sss)  = 100*std(diff(log(inv_series)));
        table_2_moms_tmp(6,sss)  = 4*100*mean(rfh_series);
        table_2_moms_tmp(7,sss)  = 100*mean(nfa_rel_series/4);
        table_2_moms_tmp(8,sss)  = corr(resid_dp, resid_carry);
        table_2_moms_tmp(9,sss)  = 4*100*mean((rA_series(2:end-1)) - mean((rfh_series(2:end-1))));
        table_2_moms_tmp(10,sss) = reg_coeffs1(2);
        table_2_moms_tmp(11,sss) = 100*trg_prms.bg_yss*mean((c_f_series./qx_series)./(4* yh_series));
        table_2_moms_tmp(12,sss) = mean(lh_series);
        table_2_moms_tmp(13,sss) = mean(lf_series);        
        table_2_moms_tmp(14,sss) = 100*mean(log(infl_h_series));
        table_2_moms_tmp(15,sss) = 100*mean(log(infl_f_series));    

    end

%% Table 3: comovements 
    if ~dis

       y_vec  = uip_pvt_series(5+1:end-1);        
       x_vec  = y_growth_series(5:end-2); 
       reg_coeffs2 = regress(y_vec,[ ones(size(x_vec)), x_vec]);
       
       table_3_moms_tmp(1,sss) = reg_coeffs2(2);
       
       y_vec  = uip_pvt_series(3:end-1);        
       x_vec  = 100*exc_retA_series(3:end-1); 
       reg_coeffs3 = regress(y_vec,[ ones(size(x_vec)), x_vec]);
       
       table_3_moms_tmp(2,sss) = reg_coeffs3(2);
       
       % from table 2 regression
       table_3_moms_tmp(3,sss) = reg_coeffs1(3);
              
       table_3_moms_tmp(4,sss) = 100*mean((h_ksav_series - h_kap_series)./(4*yh_series));
       table_3_moms_tmp(5,sss) = 100*mean((h_bh_sav_series)./(4*yh_series));
       table_3_moms_tmp(6,sss) = 100*mean((h_bf_sav_series)./(4*yh_series));

    end


%% Table 4: additional second moments 
    if ~dis

       table_4_moms_tmp(1,sss) = 4*100*std(rfh_series(2:end-1));
       table_4_moms_tmp(2,sss) = 4*100*std(exc_retA_series(2:end-1));
       table_4_moms_tmp(3,sss) = 4*100*std(rff_h_series(2:end-1) - rfh_series(2:end-1));
       table_4_moms_tmp(4,sss) = 100*std(qx_change_series(1:end-1));
       table_4_moms_tmp(5,sss) = 100*std(E_change_series(2:end-1));
       table_4_moms_tmp(6,sss) = corr( log(qx_series(3:end-1)./qx_series(2:end-2)), log(c_f_series(3:end-1)./c_f_series(2:end-2)) - log(c_h_series(3:end-1)./c_h_series(2:end-2)));
    
    end

%% Table 5: output volatility 
    
    if ~dis 

       table_5_moms_tmp(1,sss) = 100*std(log(yh_series(3:end-1)./  yh_series(2:end-2)));
       table_5_moms_tmp(2,sss) = 100*std(log(yf_series(3:end-1)./  yf_series(2:end-2)));

    end

%% Table 6: net foreign assets 

    if ~dis 
        
       table_6_moms_tmp(1,sss) = 100*std(nfa_rel_growth_series(2:end-1));
       table_6_moms_tmp(2,sss) = 100*std(nx_rely_series(2:end-1));
       table_6_moms_tmp(3,sss) = 100*std(nfa_rel_growth_series(2:end-1) - nx_rely_series(2:end-1));
       table_6_moms_tmp(4,sss) = 100*mean(nfa_rel_growth_series(2:end-1));
       table_6_moms_tmp(5,sss) = 100*mean(nx_rely_series(2:end-1));
       table_6_moms_tmp(6,sss) = 100*mean(nfa_rel_growth_series(2:end-1) - nx_rely_series(2:end-1));

    end

%% Table 7: external adjustment 

    if dis 
        
        nfa_innovations       = temp_series(2:end,jx_val_nfa)     - temp_series(1:end-1,jx_val_Enfa);
        val_innovations       = temp_series(2:end,jx_val_Eval)    - temp_series(1:end-1,jx_val_EEval);
        val_rf_innovations    = temp_series(2:end,jx_val_Evalrf)  - temp_series(1:end-1,jx_val_EEvalrf);
        val_rff_innovations   = temp_series(2:end,jx_val_Evalrff) - temp_series(1:end-1,jx_val_EEvalrff);
        val_seign_innovations = temp_series(2:end,jx_val_Evalseign) - temp_series(1:end-1,jx_val_EEvalseign);
        nx_innovations        = temp_series(2:end,jx_val_Enx)     - temp_series(1:end-1,jx_val_EEnx);
        nfa_resid_innovations = temp_series(2:end,jx_val_EnfaEnd) - temp_series(1:end-1,jx_val_EEnfaEnd);
        
        varnfa(sss) = var(nfa_innovations);
        tmp = cov(-nx_innovations,nfa_innovations);         cov1  = tmp(1,2);
        tmp = cov(-val_innovations,nfa_innovations);        cov2  = tmp(1,2);
        tmp = cov(-val_rf_innovations,nfa_innovations);     cov2A = tmp(1,2);
        tmp = cov(-val_rff_innovations,nfa_innovations);    cov2B = tmp(1,2);
        tmp = cov(-val_seign_innovations,nfa_innovations);  cov2C = tmp(1,2);
        tmp = cov(nfa_resid_innovations,nfa_innovations);   cov3  = tmp(1,2);
        
        % see collect_moments for calculation of table moments 
        table_7_moms_tmp(1,sss) =  cov1;
        table_7_moms_tmp(2,sss) =  cov2;
        table_7_moms_tmp(3,sss) =  cov2A;
        table_7_moms_tmp(4,sss) =  cov2B;
        table_7_moms_tmp(5,sss) =  cov2C;
        table_7_moms_tmp(6,sss) =  cov3;

    end

%% Table 8: swap line effects 
    
    % see collect_swap_moments.m

%% Table 9: comovements under alternative policy rules 

    % same moments as table 3 

%% Table 10: portfolios and risk premia

        table_10_moms_tmp(1,sss) = 100*mean((h_ksav_series )./(h_sav_series));
        table_10_moms_tmp(2,sss) = 100*mean((h_bh_sav_series )./(h_sav_series));
        table_10_moms_tmp(3,sss) = 100*mean((h_bf_sav_series )./(h_sav_series));
        table_10_moms_tmp(4,sss) = mean(corr_rf_spread_m1B_series);
        table_10_moms_tmp(5,sss) = mean(corr_rf_spread_m2B_series);



        
       
