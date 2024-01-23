% indices for series  

jx_z_shock   = 1;
jx_zf_shock = 2;
jx_p_shock   = 3;
jx_omg_shock = 4;

jx_k       = 5;
jx_tht_h   = 6;
jx_zf     = 7;
jx_whs     = 8;
jx_wfs     = 9;
jx_p       = 10;
jx_rho_ih  = 11;
jx_rho_if  = 12;
jx_omg     = 13;

jx_yh      = 14;
jx_yf      = 15;
jx_lh      = 16;
jx_lf      = 17;
jx_ch      = 18;
jx_cf      = 19;
jx_inv     = 20;
jx_pi      = 21;
jx_pif     = 22;
jx_P_P_hh  = 23;
jx_P_P_hf  = 24;
jx_qx      = 25;
jx_s       = 26;
jx_Erfh    = 27;
jx_Erff    = 28;
jx_Erk     = 29;
jx_q       = 30;
jx_sh_h    = 31;
jx_sh_f    = 32;
jx_ssh_h   = 33;
jx_ssh_f   = 34;

jx_aggrw   = 35;
jx_wh      = 36;
jx_wf      = 37;
jx_infl_h  = 38;
jx_infl_f  = 39;
jx_nom_i_h = 40;
jx_nom_i_f = 41;
jx_h_ksav  = 42;
jx_h_kap   = 43;
jx_f_kap   = 44;
jx_Erfh_f  = 45;
jx_Erff_f  = 46;
jx_Erk_f   = 47;
jx_Ekpa    = 48;
jx_v1      = 49;
jx_v2      = 50;

jx_invh    = 51;
jx_invf    = 52;
jx_h_bhsav = 53;
jx_bg      = 54;
jx_ce1     = 55;
jx_ce2     = 56;
jx_h_sav   = 57;
jx_f_sav   = 58;
jx_Erf_spread = 59;
jx_corr_rf_spread_mh  = 60;
jx_corr_rf_spread_mf  = 61;
jx_corr_rf_spread_mhB = 62;
jx_corr_rf_spread_mfB = 63;
jx_corr_p_mhB = 64;
jx_corr_p_mfB = 65;

jx_q1h      = 66 ;
jx_q1f      = 67 ;
jx_q1hw     = 68 ;
jx_q2h      = 69 ;
jx_q2f      = 70;
jx_q2hw     = 71;
jx_q3h      = 72;
jx_q3f      = 73;
jx_q3hw     = 74;
jx_q4h      = 75;
jx_q4f      = 76;
jx_q4hw     = 77;
jx_q19h     = 78;
jx_q19f     = 79;
jx_q19hw    = 80;
jx_q20h     = 81;
jx_q20f     = 82;
jx_q20hw    = 83;

jx_Eq1h     = 84;
jx_Eq1f     = 85;
jx_Eq1hw    = 86;
jx_Eq2h     = 87;
jx_Eq2f     = 88;
jx_Eq2hw    = 89;
jx_Eq3h     = 90;
jx_Eq3f     = 91;
jx_Eq3hw    = 92;
jx_Eq4h     = 93;
jx_Eq4f     = 94;
jx_Eq4hw    = 95;
jx_Eq19h    = 96;
jx_Eq19f    = 97;
jx_Eq19hw   = 98;
jx_Eq20h    = 99;
jx_Eq20f    = 100;
jx_Eq20hw   = 101;

jx_fEq1h    = 102;
jx_fEq1f    = 103;
jx_fEq1hw   = 104;
jx_fEq2h    = 105;
jx_fEq2f    = 106;
jx_fEq2hw   = 107;
jx_fEq3h    = 108;
jx_fEq3f    = 109;
jx_fEq3hw   = 110;
jx_fEq4h    = 111;
jx_fEq4f    = 112;
jx_fEq4hw   = 113;
jx_fEq19h   = 114;
jx_fEq19f   = 115;
jx_fEq19hw  = 116;
jx_fEq20h   = 117;
jx_fEq20f   = 118;
jx_fEq20hw  = 119;

jx_Pq1h     = 120;
jx_Pq1f     = 121;
jx_Pq1hw    = 122;
jx_Pq2h     = 123;
jx_Pq2f     = 124;
jx_Pq2hw    = 125;
jx_Pq3h     = 126;
jx_Pq3f     = 127;
jx_Pq3hw    = 128;
jx_Pq4h     = 129;
jx_Pq4f     = 130;
jx_Pq4hw    = 131;
jx_Pq19h    = 132;
jx_Pq19f    = 133;
jx_Pq19hw   = 134;
jx_Pq20h    = 135;
jx_Pq20f    = 136;
jx_Pq20hw   = 137;
% valuation vec is: 1. current nfa position using expected kappa 
%                   2. discounted final nfa position (needs to be assigne from 1 in first iteration) 
%                   3. discounted valuation term  
%                   4. discounted nx term (but with current nx, needs to be shifted by one period) 
%                   5. discounted val rf term  
%                   6. discounted val rff term  
jx_val_nfa         = 138;
jx_val_EnfaEnd     = 139;
jx_val_Eval        = 140;
jx_val_Enx         = 141;
jx_val_Evalrf      = 142;
jx_val_Evalrff     = 143;
jx_val_Evalseign   = 144;
jx_val_Enfa        = 145;
jx_val_EEnfaEnd    = 146;
jx_val_EEval       = 147;
jx_val_EEnx        = 148;
jx_val_EEvalrf     = 149;
jx_val_EEvalrff    = 150;
jx_val_EEvalseign  = 151;



