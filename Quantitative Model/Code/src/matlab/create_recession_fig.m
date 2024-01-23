% -------------------------------------------------------------------------
% create_recession_fig
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/safety
% -------------------------------------------------------------------------

%% READ DATA

    data_tmp     = csvread('../src/data/modelvsdata.csv',1,0);
    data_datevec = datenum(1995,1:3:309,1);
    
    % disaster probability: remove zeros, scale as in model
    p_series_data    = data_tmp(1:94,4);
    p_series_data = data_tmp(:,4);
    p_series_data(p_series_data==0) = min(p_series_data(p_series_data>0));
    p_series_data = p_series_data/mean(p_series_data)*param_input_vec(1,1);
    
    % safety: remove mean, scale as in model
    omg_series_data  = data_tmp(1:100,5)/(100*100);
    omg_series_data = omg_series_data - mean(omg_series_data);
    omg_scaling = param_input_vec(4,1)/std(omg_series_data);
    omg_series_data = omg_series_data* param_input_vec(4,1)/std(omg_series_data);
    
    dlog_q_data  = data_tmp(:,6);
    dlog_ch_data = data_tmp(:,7);
    dlog_cf_data = data_tmp(:,8);
    dlog_yh_data = data_tmp(:,9);
    dlog_yf_data = data_tmp(:,10);
    nfay_data    = data_tmp(:,11);
    dnfay_data   = data_tmp(:,12);
    nxy_data     = data_tmp(:,13);
    nom_ih_data  = data_tmp(:,14); 
    nom_if_data  = data_tmp(:,15)+ data_tmp(:,14);
        
        
    start_idx = find(data_datevec >= datenum(2006,4,1),1);
    end_idx = find(data_datevec >= datenum(2011,09,1),1);

    y_vec = [];
        
    temp_series = samp1_series(:,:,1)';
    extract_series; 

    y_vec(:,1,1) = 100*p_series(start_idx:end_idx);
    y_vec(:,2,1) = 100*p_series(start_idx:end_idx); 
    
    y_vec(:,2,2) = 100*omg_series(start_idx:end_idx);
    y_vec(:,1,2) = 100*omg_series(start_idx:end_idx); 
    
    y_vec(:,1,3) = 100*dlog_yh_data(start_idx:end_idx);
    y_vec(:,2,3) = 100*diff(log(yh_series(start_idx-1:end_idx)));
    
    y_vec(:,1,4) = 100*dlog_yf_data(start_idx:end_idx);
    y_vec(:,2,4) = 100*diff(log(yf_series(start_idx-1:end_idx)));
    
    tmp1 = max([max(100*(nfa_rel_growth_series(start_idx:end_idx))), max(100*(dnfay_data(start_idx:end_idx))*4), ...
                max(100*(nx_rely_series(start_idx:end_idx))), max(100*nxy_data(start_idx:end_idx)*4)]);
    tmp2 = min([min(100*(nfa_rel_growth_series(start_idx:end_idx))), min(100*(dnfay_data(start_idx:end_idx))*4), ...
                min(100*(nx_rely_series(start_idx:end_idx))), min(100*nxy_data(start_idx:end_idx)*4)]);

    y_vec(:,1,5) = 100*(dnfay_data(start_idx:end_idx))*4 ; 
    y_vec(:,2,5) = 100*(nfa_rel_growth_series(start_idx:end_idx));
    
    y_vec(:,1,6) = 100*(dnfay_data(start_idx:end_idx))*4 - 100*nxy_data(start_idx:end_idx)*4;
    y_vec(:,2,6) = 100*(nfa_rel_growth_series(start_idx:end_idx) - nx_rely_series(start_idx:end_idx));
    

    y_vec(:,1,7) = nom_ih_data(start_idx:end_idx,:)/4;
    y_vec(:,2,7) = nom_ih_series(start_idx:end_idx,:)*100-100;
    y_vec(:,5,7) = nom_if_data(start_idx:end_idx,:)/4;
    y_vec(:,6,7) = nom_if_series(start_idx:end_idx,:)*100-100;
    
    y_vec(:,1,8) = nom_if_data(start_idx:end_idx,:)/4;
    y_vec(:,2,8) = nom_if_series(start_idx:end_idx,:)*100-100;
    y_vec(:,5,8) = nom_ih_data(start_idx:end_idx,:)/4;
    y_vec(:,6,8) = nom_ih_series(start_idx:end_idx,:)*100-100;
    
    

    y_vec(:,1,9)  = NaN; % infl_data(start_idx:end_idx,:);
    y_vec(:,2,9)  = log(infl_h_series(start_idx:end_idx,:))*100;
    y_vec(:,5,9)  = log(infl_f_series(start_idx:end_idx,:))*100;

    y_vec(:,1,10) = NaN;
    y_vec(:,2,10) = log(infl_f_series(start_idx:end_idx,:))*100;
    y_vec(:,5,10) = log(infl_h_series(start_idx:end_idx,:))*100;

    y_vec(:,1,11)  = NaN;
    y_vec(:,2,11)  = log(E_rfh_series(start_idx:end_idx,:))*100;
    y_vec(:,5,11)  = log(E_rff_f_series(start_idx:end_idx,:))*100;

    y_vec(:,1,12) = NaN;
    y_vec(:,2,12) = log(E_rff_f_series(start_idx:end_idx,:))*100;
    y_vec(:,5,12) = log(E_rfh_series(start_idx:end_idx,:))*100;
    
    y_vec(:,1,13)  = NaN;
    y_vec(:,2,13)  = (yield4_hw_series(start_idx:end_idx,:))*100;
    y_vec(:,5,13)  = (yield4_f_series(start_idx:end_idx,:))*100;

    y_vec(:,1,14) = NaN;
    y_vec(:,2,14) = (yield4_f_series(start_idx:end_idx,:))*100;
    y_vec(:,5,14) = (yield4_hw_series(start_idx:end_idx,:))*100;
                
    y_vec(:,1,15) = dlog_ch_data(start_idx:end_idx,:)*100;
    y_vec(:,2,15) = (c_h_change_series(start_idx:end_idx,:))*100;

    y_vec(:,1,16) = dlog_cf_data(start_idx:end_idx,:)*100;
    y_vec(:,2,16) = (c_f_change_series(start_idx:end_idx,:))*100;

    y_vec(:,1,17) = (dlog_ch_data(start_idx:end_idx,:) + dlog_q_data(start_idx:end_idx,:) - dlog_cf_data(start_idx:end_idx,:))*100;
    y_vec(:,2,17) = (c_h_change_series(start_idx:end_idx,:) + qx_change_series(start_idx:end_idx,:) - c_f_change_series(start_idx:end_idx,:))*100;
    
    y_vec(:,1,18) = NaN;
    y_vec(:,2,18) = (tht_h_series(start_idx:end_idx,:))*100;
    y_vec(:,5,18) = (tht_h_series(start_idx:end_idx,:))*100;
    y_vec(:,6,18) = (tht_h_series(start_idx:end_idx,:))*100;

    y_vec(:,1,19) = cumsum(dlog_ch_data(start_idx:end_idx,:) + dlog_q_data(start_idx:end_idx,:) - dlog_cf_data(start_idx:end_idx,:))*100;
    y_vec(:,1,19) = y_vec(:,1,19) - y_vec(1,1,19);
    y_vec(:,2,19) = cumsum(c_h_change_series(start_idx:end_idx,:) + qx_change_series(start_idx:end_idx,:) - c_f_change_series(start_idx:end_idx,:))*100;
    y_vec(:,2,19) = y_vec(:,2,19) - y_vec(1,2,19);

    temp_series = samp2_series';
    extract_series; 
    y_vec(:,3,1) = 100*p_series(start_idx:end_idx); 
    y_vec(:,3,2) = 100*omg_series(start_idx:end_idx); 
    y_vec(:,3,3) = 100*diff(log(yh_series(start_idx-1:end_idx)));
    y_vec(:,3,4) = 100*diff(log(yf_series(start_idx-1:end_idx)));
    y_vec(:,3,5) = 100*(nfa_rel_growth_series(start_idx:end_idx));
    y_vec(:,3,6) = 100*(nfa_rel_growth_series(start_idx:end_idx) - nx_rely_series(start_idx:end_idx));
    y_vec(:,3,7) = nom_ih_series(start_idx:end_idx,:)*100-100;
    y_vec(:,3,8) = nom_if_series(start_idx:end_idx,:)*100-100;
    y_vec(:,3,9)  = log(infl_h_series(start_idx:end_idx,:))*100;
    y_vec(:,3,10) = log(infl_f_series(start_idx:end_idx,:))*100;       
    y_vec(:,3,11) = log(E_rfh_series(start_idx:end_idx,:))*100;     
    y_vec(:,3,12) = log(E_rff_f_series(start_idx:end_idx,:))*100;     
    y_vec(:,3,13) = (yield4_hw_series(start_idx:end_idx,:))*100;     
    y_vec(:,3,14) = (yield4_f_series(start_idx:end_idx,:))*100;     
    
    temp_series = samp3_series';
    extract_series; 
    y_vec(:,4,1) = 100*p_series(start_idx:end_idx); 
    y_vec(:,4,2) = 100*omg_series(start_idx:end_idx); 
    y_vec(:,4,3) = 100*diff(log(yh_series(start_idx-1:end_idx)));
    y_vec(:,4,4) = 100*diff(log(yf_series(start_idx-1:end_idx)));
    y_vec(:,4,5) = 100*(nfa_rel_growth_series(start_idx:end_idx));
    y_vec(:,4,6) = 100*(nfa_rel_growth_series(start_idx:end_idx) - nx_rely_series(start_idx:end_idx));
    y_vec(:,4,7) = nom_ih_series(start_idx:end_idx,:)*100-100;
    y_vec(:,4,8) = nom_if_series(start_idx:end_idx,:)*100-100;
    y_vec(:,4,9)  = log(infl_h_series(start_idx:end_idx,:))*100;
    y_vec(:,4,10) = log(infl_f_series(start_idx:end_idx,:))*100;                
    y_vec(:,4,11) = log(E_rfh_series(start_idx:end_idx,:))*100;     
    y_vec(:,4,12) = log(E_rff_f_series(start_idx:end_idx,:))*100;     
    y_vec(:,4,13) = (yield4_hw_series(start_idx:end_idx,:))*100;     
    y_vec(:,4,14) = (yield4_f_series(start_idx:end_idx,:))*100;                
            
    y_array = [];
    for kkk = 1:19
        y_array{kkk} = squeeze(y_vec(:,:,kkk));
    end

    x_vec      = []; x_array = [];
    x_vec(:,1) = data_datevec(start_idx:end_idx)';
    x_vec(:,2) = data_datevec(start_idx:end_idx)';
    x_vec(:,3) = data_datevec(start_idx:end_idx)';
    x_vec(:,4) = data_datevec(start_idx:end_idx)';
    x_vec(:,5) = data_datevec(start_idx:end_idx)';
    x_vec(:,6) = data_datevec(start_idx:end_idx)';
    for iii = 1:19 
    x_array{iii} =  x_vec;
    end

    title_array  = {'$p$', '$\omega$', '$ \Delta \log y$','$ \Delta \log y^\ast$', '$ \Delta nfa/y $', '$(\Delta nfa - nx)/y$'};
    ylabel_array = {'$\%$', '$\%$', '$\%$', '$\%$', '$\%$', '$\%$'};
     legend_array = {'Model','Data'}; legend_loc1 = 6; legend_loc2 = 'se';          
    fig_handle = make_figure([2,3],{2,1},{2,1,3,4,5,6},x_array,y_array, title_array, ylabel_array,legend_array, legend_loc1, legend_loc2, 1,2,'modelvsdata');
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0], 'PaperSize', [8.5 5.0]);
    print(fig_handle,[fig_path, '/fig_4'],'-dpdf', '-r600')
 
    legend_array = {'only $p$','only $\omega$','Model'}; legend_loc1 = 6; legend_loc2 = 'se';          
    fig_handle = make_figure([2,3],{3,4,2},{3,4,2,5,6},x_array,y_array, title_array, ylabel_array,legend_array, legend_loc1, legend_loc2, 1,2,'modelvsdata');
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0], 'PaperSize', [8.5 5.0]);
    print(fig_handle,[fig_path, '/fig_6'],'-dpdf', '-r600')

    title_array  = {'$i$', '$i^\ast$', '$\theta$'};
    ylabel_array = {'$\%$', '$\%$', '$\%$'};
    legend_array = {'Model', 'Data'}; legend_loc1 = 3; legend_loc2 = 'ne';
    fig_handle = make_figure([1,3], {2,1}, {2,1,5,6}, x_array([7,8,18]),y_array([7,8,18]),title_array, ylabel_array, legend_array, legend_loc1, legend_loc2, 1, 2, 'modelvsdata');
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 2.25], 'PaperSize', [8.5 2.25]);
    print(fig_handle,[fig_path,'/fig_7'],'-dpdf', '-r600')
    
    %% create recession table
    temp_series = samp1_series';
    extract_series; 

    yh_data = cumsum(dlog_yh_data);
    yf_data = cumsum(dlog_yf_data);
    
    loc_fileID = fopen([tab_path, '/table_A1.tex'],'w');
    fprintf(loc_fileID, " \\emph{Q3 07 - Q3 09} & Data & Model \\\\ \\hline \n");
    
    date_start = datenum(2007,7,1);
    date_end = datenum(2009,7,1);
    start_idx = find(data_datevec>=date_start,1);
    end_idx = find(data_datevec>=date_end,1);
    
    tmp1 = 100*(log((yh_series(end_idx))) - log(yh_series(start_idx)));
    tmp2 = 100*(((yh_data(end_idx))) - (yh_data(start_idx)));    
    fprintf(loc_fileID,'$\\Delta \\log y$ & %10.4f\\%% & %10.4f\\%%  \\\\ \n', tmp2, tmp1);

    tmp1 = 100*(log((yf_series(end_idx))) - log(yf_series(start_idx)));
    tmp2 = 100*(((yf_data(end_idx))) - (yf_data(start_idx)));
    fprintf(loc_fileID,'$\\Delta \\log y^*$ & %10.4f\\%% & %10.4f\\%%  \\\\ \n', tmp2, tmp1);

    tmp1 = 100*(((nfa_rel_series(end_idx))) - (nfa_rel_series(start_idx)));
    tmp2 = 100*(((nfay_data(end_idx))) - (nfay_data(start_idx)));
    fprintf(loc_fileID,'$\\Delta nfa/y$ & %10.4f\\%% & %10.4f\\%%  \\\\  \\hline \n', tmp2, tmp1);
    
    fclose(loc_fileID);
    
    
