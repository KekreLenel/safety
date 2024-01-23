%% State variables for numerical check
    
    % means 
        mean_k_unadj       = mean(k_series_unadj);
        mean_tht_h         = mean(tht_h_series);
        mean_wh_state      = mean(wh_state_series);
        mean_wf_state      = mean(wf_state_series);
        mean_rho_ih        = mean(rho_ih_series);
        mean_rho_if        = mean(rho_if_series);     
        mean_omg           = mean(omg_series);
        
    % std
        std_k_unadj       = std(k_series_unadj);
        std_tht_h_series  = std(tht_h_series);
        std_wh_state      = std(wh_state_series);
        std_wf_state      = std(wf_state_series);
        std_rho_ih        = std(rho_ih_series);
        std_rho_if        = std(rho_if_series);     
        std_omg           = std(omg_series);
        
        
       fprintf(fileID,'\nNUMERICAL VERIFICATION \n');
       fprintf(fileID,'-----------------------\n\n');
       fprintf(fileID,'STATES MEAN AND STD \n');
       fprintf(fileID,'-----------------------\n');
       fprintf(fileID,'        avrg.     grid_mean  std.  grid_dev\n');
       fprintf(fileID,'k       %6.2f    %6.2f   %6.2f    %6.2f\n',          mean_k_unadj   , k_grid_mean,                std_k_unadj     ,      k_grid_dev);
       fprintf(fileID,'tht_h   %6.2f%%   %6.2f%%  %6.2f%%   %6.2f%%\n', 100*mean_tht_h     , 100*tht_h_grid_mean,    100*std_tht_h_series,  100*tht_h_grid_dev);
       fprintf(fileID,'w_h     %6.2f    %6.2f   %6.2f    %6.2f\n',          mean_wh_state  , wh_grid_mean,               std_wh_state    ,       wh_grid_dev);
       fprintf(fileID,'w_f     %6.2f    %6.2f   %6.2f    %6.2f\n',          mean_wf_state  , wf_grid_mean,               std_wf_state    ,       wf_grid_dev);
       fprintf(fileID,'rho*ih  %6.2f%%   %6.2f%%  %6.2f%%   %6.2f%%\n', 100*mean_rho_ih    , ih_grid_mean*100,       100*std_rho_ih      ,    ih_grid_dev*100);
       fprintf(fileID,'rho*if  %6.2f%%   %6.2f%%  %6.2f%%   %6.2f%%\n', 100*mean_rho_if    , if_grid_mean*100,       100*std_rho_if      ,    if_grid_dev*100);
       fprintf(fileID,'omg     %6.2f%%   %6.2f%%  %6.2f%%   %6.2f%%\n', 100*mean_omg       , exp(omg_grid_mean)*100, 100*std_omg         ,    (exp(omg_grid_mean + omg_grid_dev) - exp(omg_grid_mean))*100);
       
