 fprintf(fileID,'PARAMETRIZATION\n');
fprintf(fileID,'-----------------------\n\n');
for ppp = 1:n_params
   fprintf(fileID,'%-16s   %8.4f\n', targ_prm_names{ppp}, targ_prm_table(ppp));
end

% Add parameters calibrated in steady state calculation
fprintf(fileID,'\nCALIBRATED PARAMETERS\n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'Labor disutility 1   %8.4f\n', chi0_vec(1));
fprintf(fileID,'Labor disutility 2   %8.4f\n', chi0_vec(2));

trg_prms.chi0_h = chi0_vec(1);
trg_prms.chi0_f = chi0_vec(2);
trg_prms.yh_ss   = yh_ss;


% add relevant steady state moments
fprintf(fileID,'\nSTEADY STATE VALUES\n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'yh_ss       %8.4f\n', yh_ss);
fprintf(fileID,'aggr_css   %8.4f\n', aggr_css);

% add parameters of stoch processes
fprintf(fileID,'\nLOG NORMAL PROCESSES\n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'log(omg) innovation std   %8.4f\n', sig_omg);
fprintf(fileID,'log(omg) mean             %8.4f\n', omg_grid_mean);
fprintf(fileID,'log(p)   innovation std   %8.4f\n', sig_dis);
fprintf(fileID,'log(p)   mean             %8.4f\n', dis_grid_mean);

% add target values for stoch process
fprintf(fileID,'\nLOG NORMAL PROCESSES INPUTS\n');
fprintf(fileID,'-----------------------\n\n');
fprintf(fileID,'mean p   %8.4f\n', param_input_vec(1,ccc));
fprintf(fileID,'std  p   %8.4f\n', param_input_vec(2,ccc));
fprintf(fileID,'mean wdg %8.4f\n', param_input_vec(3,ccc));
fprintf(fileID,'std  wdg %8.4f\n', param_input_vec(4,ccc));
fprintf(fileID,'skew wdg %8.4f\n', param_input_vec(5,ccc));
