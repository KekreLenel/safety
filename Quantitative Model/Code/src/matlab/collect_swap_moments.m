
% moments with fixed nominal rate 
tmp_irfs = collected_irfs(:,:,7,ccc);
scale_irfs = -14/tmp_irfs(1,irf_idxs.omg); % scale so that initial omega response is 14bp
tmp_irfs = scale_irfs*tmp_irfs;
table_8_moms(1,1,ccc) = tmp_irfs(1, irf_idxs.log_E);
table_8_moms(2,1,ccc) = tmp_irfs(1, irf_idxs.retA)+tmp_irfs(1, irf_idxs.log_infl_h);
table_8_moms(3,1,ccc) = max(tmp_irfs(:, irf_idxs.log_yh));
table_8_moms(4,1,ccc) = max(tmp_irfs(:, irf_idxs.log_yf));
table_8_moms(5,1,ccc) = tmp_irfs(1, irf_idxs.log_qx);
table_8_moms(6,1,ccc) = tmp_irfs(1, irf_idxs.excA);
table_8_moms(7,1,ccc) = tmp_irfs(1, irf_idxs.retA);
table_8_moms(8,1,ccc) = tmp_irfs(1, irf_idxs.nfa_y);
table_8_moms(9,1,ccc) = tmp_irfs(1, irf_idxs.nfa_y)-tmp_irfs(1, irf_idxs.nx_y);

% moments without fixed nominal rate 
tmp_irfs = collected_irfs(:,:,6,ccc);
scale_irfs = -14/tmp_irfs(1,irf_idxs.omg); % scale so that initial omega response is 14bp
tmp_irfs = scale_irfs*tmp_irfs;
table_8_moms(1,2,ccc) = tmp_irfs(1, irf_idxs.log_E);
table_8_moms(2,2,ccc) = tmp_irfs(1, irf_idxs.retA)+tmp_irfs(1, irf_idxs.log_infl_h);
table_8_moms(3,2,ccc) = max(tmp_irfs(:, irf_idxs.log_yh));
table_8_moms(4,2,ccc) = max(tmp_irfs(:, irf_idxs.log_yf));
table_8_moms(5,2,ccc) = tmp_irfs(1, irf_idxs.log_qx);
table_8_moms(6,2,ccc) = tmp_irfs(1, irf_idxs.excA);
table_8_moms(7,2,ccc) = tmp_irfs(1, irf_idxs.retA);
table_8_moms(8,2,ccc) = tmp_irfs(1, irf_idxs.nfa_y);
table_8_moms(9,2,ccc) = tmp_irfs(1, irf_idxs.nfa_y)-tmp_irfs(1, irf_idxs.nx_y);

