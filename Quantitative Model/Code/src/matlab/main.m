% -------------------------------------------------------------------------
% main.m: main result file, loops over comparative calibrations and 
% calls functions to read results and to produce tables and figures
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/safety
% -------------------------------------------------------------------------

clear all; close all;
tic
new_results = 1; % refresh results 

% numerical parameters
n_params         = 65;
n_sims           = 100;
n_sim_periods    = 400;
n_irf_periods    = 200;
n_sample_periods = 94;
n_irf            = 8;

% indices of calibration
ix_bm           = 1; 
ix_no_omg       = 2; 
ix_symm         = 3;
ix_symm_flex    = 4;
ix_tayl_y0      = 5;
ix_tayl_rho     = 6;
ix_no_omg_symm  = 7;        
ix_nocorr_nobg  = 8;
ix_nocorr       = 9;

% results folder 
addpath ../src/matlab
fig_path  = ['..', filesep, 'output', filesep, 'figures', filesep];
tab_path  = ['..', filesep, 'output', filesep, 'tables'];

% output files folder 
fid = fopen(['..',  filesep, 'output', filesep, 'tmp', filesep, 'n_comp.txt'], 'r');
n_comp = fscanf(fid,'%u');

% preallocate storage for tables in main text and appendix 
    table_1_params  = zeros(16,n_comp);   % externally set parameters 
    table_2_moms    = zeros(15,n_comp);   % targeted moments
    table_2_params  = zeros(15,n_comp);   % calibrated parameters 
    table_3_moms    = zeros(6,n_comp);    % comovements  
    table_4_moms    = zeros(6,n_comp);    % additional second moments 
    table_5_moms    = zeros(2,n_comp);    % output volatility 
    table_6_moms    = zeros(6,n_comp);    % net foreign assets
    table_7_moms    = zeros(6,n_comp);    % external adjustment 
    table_8_moms    = zeros(9,2,n_comp);  % dollar swap lines 
%   table 9: same moments as Table 3      % comovements under alternative policy rules 
    table_10_moms   = zeros(5,n_comp);    % portfolios and risk premia 

% initialize output series 
sim_series       = zeros(151,n_sim_periods,n_sims); % simulation 
sim_dis_series   = zeros(151,n_sim_periods,n_sims); % simulation with disaster 
p_irf_series     = zeros(151,n_irf_periods);        % p irf
omg_irf_series   = zeros(151,n_irf_periods);        % omg irf 
zf_irf_series   = zeros(151,n_irf_periods);        % zf irf
z_irf_series     = zeros(151,n_irf_periods);        % z irf
dis_irf_series   = zeros(151,n_irf_periods);        % dis irf
samp1_irf_series = zeros(137,n_irf_periods);        % no shock, fixed nominal rate
samp2_irf_series = zeros(137,n_irf_periods);        % no shock, fixed nominal rate
samp3_irf_series = zeros(137,n_irf_periods);        % no shock, fixed nominal rate
bg1_irf_series   = zeros(137,n_irf_periods);        % swp line response
fix_irf_series   = zeros(137,n_irf_periods);        % swp line, fixed nominal rate
fx0_irf_series   = zeros(137,n_irf_periods);        % no shock, fixed nominal rate

% impulse responses
collected_irfs  = zeros(n_irf_periods-2,37,n_irf,n_comp);

for ccc = 1:n_comp
     
    
    % output and parameter file paths
    data_path   = ['..', filesep, 'output', filesep, 'tmp', filesep, 'res_', num2str(ccc), filesep];
    param_file  = ['..', filesep, 'src', filesep, 'params', filesep, 'param_file_', num2str(ccc), '.csv'];

    % read results
    if new_results == 1
        read_results;
        save([data_path, 'data.mat'])
    else
        load([data_path, 'data.mat'])
    end

    % open output file
    fileID = fopen([tab_path, filesep, 'num_checks/results_', num2str(ccc), '.txt'],'w');
    fprintf(fileID,['RUN ', num2str(ccc), ' ', date, ' ', datestr(now, 'HH:MM:SS'), ' \n\n']);
    
    % write parameters into results file
    print_params; 
    
    %% business cycle moments and tables
    fprintf(fileID,'\n\nWITHOUT DISASTER REALIZATION IN SAMPLE\n');
    fprintf(fileID,'==================================================\n\n');
    
        dis = false;

        % numerical check
        temp_series = reshape(sim_series(:,:,:),size(sim_series,1),n_sims*n_sim_periods)';
        extract_series;
        numerical_check;

        % loop over samples, calc moments
        for sss = 1:n_sims 
            temp_series = sim_series(:,:,sss)';
            extract_series;
            calc_moments;
        end
        
        fprintf(fileID,'\n\nWITH DISASTER REALIZATION IN SAMPLE\n');
        fprintf(fileID,'==================================================\n\n');
    
        dis = true;
        temp_series = reshape(sim_dis_series,size(sim_series,1),n_sims*n_sim_periods)';
        extract_series;
        numerical_check;
        
        for sss = 1:n_sims 
            temp_series = sim_dis_series(:,:,sss)';
            extract_series;
            calc_moments;
        end

        % average moments 
        collect_moments;

    fclose(fileID);

    % Impulse responses 
        
    [collected_irfs(:,:,1,ccc), irf_idxs, irf_titles] = extract_irfs(p_irf_series',    trg_prms);
    [collected_irfs(:,:,2,ccc), ~, ~]                 = extract_irfs(omg_irf_series',  trg_prms); 
    [collected_irfs(:,:,3,ccc), ~, ~]                 = extract_irfs(g_irf_series',  trg_prms);  
    [collected_irfs(:,:,4,ccc), ~, ~]                 = extract_irfs(dis1_irf_series',  trg_prms); 
    [collected_irfs(:,:,5,ccc), ~, ~]                 = extract_irfs(zf_irf_series',  trg_prms);  

    [collected_irfs(:,:,6,ccc), ~, ~]                 = extract_irfs(bg1_irf_series',  trg_prms);
    [collected_irfs(:,:,7,ccc), ~, ~]                 = extract_irfs(fix_irf_series',  trg_prms);
    [collected_irfs(:,:,8,ccc), ~, ~]                 = extract_irfs(fx0_irf_series',  trg_prms);
     collected_irfs(:,:,7,ccc)                        = collected_irfs(:,:,7,ccc) - collected_irfs(:,:,8,ccc);
    
    if ccc == ix_bm
        collect_swap_moments; 
        create_recession_fig;
        create_safety_fig;
    end


end

create_tables;
create_figures;

toc
