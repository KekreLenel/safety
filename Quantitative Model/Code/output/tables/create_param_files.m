%% Script creates parameter files & prepares folder structure

    clear all;
    addpath ../src/params/
    tic 
%% benchmark parameters    
    
% relative country sizes foreign to home 

    prm.zeta           = 1.6;

% preference parameters
    
    prm.bbeta_h        = 0.98925*0.999452; 
    prm.bbeta_f        = 0.98925;

    prm.gma_h          = 21.0;
    prm.gma_f          = 23.8; 

    prm.ies_h          = 0.75;
    prm.ies_f          = 0.75;

    prm.chi            = 0.75;

    prm.sigma           = 1.5;
    home_bias           = 0.4;
    prm.varsigma_h      = 1.0/(1.0+prm.zeta) + home_bias;
    prm.varsigma_f      = 1.0/(1.0+prm.zeta) - home_bias/prm.zeta;

    prm.l_target        = 1.0;

% production parameters 

    prm.delta           = 0.025;
    prm.aalpha          = 0.33;
    prm.chiX            = 3.0;
    prm.inv_share_h     = 1.0/(1.0+prm.zeta); 

% wage stickiness 

    prm.varphi_w        = 400;
    prm.vareps_w        = 20;

% monetary policy rule 

    prm.phi_pi_h        = 1.5;
    prm.phi_pi_f        = 1.5;
    prm.phi_y_h		    = 0.0;
    prm.phi_y_f		    = 0.0;
    prm.tayl_ic_h       = 0.0049;
    prm.tayl_ic_f       = 0.0047;
    prm.rho_i           = 0.0;

% government debt parameters 

    prm.b_lmbd          = 14.0/50.0; 
    prm.bg_yss          = 0.127;

% stochastic processes 

    prm.sig_z           = 0.00165;

    prm.std_zf         = 0.0114; % this is the unconditional std, sig_zf = std_zf*sqrt(1-rho_zf^2) 
    prm.rho_zf         = 0.90;
    
    % p process
    avrg_p_trgt = 0.5/100;
    std_p_trgt  = avrg_p_trgt;
    [X, ~, FLAG]     = fsolve(@(x) 1000*[ avrg_p_trgt - exp(x(1) + x(2)/2) ; std_p_trgt^2 - (exp(x(2))-1)*exp(2*x(1) + x(2)) ], [-1; 1]);
    if (FLAG ~= 1)
        error('no solution for p estimation');
    end
    mu_log_p  = X(1);
    std_log_p = sqrt(X(2));    
    
    prm.disast_p        = mu_log_p;
    prm.disast_shock    = 0.1;
    prm.p_rho           = 0.75;
    prm.p_std           = std_log_p;
   
    % omega process  
    avrg_omg_trgt = 0.0;
    std_omg_trgt  = 0.0048;
    skew_omg_trgt = 6.1;
    [X, ~, FLAG]  = fsolve(@(x) 10000*[ skew_omg_trgt - (exp(x(2)) + 2)*sqrt(exp(x(2)) - 1) ; std_omg_trgt^2 - (exp(x(2))-1)*exp(2*x(1) + x(2)) ], [-1; 1]);
    if (FLAG ~= 1)
        error('no solution for omg estimation');
    end    
    mu_log_omg  = X(1);
    std_log_omg = sqrt(X(2));    
    omg_shift   = avrg_omg_trgt - exp(X(1) + X(2)/2);
    
    prm.omg_mean      = mu_log_omg;
    prm.omg_std       = std_log_omg;
    prm.omg_shift     = omg_shift;
    prm.omg_rho       = 0.4;
    
    prm.corr_omg_dis    = 0.5; 

% grid scale and position parameters     
    prm.tht_trgt_h      = 1.0/(1 + prm.zeta) - 0.035; 
    prm.theta_h_dev     = 0.10;
    prm.k_dev_param     = 0.1; 
    prm.k_grid_adj      = 1.025;
    prm.w_dev_param     = 0.05;     
    prm.w_grid_adj      = 1.01;
    prm.ih_grid_mean    = 0.005;
    prm.if_grid_mean    = 0.005;
    prm.ih_grid_dev     = 0.05;
    prm.if_grid_dev     = 0.05;

% Smolyak grid density and quadrature nodes 
    prm.mu_k   = 3;
    prm.mu_tht = 3; 
    prm.mu_zf = 3;
    prm.mu_wh  = 3; 
    prm.mu_wf  = 3; 
    prm.mu_p   = 3; 
    prm.mu_ih  = 0;
    prm.mu_if  = 0;
    prm.mu_omg = 3; 

% quadrature parameters 
    prm.n_quad_z   = 3;  
    prm.n_quad_zf = 3;  
    prm.n_quad_p   = 5;  
    prm.n_quad_omg = 5;  

% foreign bond trading active? 
    prm.foreign_trading = 1.0;
    
% run simulation and swap line IRFs
    prm.run_bg   = 0.0; 
    prm.run_samp = 0.0;

%% Prepare comparative statics
   
    % how many comparative statics
    n_comp      = 9; 
    
    % number of parameters
    n_params = size(fieldnames(prm),1);

    % additional inputs stored
    param_input_vec = zeros(5,n_comp);
    param_input_vec(:,1) = [avrg_p_trgt, std_p_trgt, avrg_omg_trgt, std_omg_trgt, skew_omg_trgt];


    fid = fopen('../output/tmp/n_comp.txt', 'w');
    fprintf(fid, '%i', n_comp);
    fclose(fid);


% Which parameters to vary

% Benchmark calibration 
prm_list{1} = prm;
prm_list{1}.run_samp = 1.0;
prm_list{1}.run_bg   = 1.0;

% No omg shock 
prm_list{2} = prm;
prm_list{2}.n_quad_omg = 1;
prm_list{2}.mu_omg = 0;
prm_list{2}.omg_mean    = log(0.01); % omg mean 
prm_list{2}.omg_shift   = -0.01;     % omg shift 
prm_list{2}.omg_std     = 0.0;       % omg_std 
prm_list{2}.corr_omg_dis= 1.0;       % corr omg dis (only for numerical reasons = 1) 
prm_list{2}.bbeta_h     = prm.bbeta_f*0.999537; % % 0.999557 20.7 0.999558 20.3
param_input_vec(:,2) = [avrg_p_trgt, std_p_trgt, 0.0, 0.0, 0.0];

% Symmetric case 
prm_list{3} = prm;
prm_list{3}.gma_h = (prm.tht_trgt_h/prm.gma_h  +(1 - prm.tht_trgt_h)/prm.gma_f)^(-1);
prm_list{3}.gma_f = (prm.tht_trgt_h/prm.gma_h  +(1 - prm.tht_trgt_h)/prm.gma_f)^(-1);
prm_list{3}.tht_trgt_h = 1.0/(1.0+prm.zeta);
prm_list{3}.bbeta_h = prm.bbeta_f;
param_input_vec(:,3) = param_input_vec(:,1);

% Calculate symmetric case, no stickiness
prm_list{4} = prm;
prm_list{4}.varphi_w    = 0.0;
prm_list{4}.gma_h = (prm.tht_trgt_h/prm.gma_h  +(1 - prm.tht_trgt_h)/prm.gma_f)^(-1);
prm_list{4}.gma_f = (prm.tht_trgt_h/prm.gma_h  +(1 - prm.tht_trgt_h)/prm.gma_f)^(-1);
prm_list{4}.tht_trgt_h = 1.0/(1.0+prm.zeta);
prm_list{4}.bbeta_h = prm.bbeta_f;
param_input_vec(:,4) = param_input_vec(:,1);

% output gap in the taylor rule
prm_list{5} = prm;
prm_list{5}.sig_z   = 0.0033;
prm_list{5}.std_zf = 0.0098; % 0.006*2.6/1.6;
prm_list{5}.bbeta_h = prm.bbeta_f*0.9995050;  
prm_list{5}.gma_h  = 21.2;
prm_list{5}.gma_f  = 23.8;
prm_list{5}.tayl_ic_h = -0.1312; % -0.0118;
prm_list{5}.tayl_ic_f = -0.1892; % -0.0698;
prm_list{5}.phi_y_h = 0.125;
prm_list{5}.phi_y_f = 0.125; 
param_input_vec(:, 5) = param_input_vec(:, 1);

% inertial taylor rule and output gap
prm_list{6} = prm;
prm_list{6}.sig_z   = 0.0001;
prm_list{6}.std_zf = 0.0081; %0.005*2.6/1.6;
prm_list{6}.bbeta_f = prm.bbeta_f*0.9997; 
prm_list{6}.bbeta_h = prm.bbeta_f*0.9997*0.9995400; %  99956 -34 % 600 0 % 5500 -34 % -35 5450 % -9 5350 % 2 0.9995650 % -32 0.9995050
prm_list{6}.gma_h  = 21.8;
prm_list{6}.gma_f  = 24.0;
prm_list{6}.chiX   = 3.0;
prm_list{6}.tayl_ic_h = -0.1302; % -0.0108;
prm_list{6}.tayl_ic_f = -0.1877; % -0.0683; 
prm_list{6}.phi_y_h = 0.125;
prm_list{6}.phi_y_f = 0.125; 
prm_list{6}.rho_i = 0.5; 
prm_list{6}.mu_ih  = 3;
prm_list{6}.mu_if  = 3;
param_input_vec(:, 6) = param_input_vec(:, 1);


% calculate case with no omega shocks symmetry 
prm_list{7} = prm;
prm_list{7}.n_quad_omg = 1;
prm_list{7}.mu_omg = 0;
prm_list{7}.omg_mean    = log(0.01); % omg mean 
prm_list{7}.omg_shift   = -0.01;     % omg shift 
prm_list{7}.omg_std     = 0.0;       % omg_std 
prm_list{7}.corr_omg_dis= 1.0;       % corr omg dis (only for numerical reasons = 1) 
prm_list{7}.gma_h = (prm.tht_trgt_h/prm.gma_h  +(1 - prm.tht_trgt_h)/prm.gma_f)^(-1);
prm_list{7}.gma_f = (prm.tht_trgt_h/prm.gma_h  +(1 - prm.tht_trgt_h)/prm.gma_f)^(-1);
prm_list{7}.tht_trgt_h = 1.0/(1.0+prm.zeta);
prm_list{7}.bbeta_h     = prm.bbeta_f; % *0.999495;
param_input_vec(:,7) = [avrg_p_trgt, std_p_trgt, 0.0, 0.0, 0.0];

% no correlation between omg and disaster shock  no government debt 
prm_list{8} = prm;
prm_list{8}.bbeta_h = prm.bbeta_f*0.999541; %  
prm_list{8}.corr_omg_dis = 0.0;       
prm_list{8}.bg_yss       = 0.0; 
param_input_vec(:,8) = param_input_vec(:,1);

% no correlation between shocks  
prm_list{9} = prm;
prm_list{9}.bbeta_h = prm.bbeta_f*0.999541; 
prm_list{9}.corr_omg_dis = 0.0;       
param_input_vec(:,9) = param_input_vec(:,1);


%% Store parameter files

    % Now generate parameter files
    for nnn = 1:n_comp

        name_list = string(fieldnames(prm_list{nnn}));

        prm_vec = cell2mat(struct2cell(prm_list{nnn}));

        file_str = sprintf('../src/params/param_file_%i.csv',nnn);
        fid = fopen(file_str, 'w');
        fprintf(fid, [ repmat('%s, ', [1, n_params-1]  ), ' %s\n'], name_list(1:end-1), name_list(end));
        fprintf(fid, [ repmat('%.6f, ', [1, n_params-1]  ), ' %.6f\n'], prm_vec(1:end-1), prm_vec(end));
        fclose(fid);

    end
    
    param_input_vec = param_input_vec(:,1:n_comp);
    save('../src/params/additional_params','param_input_vec')



%% store shock series for recession simulation
    % load raw series 
    data_tmp = csvread('../src/data/modelvsdata.csv',1,0);
    p_data    = data_tmp(1:94,4);
    omg_data  = data_tmp(1:100,5)/(100*100);
    
    % scale to target 
    p_data = p_data/mean(p_data)*avrg_p_trgt;

    trgt_moments    = [];
    trgt_moments(1) = corr(p_data(2:end), p_data(1:end-1));
    trgt_moments(2) = corr(omg_data(2:end), omg_data(1:end-1));
    trgt_moments(3) = corr(omg_data(1:94),p_data);    
    
    % remove zeros in p data (for feeding into model)
    p_data(p_data<0.05*mean(p_data)) = 0.05*mean(p_data);
    log_p_data = log(p_data);
    log_p_data = log_p_data - mu_log_p;
    p_shocks = log_p_data(2:end) - prm.p_rho*log_p_data(1:end-1);
    p_shocks = [log_p_data(1)+mu_log_p;p_shocks];
    
    log_p_avrg = log(ones(size(log_p_data))*avrg_p_trgt)  - mu_log_p;
    p_shocks_0 = log_p_avrg(2:end) - prm.p_rho*log_p_avrg(1:end-1);
    p_shocks_0 = [log_p_avrg(1)+mu_log_p;p_shocks_0];    
    
    omg_data = omg_data(1:94) - mean(omg_data);
    omg_data = omg_data*std_omg_trgt/std(omg_data);
    omg_data = omg_data - omg_shift;
    omg_data(omg_data < 0.0) = 0.0005;
    log_omg_data = log(omg_data);
    log_omg_data = log_omg_data - mu_log_omg;
    omg_shocks = log_omg_data(2:end) - prm.omg_rho*log_omg_data(1:end-1);
    omg_shocks = [log_omg_data(1)+mu_log_omg; omg_shocks];
    
    log_omg_avrg = log(avrg_omg_trgt*ones(size(omg_data)) - omg_shift)- mu_log_omg;
    omg_shocks_0 = log_omg_avrg(2:end) - prm.omg_rho*log_omg_avrg(1:end-1);
    omg_shocks_0 = [log_omg_avrg(1)+mu_log_omg; omg_shocks_0];    
    
    n_sample = size(omg_shocks,1);
    n_shocks = 4; 
    
    sample_mat = zeros(n_shocks,n_sample);
    sample_mat(3,:) = p_shocks;
    sample_mat(4,:) = omg_shocks;
    file_str = sprintf('../src/params/sample1_mat.txt');
    fid = fopen(file_str, 'w');
    fprintf(fid, [ repmat('%10.6f,  ', [1, n_sample*n_shocks]  )], sample_mat(:));
    fclose(fid);    
    
    n_sample   = size(omg_shocks,1);
    sample_mat = zeros(n_shocks,n_sample);
    sample_mat(3,:) = p_shocks;
    sample_mat(4,:) = omg_shocks_0;
    file_str = sprintf('../src/params/sample2_mat.txt');
    fid = fopen(file_str, 'w');
    fprintf(fid, [ repmat('%10.6f,  ', [1, n_sample*n_shocks]  )], sample_mat(:));
    fclose(fid);    
    
    n_sample   = size(omg_shocks,1);
    sample_mat = zeros(n_shocks,n_sample);
    sample_mat(3,:) = p_shocks_0;
    sample_mat(4,:) = omg_shocks;
    file_str = sprintf('../src/params/sample3_mat.txt');
    fid = fopen(file_str, 'w');
    fprintf(fid, [ repmat('%10.6f,  ', [1, n_sample*n_shocks]  )], sample_mat(:));
    fclose(fid);        
    
    
%% test accuracy of calibration 

    test_calib = 1;
    if test_calib == 1 
        
        n_sim = 1000;
        n_sim_periods = 100;
        
        sig_p     = std_log_p   *sqrt(1 -  prm.p_rho^2);
        sig_omg   = std_log_omg *sqrt(1 -  prm.omg_rho^2);
        cov_p_omg  = prm.corr_omg_dis*sig_p*sig_omg;
        chol_mat = chol([sig_p^2, cov_p_omg; cov_p_omg, sig_omg^2]);
        
        rng('default');
        sim_moments = zeros(3,n_sim);
        for sss = 1:n_sim

            sim_shock_mat = randn(n_sim_periods,2)*chol_mat;

            a = [1, -prm.p_rho];                   
            b = 1;                       
            p_sim = filter(b,a, sim_shock_mat(:,1)); 
            a = [1, -prm.omg_rho];                   
            b = 1;                       
            omg_sim = filter(b,a, sim_shock_mat(:,2)); 

            sim_moments(1,sss) = corr(exp(p_sim(2:end) + mu_log_p), exp(p_sim(1:end-1) + mu_log_p));
            sim_moments(2,sss) = corr(exp(omg_sim(2:end) + mu_log_omg), exp(omg_sim(1:end-1) + mu_log_omg));
            sim_moments(3,sss) = corr(exp(p_sim + mu_log_p), exp(omg_sim + mu_log_omg));

        end

        moms = mean(sim_moments,2)';        
        
    end
    toc
    
