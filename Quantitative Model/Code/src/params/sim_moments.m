function moms = sim_moments(params)

    n_sim         = 1000; % number of process samples.
    n_sim_periods = 100;  % number of periods per sample
    
    mu_log_p    = params(1); 
    std_log_p   = params(2);
    p_rho       = params(3); 
    mu_log_omg  = params(4); 
    std_log_omg = params(5); 
    omg_rho     = params(6);
    corr_omg_p  = params(7);
    
    sig_p     = std_log_p   *sqrt(1 - p_rho^2);
    sig_omg   = std_log_omg *sqrt(1 - omg_rho^2);
    cov_p_omg  = corr_omg_p*sig_p*sig_omg;
    chol_mat = chol([sig_p^2, cov_p_omg; cov_p_omg, sig_omg^2]);
    
    rng('default');
    sim_moments = zeros(3,n_sim);
    for sss = 1:n_sim
        
        sim_shock_mat = randn(n_sim_periods,2)*chol_mat;

        a = [1, -p_rho];                   
        b = 1;                       
        p_sim = filter(b,a, sim_shock_mat(:,1)); 
        a = [1, -omg_rho];                   
        b = 1;                       
        omg_sim = filter(b,a, sim_shock_mat(:,2)); 

        sim_moments(1,sss) = corr(exp(p_sim(2:end) + mu_log_p), exp(p_sim(1:end-1) + mu_log_p));
        sim_moments(2,sss) = corr(exp(omg_sim(2:end) + mu_log_omg), exp(omg_sim(1:end-1) + mu_log_omg));
        sim_moments(3,sss) = corr(exp(p_sim + mu_log_p), exp(omg_sim + mu_log_omg));
            
    end
    
    moms = mean(sim_moments,2)';

end