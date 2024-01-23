! -------------------------------------------------------------------------
! mod_param.f90: declares parameters and global variables
! - init_setup: initial setup of model parameters, grids, quadrature shocks, ...
! - grid_setup: Smolyak grid setup
! -------------------------------------------------------------------------
! authors:         Rohan Kekre and Moritz Lenel
! for updates see: https://github.com/KekreLenel/safety
! -------------------------------------------------------------------------

        
module mod_param

use base_lib
use omp_lib

implicit none
private

public  ::  init_setup, grid_setup, aalpha, bbeta_vec, ies_vec, ddelta, tht_trgt_vec, sigma, varsigma_vec,  &  
            zeta, l_target, n_I, theta_target, chi, chi0_vec, results_folder, k_dev_param, w_dev_param, &
            k_grid_dev,  wh_grid_dev, wf_grid_dev, k_grid_mean, wh_grid_mean, wf_grid_mean, n_states,       &
            idx_k, idx_thh, idx_zf, idx_wh, idx_wf, n_quad, state_grid,  sidx_z, sidx_zf, next_zf_mat,   & 
            smolyak_d, smol_polynom, max_smol_level, smol_elem_ani, tht_h_grid_mean, tht_h_grid_dev,        & 
            zf_grid_mean, zf_grid_dev, dz_vec, inv_share_h,  dz_vec_adj, gma_vec, v_normalization_vec,      & 
            quad_weight_vec, shock_grid, n_interp, w_grid_adj, low_guess_fixed, high_guess_fixed, chiX,     & 
            k_grid_adj,  wealth_share_grid, vareps_w, phi_w, n_shocks, n_bond, phi_h, phi_f, tayl_ic_h,     & 
            tayl_ic_f, no_shock_idx, n_active_dims, vector_mus_dimensions,  next_dis_mat, idx_dis,          &
            dis_grid_dev, dis_grid_mean, smol_grid, rho_i, ih_grid_mean, if_grid_mean, ih_grid_dev,         & 
            if_grid_dev, idx_ih, idx_if, disast_std, disast_p, next_omg_mat, idx_omg, omg_grid_dev,         & 
            omg_grid_mean, sidx_p, irf_indices, irf_shock_sizes, sidx_dis, sidx_omg, sig_z, sig_omg,        & 
            sig_dis, safe_low_guess_fixed, safe_high_guess_fixed, foreign_trading, b_lmbd, bg_yss,          &  
            irf_indices_ll, irf_ll_shock_sizes, sidx_bg, irf_ll_length, safe_share_base, run_bg,            & 
            omg_shift, n_uni_quad_p, n_uni_quad_omg, n_uni_quad_z, n_uni_quad_zf, omg_rho,                 & 
            bbeta_coeff, n_sample, n_sims, n_irf_periods, n_sim_periods, n_burn, n_sample_periods, phi_yh,  & 
            phi_yf, run_samp, n_valuation, max_iter, max_iter_irf, irf_ll_hold_fixed


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! states: 1. aggr_capital, 2. wealth_share_h, 3. relative_tfp, 4. wage_h, 
    !         5. wage_f 6. disaster_state 7. nom_ih_last 8. nom_if_last 9. omg 
    integer, parameter  ::  smolyak_d = 9, n_I = 2, n_bond = 20, n_interp = 52+6*12
    integer, parameter  ::  n_sample = 100, n_sims = 100, n_irf_periods = 200, n_sim_periods = 400, n_burn = 10000, & 
                            n_sample_periods = 94, n_valuation = 500, max_iter = 5000, max_iter_irf = 500

    integer, save       ::  foreign_trading, vector_mus_dimensions(smolyak_d)    

    ! which shoks to do IRFs for
    integer, parameter  ::  sidx_z = 1, sidx_zf = 2, sidx_p = 3, sidx_omg = 4, sidx_dis = 5, n_irf = 6,   & 
                            sidx_bg = 1, n_irf_ll = 3, & 
                            irf_indices(n_irf) = [0, 1, 2, 3, 4, 5], irf_indices_ll(n_irf_ll) = [1, 1, 6], & 
                            irf_ll_length = 40, irf_ll_hold_fixed(n_irf_ll) = [0, 1, 1] 
    real(dp), save      ::  irf_shock_sizes(n_irf), irf_ll_shock_sizes(n_irf_ll)

    integer, parameter  ::  n_shocks  = 4
    integer, save       ::  n_uni_quad_z, n_uni_quad_zf, n_uni_quad_p, n_uni_quad_omg, max_smol_level 

    integer, save       ::  n_active_dims
    integer, save       ::  n_states, no_shock_idx(1)

    ! PARAMETERS
    real(dp), save      ::  bbeta_h, bbeta_f, gma_h, gma_f, ies_h, ies_f, chi, ddelta, aalpha, sig_zf, rho_zf, zf_std,   & 
                            sig_z, phi_h, phi_f, tayl_ic_h, tayl_ic_f, disast_p, disast_p_in, disast_shock, disast_rho,     & 
                            disast_std, disast_std_in, sig_dis, sig_omg, vareps_w, phi_w, chiX, tht_trgt_vec(n_I),          & 
                            omg_shift, bbeta_coeff, sigma, zeta, l_target, theta_target, chi0_vec(n_I),               & 
                            v_normalization_vec(n_I), safe_low_guess_fixed, safe_high_guess_fixed, low_guess_fixed,         & 
                            high_guess_fixed, inv_share_h, k_dev_param, w_dev_param, w_grid_adj, k_grid_adj, rho_i,         & 
                            ih_grid_mean, if_grid_mean, ih_grid_dev, if_grid_dev, omg_std, omg_mean, omg_rho,         & 
                            safe_share_base, b_lmbd, bg_yss, phi_yh, phi_yf 

    ! grid indices
    integer, parameter ::   idx_k = 1, idx_thh = 2, idx_zf = 3, idx_wh = 4, idx_wf = 5,     &
                            idx_dis = 6, idx_ih = 7, idx_if = 8, idx_omg = 9 

    integer, save       ::  n_GH_points, run_bg, run_samp,  n_quad, n_omg_nodes, n_dis_nodes

    real(dp), save      ::  bbeta_vec(n_I), gma_vec(n_I), ies_vec(n_I), varsigma_vec(n_I)


    real(dp), allocatable, save :: smol_grid(:,:), state_grid(:,:), smol_polynom(:,:), wealth_share_grid(:,:)
    integer, allocatable, save  ::  smolyak_elem_iso(:,:), smol_elem_ani(:,:)

    real(dp), save :: k_grid_dev, tht_h_grid_dev, zf_grid_dev, wh_grid_dev, wf_grid_dev, &   
                      k_grid_mean, tht_h_grid_mean, zf_grid_mean, wh_grid_mean, wf_grid_mean, &
                      dis_grid_dev, dis_grid_mean, omg_grid_dev, omg_grid_mean

    real(dp), allocatable, save :: shock_grid(:,:), quad_weight_vec(:), dz_vec(:), dz_vec_adj(:)
    real(dp), allocatable, save :: next_zf_mat(:,:),  next_dis_mat(:,:), next_omg_mat(:,:) 

    character, save ::    results_folder*100                      

    contains

    ! --------------------------------------------------------------------------
    ! init_setup:     creates variables needed
    !                   called once
    ! --------------------------------------------------------------------------
    subroutine init_setup

        use mod_smolyak
        use mod_normal
        use mod_markov
        ! use mod_sparsequad

        integer                         ::  ddd, ggg, zzz, hhh, fff, iii, www, eee

        ! parameter input
        integer, parameter              ::  n_prm = 65
        integer                         ::  command_len, narg
        character                       ::  command_input*2, param_file*100, n_prm_str*2
        character(12), dimension(n_prm) ::  paramnames
        real(dp), save                 ::  param_input_vec(n_prm), corr_omg_dis, cov_mat(n_shocks, n_shocks), &
                                            quad_vec_temp(n_shocks) 

        real(dp), allocatable, save     ::  uni_weight_vec_omg(:), quad_vec_omg(:), uni_weight_vec_z(:), & 
                                            uni_weight_vec_zf(:), uni_weight_vec_dis(:), quad_vec_z(:), & 
                                            quad_vec_zf(:), quad_vec_dis(:)

        integer, allocatable            ::  mus_dimensions_redux(:)

        integer                         ::  counter, info


        ! read command line input
        narg=command_argument_count()
        if (narg>0) then
            call get_command_argument(1,command_input,command_len)
            param_file = '../src/params/param_file_' // command_input(1:command_len) // '.csv'
            results_folder = '../output/tmp/res_' // command_input(1:command_len) // '/'
        else
            write(*,*) 'ERROR: No run index specified.'
            stop
        endif
        
        ! open file
        open (unit = 10, file = param_file, ACTION="read",  &
                & FORM="formatted", ACCESS="sequential")


        ! read parameters
        read(10,*)            paramnames
        write(n_prm_str,'(i2)') n_prm
        read(10,'(' // n_prm_str // 'f100.0)') param_input_vec
        close(10)
        
        ! relative country size foreign to home 
        zeta         = param_input_vec(1)

        ! preferences 
        bbeta_vec(1) = param_input_vec(2)
        bbeta_vec(2) = param_input_vec(3)

        gma_vec(1)   = param_input_vec(4)
        gma_vec(2)   = param_input_vec(5)

        ies_vec(1)   = param_input_vec(6)
        ies_vec(2)   = param_input_vec(7)
        
        chi          = param_input_vec(8)

        sigma           = param_input_vec(9)
        varsigma_vec(1) = param_input_vec(10)
        varsigma_vec(2) = param_input_vec(11)

        l_target     = param_input_vec(12)
        
        ! production
        ddelta       = param_input_vec(13)
        aalpha       = param_input_vec(14)
        chiX         = param_input_vec(15)
        inv_share_h  = param_input_vec(16)

        ! wage stickiness
        phi_w        = param_input_vec(17)
        vareps_w     = param_input_vec(18)
    
        ! monetary policy 
        phi_h        = param_input_vec(19)
        phi_f        = param_input_vec(20)
        phi_yh       = param_input_vec(21)
        phi_yf       = param_input_vec(22)
        tayl_ic_h    = param_input_vec(23)
        tayl_ic_f    = param_input_vec(24)
        rho_i        = param_input_vec(25)

        ! fiscal policy  
        b_lmbd       = param_input_vec(26)
        bg_yss       = param_input_vec(27)
        
        ! stochastic processes 
        sig_z        = maxval([sqrt_eps,param_input_vec(28)])

        zf_std      = maxval([sqrt_eps,param_input_vec(29)])
        rho_zf      = param_input_vec(30)
        
        disast_p_in  = param_input_vec(31)
        disast_shock = param_input_vec(32)
        disast_rho   = param_input_vec(33)
        disast_std_in= maxval([sqrt_eps,param_input_vec(34)])
        disast_p     = exp(disast_p_in + (disast_std_in**2)/2)
        disast_std   = sqrt((exp(disast_std_in**2) - 1)*exp(2*disast_p_in + disast_std_in**2))


        omg_mean   = param_input_vec(35)
        omg_std    = param_input_vec(36)
        omg_shift  = param_input_vec(37)
        omg_rho    = param_input_vec(38)

        corr_omg_dis = param_input_vec(39)


        tht_trgt_vec(1) = param_input_vec(40)
        tht_trgt_vec(2) = 1.0_dp - tht_trgt_vec(1)
        tht_h_grid_dev  = param_input_vec(41);
        k_dev_param     = param_input_vec(42); 
        k_grid_adj   = param_input_vec(43)
        w_dev_param     = param_input_vec(44); 
        w_grid_adj   = param_input_vec(45);
        ih_grid_mean = param_input_vec(46);
        if_grid_mean = param_input_vec(47);
        ih_grid_dev  = param_input_vec(48);
        if_grid_dev  = param_input_vec(49);

        vector_mus_dimensions(1)  = int(param_input_vec(50));
        vector_mus_dimensions(2)  = int(param_input_vec(51));
        vector_mus_dimensions(3)  = int(param_input_vec(52));
        vector_mus_dimensions(4)  = int(param_input_vec(53));
        vector_mus_dimensions(5)  = int(param_input_vec(54));
        vector_mus_dimensions(6)  = int(param_input_vec(55));
        vector_mus_dimensions(7)  = int(param_input_vec(56))
        vector_mus_dimensions(8)  = int(param_input_vec(57))
        vector_mus_dimensions(9)  = int(param_input_vec(58))
        max_smol_level = maxval(vector_mus_dimensions)
        
        n_uni_quad_z   = int(param_input_vec(59))
        n_uni_quad_zf = int(param_input_vec(60))
        n_uni_quad_p   = int(param_input_vec(61))
        n_uni_quad_omg = int(param_input_vec(62))

        foreign_trading  = int(param_input_vec(63))
        run_bg           = int(param_input_vec(64))
        run_samp         = int(param_input_vec(65))

        bbeta_coeff   = 0.001
        
        sig_dis = disast_std_in*sqrt(1.0-disast_rho**2)
        sig_omg = omg_std    *sqrt(1.0-omg_rho**2)
        sig_zf = zf_std      *sqrt(1.0-rho_zf**2)
        
        if (omg_std == 0.0 .and. corr_omg_dis == 1.0_dp) then
                sig_omg = 1.0_dp
        endif
        
        irf_shock_sizes = [0.0_dp, -5.0_dp*maxval([sqrt_eps,sig_z]), -2.0_dp*maxval([sqrt_eps,sig_zf]),  & 
                           2*maxval([sqrt_eps,sig_dis]),     2.0*sig_omg,    0.0_dp] 
        
        irf_ll_shock_sizes = [-1.0*sig_omg, -1.0*sig_omg, 0.0_dp]
        
        if (corr_omg_dis > 1.0-sqrt_eps) then 
            n_uni_quad_omg = 1
        endif
        
        n_GH_points = n_uni_quad_z*n_uni_quad_zf*n_uni_quad_p*n_uni_quad_omg
        
        n_quad = n_GH_points + 2
       
        allocate(uni_weight_vec_z(n_uni_quad_z), quad_vec_z(n_uni_quad_z), & 
                 uni_weight_vec_zf(n_uni_quad_zf), quad_vec_zf(n_uni_quad_zf), &
                 uni_weight_vec_dis(n_uni_quad_p), quad_vec_dis(n_uni_quad_p), &
                 uni_weight_vec_omg(n_uni_quad_omg), quad_vec_omg(n_uni_quad_omg))

        allocate(shock_grid(n_quad,n_shocks), quad_weight_vec(n_GH_points+1), dz_vec(n_quad), dz_vec_adj(n_quad))

        low_guess_fixed       = -10.0_dp
        high_guess_fixed      =   5.0_dp
        safe_low_guess_fixed  = -10.0_dp
        safe_high_guess_fixed =  10.0_dp
        call print_paramters()
      

        ! Covariance matrix for shocks
        !  sidx_z = 1, sidx_zf = 2, sidx_p = 3, sidx_omg = 4, & 
        cov_mat = 0.0_dp
        cov_mat(1,1) = sig_z**2;   cov_mat(2,2) = sig_zf**2
        cov_mat(3,3) = sig_dis**2; cov_mat(4,4) = sig_omg**2

        ! ONLY STORE UPPER TRIANGULAR 
        cov_mat(1,2) = 0.0   
        cov_mat(1,3) = 0.0   
        cov_mat(1,4) = 0.0  
        
        if (corr_omg_dis < 1.0_dp - sqrt_eps .and. omg_std > sqrt_eps) then
        cov_mat(3,4) = sig_omg*sig_dis*corr_omg_dis 
        else
        cov_mat(3,4) = 0.0_dp
        endif

        if (sig_z == 0.0) then
                cov_mat(1,1) = 1.0;
        endif
        if(sig_zf == 0.0) then
                cov_mat(2,2) = 1.0;
        endif
        if(sig_dis == 0.0) then
                cov_mat(3,3) = 1.0;
        endif

        call F07FDF('U', n_shocks, cov_mat, n_shocks, INFO)
        if (info /= 0) then 
            write(*,*) 'ERROR IN FACTOR DECOMPOSITION', info
            stop
        endif
        if (omg_std == 0.0 .and. corr_omg_dis == 1.0_dp) then
                sig_omg = 0.0_dp
                cov_mat(4,4) = 0.0_dp
        endif
       
       
        write(*,*) 
        write(*,*) 'CREATE SMOLYAK GRID'
        write(*,*) '---------------------------------------------'
        write(*,*) 
        ! adjust grid dimensions based on parameters
        if (rho_i > sqrt_eps) then 
            vector_mus_dimensions(idx_ih) = 3       
            vector_mus_dimensions(idx_if) = 3       
        endif

        n_active_dims = count(vector_mus_dimensions > 0)
        write(*,"(A19,I2,A1,I2)") ' Active dimensions:',n_active_dims,'/',smolyak_d
        
        allocate(mus_dimensions_redux(n_active_dims))
        counter = 0
        do ddd = 1,smolyak_d 
        if (vector_mus_dimensions(ddd) > 0) then
            counter = counter + 1
            mus_dimensions_redux(counter) = vector_mus_dimensions(ddd)
        endif
        enddo


        ! create grid 
        smolyak_elem_iso = Smolyak_Elem_Isotrop(n_active_dims, maxval(mus_dimensions_redux))
        smol_elem_ani    = Smolyak_Elem_Anisotrop(smolyak_elem_iso, n_active_dims, mus_dimensions_redux)
        smol_grid        = Smolyak_Grid(n_active_dims,maxval(mus_dimensions_redux), smol_elem_ani)
        smol_polynom     = Smolyak_Polynomial(smol_grid,n_active_dims,maxval(mus_dimensions_redux), smol_elem_ani); 
        n_states         = size(smol_grid,1)

        write(*,"(A19,I5)") ' Number of states:  ', n_states
        write(*,*)

        ! allocate wealth share grid
        allocate(wealth_share_grid(n_I,n_states), state_grid(smolyak_d, n_states))

        
        ! prepare quadrature nodes
        if (n_uni_quad_z > 1) then
            call get_quadrature_points(0.0_dp, 1.0_dp, n_uni_quad_z, quad_vec_z, uni_weight_vec_z)
        else
            quad_vec_z = 0.0_dp
            uni_weight_vec_z = 1.0_dp
        endif

        ! prepare quadrature nodes
        if (n_uni_quad_zf > 1) then
            call get_quadrature_points(0.0_dp, 1.0_dp, n_uni_quad_zf, quad_vec_zf, uni_weight_vec_zf)
        else
            quad_vec_zf       = 0.0_dp
            uni_weight_vec_zf = 1.0_dp
        endif
        
        ! prepare quadrature nodes
        if (n_uni_quad_p > 1) then
            call get_quadrature_points(0.0_dp, 1.0_dp, n_uni_quad_p, quad_vec_dis, uni_weight_vec_dis)
        else
            quad_vec_dis       = 0.0_dp
            uni_weight_vec_dis = 1.0_dp
        endif
        
        if (n_uni_quad_omg > 1) then
            call get_quadrature_points(0.0_dp, 1.0_dp, n_uni_quad_omg, quad_vec_omg, uni_weight_vec_omg)
        else
            quad_vec_omg       = 0.0_dp
            uni_weight_vec_omg = 1.0_dp
        endif

       

        counter = 0
        !DIR$ NOUNROLL
        do ggg = 1, n_uni_quad_z
        do zzz = 1, n_uni_quad_zf
        do ddd = 1, n_uni_quad_p
        do www = 1, n_uni_quad_omg
            counter = counter + 1
            
            if (corr_omg_dis > 1.0-sqrt_eps) then 
            quad_vec_temp  = [ quad_vec_z(ggg), quad_vec_zf(zzz), quad_vec_dis(ddd), quad_vec_dis(ddd)]
            else
            quad_vec_temp  = [ quad_vec_z(ggg), quad_vec_zf(zzz), quad_vec_dis(ddd), quad_vec_omg(www)]
            endif
            
            do iii = 1,n_shocks
            shock_grid( counter, iii) = sum(cov_mat(iii,:)*quad_vec_temp) 
            enddo

            quad_weight_vec(counter) =  uni_weight_vec_z(ggg)*uni_weight_vec_zf(zzz)*  & 
                                        uni_weight_vec_dis(ddd)*uni_weight_vec_omg(www)


        end do
        end do
        end do
        end do


        ! no shock node 
        quad_weight_vec(n_GH_points + 1) = 0.0_dp
        shock_grid(n_quad-1,:)           = 0.0_dp
        
        ! disast shock node 
        shock_grid(n_quad, 1)            = - 1.0*disast_shock
        shock_grid(n_quad, 2:n_shocks)   = 0.0_dp
        
        no_shock_idx =  n_quad-1        
        
        ! extract dz_vec
        dz_vec             = shock_grid(:,1)
        dz_vec_adj         = dz_vec
        dz_vec_adj(n_quad) = 0.0_dp

        ! smolyak grid deviation
        zf_grid_dev    = maxval([sqrt_eps, 3*zf_std])
        dis_grid_dev    = 3*disast_std_in
        omg_grid_dev    = maxval([sqrt_eps, 3*omg_std])

        tht_h_grid_mean  = tht_trgt_vec(1)
        zf_grid_mean  = 0.0 
        dis_grid_mean = disast_p_in
        omg_grid_mean = omg_mean

        allocate( next_zf_mat(n_quad, n_states), next_dis_mat(n_quad, n_states), next_omg_mat(n_quad, n_states))
      

        write(*,*) 
        write(*,*) 'FINISHED SETUP'
        write(*,*) '---------------------------------------------'
        write(*,*) 

    end subroutine init_setup


    subroutine grid_setup()
    integer  :: iii, counter, sss 
    
    ! ! Store grid indices

    counter = 0
    
    if (vector_mus_dimensions(1)>0) then 
    counter = counter + 1
    state_grid(1,:) = smol_grid(:,counter)*k_grid_dev      + k_grid_mean
    else 
    state_grid(1,:) = k_grid_mean
    endif

    if (vector_mus_dimensions(2)>0) then 
    counter = counter + 1
    state_grid(2,:) = smol_grid(:,counter)*tht_h_grid_dev  + tht_h_grid_mean
    else
    state_grid(2,:) = tht_h_grid_mean
    endif

    if (vector_mus_dimensions(3)>0) then 
    counter = counter + 1
    state_grid(3,:) = smol_grid(:,counter)*zf_grid_dev  + zf_grid_mean
    else 
    state_grid(3,:) = zf_grid_mean
    endif
    
    if (vector_mus_dimensions(4)>0) then 
    counter = counter + 1
    state_grid(4,:) = smol_grid(:,counter)*wh_grid_dev   + wh_grid_mean
    else 
    state_grid(4,:) = wh_grid_mean
    endif


    if (vector_mus_dimensions(5)>0) then 
    counter = counter + 1
    state_grid(5,:) = smol_grid(:,counter)*wf_grid_dev   + wf_grid_mean
    else 
    state_grid(5,:) = wf_grid_mean
    endif

    if (vector_mus_dimensions(6)>0) then 
    counter = counter + 1
    state_grid(6,:) = smol_grid(:,counter)*dis_grid_dev   + dis_grid_mean
    else 
    state_grid(6,:) = dis_grid_mean
    endif

    if (vector_mus_dimensions(7)>0) then 
    counter = counter + 1
    state_grid(7,:) = smol_grid(:,counter)*ih_grid_dev   + ih_grid_mean
    else 
    state_grid(7,:) = ih_grid_mean
    endif
    
    if (vector_mus_dimensions(8)>0) then 
    counter = counter + 1
    state_grid(8,:) = smol_grid(:,counter)*if_grid_dev   + if_grid_mean
    else 
    state_grid(8,:) = if_grid_mean
    endif
    
    if (vector_mus_dimensions(9)>0) then 
    counter = counter + 1
    state_grid(9,:) = smol_grid(:,counter)*omg_grid_dev   + omg_grid_mean
    else 
    state_grid(9,:) = omg_grid_mean
    endif
    

    wealth_share_grid(1,:) = state_grid(idx_thh,:)
    wealth_share_grid(2,:) = (1.0-state_grid(idx_thh,:))
    
 
    do sss = 1,n_states
        next_zf_mat(:,sss)  = (1.0 - rho_zf)*zf_grid_mean + rho_zf*state_grid(idx_zf,sss) &
                               + shock_grid(:,sidx_zf)
        next_dis_mat(:,sss)  = (1.0 - disast_rho)*disast_p_in + disast_rho*state_grid(idx_dis,sss) & 
                               + shock_grid(:,sidx_p)
        
        next_omg_mat(:,sss)  = (1.0 - omg_rho)*omg_grid_mean + omg_rho*state_grid(idx_omg,sss) & 
                               + shock_grid(:,sidx_omg)
        
    enddo



    end subroutine grid_setup


subroutine print_paramters()

        write(*,*) 'PARAMETRIZATION'
        write(*,*) '---------------------------------------------'
        write(*,*)
        write(*,"(A21, F10.4, F10.4)") ' bbeta        = ' , bbeta_vec
        write(*,"(A21, F10.4, F10.4)") ' gma          = ' , gma_vec
        write(*,"(A21, F10.4, F10.4)") ' ies          = ' , ies_vec
        write(*,"(A21, F10.4)")        ' chi          = ' , chi
        write(*,"(A21, F10.4)")        ' sigma        = ' , sigma
        write(*,"(A21, F10.4, F10.4)") ' varsigma     = ' , varsigma_vec
        write(*,"(A21, F10.4)")        ' ddelta       = ' , ddelta
        write(*,"(A21, F10.4)")        ' aalpha       = ' , aalpha
        write(*,"(A21, F10.4)")        ' chiX         = ' , chiX 
        write(*,"(A21, F10.4)")        ' sig_zf      = ' , sig_zf  
        write(*,"(A21, F10.4)")        ' rho_zf      = ' , rho_zf  
        write(*,"(A21, F10.4)")        ' sig_z        = ' , sig_z  
        write(*,"(A21, F10.4)")        ' phi_h        = ' , phi_h  
        write(*,"(A21, F10.4)")        ' phi_f        = ' , phi_f  
        write(*,"(A21, F10.4)")        ' tayl_ic_h    = ' , tayl_ic_h  
        write(*,"(A21, F10.4)")        ' tayl_ic_f    = ' , tayl_ic_f  
        write(*,"(A21, F10.4)")        ' disast_p     = ' , disast_p  
        write(*,"(A21, F10.4)")        ' disast_shock = ' , disast_shock  
        write(*,"(A21, F10.4)")        ' disast_rho   = ' , disast_rho  
        write(*,"(A21, F10.4)")        ' disast_std   = ' , disast_std  
        write(*,"(A21, F10.4)")        ' vareps_w     = ' , vareps_w  
        write(*,"(A21, F10.4)")        ' phi_w        = ' , phi_w  
        write(*,"(A21, F10.4)")        ' sigma        = ' , sigma        
        write(*,"(A21, F10.4)")        ' zeta         = ' , zeta
        write(*,"(A21, F10.4)")        ' l_target     = ' , l_target     
        write(*,"(A21, F10.4, F10.4)") ' tht_target   = ' , tht_trgt_vec
        write(*,*)
        write(*,*) '---------------------------------------------'
        write(*,*)

end subroutine print_paramters


end module mod_param


