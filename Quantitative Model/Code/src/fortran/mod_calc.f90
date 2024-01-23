! -------------------------------------------------------------------------
! mod_calc.f90: solution module, public subroutines:
! - calc_steady: solves for steady state to center smolyak grid 
! - calc_sol: solves for stationary equilibrium and creates impulse respones
! -------------------------------------------------------------------------
! authors:         Rohan Kekre and Moritz Lenel
! for updates see: https://github.com/KekreLenel/safety
! -------------------------------------------------------------------------

module mod_calc

use omp_lib
use base_lib,  only: dp, sqrt_eps, eps, initBrent, updateBrent
use mod_param, only: n_I

implicit none

private

public  :: calc_steady, calc_sol 

real(dp) :: v_ss_vec(n_I), mc_ss_vec(n_I), k_ss, s_ss, rf_ss, l_ss_vec(2), q_ss, & 
            tot_wealth_ss_vec(n_I), q_l_ss_vec(n_I), c_cost_ss_vec(n_I), w_ss_vec(2)
real(dp), parameter :: conv_crit = 1E-08_dp ! numerical concergence criterion

contains

subroutine calc_steady()
! calculates steady state values, prepares state grid construction

    use mod_param, only :   bbeta_vec, IES_vec, ddelta, aalpha, tht_trgt_vec, inv_share_h, sigma, &
                            varsigma_vec, zeta, l_target, n_I, chi, chi0_vec, results_folder, &
                            k_grid_dev, wh_grid_dev, wf_grid_dev, k_grid_mean, wh_grid_mean, &
                            v_normalization_vec, k_dev_param, w_dev_param, w_grid_adj, k_grid_adj, & 
                            wf_grid_mean, sig_dis, sig_omg, omg_grid_mean, dis_grid_mean

    real(dp) :: bbeta_avrg, ies_avrg, rk_ss, kappa_h_ss, kappa_f_ss,  & 
                c_h_aggr_ss, c_f_aggr_ss, c_h_ss_vec(n_I), aggr_css, y_h_ss, y_f_ss, &
                c_f_ss_vec(n_I), theta_ss_vec(n_I), zf_ss(2), pi_ss, aggr_wealth, diff, & 
                s_ss_new, chi0, chi0_new, c_bundle, c_bundle_ss_vec(n_I), inv_h, inv_f

    integer  :: iter, iii
    
    ! average productivity level
    zf_ss = [1.0_dp, zeta]

    ! average discount factor and IES
    bbeta_avrg = sum(tht_trgt_vec*bbeta_vec) 
    IES_avrg   = sum(tht_trgt_vec*IES_vec)

    ! steady state returns
    rk_ss = 1.0_dp/bbeta_avrg-1.0_dp
    rf_ss = 1.0_dp/bbeta_avrg-1.0_dp
    
    ! labor in steady state (targeted) 
    l_ss_vec = l_target

    ! initalize home good consumption of foreign
    s_ss = 1.0_dp    
    
    ! wealth shares (targeted)
    theta_ss_vec = tht_trgt_vec 
    

    diff = 1.0_dp; iter = 0
    do while (diff > sqrt_eps)
    
        iter = iter+1

        q_ss  =inv_share_h**(1.0/(1.0-sigma)) * (1 + &
              (1.0-inv_share_h)/inv_share_h*s_ss**(sigma-1.0))**(1.0/(1.0-sigma))         
        

        pi_ss = ( rk_ss + ddelta ) * q_ss 

        kappa_h_ss = zf_ss(1)*l_ss_vec(1)*(pi_ss/aalpha)**(1.0_dp/(aalpha-1.0_dp))
        kappa_f_ss = zf_ss(2)*l_ss_vec(2)*(s_ss*pi_ss/aalpha)**(1.0_dp/(aalpha-1.0_dp))

        k_ss = kappa_h_ss + kappa_f_ss

        inv_h = ddelta*k_ss*q_ss/(1.0+(s_ss**(sigma-1.0))*(1.0-inv_share_h)/(inv_share_h))
        inv_f = inv_h * s_ss**(sigma)*(1.0-inv_share_h)/inv_share_h 

        ! aggregate output
        y_h_ss  =  (kappa_h_ss**aalpha)*( (zf_ss(1)*l_ss_vec(1))**(1.0_dp-aalpha) )
        y_f_ss  =  (kappa_f_ss**aalpha)*( (zf_ss(2)*l_ss_vec(2))**(1.0_dp-aalpha) )
     
        c_h_aggr_ss = (zf_ss(1)*l_ss_vec(1))**(1.0-aalpha) * kappa_h_ss**aalpha - inv_h 
        c_f_aggr_ss = (zf_ss(2)*l_ss_vec(2))**(1.0-aalpha) * kappa_f_ss**aalpha - inv_f 

        aggr_wealth = (pi_ss + (1.0-ddelta)*q_ss)*(kappa_f_ss+kappa_h_ss) 

        ! find real wages, expressed in home goods
        w_ss_vec(1)   = (1.0_dp - aalpha)*(kappa_h_ss**aalpha)*l_ss_vec(1)**(-aalpha) & 
                        *(zf_ss(1)**(1.0_dp-aalpha))
        w_ss_vec(2)   = (1.0_dp - aalpha)*(kappa_f_ss**aalpha)*l_ss_vec(2)**(-aalpha) & 
                        *s_ss**(-1)*(zf_ss(2)**(1.0_dp-aalpha))


        do iii = 1,n_I 
            c_h_ss_vec(iii) = ( theta_ss_vec(iii)*q_ss*(kappa_f_ss + kappa_h_ss)*rf_ss + w_ss_vec(iii) *l_ss_vec(iii) ) &
                              /(1.0+(s_ss**(sigma-1.0))*(1.0-varsigma_vec(iii))/(varsigma_vec(iii)))
        enddo
         
        c_f_ss_vec    = c_h_ss_vec*s_ss**(sigma)*(1.0-varsigma_vec)/(varsigma_vec)

        s_ss_new      = ( c_f_aggr_ss / sum(c_h_ss_vec*(1.0-varsigma_vec)/varsigma_vec) )**(1.0_dp/sigma)

        diff  = abs(s_ss_new - s_ss) 
        
        s_ss = s_ss + 0.01*(s_ss_new-s_ss)

        if (iter>10000) then; write(*,*) 'no convergence s_ss'; stop; endif
    enddo

    ! c_cost_vec 
    c_cost_ss_vec =  varsigma_vec**(1.0/(1.0-sigma)) * (1 + &
          (1.0-varsigma_vec)/varsigma_vec*s_ss**(sigma-1))**(1.0/(1.0-sigma))         


    ! now find disutility parameter to get labor right
    do iii = 1,n_I
        c_bundle = calc_bundle(c_h_ss_vec(iii), c_f_ss_vec(iii), varsigma_vec(iii) )

            chi0     = 0.5_dp
            diff     = 1.0_dp
            iter = 0
            do while (diff>1E-12_dp)

                iter = iter + 1
               
                chi0_new = IES_vec(iii)*w_ss_vec(iii)/c_cost_ss_vec(iii) /c_bundle  * l_ss_vec(iii)**(-1.0_dp) * &
                         (1.0_dp + (1.0_dp/IES_vec(iii) - 1.0_dp) *chi0 * chi / (1.0_dp+chi)  )
               
                diff = abs(chi0-chi0_new)
                chi0 = chi0 + 0.5*(chi0_new-chi0)

                if (iter > 1000) then
                    write(*,*) 'ERROR: No convergence in finding chi0.'
                    stop
                end if

            end do


        chi0_vec(iii) = chi0


        v_normalization_vec(iii) = (1-bbeta_vec(iii)) * ( c_bundle**(1.0-1/ies_vec(iii))  *(1.0+(1.0/IES_vec(iii)-1.0)*chi0/(1.0 + 1.0/chi) )**( & 
                                        1.0/ies_vec(iii) ) )**(-1)
        v_ss_vec(iii) = ( v_normalization_vec(iii)/(1.0-bbeta_vec(iii)) )**(1.0/(1.0 - 1.0/IES_vec(iii))) * & 
                        c_bundle*(1+(1/IES_vec(iii)-1)*chi0/(1.0 + 1.0/chi))**( (1.0/ies_vec(iii))/(1.0-1.0/ies_vec(iii)) )
        mc_ss_vec(iii) = ((1+(1/IES_vec(iii)-1)*chi0/(1.0 + 1.0/chi)) / c_bundle)**(1/IES_vec(iii))
        c_bundle_ss_vec(iii) = c_bundle

        
    enddo

    q_l_ss_vec = w_ss_vec/(1.0-bbeta_vec)

    tot_wealth_ss_vec = q_l_ss_vec + theta_ss_vec*k_ss*(1.0+rk_ss)

    write(*,*) 'STEADY STATE'
    write(*,*) '---------------------------------------------'
    write(*,"(A42)")                   '                        HOME       FOREIGN'
    write(*,*)
    write(*,"(A21, F10.4, A1, F10.4)") ' k_ss   =            ' , kappa_h_ss,    ' ', kappa_f_ss
    write(*,"(A21, F10.4, A1, F10.4)") ' y_ss   =            ' , y_h_ss,        ' ', y_f_ss
    write(*,"(A21, F10.4, A1, F10.4)") ' w_ss   =            ' , w_ss_vec(1),   ' ', w_ss_vec(2)
    write(*,"(A21, F10.4, A1, F10.4)") ' l_ss   =            ' , l_ss_vec(1),   ' ', l_ss_vec(2)
    write(*,"(A21, F10.4)")            ' x_ss   =            ' , k_ss * ddelta
    write(*,"(A21, F10.4)")            ' s_ss   =            ' , s_ss 
    write(*,"(A21, F10.4, A1, F10.4)") ' q_l_ss =            ' , q_l_ss_vec(1), ' ', q_l_ss_vec(2)
    write(*,"(A21, F10.4, A1, F10.4)") ' chi0   =            ' , chi0_vec(1), ' ', chi0_vec(2)
    write(*,"(A21, F10.4, A1, F10.4)") ' c_h_ss =            ' , c_h_ss_vec(1), ' ', c_h_ss_vec(2)
    write(*,"(A21, F10.4, A1, F10.4)") ' c_f_ss =            ' , c_f_ss_vec(1), ' ', c_f_ss_vec(2)
    write(*,"(A21, F10.4, A1, F10.4)") ' v_ss   =            ' , v_ss_vec(1), ' ', v_ss_vec(2)
    write(*,*)
    

    k_grid_mean   = k_grid_adj*(kappa_h_ss + kappa_f_ss)
    wh_grid_mean  = w_grid_adj*w_ss_vec(1)
    wf_grid_mean  = w_grid_adj*w_ss_vec(2)

    k_grid_dev   = k_dev_param*k_grid_mean
    wh_grid_dev  = w_dev_param*wh_grid_mean
    wf_grid_dev  = w_dev_param*wf_grid_mean
 

    aggr_css = sum(c_cost_ss_vec*c_bundle_ss_vec)
    

    open (unit = 10, file = trim(results_folder) // 'extra_data.csv', ACTION="write",  &
            & FORM="formatted", ACCESS="sequential")

    write(10,'(8F10.4)') chi0_vec(1), chi0_vec(2), sig_dis, sig_omg, dis_grid_mean, omg_grid_mean, y_h_ss, aggr_css

    close(10)

end subroutine calc_steady



subroutine calc_sol() 
! outer solution routine

use mod_smolyak, only : Smolyak_Polynomial2
use mod_param, only:    n_states, idx_k, no_shock_idx, shock_grid, irf_indices_ll, idx_zf, & 
                        idx_wf, inv_share_h, ddelta, idx_omg, idx_dis, state_grid, n_quad,  &
                        aalpha, smolyak_d, n_active_dims, next_zf_mat, smol_polynom,       & 
                        disast_p, max_smol_level, smol_elem_ani, k_grid_mean, k_grid_dev,   & 
                        tht_h_grid_mean,  tht_h_grid_dev, wh_grid_mean, wh_grid_dev,        & 
                        wf_grid_mean, wf_grid_dev, zf_grid_mean, zf_grid_dev, zeta,   & 
                        n_shocks, dz_vec, dz_vec_adj, bbeta_vec, wealth_share_grid,           & 
                        varsigma_vec, sigma, results_folder, n_interp, run_bg,              & 
                        vector_mus_dimensions, next_dis_mat, dis_grid_mean, dis_grid_dev,   & 
                        n_bond, ih_grid_mean, ih_grid_dev, if_grid_mean, if_grid_dev,       & 
                        omg_shift, next_omg_mat, omg_grid_mean, omg_grid_dev, irf_indices,  &
                        irf_shock_sizes, sig_omg, sig_dis, sidx_zf, sidx_z,         & 
                        sidx_dis, sidx_p, sidx_omg, sidx_bg,                                &  
                        irf_ll_hold_fixed, irf_ll_shock_sizes, irf_ll_length, b_lmbd,       & 
                        omg_rho, quad_weight_vec, n_valuation, max_iter, max_iter_irf

integer, parameter  :: n_nxt = 8 + 2*n_I   ! number of imterpolated next period variables

real(dp) :: v_mat(n_I, n_states), mc_mat(n_I, n_states), k_next_mat(n_quad, n_states),  q_mat(n_states),        & 
            kappa_mat(2,n_states), l_aggr_mat(2, n_states), next_state_mat(n_quad, n_active_dims, n_states),    &
            next_state_mat_new(n_quad, n_active_dims, n_states), nxt_mat(n_quad, n_nxt), nxt_mat_2(n_quad,3),   & 
            nxt_bond_prices(n_quad, 3), nxt_mat_trans(1, n_nxt), nxt_mat_2_trans(1,3), p_dis,                   & 
            nxt_valuation(n_quad, 7), interp_input_mat(n_states, n_nxt), nom_i_mat(2,n_states),                 & 
            s_mat(n_states), nom_i_vec(2), s, big_weight_vec(n_quad), w_choice_mat(2, n_states), q_current,     & 
            infl_mat(2, n_states), c_spending_vec(n_I),  c_spending_mat(n_I, n_states),                         & 
            share_mat(n_I, n_states), share_vec(n_I),  bF_share_mat(n_I, n_states), bF_share_vec(n_I),      & 
            infl_vec(2), l_aggr_vec(2), q_new, s_new, k_next_new, c_spending_vec_new(n_I), infl_new(2),         & 
            w_choice_new(2), v_new(n_I), mc_new(n_I), theta_nxt(n_quad), l_aggr_new(2),                         & 
            v_mat_new(n_I, n_states), mc_mat_new(n_I, n_states), l_aggr_mat_new(2,n_states),                    & 
            q_mat_new(n_states), nom_i_mat_new(2, n_states), inv_f_mat(n_states), s_mat_new(n_states),          & 
            c_spending_mat_new(n_I, n_states), inv_h_mat(n_states), w_choice_mat_new(2, n_states),              & 
            k_next_mat_new(n_quad, n_states), inv_f, inv_h, bF_share_mat_new(n_I, n_states),                  & 
            share_mat_new(n_I, n_states), infl_mat_new(2, n_states),                                            & 
            results_vec(n_interp-6*12), results_mat(n_states, n_interp), q_bond_mat(n_states, 3, n_bond+1),     & 
            valuation_vec(7), results_irf_mat(n_states, n_interp), q_bond_irf_mat(n_states, 3, n_bond+1),       & 
            bond_supply_shock, bond_supply_shock_nxt, valuation_mat_store(14,n_states),                         &
            q_bond_irf_mat_tmp(n_states, 3, n_bond+1), valuation_mat(n_states,7),                               & 
            valuation_mat_tmp(n_states,7), valuation_irf_mat(n_states,7)


real(dp) :: next_state_irf_mat(n_quad, n_active_dims, n_states),  nom_i_irf_mat(2,n_states),            & 
            s_irf_mat(n_states), q_irf_mat(n_states), c_spending_irf_mat(n_I,n_states),                 & 
            share_irf_mat(n_I,n_states), infl_irf_mat(2,n_states), bF_share_irf_mat(n_I,n_states),    &
            l_aggr_irf_mat(2,n_states), next_state_irf_mat_new(n_quad, n_active_dims, n_states),        & 
            M_vec_mat(n_quad, n_I, n_states), inv_h_irf_mat(n_states), inv_f_irf_mat(n_states),         &
            bF_share_irf_mat_new(n_I,n_states), share_irf_mat_new(n_I,n_states),                      &
            w_choice_irf_mat_new(2,n_states), infl_irf_mat_new(2,n_states), v_irf_mat(n_I,n_states),    & 
            mc_irf_mat(n_I, n_states), q_irf_mat_new(n_states), l_aggr_irf_mat_new(2,n_states),         & 
            v_irf_mat_new(n_I,n_states), mc_irf_mat_new(n_I,n_states), s_irf_mat_new(n_states),         & 
            nom_i_irf_mat_new(2,n_states), c_spending_irf_mat_new(n_I,n_states), E_rb(6),               &
            k_next_irf_mat_new(n_quad,n_states), zf_tmp(n_states,2), who_prices_bond(3),               & 
            next_state_irftransition_mat_new(1,n_active_dims, n_states), w_choice_irf_mat(2,n_states),  & 
            next_state_irftransition_mat(1,n_active_dims, n_states), k_next_irf_mat(n_quad, n_states),  & 
            trans_shock_vec(4)

character :: filename*100, n_shock_str*20
! interpolation variables
integer, parameter :: nrhs = n_nxt

integer  :: lda, ldaf, ldb, ldx
integer  :: ipiv(n_states), iwork(n_states)
integer  :: info
real(dp) :: af(n_states,n_states), r(n_states), & 
            c(n_states), smol_coeffs(n_states,nrhs), work(4*n_states), ferr(nrhs), &  
            berr(nrhs), polyn_points(n_quad,n_states), polyn_points_trans(1,n_states)
real(dp) ::  rcond
character(1) :: equed

character :: shock_char*1, shock_char2*2, iter_char*4

integer, parameter :: nrhs2 = 3
real(dp) :: bond_coeffs(n_states,nrhs2), ferr2(nrhs2), &  
            berr2(nrhs2)
integer, parameter :: nrhs3 = 7
real(dp) :: val_coeffs(n_states,nrhs3), ferr3(nrhs3), &  
            berr3(nrhs3)

real(dp) :: diff, timer_start, timer_end
integer  :: outer_iter, ifail
integer  :: sss, iii, counter, bbb, fff, ttt, bbb_counter
real(dp) :: next_state(n_quad, n_active_dims), excess_b, k_next(n_quad), nxt_omg(n_quad),  next_zf(n_quad), & 
            infl_h, ih_fixed, infl_A, infl_B, step_sign, step_sign_prev, excess_b_B
real(dp) :: brnt_xbracket(2), brnt_ybracket(2), brnt_xnew, brnt_ynew, brnt_xvec(4), brnt_yvec(3)
logical  :: brnt_mflag
integer  :: brnt_sflag, brnt_iter, bracket_success

    ! ==================================================================
    ! SETUP
    ! ==================================================================
    lda  = n_states; ldaf = n_states 
    ldb  = n_states; ldx  = n_states
   
    ! ==================================================================
    ! INITIALIZATION
    ! ==================================================================

    ! initialize various matrices 
    l_aggr_mat      = 1.0_dp 
    q_mat           = q_ss
    nom_i_mat(1,:)  = (1.0 + rf_ss)/(1+exp(state_grid(idx_omg,:)))
    nom_i_mat(2,:)  = 1.0 + rf_ss - nom_i_mat(1,:)        ! nom_i_mat(2,:) is difference of i* to i
    s_mat           = s_ss ! exp(state_grid(idx_eps,:)) ! s_ss
    infl_mat        = 1.0_dp
    q_bond_mat      = 1.0_dp
    share_mat       = 0.0_dp
    inv_h_mat       = ddelta*state_grid(idx_k,:)*inv_share_h
    inv_f_mat       = ddelta*state_grid(idx_k,:)*(1.0-inv_share_h)
    
    ! start from only dollar holdings
    bF_share_mat = 0.0_dp


    ! initialize value function
    do iii = 1,n_I
        v_mat(iii,:) = v_ss_vec(iii)
        mc_mat(iii,:) = mc_ss_vec(iii)
    enddo

    
    ! collect zeta*zf
    zf_tmp(:,1) = 1.0
    zf_tmp(:,2) = zeta*exp(state_grid(idx_zf,:))

    ! initialize deployed capital
    kappa_mat(1,:) = state_grid(idx_k,:) * (s_mat**(1.0/(1.0-aalpha))) * &
                     zf_tmp(:,1)/zf_tmp(:,2) * l_aggr_mat(1,:)/l_aggr_mat(2,:) / & 
                     (1.0 +   (s_mat**(1.0/(1.0-aalpha))) * &
                     zf_tmp(:,1)/zf_tmp(:,2) * l_aggr_mat(1,:)/l_aggr_mat(2,:) )
    kappa_mat(2,:) = state_grid(idx_k,:) - kappa_mat(1,:)
    
    M_vec_mat    = 0.0_dp

    ! initialize wages (expressed in home goods)
    w_choice_mat(1,:) = (1.0-aalpha)* kappa_mat(1,:)**aalpha * l_ss_vec(1)**(-aalpha)  & 
                        * ( zf_tmp(:,1)**(1.0 -aalpha) ) 
    w_choice_mat(2,:) = (1.0-aalpha)* kappa_mat(2,:)**aalpha * l_ss_vec(2)**(-aalpha) &
                        *s_ss**(-1) * ( zf_tmp(:,2)**(1.0-aalpha) )
        
        
    ! initialize consumption guess 
    ! (this is aggregate consumption spending in terms of home produced good)
    do iii = 1,n_I
            c_spending_mat(iii, :) = w_choice_mat(iii,:)*l_ss_vec(iii) + & 
                                     state_grid(idx_k,:)*wealth_share_grid(iii,:)*(1.0_dp/bbeta_vec(iii)-1.0_dp)
    enddo

    ! initialize next state transition matrix for active states
    next_state_mat = 0.0_dp
    do sss = 1, n_states
    ! states: 1. aggr_capital, 2. wealth_share_h, 3. relative_tfp,  4. wage_h, 
    !         5. wage_f 6. disaster_state 7. i_h 8. i_f 9. omg 
        
        counter = 0
        
        k_next_mat(:, sss) = state_grid(idx_k,sss)*exp(dz_vec)/exp(dz_vec_adj)

        if (vector_mus_dimensions(1) > 0) then
        counter = counter + 1
        next_state_mat(:,counter,sss)  = (state_grid(idx_k,sss)/exp(dz_vec_adj) - k_grid_mean ) /k_grid_dev
        endif

        if (vector_mus_dimensions(2) > 0) then
        counter = counter + 1
        next_state_mat(:,counter,sss)  = 0.0_dp 
        endif

        if (vector_mus_dimensions(3) > 0) then
        counter = counter + 1
        next_state_mat(:,counter,sss) = (next_zf_mat(:,sss)       - zf_grid_mean) /zf_grid_dev
        endif

        if (vector_mus_dimensions(4) > 0) then
        counter = counter + 1
        next_state_mat(:,counter,sss) = (w_choice_mat(1,sss) - wh_grid_mean) /wh_grid_dev
        endif

        if (vector_mus_dimensions(5) > 0) then
        counter = counter + 1
        next_state_mat(:,counter,sss) = (w_choice_mat(2,sss) - wf_grid_mean) /wf_grid_dev
        endif

        if (vector_mus_dimensions(6) > 0) then
        counter = counter + 1
        if (dis_grid_dev < sqrt_eps) then
        next_state_mat(:,counter,sss) = 0.0_dp 
        else
        next_state_mat(:,counter,sss) = (next_dis_mat(:,sss)   - dis_grid_mean) /dis_grid_dev
        endif
        endif

        if (vector_mus_dimensions(7) > 0) then
        counter = counter + 1
        next_state_mat(:,counter,sss) = 0.0_dp
        endif
        
        if (vector_mus_dimensions(8) > 0) then
        counter = counter + 1
        next_state_mat(:,counter,sss) = 0.0_dp 
        endif

        if (vector_mus_dimensions(9) > 0) then
            counter = counter + 1
            if (omg_grid_dev < sqrt_eps) then
                next_state_mat(:,counter,sss) = 0.0_dp 
            else
                next_state_mat(:,counter,sss) = (next_omg_mat(:,sss)   - omg_grid_mean) /omg_grid_dev
            endif
        endif

        

    enddo
    
    ! make sure transition matrix is within bounds
    where (next_state_mat > 1.0_dp)
        next_state_mat = 1.0_dp 
    elsewhere (next_state_mat < -1.0_dp)
        next_state_mat = -1.0_dp 
    endwhere

    ! Just for preallocation
    next_state_mat_new = next_state_mat

    ! ==================================================================
    ! ITERATION
    ! ==================================================================
    
    ! start loop, solve until convergence
    outer_iter = 0; diff = 1.0_dp
    do while ( diff > conv_crit .and. outer_iter < max_iter) 
    
        
        ! add to iteration counter
        outer_iter = outer_iter + 1

        ! initialize OMP timer
        timer_start = omp_get_wtime()

        ! variables that will be interpolated when 
        ! forming next period expectations:
        ! if (outer_iter < 1000) then
        do iii = 1,n_I
            interp_input_mat(:,iii)  = v_mat(iii,:)
            interp_input_mat(:,n_I + iii)  = mc_mat(iii,:)
        enddo
        interp_input_mat(:,2*n_I + 1)  = s_mat
        interp_input_mat(:,2*n_I + 2)  = q_mat
        interp_input_mat(:,2*n_I + 3)  = l_aggr_mat(1,:)
        interp_input_mat(:,2*n_I + 4)  = l_aggr_mat(2,:)
        interp_input_mat(:,2*n_I + 5)  = infl_mat(1,:)
        interp_input_mat(:,2*n_I + 6)  = infl_mat(2,:)
        interp_input_mat(:,2*n_I + 7) = c_spending_mat(1,:) 
        interp_input_mat(:,2*n_I + 8) = c_spending_mat(2,:) 

        ! endif

        ! solve for polynomial coefficients for interpolation
        call F07ABF( 'Equilibration','No transpose',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,    &
                      equed,r,c,interp_input_mat,ldb,smol_coeffs,ldx,rcond,ferr,berr,work,iwork,info )
        
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED( & 
        !$OMP next_state_mat, smol_elem_ani, smol_coeffs, state_grid, &
        !$OMP nom_i_mat, s_mat, w_choice_mat, n_states, dz_vec, next_omg_mat, & 
        !$OMP k_grid_dev, k_grid_mean, zf_grid_dev, zf_grid_mean, q_mat, &
        !$OMP c_spending_mat, share_mat, infl_mat, varsigma_vec, sigma, & 
        !$OMP bF_share_mat, l_aggr_mat, next_state_mat_new, inv_h_mat, inv_f_mat, &
        !$OMP bF_share_mat_new, share_mat_new, disast_p, vector_mus_dimensions, &
        !$OMP wh_grid_mean, wh_grid_dev, wf_grid_mean, wf_grid_dev, n_active_dims, M_vec_mat, &
        !$OMP tht_h_grid_mean, tht_h_grid_dev, mc_mat_new, &
        !$OMP w_choice_mat_new, infl_mat_new, n_quad, max_smol_level, &
        !$OMP q_mat_new, l_aggr_mat_new, v_mat_new, s_mat_new, nom_i_mat_new, outer_iter, &
        !$OMP c_spending_mat_new, k_next_mat_new, dz_vec_adj, next_zf_mat, &
        !$OMP k_next_mat, ih_grid_mean, ih_grid_dev, if_grid_mean, if_grid_dev) &
        !$OMP PRIVATE( &
        !$OMP polyn_points, nxt_mat, inv_h, inv_f, mc_new, &
        !$OMP nom_i_vec, s, nxt_mat_2, counter, &
        !$OMP q_current, c_spending_vec, share_vec, &
        !$OMP bF_share_vec, infl_vec , l_aggr_vec, &
        !$OMP q_new, s_new, k_next_new, c_spending_vec_new, infl_new, &
        !$OMP sss, v_new, l_aggr_new, w_choice_new, theta_nxt)
        !$OMP DO SCHEDULE(static)
        do sss = 1, n_states 

            ! interpolate future values 
            polyn_points = Smolyak_Polynomial2(next_state_mat(:,:,sss), n_active_dims, &
                                      n_quad, n_states, max_smol_level, smol_elem_ani)

            ! use MKL Library for multiplication, same as
            ! nxt_mat        = matmul(polyn_points, smol_coeffs)
            CALL DGEMM('N','N', n_quad, nrhs, n_states, 1.0_dp, polyn_points, n_quad, & 
                               smol_coeffs, n_states, 0.0_dp, nxt_mat, n_quad)   
            
            ! non interpolated values
            nxt_mat_2(:, 1) = k_next_mat(:,sss) 
            nxt_mat_2(:, 2) = next_zf_mat(:,sss) 
            nxt_mat_2(:, 3) = next_omg_mat(:,sss)


            ! current state and guesses 
            nom_i_vec         = nom_i_mat(:,sss)
            s                 = s_mat(sss)
            q_current         = q_mat(sss)
            c_spending_vec    = c_spending_mat(:,sss)
            share_vec         = share_mat(:,sss)
            bF_share_vec    = bF_share_mat(:,sss)
            infl_vec          = infl_mat(:,sss)
            l_aggr_vec        = l_aggr_mat(:,sss)
            
            ! calculate current state equilibrium 
            call calc_equilibrium_and_update( nxt_mat, n_nxt, nxt_mat_2, sss, outer_iter, &
                        c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                        nom_i_vec, share_vec, bF_share_vec, inv_h, inv_f,  0.0_dp, 0.0_dp, &
                        q_new, s_new, k_next_new, c_spending_vec_new, & 
                        infl_new, v_new, mc_new, l_aggr_new, w_choice_new, theta_nxt, M_vec_mat(:, :, sss))
                
            v_mat_new(:,sss)           = v_new
            mc_mat_new(:, sss)         = mc_new
            s_mat_new(sss)             = s_new
            c_spending_mat_new(:,sss)  = c_spending_vec_new
            q_mat_new(sss)             = q_new 
            inv_f_mat(sss)             = inv_f 
            inv_h_mat(sss)             = inv_h 
            l_aggr_mat_new(:,sss)      = l_aggr_new
            nom_i_mat_new(:,sss)       = nom_i_vec
            w_choice_mat_new(:,sss)    = w_choice_new
            share_mat_new(:,sss)       = share_vec          
            bF_share_mat_new(:,sss)  = bF_share_vec     
            infl_mat_new(:,sss)        = infl_new
            ! adjustment k_next: depreciation in disaster states
            k_next_mat_new(:,sss)      = k_next_new/exp(dz_vec_adj)*exp(dz_vec)

        counter = 0
        
        if (vector_mus_dimensions(1) > 0) then
            counter = counter + 1
            next_state_mat_new(:, counter, sss)   = (k_next_new/exp(dz_vec_adj) -k_grid_mean )/k_grid_dev
        endif

        if (vector_mus_dimensions(2) > 0) then
            counter = counter + 1
            next_state_mat_new(:, counter, sss) = (theta_nxt - tht_h_grid_mean)/tht_h_grid_dev
        endif

        if (vector_mus_dimensions(3) > 0) then
        counter = counter + 1
        endif
       

        if (vector_mus_dimensions(4) > 0) then
            counter = counter + 1
            next_state_mat_new(:, counter, sss)  = (w_choice_new(1)/exp(dz_vec_adj)/ & 
                                                          (varsigma_vec(1) +  & 
                                                          (1.0-varsigma_vec(1) )*s**(sigma-1.0))**(1.0/(sigma-1.0)) & 
                                                          -wh_grid_mean) /wh_grid_dev
        endif

        if (vector_mus_dimensions(5) > 0) then
            counter = counter + 1
            next_state_mat_new(:, counter, sss)  = (w_choice_new(2)/exp(dz_vec_adj)*s/ & 
                                                          (varsigma_vec(2)*s**(1.0_dp-sigma) + (1.0-varsigma_vec(2)))**(1.0/(sigma-1.0)) & 
                                                          -wf_grid_mean) /wf_grid_dev
        endif

        if (vector_mus_dimensions(6) > 0) then
        counter = counter + 1
        endif
        

        if (vector_mus_dimensions(7) > 0) then
            counter = counter + 1
            next_state_mat_new(:, counter, sss)  = ((nom_i_vec(1)-1.0) - ih_grid_mean) /ih_grid_dev 
        endif

        if (vector_mus_dimensions(8) > 0) then
            counter = counter + 1
            next_state_mat_new(:, counter, sss)  = ((nom_i_vec(1) + nom_i_vec(2)-1.0) - if_grid_mean) /if_grid_dev 
        endif
        
        if (vector_mus_dimensions(9) > 0) then
        counter = counter + 1
        endif

        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        
    ! constrain transition matrix
    where(next_state_mat_new > 1.0_dp) 
        next_state_mat_new   = 1.0_dp 
    endwhere
    where(next_state_mat_new < -1.0_dp)
        next_state_mat_new   = -1.0_dp
    endwhere

    ! check convergence
    diff =  maxval([ maxval(abs(log(v_mat_new/v_mat))), &
                     maxval(abs(log(mc_mat_new/mc_mat))), &
                     maxval(abs(log(c_spending_mat_new/c_spending_mat))), &
                     maxval(abs((share_mat_new-share_mat))), &
                     maxval(abs(log(k_next_mat_new/k_next_mat))), &
                     maxval(abs(log(l_aggr_mat_new/l_aggr_mat))), &
                     maxval(abs((nom_i_mat_new-nom_i_mat))), &
                     maxval(abs(log(w_choice_mat_new/w_choice_mat))), &
                     maxval(abs(log(q_mat_new/q_mat))), &
                     maxval(abs(log(infl_mat_new/infl_mat))), &
                     maxval(abs(next_state_mat - next_state_mat_new)) ])

    timer_end = omp_get_wtime()


    write(*,*) ''
    write(*,"(A11, i5)")    ' Iteration ', outer_iter
    write(*,"(A11, f8.2)")  ' Calc. time', timer_end-timer_start
    write(*,*) ''

    write(*,"(A19,e12.4,A1)") ' Changes (max diff: ', diff , ')'
    write(*,*) '------------------------------------------------------------------------'
   
    write(*,"(A13, e12.4, e12.4)")     ' v       = ' , &  
    maxval(abs(log(v_mat_new(1,:)/v_mat(1,:)))), & 
    maxval(abs(log(v_mat_new(2,:)/v_mat(2,:))))
    write(*,"(A13, e12.4, e12.4)")     ' mc       = ' , &  
    maxval(abs(log(mc_mat_new(1,:)/mc_mat(1,:)))), & 
    maxval(abs(log(mc_mat_new(2,:)/mc_mat(2,:))))
    write(*,"(A13, e12.4, e12.4)")     ' c_spend   = ' , &
    maxval(abs(log(c_spending_mat_new(1,:)/c_spending_mat(1,:)))), & 
    maxval(abs(log(c_spending_mat_new(2,:)/c_spending_mat(2,:))))
    write(*,"(A13, e12.4, e12.4)")     ' share_mat = ' , &
    maxval(abs((share_mat_new(1,:)-share_mat(1,:)))), & 
    maxval(abs((share_mat_new(2,:)-share_mat(2,:))))
    write(*,"(A13, e12.4, e12.4)")     ' safes_mat = ' , &
    maxval(abs((bF_share_mat_new(1,:)-bF_share_mat(1,:)))), & 
    maxval(abs((bF_share_mat_new(2,:)-bF_share_mat(2,:))))
    write(*,"(A13, e12.4)")     ' k         = ' , &
    maxval(abs(log(k_next_mat_new/k_next_mat)))
    write(*,"(A13, e12.4)")     ' s         = ' , &
    maxval(abs(log(s_mat_new/s_mat)))
    do sss = 1,n_active_dims 
    write(*,"(A13, e12.4)")     ' next_state= ' , &
    maxval(abs((next_state_mat_new(:,sss,:)-next_state_mat(:,sss,:))))
    enddo
    write(*,"(A13, e12.4)")     ' q         = ' , &
    maxval(abs(log(q_mat_new/q_mat)))
    write(*,"(A13, e12.4, e12.4)")     ' nom_i     = ' , &
    maxval(abs((nom_i_mat_new(1,:)-nom_i_mat(1,:)))), & 
    maxval(abs((nom_i_mat_new(2,:)-nom_i_mat(2,:))))
    write(*,"(A13, e12.4, e12.4)")     ' infl      = ' , &
    maxval(abs((infl_mat_new(1,:)-infl_mat(1,:)))), & 
    maxval(abs((infl_mat_new(2,:)-infl_mat(2,:))))
    write(*,"(A13, e12.4, e12.4)")     ' l         = ' , &
    maxval(abs(log(l_aggr_mat_new(1,:)/l_aggr_mat(1,:)))), & 
    maxval(abs(log(l_aggr_mat_new(2,:)/l_aggr_mat(2,:))))
    write(*,"(A13, e12.4, e12.4)")     ' w         = ' , &
    maxval(abs(log(w_choice_mat_new(1,:)/w_choice_mat(1,:)))), & 
    maxval(abs(log(w_choice_mat_new(2,:)/w_choice_mat(2,:))))


    write(*,*)
    write(*,"(A13, e12.4, e12.4)") ' max/min v1 ', minval(v_mat_new(1,:)), maxval(v_mat_new(1,:))
    write(*,"(A13, e12.4, e12.4)") ' max/min v2 ', minval(v_mat_new(2,:)), maxval(v_mat_new(2,:))
    write(*,"(A13, e12.4, e12.4)") ' max/min m1 ', minval(mc_mat_new(1,:)), maxval(mc_mat_new(1,:))
    write(*,"(A13, e12.4, e12.4)") ' max/min m2 ', minval(mc_mat_new(2,:)), maxval(mc_mat_new(2,:))
    write(*,"(A13, e12.4, e12.4)") ' max/min s1 ', minval(bF_share_mat_new(1,:)), maxval(bF_share_mat_new(1,:))
    write(*,"(A13, e12.4, e12.4)") ' max/min s2 ', minval(bF_share_mat_new(2,:)), maxval(bF_share_mat_new(2,:))
    write(*,*)
    
    ! update value function
    v_mat           = v_mat + 1.0*(v_mat_new-v_mat )
    mc_mat          = mc_mat + 1.0*(mc_mat_new - mc_mat)
    
    next_state_mat = next_state_mat +  0.1_dp*(next_state_mat_new-next_state_mat)
      
     q_mat = q_mat + 0.2_dp*(q_mat_new-q_mat)
 
     s_mat = s_mat + 0.1_dp*(s_mat_new-s_mat)
 
       l_aggr_mat = l_aggr_mat + 0.2_dp*(l_aggr_mat_new-l_aggr_mat)
 
       w_choice_mat = w_choice_mat*(1 + 0.2_dp*log(w_choice_mat_new/w_choice_mat))


    c_spending_mat = c_spending_mat      + 0.2*(c_spending_mat_new - c_spending_mat)
    nom_i_mat(1,:) = nom_i_mat(1,:)      + 0.1*(nom_i_mat_new(1,:) - nom_i_mat(1,:))
    nom_i_mat(2,:) = nom_i_mat(2,:)      + 0.1*(nom_i_mat_new(2,:) - nom_i_mat(2,:))
   
    
    infl_mat = infl_mat          + 0.05*(infl_mat_new - infl_mat)
 
    k_next_mat  = k_next_mat          + 0.2*(k_next_mat_new - k_next_mat)
    share_mat   = share_mat           + 0.1*(share_mat_new - share_mat) 
    
    if (outer_iter < 101) then
        bF_share_mat = bF_share_mat + 0.0_dp*(bF_share_mat_new-bF_share_mat)
    else 
        bF_share_mat = bF_share_mat + 0.05_dp*(bF_share_mat_new-bF_share_mat)
    endif
     
    
     
           
    enddo

! ============================================================================================================
! COLLECT RESULTS AND PRICE BONDS
! ============================================================================================================

    ! Recursive calculation: last period payoff = 1
    q_bond_mat(:,:,1) = 1.0_dp
    bbb_counter = 0   
    do bbb = 1,n_bond
        if (bbb <= 4 .or. bbb > n_bond-2) then 
            bbb_counter = bbb_counter + 1
        endif 
    
    ! solve for polynomial coefficients
    ! call omp_set_num_threads(4)
    call F07ABF( 'Equilibration','No transpose',n_states, nrhs2, smol_polynom, lda,af,ldaf,ipiv,    &
                  equed,r,c, q_bond_mat(:,:,bbb),ldb,bond_coeffs,ldx,rcond,ferr2,berr2,work,iwork,info)
    
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED( & 
    !$OMP next_state_mat, smol_elem_ani, smol_coeffs, state_grid, valuation_mat, &
    !$OMP nom_i_mat, s_mat, w_choice_mat, n_states, dz_vec, & 
    !$OMP k_grid_dev, k_grid_mean, zf_grid_dev, zf_grid_mean, q_mat, &
    !$OMP c_spending_mat, share_mat, infl_mat, varsigma_vec, sigma, & 
    !$OMP bF_share_mat, l_aggr_mat, next_state_mat_new, inv_h_mat, inv_f_mat, &
    !$OMP bF_share_mat_new, share_mat_new, bbb_counter, next_omg_mat, &
    !$OMP wh_grid_mean, wh_grid_dev, wf_grid_mean, wf_grid_dev, n_active_dims, &
    !$OMP tht_h_grid_mean, tht_h_grid_dev, &
    !$OMP w_choice_mat_new, infl_mat_new, max_smol_level, &
    !$OMP q_mat_new, l_aggr_mat_new, v_mat_new, s_mat_new, nom_i_mat_new, &
    !$OMP c_spending_mat_new, k_next_mat_new, dz_vec_adj, next_zf_mat, &
    !$OMP results_mat, k_next_mat, bbb, bond_coeffs, q_bond_mat, n_quad) &
    !$OMP PRIVATE( &
    !$OMP polyn_points, nxt_mat, inv_h, inv_f, &
    !$OMP nom_i_vec, s, nxt_mat_2, counter, who_prices_bond, &
    !$OMP q_current, c_spending_vec, share_vec, nxt_bond_prices, &
    !$OMP bF_share_vec, infl_vec , l_aggr_vec,E_rb,  &
    !$OMP q_new, s_new, k_next_new, c_spending_vec_new, infl_new, &
    !$OMP sss, v_new, l_aggr_new, w_choice_new, results_vec)
    !$OMP DO SCHEDULE(static)
    do sss = 1, n_states 
    
        ! interpolate future values 
        polyn_points = Smolyak_Polynomial2(next_state_mat(:,:,sss), n_active_dims, &
                                  n_quad, n_states, max_smol_level, smol_elem_ani)

        ! use MKL Library for multiplication, same as
        ! nxt_mat = matmul(polyn_points, smol_coeffs)
        CALL DGEMM('N','N', n_quad, nrhs, n_states, 1.0_dp, polyn_points, n_quad, & 
                           smol_coeffs, n_states, 0.0_dp, & 
                           nxt_mat, n_quad)   

        CALL DGEMM('N','N', n_quad, nrhs2, n_states, 1.0_dp, polyn_points, n_quad, & 
                           bond_coeffs, n_states, 0.0_dp, & 
                           nxt_bond_prices, n_quad)   

        ! non interpolated values
        nxt_mat_2(:, 1) = k_next_mat(:,sss) 
        nxt_mat_2(:, 2) = next_zf_mat(:,sss) 
        nxt_mat_2(:, 3) = next_omg_mat(:,sss) 

        ! current state and guesses 
        nom_i_vec         = nom_i_mat(:,sss)
        s                 = s_mat(sss)
        q_current         = q_mat(sss)
        c_spending_vec    = c_spending_mat(:,sss)
        share_vec         = share_mat(:,sss)
        bF_share_vec    = bF_share_mat(:,sss)
        infl_vec          = infl_mat(:,sss)
        l_aggr_vec        = l_aggr_mat(:,sss)
        inv_h             = inv_h_mat(sss)
        inv_f             = inv_f_mat(sss)
        
        call calc_bond_prices( nxt_mat, n_nxt, nxt_mat_2, nxt_bond_prices, sss, &
                    c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                    nom_i_vec, share_vec, bF_share_vec, inv_h, inv_f, 0.0_dp, 0.0_dp, & 
                    who_prices_bond, q_bond_mat(sss,:,bbb+1), E_rb, results_vec)
    
        ! stuff stored only for results
        results_mat(sss,1:n_interp-6*12) = results_vec
       
        if (bbb <= 4 .or. bbb > n_bond-2) then 
        results_mat(sss,n_interp-6*12 + 1 + (bbb_counter-1)*3:n_interp-6*12+ 3*bbb_counter) = q_bond_mat(sss,:,bbb+1)
        results_mat(sss,n_interp-6*9  + 1 + (bbb_counter-1)*6:n_interp-6*9 + 6*bbb_counter) = E_rb
        results_mat(sss,n_interp-6*3  + 1 + (bbb_counter-1)*3:n_interp-6*3 + 3*bbb_counter) = who_prices_bond
        endif 


        enddo
    !$OMP END DO
    !$OMP END PARALLEL

    enddo
    
    !! CALC VALUATION
    valuation_mat = 0.0_dp

    do bbb = 1,n_valuation - 1
    
    ! solve for polynomial coefficients
    call F07ABF( 'Equilibration','No transpose',n_states, nrhs3, smol_polynom, lda,af,ldaf,ipiv,    &
                  equed,r,c, valuation_mat,ldb, val_coeffs,ldx,rcond,ferr3,berr3,work,iwork,info)
    
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED( & 
    !$OMP next_state_mat, smol_elem_ani, smol_coeffs, state_grid, valuation_mat, &
    !$OMP nom_i_mat, s_mat, w_choice_mat, n_states, dz_vec, & 
    !$OMP k_grid_dev, k_grid_mean, zf_grid_dev, zf_grid_mean, q_mat, &
    !$OMP c_spending_mat, share_mat, infl_mat, varsigma_vec, sigma, & 
    !$OMP bF_share_mat, l_aggr_mat, next_state_mat_new, inv_h_mat, inv_f_mat, &
    !$OMP bF_share_mat_new, share_mat_new, valuation_mat_tmp, next_omg_mat, &
    !$OMP wh_grid_mean, wh_grid_dev, wf_grid_mean, wf_grid_dev, n_active_dims, &
    !$OMP w_choice_mat_new, infl_mat_new, max_smol_level, &
    !$OMP q_mat_new, l_aggr_mat_new, v_mat_new, s_mat_new, nom_i_mat_new, &
    !$OMP c_spending_mat_new, k_next_mat_new, dz_vec_adj, next_zf_mat, &
    !$OMP results_mat, k_next_mat, bbb, bond_coeffs, q_bond_mat, n_quad, val_coeffs) &
    !$OMP PRIVATE( &
    !$OMP polyn_points, nxt_mat, inv_h, inv_f, &
    !$OMP nom_i_vec, s, nxt_mat_2, counter, who_prices_bond, &
    !$OMP q_current, c_spending_vec, share_vec, nxt_bond_prices, &
    !$OMP bF_share_vec, infl_vec , l_aggr_vec, &
    !$OMP q_new, s_new, k_next_new, c_spending_vec_new, infl_new, &
    !$OMP sss, v_new, l_aggr_new, w_choice_new, results_vec, nxt_valuation, valuation_vec)
    !$OMP DO SCHEDULE(static)
    do sss = 1, n_states 
    
        ! interpolate future values 
        polyn_points = Smolyak_Polynomial2(next_state_mat(:,:,sss), n_active_dims, &
                                  n_quad, n_states, max_smol_level, smol_elem_ani)

        ! use MKL Library for multiplication, same as
        ! nxt_mat = matmul(polyn_points, smol_coeffs)
        CALL DGEMM('N','N', n_quad, nrhs, n_states, 1.0_dp, polyn_points, n_quad, & 
                           smol_coeffs, n_states, 0.0_dp, & 
                           nxt_mat, n_quad)   

        CALL DGEMM('N','N', n_quad, nrhs3, n_states, 1.0_dp, polyn_points, n_quad, & 
                    val_coeffs, n_states, 0.0_dp, nxt_valuation, n_quad)   

        ! non interpolated values
        nxt_mat_2(:, 1) = k_next_mat(:,sss) 
        nxt_mat_2(:, 2) = next_zf_mat(:,sss) 
        nxt_mat_2(:, 3) = next_omg_mat(:,sss) 

        ! current state and guesses 
        nom_i_vec         = nom_i_mat(:,sss)
        s                 = s_mat(sss)
        q_current         = q_mat(sss)
        c_spending_vec    = c_spending_mat(:,sss)
        share_vec         = share_mat(:,sss)
        bF_share_vec    = bF_share_mat(:,sss)
        infl_vec          = infl_mat(:,sss)
        l_aggr_vec        = l_aggr_mat(:,sss)
        inv_h             = inv_h_mat(sss)
        inv_f             = inv_f_mat(sss)
        
        call calc_valuation( nxt_mat, n_nxt, nxt_mat_2, nxt_valuation, sss, &
                    c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                    nom_i_vec, share_vec, bF_share_vec, inv_h, inv_f, 0.0_dp, 0.0_dp, valuation_vec)
        
        if (bbb == 1) then ! in run 1 assign period H nfa term 
            valuation_vec(2) = valuation_vec(1)
        endif

        valuation_mat(sss,:) = valuation_vec
         
        if (bbb == n_valuation - 2) then 
            valuation_mat_tmp(sss,:) = valuation_mat(sss,:)
        endif
    

        enddo
    !$OMP END DO
    !$OMP END PARALLEL
    enddo

    call F07ABF( 'Equilibration','No transpose',n_states, nrhs3, smol_polynom, lda,af,ldaf,ipiv,    &
                  equed,r,c, valuation_mat,ldb, val_coeffs,ldx,rcond,ferr3,berr3,work,iwork,info)

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED( & 
    !$OMP next_state_mat, smol_elem_ani, smol_coeffs, max_smol_level,  &
    !$OMP n_states, n_active_dims, quad_weight_vec, state_grid, &
    !$OMP bbb, val_coeffs, n_quad, valuation_mat) &
    !$OMP PRIVATE( &
    !$OMP polyn_points, p_dis, nxt_valuation, big_weight_vec,&
    !$OMP sss)
    !$OMP DO SCHEDULE(static)
    do sss = 1, n_states 
    
        ! interpolate future values 
        polyn_points = Smolyak_Polynomial2(next_state_mat(:,:,sss), n_active_dims, &
                                  n_quad, n_states, max_smol_level, smol_elem_ani)

        CALL DGEMM('N','N', n_quad, nrhs3, n_states, 1.0_dp, polyn_points, n_quad, & 
                    val_coeffs, n_states, 0.0_dp, nxt_valuation, n_quad)   
        p_dis       = exp(state_grid(idx_dis, sss))

        ! current weight vec
        big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
        big_weight_vec(n_quad)     = p_dis

        ! shift nx term to expected nx 
        valuation_mat(sss,4)       = sum(big_weight_vec*nxt_valuation(:,4))
        
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call F07ABF( 'Equilibration','No transpose',n_states, nrhs3, smol_polynom, lda,af,ldaf,ipiv,    &
                  equed,r,c, valuation_mat,ldb, val_coeffs,ldx,rcond,ferr3,berr3,work,iwork,info)

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED( & 
    !$OMP next_state_mat, smol_elem_ani, smol_coeffs, valuation_mat_store,  &
    !$OMP n_states, n_active_dims, quad_weight_vec, max_smol_level, state_grid, &
    !$OMP bbb, val_coeffs, n_quad, valuation_mat) &
    !$OMP PRIVATE( &
    !$OMP polyn_points, p_dis, nxt_valuation, big_weight_vec,&
    !$OMP sss)
    !$OMP DO SCHEDULE(static)
    do sss = 1, n_states 
    
        ! interpolate future values 
        polyn_points = Smolyak_Polynomial2(next_state_mat(:,:,sss), n_active_dims, &
                                  n_quad, n_states, max_smol_level, smol_elem_ani)
        CALL DGEMM('N','N', n_quad, nrhs3, n_states, 1.0_dp, polyn_points, n_quad, & 
                    val_coeffs, n_states, 0.0_dp, nxt_valuation, n_quad)   

        p_dis       = exp(state_grid(idx_dis, sss))

        ! current weight vec
        big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
        big_weight_vec(n_quad)     = p_dis
        
        ! calculated one period ahead expectations of above terms 
        valuation_mat_store(1  ,sss) = valuation_mat(sss,1)
        valuation_mat_store(2  ,sss) = valuation_mat(sss,2)
        valuation_mat_store(3  ,sss) = valuation_mat(sss,3)
        valuation_mat_store(4  ,sss) = valuation_mat(sss,4)
        valuation_mat_store(5  ,sss) = valuation_mat(sss,5)
        valuation_mat_store(6  ,sss) = valuation_mat(sss,6)
        valuation_mat_store(7  ,sss) = valuation_mat(sss,7)
        valuation_mat_store(7+1,sss) = sum(big_weight_vec*nxt_valuation(:,1))
        valuation_mat_store(7+2,sss) = sum(big_weight_vec*nxt_valuation(:,2))
        valuation_mat_store(7+3,sss) = sum(big_weight_vec*nxt_valuation(:,3))
        valuation_mat_store(7+4,sss) = sum(big_weight_vec*nxt_valuation(:,4))
        valuation_mat_store(7+5,sss) = sum(big_weight_vec*nxt_valuation(:,5))
        valuation_mat_store(7+6,sss) = sum(big_weight_vec*nxt_valuation(:,6))
        valuation_mat_store(7+7,sss) = sum(big_weight_vec*nxt_valuation(:,7))
        
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    write(*,*) 'MODEL SOLUTION DONE - STORING SOLUTION FILES'
    write(*,*) 
    flush(6)

    ! print all results
    open (unit = 10, file = trim(results_folder) // 'next_state_mat.dat', ACTION="write", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    write(10) next_state_mat; close(10)

    ! print all results
    open (unit = 10, file = trim(results_folder) // 'results_mat.dat', ACTION="write", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    write(10) results_mat; close(10)

    ! print all results
    open (unit = 10, file = trim(results_folder) // 'valuation_mat.dat', ACTION="write", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    write(10) valuation_mat_store; close(10)


! ============================================================================================================
! IRF CALCULATIONS 
! ============================================================================================================
write(*,*) 'IRF CALCULATIONS'
do fff = 1,size(irf_indices,1)

     ! for standard shocks: next period values = model solution 
     next_state_irf_mat     = next_state_mat
     next_state_irf_mat_new = next_state_mat
     k_next_irf_mat         = k_next_mat
     nom_i_irf_mat          = nom_i_mat 
     s_irf_mat              = s_mat
     q_irf_mat              = q_mat
     c_spending_irf_mat     = c_spending_mat 
     share_irf_mat          = share_mat 
     bF_share_irf_mat     = bF_share_mat 
     infl_irf_mat           = infl_mat 
     l_aggr_irf_mat         = l_aggr_mat 
     inv_h_irf_mat          = inv_h_mat
     inv_f_irf_mat          = inv_f_mat
     w_choice_irf_mat       = w_choice_mat 
     v_irf_mat              = v_mat 
     mc_irf_mat             = mc_mat
  

    ! calculation transition 
    if (irf_indices(fff) == sidx_dis) then 
        next_state_irftransition_mat(1,:,:) = next_state_mat(n_quad,:,:)
    else

        trans_shock_vec = shock_grid(no_shock_idx(1),:)
        if (irf_indices(fff) == sidx_z) then 
            trans_shock_vec(irf_indices(fff)) = irf_shock_sizes(fff) 
        elseif (irf_indices(fff) == sidx_p) then 
            trans_shock_vec(irf_indices(fff)) = irf_shock_sizes(fff)
        elseif (irf_indices(fff) == sidx_omg) then 
            trans_shock_vec(irf_indices(fff)) = irf_shock_sizes(fff)
        elseif (irf_indices(fff) == sidx_zf) then 
            trans_shock_vec(irf_indices(fff)) = irf_shock_sizes(fff)
        endif

        next_state_irftransition_mat     = next_state_mat(no_shock_idx,:,:)
        next_state_irftransition_mat_new = next_state_mat(no_shock_idx,:,:)

        do iii = 1,n_I
            interp_input_mat(:,iii)        = v_irf_mat(iii,:)
            interp_input_mat(:,n_I + iii)  = mc_irf_mat(iii,:)
        enddo

        interp_input_mat(:,2*n_I + 1)  = s_irf_mat
        interp_input_mat(:,2*n_I + 2)  = q_irf_mat
        interp_input_mat(:,2*n_I + 3)  = l_aggr_irf_mat(1,:)
        interp_input_mat(:,2*n_I + 4)  = l_aggr_irf_mat(2,:)
        interp_input_mat(:,2*n_I + 5)  = infl_irf_mat(1,:)
        interp_input_mat(:,2*n_I + 6)  = infl_irf_mat(2,:)
        interp_input_mat(:,2*n_I + 7) = c_spending_irf_mat(1,:)
        interp_input_mat(:,2*n_I + 8) = c_spending_irf_mat(2,:)
        
        ! solve for polynomial coefficients
        call F07ABF( 'Equilibration','No transpose',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,    &
                      equed,r,c,interp_input_mat,ldb,smol_coeffs,ldx,rcond,ferr,berr,work,iwork,info)
        
        diff = 1.0_dp
        outer_iter = 0
        do while (diff > conv_crit)
        outer_iter = outer_iter + 1

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED( & 
        !$OMP next_state_irftransition_mat, smol_elem_ani, smol_coeffs, state_grid, &
        !$OMP nom_i_mat, s_mat, w_choice_mat, n_states, dz_vec, vector_mus_dimensions, & 
        !$OMP k_grid_dev, k_grid_mean, zf_grid_dev, zf_grid_mean, q_mat, next_omg_mat, &
        !$OMP c_spending_mat, share_mat, infl_mat, varsigma_vec, sigma, & 
        !$OMP bF_share_mat, l_aggr_mat, next_state_irftransition_mat_new, inv_h_mat, inv_f_mat, &
        !$OMP bF_share_mat_new, share_mat_new, disast_p, &
        !$OMP wh_grid_mean, wh_grid_dev, wf_grid_mean, wf_grid_dev, n_active_dims, next_dis_mat, &
        !$OMP omg_grid_mean, omg_grid_dev, dis_grid_mean, dis_grid_dev, no_shock_idx, &
        !$OMP tht_h_grid_mean, tht_h_grid_dev, max_smol_level, &
        !$OMP w_choice_mat_new, infl_mat_new, &
        !$OMP q_mat_new, l_aggr_mat_new, v_mat_new, s_mat_new, nom_i_mat_new, trans_shock_vec, &
        !$OMP c_spending_mat_new, k_next_mat_new, dz_vec_adj, next_zf_mat, &
        !$OMP k_next_mat) &
        !$OMP PRIVATE( &
        !$OMP polyn_points_trans, nxt_mat_trans, inv_h, inv_f, &
        !$OMP nom_i_vec, s, nxt_mat_2_trans, counter, &
        !$OMP q_current, c_spending_vec, share_vec, &
        !$OMP bF_share_vec, infl_vec , l_aggr_vec, &
        !$OMP q_new, s_new, k_next_new, c_spending_vec_new, infl_new, &
        !$OMP sss, v_new, l_aggr_new, w_choice_new, theta_nxt)
        !$OMP DO SCHEDULE(static)
        do sss = 1, n_states 
        

            ! interpolate future values 
            polyn_points_trans = Smolyak_Polynomial2(next_state_irftransition_mat(1,:,sss), n_active_dims, &
                                      1, n_states, max_smol_level, smol_elem_ani)
                                                                
            ! use MKL Library for multiplication, same as
            ! nxt_mat        = matmul(polyn_points, smol_coeffs)
            CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points_trans, 1, & 
                               smol_coeffs, n_states, 0.0_dp, & 
                               nxt_mat_trans, 1)   

            ! non interpolated values
            nxt_mat_2_trans(1, 1) = k_next_mat(no_shock_idx(1),sss) 
            nxt_mat_2_trans(1, 2) = next_zf_mat(no_shock_idx(1),sss) 
            nxt_mat_2_trans(1, 3) = next_omg_mat(no_shock_idx(1),sss) 

            ! current state and guesses 
            nom_i_vec         = nom_i_mat(:,sss)
            s                 = s_mat(sss)
            q_current         = q_mat(sss)
            c_spending_vec    = c_spending_mat(:,sss)
            share_vec         = share_mat(:,sss)
            bF_share_vec    = bF_share_mat(:,sss)
            l_aggr_vec        = l_aggr_mat(:,sss)
            infl_vec          = infl_mat(:,sss)
            
            call calc_unexpected_transition( nxt_mat_trans, n_nxt, nxt_mat_2_trans, sss, &
                        c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                        nom_i_vec, share_vec, bF_share_vec, trans_shock_vec,  &
                        0.0_dp, 0.0_dp, theta_nxt(1))
            
            w_choice_new = w_choice_mat(:,sss)

        counter = 0
        
        if (vector_mus_dimensions(1) > 0) then
            counter = counter + 1
            next_state_irftransition_mat_new(1, counter, sss)   = (k_next_mat(1,sss)/exp(trans_shock_vec(sidx_z)) -k_grid_mean )/k_grid_dev
        endif
        
        if (vector_mus_dimensions(2) > 0) then
            counter = counter + 1
            next_state_irftransition_mat_new(1, counter, sss) = (theta_nxt(1)  &
                                                           - tht_h_grid_mean)/tht_h_grid_dev
        endif

    
        if (vector_mus_dimensions(3) > 0) then
        counter = counter + 1
        next_state_irftransition_mat_new(1,counter,sss) =  (next_zf_mat(no_shock_idx(1),sss) + trans_shock_vec(sidx_zf)  - zf_grid_mean) /zf_grid_dev
        endif
    
        
    
        if (vector_mus_dimensions(4) > 0) then
            counter = counter + 1
            next_state_irftransition_mat_new(1, counter, sss)  = (w_choice_new(1)/exp(trans_shock_vec(sidx_z))/ & 
                                                          (varsigma_vec(1)    & 
                                                          + (1.0-varsigma_vec(1))*s**(sigma-1.0))**(1.0/(sigma-1.0)) & 
                                                          -wh_grid_mean) /wh_grid_dev
        endif
    
        if (vector_mus_dimensions(5) > 0) then
            counter = counter + 1
            next_state_irftransition_mat_new(1, counter, sss)  = (w_choice_new(2)/exp(trans_shock_vec(sidx_z))*s/ & 
                                                          (varsigma_vec(2)*s**(1.0_dp-sigma) + (1.0-varsigma_vec(2)))**(1.0/(sigma-1.0)) & 
                                                          -wf_grid_mean) /wf_grid_dev
        endif
    
        if (vector_mus_dimensions(6) > 0) then
            counter = counter + 1
            if (dis_grid_dev < sqrt_eps) then
            next_state_irftransition_mat_new(1,counter,sss) = 0.0_dp 
            else
            next_state_irftransition_mat_new(:,counter,sss) = (next_dis_mat(no_shock_idx(1),sss) + trans_shock_vec(sidx_p)   - dis_grid_mean) /dis_grid_dev
            endif
        endif
    
        if (vector_mus_dimensions(7) > 0) then
            counter = counter + 1
        endif
    
        if (vector_mus_dimensions(8) > 0) then
            counter = counter + 1
        endif
        
        if (vector_mus_dimensions(9) > 0) then
            counter = counter + 1
            if (omg_grid_dev < sqrt_eps) then
                next_state_irftransition_mat_new(1,counter,sss) = 0.0_dp 
            else
                next_state_irftransition_mat_new(:,counter,sss) = (next_omg_mat(no_shock_idx(1),sss) + trans_shock_vec(sidx_omg)   - omg_grid_mean) /omg_grid_dev
            endif
        endif

        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    
    where(next_state_irftransition_mat_new> 1.0_dp) 
        next_state_irftransition_mat_new= 1.0_dp 
    endwhere
    where(next_state_irftransition_mat_new < -1.0_dp)
        next_state_irftransition_mat_new   = -1.0_dp
    endwhere
    
    if (outer_iter < 10) then
        diff = 1.0_dp
    else
        diff =  maxval(abs(next_state_irftransition_mat - next_state_irftransition_mat_new)) 
    endif

    next_state_irftransition_mat = next_state_irftransition_mat +  0.2_dp*(next_state_irftransition_mat_new-next_state_irftransition_mat)

    enddo
    endif
    
        write(shock_char,'(i1)') irf_indices(fff)
        
        ! store transition matrix
        open (unit = 10, file = trim(results_folder) // 'next_state_shocktrans_mat_' // shock_char // '.dat', ACTION="write", STATUS="replace", &
                & FORM="unformatted", ACCESS="STREAM")
        write(10) next_state_irftransition_mat; close(10)

enddo
write(*,*) 'DONE'
write(*,*) 
flush(6)


if (run_bg == 1) then
! ============================================================================================================
! LONG-LIVED IRF CALCULATIONS 
! ============================================================================================================
! patch long-lived impulse response to model solution 
! count backwards

write(*,*) 'LONG-LIVED IRF CALCULATIONS'
    do fff = 1,size(irf_indices_ll,1)

        ! initially next period values as standard model solution
        next_state_irf_mat     = next_state_mat
        next_state_irf_mat_new = next_state_mat
        k_next_irf_mat         = k_next_mat
        nom_i_irf_mat          = nom_i_mat 
        s_irf_mat              = s_mat
        q_irf_mat              = q_mat
        c_spending_irf_mat     = c_spending_mat 
        share_irf_mat          = share_mat 
        bF_share_irf_mat     = bF_share_mat 
        infl_irf_mat           = infl_mat 
        l_aggr_irf_mat         = l_aggr_mat 
        inv_h_irf_mat          = inv_h_mat
        inv_f_irf_mat          = inv_f_mat
        w_choice_irf_mat       = w_choice_mat 
        v_irf_mat              = v_mat 
        mc_irf_mat             = mc_mat
 
        q_bond_irf_mat         = q_bond_mat


        ! count time backwards
        do ttt = 1000 + irf_ll_length, 1001, -1

            bond_supply_shock = 0.0 
            bond_supply_shock_nxt = 0.0 

        
        ! next period values
        do iii = 1,n_I
            interp_input_mat(:,iii)  = v_irf_mat(iii,:)
            interp_input_mat(:,n_I+iii)  = mc_irf_mat(iii,:)
        enddo
            
        interp_input_mat(:,2*n_I + 1)  = s_irf_mat
        interp_input_mat(:,2*n_I + 2)  = q_irf_mat
        interp_input_mat(:,2*n_I + 3)  = l_aggr_irf_mat(1,:)
        interp_input_mat(:,2*n_I + 4)  = l_aggr_irf_mat(2,:)
        interp_input_mat(:,2*n_I + 5)  = infl_irf_mat(1,:)
        interp_input_mat(:,2*n_I + 6)  = infl_irf_mat(2,:)
        interp_input_mat(:,2*n_I + 7) = c_spending_irf_mat(1,:)
        interp_input_mat(:,2*n_I + 8) = c_spending_irf_mat(2,:)
            
        ! solve for polynomial coefficients
        call F07ABF( 'Equilibration','No transpose',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,    &
                      equed,r,c,interp_input_mat,ldb,smol_coeffs,ldx,rcond,ferr,berr,work,iwork,info)
        
        if (ttt > 1000 + irf_ll_hold_fixed(fff)) then 
        outer_iter = 0
        diff = 1.0_dp
        do while (diff > conv_crit .and. outer_iter < max_iter_irf) 
            outer_iter = outer_iter + 1
    
            !$OMP PARALLEL DEFAULT(NONE) &
            !$OMP SHARED( & 
            !$OMP next_state_irf_mat, smol_elem_ani, smol_coeffs, state_grid, vector_mus_dimensions, &
            !$OMP nom_i_irf_mat, s_irf_mat, n_states, dz_vec, disast_p, & 
            !$OMP k_grid_dev, k_grid_mean, zf_grid_dev, zf_grid_mean, q_irf_mat, &
            !$OMP c_spending_irf_mat, share_irf_mat, infl_irf_mat, & 
            !$OMP varsigma_vec, sigma, M_vec_mat, & 
            !$OMP bF_share_irf_mat, l_aggr_irf_mat, next_state_irf_mat_new, & 
            !$OMP inv_h_irf_mat, inv_f_irf_mat, &
            !$OMP bF_share_irf_mat_new, share_irf_mat_new, next_omg_mat, &
            !$OMP wh_grid_mean, wh_grid_dev, wf_grid_mean, wf_grid_dev, n_active_dims, mc_irf_mat_new, &
            !$OMP tht_h_grid_mean, tht_h_grid_dev, max_smol_level, &
            !$OMP w_choice_irf_mat_new, infl_irf_mat_new, &
            !$OMP q_irf_mat_new, l_aggr_irf_mat_new, v_irf_mat_new, s_irf_mat_new, & 
            !$OMP nom_i_irf_mat_new, w_choice_irf_mat, n_quad, bond_supply_shock, bond_supply_shock_nxt, &
            !$OMP c_spending_irf_mat_new, k_next_irf_mat_new, dz_vec_adj, next_zf_mat, &
            !$OMP k_next_irf_mat, ih_grid_mean, ih_grid_dev, if_grid_mean, if_grid_dev) &
            !$OMP PRIVATE( &
            !$OMP polyn_points, nxt_mat, inv_h, inv_f, mc_new, &
            !$OMP nom_i_vec, s, nxt_mat_2, counter, &
            !$OMP q_current, c_spending_vec, share_vec, &
            !$OMP bF_share_vec, infl_vec , l_aggr_vec, &
            !$OMP q_new, s_new, k_next_new, c_spending_vec_new, infl_new, &
            !$OMP sss, v_new, l_aggr_new, w_choice_new, theta_nxt)
            !$OMP DO SCHEDULE(static)
            do sss = 1, n_states 
           
                ! interpolate future values 
                polyn_points = Smolyak_Polynomial2(next_state_irf_mat(:,:,sss), n_active_dims, &
                                          n_quad, n_states, max_smol_level, smol_elem_ani)
    
                ! use MKL Library for multiplication, same as
                ! nxt_mat        = matmul(polyn_points, smol_coeffs)
                CALL DGEMM('N','N', n_quad, nrhs, n_states, 1.0_dp, polyn_points, n_quad, & 
                                   smol_coeffs, n_states, 0.0_dp, nxt_mat, n_quad)   
    
                ! non interpolated values
                nxt_mat_2(:, 1) = k_next_irf_mat(:,sss) 
                nxt_mat_2(:, 2) = next_zf_mat(:,sss) 
                nxt_mat_2(:, 3) = next_omg_mat(:,sss) 
    
                ! current state and guesses 
                nom_i_vec         = nom_i_irf_mat(:,sss)
                s                 = s_irf_mat(sss)
                q_current         = q_irf_mat(sss)
                c_spending_vec    = c_spending_irf_mat(:,sss)
                share_vec         = share_irf_mat(:,sss)
                bF_share_vec    = bF_share_irf_mat(:,sss)
                infl_vec          = infl_irf_mat(:,sss)
                l_aggr_vec        = l_aggr_irf_mat(:,sss)

                call calc_equilibrium_and_update( nxt_mat, n_nxt, nxt_mat_2, sss, 10000, &
                            c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                            nom_i_vec, share_vec, bF_share_vec, inv_h, inv_f,  & 
                            bond_supply_shock, bond_supply_shock_nxt, &
                            q_new, s_new, k_next_new, c_spending_vec_new, & 
                            infl_new, v_new, mc_new, l_aggr_new, w_choice_new, theta_nxt, M_vec_mat(:, :, sss))
    
                v_irf_mat_new(:,sss)           = v_new 
                mc_irf_mat_new(:, sss)         = mc_new
                s_irf_mat_new(sss)             = s_new
                c_spending_irf_mat_new(:,sss)  = c_spending_vec_new
                k_next_irf_mat_new(:,sss)      = k_next_new/exp(dz_vec_adj)*exp(dz_vec)
                q_irf_mat_new(sss)             = q_new 
                inv_f_irf_mat(sss)         = inv_f 
                inv_h_irf_mat(sss)         = inv_h 
                l_aggr_irf_mat_new(:,sss)      = l_aggr_new
                nom_i_irf_mat_new(:,sss)       = nom_i_vec
                w_choice_irf_mat_new(:,sss)    = w_choice_new
                share_irf_mat_new(:,sss)       = share_vec          
                bF_share_irf_mat_new(:,sss)  = bF_share_vec     
                infl_irf_mat_new(:,sss)        = infl_new
                
    
                counter = 0
                
                if (vector_mus_dimensions(1) > 0) then
                    counter = counter + 1
                    next_state_irf_mat_new(:, counter, sss)   = (k_next_new/exp(dz_vec_adj) -k_grid_mean )/k_grid_dev
                endif
    
                if (vector_mus_dimensions(2) > 0) then
                    counter = counter + 1
                    next_state_irf_mat_new(:, counter, sss) = (theta_nxt - tht_h_grid_mean)/tht_h_grid_dev
                endif
    
                if (vector_mus_dimensions(3) > 0) then
                counter = counter + 1
                endif
    
                if (vector_mus_dimensions(4) > 0) then
                    counter = counter + 1
                    next_state_irf_mat_new(:, counter, sss)  = (w_choice_new(1)/exp(dz_vec_adj)/ & 
                                                                  (varsigma_vec(1)  + & 
                                                                  (1.0-varsigma_vec(1))*s**(sigma-1.0))**(1.0/(sigma-1.0)) & 
                                                                  -wh_grid_mean) /wh_grid_dev
                endif
    
                if (vector_mus_dimensions(5) > 0) then
                    counter = counter + 1
                    next_state_irf_mat_new(:, counter, sss)  = (w_choice_new(2)/exp(dz_vec_adj)*s/ & 
                                                                  (varsigma_vec(2)*s**(1.0_dp-sigma) + (1.0-varsigma_vec(2)))**(1.0/(sigma-1.0)) & 
                                                                  -wf_grid_mean) /wf_grid_dev
                endif
    
                if (vector_mus_dimensions(6) > 0) then
                counter = counter + 1
                endif
    
                if (vector_mus_dimensions(7) > 0) then
                    counter = counter + 1
                    next_state_irf_mat_new(:, counter, sss)  = ((nom_i_vec(1)-1.0) - ih_grid_mean) /ih_grid_dev 
                endif
    
                if (vector_mus_dimensions(8) > 0) then
                    counter = counter + 1
                    next_state_irf_mat_new(:, counter, sss)  = ((nom_i_vec(1) + nom_i_vec(2)-1.0) - if_grid_mean) /if_grid_dev 
                endif

                if (vector_mus_dimensions(9) > 0) then
                counter = counter + 1
                endif
    
    

            enddo
            !$OMP END DO
            !$OMP END PARALLEL
        
            where(next_state_irf_mat_new > 1.0_dp) 
                next_state_irf_mat_new   = 1.0_dp 
            endwhere
            where(next_state_irf_mat_new < -1.0_dp)
                next_state_irf_mat_new   = -1.0_dp
            endwhere
    
        ! check convergence
        if (outer_iter < 10) then
        diff = 1.0_dp
        else
        ! check convergence
        diff =  maxval([ maxval(abs(log(v_irf_mat_new/v_irf_mat))), &
                         maxval(abs(log(mc_irf_mat_new/mc_irf_mat))), &
                         maxval(abs((c_spending_irf_mat_new-c_spending_irf_mat))), &
                         maxval(abs((share_irf_mat_new-share_irf_mat))), &
                         maxval(abs(log(k_next_irf_mat_new/k_next_irf_mat))), &
                         maxval(abs(log(l_aggr_irf_mat_new/l_aggr_irf_mat))), &
                         maxval(abs((nom_i_irf_mat_new-nom_i_irf_mat))), &
                         maxval(abs(log(w_choice_irf_mat_new/w_choice_irf_mat))), &
                         maxval(abs(log(q_irf_mat_new/q_irf_mat))), &
                         maxval(abs(log(infl_irf_mat_new/infl_irf_mat))), &
                         maxval(abs(next_state_irf_mat - next_state_irf_mat_new)) ])
    
        endif
    
        timer_end = omp_get_wtime()
    
        write(*,*) ''
        write(*,"(A15, i5, i5, i5)") ' IRF iteration ', fff, outer_iter, ttt
        write(*,"(A11, f8.2)") ' Calc. time', timer_end-timer_start
        write(*,*) ''
    
        write(*,"(A19,e12.4,A1)") 'Changes (max diff: ', diff , ')'
        write(*,*) '------------------------------------------------------------------------'
        write(*,"(A13, e12.4, e12.4)")     ' v         = ' , &
        maxval(abs(log(v_irf_mat_new(1,:)/v_irf_mat(1,:)))), & 
        maxval(abs(log(v_irf_mat_new(2,:)/v_irf_mat(2,:))))
        write(*,"(A13, e12.4, e12.4)")     ' mc         = ' , &
        maxval(abs(log(mc_irf_mat_new(1,:)/mc_irf_mat(1,:)))), & 
        maxval(abs(log(mc_irf_mat_new(2,:)/mc_irf_mat(2,:))))
        write(*,"(A13, e12.4, e12.4)")     ' c_spend   = ' , &
        maxval(abs(log(c_spending_irf_mat_new(1,:)/c_spending_irf_mat(1,:)))), & 
        maxval(abs(log(c_spending_irf_mat_new(2,:)/c_spending_irf_mat(2,:))))
        write(*,"(A13, e12.4, e12.4)")     ' share_mat = ' , &
        maxval(abs((share_irf_mat_new(1,:)-share_irf_mat(1,:)))), & 
        maxval(abs((share_irf_mat_new(2,:)-share_irf_mat(2,:))))
        write(*,"(A13, e12.4)")     ' k         = ' , &
        maxval(abs(log(k_next_irf_mat_new/k_next_irf_mat)))
        write(*,"(A13, e12.4)")     ' s         = ' , &
        maxval(abs(log(s_irf_mat_new/s_irf_mat)))
        do sss = 1,n_active_dims 
        write(*,"(A13, e12.4)")     ' next_state= ' , &
        maxval(abs((next_state_irf_mat_new(:,sss,:)-next_state_irf_mat(:,sss,:))))
        enddo
        write(*,"(A13, e12.4)")     ' q         = ' , &
        maxval(abs(log(q_irf_mat_new/q_irf_mat)))
        write(*,"(A13, e12.4, e12.4)")     ' nom_i     = ' , &
        maxval(abs((nom_i_irf_mat_new(1,:)-nom_i_irf_mat(1,:)))), & 
        maxval(abs((nom_i_irf_mat_new(2,:)-nom_i_irf_mat(2,:))))
        write(*,"(A13, e12.4, e12.4)")     ' l         = ' , &
        maxval(abs(log(l_aggr_irf_mat_new(1,:)/l_aggr_irf_mat(1,:)))), & 
        maxval(abs(log(l_aggr_irf_mat_new(2,:)/l_aggr_irf_mat(2,:))))
        write(*,"(A13, e12.4, e12.4)")     ' w         = ' , &
        maxval(abs(log(w_choice_irf_mat_new(1,:)/w_choice_irf_mat(1,:)))), & 
        maxval(abs(log(w_choice_irf_mat_new(2,:)/w_choice_irf_mat(2,:))))
        
        ! update value function
        v_irf_mat           = v_irf_mat + 1.0*(v_irf_mat_new-v_irf_mat )
        mc_irf_mat          = mc_irf_mat + 1.0*(mc_irf_mat_new - mc_irf_mat)
        next_state_irf_mat = next_state_irf_mat +  0.2_dp*(next_state_irf_mat_new-next_state_irf_mat)
         
        where (abs((q_irf_mat-q_irf_mat_new)) > 0.01)
           q_irf_mat = q_irf_mat + 0.002_dp*sign(1.0_dp, q_irf_mat_new - q_irf_mat)
        elsewhere
           q_irf_mat = q_irf_mat + 0.1_dp*(q_irf_mat_new-q_irf_mat)
        endwhere
    
        where (abs((s_irf_mat-s_irf_mat_new)) > 0.05)
           s_irf_mat = s_irf_mat + 0.025_dp*sign(1.0_dp, s_irf_mat_new - s_irf_mat)
        elsewhere
           s_irf_mat = s_irf_mat + 0.5_dp*(s_irf_mat_new-s_irf_mat)
        endwhere
    
        where (abs(l_aggr_irf_mat-l_aggr_irf_mat_new) > 0.01)
           l_aggr_irf_mat = l_aggr_irf_mat + 0.002_dp*sign(1.0_dp, l_aggr_irf_mat_new - l_aggr_irf_mat)
        elsewhere
           l_aggr_irf_mat = l_aggr_irf_mat + 0.2_dp*(l_aggr_irf_mat_new-l_aggr_irf_mat)
        endwhere
    
        where (abs(log(w_choice_irf_mat/w_choice_irf_mat_new)) > 0.01)
           w_choice_irf_mat = w_choice_irf_mat*(1 + 0.002_dp*sign(1.0_dp, w_choice_irf_mat_new - w_choice_irf_mat))
        elsewhere
           w_choice_irf_mat = w_choice_irf_mat*(1 + 0.2_dp*log(w_choice_irf_mat_new/w_choice_irf_mat))
        endwhere
    
    
        c_spending_irf_mat = c_spending_irf_mat      + 0.2*(c_spending_irf_mat_new - c_spending_irf_mat)
        nom_i_irf_mat      = nom_i_irf_mat           + 0.5*(nom_i_irf_mat_new - nom_i_irf_mat)
        infl_irf_mat       = infl_irf_mat            + 0.2*(infl_irf_mat_new - infl_irf_mat)
        k_next_irf_mat     = k_next_irf_mat          + 0.2*(k_next_irf_mat_new - k_next_irf_mat)
        share_irf_mat      = share_irf_mat           + 0.5*(share_irf_mat_new - share_irf_mat) 
        bF_share_irf_mat = bF_share_irf_mat      + 0.05_dp*(bF_share_irf_mat_new-bF_share_irf_mat)
    
        enddo
        else  

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED( & 
        !$OMP next_state_irf_mat, smol_elem_ani, smol_coeffs, state_grid, &
        !$OMP nom_i_irf_mat, s_irf_mat, w_choice_irf_mat, n_states, dz_vec, & 
        !$OMP k_grid_dev, k_grid_mean, zf_grid_dev, zf_grid_mean, q_irf_mat, &
        !$OMP c_spending_irf_mat, share_irf_mat, infl_irf_mat, varsigma_vec, &  
        !$OMP sigma, bbb_counter, & 
        !$OMP bF_share_irf_mat, l_aggr_irf_mat, next_state_irf_mat_new, & 
        !$OMP inv_h_irf_mat, inv_f_irf_mat, &
        !$OMP bF_share_irf_mat_new, share_irf_mat_new, bond_supply_shock, bond_supply_shock_nxt, &
        !$OMP wh_grid_mean, wh_grid_dev, wf_grid_mean, wf_grid_dev, n_active_dims, v_irf_mat, mc_irf_mat, &
        !$OMP n_quad, next_omg_mat, &
        !$OMP w_choice_irf_mat_new, infl_irf_mat_new, &
        !$OMP q_irf_mat_new, l_aggr_irf_mat_new, v_irf_mat_new, s_irf_mat_new, nom_i_irf_mat_new, &
        !$OMP c_spending_irf_mat_new, k_next_irf_mat_new, dz_vec_adj, next_zf_mat, &
        !$OMP results_irf_mat, k_next_irf_mat, bbb, bond_coeffs, q_bond_irf_mat) &
        !$OMP PRIVATE( &
        !$OMP ih_fixed, infl_A, infl_B, inv_h, inv_f, mc_new, next_state, excess_b,  &
        !$OMP nom_i_vec, s, brnt_xbracket, brnt_ybracket, brnt_xnew, brnt_ynew, brnt_xvec, brnt_yvec,  &
        !$OMP q_current, c_spending_vec, share_vec, brnt_mflag, brnt_sflag, brnt_iter, &
        !$OMP bF_share_vec, infl_vec , l_aggr_vec, k_next, next_zf, step_sign, step_sign_prev,  &
        !$OMP q_new, s_new, k_next_new, c_spending_vec_new, infl_new, nxt_omg, &
        !$OMP sss, v_new, l_aggr_new, w_choice_new, bracket_success, excess_b_B)
        !$OMP DO SCHEDULE(static)
        do sss = 1,n_states 

            
            
            ! get initial guess 
            ih_fixed = 1.01 
            infl_A   = 1.0
            bracket_success = 0 
            step_sign_prev = 0.0_dp
            do while (bracket_success == 0)

                next_state        = next_state_irf_mat(:,:,sss)
                k_next            = k_next_irf_mat(:,sss) 
                next_zf          = next_zf_mat(:,sss) 
                nxt_omg           = next_omg_mat(:,sss) 
                nom_i_vec         = nom_i_irf_mat(:,sss)
                s                 = s_irf_mat(sss)
                q_current         = q_irf_mat(sss)
                c_spending_vec    = c_spending_irf_mat(:,sss)
                share_vec         = share_irf_mat(:,sss)
                bF_share_vec    = bF_share_irf_mat(:,sss)
                infl_vec          = infl_irf_mat(:,sss)
                l_aggr_vec        = l_aggr_irf_mat(:,sss)
                inv_h             = inv_h_irf_mat(sss)
                inv_f             = inv_f_irf_mat(sss)

                call update_ifixed(infl_A, ih_fixed, next_state, n_nxt, sss, smol_coeffs, k_next, next_zf, nxt_omg, nom_i_vec, s, q_current, & 
                             c_spending_vec, share_vec, bF_share_vec, infl_vec, l_aggr_vec, inv_h, inv_f, & 
                             bond_supply_shock, bond_supply_shock_nxt, &
                             excess_b, v_new, mc_new, w_choice_new &
                             )
                if (excess_b > 0) then 
                    step_sign = -1.0_dp
                else 
                    step_sign = 1.0_dp
                end if
                
                if (step_sign == step_sign_prev .or. step_sign_prev == 0.0_dp) then 
                    infl_B = infl_A 
                    excess_b_B = excess_b
                    infl_A = infl_A + step_sign*0.01
                    step_sign_prev = step_sign
                    if (infl_A > 1.1) then 
                       infl_A = 1.1 
                       ih_fixed = ih_fixed + 0.01
                       step_sign_prev = 0.0_dp
                    elseif (infl_A < 0.9) then 
                       infl_A = 0.9 
                       ih_fixed = ih_fixed - 0.01
                       step_sign_prev = 0.0_dp
                    end if
                elseif (step_sign_prev == -step_sign) then 
                    bracket_success = 1 
                    brnt_xbracket = [infl_A, infl_B] 
                    brnt_ybracket = [excess_b, excess_b_B]
                endif
                    
            enddo 
             
            call initBrent(brnt_xbracket, brnt_ybracket, brnt_xnew, brnt_ynew, brnt_xvec, brnt_yvec, brnt_mflag, brnt_sflag, brnt_iter)
            
            do while (brnt_sflag ==0) 

            call updateBrent(brnt_xnew, brnt_ynew, brnt_xvec, brnt_yvec, brnt_mflag, brnt_sflag, brnt_iter, 1E-14_dp, 1E-14_dp, 1000)

            next_state        = next_state_irf_mat(:,:,sss)
            k_next            = k_next_irf_mat(:,sss) 
            next_zf          = next_zf_mat(:,sss) 
            nxt_omg           = next_omg_mat(:,sss) 
            nom_i_vec         = nom_i_irf_mat(:,sss)
            s                 = s_irf_mat(sss)
            q_current         = q_irf_mat(sss)
            c_spending_vec    = c_spending_irf_mat(:,sss)
            share_vec         = share_irf_mat(:,sss)
            bF_share_vec    = bF_share_irf_mat(:,sss)
            infl_vec          = infl_irf_mat(:,sss)
            l_aggr_vec        = l_aggr_irf_mat(:,sss)
            inv_h             = inv_h_irf_mat(sss)
            inv_f             = inv_f_irf_mat(sss)

            call update_ifixed(brnt_xnew, ih_fixed, next_state, n_nxt, sss, smol_coeffs, k_next, next_zf, nxt_omg, nom_i_vec, s, q_current, & 
                         c_spending_vec, share_vec, bF_share_vec, infl_vec, l_aggr_vec, inv_h, inv_f, & 
                         bond_supply_shock, bond_supply_shock_nxt, & 
                         brnt_ynew, v_new, mc_new, w_choice_new & 
                         )
            enddo 
            
            if (brnt_sflag > 2) then 
                write(*,*) 'Brent failed to converge'
                stop
            endif


            ! store solution 
            v_irf_mat(:,sss)           = v_new 
            mc_irf_mat(:, sss)         = mc_new
            s_irf_mat(sss)             = s
            c_spending_irf_mat(:,sss)  = c_spending_vec
            k_next_irf_mat(:,sss)      = k_next/exp(dz_vec_adj)*exp(dz_vec)
            q_irf_mat(sss)             = q_current
            inv_f_irf_mat(sss)         = inv_f 
            inv_h_irf_mat(sss)         = inv_h 
            l_aggr_irf_mat(:,sss)      = l_aggr_vec
            nom_i_irf_mat(:,sss)       = [ih_fixed, nom_i_vec(2)]
            w_choice_irf_mat(:,sss)    = w_choice_new
            share_irf_mat(:,sss)       = share_vec          
            bF_share_irf_mat(:,sss)  = bF_share_vec     
            infl_irf_mat(:,sss)        = [brnt_xnew, infl_vec(2)]
            next_state_irf_mat(:,:,sss)= next_state

        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        endif

        q_bond_irf_mat_tmp = q_bond_irf_mat
        
        bbb_counter = 0   
        do bbb = 1,n_bond
            if (bbb <= 4 .or. bbb > n_bond-2) then 
                bbb_counter = bbb_counter + 1
            endif 
        
        ! solve for polynomial coefficients
        ! call omp_set_num_threads(4)
        call F07ABF( 'Equilibration','No transpose',n_states, nrhs2, smol_polynom, lda,af,ldaf,ipiv,    &
                      equed,r,c, q_bond_irf_mat_tmp(:,:,bbb),ldb,bond_coeffs,ldx,rcond,ferr2,berr2,work,iwork,info)
        
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED( & 
        !$OMP next_state_irf_mat, smol_elem_ani, smol_coeffs, state_grid, &
        !$OMP nom_i_irf_mat, s_irf_mat, w_choice_irf_mat, n_states, dz_vec, & 
        !$OMP k_grid_dev, k_grid_mean, zf_grid_dev, zf_grid_mean, q_irf_mat, &
        !$OMP c_spending_irf_mat, share_irf_mat, infl_irf_mat, varsigma_vec, &  
        !$OMP sigma, bbb_counter, & 
        !$OMP bF_share_irf_mat, l_aggr_irf_mat, next_state_irf_mat_new, & 
        !$OMP inv_h_irf_mat, inv_f_irf_mat, next_omg_mat, &
        !$OMP bF_share_irf_mat_new, share_irf_mat_new, bond_supply_shock, bond_supply_shock_nxt, &
        !$OMP wh_grid_mean, wh_grid_dev, wf_grid_mean, wf_grid_dev, n_active_dims, &
        !$OMP n_quad, max_smol_level, &
        !$OMP w_choice_irf_mat_new, infl_irf_mat_new, &
        !$OMP q_irf_mat_new, l_aggr_irf_mat_new, v_irf_mat_new, s_irf_mat_new, nom_i_irf_mat_new, &
        !$OMP c_spending_irf_mat_new, k_next_irf_mat_new, dz_vec_adj, next_zf_mat, &
        !$OMP results_irf_mat, k_next_irf_mat, bbb, bond_coeffs, q_bond_irf_mat) &
        !$OMP PRIVATE( &
        !$OMP polyn_points, nxt_mat, inv_h, inv_f, E_rb, &
        !$OMP nom_i_vec, s, nxt_mat_2, counter, who_prices_bond, &
        !$OMP q_current, c_spending_vec, share_vec, nxt_bond_prices, &
        !$OMP bF_share_vec, infl_vec , l_aggr_vec, &
        !$OMP q_new, s_new, k_next_new, c_spending_vec_new, infl_new, &
        !$OMP sss, v_new, l_aggr_new, w_choice_new, results_vec, valuation_vec)
        !$OMP DO SCHEDULE(static)
        do sss = 1, n_states 
        

            ! interpolate future values 
            polyn_points = Smolyak_Polynomial2(next_state_irf_mat(:,:,sss), n_active_dims, &
                                      n_quad, n_states, max_smol_level, smol_elem_ani)
                                                                                         
            ! use MKL Library for multiplication, same as
            ! nxt_mat        = matmul(polyn_points, smol_coeffs)
            CALL DGEMM('N','N', n_quad, nrhs, n_states, 1.0_dp, polyn_points, n_quad, & 
                               smol_coeffs, n_states, 0.0_dp, & 
                               nxt_mat, n_quad)   

            CALL DGEMM('N','N', n_quad, nrhs2, n_states, 1.0_dp, polyn_points, n_quad, & 
                               bond_coeffs, n_states, 0.0_dp, & 
                               nxt_bond_prices, n_quad)   

            ! non interpolated values
            nxt_mat_2(:, 1) = k_next_irf_mat(:,sss) 
            nxt_mat_2(:, 2) = next_zf_mat(:,sss) 
            nxt_mat_2(:, 3) = next_omg_mat(:,sss) 

            ! current state and guesses 
            nom_i_vec         = nom_i_irf_mat(:,sss)
            s                 = s_irf_mat(sss)
            q_current         = q_irf_mat(sss)
            c_spending_vec    = c_spending_irf_mat(:,sss)
            share_vec         = share_irf_mat(:,sss)
            bF_share_vec    = bF_share_irf_mat(:,sss)
            infl_vec          = infl_irf_mat(:,sss)
            l_aggr_vec        = l_aggr_irf_mat(:,sss)
            inv_h             = inv_h_irf_mat(sss)
            inv_f             = inv_f_irf_mat(sss)
            
            call calc_bond_prices( nxt_mat, n_nxt, nxt_mat_2, nxt_bond_prices, sss, &
                        c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                        nom_i_vec, share_vec, bF_share_vec, inv_h, inv_f, & 
                        bond_supply_shock, bond_supply_shock_nxt, & 
                        who_prices_bond, q_bond_irf_mat(sss,:,bbb+1), E_rb, results_vec)
        
            ! stuff stored only for results
            results_irf_mat(sss,1:n_interp-6*12)  = results_vec
            if (bbb <= 4 .or. bbb > n_bond-2) then 
            results_irf_mat(sss,n_interp-6*12 + 1 + (bbb_counter-1)*3:n_interp-6*12+ 3*bbb_counter) = q_bond_irf_mat(sss,:,bbb+1)
            results_irf_mat(sss,n_interp-6*9  + 1 + (bbb_counter-1)*6:n_interp-6*9 + 6*bbb_counter) = E_rb
            results_irf_mat(sss,n_interp-6*3  + 1 + (bbb_counter-1)*3:n_interp-6*3 + 3*bbb_counter) = who_prices_bond
            endif 
        
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    enddo    

    if (ttt == 1001) then

        bond_supply_shock = 0.0_dp
        bond_supply_shock_nxt = 0.0_dp

        trans_shock_vec = shock_grid(no_shock_idx(1),:)
        if (irf_indices_ll(fff) == sidx_bg) then 
            trans_shock_vec(sidx_omg) = irf_ll_shock_sizes(fff)
        endif
        next_state_irftransition_mat     = next_state_mat(no_shock_idx,:,:)
        next_state_irftransition_mat_new = next_state_mat(no_shock_idx,:,:)

        do iii = 1,n_I
        interp_input_mat(:,iii)     = v_irf_mat(iii,:)
        interp_input_mat(:,n_I + iii)     = mc_irf_mat(iii,:)
        enddo

        interp_input_mat(:,2*n_I + 1)  = s_irf_mat
        interp_input_mat(:,2*n_I + 2)  = q_irf_mat
        interp_input_mat(:,2*n_I + 3)  = l_aggr_irf_mat(1,:)
        interp_input_mat(:,2*n_I + 4)  = l_aggr_irf_mat(2,:)
        interp_input_mat(:,2*n_I + 5)  = infl_irf_mat(1,:)
        interp_input_mat(:,2*n_I + 6)  = infl_irf_mat(2,:)
        interp_input_mat(:,2*n_I + 7) = c_spending_irf_mat(1,:)
        interp_input_mat(:,2*n_I + 8) = c_spending_irf_mat(2,:)
        
        ! solve for polynomial coefficients
        call F07ABF( 'Equilibration','No transpose',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,    &
                      equed,r,c,interp_input_mat,ldb,smol_coeffs,ldx,rcond,ferr,berr,work,iwork,info)
        
        diff = 1.0_dp
        outer_iter = 0
        do while (diff > conv_crit)
        outer_iter = outer_iter + 1

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED( & 
        !$OMP next_state_irftransition_mat, smol_elem_ani, smol_coeffs, state_grid, vector_mus_dimensions, &
        !$OMP nom_i_mat, s_mat, w_choice_mat, n_states, dz_vec, bond_supply_shock, & 
        !$OMP k_grid_dev, k_grid_mean, zf_grid_dev, zf_grid_mean, q_mat, next_omg_mat, &
        !$OMP c_spending_mat, share_mat, infl_mat, varsigma_vec, sigma, & 
        !$OMP bF_share_mat, l_aggr_mat, next_state_irftransition_mat_new, inv_h_mat, inv_f_mat, &
        !$OMP bF_share_mat_new, share_mat_new, disast_p, irf_ll_shock_sizes, &
        !$OMP wh_grid_mean, wh_grid_dev, wf_grid_mean, wf_grid_dev, n_active_dims, next_dis_mat, &
        !$OMP omg_grid_mean, omg_grid_dev, dis_grid_mean, dis_grid_dev, no_shock_idx, bond_supply_shock_nxt, &
        !$OMP tht_h_grid_mean, tht_h_grid_dev, max_smol_level, &
        !$OMP w_choice_mat_new, infl_mat_new, fff, omg_shift, &
        !$OMP q_mat_new, l_aggr_mat_new, v_mat_new, s_mat_new, nom_i_mat_new, trans_shock_vec, &
        !$OMP c_spending_mat_new, k_next_mat_new, dz_vec_adj, next_zf_mat, &
        !$OMP k_next_mat) &
        !$OMP PRIVATE( &
        !$OMP polyn_points_trans, nxt_mat_trans, inv_h, inv_f, &
        !$OMP nom_i_vec, s, nxt_mat_2_trans, counter, &
        !$OMP q_current, c_spending_vec, share_vec, &
        !$OMP bF_share_vec, infl_vec , l_aggr_vec, &
        !$OMP q_new, s_new, k_next_new, c_spending_vec_new, infl_new, &
        !$OMP sss, v_new, l_aggr_new, w_choice_new, theta_nxt)
        !$OMP DO SCHEDULE(static)
        do sss = 1, n_states 
        

            ! interpolate future values 
            polyn_points_trans = Smolyak_Polynomial2(next_state_irftransition_mat(1,:,sss), n_active_dims, &
                                      1, n_states, max_smol_level, smol_elem_ani)
                                                                
            ! use MKL Library for multiplication, same as
            ! nxt_mat        = matmul(polyn_points, smol_coeffs)
            CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points_trans, 1, & 
                               smol_coeffs, n_states, 0.0_dp, & 
                               nxt_mat_trans, 1)   

            ! non interpolated values 
            nxt_mat_2_trans(1, 1) = k_next_mat(no_shock_idx(1),sss) 
            nxt_mat_2_trans(1, 2) = next_zf_mat(no_shock_idx(1),sss) 
            nxt_mat_2_trans(1, 3) = next_omg_mat(no_shock_idx(1),sss) 

            ! current state and guesses 
            nom_i_vec         = nom_i_mat(:,sss)
            s                 = s_mat(sss)
            q_current         = q_mat(sss)
            c_spending_vec    = c_spending_mat(:,sss)
            share_vec         = share_mat(:,sss)
            bF_share_vec    = bF_share_mat(:,sss)
            l_aggr_vec        = l_aggr_mat(:,sss)
            infl_vec          = infl_mat(:,sss)
            
            call calc_unexpected_transition( nxt_mat_trans, n_nxt, nxt_mat_2_trans, sss, &
                        c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                        nom_i_vec, share_vec, bF_share_vec, trans_shock_vec,   &
                        bond_supply_shock, bond_supply_shock_nxt, theta_nxt(1))

            w_choice_new = w_choice_mat(:,sss)

        counter = 0
        
        if (vector_mus_dimensions(1) > 0) then
            counter = counter + 1
            next_state_irftransition_mat_new(1, counter, sss)   = (k_next_mat(1,sss)/exp(trans_shock_vec(sidx_z)) -k_grid_mean )/k_grid_dev
        endif
        
        if (vector_mus_dimensions(2) > 0) then
            counter = counter + 1
            next_state_irftransition_mat_new(:, counter, sss) = (theta_nxt(1) - tht_h_grid_mean)/tht_h_grid_dev
        endif

        if (vector_mus_dimensions(3) > 0) then
        counter = counter + 1
        endif
    
        if (vector_mus_dimensions(4) > 0) then
            counter = counter + 1
            next_state_irftransition_mat_new(1, counter, sss)  = (w_choice_new(1)/exp(trans_shock_vec(sidx_z))/ & 
                                                          (varsigma_vec(1)  & 
                                                          + (1.0-varsigma_vec(1) )*s**(sigma-1.0))**(1.0/(sigma-1.0)) & 
                                                          -wh_grid_mean) /wh_grid_dev
        endif
    
        if (vector_mus_dimensions(5) > 0) then
            counter = counter + 1
            next_state_irftransition_mat_new(1, counter, sss)  = (w_choice_new(2)/exp(trans_shock_vec(sidx_z))*s/ & 
                                                          (varsigma_vec(2)*s**(1.0_dp-sigma) + (1.0-varsigma_vec(2)))**(1.0/(sigma-1.0)) & 
                                                          -wf_grid_mean) /wf_grid_dev
        endif
    
        if (vector_mus_dimensions(6) > 0) then
            counter = counter + 1
            if (dis_grid_dev < sqrt_eps) then
            next_state_irftransition_mat_new(1,counter,sss) = 0.0_dp  ! + trans_shock_vec(sidx_dis)/dis_grid_dev
            else
            next_state_irftransition_mat_new(:,counter,sss) = (next_dis_mat(no_shock_idx,sss) + trans_shock_vec(sidx_p)   - dis_grid_mean) /dis_grid_dev
            endif
        endif
    
        if (vector_mus_dimensions(7) > 0) then
            counter = counter + 1
        endif
    
        if (vector_mus_dimensions(8) > 0) then
            counter = counter + 1
        endif
       
        if (vector_mus_dimensions(9) > 0) then
            counter = counter + 1
            if (omg_grid_dev < sqrt_eps) then
                next_state_irftransition_mat_new(1,counter,sss) = 0.0_dp ! trans_shock_vec(sidx_omg)/omg_grid_dev
            else
                next_state_irftransition_mat_new(:,counter,sss) = (next_omg_mat(no_shock_idx,sss) + trans_shock_vec(sidx_omg)   - omg_grid_mean) /omg_grid_dev
            endif
        endif
    

        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    
    where(next_state_irftransition_mat_new> 1.0_dp) 
        next_state_irftransition_mat_new= 1.0_dp 
    endwhere
    where(next_state_irftransition_mat_new < -1.0_dp)
        next_state_irftransition_mat_new   = -1.0_dp
    endwhere
    
    if (outer_iter < 10) then
        diff = 1.0_dp
    else
        diff =  maxval(abs(next_state_irftransition_mat - next_state_irftransition_mat_new)) 
    endif

    next_state_irftransition_mat = next_state_irftransition_mat +  0.2_dp*(next_state_irftransition_mat_new-next_state_irftransition_mat)

    enddo
    endif


    write(iter_char, '(i4)') ttt
    write(shock_char2,'(i2)') 10+fff 
        
    ! print all results
    open (unit = 10, file = trim(results_folder) // 'next_state_shock_mat_' // shock_char2 // '_'  // iter_char  // '.dat', ACTION="write", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    write(10) next_state_irf_mat; close(10)

    if (ttt == 1001) then
    ! print all results
    open (unit = 10, file = trim(results_folder) // 'next_state_shocktrans_mat_' // shock_char2  // '_'  // iter_char  // '.dat', ACTION="write", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    write(10) next_state_irftransition_mat; close(10)
    endif

    ! print all results
    open (unit = 10, file = trim(results_folder) // 'results_shock_mat_' // shock_char2  // '_'  // iter_char  // '.dat', ACTION="write", STATUS="replace", &
            & FORM="unformatted", ACCESS="STREAM")
    write(10) results_irf_mat; close(10)



    enddo ! ttt
enddo
write(*,*) 'DONE'
flush(6)
endif
end subroutine calc_sol


subroutine calc_unexpected_transition( nxt_mat, n_nxt, nxt_mat_2, sss, &
                                       c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                                       nom_i_vec, share_vec, bF_share_vec, trans_shock_vec, &
                                       bond_supply_shock, bond_supply_shock_nxt, theta_nxt)
    
    use mod_param, only: n_I, idx_zf, idx_k, aalpha, smolyak_d, quad_weight_vec, &
                         ddelta, state_grid, varsigma_vec, sigma, idx_dis, b_lmbd, omg_shift, &
                         chix, ies_vec, gma_vec, disast_p, dz_vec, &
                         bbeta_vec, v_normalization_vec, chi0_vec, chi, wealth_share_grid, idx_omg, &
                         n_interp,  idx_wh, idx_wf, phi_h, phi_f, &
                         sidx_z, n_shocks, bg_yss,omg_grid_mean, & 
                         zeta, bbeta_coeff, idx_thh, tht_h_grid_mean, phi_yh, phi_yf

                        
    real(dp), intent(in) :: nxt_mat_2(1, 3), l_aggr_vec(2), trans_shock_vec(n_shocks), &
                            nxt_mat(1, n_nxt), s, bond_supply_shock, bond_supply_shock_nxt, &
                            q_current, c_spending_vec(n_I) , infl_vec(2)
    real(dp), intent(inout) :: bF_share_vec(n_I), share_vec(n_I), nom_i_vec(2)
    real(dp), intent(out)   ::  theta_nxt(1)

    integer, intent(in)  :: sss, n_nxt

    integer :: iii, ggg, ifail, nnn

    real(dp) :: s_vec(2), zf_vec(2), k_aggr, kappa_vec(2), & 
                q_nxt(1), rk_vec(1), dz_temp, aggr_wealth, &
                w_next_choice(1, 2), pi_nxt(1), y_next(1, 2), &
                kappa_nxt(1,2), homegood_infl_nxt(1, 2), &
                k_nxt(1), s_nxt(1), zf_nxt(1,2), omg_nxt(1), &
                rf_vec(1, 2), pi_current, y_current(2), &
                excess_b, k_temp_vec(n_I), b_temp_vec(n_I), c_temp(n_I), & 
                share_temp(n_I), bF_share_temp(n_I), infl_nxt(1,2), &
                c_vec(n_I), ch_vec(n_I), chi0, savings_vec(n_I)

    real(dp) :: next_k

    real(dp), parameter :: brent_delta = 5E-16_dp, min_cons_sav = 1E-08_dp
 
    real(dp) :: &
             w_choice_vec(2), &
             wealth, consumption, r_alpha_vec(1), &
             diff_inner, savings, objf, objgrd(2), cons_update, IES, bbeta, &
             v_normalization, v_vec(1), EV, &
             M_vec_mat(1, n_I), w_temp, w_diff, w_temp_new, temp_vec(1), &
             share_low, share_high, real_bond, c_spend_temp(n_I), &
             mid_point, k_param, derivative_consumption, M_weighted_vec(1), &
             share_guess, gov_debt, wealth_tax_income, gov_surplus, lump_sum_taxes(n_I), gma, &
             ch, cf, c_spend, labor, c_cost, c_cost_nxt(1), bbeta_adj(n_I), & 
             v_vec_twisted(1), iota(1), & 
             l_aggr_nxt(1,2), nxt_wealth_vec(1,n_I), seigniorage_nxt(1),&
             s_update, P_div_P_h(2), P_div_P_h_nxt(1,2), w_current_next(1), seignorage_transfer(1,n_I), &
             w_current(2), w_choice_new_vec(2), kappa_guess, kappa_share, util, util_c_deriv, labor_part, &
             l_aggr_vec_temp(2), l_aggr_temp(2), l_aggr_temp_new(2), w_choice_temp(2), w_choice_temp_new(2), &
             kappa_vec_new(2), p_dis,  ch_spending_nxt(1), cf_spending_nxt(1),  omg, & 
             E_rfh, E_rff, E_rk, investment, q_temp, test_ratio, tax_vec(1)
 
 integer :: iter_inner, mflag, ix1, ix2

! =============================================================================
! CURRENT STATE ===============================================================
! =============================================================================

    ! current capital stock
    k_aggr = state_grid(idx_k, sss)

    ! store relative productivity in vector
    zf_vec = [ 1.0_dp , zeta*exp(state_grid(idx_zf,sss)) ]

    ! current disaster probability
    p_dis   = exp(state_grid(idx_dis, sss))

    ! w_current is stored as w*P_h/P, so in bundle units
    ! allows to calculate w'/w*local_infl = w'/w_current*infl*P'_h/P'
    w_current = [state_grid(idx_wh, sss), state_grid(idx_wf, sss)]

! =============================================================================
! CURRENT PERIOD GUESSES ======================================================
! =============================================================================

    ! store relative prices in vector
    s_vec   = [1.0_dp, s]

    ! given current guess of aggregate labor and wage, allocate capital across countries
    ! allows to calculate initial wealth
    kappa_vec(1) = k_aggr * ( (s_vec(2) )**(1.0/(1.0-aalpha))) * &
                      zf_vec(1) / zf_vec(2) * l_aggr_vec(1)/l_aggr_vec(2) / & 
                   (1.0 +   ( (s_vec(2) )**(1.0/(1.0-aalpha))) * &
                      zf_vec(1) / zf_vec(2) * l_aggr_vec(1)/l_aggr_vec(2) )

    kappa_vec(2) = k_aggr - kappa_vec(1)

    ! given this calculate profits and aggregate wealth
    pi_current   = aalpha * (zf_vec(1)*l_aggr_vec(1)/kappa_vec(1))**(1.0-aalpha)
    
    
    y_current    = (zf_vec**(1.0-aalpha)) * kappa_vec**aalpha * l_aggr_vec**(1.0-aalpha)
    
    w_choice_vec(1) = (1.0-aalpha)*y_current(1)/l_aggr_vec(1)
    w_choice_vec(2) = s**(-1)*(1.0-aalpha)*y_current(2)/l_aggr_vec(2)

    ! cost of unit of bundle in terms of ch
    P_div_P_h(1) = (varsigma_vec(1) + (1.0 -varsigma_vec(1) )*s**(sigma-1))**(1/(1-sigma))
    P_div_P_h(2) = s**(-1) * (varsigma_vec(2)*s**(1.0-sigma) + (1.0 -varsigma_vec(2)))**(1/(1-sigma))
    
    aggr_wealth  = k_aggr*((1-ddelta)*q_current + pi_current) 

   
    ! set the omg - with rebating only appears in taylor rule
    omg      = exp(state_grid(idx_omg,sss)) + omg_shift  - b_lmbd*bond_supply_shock

    
    bbeta_adj(1) = - bbeta_coeff*(state_grid(idx_thh,sss) - tht_h_grid_mean)
    bbeta_adj(2) = + bbeta_coeff*(state_grid(idx_thh,sss) - tht_h_grid_mean)
       
! =============================================================================
! NEXT PERIOD =================================================================
! =============================================================================
    
    dz_temp = trans_shock_vec(sidx_z)

    ! tomorrow's price, labor choice, infl and price of time endowment
    s_nxt               = nxt_mat(:,2*n_I + 1)
    q_nxt               = nxt_mat(:,2*n_I + 2) 
    l_aggr_nxt(:,1)     = nxt_mat(:,2*n_I + 3)
    l_aggr_nxt(:,2)     = nxt_mat(:,2*n_I + 4)
    infl_nxt(:,1)       = nxt_mat(:,2*n_I + 5)
    infl_nxt(:,2)       = nxt_mat(:,2*n_I + 6)
    ch_spending_nxt     = nxt_mat(:,2*n_I + 7)*exp(dz_temp)
    cf_spending_nxt     = nxt_mat(:,2*n_I + 8)*exp(dz_temp)

    ! currently assumed capital in the next period
    ! vector because depreciation shock in disaster state 
    k_nxt      = nxt_mat_2(:,1)

    ! relative productivity next period
    zf_nxt(:,1)   = 1.0_dp 
    zf_nxt(:,2)   = zeta*exp(nxt_mat_2(:,2))
    
    omg_nxt      = exp(nxt_mat_2(:,3)) + omg_shift  - b_lmbd*bond_supply_shock_nxt

    homegood_infl_nxt(:,1) = infl_nxt(:,1)*( ( varsigma_vec(1)  & 
                                        + (1.0-varsigma_vec(1) )*    s**(sigma-1.0) ) / & 
                                          ( varsigma_vec(1)  & 
                                          + (1.0-varsigma_vec(1) )*s_nxt**(sigma-1.0) ) )**(1.0/(1.0-sigma)) 
    homegood_infl_nxt(:,2) = infl_nxt(:,2)* s_nxt/s *( ( varsigma_vec(2)*    s**(1.0-sigma)+ (1.0-varsigma_vec(2)) ) / & 
                                          ( varsigma_vec(2)*s_nxt**(1.0-sigma)+ (1.0-varsigma_vec(2)) ) )**(1.0/(1.0-sigma)) 

    ! given employed labor, allocate capital to equalize capital returns across countries
    kappa_nxt(:,1) = k_nxt *  ((s_nxt)**(1.0/(1.0-aalpha))) * zf_nxt(:,1) / zf_nxt(:,2) * l_aggr_nxt(:,1)/l_aggr_nxt(:,2) / & 
                       (1.0 + ((s_nxt)**(1.0/(1.0-aalpha))) * zf_nxt(:,1) / zf_nxt(:,2) * l_aggr_nxt(:,1)/l_aggr_nxt(:,2) )
    
    kappa_nxt(:,2) = k_nxt - kappa_nxt(:,1)

    next_k       = k_nxt(1)

    investment = next_k - k_aggr*(1.0-ddelta)

    ! output in next period
    y_next(:,1) = exp((log(zf_nxt(:,1)) + dz_temp)*(1.0_dp-aalpha))*l_aggr_nxt(:,1)* & 
                  ( (kappa_nxt(:,1)/l_aggr_nxt(:,1))**aalpha)
    y_next(:,2) = exp((log(zf_nxt(:,2)) + dz_temp)*(1.0_dp-aalpha))* l_aggr_nxt(:,2)* & 
                  ( (kappa_nxt(:,2)/l_aggr_nxt(:,2))**aalpha)
    
    w_next_choice(:,1) = (1-aalpha)*y_next(:,1)/l_aggr_nxt(:,1)
    w_next_choice(:,2) = s_nxt**(-1) * (1-aalpha)*y_next(:,2)/l_aggr_nxt(:,2)

    ! profit next period per unit of capital
    pi_nxt = aalpha*exp( (log(zf_nxt(:,1)) + dz_vec)*(1.0_dp-aalpha)) &
           *(l_aggr_nxt(:,1)**(1.0_dp-aalpha))*kappa_nxt(:,1)**(aalpha-1.0)
   
    ! capital return, adjustment in disaster state
    rk_vec = ( (1.0_dp-ddelta)*q_nxt + pi_nxt)/q_current
    rk_vec = ( (1.0_dp-ddelta)*q_nxt)/(q_current - pi_current)

    P_div_P_h_nxt(:,1) = (varsigma_vec(1)  + & 
                         (1.0 -varsigma_vec(1) )*s_nxt**(sigma-1))**(1.0/(1.0-sigma))
    P_div_P_h_nxt(:,2) = s_nxt**(-1) * (varsigma_vec(2)*s_nxt**(1.0-sigma) + (1.0 -varsigma_vec(2)))**(1.0/(1.0-sigma))

    ! calculate risk-free returns, given nominal rate guess 
    ! transform everything in units of home good
    rf_vec(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1) 
    rf_vec(:,2) = (nom_i_vec(1) + nom_i_vec(2))/homegood_infl_nxt(:,2) 

    ! seignorage: next period's omega creates seigniorage payout next period
    if (bg_yss > 0) then 
        seigniorage_nxt = cf_spending_nxt/(ch_spending_nxt + cf_spending_nxt)*(bg_yss + bond_supply_shock_nxt)*(ch_spending_nxt + cf_spending_nxt)*omg_nxt
        seignorage_transfer(:,1) =  seigniorage_nxt 
        seignorage_transfer(:,2) = -seigniorage_nxt
    else
    seignorage_transfer = 0.0_dp
    endif

!! =============================================================================
!! STORE VALUE FUNCTION, NEXT CAPITAL, NEXT WEALTH SHARES ======================
!! =============================================================================

    c_spend_temp = c_spending_vec
    do iii = 1, n_I

            ! read in preferece parameters for current agent
            IES   = IES_vec(iii)
            bbeta = bbeta_vec(iii) + bbeta_adj(iii)
            chi0   = chi0_vec(iii)
            v_normalization = v_normalization_vec(iii)
          
            c_cost     = P_div_P_h(iii)
            c_cost_nxt = P_div_P_h_nxt(:,iii)
          
            wealth     = aggr_wealth*wealth_share_grid(iii,sss) 
            
            labor      = l_aggr_vec(iii)


            c_spend =  minval( [maxval([c_spend_temp(iii), min_cons_sav]), & 
                       wealth + w_choice_vec(iii)*labor-min_cons_sav] )
            
            savings     = wealth + w_choice_vec(iii)*labor - c_spend
          
            consumption = c_spend/c_cost
            c_vec(iii) = consumption

            r_alpha_vec = (share_vec(iii)-bF_share_vec(iii))*rf_vec(:,1)/(1.0_dp - omg) & 
                                               + bF_share_vec(iii)*rf_vec(:,2) & 
                              + (1.0_dp - share_vec(iii))*rk_vec

            labor_part =  1.0_dp + (1.0_dp/IES - 1.0_dp) * chi0*chi/(1.0_dp + chi) * labor**( (1.0_dp+chi)/chi )
            nxt_wealth_vec(:,iii) = (savings*r_alpha_vec + seignorage_transfer(:,iii))

            savings_vec(iii) = savings
       
    enddo
      
    theta_nxt = nxt_wealth_vec(:,1)/sum(nxt_wealth_vec,2)
    
end subroutine calc_unexpected_transition

subroutine calc_equilibrium_and_update( nxt_mat, n_nxt, nxt_mat_2, sss, outer_iter, &
                            c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                            nom_i_vec, share_vec, bF_share_vec, inv_h, inv_f, bond_supply_shock, bond_supply_shock_nxt,  &
                            q_new, s_new, k_next_new, c_spending_vec_new, infl_new, &
                            v_new, mc_new, l_aggr_new, w_choice_new, theta_nxt, M_vec_mat)
    
    use mod_param, only:    n_I, n_quad,  idx_zf, idx_k, aalpha, smolyak_d, quad_weight_vec, zeta, bbeta_coeff, &
                            dz_vec, ddelta, state_grid, varsigma_vec, sigma, idx_dis, inv_share_h, omg_rho, idx_thh, &
                            chix, ies_vec, gma_vec, omg_shift, b_lmbd, bg_yss,&  
                            bbeta_vec, v_normalization_vec, chi0_vec, chi, wealth_share_grid, &
                            vareps_w, n_interp, &
                            phi_w, idx_wh, idx_wf, phi_h, phi_f, tayl_ic_h, tayl_ic_f, tht_h_grid_mean, & 
                            idx_ih, idx_if, rho_i, disast_p, disast_std, idx_omg, &
                            foreign_trading, omg_grid_mean, phi_yh, phi_yf 

    real(dp), intent(in) :: nxt_mat_2(n_quad, 3), l_aggr_vec(2), nxt_mat(n_quad, n_nxt), s, &
                            q_current, c_spending_vec(n_I), infl_vec(2), bond_supply_shock, bond_supply_shock_nxt 
    real(dp), intent(inout) :: bF_share_vec(n_I), share_vec(n_I), nom_i_vec(2), M_vec_mat(n_quad, n_I)
    real(dp), intent(out)   :: q_new, s_new, k_next_new, c_spending_vec_new(n_I), infl_new(2), inv_f, inv_h, &
                               v_new(n_I), l_aggr_new(2), w_choice_new(2), theta_nxt(n_quad), mc_new(n_I)

    integer, intent(in)  :: sss, n_nxt, outer_iter

    integer :: iii, ggg, ifail, nnn

    real(dp) :: s_vec(2), zf_vec(2), k_aggr, kappa_vec(2), q_nxt(n_quad), rk_vec(n_quad), &
                aggr_wealth, v_temp(n_quad, n_I),  mc_temp(n_quad, n_I), &
                w_next_choice(n_quad, 2), pi_nxt(n_quad), y_next(n_quad, 2), &
                kappa_nxt(n_quad,2), homegood_infl_nxt(n_quad, 2), &
                k_nxt(n_quad), s_nxt(n_quad), zf_nxt(n_quad,2), omg_nxt(n_quad), &
                rf_vec(n_quad, 2), pi_current(2), y_current(2), &
                excess_b, k_temp_vec(n_I), b_temp_vec(n_I), c_temp(n_I), & 
                excess_b_home, bf_temp_vec(n_I), bh_temp_vec(n_I), & 
                share_temp(n_I), bF_share_temp(n_I), infl_nxt(n_quad,2), &
                c_vec(n_I), ch_vec(n_I), chi0, savings_vec(n_I), &
                labor_nxt(n_quad)

    real(dp) :: next_k

    real(dp), parameter :: brent_delta = 5E-16_dp, min_cons_sav = 1E-08_dp
 
    real(dp) :: &
             nom_i_a, nom_i_b, rf_vec_temp(n_quad,2), excess_b_A, excess_b_B, w_choice_vec(2), &
             nom_i_C, temp, excess_b_C, excess_b_S, nom_i_S, nom_i_D, wealth, consumption, r_alpha_vec(n_quad), &
             diff_inner, iter_inner, savings, objf, objgrd(2), cons_update, IES, bbeta, &
             v_normalization, next_period_share(n_quad), v_vec(n_quad), EV, &
             w_temp, w_diff, w_temp_new, temp_vec(n_quad), share_low, share_high, real_bond, c_spend_temp(n_I), &
             mid_point, k_param, derivative_consumption, M_weighted_vec(n_quad), nom_i_guess, &
             share_guess, gov_debt,  gma, bbeta_adj(n_I), ch, cf, c_spend, labor, c_cost, c_cost_nxt(n_quad), & 
             v_vec_twisted(n_quad), ih_last, if_last, iota(n_quad), l_aggr_nxt(n_quad,2), nxt_wealth_vec(n_quad,n_I), &
             s_update, diff, P_div_P_h(2), P_div_P_h_nxt(n_quad,2), w_current_next(n_quad), &
             w_current(2), w_choice_new_vec(2), kappa_guess, kappa_share, util, util_c_deriv, labor_part, &
             l_aggr_vec_temp(2), l_aggr_temp(2), l_aggr_temp_new(2), w_choice_temp(2), w_choice_temp_new(2), &
             kappa_vec_new(2), p_dis, big_weight_vec(n_quad), deriv_helper(n_quad), labor_new, &
             E_rfh, E_rff, E_rk, investment, q_temp, test_ratio, l_temp_vec(n_I), & 
             temp_bond_holdings, idiff, nom_i_vec_in(2), omg, r_alpha_omg_vec(n_quad), &
             ch_spending_nxt(n_quad), cf_spending_nxt(n_quad), seigniorage_nxt(n_quad),& 
             tax_vec(n_quad), seignorage_transfer(n_quad, n_I), wealth_vec(n_I)

 
 integer :: iter, mflag, ix1, ix2, safe_constraint_binding_vec(n_I), constraint_binding_vec(n_I), bond_pricer, iiter


    nom_i_vec_in = nom_i_vec ! store incoming nominal rate guess for Taylor rule

! =============================================================================
! CURRENT STATE ===============================================================
! =============================================================================

    ! current capital stock
    k_aggr        = state_grid(idx_k,   sss)
    
    ! store relative productivity in vector, foreign with relative productivity shock 
    zf_vec = [ 1.0_dp, zeta*exp(state_grid(idx_zf, sss))]

    ! current disaster probability
    p_dis       = exp(state_grid(idx_dis, sss))
    

    ! current weight vec, given disaster probability
    big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
    big_weight_vec(n_quad)     = p_dis

    ! w_current is stored as w / (P /P_h) and w / (P* / P_h*), so in bundle units
    ! allows to calculate w'/w*P'_h/P_h = w'/w_current*infl*P'_h/P'
    w_current = [state_grid(idx_wh, sss), state_grid(idx_wf, sss)]

    ih_last = state_grid(idx_ih,sss)
    if_last = state_grid(idx_if,sss)
   

    ! theta directly used in next step 

! =============================================================================
! STEP 1: CURRENT PERIOD GUESSES ==============================================
! =============================================================================

    ! store relative prices home/foreign good in vector
    s_vec   = [1.0_dp, s]

    ! given current guess of aggregate labor and wage, allocate capital across countries
    ! allows to calculate initial wealth 
    kappa_vec(1) = k_aggr * ( (s_vec(2) )**(1.0/(1.0-aalpha))) * &
                      zf_vec(1) / zf_vec(2) * l_aggr_vec(1)/l_aggr_vec(2) / & 
                   (1.0 +   ( (s_vec(2) )**(1.0/(1.0-aalpha))) * &
                      zf_vec(1) / zf_vec(2) * l_aggr_vec(1)/l_aggr_vec(2) )

    kappa_vec(2) = k_aggr - kappa_vec(1)

    ! given this calculate profits and aggregate wealth
    pi_current   = (s_vec**(-1))*aalpha * (zf_vec*l_aggr_vec/kappa_vec)**(1.0-aalpha)
   
    ! current output 
    y_current    = (zf_vec**(1.0-aalpha)) * kappa_vec**aalpha * l_aggr_vec**(1.0-aalpha)
   
    ! wages calculated in home consumption goods
    w_choice_vec = (s_vec**(-1))*(1.0-aalpha)*y_current/l_aggr_vec

    ! cost of unit of bundle in terms of home good (P/P_h and P*/P^*_h 
    ! (one bundle unit valued at P, divided by cost of home good -> units of home good per bundle
    P_div_P_h(1) = (    varsigma_vec(1)  & 
                        + (1.0 - varsigma_vec(1) )*(s**(sigma-1.0))  )**(1.0/(1.0-sigma))
    P_div_P_h(2) = s**(-1) * (varsigma_vec(2)*s**(1.0-sigma) + (1.0 -varsigma_vec(2)))**(1.0/(1.0-sigma))
    
    aggr_wealth  = kappa_vec(1)*( (1-ddelta)*q_current + pi_current(1) ) +  &
                   kappa_vec(2)*( (1-ddelta)*q_current + pi_current(2) )
    wealth_vec = aggr_wealth * wealth_share_grid(:, sss)               


    ! set the omg - with rebating only appears in taylor rule
    omg      = exp(state_grid(idx_omg,sss)) + omg_shift  - b_lmbd*bond_supply_shock
    
    ! beta adjustment for stationarity if bbeta_coeff > 0 
    bbeta_adj(1) = - bbeta_coeff*(state_grid(idx_thh,sss) - tht_h_grid_mean)
    bbeta_adj(2) = + bbeta_coeff*(state_grid(idx_thh,sss) - tht_h_grid_mean)

! =============================================================================
! STEP 2: NEXT PERIOD =========================================================
! =============================================================================

    ! value functions
    do iii = 1,n_I
        v_temp(:,iii)  = nxt_mat(:, iii )
        mc_temp(:,iii) = nxt_mat(:, n_I + iii)
    enddo

    ! tomorrow's price, labor choice, infl and price of time endowment
    s_nxt               = nxt_mat(:,2*n_I + 1)
    q_nxt               = nxt_mat(:,2*n_I + 2) 
    l_aggr_nxt(:,1)     = nxt_mat(:,2*n_I + 3)
    l_aggr_nxt(:,2)     = nxt_mat(:,2*n_I + 4)
    infl_nxt(:,1)       = nxt_mat(:,2*n_I + 5)
    infl_nxt(:,2)       = nxt_mat(:,2*n_I + 6)
    ch_spending_nxt     = nxt_mat(:,2*n_I + 7)*exp(dz_vec)
    cf_spending_nxt     = nxt_mat(:,2*n_I + 8)*exp(dz_vec)

    ! currently assumed capital in the next period (tilde bar k * exp(varphi_+1)
    ! vector because depreciation shock in disaster state 
    k_nxt      = nxt_mat_2(:,1)

    ! chosen capital, absent depreciation shock (tilde bar k)
    next_k         = k_nxt(1) 

    ! relative productivity next period
    zf_nxt(:,1)   = 1.0_dp 
    zf_nxt(:,2)   = zeta*exp(nxt_mat_2(:,2))

    omg_nxt      = exp(nxt_mat_2(:,3)) + omg_shift  - b_lmbd*bond_supply_shock_nxt

    ! local inflation P'_h/P_h = (P/P_h)/(P'/P'_h) * P'/P and equivalently P'_f/P_f
    homegood_infl_nxt(:,1) = infl_nxt(:,1)*( &  
                           ( varsigma_vec(1)  + (1.0-varsigma_vec(1))*    s**(sigma-1.0) ) / & 
                           ( varsigma_vec(1)  + (1.0-varsigma_vec(1))*s_nxt**(sigma-1.0) )  )**(1.0/(1.0-sigma)) 
    homegood_infl_nxt(:,2) = infl_nxt(:,2)* s_nxt/s *( ( varsigma_vec(2)*    s**(1.0-sigma)+ (1.0-varsigma_vec(2)) ) / & 
                                           ( varsigma_vec(2)*s_nxt**(1.0-sigma)+ (1.0-varsigma_vec(2)) ) )**(1.0/(1.0-sigma)) 

    ! given employed labor, allocate capital to equalize capital returns across countries
    kappa_nxt(:,1) = k_nxt *  ((s_nxt)**(1.0/(1.0-aalpha))) * zf_nxt(:,1) / zf_nxt(:,2) * l_aggr_nxt(:,1)/l_aggr_nxt(:,2) / & 
                       (1.0 + ((s_nxt)**(1.0/(1.0-aalpha))) * zf_nxt(:,1) / zf_nxt(:,2) * l_aggr_nxt(:,1)/l_aggr_nxt(:,2) )

    kappa_nxt(:,2) = k_nxt - kappa_nxt(:,1)



    ! output in next period
    y_next(:,1) = exp((log(zf_nxt(:,1)) + dz_vec)*(1.0_dp-aalpha))*l_aggr_nxt(:,1)* & 
                  ( (kappa_nxt(:,1)/l_aggr_nxt(:,1))**aalpha)
    y_next(:,2) = exp((log(zf_nxt(:,2)) + dz_vec)*(1.0_dp-aalpha))* l_aggr_nxt(:,2)* & 
                  ( (kappa_nxt(:,2)/l_aggr_nxt(:,2))**aalpha)
    
   
    w_next_choice(:,1) = (1-aalpha)*y_next(:,1)/l_aggr_nxt(:,1)
    w_next_choice(:,2) = s_nxt**(-1) * (1-aalpha)*y_next(:,2)/l_aggr_nxt(:,2)

    ! profit next period per unit of capital
        pi_nxt = aalpha*exp( (log(zf_nxt(:,1)) + dz_vec)*(1.0_dp-aalpha)) &
               *(l_aggr_nxt(:,1)**(1.0_dp-aalpha))*kappa_nxt(:,1)**(aalpha-1.0)

    ! capital return, adjustment in disaster state
    rk_vec = ( (1.0_dp-ddelta)*q_nxt + pi_nxt)/q_current
    rk_vec(n_quad) =exp(dz_vec(n_quad))*rk_vec(n_quad) 

    ! cost of unit of bundle in terms of home good (P/P_h and P*/P^*_h 
    ! (one bundle unit valued at P, divided by cost of home good -> units of home good per bundle
    P_div_P_h_nxt(:,1) = (varsigma_vec(1)  & 
                       + (1.0 -varsigma_vec(1)  )*s_nxt**(sigma-1))**(1.0/(1.0-sigma))
    P_div_P_h_nxt(:,2) = s_nxt**(-1) * (varsigma_vec(2)*s_nxt**(1.0-sigma) + (1.0 -varsigma_vec(2)))**(1.0/(1.0-sigma))

    ! calculate risk-free returns, given nominal rate guess 
    ! transform everything in units of home good
    rf_vec(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1) 
    rf_vec(:,2) = (nom_i_vec(1) + nom_i_vec(2))/homegood_infl_nxt(:,2) 
   
    ! seignorage
    if (bg_yss > 0) then 
        seigniorage_nxt = cf_spending_nxt/(ch_spending_nxt + cf_spending_nxt)*(bg_yss + bond_supply_shock_nxt)*(ch_spending_nxt + cf_spending_nxt)*omg_nxt
        seignorage_transfer(:,1) =  seigniorage_nxt 
        seignorage_transfer(:,2) = -seigniorage_nxt
    else
        seignorage_transfer = 0.0_dp
    endif

! =============================================================================
! STEP 3: UPDATE THE PRICE OF CAPITAL 
! =============================================================================

    ! update price of capital but only for next iteration
    q_new =  inv_share_h**(1.0/(1.0-sigma)) * (1 + &
              (1.0-inv_share_h)/inv_share_h*(s)**(sigma-1.0))**(1.0/(1.0-sigma)) & 
              * (next_k/state_grid(idx_k,sss))**chiX

! =============================================================================
! STEP 4: SOLVE FOR NOMINAL RATE IN EQUILIBRIUM ===============================
! =============================================================================
    
    ! at this point inflation and consumption choices are held fixed
    ! excess demand is solved for aggregate bond portfolios for both countries
    ! holding fixed within country composition, which is adjusted after
    if (outer_iter > 10) then

    ! check whether it generates excess demand or not
    c_spend_temp    = c_spending_vec
    c_temp          = c_spend_temp/P_div_P_h
    
    bF_share_temp = bF_share_vec
    share_temp      = share_vec
    
    rf_vec_temp(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1) 
    rf_vec_temp(:,2) = (nom_i_vec(1) + nom_i_vec(2))/homegood_infl_nxt(:,2) 

    call calc_excess_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, &
                              bbeta_adj,  P_div_P_h, P_div_P_h_nxt,   &  
                              q_current, w_choice_vec, l_aggr_vec, sss, big_weight_vec, &  
                              c_temp, c_spend_temp, share_temp, bF_share_temp, k_temp_vec, b_temp_vec, excess_b, constraint_binding_vec, mc_temp)

                                
    nom_i_guess = nom_i_vec(1)
    
    ! store
    nom_i_A    = nom_i_guess
    excess_b_A = excess_b

    ! find bound between solver finds solution for nominal rate
    if (excess_b == 0) then

      excess_b_S = excess_b

    else
      if (excess_b > 0.0_dp) then

        iter = 0
        do while (excess_b > 0.0_dp .and. iter <= 5000)
            iter = iter + 1

            ! initial guess - last period nominal rate
            nom_i_guess = nom_i_guess - 1E-02_dp

            ! check whether it generates excess demand or not
            rf_vec_temp(:,1) = nom_i_guess/homegood_infl_nxt(:,1)
            rf_vec_temp(:,2) = (nom_i_guess + nom_i_vec(2))/homegood_infl_nxt(:,2) 

            share_temp = share_vec
            call calc_excess_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, &
                            bbeta_adj, P_div_P_h, P_div_P_h_nxt,   &  
                            q_current, w_choice_vec, l_aggr_vec, sss, big_weight_vec, &  
                            c_temp, c_spend_temp, share_temp, bF_share_temp, k_temp_vec, b_temp_vec, excess_b, constraint_binding_vec, mc_temp)

        enddo
      
        if (iter >= 999) then
          write(*,*) 'No convergence in finding nom_i_low', sss
          stop
        end if

        nom_i_B = nom_i_guess
        excess_b_B = excess_b

      else

        iter = 0
        do while (excess_b < 0.0_dp .and. iter <= 5000)
            
            iter = iter + 1

            ! initial guess - last period nominal rate
            nom_i_guess  = nom_i_guess + 1E-02_dp

            rf_vec_temp(:,1) = nom_i_guess/homegood_infl_nxt(:,1) 
            rf_vec_temp(:,2) = (nom_i_guess + nom_i_vec(2))/homegood_infl_nxt(:,2) 

            ! check whether it generates excess demand or not
            share_temp = share_vec
            call calc_excess_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, & 
                                        seignorage_transfer, &
                                        bbeta_adj, P_div_P_h, P_div_P_h_nxt,   &  
                                        q_current, w_choice_vec, l_aggr_vec, sss, big_weight_vec, &  
                                        c_temp, c_spend_temp, share_temp, bF_share_temp, &
                                        k_temp_vec, b_temp_vec, excess_b, constraint_binding_vec, mc_temp)

        enddo

        if (iter >= 999) then
            write(*,*) 'No convergence in finding nom_i_high', sss
         !   stop
        end if
        nom_i_B = nom_i_guess
        excess_b_B = excess_b

      endif

      ! check that root is bracketed
      if ((excess_b_A*excess_b_B) > 0) then
         write(*,*) 'ERROR: Initial bracket does not contain root. A'
         stop
      endif

      ! swap a and b
      if ( abs(excess_b_A) < abs(excess_b_B) )then
         temp = nom_i_A
         nom_i_A = nom_i_B
         nom_i_B = temp

         temp = excess_b_A
         excess_b_A = excess_b_B
         excess_b_B = temp

      endif

      nom_i_C    = nom_i_A
      excess_b_C = excess_b_A
      excess_b_S = excess_b_B
      mflag = 1 ! set mflag

      iiter = 0
      do while ( excess_b_S /= 0 .and. abs(nom_i_A - nom_i_B) > 1E-15_dp)
        
        iiter = iiter + 1
        
        if ( (excess_b_A /= excess_b_C) .and. (excess_b_C /= excess_b_B) ) then ! inverse quadratic interpolation

          nom_i_S = nom_i_A * excess_b_B * excess_b_C / ( (excess_b_A - excess_b_B) * (excess_b_A - excess_b_C) ) + &
                    nom_i_B * excess_b_A * excess_b_C / ( (excess_b_B - excess_b_A) * (excess_b_B - excess_b_C) ) + &
                    nom_i_C * excess_b_A * excess_b_B / ( (excess_b_C - excess_b_A) * (excess_b_C - excess_b_B) )

        else ! secant method

          nom_i_S = nom_i_B - excess_b_B * (nom_i_B - nom_i_A) /   (excess_b_B - excess_b_A)

        endif

        if ( ( (  ( nom_i_S > ( 3*nom_i_A + nom_i_B )/4 .and. nom_i_S < nom_i_B) .or. &
                  ( nom_i_S < ( 3*nom_i_A + nom_i_B )/4 .and. nom_i_S > nom_i_B)  ) == .FALSE. ) .or. &
             (mflag == 1 .and. abs(nom_i_S - nom_i_B) >= abs(nom_i_B - nom_i_C)/2  )             .or. &
             (mflag == 0 .and. abs(nom_i_S - nom_i_B) >= abs(nom_i_B - nom_i_D)/2  )             .or. &
             (mflag == 1 .and. abs(nom_i_B - nom_i_C) <  abs(brent_delta)  )                     .or. &
             (mflag == 0 .and. abs(nom_i_B - nom_i_D) <  abs(brent_delta)  )                          &
           ) then

          nom_i_S = (nom_i_A + nom_i_B )/ 2
          mflag = 1
        else
          mflag = 0
        endif

          rf_vec_temp(:,1) = nom_i_S/homegood_infl_nxt(:,1) 
          rf_vec_temp(:,2) = (nom_i_S + nom_i_vec(2))/homegood_infl_nxt(:,2) 

            ! check whether it generates excess demand or not
            share_temp = share_vec
            call calc_excess_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, & 
                bbeta_adj, P_div_P_h, P_div_P_h_nxt,   &  
                            q_current, w_choice_vec, l_aggr_vec, sss, big_weight_vec, &  
                            c_temp, c_spend_temp, share_temp, bF_share_temp, k_temp_vec, b_temp_vec, excess_b_S, constraint_binding_vec, mc_temp)
                           
        nom_i_D = nom_i_C
        nom_i_C = nom_i_B
        excess_b_C = excess_b_B

        if ((excess_b_A*excess_b_S) < 0) then
          nom_i_B = nom_i_S
          excess_b_B = excess_b_S
        else
          nom_i_A = nom_i_S
          excess_b_A = excess_b_S
        endif

        ! swap a and b
        if ( abs(excess_b_A) < abs(excess_b_B) )then
           temp    = nom_i_A
           nom_i_A = nom_i_B
           nom_i_B = temp

           temp       = excess_b_A
           excess_b_A = excess_b_B
           excess_b_B = temp

        endif

            if (iiter > 2000) then 
                write(*,*) 'no convergence nom_i'
                stop
            endif
      enddo

      nom_i_guess = nom_i_B
    endif

    if (abs(excess_b_S)>1E-02_dp) then
        write(*,*) 'ERROR: EXCESS AGGREGATE BONDS HOLDINGS.'
        stop
    endif
    

    ! FOUND EQUILIBRIUM -------------------------------------------------
    ! store these choices
    nom_i_vec(1) = nom_i_guess
    endif


    rf_vec(:,1)  = nom_i_vec(1)/homegood_infl_nxt(:,1) 
    rf_vec(:,2)  = (nom_i_vec(1) + nom_i_vec(2))/homegood_infl_nxt(:,2)

! =============================================================================
! 5. GIVEN CURRENT OVERALL BOND HOLDINGS, UPDATE NOMINAL RATE SPREAD
! =============================================================================
    if (foreign_trading == 1 .and. outer_iter > 100) then 
    
    bF_share_temp = bF_share_vec 
    
    rf_vec_temp(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1) 
    rf_vec_temp(:,2) = (nom_i_vec(1) + nom_i_vec(2))/homegood_infl_nxt(:,2) 


    call calc_excess_foreign_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, bbeta_adj,  P_div_P_h, P_div_P_h_nxt, &  
                                       w_choice_vec, l_aggr_vec, sss, big_weight_vec, mc_temp, &  
                                       c_temp, c_spend_temp, share_vec, bF_share_temp, bh_temp_vec, bf_temp_vec, & 
                                       excess_b, safe_constraint_binding_vec)

                            
    nom_i_guess = nom_i_vec(2)
    
    ! store
    nom_i_A    = nom_i_guess
    excess_b_A = excess_b

      ! find bound between solver finds solution for nominal rate
      if (excess_b == 0) then

        excess_b_S = excess_b

      else
        if (excess_b > 0.0_dp) then

          iter = 0
          do while (excess_b > 0.0_dp .and. iter <= 5000)
              iter = iter + 1

              ! initial guess - last period nominal rate
              nom_i_guess = nom_i_guess - 1E-03_dp

              ! check whether it generates excess demand or not
              rf_vec_temp(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1)
              rf_vec_temp(:,2) = (nom_i_guess + nom_i_vec(1))/homegood_infl_nxt(:,2) 

              bF_share_temp = bF_share_vec
              call calc_excess_foreign_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, bbeta_adj,  P_div_P_h, P_div_P_h_nxt, &  
                                          w_choice_vec, l_aggr_vec, sss, big_weight_vec, mc_temp, &  
                                          c_temp, c_spend_temp, share_vec, bF_share_temp, bh_temp_vec, bf_temp_vec,&
                                          excess_b, safe_constraint_binding_vec)

          enddo
        
          if (iter >= 999) then
            write(*,*) 'No convergence in finding nom_i_low f'
            stop
          end if

          nom_i_B = nom_i_guess
          excess_b_B = excess_b

        else

          iter = 0
          do while (excess_b < 0.0_dp .and. iter <= 5000)
              
              iter = iter + 1

                ! initial guess - last period nominal rate
                nom_i_guess  = nom_i_guess + 1E-03_dp

                rf_vec_temp(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1) 
                rf_vec_temp(:,2) = (nom_i_guess + nom_i_vec(1))/homegood_infl_nxt(:,2) 

                ! check whether it generates excess demand or not
                bF_share_temp = bF_share_vec
                call calc_excess_foreign_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, bbeta_adj, P_div_P_h, P_div_P_h_nxt,&  
                                          w_choice_vec, l_aggr_vec, sss, big_weight_vec, mc_temp, &  
                                          c_temp, c_spend_temp, share_vec, bF_share_temp, bh_temp_vec, bf_temp_vec, & 
                                          excess_b, safe_constraint_binding_vec)
          enddo

          if (iter >= 999) then
              write(*,*) 'No convergence in finding nom_i_high foreign'
              stop
          end if
          nom_i_B = nom_i_guess
          excess_b_B = excess_b

        endif

        ! check that root is bracketed
        if ((excess_b_A*excess_b_B) > 0) then
           write(*,*) 'ERROR: Initial bracket does not contain root. B'
           stop
        endif

        ! swap a and b
        if ( abs(excess_b_A) < abs(excess_b_B) )then
           temp = nom_i_A
           nom_i_A = nom_i_B
           nom_i_B = temp

           temp = excess_b_A
           excess_b_A = excess_b_B
           excess_b_B = temp

        endif

        nom_i_C    = nom_i_A
        excess_b_C = excess_b_A
        excess_b_S = excess_b_B
        mflag = 1 ! set mflag

        iiter = 0 
        do while ( excess_b_S /= 0 .and. abs(nom_i_A - nom_i_B) > 1E-15_dp)
            
            iiter = iiter + 1

            if ( (excess_b_A /= excess_b_C) .and. (excess_b_C /= excess_b_B) ) then ! inverse quadratic interpolation

            nom_i_S = nom_i_A * excess_b_B * excess_b_C / ( (excess_b_A - excess_b_B) * (excess_b_A - excess_b_C) ) + &
                      nom_i_B * excess_b_A * excess_b_C / ( (excess_b_B - excess_b_A) * (excess_b_B - excess_b_C) ) + &
                      nom_i_C * excess_b_A * excess_b_B / ( (excess_b_C - excess_b_A) * (excess_b_C - excess_b_B) )

          else ! secant method

            nom_i_S = nom_i_B - excess_b_B * (nom_i_B - nom_i_A) /   (excess_b_B - excess_b_A)

          endif

          if ( ( (  ( nom_i_S > ( 3*nom_i_A + nom_i_B )/4 .and. nom_i_S < nom_i_B) .or. &
                    ( nom_i_S < ( 3*nom_i_A + nom_i_B )/4 .and. nom_i_S > nom_i_B)  ) == .FALSE. ) .or. &
               (mflag == 1 .and. abs(nom_i_S - nom_i_B) >= abs(nom_i_B - nom_i_C)/2  )             .or. &
               (mflag == 0 .and. abs(nom_i_S - nom_i_B) >= abs(nom_i_B - nom_i_D)/2  )             .or. &
               (mflag == 1 .and. abs(nom_i_B - nom_i_C) <  abs(brent_delta)  )                     .or. &
               (mflag == 0 .and. abs(nom_i_B - nom_i_D) <  abs(brent_delta)  )                          &
             ) then

            nom_i_S = (nom_i_A + nom_i_B )/ 2
            mflag = 1
          else
            mflag = 0
          endif

              rf_vec_temp(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1) 
              rf_vec_temp(:,2) = (nom_i_S + nom_i_vec(1))/homegood_infl_nxt(:,2) 

            ! check whether it generates excess demand or not
            bF_share_temp = bF_share_vec
            call calc_excess_foreign_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, & 
                                              bbeta_adj, P_div_P_h, P_div_P_h_nxt,   &  
                                              w_choice_vec, l_aggr_vec, sss, big_weight_vec, mc_temp,&  
                                              c_temp, c_spend_temp, share_vec, bF_share_temp, bh_temp_vec, bf_temp_vec,& 
                                              excess_b_S, safe_constraint_binding_vec)

          nom_i_D = nom_i_C
          nom_i_C = nom_i_B
          excess_b_C = excess_b_B

          if ((excess_b_A*excess_b_S) < 0) then
            nom_i_B = nom_i_S
            excess_b_B = excess_b_S
          else
            nom_i_A = nom_i_S
            excess_b_A = excess_b_S
          endif

          ! swap a and b
          if ( abs(excess_b_A) < abs(excess_b_B) )then
             temp    = nom_i_A
             nom_i_A = nom_i_B
             nom_i_B = temp

             temp       = excess_b_A
             excess_b_A = excess_b_B
             excess_b_B = temp

          endif

            if (iiter > 2000) then 
                write(*,*) 'no convergence nom_i f'
                stop
            endif


        enddo
        
        nom_i_guess = nom_i_B
      endif

      if (abs(excess_b_B)>1E-02_dp) then
          write(*,*) 'WARNING: EXCESS AGGREGATE BONDS HOLDINGS (FOREIGN).'
      endif

    nom_i_vec(2) = nom_i_B

endif

! =============================================================================
! 6. STORE VALUE FUNCTION, NEXT CAPITAL, NEXT WEALTH SHARES ======================
! =============================================================================

    c_spend_temp = c_spending_vec
    do iii = 1, n_I

            ! read in preferece parameters for current agent
            IES   = IES_vec(iii)
            bbeta = bbeta_vec(iii) + bbeta_adj(iii)
            gma   = gma_vec(iii)
            chi0   = chi0_vec(iii)
            v_normalization = v_normalization_vec(iii)
          
            c_cost     = P_div_P_h(iii)
            c_cost_nxt = P_div_P_h_nxt(:,iii)
          
            wealth     = aggr_wealth*wealth_share_grid(iii,sss) 

            labor      = l_aggr_vec(iii)

            c_spend =  minval( [maxval([c_spend_temp(iii), min_cons_sav]), & 
                       wealth + w_choice_vec(iii)*labor-min_cons_sav] )
          
            consumption = c_spend/c_cost
            c_vec(iii) = consumption
            
            r_alpha_vec = (share_vec(iii)-bF_share_vec(iii))*rf_vec(:,1)/(1.0_dp - omg) & 
                                           + bF_share_vec(iii)*rf_vec(:,2) & 
                          + (1.0_dp - share_vec(iii))*rk_vec

            labor_part =  1.0_dp + (1.0_dp/IES - 1.0_dp) * chi0*chi/(1.0_dp + chi) * labor**( (1.0_dp+chi)/chi )

            savings     = wealth + w_choice_vec(iii)*labor - consumption*c_cost
            
            ! calculate next period's wealth shares
            next_period_share = ( savings*r_alpha_vec + exp(dz_vec)*q_l_ss_vec(iii)  + seignorage_transfer(:,iii))/(tot_wealth_ss_vec(iii))

            v_vec_twisted = (v_temp(:,iii)*next_period_share)**(1.0-gma)
            EV = sum(big_weight_vec * v_vec_twisted) **(1.0 / (1.0 - gma))
            ! update wealth transitions accordingly
            savings_vec(iii) = savings
            nxt_wealth_vec(:,iii) = (savings*r_alpha_vec + seignorage_transfer(:,iii))
            
            objf = (v_normalization * labor_part**(1.0_dp/ies)*(consumption)**((ies - 1.0_dp)/ies) + bbeta *EV**((IES - 1.0_dp)/IES))**(ies / (ies - 1.0_dp))

            v_new(iii) = objf*( (wealth + q_l_ss_vec(iii))/tot_wealth_ss_vec(iii))**(-1)

    enddo

    theta_nxt = nxt_wealth_vec(:,1)/(sum(nxt_wealth_vec,2) )

    ! given new capital choices: next period capital
    k_next_new = sum(savings_vec*(1.0-share_vec)/q_current)

! =============================================================================
! 7. UPDATE TERMS OF TRADE ====================================================
! =============================================================================
      
      investment = next_k - k_aggr*(1.0-ddelta)
      
      ! check these derivations 
      s_new = s
      diff = 1.0_dp
      iter_inner = 0 
      do while (diff > conv_crit) 
          iter_inner = iter_inner + 1
          
          inv_h = investment*q_current/(1.0+((s_new)**(sigma-1.0))*(1.0-inv_share_h)/(inv_share_h))
          inv_f = inv_h * (s_new)**(sigma)*(1.0-inv_share_h)/inv_share_h 

          ! consumption in home units, cancels out below
          ch_vec(1) = c_spending_vec(1)/ ( (1.0-varsigma_vec(1) ) & 
              /(varsigma_vec(1) ) * s_new**(sigma-1.0_dp) + 1.0_dp) 
          ch_vec(2) = c_spending_vec(2)/ ( (1.0-varsigma_vec(2))/varsigma_vec(2) * s_new**(sigma-1.0_dp) + 1.0_dp) 
          s_update = ( sum((c_spending_vec - ch_vec))/(sum(ch_vec)) * (y_current(1)-inv_h)/( (y_current(2)-inv_f)) )**(-1)
          
          diff = abs(s_update-s_new)
          if (isnan(s_new)) then 
              write(*,*) 'NaN s'
              stop
          endif
          
          s_new = s_new + 0.1*(s_update-s_new)
         
          if (iter_inner > 10000) then 
             write(*,*) 'no convergence s', diff
             stop
         endif
      enddo
        
! =============================================================================
! 8. UPDATE CONSUMPTION CHOICE ===================================================
! =============================================================================

    c_spend_temp = c_spending_vec
    do iii = 1, n_I

            ! read in preferece parameters for current agent
            IES   = IES_vec(iii)
            bbeta = bbeta_vec(iii) + bbeta_adj(iii)
                gma   = gma_vec(iii)
            
            chi0   = chi0_vec(iii)
            v_normalization = v_normalization_vec(iii)
          
            c_cost     = P_div_P_h(iii)
            c_cost_nxt = P_div_P_h_nxt(:,iii)
          
            wealth     = aggr_wealth*wealth_share_grid(iii,sss) 

                labor = l_aggr_vec(iii)
                labor_nxt = l_aggr_nxt(:, iii)
           
            c_spend =  minval( [maxval([c_spend_temp(iii), min_cons_sav]), & 
                       wealth + w_choice_vec(iii)*labor-min_cons_sav] )
          
            consumption = c_spend/c_cost
           
            ! returns for resource accumulation
            r_alpha_vec = (share_vec(iii)-bF_share_vec(iii))*rf_vec(:,1)/(1.0_dp - omg) & 
                                           + bF_share_vec(iii)*rf_vec(:,2) & 
                          + (1.0_dp - share_vec(iii))*rk_vec

            labor_part =  1.0_dp + (1.0_dp/IES - 1.0_dp) * chi0*chi/(1.0_dp + chi) * labor**( (1.0_dp+chi)/chi )
            
            diff_inner = 1.0_dp
            iter_inner = 0

            do while (diff_inner>conv_crit)
             
                iter_inner  = iter_inner+1
                savings     = wealth + w_choice_vec(iii)*labor - consumption*c_cost
             
                ! calculate next period's wealth shares
                next_period_share = ( savings*r_alpha_vec + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer(:,iii))/tot_wealth_ss_vec(iii)
                
                ! returns for FOC (always include omg)
                r_alpha_omg_vec = (share_vec(iii)-bF_share_vec(iii))*rf_vec(:,1)/(1.0_dp - omg) & 
                                           + bF_share_vec(iii)*rf_vec(:,2) & 
                          + (1.0_dp - share_vec(iii))*rk_vec
                
                util_c_deriv = (labor_part / consumption)**(1/IES)
                mc_new(iii) = util_c_deriv*((q_l_ss_vec(iii) + wealth)/tot_wealth_ss_vec(iii))**(1/IES)
                v_vec_twisted = (v_temp(:,iii)*next_period_share)**(1.0_dp-gma)
             
                EV = sum(big_weight_vec*v_vec_twisted)**( 1.0_dp / (1.0_dp - gma))
                M_vec_mat(:, iii) = bbeta * mc_temp(:, iii)/util_c_deriv * ((v_temp(:, iii)/EV) **(1.0_dp/IES - gma))  * (next_period_share **(-gma))
       
                temp = sum(M_vec_mat(:, iii)*r_alpha_omg_vec*c_cost/c_cost_nxt*big_weight_vec)

                cons_update = consumption / (temp)**IES 

                cons_update = minval([maxval([cons_update, min_cons_sav]), (wealth + w_choice_vec(iii)*labor-min_cons_sav)/c_cost])
                
                diff_inner = abs(consumption-cons_update)
                 
                consumption = consumption + 0.5*(cons_update-consumption)
                
                if (phi_w == 0) then
                    labor_new = ( w_choice_vec(iii) / P_div_P_h(iii)    *   &
                                (1.0_dp + (1.0_dp/IES_vec(iii) - 1.0_dp) * chi0_vec(iii)*chi/(1.0_dp + chi) * &
                                labor**( (1.0_dp+chi)/chi ) ) /    &
                                ( 1.0_dp/IES_vec(iii) * consumption * chi0_vec(iii) )  )**(chi) 
                    labor = labor + 0.2*(labor_new-labor)
                endif

                if (iter_inner > 1000) then
                    if (diff_inner > 0.1_dp) then
                        write(*,*) 'NO CONVERGENCE 1.'
                        stop
                    else
                        write(*,*) 'WARNING: NO CONVERGENCE 1.'
                        diff_inner = 0.0_dp
                    endif
                endif
            enddo
        
        c_spending_vec_new(iii) = consumption*c_cost

        savings_vec(iii) = savings
        
        if (phi_w == 0) then 
            l_temp_vec(iii) = labor
        endif

    enddo
    
! =============================================================================
! 9. UPDATE LABOR ALLOCATION =====================================
! =============================================================================
    
    if (phi_w > 0.0_dp) then
    l_aggr_temp   = l_aggr_vec 
    w_choice_temp = w_choice_vec

    ! update wages in both countries
    do iii = 1,n_I

          w_temp = w_choice_temp(iii)
          w_diff = 1.0_dp
          iter = 0
          do while ( (w_diff > conv_crit .and. iter<100) .or. iter < 10 )

            iter = iter + 1

            w_current_next = w_temp 
            w_current_next(n_quad) = w_temp*exp(dz_vec(n_quad))
           
                    gma   = gma_vec(iii)
              w_temp_new = ( &
                            vareps_w/(1.0_dp-vareps_w)    * P_div_P_h(iii)   &
                            * ( -  1.0_dp/IES_vec(iii) * c_vec(iii) * chi0_vec(iii) * (l_aggr_temp(iii)**(1.0/chi)) / &
                               (1.0_dp + (1.0_dp/IES_vec(iii) - 1.0_dp) * chi0_vec(iii)*chi/(1.0_dp + chi) * &
                           l_aggr_temp(iii)**( (1.0_dp+chi)/chi ) ) )    &
                        + w_temp*phi_w/(1.0_dp-vareps_w)  * ( ( &
                                w_temp/w_current(iii)*infl_vec(iii)/P_div_P_h(iii) - 1.0 )   &
                               *w_temp/w_current(iii)*infl_vec(iii)/P_div_P_h(iii)  &
                            - sum( M_vec_mat(:,iii) * big_weight_vec * & 
                             homegood_infl_nxt(:,iii)/infl_nxt(:,iii)* &
                              ( w_next_choice(:,iii)/w_current_next*homegood_infl_nxt(:,iii) - 1.0 ) * &
                             (( w_next_choice(:,iii)/w_current_next)**2) * homegood_infl_nxt(:,iii)* l_aggr_nxt(:,iii)/l_aggr_temp(iii) ) &
                             ) & 
                             )

           w_diff = abs(w_temp_new-w_temp)
           w_temp = w_temp + 0.005_dp*(w_temp_new-w_temp)


          end do

          if (isnan(w_temp)) then
            write(*,*) 'error w isnan', M_vec_mat
            stop
          elseif (w_temp<=0) then
            write(*,*) 'error w<0'
            stop
          endif

          l_aggr_temp(iii) = kappa_vec(iii)*(( ( zf_vec(iii)**(1.0_dp-aalpha) )*(1.0_dp-aalpha)/w_temp/s_vec(iii))**(1.0_dp/aalpha))
          w_choice_temp(iii) =  w_temp

    enddo

    l_aggr_new = l_aggr_temp
    else
        
        l_aggr_new(1) = l_temp_vec(1)
        l_aggr_new(2) = l_temp_vec(2)

    endif

    w_choice_new = w_choice_vec

! =============================================================================
! 10. UPDATE INFLATION RATE =======================================================
! =============================================================================
         
        infl_new(1) = ((nom_i_vec_in(1)/((1.0+ih_last)**rho_i))**(1.0/(1.0-rho_i)) * exp(-tayl_ic_h)/(y_current(1)**phi_yh))**(1.0_dp/phi_h)
            
        infl_new(2) = (((nom_i_vec_in(1)+nom_i_vec_in(2))/((1.0+if_last)**rho_i))**(1.0/(1.0-rho_i)) *  exp(-tayl_ic_f)/(y_current(2)**phi_yf))**(1.0_dp/phi_f)
            
! =============================================================================
! UPDATE BOND PRICING IF NOT CALCULATED IN EQUILIBRIUM (INITIAL ITERATIONS ONLY)
! =============================================================================
    if (outer_iter <= 10) then   
        

        share_vec = 0.0_dp
        nom_i_vec(1) = maxval( [ 1*sum(big_weight_vec*M_vec_mat(:,1)/(1.0-omg)*(1.0_dp/homegood_infl_nxt(:,1))* P_div_P_h(1)/P_div_P_h_nxt(:,1)), & 
                                 1*sum(big_weight_vec*M_vec_mat(:,2)/(1.0-omg)*(1.0_dp/homegood_infl_nxt(:,1))* P_div_P_h(2)/P_div_P_h_nxt(:,2))  ] )**(-1) 

    else
        share_vec = share_temp
    endif

    if (foreign_trading == 1 .and. outer_iter > 100) then     
        bF_share_vec  = bF_share_temp 
    else
        bF_share_vec = 0.0
        nom_i_vec(2) = maxval( [ 1*sum(big_weight_vec*M_vec_mat(:,1)*(1.0_dp/homegood_infl_nxt(:,2))* P_div_P_h(1)/P_div_P_h_nxt(:,1)), & 
                                 1*sum(big_weight_vec*M_vec_mat(:,2)*(1.0_dp/homegood_infl_nxt(:,2))* P_div_P_h(2)/P_div_P_h_nxt(:,2))  ] )**(-1) - nom_i_vec(1)

    endif

end subroutine calc_equilibrium_and_update

subroutine calc_bond_prices( nxt_mat, n_nxt, nxt_mat_2, nxt_bond_prices,  sss, &
                               c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                               nom_i_vec, share_vec, bF_share_vec, inv_h, inv_f, bond_supply_shock, bond_supply_shock_nxt, &
                               who_prices_bond, bond_prices_out, E_rb, results_vec)

    
    use mod_param, only: n_I, n_quad,  idx_zf, idx_k, aalpha, smolyak_d, quad_weight_vec, omg_rho, &
                         dz_vec, ddelta, state_grid, varsigma_vec, sigma, idx_dis, tht_h_grid_mean, &
                         chix, ies_vec, gma_vec, omg_shift, b_lmbd, &
                         bbeta_vec, v_normalization_vec, chi0_vec, chi, wealth_share_grid, n_bond, &
                         n_interp, idx_wh, idx_wf, phi_h, phi_f, next_dis_mat, &
                         disast_p, disast_std, idx_omg, zeta, bg_yss, bbeta_coeff, idx_thh, omg_grid_mean

    real(dp), intent(in) :: nxt_mat_2(n_quad, 3), l_aggr_vec(2), nxt_bond_prices(n_quad,3), &
                            nxt_mat(n_quad, n_nxt), s, inv_h, inv_f, bond_supply_shock, bond_supply_shock_nxt, &
                            q_current, c_spending_vec(n_I), infl_vec(2) 
    real(dp), intent(in) :: bF_share_vec(n_I), share_vec(n_I), nom_i_vec(2)
    real(dp), intent(out)   :: bond_prices_out(3), results_vec(n_interp - 6*12), who_prices_bond(3), E_rb(6)

    integer, intent(in)  :: sss, n_nxt

    integer :: iii, ggg, ifail, nnn

    real(dp) :: s_vec(2), zf_vec(2), k_aggr, kappa_vec(2), & 
                q_nxt(n_quad), rk_vec(n_quad), &
                aggr_wealth, v_temp(n_quad, n_I), mc_temp(n_quad, n_I),  &
                w_next_choice(n_quad, 2), & 
                pi_nxt(n_quad), y_next(n_quad, 2), &
                kappa_nxt(n_quad,2), homegood_infl_nxt(n_quad, 2), &
                k_nxt(n_quad), s_nxt(n_quad), zf_nxt(n_quad,2), &
                rf_vec(n_quad, 2), pi_current(2), y_current(2), &
                excess_b, k_temp_vec(n_I), b_temp_vec(n_I), c_temp(n_I), & 
                share_temp(n_I), bF_share_temp(n_I), infl_nxt(n_quad,2), &
                c_vec(n_I), ch_vec(n_I), chi0, savings_vec(n_I)

    real(dp) :: next_k

    real(dp), parameter :: brent_delta = 5E-16_dp, min_cons_sav = 1E-08_dp
 
    real(dp) :: &
             nom_i_a, nom_i_b, excess_b_A, excess_b_B, w_choice_vec(2), &
             nom_i_C, temp, excess_b_C, excess_b_S, nom_i_S, nom_i_D, wealth, consumption, r_alpha_vec(n_quad), &
             diff_inner, iter_inner, savings, objf, objgrd(2), cons_update, IES, bbeta, &
             v_normalization, next_period_share(n_quad), v_vec(n_quad), EV, &
             M_vec_mat(n_quad, n_I), weight_vec(n_I), w_diff, temp_vec(n_quad), &
             share_low, share_high, real_bond, c_spend_temp(n_I), &
             mid_point, k_param, derivative_consumption, M_weighted_vec(n_quad), nom_i_guess, &
             share_guess, gov_debt,  gma, seignorage_transfer(n_quad, n_I), bbeta_adj(n_I), &
             ch, cf, c_spend, labor, c_cost, c_cost_nxt(n_quad), & 
             v_vec_twisted(n_quad), iota(n_quad), seigniorage_nxt(n_quad),& 
             l_aggr_nxt(n_quad,2), &
             s_update, P_div_P_h(2), P_div_P_h_nxt(n_quad,2), w_current_next(n_quad), &
             w_current(2), w_choice_new_vec(2), kappa_guess, kappa_share, util, util_c_deriv, labor_part, &
             l_aggr_vec_temp(2), l_aggr_temp(2), l_aggr_temp_new(2), w_choice_temp(2), w_choice_temp_new(2), &
             kappa_vec_new(2), p_dis, big_weight_vec(n_quad), deriv_helper(n_quad),  omg, & 
             E_rfh, E_rff, E_rk, investment, q_temp, test_ratio,  v_store(n_I), & 
             ce_store(n_I), omg_nxt(n_quad), ch_spending_nxt(n_quad), cf_spending_nxt(n_quad), & 
             consumption_imports, consumption_exports, nfa_tmp, E_val, nx_tmp, E_val_rk, E_val_rfh, E_val_rff, E_rf_spread,&
             cov_rf_spread_m(n_I), var_rf_spread, var_m(n_I), corr_rf_spread_m(12)
 
 
 integer :: mflag, ix1, ix2

! =============================================================================
! CURRENT STATE ===============================================================
! =============================================================================
    

    
    ! current capital stock
    k_aggr = state_grid(idx_k, sss)

    ! store relative productivity in vector
    zf_vec = [ 1.0_dp, zeta*exp(state_grid(idx_zf, sss)) ]

    ! current disaster probability
    p_dis       = exp(state_grid(idx_dis, sss))


    ! current weight vec
    big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
    big_weight_vec(n_quad)     = p_dis

    ! w_current is stored as w*P_h/P, so in bundle units
    ! allows to calculate w'/w*local_infl = w'/w_current*infl*P'_h/P'
    w_current = [state_grid(idx_wh, sss), state_grid(idx_wf, sss)]

! =============================================================================
! CURRENT PERIOD GUESSES ======================================================
! =============================================================================

    ! store relative prices in vector
    s_vec   = [1.0_dp, s]

    ! given current guess of aggregate labor and wage, allocate capital across countries
    ! allows to calculate initial wealth
    kappa_vec(1) = k_aggr * ( (s_vec(2) )**(1.0/(1.0-aalpha))) * &
                      zf_vec(1) / zf_vec(2) * l_aggr_vec(1)/l_aggr_vec(2) / & 
                   (1.0 +   ( (s_vec(2) )**(1.0/(1.0-aalpha))) * &
                      zf_vec(1) / zf_vec(2) * l_aggr_vec(1)/l_aggr_vec(2) )

    kappa_vec(2) = k_aggr - kappa_vec(1)

    ! given this calculate profits and aggregate wealth
    pi_current   = (s_vec**(-1))*aalpha * (zf_vec*l_aggr_vec/kappa_vec)**(1.0-aalpha)
    
    aggr_wealth  = kappa_vec(1)*( (1-ddelta)*q_current + pi_current(1) ) +  &
                   kappa_vec(2)*( (1-ddelta)*q_current + pi_current(2) )
    
    y_current    = (zf_vec**(1.0-aalpha)) * kappa_vec**aalpha * l_aggr_vec**(1.0-aalpha)
    
    w_choice_vec = (s_vec**(-1))*(1.0-aalpha)*y_current/l_aggr_vec

    ! cost of unit of bundle in terms of ch
    P_div_P_h(1) = (varsigma_vec(1)  & 
        + (1.0 -varsigma_vec(1) )*s**(sigma-1))**(1/(1-sigma))
    P_div_P_h(2) = s**(-1) * (varsigma_vec(2)*s**(1.0-sigma) + (1.0 -varsigma_vec(2)))**(1/(1-sigma))
    
    omg      = exp(state_grid(idx_omg,sss)) + omg_shift  - b_lmbd*bond_supply_shock

    
    bbeta_adj(1) = - bbeta_coeff*(state_grid(idx_thh,sss) - tht_h_grid_mean)
    bbeta_adj(2) = + bbeta_coeff*(state_grid(idx_thh,sss) - tht_h_grid_mean)
! =============================================================================
! NEXT PERIOD =================================================================
! =============================================================================


    ! currently assumed capital in the next period
    ! vector because depreciation shock in disaster state 
    k_nxt      = nxt_mat_2(:,1)

    ! relative productivity next period
    zf_nxt(:,1)   = 1.0_dp 
    zf_nxt(:,2)   = zeta*exp(nxt_mat_2(:,2))
    omg_nxt      = exp(nxt_mat_2(:,3)) + omg_shift  - b_lmbd*bond_supply_shock_nxt

    ! value functions
    do iii = 1,n_I
        v_temp(:,iii) = nxt_mat(:, iii )
        mc_temp(:, iii) = nxt_mat(:, n_I + iii)
    enddo

    ! tomorrow's price, labor choice, infl and price of time endowment
    s_nxt               = nxt_mat(:,2*n_I + 1)
    q_nxt               = nxt_mat(:,2*n_I + 2) 
    l_aggr_nxt(:,1)     = nxt_mat(:,2*n_I + 3)
    l_aggr_nxt(:,2)     = nxt_mat(:,2*n_I + 4)
    infl_nxt(:,1)       = nxt_mat(:,2*n_I + 5)
    infl_nxt(:,2)       = nxt_mat(:,2*n_I + 6)
    ch_spending_nxt     = nxt_mat(:,2*n_I + 7)*exp(dz_vec)
    cf_spending_nxt     = nxt_mat(:,2*n_I + 8)*exp(dz_vec)

    homegood_infl_nxt(:,1) = infl_nxt(:,1)*( ( varsigma_vec(1)  + (1.0-varsigma_vec(1) )*    s**(sigma-1.0) ) / & 
                                          ( varsigma_vec(1) + (1.0-varsigma_vec(1) )*s_nxt**(sigma-1.0) ) )**(1.0/(1.0-sigma)) 
    homegood_infl_nxt(:,2) = infl_nxt(:,2)* s_nxt/s*( ( varsigma_vec(2)*    s**(1.0-sigma)+ (1.0-varsigma_vec(2)) ) / & 
                                          ( varsigma_vec(2)*s_nxt**(1.0-sigma)+ (1.0-varsigma_vec(2)) ) )**(1.0/(1.0-sigma)) 

    ! given employed labor, allocate capital to equalize capital returns across countries
    kappa_nxt(:,1) = k_nxt *  ((s_nxt)**(1.0/(1.0-aalpha))) * zf_nxt(:,1) / zf_nxt(:,2) * l_aggr_nxt(:,1)/l_aggr_nxt(:,2) / & 
                       (1.0 + ((s_nxt)**(1.0/(1.0-aalpha))) * zf_nxt(:,1) / zf_nxt(:,2) * l_aggr_nxt(:,1)/l_aggr_nxt(:,2) )
    
    kappa_nxt(:,2) = k_nxt - kappa_nxt(:,1)

    next_k       = k_nxt(1)

    investment = next_k - k_aggr*(1.0-ddelta)

    ! output in next period
    y_next(:,1) = exp((log(zf_nxt(:,1)) + dz_vec)*(1.0_dp-aalpha))*l_aggr_nxt(:,1)* & 
                  ( (kappa_nxt(:,1)/l_aggr_nxt(:,1))**aalpha)
    y_next(:,2) = exp((log(zf_nxt(:,2)) + dz_vec)*(1.0_dp-aalpha))* l_aggr_nxt(:,2)* & 
                  ( (kappa_nxt(:,2)/l_aggr_nxt(:,2))**aalpha)
    
    w_next_choice(:,1) = (1-aalpha)*y_next(:,1)/l_aggr_nxt(:,1)
    w_next_choice(:,2) = s_nxt**(-1) * (1-aalpha)*y_next(:,2)/l_aggr_nxt(:,2)

    ! profit next period per unit of capital
    pi_nxt = aalpha*exp( (log(zf_nxt(:,1)) + dz_vec)*(1.0_dp-aalpha)) &
           *(l_aggr_nxt(:,1)**(1.0_dp-aalpha))*kappa_nxt(:,1)**(aalpha-1.0)
   
    ! capital return, adjustment in disaster state
    rk_vec = ( (1.0_dp-ddelta)*q_nxt + pi_nxt)/q_current
    rk_vec(n_quad) =exp(dz_vec(n_quad))*rk_vec(n_quad) 

    P_div_P_h_nxt(:,1) = (varsigma_vec(1)  + (1.0 -varsigma_vec(1)  )*s_nxt**(sigma-1))**(1.0/(1.0-sigma)) 
    P_div_P_h_nxt(:,2) = s_nxt**(-1) * (varsigma_vec(2)*s_nxt**(1.0-sigma) + (1.0 -varsigma_vec(2)))**(1.0/(1.0-sigma))

    ! calculate risk-free returns, given nominal rate guess 
    ! transform everything in units of home good
    rf_vec(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1)
    rf_vec(:,2) = (nom_i_vec(1) + nom_i_vec(2))/homegood_infl_nxt(:,2)
    
    ! seignorage
    if (bg_yss > 0) then 
        seigniorage_nxt = cf_spending_nxt/(ch_spending_nxt + cf_spending_nxt)*(bg_yss + bond_supply_shock_nxt)*(ch_spending_nxt + cf_spending_nxt)*omg_nxt
        seignorage_transfer(:,1) =  seigniorage_nxt 
        seignorage_transfer(:,2) = -seigniorage_nxt
    else
    seignorage_transfer = 0.0_dp
    endif

! =============================================================================
! CALCULATE MARG. UTIL AND PRICING KERNEL =====================================
! =============================================================================

    M_vec_mat = 0.0_dp
    do iii = 1, n_I

          IES             = IES_vec(iii)
          bbeta           = bbeta_vec(iii) + bbeta_adj(iii)
          gma             = gma_vec(iii)
          v_normalization = v_normalization_vec(iii)
          chi0            = chi0_vec(iii) 

          c_cost     = P_div_P_h(iii)
          c_cost_nxt = P_div_P_h_nxt(:,iii)

          wealth = aggr_wealth*wealth_share_grid(iii,sss) 

          labor      = l_aggr_vec(iii)
          
          c_spend =  minval([maxval([c_spending_vec(iii), min_cons_sav]), & 
                     wealth + w_choice_vec(iii)*labor-min_cons_sav])
          
          consumption = c_spend/c_cost
          c_vec(iii)  = consumption

          r_alpha_vec = (share_vec(iii)-bF_share_vec(iii))*rf_vec(:,1)/(1.0_dp - omg) & 
                                           + bF_share_vec(iii)*rf_vec(:,2) & 
                          + (1.0_dp - share_vec(iii))*rk_vec

          savings     = wealth + w_choice_vec(iii)*labor - c_spend 
          if (savings < 0.0) then 
              write(*,*) 'Savings:', savings
              stop
          endif
          savings_vec(iii) = savings

          ! calculate next period's wealth shares
          next_period_share = (savings*r_alpha_vec + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer(:,iii) )/(tot_wealth_ss_vec(iii))

          v_vec_twisted = (v_temp(:,iii)*next_period_share)**(1.0-gma)

          EV = sum(big_weight_vec*v_vec_twisted)

          call util_fun(consumption, labor, iii, util, util_c_deriv, labor_part)
          if (IES > 1.0_dp-sqrt_eps .and. IES < 1.0_dp + sqrt_eps) then 
            objf = exp( v_normalization*util +  bbeta/(1.0_dp-gma)*log(EV) )
          else
            objf = ( v_normalization*util +  bbeta*  (EV**( (1.0_dp - (1.0_dp/IES))/(1.0_dp-gma) )) )**(1.0_dp/(1.0_dp-1.0_dp/IES))
          endif  

          v_store(iii)  = objf
          ce_store(iii) = EV**(1.0_dp/(1.0_dp - gma))

          M_vec_mat(:,iii) =  bbeta * mc_temp(:, iii) / util_c_deriv * (v_temp(:, iii) / ce_store(iii))**(1/ies - gma) * next_period_share**(-gma)
    enddo

    if (any(isnan(M_vec_mat))) then
        write(*,*) 'M_vec_mat is NaN.'
        stop
    endif


     ! Home bond pricing
     bond_prices_out(1) = maxval([sum(big_weight_vec*M_vec_mat(:,1)*(nxt_bond_prices(:,1)/homegood_infl_nxt(:,1))* P_div_P_h(1)/P_div_P_h_nxt(:,1)), &
                                  sum(big_weight_vec*M_vec_mat(:,2)*(nxt_bond_prices(:,1)/homegood_infl_nxt(:,1))* P_div_P_h(2)/P_div_P_h_nxt(:,2)) &
                                      ])

     ! Foreign bond pricing 
     bond_prices_out(2) = maxval([ sum(big_weight_vec*M_vec_mat(:,1)*nxt_bond_prices(:,2)/homegood_infl_nxt(:,2)* P_div_P_h(1)/P_div_P_h_nxt(:,1)), & 
                                   sum(big_weight_vec*M_vec_mat(:,2)*nxt_bond_prices(:,2)/homegood_infl_nxt(:,2)* P_div_P_h(2)/P_div_P_h_nxt(:,2))])

     bond_prices_out(3) = maxval([sum(big_weight_vec*M_vec_mat(:,1)*(nxt_bond_prices(:,3)/homegood_infl_nxt(:,1))/(1.0-omg)* P_div_P_h(1)/P_div_P_h_nxt(:,1)), &
                                  sum(big_weight_vec*M_vec_mat(:,2)*(nxt_bond_prices(:,3)/homegood_infl_nxt(:,1))/(1.0-omg)* P_div_P_h(2)/P_div_P_h_nxt(:,2)) &
                                      ])
    
    who_prices_bond(1) = real( maxloc([ sum(big_weight_vec*M_vec_mat(:,1)*(nxt_bond_prices(:,1)/homegood_infl_nxt(:,1))* P_div_P_h(1)/P_div_P_h_nxt(:,1) ), & 
                                        sum(big_weight_vec*M_vec_mat(:,2)*(nxt_bond_prices(:,1)/homegood_infl_nxt(:,1))* P_div_P_h(2)/P_div_P_h_nxt(:,2)) ], 1, .true.), dp)

    who_prices_bond(2) = real( maxloc([ sum(big_weight_vec*M_vec_mat(:,1)*nxt_bond_prices(:,2)/homegood_infl_nxt(:,2)* P_div_P_h(1)/P_div_P_h_nxt(:,1)), & 
                                        sum(big_weight_vec*M_vec_mat(:,2)*nxt_bond_prices(:,2)/homegood_infl_nxt(:,2)* P_div_P_h(2)/P_div_P_h_nxt(:,2))], 1, .true.), dp)

    who_prices_bond(3) = real( maxloc([ sum(big_weight_vec*M_vec_mat(:,1)*(nxt_bond_prices(:,3)/homegood_infl_nxt(:,1))/(1.0-omg)* P_div_P_h(1)/P_div_P_h_nxt(:,1) ), & 
                                        sum(big_weight_vec*M_vec_mat(:,2)*(nxt_bond_prices(:,3)/homegood_infl_nxt(:,1))/(1.0-omg)* P_div_P_h(2)/P_div_P_h_nxt(:,2)) ], 1, .true.), dp)


! =============================================================================
! CALCULATE EXPECTED RETURNS ===================================================
! =============================================================================
E_rfh     = sum(big_weight_vec*((rf_vec(:,1))*P_div_P_h(1)/P_div_P_h_nxt(:,1))) 

E_rff     = sum(big_weight_vec*(rf_vec(:,2)*P_div_P_h(1)/P_div_P_h_nxt(:,1))) 
E_rk      = sum(big_weight_vec*(rk_vec*P_div_P_h(1)/P_div_P_h_nxt(:,1))) 

E_rb(1)   = sum(big_weight_vec*(nxt_bond_prices(:,1)/infl_nxt(:,1))/bond_prices_out(1))
E_rb(2)   = sum(big_weight_vec*(nxt_bond_prices(:,2)/infl_nxt(:,2))/bond_prices_out(2)*P_div_P_h_nxt(:,2)/P_div_P_h(2)*P_div_P_h(1)/P_div_P_h_nxt(:,1))
E_rb(3)   = sum(big_weight_vec*(nxt_bond_prices(:,3)/infl_nxt(:,1))/bond_prices_out(3))
E_rb(4)   = sum(big_weight_vec*(nxt_bond_prices(:,1)/infl_nxt(:,1))/bond_prices_out(1)*P_div_P_h_nxt(:,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(:,2))
E_rb(5)   = sum(big_weight_vec*(nxt_bond_prices(:,2)/infl_nxt(:,2))/bond_prices_out(2))
E_rb(6)   = sum(big_weight_vec*(nxt_bond_prices(:,3)/infl_nxt(:,1))/bond_prices_out(3)*P_div_P_h_nxt(:,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(:,2))



E_rf_spread = sum(big_weight_vec*(log(rf_vec(:,2)/rf_vec(:,1))))
var_rf_spread = sum(big_weight_vec*(log(rf_vec(:,2)/rf_vec(:,1))**2)) - E_rf_spread**2
var_m(1) = sum(big_weight_vec*(log(M_vec_mat(:,1)))**2) - sum(big_weight_vec*(log(M_vec_mat(:,1))))**2
var_m(2) = sum(big_weight_vec*(log(M_vec_mat(:,2)*P_div_P_h_nxt(:,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(:,2)))**2) - sum(big_weight_vec*(log(M_vec_mat(:,2)*P_div_P_h_nxt(:,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(:,2))))**2

cov_rf_spread_m(1) = sum(big_weight_vec*((log(rf_vec(:, 2)/rf_vec(:,1)))*log(M_vec_mat(:,1)))) - sum(big_weight_vec*(log(rf_vec(:,2)/rf_vec(:,1))))*sum(big_weight_vec*(log(M_vec_mat(:,1))))
cov_rf_spread_m(2) = sum(big_weight_vec*((log(rf_vec(:, 2)/rf_vec(:,1)))*log(M_vec_mat(:,2)*P_div_P_h_nxt(:,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(:,2)))) - sum(big_weight_vec*(log(rf_vec(:,2)/rf_vec(:,1))))*sum(big_weight_vec*(log(M_vec_mat(:,2)*P_div_P_h_nxt(:,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(:,2))))

corr_rf_spread_m(1) = cov_rf_spread_m(1)/sqrt(var_rf_spread*var_m(1))
corr_rf_spread_m(2) = cov_rf_spread_m(2)/sqrt(var_rf_spread*var_m(2))


E_rf_spread = sum(quad_weight_vec*(log(rf_vec(1:n_quad-1,2)/rf_vec(1:n_quad-1,1))))
var_rf_spread = sum(quad_weight_vec*(log(rf_vec(1:n_quad-1,2)/rf_vec(1:n_quad-1,1))**2)) - E_rf_spread**2

var_m(1) = sum(quad_weight_vec*(log(M_vec_mat(1:n_quad-1,1)))**2) - sum(big_weight_vec*(log(M_vec_mat(1:n_quad-1,1))))**2
var_m(2) = sum(quad_weight_vec*(log(M_vec_mat(1:n_quad-1,2)*P_div_P_h_nxt(1:n_quad-1,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(1:n_quad-1,2)))**2)  & 
         - sum(quad_weight_vec*(log(M_vec_mat(1:n_quad-1,2)*P_div_P_h_nxt(1:n_quad-1,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(1:n_quad-1,2))))**2

cov_rf_spread_m(1) = sum(quad_weight_vec*((log(rf_vec(1:n_quad-1, 2)/rf_vec(1:n_quad-1,1)))*log(M_vec_mat(1:n_quad-1,1)))) &
                   - sum(quad_weight_vec*( log(rf_vec(1:n_quad-1,2)/ rf_vec(1:n_quad-1,1))))*sum(quad_weight_vec*(log(M_vec_mat(1:n_quad-1,1))))
cov_rf_spread_m(2) = sum(quad_weight_vec*((log(rf_vec(1:n_quad-1, 2)/rf_vec(1:n_quad-1,1)))*log(M_vec_mat(1:n_quad-1,2)*P_div_P_h_nxt(1:n_quad-1,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(1:n_quad-1,2)))) & 
                   - sum(quad_weight_vec*( log(rf_vec(1:n_quad-1,2)/ rf_vec(1:n_quad-1,1))))*sum(quad_weight_vec*(log(M_vec_mat(1:n_quad-1,2)*P_div_P_h_nxt(1:n_quad-1,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(1:n_quad-1,2))))

corr_rf_spread_m(5) = cov_rf_spread_m(1)/sqrt(var_rf_spread*var_m(1))
corr_rf_spread_m(6) = cov_rf_spread_m(2)/sqrt(var_rf_spread*var_m(2))

! correlation with disaster probability changes
E_rf_spread = sum(quad_weight_vec*next_dis_mat(:,sss))
var_rf_spread = sum(quad_weight_vec*((next_dis_mat(:,sss))**2)) - E_rf_spread**2

cov_rf_spread_m(1) = sum(quad_weight_vec*((next_dis_mat(:,sss))*log(M_vec_mat(1:n_quad-1,1)))) &
    - sum(quad_weight_vec*(next_dis_mat(:,sss)))*sum(quad_weight_vec*(log(M_vec_mat(1:n_quad-1,1))))
cov_rf_spread_m(2) = sum(quad_weight_vec*((next_dis_mat(:,sss))*log(M_vec_mat(1:n_quad-1,2)*P_div_P_h_nxt(1:n_quad-1,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(1:n_quad-1,2)))) & 
    - sum(quad_weight_vec*(next_dis_mat(:,sss)))*sum(quad_weight_vec*(log(M_vec_mat(1:n_quad-1,2)*P_div_P_h_nxt(1:n_quad-1,1)/P_div_P_h(1)*P_div_P_h(2)/P_div_P_h_nxt(1:n_quad-1,2))))

corr_rf_spread_m(9) = cov_rf_spread_m(1)/sqrt(var_rf_spread*var_m(1))
corr_rf_spread_m(10) = cov_rf_spread_m(2)/sqrt(var_rf_spread*var_m(2))



E_rf_spread = sum(big_weight_vec*(log(rf_vec(:,2)/rf_vec(:,1))))

results_vec(1) = y_current(1)
results_vec(2) = y_current(2)
results_vec(3) = l_aggr_vec(1)
results_vec(4) = l_aggr_vec(2)
results_vec(5) = c_vec(1) 
results_vec(6) = c_vec(2) 
results_vec(7) = k_nxt(1) - (1.0-ddelta)*k_aggr
results_vec(8) = pi_current(1)/P_div_P_h(1) 
results_vec(9) = pi_current(2)/P_div_P_h(1) 
results_vec(10)= P_div_P_h(1) 
results_vec(11)= P_div_P_h(2) 
results_vec(12)= P_div_P_h(1) / P_div_P_h(2) 
results_vec(13)= s
results_vec(14)= E_rfh
results_vec(15)= E_rff
results_vec(16)= E_rk
results_vec(17)= q_current/P_div_P_h(1)
do iii = 1,n_I
results_vec(17+iii)= share_vec(iii)
results_vec(19+iii)= share_vec(iii) - bF_share_vec(iii)
enddo
results_vec(22)= aggr_wealth/P_div_P_h(1)
results_vec(23)= w_choice_vec(1)/P_div_P_h(1)
results_vec(24)= w_choice_vec(2)/P_div_P_h(2)
results_vec(25)= infl_vec(1)
results_vec(26)= infl_vec(2)
results_vec(27)= nom_i_vec(1) 
results_vec(28)= nom_i_vec(1) + nom_i_vec(2) 
results_vec(29)= (savings_vec(1)*(1.0-share_vec(1)))/P_div_P_h(1)
results_vec(30)= kappa_vec(1)
results_vec(31)= kappa_vec(2)
results_vec(32)= sum(big_weight_vec*(rf_vec(:,1)*P_div_P_h(2)/P_div_P_h_nxt(:,2)))
results_vec(33)= sum(big_weight_vec*(rf_vec(:,2)*P_div_P_h(2)/P_div_P_h_nxt(:,2)))
results_vec(34)= sum(big_weight_vec*(rk_vec*P_div_P_h(2)/P_div_P_h_nxt(:,2)))
results_vec(35)= sum(big_weight_vec*kappa_nxt(:,1)/k_nxt*k_nxt(1))
results_vec(36)= v_store(1) 
results_vec(37)= v_store(2) 
results_vec(38)= inv_h/P_div_P_h(1)
results_vec(39)= s**(-1) * inv_f/P_div_P_h(2)
results_vec(40)= savings_vec(1)*(share_vec(1) - bF_share_vec(1))/P_div_P_h(1)
results_vec(41)= bond_supply_shock*sum(c_spending_vec)/P_div_P_h(1)
results_vec(42)= ce_store(1) 
results_vec(43)= ce_store(2) 
results_vec(44)= savings_vec(1)/P_div_P_h(1) 
results_vec(45)= savings_vec(2)/P_div_P_h(1) 
results_vec(46) = E_rf_spread
results_vec(47) = corr_rf_spread_m(1)
results_vec(48) = corr_rf_spread_m(2)
results_vec(49) = corr_rf_spread_m(5)
results_vec(50) = corr_rf_spread_m(6)
results_vec(51) = corr_rf_spread_m(9)
results_vec(52) = corr_rf_spread_m(10)
end subroutine calc_bond_prices

subroutine calc_valuation( nxt_mat, n_nxt, nxt_mat_2, nxt_valuation,  sss, &
                               c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                               nom_i_vec, share_vec, bF_share_vec, inv_h, inv_f, bond_supply_shock, bond_supply_shock_nxt, &
                               valuation_vec)

    
    use mod_param, only: n_I, n_quad,  idx_zf, idx_k, aalpha, smolyak_d, quad_weight_vec, &
                         dz_vec, ddelta, state_grid, varsigma_vec, sigma, idx_dis, bg_yss, &
                         chix, ies_vec, gma_vec, omg_shift, b_lmbd, &
                         bbeta_vec, v_normalization_vec, chi0_vec, chi, wealth_share_grid, n_bond, &
                         n_interp, &
                         phi_w, idx_wh, idx_wf, phi_h, phi_f, &
                         disast_p, disast_std, idx_omg, & 
                         zeta, bg_yss, omg_rho, omg_grid_mean

    real(dp), intent(in) :: nxt_mat_2(n_quad, 3), l_aggr_vec(2), nxt_valuation(n_quad,7), &
                            nxt_mat(n_quad, n_nxt), s, inv_h, inv_f, bond_supply_shock, bond_supply_shock_nxt, &  
                            q_current, c_spending_vec(n_I), infl_vec(2) 
    real(dp), intent(in) :: bF_share_vec(n_I), share_vec(n_I), nom_i_vec(2)
    real(dp), intent(out)   :: valuation_vec(7)

    integer, intent(in)  :: sss, n_nxt

    integer :: iii, ggg, ifail, nnn

    real(dp) :: s_vec(2), zf_vec(2), k_aggr, kappa_vec(2), & 
                q_nxt(n_quad), rk_vec(n_quad), &
                aggr_wealth, v_temp(n_quad, n_I),  &
                w_next_choice(n_quad, 2), & 
                pi_nxt(n_quad), y_next(n_quad, 2), &
                kappa_nxt(n_quad,2), homegood_infl_nxt(n_quad, 2), &
                k_nxt(n_quad), s_nxt(n_quad), zf_nxt(n_quad,2), &
                rf_vec(n_quad, 2), pi_current(2), y_current(2), &
                excess_b, k_temp_vec(n_I), b_temp_vec(n_I), c_temp(n_I), & 
                share_temp(n_I), bF_share_temp(n_I), infl_nxt(n_quad,2), &
                c_vec(n_I), ch_vec(n_I), chi0, savings_vec(n_I)

    real(dp) :: next_k

    real(dp), parameter :: brent_delta = 5E-16_dp, min_cons_sav = 1E-08_dp
 
    real(dp) :: &
             nom_i_a, nom_i_b, excess_b_A, excess_b_B, w_choice_vec(2), &
             nom_i_C, temp, excess_b_C, excess_b_S, nom_i_S, nom_i_D, wealth, consumption, r_alpha_vec(n_quad), &
             diff_inner, iter_inner, savings, objf, objgrd(2), cons_update, IES, bbeta, &
             v_normalization, v_vec(n_quad), EV, &
             M_vec_mat(n_quad, n_I), weight_vec(n_I), w_diff, temp_vec(n_quad), &
             share_low, share_high, real_bond, c_spend_temp(n_I), &
             mid_point, k_param, derivative_consumption, M_weighted_vec(n_quad), &
             share_guess, gov_debt,  gma, seignorage_transfer(n_quad,n_I), seigniorage_nxt(n_quad),  &
             ch, cf, c_spend, labor, c_cost, c_cost_nxt(n_quad), ch_spending_nxt(n_quad), cf_spending_nxt(n_quad),  & 
             v_vec_twisted(n_quad), iota(n_quad), omg_nxt(n_quad), & 
             l_aggr_nxt(n_quad,2), &
             s_update, P_div_P_h(2), P_div_P_h_nxt(n_quad,2), w_current_next(n_quad), &
             w_current(2), w_choice_new_vec(2), kappa_guess, kappa_share, util, util_c_deriv, labor_part, &
             l_aggr_vec_temp(2), l_aggr_temp(2), l_aggr_temp_new(2), w_choice_temp(2), w_choice_temp_new(2), &
             kappa_vec_new(2), p_dis, big_weight_vec(n_quad), deriv_helper(n_quad),  omg, & 
             E_rfh, E_rff, E_rk, investment, q_temp, test_ratio,  v_store(n_I), & 
             ce_store(n_I), & 
             consumption_imports, consumption_exports, E_kpa, nfa_tmp, E_val, nx_tmp, E_val_seign, E_val_rk, E_val_rfh, E_val_rff, & 
         E_nfa_end, val_tmp(n_quad), E_nx
 
 
 integer :: mflag, ix1, ix2

! =============================================================================
! CURRENT STATE ===============================================================
! =============================================================================

    
    ! current capital stock
    k_aggr = state_grid(idx_k, sss)

    ! store relative productivity in vector
    zf_vec = [ 1.0_dp , zeta*exp(state_grid(idx_zf, sss)) ]

    ! current disaster probability
    p_dis       = exp(state_grid(idx_dis, sss))


    ! current weight vec
    big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
    big_weight_vec(n_quad)     = p_dis

    ! w_current is stored as w*P_h/P, so in bundle units
    ! allows to calculate w'/w*local_infl = w'/w_current*infl*P'_h/P'
    w_current = [state_grid(idx_wh, sss), state_grid(idx_wf, sss)]

! =============================================================================
! CURRENT PERIOD GUESSES ======================================================
! =============================================================================

    ! store relative prices in vector
    s_vec   = [1.0_dp, s]

    ! given current guess of aggregate labor and wage, allocate capital across countries
    ! allows to calculate initial wealth
    kappa_vec(1) = k_aggr * ( (s_vec(2) )**(1.0/(1.0-aalpha))) * &
                      zf_vec(1) / zf_vec(2) * l_aggr_vec(1)/l_aggr_vec(2) / & 
                   (1.0 +   ( (s_vec(2) )**(1.0/(1.0-aalpha))) * &
                      zf_vec(1) / zf_vec(2) * l_aggr_vec(1)/l_aggr_vec(2) )

    kappa_vec(2) = k_aggr - kappa_vec(1)

    ! given this calculate profits and aggregate wealth
    pi_current   = (s_vec**(-1))*aalpha * (zf_vec*l_aggr_vec/kappa_vec)**(1.0-aalpha)
    
    aggr_wealth  = kappa_vec(1)*( (1-ddelta)*q_current + pi_current(1) ) +  &
                   kappa_vec(2)*( (1-ddelta)*q_current + pi_current(2) )
    
    y_current    = (zf_vec**(1.0-aalpha)) * kappa_vec**aalpha * l_aggr_vec**(1.0-aalpha)
    
    w_choice_vec = (s_vec**(-1))*(1.0-aalpha)*y_current/l_aggr_vec

    ! cost of unit of bundle in terms of ch
    P_div_P_h(1) = (varsigma_vec(1)  & 
        + (1.0 -varsigma_vec(1) )*s**(sigma-1))**(1/(1-sigma))
    P_div_P_h(2) = s**(-1) * (varsigma_vec(2)*s**(1.0-sigma) + (1.0 -varsigma_vec(2)))**(1/(1-sigma))
    
    omg      = exp(state_grid(idx_omg,sss)) + omg_shift  - b_lmbd*bond_supply_shock

! =============================================================================
! NEXT PERIOD =================================================================
! =============================================================================

    ! currently assumed capital in the next period
    ! vector because depreciation shock in disaster state 
    k_nxt      = nxt_mat_2(:,1)

    ! relative productivity next period
    zf_nxt(:,1)   = 1.0_dp 
    zf_nxt(:,2)   = zeta*exp(nxt_mat_2(:,2))
    omg_nxt      = exp(nxt_mat_2(:,3)) + omg_shift  - b_lmbd*bond_supply_shock_nxt


    ! tomorrow's price, labor choice, infl and price of time endowment
    s_nxt               = nxt_mat(:,2*n_I + 1)
    q_nxt               = nxt_mat(:,2*n_I + 2) 
    l_aggr_nxt(:,1)     = nxt_mat(:,2*n_I + 3)
    l_aggr_nxt(:,2)     = nxt_mat(:,2*n_I + 4)
    infl_nxt(:,1)       = nxt_mat(:,2*n_I + 5)
    infl_nxt(:,2)       = nxt_mat(:,2*n_I + 6)
    ch_spending_nxt     = nxt_mat(:,2*n_I + 7)*exp(dz_vec)
    cf_spending_nxt     = nxt_mat(:,2*n_I + 8)*exp(dz_vec)

    homegood_infl_nxt(:,1) = infl_nxt(:,1)*( ( varsigma_vec(1) +   & 
    (1.0-varsigma_vec(1) )*    s**(sigma-1.0) ) / & 
    ( varsigma_vec(1)  + & 
    (1.0-varsigma_vec(1))*s_nxt**(sigma-1.0) ) )**(1.0/(1.0-sigma)) 
    homegood_infl_nxt(:,2) = infl_nxt(:,2)* s_nxt/s*( ( varsigma_vec(2)*    s**(1.0-sigma)+ (1.0-varsigma_vec(2)) ) / & 
                                          ( varsigma_vec(2)*s_nxt**(1.0-sigma)+ (1.0-varsigma_vec(2)) ) )**(1.0/(1.0-sigma)) 

    ! given employed labor, allocate capital to equalize capital returns across countries
    kappa_nxt(:,1) = k_nxt *  ((s_nxt)**(1.0/(1.0-aalpha))) * zf_nxt(:,1) / zf_nxt(:,2) * l_aggr_nxt(:,1)/l_aggr_nxt(:,2) / & 
                       (1.0 + ((s_nxt)**(1.0/(1.0-aalpha))) * zf_nxt(:,1) / zf_nxt(:,2) * l_aggr_nxt(:,1)/l_aggr_nxt(:,2) )
    

    kappa_nxt(:,2) = k_nxt - kappa_nxt(:,1)

    next_k       = k_nxt(1)

    investment = next_k - k_aggr*(1.0-ddelta)

    ! output in next period
    y_next(:,1) = exp((log(zf_nxt(:,1)) + dz_vec)*(1.0_dp-aalpha))*l_aggr_nxt(:,1)* & 
                  ( (kappa_nxt(:,1)/l_aggr_nxt(:,1))**aalpha)
    y_next(:,2) = exp((log(zf_nxt(:,2)) + dz_vec)*(1.0_dp-aalpha))* l_aggr_nxt(:,2)* & 
                  ( (kappa_nxt(:,2)/l_aggr_nxt(:,2))**aalpha)
    
    w_next_choice(:,1) = (1-aalpha)*y_next(:,1)/l_aggr_nxt(:,1)
    w_next_choice(:,2) = s_nxt**(-1) * (1-aalpha)*y_next(:,2)/l_aggr_nxt(:,2)

    ! profit next period per unit of capital
    pi_nxt = aalpha*exp( (log(zf_nxt(:,1)) + dz_vec)*(1.0_dp-aalpha)) &
           *(l_aggr_nxt(:,1)**(1.0_dp-aalpha))*kappa_nxt(:,1)**(aalpha-1.0)
   
    ! capital return, adjustment in disaster state
    rk_vec = ( (1.0_dp-ddelta)*q_nxt + pi_nxt)/q_current
    rk_vec(n_quad) =exp(dz_vec(n_quad))*rk_vec(n_quad) 

    P_div_P_h_nxt(:,1) = (varsigma_vec(1)  + (1.0 -varsigma_vec(1) )*s_nxt**(sigma-1))**(1.0/(1.0-sigma))
    P_div_P_h_nxt(:,2) = s_nxt**(-1) * (varsigma_vec(2)*s_nxt**(1.0-sigma) + (1.0 -varsigma_vec(2)))**(1.0/(1.0-sigma))

    ! calculate risk-free returns, given nominal rate guess 
    ! transform everything in units of home good
    rf_vec(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1)
    rf_vec(:,2) = (nom_i_vec(1) + nom_i_vec(2))/homegood_infl_nxt(:,2)
    
    if (bg_yss > 0) then 
        seigniorage_nxt = cf_spending_nxt/(ch_spending_nxt + cf_spending_nxt)*(bg_yss + bond_supply_shock_nxt)*(ch_spending_nxt + cf_spending_nxt)*omg_nxt
        seignorage_transfer(:,1) =  seigniorage_nxt 
        seignorage_transfer(:,2) = -seigniorage_nxt
    else
    seignorage_transfer = 0.0_dp
    endif

! =============================================================================
! CALCULATE MARG. UTIL AND PRICING KERNEL =====================================
! =============================================================================

    do iii = 1, n_I


          c_cost     = P_div_P_h(iii)
          c_cost_nxt = P_div_P_h_nxt(:,iii)

          wealth = aggr_wealth*wealth_share_grid(iii,sss) 

          labor      = l_aggr_vec(iii)
          
          c_spend =  minval([maxval([c_spending_vec(iii), min_cons_sav]), & 
                     wealth + w_choice_vec(iii)*labor-min_cons_sav])
          
          consumption = c_spend/c_cost
          c_vec(iii)  = consumption

          savings     = wealth + w_choice_vec(iii)*labor - c_spend 
          if (savings < 0.0) then 
              write(*,*) 'Savings:', savings
              stop
          endif
          savings_vec(iii) = savings


    enddo


! =============================================================================
! CALCULATE EXPECTED RETURNS ===================================================
! =============================================================================
E_rfh     = sum(big_weight_vec*((rf_vec(:,1))*P_div_P_h(1)/P_div_P_h_nxt(:,1))) 

E_rff     = sum(big_weight_vec*(rf_vec(:,2)*P_div_P_h(1)/P_div_P_h_nxt(:,1))) 
E_rk      = sum(big_weight_vec*(rk_vec*P_div_P_h(1)/P_div_P_h_nxt(:,1))) 
ch_vec(2) = c_spending_vec(1)/ ( (1.0-varsigma_vec(1) ) & 
             /(varsigma_vec(1) ) * s**(sigma-1.0_dp) + 1.0_dp) 
ch_vec(2) = c_spending_vec(2)/ ( (1.0-varsigma_vec(2))/varsigma_vec(2) * s**(sigma-1.0_dp) + 1.0_dp) 

consumption_exports = ch_vec(2)/P_div_P_h(1)
consumption_imports = (c_spending_vec(1) - ch_vec(1))/P_div_P_h(1)

! valuation: 1: nfa
E_kpa     = sum(big_weight_vec*kappa_nxt(:,1)/k_nxt*k_nxt(1))
nfa_tmp   = ((savings_vec(1)) - q_current*E_kpa)/P_div_P_h(1)

E_nfa_end = sum(big_weight_vec*nxt_valuation(:,2)/(rk_vec*P_div_P_h(1)/P_div_P_h_nxt(:,1)))

val_tmp  = (savings_vec(1)*share_vec(1))*((rf_vec(:,1)/(1-omg) - rk_vec)*P_div_P_h(1)/P_div_P_h_nxt(:,1)) & 
         + (savings_vec(1)*bF_share_vec(1))*(rf_vec(:,2) - rf_vec(:,1)/(1-omg))*P_div_P_h(1)/P_div_P_h_nxt(:,1) &  
         + seignorage_transfer(:,1)*P_div_P_h(1)/P_div_P_h_nxt(:,1)

E_val    = sum(big_weight_vec/(rk_vec*P_div_P_h(1)/P_div_P_h_nxt(:,1))*(val_tmp + nxt_valuation(:,3)))
nx_tmp   = consumption_exports + inv_h/P_div_P_h(1)  + q_current*((1-ddelta)*kappa_vec(1) - E_kpa)/P_div_P_h(1) - consumption_imports
E_nx     = nx_tmp +  sum(big_weight_vec/(rk_vec*P_div_P_h(1)/P_div_P_h_nxt(:,1))*nxt_valuation(:,4))
val_tmp  = (savings_vec(1)*share_vec(1))*((rf_vec(:,1)/(1-omg) - rk_vec)*P_div_P_h(1)/P_div_P_h_nxt(:,1))
E_val_rfh= sum(big_weight_vec/(rk_vec*P_div_P_h(1)/P_div_P_h_nxt(:,1))*(val_tmp + nxt_valuation(:,5)))
val_tmp  = (savings_vec(1)*bF_share_vec(1))*(rf_vec(:,2) - rf_vec(:,1)/(1-omg))*P_div_P_h(1)/P_div_P_h_nxt(:,1)
E_val_rff= sum(big_weight_vec/(rk_vec*P_div_P_h(1)/P_div_P_h_nxt(:,1))*(val_tmp + nxt_valuation(:,6)))
val_tmp  = seignorage_transfer(:,1) *P_div_P_h(1)/P_div_P_h_nxt(:,1)

E_val_seign= sum(big_weight_vec/(rk_vec*P_div_P_h(1)/P_div_P_h_nxt(:,1))*(val_tmp + nxt_valuation(:,7)))

valuation_vec = [nfa_tmp, E_nfa_end, E_val, E_nx, E_val_rfh, E_val_rff, E_val_seign]
! valuation vec is: 1. current nfa position using expected kappa 
!                   2. discounted final nfa position (needs to be assigne from 1 in first iteration) 
!                   3. discounted valuation term  
!                   4. discounted nx term (but with current nx, needs to be shifted by one period) 
!                   5. discounted val rf term  
!                   6. discounted val rff term  


end subroutine calc_valuation 

subroutine calc_excess_bond_nom( aggr_wealth, omg, v_vecs, rf_vec, rk, seignorage_transfer, & 
                                 bbeta_adj, c_cost_vec, c_cost_nxt_vec,  &  
                                 q_current, w_choice_vec, l_vec,  sss, big_weight_vec, &
                                 c_temp, c_spend_temp,  share_guess_vec, bF_share_guess_vec, & 
                                 k_temp_vec, b_temp_vec, excess_b, constraint_binding_vec, mc_nxt_mat)
                                 
use mod_param, only: n_quad, n_I, low_guess_fixed, high_guess_fixed, omg_shift, b_lmbd, state_grid, &
                     wealth_share_grid, bg_yss, idx_omg, &
                     foreign_trading
                     

integer, intent(in)  :: sss
real(dp), intent(in) :: aggr_wealth, omg, v_vecs(n_quad, n_I),  bbeta_adj(n_I), & 
                        rf_vec(n_quad, 2), rk(n_quad), big_weight_vec(n_quad),   &  
                        q_current, w_choice_vec(2), l_vec(2),  seignorage_transfer(n_quad, n_I), &
                        c_cost_vec(2), c_cost_nxt_vec(n_quad, 2), c_temp(n_I), c_spend_temp(n_I), mc_nxt_mat(n_quad, n_I) 
real(dp), intent(inout) ::  share_guess_vec(n_I), bF_share_guess_vec(n_I)
real(dp), intent(out)   :: k_temp_vec(n_I), b_temp_vec(n_I), excess_b
integer, intent(out)    :: constraint_binding_vec(n_I)

real(dp) :: rf(n_quad)
real(dp) :: share_guess, share_low, share_high, wealth, temp(n_quad), savings, objf, &
            labor, derivative_consumption, wage
integer  :: iii


 
! for each agent calculate portfolio choice between bonds and capital 
do iii = 1,n_I

    
    wage = w_choice_vec(iii)

    ! marginal portfolio choice over home bond 
    rf = rf_vec(:,1)/(1.0-omg)

    wealth = aggr_wealth*wealth_share_grid(iii,sss) 

    ! agent works
        labor  = l_vec(iii)
 
    ! savings, consumption choice held fixed 
    savings =  wealth + wage*labor - c_spend_temp(iii)

    ! check for negative savings 
    if (savings <0.0_dp) then
       write(*,*) 'ERROR: negative savings in calc_excess_bond_nom'
       stop
    endif

        ! find leverage lower limit to prevent negative payoffs: (a-b)*r1 + b*r2 + (1-a)*rk >0 -rk 
        temp = (-rk - bF_share_guess_vec(iii)*(rf_vec(:,2)-rf_vec(:,1)))/(rf-rk + sqrt_eps)
        where (temp >= 0.0_dp)
            temp = -10000.0_dp
        end where
        
        ! max of limit and imposed limit
        share_low = maxval([low_guess_fixed, maxval(temp)+ sqrt_eps])

        ! find upper natural leverage limit to prevent negative payoffs
        temp = (-rk - bF_share_guess_vec(iii)*(rf_vec(:,2)-rf_vec(:,1)))/(rf-rk + sqrt_eps)
        where (temp <= 0.0_dp)
            temp = 10000.0_dp
        end where

        ! min of limit and imposed limit
        share_high = minval([high_guess_fixed, minval(temp) - sqrt_eps ])
        
        ! ensure initial guess within boundaries
        share_guess = maxval([minval([share_guess_vec(iii), share_high]), share_low]) ! start from share guess but stay in limits
   
            ! calculate portfolio share 

        call calc_portfolio_share(share_low, share_high, iii, sss, big_weight_vec, savings, labor, c_spend_temp(iii) / c_cost_vec(iii), mc_nxt_mat(:, iii),  &
        v_vecs(:,iii), rf_vec, rk, seignorage_transfer(:,iii), omg,  & 
        bbeta_adj, c_cost_vec(iii), c_cost_nxt_vec(:, iii), &
        share_guess, bF_share_guess_vec(iii), constraint_binding_vec(iii))
    
        ! store results: can everyone adjust portfolio?
        share_guess_vec(iii) = share_guess 

        k_temp_vec(iii) = savings*(1.0_dp-share_guess)/q_current
        b_temp_vec(iii) = savings*(( (share_guess-bF_share_guess_vec(iii)) & 
                                   + bF_share_guess_vec(iii) ) )
           
 

enddo

! calculate excess bond demand
excess_b =  sum(b_temp_vec) 

end subroutine calc_excess_bond_nom

! calculates the optimal portfolio between capital and bonds given a savings decision
subroutine calc_portfolio_share(share_low_in, share_high_in, iii, zzz, big_weight_vec, savings, labor, consumption, mc_nxt, &
    v, rf_vec, rk_vec, seignorage_transfer, omg, bbeta_adj, c_cost, c_cost_nxt, share_guess, bF_share, constraint_binding)

    use mod_param, only: gma_vec, IES_vec, idx_dis, disast_p, &
              bbeta_vec, n_quad, dz_vec, state_grid, disast_std, disast_p, idx_omg,  &
              omg_shift, b_lmbd, n_I
              

    ! worker_capital_SSfrac, 
    real(dp), intent(in)    :: v(n_quad), share_low_in, share_high_in, &
                    rf_vec(n_quad,2), rk_vec(n_quad), seignorage_transfer(n_quad), savings, labor, c_cost, omg, & 
                    c_cost_nxt(n_quad), big_weight_vec(n_quad), consumption, bF_share, bbeta_adj(n_I), mc_nxt(n_quad)
    integer, intent(in)     :: iii, zzz
    real(dp), intent(inout) :: share_guess
    integer, intent(out)    :: constraint_binding

    real(dp) :: IES, share_FOC_low, share_FOC_hi, share_guess_A, share_guess_B, share_guess_C, share_guess_S, share_guess_D, share_FOC, &
     share_FOC_A, share_FOC_B, share_FOC_C, share_FOC_S, temp_vec(n_quad), temp, brent_delta, rk(n_quad), &
     objgrd(2), r_alpha(n_quad), share_temp, share_low, share_high, cons_update, gma, rk_FOC(n_quad), rf_home(n_quad), rf_foreign(n_quad), &
     next_period_share(n_quad), v_vec_twisted(n_quad), objf, rf(n_quad), bbeta, chi0, EV, c_cost_growth(n_quad)

    integer :: mflag, s_iter   

        share_low  = share_low_in 

    constraint_binding = 0

    share_high = share_high_in

    ! for risk on calculations
    gma   = gma_vec(iii)
    IES = IES_vec(iii)
    
    c_cost_growth = c_cost_nxt/c_cost 

    rf = rf_vec(:,1)/(1-omg)
    rk = rk_vec
    rk_FOC = rk
    rf_home = rf_vec(:,1)/(1.0 - omg)
    rf_foreign = rf_vec(:,2)

     if (share_low<-20.0_dp) then

        share_temp = -20.0_dp
        r_alpha = portfolio_return(share_temp, bF_share, rf_home, rf_foreign, rk) 
        next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                            (tot_wealth_ss_vec(iii))
        share_FOC_low = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf, rk_FOC, c_cost_growth, gma, IES) 
    
        if (share_FOC_low > 0.0_dp) then
            share_low = -20.0_dp
        else
            share_temp = share_low
        endif

    else
        share_temp = share_low
    endif

    r_alpha = portfolio_return(share_temp, bF_share, rf_home, rf_foreign, rk) 
    next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                        (tot_wealth_ss_vec(iii))
    share_FOC_low = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf, rk_FOC, c_cost_growth, gma, IES) 


    if (share_FOC_low <= 0.0_dp) then

        share_guess = share_low
        constraint_binding = 1.0

    else

        share_temp = share_high
        r_alpha = portfolio_return(share_temp, bF_share, rf_home, rf_foreign, rk) 
        next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                            (tot_wealth_ss_vec(iii))
        share_FOC_hi = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf, rk_FOC, c_cost_growth, gma, IES) 

        if (share_FOC_hi >= 0.0_dp) then
            share_guess = share_high
            constraint_binding = -1
        else
            share_temp = share_guess
            r_alpha = portfolio_return(share_temp, bF_share, rf_home, rf_foreign, rk) 
            next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                                (tot_wealth_ss_vec(iii))
            share_FOC = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf, rk_FOC, c_cost_growth, gma, IES) 

            if (share_FOC /= 0.0_dp) then

                share_guess_B = share_guess ! make this highest value
                share_FOC_B   = share_FOC

                if (share_FOC <= 0.0_dp) then ! share_FOC<0 means optimum lies to the left

                    share_guess_A = share_low
                    share_temp = share_guess_A
                    r_alpha = portfolio_return(share_temp, bF_share, rf_home, rf_foreign, rk) 
                    next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                                        (tot_wealth_ss_vec(iii))
                    share_FOC_A = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf, rk_FOC, c_cost_growth, gma, IES) 
                 else ! share_FOC > 0 so optimum lies to the right

                    share_guess_A = share_high
                    share_temp = share_guess_A
                    r_alpha = portfolio_return(share_temp, bF_share, rf_home, rf_foreign, rk) 
                    next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                                        (tot_wealth_ss_vec(iii))
                    share_FOC_A = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf, rk_FOC, c_cost_growth, gma, IES) 

                endif


                ! check that root is bracketed
                if ((share_FOC_A*share_FOC_B) > 0) then
                   write(*,*) 'ERROR: Initial bracket does not contain root. D'
                   stop
                endif

                ! swap a and b
                if ( abs(share_FOC_A) < abs(share_FOC_B) )then
                    temp = share_guess_A
                    share_guess_A = share_guess_B
                    share_guess_B = temp

                    temp = share_FOC_A
                    share_FOC_A = share_FOC_B
                    share_FOC_B = temp
                endif


                share_guess_C    = share_guess_A
                share_FOC_C      = share_FOC_A
                share_FOC_S      = share_FOC_B
                mflag = 1 ! set mflag
                brent_delta = 5E-16_dp
                s_iter = 1 
                do while ( share_FOC_S /= 0 .and. abs(share_guess_A - share_guess_B) > 1E-14_dp)
                    s_iter = s_iter+1 

                    if ( (share_FOC_A /= share_FOC_C) .and. (share_FOC_C /= share_FOC_B) ) then ! inverse quadratic interpolation

                        share_guess_S = share_guess_A * share_FOC_B * share_FOC_C / ( (share_FOC_A - share_FOC_B) * (share_FOC_A - share_FOC_C) ) + &
                                        share_guess_B * share_FOC_A * share_FOC_C / ( (share_FOC_B - share_FOC_A) * (share_FOC_B - share_FOC_C) ) + &
                                        share_guess_C * share_FOC_A * share_FOC_B / ( (share_FOC_C - share_FOC_A) * (share_FOC_C - share_FOC_B) )

                    else ! secant method

                        share_guess_S = share_guess_B - share_FOC_B * (share_guess_B - share_guess_A) /   (share_FOC_B - share_FOC_A)

                    endif

                    if ( ( (  ( share_guess_S > ( 3*share_guess_A + share_guess_B )/4 .and. share_guess_S < share_guess_B) .or. &
                               ( share_guess_S < ( 3*share_guess_A + share_guess_B )/4 .and. share_guess_S > share_guess_B)  ) == .FALSE. ) .or. &
                          (mflag == 1 .and. abs(share_guess_S - share_guess_B) >= abs(share_guess_B - share_guess_C)/2  )             .or. &
                          (mflag == 0 .and. abs(share_guess_S - share_guess_B) >= abs(share_guess_B - share_guess_D)/2  )             .or. &
                          (mflag == 1 .and. abs(share_guess_B - share_guess_C) <  abs(brent_delta)  )                     .or. &
                          (mflag == 0 .and. abs(share_guess_B - share_guess_D) <  abs(brent_delta)  )                          &
                        ) then

                       share_guess_S = (share_guess_A + share_guess_B )/ 2
                       mflag = 1
                    else
                        mflag = 0
                    endif

                    share_temp = share_guess_S
                    r_alpha = portfolio_return(share_temp, bF_share, rf_home, rf_foreign, rk) 
                    next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                                        (tot_wealth_ss_vec(iii))
                    share_FOC_S = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf, rk_FOC, c_cost_growth, gma, IES) 
                             

                    share_guess_D = share_guess_C
                    share_guess_C = share_guess_B
                    share_FOC_C = share_FOC_B

                    if ((share_FOC_A*share_FOC_S) < 0) then
                       share_guess_B = share_guess_S
                       share_FOC_B = share_FOC_S
                    else
                       share_guess_A = share_guess_S
                       share_FOC_A = share_FOC_S
                    endif

                    ! swap a and b
                    if ( abs(share_FOC_A) < abs(share_FOC_B) ) then

                        temp = share_guess_A
                        share_guess_A = share_guess_B
                        share_guess_B = temp

                        temp = share_FOC_A
                        share_FOC_A = share_FOC_B
                        share_FOC_B = temp

                    endif

                    if (s_iter > 1000) then 
                         write(*,*) 'no share convergence'
                         stop
                    endif
                enddo

            share_guess = share_guess_B

        endif

    endif
endif

end subroutine calc_portfolio_share



subroutine calc_excess_foreign_bond_nom(aggr_wealth, omg, v_vecs, rf_vec, rk, seignorage_transfer, & 
                                 bbeta_adj, c_cost_vec, c_cost_nxt_vec,  &  
                                 w_choice_vec, l_vec, sss, big_weight_vec, mc_nxt_mat, &
                                 c_temp, c_spend_temp, share_guess_vec, bF_share_guess_vec, & 
                                 bh_temp_vec, bf_temp_vec, excess_b_home, constraint_binding_vec)

use mod_param, only: n_quad, n_I, safe_low_guess_fixed, safe_high_guess_fixed, &
                     wealth_share_grid, bg_yss, state_grid, idx_omg, &
                     omg_shift, b_lmbd

integer, intent(in)  :: sss
real(dp), intent(in) :: aggr_wealth, omg, v_vecs(n_quad, n_I), & 
                        rf_vec(n_quad, 2), rk(n_quad), big_weight_vec(n_quad), &  
                        w_choice_vec(2), l_vec(2), share_guess_vec(n_I), seignorage_transfer(n_quad, n_I), &
                        c_cost_vec(2), c_cost_nxt_vec(n_quad, 2), c_temp(n_I), c_spend_temp(n_I), bbeta_adj(n_I), mc_nxt_mat(n_quad, n_I)
real(dp), intent(inout) ::   bF_share_guess_vec(n_I)
real(dp), intent(out)   :: bh_temp_vec(n_I), bf_temp_vec(n_I), excess_b_home
integer, intent(out)    :: constraint_binding_vec(n_I)

real(dp) :: rf(n_quad), temp_bound(n_quad)
real(dp) :: share_guess, share_low, share_high, wealth, temp(n_quad), savings, objf, &
            labor, derivative_consumption, wage, temp_r(n_quad), bound1, bound2
integer  :: iii

    
do iii = 1,n_I

    wage = w_choice_vec(iii)

    wealth = aggr_wealth*wealth_share_grid(iii,sss)

    ! agent works?
        labor  = l_vec(iii)

    ! savings (consumption decision held fixed)
    savings =  wealth + wage*labor - c_spend_temp(iii)

    ! check if positive savings  
    if (savings <0.0_dp) then
       write(*,*) 'ERROR: negative savings calc_excess_foreign_bond_nom'
       write(*,*) savings, sss
       stop
    endif
    
        temp_r = rf_vec(:,2) - rf_vec(:,1)/(1.0 - omg)
       
        temp_bound = -1000.0_dp
        where (temp_r > 0)
        temp_bound  =  (-(1-share_guess_vec(iii))*rk - share_guess_vec(iii)*rf_vec(:,1))/&
                              temp_r
        endwhere



        bound1  =  maxval( temp_bound )

        temp_bound = 1000.0_dp
        where (temp_r < 0)
        temp_bound  =  (-(1-share_guess_vec(iii))*rk - share_guess_vec(iii)*rf_vec(:,1))/&
                              temp_r
        endwhere
        
        bound2  =  minval( temp_bound )

     share_low = maxval([safe_low_guess_fixed, bound1 + 0.01])

     share_high = minval([safe_high_guess_fixed, bound2 - 0.01])

      share_guess =  maxval([minval([bF_share_guess_vec(iii), share_high]), share_low]) ! start from share guess but stay in limits
    
        call calc_bond_portfolio_share(share_low, share_high, iii, sss, big_weight_vec, savings, labor, c_spend_temp(iii) / c_cost_vec(iii), mc_nxt_mat(:, iii),  &
        v_vecs(:,iii), rf_vec(:,1)/(1.0-omg),   rf_vec(:,2), &
        rk, seignorage_transfer(:,iii), omg, bbeta_adj, c_cost_vec(iii), c_cost_nxt_vec(:,iii), &
        share_guess, share_guess_vec(iii), constraint_binding_vec(iii) )

        
     


       ! store results
       bF_share_guess_vec(iii) = share_guess
       bf_temp_vec(iii) = savings*share_guess
       bh_temp_vec(iii) = (share_guess_vec(iii) - share_guess)*savings

enddo

! calculate excess bond demand
excess_b_home =   -1.0_dp*sum(bh_temp_vec)

end subroutine calc_excess_foreign_bond_nom

! calculates the optimal portfolio between capital and bonds given a savings decision
subroutine calc_bond_portfolio_share(share_low_in, share_high_in, iii, zzz, big_weight_vec, savings, labor, consumption, mc_nxt,  &
           v, rf_home, rf_foreign, rk, seignorage_transfer, omg, bbeta_adj, c_cost, c_cost_nxt, share_guess, bond_share, constraint_binding)
 
use mod_param, only: gma_vec, IES_vec, v_normalization_vec, idx_dis, idx_dis, &
                     bbeta_vec, n_quad, dz_vec, state_grid, disast_std, disast_p, idx_omg, omg_shift, b_lmbd

! worker_capital_SSfrac, 
real(dp), intent(in)    :: v(n_quad), share_low_in, share_high_in, &
                           rf_home(n_quad), rf_foreign(n_quad), rk(n_quad), seignorage_transfer(n_quad), savings, labor, c_cost, & 
                           c_cost_nxt(n_quad), big_weight_vec(n_quad), consumption, bond_share, omg, bbeta_adj(n_I), mc_nxt(n_quad)
integer, intent(in)     :: iii, zzz
real(dp), intent(inout) :: share_guess
integer, intent(out)    :: constraint_binding

real(dp) :: IES, share_FOC_low, share_FOC_hi, share_guess_A, share_guess_B, share_guess_C, share_guess_S, share_guess_D, share_FOC, &
            share_FOC_A, share_FOC_B, share_FOC_C, share_FOC_S, temp_vec(n_quad), temp, brent_delta, &
            objgrd(2), r_alpha(n_quad), share_temp, share_low, share_high, cons_update, gma, &
            next_period_share(n_quad), v_vec_twisted(n_quad), objf, EV, c_cost_growth(n_quad)

integer :: mflag, s_iter
    
    share_low  = share_low_in 
    share_high = share_high_in
    constraint_binding = 0
 
        gma   = gma_vec(iii)
    IES = IES_vec(iii)

    c_cost_growth = c_cost_nxt/c_cost 

    share_temp = share_low
    r_alpha = portfolio_return(bond_share, share_temp, rf_home, rf_foreign, rk)
    next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                        (tot_wealth_ss_vec(iii))
    share_FOC_low = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf_foreign, rf_home, c_cost_growth, gma, IES) 
     
    if (share_FOC_low == 0.0_dp) then            
        share_guess = share_low                  
    elseif (share_FOC_low <= 0.0_dp) then
        share_guess = share_low
        constraint_binding = -1
    else
        share_temp = share_high
        r_alpha = portfolio_return(bond_share, share_temp, rf_home, rf_foreign, rk)
        next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                            (tot_wealth_ss_vec(iii))
        share_FOC_hi = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf_foreign, rf_home, c_cost_growth, gma, IES) 
 
        if (share_FOC_hi == 0.0_dp) then
            share_guess = share_high
        elseif (share_FOC_hi >= 0.0_dp) then
            share_guess = share_high
            constraint_binding = 1
        else

            share_temp = share_guess
            r_alpha = portfolio_return(bond_share, share_temp, rf_home, rf_foreign, rk)
            next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                                (tot_wealth_ss_vec(iii))
            share_FOC = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf_foreign, rf_home, c_cost_growth, gma, IES) 

 
            if (share_FOC /= 0.0_dp) then
 
                share_guess_B = share_guess ! make this highest value
                share_FOC_B   = share_FOC
 
                if (share_FOC <= 0.0_dp) then ! share_FOC<0 means optimum lies to the left

                    share_guess_A = share_low
                else ! share_FOC > 0 so optimum lies to the right
                     share_guess_A = share_high
                endif

                share_temp = share_guess_A
                r_alpha = portfolio_return(bond_share, share_temp, rf_home, rf_foreign, rk)
                next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                                    (tot_wealth_ss_vec(iii))
                share_FOC_A = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf_foreign, rf_home, c_cost_growth, gma, IES) 
 
                ! check that root is bracketed
                if ((share_FOC_A*share_FOC_B) > 0) then
                    write(*,*) 'ERROR: Initial bracket does not contain root. C'
                    stop
                endif
 
                ! swap a and b
                if ( abs(share_FOC_A) < abs(share_FOC_B) )then
                    temp = share_guess_A
                    share_guess_A = share_guess_B
                    share_guess_B = temp
 
                    temp = share_FOC_A
                    share_FOC_A = share_FOC_B
                    share_FOC_B = temp
                endif
 
 
               share_guess_C    = share_guess_A
               share_FOC_C      = share_FOC_A
               share_FOC_S      = share_FOC_B
               mflag = 1 ! set mflag
               brent_delta = 5E-16_dp
               s_iter = 0  
               do while ( share_FOC_S /= 0 .and. abs(share_guess_A - share_guess_B) > 1E-14_dp)
                    s_iter = s_iter + 1

                 if ( (share_FOC_A /= share_FOC_C) .and. (share_FOC_C /= share_FOC_B) ) then ! inverse quadratic interpolation
 
                   share_guess_S = share_guess_A * share_FOC_B * share_FOC_C / ( (share_FOC_A - share_FOC_B) * (share_FOC_A - share_FOC_C) ) + &
                             share_guess_B * share_FOC_A * share_FOC_C / ( (share_FOC_B - share_FOC_A) * (share_FOC_B - share_FOC_C) ) + &
                             share_guess_C * share_FOC_A * share_FOC_B / ( (share_FOC_C - share_FOC_A) * (share_FOC_C - share_FOC_B) )
 
                 else ! secant method
 
                   share_guess_S = share_guess_B - share_FOC_B * (share_guess_B - share_guess_A) /   (share_FOC_B - share_FOC_A)
 
                 endif
 
                 if ( ( (  ( share_guess_S > ( 3*share_guess_A + share_guess_B )/4 .and. share_guess_S < share_guess_B) .or. &
                           ( share_guess_S < ( 3*share_guess_A + share_guess_B )/4 .and. share_guess_S > share_guess_B)  ) == .FALSE. ) .or. &
                      (mflag == 1 .and. abs(share_guess_S - share_guess_B) >= abs(share_guess_B - share_guess_C)/2  )             .or. &
                      (mflag == 0 .and. abs(share_guess_S - share_guess_B) >= abs(share_guess_B - share_guess_D)/2  )             .or. &
                      (mflag == 1 .and. abs(share_guess_B - share_guess_C) <  abs(brent_delta)  )                     .or. &
                      (mflag == 0 .and. abs(share_guess_B - share_guess_D) <  abs(brent_delta)  )                          &
                    ) then
 
                   share_guess_S = (share_guess_A + share_guess_B )/ 2
                   mflag = 1
                 else
                   mflag = 0
                 endif

                 share_temp = share_guess_S
                 r_alpha = portfolio_return(bond_share, share_temp, rf_home, rf_foreign, rk)
                 next_period_share = (savings*r_alpha   + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer )/ & 
                                     (tot_wealth_ss_vec(iii))
                 share_FOC_S = portfolio_foc(next_period_share, v, mc_nxt, big_weight_vec, rf_foreign, rf_home, c_cost_growth, gma, IES) 
                 
                 share_guess_D = share_guess_C
                 share_guess_C = share_guess_B
                 share_FOC_C = share_FOC_B
 
                 if ((share_FOC_A*share_FOC_S) < 0) then
                   share_guess_B = share_guess_S
                   share_FOC_B = share_FOC_S
                 else
                   share_guess_A = share_guess_S
                   share_FOC_A = share_FOC_S
                 endif
 
                 ! swap a and b
                 if ( abs(share_FOC_A) < abs(share_FOC_B) ) then
 
                    temp = share_guess_A
                    share_guess_A = share_guess_B
                    share_guess_B = temp
 
                    temp = share_FOC_A
                    share_FOC_A = share_FOC_B
                    share_FOC_B = temp
 
                 endif

                    if (s_iter > 1000) then 
                        write(*,*) 'no share convergence b'
                        write(*,*) iii
                        write(*,*) bond_share 
                        write(*,*) share_FOC_A, share_FOC_B
                        write(*,*) share_guess_A, share_guess_B, share_guess_S, share_guess_C
                        stop
                    endif

               enddo
 
           share_guess = share_guess_B
 
         endif
 
       endif
     endif
 
end subroutine calc_bond_portfolio_share


! utility function and derivative w.r.t. to consumption
subroutine util_fun(consumption, labor, iii, util, util_c_deriv, labor_part) 
use mod_param, only: chi0_vec, ies_vec, chi

integer, intent(in)   :: iii
real(dp), intent(in)  :: consumption, labor
real(dp), intent(out) :: util, util_c_deriv, labor_part

real(dp) :: chi0, ies 

    chi0 = chi0_vec(iii)
    ies  = ies_vec(iii)

    if (IES > 1.0_dp-sqrt_eps .and. IES < 1.0_dp + sqrt_eps) then 
        labor_part = 1.0_dp
        util = log(consumption) -  chi0*chi/(1.0_dp + chi) * labor**( (1.0_dp+chi)/chi ) 
        util_c_deriv = 1.0_dp/consumption    
    else 
        labor_part =  (1.0_dp + (1.0_dp/IES - 1.0_dp) * chi0*chi/(1.0_dp + chi) * labor**( (1.0_dp+chi)/chi ) )**(1.0_dp/IES)
        util = (consumption**(1.0_dp-1.0_dp/IES)) * labor_part
        util_c_deriv = util/consumption    
    endif

end subroutine util_fun

function portfolio_foc(next_period_wealth, v, mc, big_weight_vec, r1, r2, c_cost_growth, gma, ies) result(foc)
    
    use mod_param, only: n_quad

    real(dp), intent(in)  :: next_period_wealth(n_quad), big_weight_vec(n_quad), v(n_quad), mc(n_quad),  &
                             r1(n_quad), r2(n_quad), c_cost_growth(n_quad), gma, ies
    real(dp) :: foc, sdf_component(n_quad)
    
    sdf_component = (v**(1.0/ies - gma)) * mc
    foc = sum(big_weight_vec * sdf_component * (next_period_wealth**(-gma)) * (r1 - r2) / c_cost_growth) &
              /abs(sum(big_weight_vec*sdf_component)) ! just for scaling 

end function

function portfolio_return(bond_share, foreign_share, rh, rf, rk) result(r_alpha)

    use mod_param, only: n_quad
    real(dp), intent(in)  :: bond_share, foreign_share, rh(n_quad), rf(n_quad), rk(n_quad)
    real(dp) :: r_alpha(n_quad)

    r_alpha = (bond_share - foreign_share)*rh + foreign_share*rf + (1.0_dp - bond_share)*rk

end function


! calculates consumption bundle 
function calc_bundle(c_h, c_f, varsigma) result(c_bundle) 
    use mod_param, only : sigma
    real(dp) :: c_bundle
    real(dp), intent(in)  :: c_h, c_f, varsigma

    c_bundle = ( varsigma**(1.0/sigma) * c_h**((sigma-1.0)/sigma) + & 
                    (1.0-varsigma)**(1.0/sigma) * c_f**((sigma-1.0)/sigma))**(sigma/(sigma-1.0))

end function calc_bundle

subroutine calc_equilibrium_ifixed( nxt_mat, n_nxt, nxt_mat_2, sss, outer_iter, &
                            c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                            nom_i_vec, share_vec, share_vec_new, bF_share_vec, inv_h, inv_f, bond_supply_shock, bond_supply_shock_nxt, &
                            q_new, s_new, k_next_new, c_spending_vec_new, infl_new, &
                            v_new, mc_new, l_aggr_new, w_choice_new, theta_nxt, M_vec_mat, excess_b_out)
    
    use mod_param, only:    n_I, n_quad,  idx_zf, idx_k, aalpha, smolyak_d, quad_weight_vec, zeta, bbeta_coeff, &
                            dz_vec, ddelta, state_grid, varsigma_vec, sigma, idx_dis, inv_share_h, omg_rho, idx_thh, &
                            chix, ies_vec, gma_vec, omg_shift, b_lmbd, bg_yss, &
                            bbeta_vec, v_normalization_vec, chi0_vec, chi, wealth_share_grid, &
                            vareps_w, n_interp, &
                            phi_w, idx_wh, idx_wf, phi_h, phi_f, tayl_ic_h, tayl_ic_f, tht_h_grid_mean, & 
                            idx_ih, idx_if, rho_i, disast_p, disast_std, idx_omg, &
                            foreign_trading, &
                            omg_grid_mean, phi_yh, phi_yf 

    real(dp), intent(in) :: nxt_mat_2(n_quad, 3), l_aggr_vec(2), &
                            nxt_mat(n_quad, n_nxt), share_vec(n_I),s, &
                            q_current, c_spending_vec(n_I), infl_vec(2), bond_supply_shock, bond_supply_shock_nxt
    real(dp), intent(inout) :: bF_share_vec(n_I),  nom_i_vec(2), inv_f, inv_h, M_vec_mat(n_quad, n_I)
    real(dp), intent(out)   :: q_new, s_new, k_next_new, c_spending_vec_new(n_I), infl_new(2), excess_b_out, &
                               v_new(n_I), l_aggr_new(2), share_vec_new(n_I),w_choice_new(2), theta_nxt(n_quad), mc_new(n_I)

    integer, intent(in)  :: sss, n_nxt, outer_iter

    integer :: iii, ggg, ifail, nnn

    real(dp) :: s_vec(2), zf_vec(2), k_aggr, kappa_vec(2), & 
                q_nxt(n_quad), rk_vec(n_quad), &
                aggr_wealth, v_temp(n_quad, n_I),  mc_temp(n_quad, n_I), &
                w_next_choice(n_quad, 2), & 
                pi_nxt(n_quad), y_next(n_quad, 2), &
                kappa_nxt(n_quad,2), homegood_infl_nxt(n_quad, 2), &
                k_nxt(n_quad), s_nxt(n_quad), zf_nxt(n_quad,2), &
                rf_vec(n_quad, 2), pi_current(2), y_current(2), &
                k_temp_vec(n_I), b_temp_vec(n_I), c_temp(n_I), excess_b, & 
                excess_b_home, bf_temp_vec(n_I), bh_temp_vec(n_I), & 
                share_temp(n_I), bF_share_temp(n_I), infl_nxt(n_quad,2), &
                c_vec(n_I), ch_vec(n_I), chi0, savings_vec(n_I), &
                labor_nxt(n_quad)

    real(dp) :: next_k

    real(dp), parameter :: brent_delta = 5E-16_dp, min_cons_sav = 1E-08_dp
 
    real(dp) :: &
             nom_i_a, nom_i_b, rf_vec_temp(n_quad,2), excess_b_A, excess_b_B, w_choice_vec(2), &
             nom_i_C, temp, excess_b_C, excess_b_S, nom_i_S, nom_i_D, wealth, consumption, r_alpha_vec(n_quad), &
             diff_inner, iter_inner, savings, objf, objgrd(2), cons_update, IES, bbeta, &
             v_normalization, next_period_share(n_quad), v_vec(n_quad), EV, &
             weight_vec(n_I), w_temp, w_diff, w_temp_new, temp_vec(n_quad), & 
             share_low, share_high, real_bond, c_spend_temp(n_I), &
             mid_point, k_param, derivative_consumption, M_weighted_vec(n_quad), nom_i_guess, &
             share_guess, gov_debt,  gma, bbeta_adj(n_I), omg_nxt(n_quad), &
             ch, cf, c_spend, labor, c_cost, c_cost_nxt(n_quad), & 
             v_vec_twisted(n_quad), ih_last, if_last, iota(n_quad), & 
             l_aggr_nxt(n_quad,2), nxt_wealth_vec(n_quad,n_I), &
             s_update, diff, P_div_P_h(2), P_div_P_h_nxt(n_quad,2), w_current_next(n_quad), &
             w_current(2), w_choice_new_vec(2), kappa_guess, kappa_share, util, util_c_deriv, labor_part, &
             l_aggr_vec_temp(2), l_aggr_temp(2), l_aggr_temp_new(2), w_choice_temp(2), w_choice_temp_new(2), &
             kappa_vec_new(2), p_dis, big_weight_vec(n_quad), deriv_helper(n_quad), labor_new, &
             E_rfh, E_rff, E_rk, investment, q_temp, test_ratio,  l_temp_vec(n_I), & 
             temp_bond_holdings, idiff, nom_i_vec_in(2), omg, r_alpha_omg_vec(n_quad), &
             ch_spending_nxt(n_quad), cf_spending_nxt(n_quad), seigniorage_nxt(n_quad),& 
             tax_vec(n_quad), seignorage_transfer(n_quad, n_I), wealth_vec(n_I)

 
 integer :: iter, mflag, ix1, ix2, safe_constraint_binding_vec(n_I), constraint_binding_vec(n_I), bond_pricer, iiter


    nom_i_vec_in = nom_i_vec ! store incoming nominal rate guess for Taylor rule

! =============================================================================
! CURRENT STATE ===============================================================
! =============================================================================

    ! current capital stock
    k_aggr        = state_grid(idx_k,   sss)

    ! store relative productivity in vector, foreign with relative productivity shock 
    zf_vec = [ 1.0_dp, zeta*exp(state_grid(idx_zf, sss))]

    ! current disaster probability
    p_dis       = exp(state_grid(idx_dis, sss))
    

    ! current weight vec, given disaster probability
    big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
    big_weight_vec(n_quad)     = p_dis

    ! w_current is stored as w / (P /P_h) and w / (P* / P_h*), so in bundle units
    ! allows to calculate w'/w*P'_h/P_h = w'/w_current*infl*P'_h/P'
    w_current = [state_grid(idx_wh, sss), state_grid(idx_wf, sss)]

! =============================================================================
! CURRENT PERIOD GUESSES ======================================================
! =============================================================================

    ! store relative prices home/foreign good in vector
    s_vec   = [1.0_dp, s]

    ! given current guess of aggregate labor and wage, allocate capital across countries
    ! allows to calculate initial wealth 
    kappa_vec(1) = k_aggr * ( (s_vec(2) )**(1.0/(1.0-aalpha))) * &
                      zf_vec(1) / zf_vec(2) * l_aggr_vec(1)/l_aggr_vec(2) / & 
                   (1.0 +   ( (s_vec(2) )**(1.0/(1.0-aalpha))) * &
                      zf_vec(1) / zf_vec(2) * l_aggr_vec(1)/l_aggr_vec(2) )
    
    kappa_vec(2) = k_aggr - kappa_vec(1)

    ! given this calculate profits and aggregate wealth
    pi_current   = (s_vec**(-1))*aalpha * (zf_vec*l_aggr_vec/kappa_vec)**(1.0-aalpha)
   
    ! current output 
    y_current    = (zf_vec**(1.0-aalpha)) * kappa_vec**aalpha * l_aggr_vec**(1.0-aalpha)
    
    ! wages calculated in home consumption goods
    w_choice_vec = (s_vec**(-1))*(1.0-aalpha)*y_current/l_aggr_vec

    ! cost of unit of bundle in terms of home good (P/P_h and P*/P^*_h 
    ! (one bundle unit valued at P, divided by cost of home good -> units of home good per bundle
    P_div_P_h(1) = (    varsigma_vec(1)  & 
                        + (1.0 - varsigma_vec(1) )*(s**(sigma-1.0))  )**(1.0/(1.0-sigma))
    P_div_P_h(2) = s**(-1) * (varsigma_vec(2)*s**(1.0-sigma) + (1.0 -varsigma_vec(2)))**(1.0/(1.0-sigma))
    
    aggr_wealth  = kappa_vec(1)*( (1-ddelta)*q_current + pi_current(1) ) +  &
                   kappa_vec(2)*( (1-ddelta)*q_current + pi_current(2) )
    wealth_vec = aggr_wealth * wealth_share_grid(:, sss)               
    

    ! set the omg - with rebating only appears in taylor rule
    omg      = exp(state_grid(idx_omg,sss)) + omg_shift  - b_lmbd*bond_supply_shock
   
    

    
    ! beta adjustment for stationarity if bbeta_coeff > 0 
    bbeta_adj(1) = - bbeta_coeff*(state_grid(idx_thh,sss) - tht_h_grid_mean)
    bbeta_adj(2) = + bbeta_coeff*(state_grid(idx_thh,sss) - tht_h_grid_mean)

! =============================================================================
! NEXT PERIOD =================================================================
! =============================================================================

    ! currently assumed capital in the next period
    ! vector because depreciation shock in disaster state 
    k_nxt      = nxt_mat_2(:,1)

    ! relative productivity next period
    zf_nxt(:,1)   = 1.0_dp
    zf_nxt(:,2)   = zeta*exp(nxt_mat_2(:,2))
    omg_nxt      = exp(nxt_mat_2(:,3)) + omg_shift  - b_lmbd*bond_supply_shock_nxt

    ! value functions
    do iii = 1,n_I
        v_temp(:,iii)  = nxt_mat(:, iii )
        mc_temp(:,iii) = nxt_mat(:, n_I + iii)
    enddo

    ! tomorrow's price, labor choice, infl and price of time endowment
    s_nxt               = nxt_mat(:,2*n_I + 1)
    q_nxt               = nxt_mat(:,2*n_I + 2) 
    l_aggr_nxt(:,1)     = nxt_mat(:,2*n_I + 3)
    l_aggr_nxt(:,2)     = nxt_mat(:,2*n_I + 4)
    infl_nxt(:,1)       = nxt_mat(:,2*n_I + 5)
    infl_nxt(:,2)       = nxt_mat(:,2*n_I + 6)
    ch_spending_nxt     = nxt_mat(:,2*n_I + 7)*exp(dz_vec)
    cf_spending_nxt     = nxt_mat(:,2*n_I + 8)*exp(dz_vec)


    ! local inflation P'_h/P_h = (P/P_h)/(P'/P'_h) * P'/P and equivalently P'_f/P_f
    
    homegood_infl_nxt(:,1) = infl_nxt(:,1)*( &  
                           ( varsigma_vec(1)  +  (1.0-varsigma_vec(1) )*    s**(sigma-1.0) ) / & 
                           ( varsigma_vec(1)  +  (1.0-varsigma_vec(1) )*s_nxt**(sigma-1.0) )  )**(1.0/(1.0-sigma)) 
    homegood_infl_nxt(:,2) = infl_nxt(:,2)* s_nxt/s *( ( varsigma_vec(2)*    s**(1.0-sigma)+ (1.0-varsigma_vec(2)) ) / & 
                                           ( varsigma_vec(2)*s_nxt**(1.0-sigma)+ (1.0-varsigma_vec(2)) ) )**(1.0/(1.0-sigma)) 

    ! given employed labor, allocate capital to equalize capital returns across countries
        kappa_nxt(:,1) = k_nxt *  ((s_nxt)**(1.0/(1.0-aalpha))) * zf_nxt(:,1) / zf_nxt(:,2) * l_aggr_nxt(:,1)/l_aggr_nxt(:,2) / & 
                           (1.0 + ((s_nxt)**(1.0/(1.0-aalpha))) * zf_nxt(:,1) / zf_nxt(:,2) * l_aggr_nxt(:,1)/l_aggr_nxt(:,2) )

    kappa_nxt(:,2) = k_nxt - kappa_nxt(:,1)

    ! chosen capital, absent depreciation shock 
    next_k         = k_nxt(1)

    investment = next_k - k_aggr*(1.0-ddelta)

    ! update price of capital but only for next iteration
    q_new =  inv_share_h**(1.0/(1.0-sigma)) * (1 + &
              (1.0-inv_share_h)/inv_share_h*(s)**(sigma-1.0))**(1.0/(1.0-sigma)) & 
              * (next_k/state_grid(idx_k,sss))**chiX

    ! output in next period
    y_next(:,1) = exp((log(zf_nxt(:,1)) + dz_vec)*(1.0_dp-aalpha))*l_aggr_nxt(:,1)* & 
                  ( (kappa_nxt(:,1)/l_aggr_nxt(:,1))**aalpha)
    y_next(:,2) = exp((log(zf_nxt(:,2)) + dz_vec)*(1.0_dp-aalpha))* l_aggr_nxt(:,2)* & 
                  ( (kappa_nxt(:,2)/l_aggr_nxt(:,2))**aalpha)
    
   
    w_next_choice(:,1) = (1-aalpha)*y_next(:,1)/l_aggr_nxt(:,1)
    w_next_choice(:,2) = s_nxt**(-1) * (1-aalpha)*y_next(:,2)/l_aggr_nxt(:,2)

    ! profit next period per unit of capital
        pi_nxt = aalpha*exp( (log(zf_nxt(:,1)) + dz_vec)*(1.0_dp-aalpha)) &
               *(l_aggr_nxt(:,1)**(1.0_dp-aalpha))*kappa_nxt(:,1)**(aalpha-1.0)

    ! capital return, adjustment in disaster state
    rk_vec = ( (1.0_dp-ddelta)*q_nxt + pi_nxt)/q_current
    rk_vec(n_quad) =exp(dz_vec(n_quad))*rk_vec(n_quad) 

    ! cost of unit of bundle in terms of home good (P/P_h and P*/P^*_h 
    ! (one bundle unit valued at P, divided by cost of home good -> units of home good per bundle
    P_div_P_h_nxt(:,1) = (varsigma_vec(1)   & 
                       + (1.0 -varsigma_vec(1))*s_nxt**(sigma-1))**(1.0/(1.0-sigma))
    P_div_P_h_nxt(:,2) = s_nxt**(-1) * (varsigma_vec(2)*s_nxt**(1.0-sigma) + (1.0 -varsigma_vec(2)))**(1.0/(1.0-sigma))

    ! calculate risk-free returns, given nominal rate guess 
    ! transform everything in units of home good
    rf_vec(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1) 
    rf_vec(:,2) = (nom_i_vec(1) + nom_i_vec(2))/homegood_infl_nxt(:,2) 
   
    ! seignorage
    if (bg_yss > 0) then 
        seigniorage_nxt = cf_spending_nxt/(ch_spending_nxt + cf_spending_nxt)*(bg_yss + bond_supply_shock_nxt)*(ch_spending_nxt + cf_spending_nxt)*omg_nxt
        seignorage_transfer(:,1) =  seigniorage_nxt 
        seignorage_transfer(:,2) = -seigniorage_nxt
    else
        seignorage_transfer = 0.0_dp
    endif

! =============================================================================
! SOLVE FOR NOMINAL RATE IN EQUILIBRIUM =======================================
! =============================================================================

    ! check whether it generates excess demand or not
    c_spend_temp    = c_spending_vec
    c_temp          = c_spend_temp/P_div_P_h
    
    bF_share_temp = bF_share_vec
    share_temp      = share_vec
    
    rf_vec_temp(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1) 
    rf_vec_temp(:,2) = (nom_i_vec(1) + nom_i_vec(2))/homegood_infl_nxt(:,2) 

    call calc_excess_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, &
                              bbeta_adj,  P_div_P_h, P_div_P_h_nxt,   &  
                              q_current, w_choice_vec, l_aggr_vec, sss, big_weight_vec, &  
                              c_temp, c_spend_temp, share_temp, bF_share_temp, k_temp_vec, b_temp_vec, excess_b, constraint_binding_vec, mc_temp)

    excess_b_out = excess_b 

    rf_vec(:,1)  = nom_i_vec(1)/homegood_infl_nxt(:,1) 
    rf_vec(:,2)  = (nom_i_vec(1) + nom_i_vec(2))/homegood_infl_nxt(:,2)

! =============================================================================
! GIVEN CURRENT OVERALL BOND HOLDINGS, UPDATE NOMINAL RATE SPREAD
! =============================================================================
    if (foreign_trading == 1) then 
    
    bF_share_temp = bF_share_vec 
    
    rf_vec_temp(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1) 
    rf_vec_temp(:,2) = (nom_i_vec(1) + nom_i_vec(2))/homegood_infl_nxt(:,2) 


    call calc_excess_foreign_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, bbeta_adj,  P_div_P_h, P_div_P_h_nxt, &  
                                       w_choice_vec, l_aggr_vec, sss, big_weight_vec, mc_temp, &  
                                       c_temp, c_spend_temp, share_vec, bF_share_temp, bh_temp_vec, bf_temp_vec, & 
                                       excess_b, safe_constraint_binding_vec)

                            
    nom_i_guess = nom_i_vec(2)
    
    ! store
    nom_i_A    = nom_i_guess
    excess_b_A = excess_b

      ! find bound between solver finds solution for nominal rate
      if (excess_b == 0) then

        excess_b_S = excess_b

      else
        if (excess_b > 0.0_dp) then

          iter = 0
          do while (excess_b > 0.0_dp .and. iter <= 5000)
              iter = iter + 1

              ! initial guess - last period nominal rate
              nom_i_guess = nom_i_guess - 1E-03_dp

              ! check whether it generates excess demand or not
              rf_vec_temp(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1)
              rf_vec_temp(:,2) = (nom_i_guess + nom_i_vec(1))/homegood_infl_nxt(:,2) 

              bF_share_temp = bF_share_vec
              call calc_excess_foreign_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, & 
                                        bbeta_adj,  P_div_P_h, P_div_P_h_nxt, &  
                                          w_choice_vec, l_aggr_vec, sss, big_weight_vec, mc_temp, &  
                                          c_temp, c_spend_temp, share_vec, bF_share_temp, bh_temp_vec, bf_temp_vec,&
                                          excess_b, safe_constraint_binding_vec)

          enddo
        
          if (iter >= 999) then
            write(*,*) 'No convergence in finding nom_i_low f'
            stop
          end if

          nom_i_B = nom_i_guess
          excess_b_B = excess_b

        else

          iter = 0
          do while (excess_b < 0.0_dp .and. iter <= 5000)
              
              iter = iter + 1

                ! initial guess - last period nominal rate
                nom_i_guess  = nom_i_guess + 1E-03_dp

                rf_vec_temp(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1) 
                rf_vec_temp(:,2) = (nom_i_guess + nom_i_vec(1))/homegood_infl_nxt(:,2) 

                ! check whether it generates excess demand or not
                bF_share_temp = bF_share_vec
                call calc_excess_foreign_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, bbeta_adj, P_div_P_h, P_div_P_h_nxt,&  
                                          w_choice_vec, l_aggr_vec, sss, big_weight_vec, mc_temp, &  
                                          c_temp, c_spend_temp, share_vec, bF_share_temp, bh_temp_vec, bf_temp_vec, & 
                                          excess_b, safe_constraint_binding_vec)
          enddo

          if (iter >= 999) then
              write(*,*) 'No convergence in finding nom_i_high foreign'
              write(*,*) sss, excess_b
              write(*,*) share_vec 
              write(*,*) bF_share_temp 
              stop
          end if
          nom_i_B = nom_i_guess
          excess_b_B = excess_b

        endif

        ! check that root is bracketed
        if ((excess_b_A*excess_b_B) > 0) then
           write(*,*) 'ERROR: Initial bracket does not contain root. B'
           stop
        endif

        ! swap a and b
        if ( abs(excess_b_A) < abs(excess_b_B) )then
           temp = nom_i_A
           nom_i_A = nom_i_B
           nom_i_B = temp

           temp = excess_b_A
           excess_b_A = excess_b_B
           excess_b_B = temp

        endif

        nom_i_C    = nom_i_A
        excess_b_C = excess_b_A
        excess_b_S = excess_b_B
        mflag = 1 ! set mflag

        iiter = 0 
        do while ( excess_b_S /= 0 .and. abs(nom_i_A - nom_i_B) > 1E-15_dp)
            
            iiter = iiter + 1

            if ( (excess_b_A /= excess_b_C) .and. (excess_b_C /= excess_b_B) ) then ! inverse quadratic interpolation

            nom_i_S = nom_i_A * excess_b_B * excess_b_C / ( (excess_b_A - excess_b_B) * (excess_b_A - excess_b_C) ) + &
                      nom_i_B * excess_b_A * excess_b_C / ( (excess_b_B - excess_b_A) * (excess_b_B - excess_b_C) ) + &
                      nom_i_C * excess_b_A * excess_b_B / ( (excess_b_C - excess_b_A) * (excess_b_C - excess_b_B) )

          else ! secant method

            nom_i_S = nom_i_B - excess_b_B * (nom_i_B - nom_i_A) /   (excess_b_B - excess_b_A)

          endif

          if ( ( (  ( nom_i_S > ( 3*nom_i_A + nom_i_B )/4 .and. nom_i_S < nom_i_B) .or. &
                    ( nom_i_S < ( 3*nom_i_A + nom_i_B )/4 .and. nom_i_S > nom_i_B)  ) == .FALSE. ) .or. &
               (mflag == 1 .and. abs(nom_i_S - nom_i_B) >= abs(nom_i_B - nom_i_C)/2  )             .or. &
               (mflag == 0 .and. abs(nom_i_S - nom_i_B) >= abs(nom_i_B - nom_i_D)/2  )             .or. &
               (mflag == 1 .and. abs(nom_i_B - nom_i_C) <  abs(brent_delta)  )                     .or. &
               (mflag == 0 .and. abs(nom_i_B - nom_i_D) <  abs(brent_delta)  )                          &
             ) then

            nom_i_S = (nom_i_A + nom_i_B )/ 2
            mflag = 1
          else
            mflag = 0
          endif

              rf_vec_temp(:,1) = nom_i_vec(1)/homegood_infl_nxt(:,1) 
              rf_vec_temp(:,2) = (nom_i_S + nom_i_vec(1))/homegood_infl_nxt(:,2) 

            ! check whether it generates excess demand or not
            bF_share_temp = bF_share_vec
            call calc_excess_foreign_bond_nom(aggr_wealth, omg, v_temp, rf_vec_temp, rk_vec, seignorage_transfer, & 
                                              bbeta_adj, P_div_P_h, P_div_P_h_nxt,   &  
                                              w_choice_vec, l_aggr_vec, sss, big_weight_vec, mc_temp,&  
                                              c_temp, c_spend_temp, share_vec, bF_share_temp, bh_temp_vec, bf_temp_vec,& 
                                              excess_b_S, safe_constraint_binding_vec)

          nom_i_D = nom_i_C
          nom_i_C = nom_i_B
          excess_b_C = excess_b_B

          if ((excess_b_A*excess_b_S) < 0) then
            nom_i_B = nom_i_S
            excess_b_B = excess_b_S
          else
            nom_i_A = nom_i_S
            excess_b_A = excess_b_S
          endif

          ! swap a and b
          if ( abs(excess_b_A) < abs(excess_b_B) )then
             temp    = nom_i_A
             nom_i_A = nom_i_B
             nom_i_B = temp

             temp       = excess_b_A
             excess_b_A = excess_b_B
             excess_b_B = temp

          endif

            if (iiter > 2000) then 
                write(*,*) 'no convergence nom_i f'
                stop
            endif


        enddo
        
        nom_i_guess = nom_i_B
      endif

      if (abs(excess_b_B)>1E-02_dp) then
          write(*,*) 'WARNING: EXCESS AGGREGATE BONDS HOLDINGS (FOREIGN).'
          write(*,*) sss
          write(*,*) bF_share_temp
          write(*,*)
          write(*,*) rf_vec_temp(:,1) - rf_vec_temp(:,2)
          write(*,*)
          write(*,*) rf_vec_temp(:,1) - rf_vec_temp(:,2)
          ! stop
      endif

    nom_i_vec(2) = nom_i_B

endif

! =============================================================================
! STORE VALUE FUNCTION, NEXT CAPITAL, NEXT WEALTH SHARES ======================
! =============================================================================

    c_spend_temp = c_spending_vec
    do iii = 1, n_I

            ! read in preferece parameters for current agent
            IES   = IES_vec(iii)
            bbeta = bbeta_vec(iii) + bbeta_adj(iii)
            gma   = gma_vec(iii)
            chi0   = chi0_vec(iii)
            v_normalization = v_normalization_vec(iii)
          
            c_cost     = P_div_P_h(iii)
            c_cost_nxt = P_div_P_h_nxt(:,iii)
          
            wealth     = aggr_wealth*wealth_share_grid(iii,sss) 

            labor      = l_aggr_vec(iii)

            c_spend =  minval( [maxval([c_spend_temp(iii), min_cons_sav]), & 
                       wealth + w_choice_vec(iii)*labor-min_cons_sav] )
          
            consumption = c_spend/c_cost
            c_vec(iii) = consumption
            
            r_alpha_vec = (share_vec(iii)-bF_share_vec(iii))*rf_vec(:,1)/(1.0_dp - omg) & 
                                           + bF_share_vec(iii)*rf_vec(:,2) & 
                          + (1.0_dp - share_vec(iii))*rk_vec

            labor_part =  1.0_dp + (1.0_dp/IES - 1.0_dp) * chi0*chi/(1.0_dp + chi) * labor**( (1.0_dp+chi)/chi )

            savings     = wealth + w_choice_vec(iii)*labor - consumption*c_cost
            
            ! calculate next period's wealth shares
            next_period_share = ( savings*r_alpha_vec + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer(:,iii))/(tot_wealth_ss_vec(iii))

            v_vec_twisted = (v_temp(:,iii)*next_period_share)**(1.0-gma)
            EV = sum(big_weight_vec * v_vec_twisted) **(1.0 / (1.0 - gma))
            ! update wealth transitions accordingly
            savings_vec(iii) = savings
            nxt_wealth_vec(:,iii) = (savings*r_alpha_vec + seignorage_transfer(:,iii))
            
            objf = (v_normalization * labor_part**(1.0_dp/ies)*(consumption)**((ies - 1.0_dp)/ies) + bbeta *EV**((IES - 1.0_dp)/IES))**(ies / (ies - 1.0_dp))

            v_new(iii) = objf*( (wealth + q_l_ss_vec(iii))/tot_wealth_ss_vec(iii))**(-1)

    enddo

    theta_nxt = nxt_wealth_vec(:,1)/(sum(nxt_wealth_vec,2) )

    ! given new capital choices: next period capital
    k_next_new = sum(savings_vec*(1.0-share_vec)/q_current)

! =============================================================================
! UPDATE RELATIVE COST S ======================================================
! =============================================================================
      
      investment = next_k - k_aggr*(1.0-ddelta)
      
      ! check these derivations 
      s_new = s
      diff = 1.0_dp
      iter_inner = 0 
      do while (diff > conv_crit) 
          iter_inner = iter_inner + 1
          
          !q_temp =  inv_share_h**(1.0/(1.0-sigma)) * (1 + &
          !    (1.0-inv_share_h)/inv_share_h*(s_new*exp(-1*current_eps))**(sigma-1.0))**(1.0/(1.0-sigma)) & 
          !    * (next_k/state_grid(idx_k,sss))**chiX
          inv_h = investment*q_current/(1.0+((s_new)**(sigma-1.0))*(1.0-inv_share_h)/(inv_share_h))
          inv_f = inv_h * (s_new)**(sigma)*(1.0-inv_share_h)/inv_share_h 

          ! consumption in home units, cancels out below
          ch_vec(1) = c_spending_vec(1)/ ( (1.0-varsigma_vec(1) ) & 
              /(varsigma_vec(1) ) * s_new**(sigma-1.0_dp) + 1.0_dp) 
          ch_vec(2) = c_spending_vec(2)/ ( (1.0-varsigma_vec(2))/varsigma_vec(2) * s_new**(sigma-1.0_dp) + 1.0_dp) 
          s_update = ( sum((c_spending_vec - ch_vec))/(sum( ch_vec)) * (y_current(1)-inv_h)/( (y_current(2)-inv_f)) )**(-1)
          
          diff = abs(s_update-s_new)
          if (isnan(s_new)) then 
              write(*,*) 'NaN s'
              write(*,*) y_current - [inv_h, inv_f], (sum(ch_vec))
              write(*,*) sss, s, s_new, iter_inner, investment
              write(*,*) y_current, inv_h, inv_f 
              stop
          endif
          
          s_new = s_new + 0.1*(s_update-s_new)
         
          if (iter_inner > 10000) then 
             write(*,*) 'no convergence s', diff
             stop
         endif
      enddo
        
      ! write(*,*) s_new

! =============================================================================
! UPDATE CONSUMPTION CHOICE ===================================================
! =============================================================================

    c_spend_temp = c_spending_vec
    do iii = 1, n_I

            ! read in preferece parameters for current agent
            IES   = IES_vec(iii)
            bbeta = bbeta_vec(iii) + bbeta_adj(iii)
                gma   = gma_vec(iii)
            
            chi0   = chi0_vec(iii)
            v_normalization = v_normalization_vec(iii)
          
            c_cost     = P_div_P_h(iii)
            c_cost_nxt = P_div_P_h_nxt(:,iii)
          
            wealth     = aggr_wealth*wealth_share_grid(iii,sss) 

                labor = l_aggr_vec(iii)
                labor_nxt = l_aggr_nxt(:, iii)
           
            c_spend =  minval( [maxval([c_spend_temp(iii), min_cons_sav]), & 
                       wealth + w_choice_vec(iii)*labor-min_cons_sav] )
          
            consumption = c_spend/c_cost
           
            ! returns for resource accumulation
            r_alpha_vec = (share_vec(iii)-bF_share_vec(iii))*rf_vec(:,1)/(1.0_dp - omg) & 
                                           + bF_share_vec(iii)*rf_vec(:,2) & 
                          + (1.0_dp - share_vec(iii))*rk_vec

            labor_part =  1.0_dp + (1.0_dp/IES - 1.0_dp) * chi0*chi/(1.0_dp + chi) * labor**( (1.0_dp+chi)/chi )
            
            diff_inner = 1.0_dp
            iter_inner = 0

            do while (diff_inner>conv_crit)
             
                iter_inner  = iter_inner+1
                savings     = wealth + w_choice_vec(iii)*labor - consumption*c_cost
             
                ! calculate next period's wealth shares
                next_period_share = ( savings*r_alpha_vec + exp(dz_vec)*q_l_ss_vec(iii) + seignorage_transfer(:,iii))/tot_wealth_ss_vec(iii)
                
                ! returns for FOC (always include omg)
                r_alpha_omg_vec = (share_vec(iii)-bF_share_vec(iii))*rf_vec(:,1)/(1.0_dp - omg) & 
                                           + bF_share_vec(iii)*rf_vec(:,2) & 
                          + (1.0_dp - share_vec(iii))*rk_vec
                
                util_c_deriv = (labor_part / consumption)**(1/IES)
                mc_new(iii) = util_c_deriv*((q_l_ss_vec(iii) + wealth)/tot_wealth_ss_vec(iii))**(1/IES)
                v_vec_twisted = (v_temp(:,iii)*next_period_share)**(1.0_dp-gma)
             
                EV = sum(big_weight_vec*v_vec_twisted)**( 1.0_dp / (1.0_dp - gma))
                M_vec_mat(:, iii) = bbeta * mc_temp(:, iii)/util_c_deriv * ((v_temp(:, iii)/EV) **(1.0_dp/IES - gma))  * (next_period_share **(-gma))
       
                    temp = sum(M_vec_mat(:, iii)*r_alpha_omg_vec*c_cost/c_cost_nxt*big_weight_vec)

                cons_update = consumption / (temp)**IES 

                cons_update = minval([maxval([cons_update, min_cons_sav]), (wealth + w_choice_vec(iii)*labor-min_cons_sav)/c_cost])
                
                diff_inner = abs(consumption-cons_update)
                 
                consumption = consumption + 0.5*(cons_update-consumption)
                
                if (phi_w == 0) then
                    labor_new = ( w_choice_vec(iii) / P_div_P_h(iii)    *   &
                                (1.0_dp + (1.0_dp/IES_vec(iii) - 1.0_dp) * chi0_vec(iii)*chi/(1.0_dp + chi) * &
                                labor**( (1.0_dp+chi)/chi ) ) /    &
                                ( 1.0_dp/IES_vec(iii) * consumption * chi0_vec(iii) )  )**(chi) 
                    labor = labor + 0.2*(labor_new-labor)
                endif

                if (iter_inner > 1000) then
                    if (diff_inner > 0.1_dp) then
                        write(*,*) 'NO CONVERGENCE 1.'
                        stop
                    else
                        write(*,*) 'WARNING: NO CONVERGENCE 1.'
                        diff_inner = 0.0_dp
                    endif
                endif
            enddo
        
        c_spending_vec_new(iii) = consumption*c_cost
        ! c_vec(iii) = consumption 
        savings_vec(iii) = savings
        
        if (phi_w == 0) then 
            l_temp_vec(iii) = labor
        endif

    enddo
      
! =============================================================================
! UPDATE LABOR, WAGES, CAPITAL ALLOCATION =====================================
! =============================================================================
    
    if (phi_w > 0.0_dp) then
    l_aggr_temp   = l_aggr_vec 
    w_choice_temp = w_choice_vec

    ! update wages in both countries
    do iii = 1,2

          w_temp = w_choice_temp(iii)
          w_diff = 1.0_dp
          iter = 0
          do while ( (w_diff > conv_crit .and. iter<100) .or. iter < 10 )

            iter = iter + 1

            w_current_next = w_temp 
            w_current_next(n_quad) = w_temp*exp(dz_vec(n_quad))
           
                    gma   = gma_vec(iii)
              w_temp_new = ( &
                            vareps_w/(1.0_dp-vareps_w)   * P_div_P_h(iii)   &
                            * ( -  1.0_dp/IES_vec(iii) * c_vec(iii) * chi0_vec(iii) * (l_aggr_temp(iii)**(1.0/chi)) / &
                               (1.0_dp + (1.0_dp/IES_vec(iii) - 1.0_dp) * chi0_vec(iii)*chi/(1.0_dp + chi) * &
                           l_aggr_temp(iii)**( (1.0_dp+chi)/chi ) ) )    &
                        + w_temp*phi_w/(1.0_dp-vareps_w) * ( ( &
                                w_temp/w_current(iii)*infl_vec(iii)/P_div_P_h(iii) - 1.0)   &
                               *w_temp/w_current(iii)*infl_vec(iii)/P_div_P_h(iii)  &
                            - sum( M_vec_mat(:,iii) * big_weight_vec * & 
                             homegood_infl_nxt(:,iii)/infl_nxt(:,iii)* &
                              ( w_next_choice(:,iii)/w_current_next*homegood_infl_nxt(:,iii) - 1.0) * &
                             (( w_next_choice(:,iii)/w_current_next)**2) * homegood_infl_nxt(:,iii)* l_aggr_nxt(:,iii)/l_aggr_temp(iii) ) &
                             ) & 
                             )
    

           w_diff = abs(w_temp_new-w_temp)
           w_temp = w_temp + 0.005_dp*(w_temp_new-w_temp)


          end do

          if (isnan(w_temp)) then
            write(*,*) 'error w isnan', M_vec_mat
            stop
          elseif (w_temp<=0) then
            write(*,*) 'error w<0'
            stop
          endif

          l_aggr_temp(iii) = kappa_vec(iii)*(( ( zf_vec(iii)**(1.0_dp-aalpha) )*(1.0_dp-aalpha)/w_temp/s_vec(iii))**(1.0_dp/aalpha))
          w_choice_temp(iii) =  w_temp

    enddo

    l_aggr_new = l_aggr_temp
    else
        
        l_aggr_new(1) = l_temp_vec(1)
        l_aggr_new(2) = l_temp_vec(2)

    endif

    w_choice_new = w_choice_vec

! =============================================================================
! UPDATE INFLATION RATE =======================================================
! =============================================================================

    ! taylor rule: exp(i) = (rho*exp(i_last) + (1-rho)* exp(k)exp(pi)^1.5)exp(m_shock) ->  exp(pi)=(exp(i)exp(-m_shock) - rho*exp(i_last))*exp(-k)/(1-rho_i))^(1/1.5)
    ! rf fix : i up : E[pi] up, p down, current inflation down
         
            ih_last = state_grid(idx_ih,sss)
            if_last = state_grid(idx_if,sss)
    
        infl_new(1) = ((nom_i_vec_in(1)/((1.0+ih_last)**rho_i))**(1.0/(1.0-rho_i)) * exp(-tayl_ic_h)/(y_current(1)**phi_yh))**(1.0_dp/phi_h)
            
        infl_new(2) = (((nom_i_vec_in(1)+nom_i_vec_in(2))/((1.0+if_last)**rho_i))**(1.0/(1.0-rho_i)) *  exp(-tayl_ic_f)/(y_current(2)**phi_yf))**(1.0_dp/phi_f)
        share_vec_new = share_temp

    if (foreign_trading == 1 .and. outer_iter > 20) then     
        bF_share_vec  = bF_share_temp 
    else
        bF_share_vec = 0.0
        nom_i_vec(2) = maxval( [ 1*sum(big_weight_vec*M_vec_mat(:,1)*(1.0_dp/homegood_infl_nxt(:,2))* P_div_P_h(1)/P_div_P_h_nxt(:,1)), & 
                                 1*sum(big_weight_vec*M_vec_mat(:,2)*(1.0_dp/homegood_infl_nxt(:,2))* P_div_P_h(2)/P_div_P_h_nxt(:,2))  ] )**(-1) - nom_i_vec(1)

    endif


end subroutine calc_equilibrium_ifixed


subroutine update_ifixed(infl_h, ih_fixed, next_state, n_nxt, sss, smol_coeffs, k_next, next_zf, next_omg, nom_i_vec, s, q_current, & 
                         c_spending_vec, share_vec, bF_share_vec, infl_vec, l_aggr_vec, inv_h, inv_f, bond_supply_shock, bond_supply_shock_nxt, &
                         excess_b, v_new, mc_new, w_choice_new)

    use mod_smolyak, only : Smolyak_Polynomial2
    use mod_param, only: n_active_dims, n_quad, n_states, max_smol_level, smol_elem_ani , n_I, vector_mus_dimensions, dz_vec, &
    dz_vec_adj, k_grid_mean, k_grid_dev, varsigma_vec, &
    wh_grid_mean, wh_grid_dev, wf_grid_mean, wf_grid_dev, ih_grid_mean, ih_grid_dev, if_grid_mean, if_grid_dev, & 
    tht_h_grid_mean, tht_h_grid_dev, sigma, n_states
    
    integer, intent(in) :: n_nxt, sss
   real(dp), intent(in) :: infl_h, ih_fixed, bond_supply_shock, bond_supply_shock_nxt, next_zf(n_quad), next_omg(n_quad)
   real(dp), intent(inout) :: next_state(n_quad, n_active_dims), smol_coeffs(n_states,n_nxt), k_next(n_quad), &
                               nom_i_vec(2), s, q_current, & 
                              c_spending_vec(n_I), share_vec(n_I), bF_share_vec(n_I), infl_vec(2), l_aggr_vec(2), inv_h, inv_f
   real(dp), intent(out) :: mc_new(n_I), v_new(n_I), excess_b, w_choice_new(2)
   integer :: counter , outer_iter
   real(dp) :: share_vec_new(n_I), polyn_points(n_quad, n_states), nxt_mat(n_quad, n_nxt), nxt_mat_2(n_quad, 3), q_new, s_new, k_next_new, c_spending_vec_new(n_I), & 
               next_state_new(n_quad, n_active_dims), infl_new(2), l_aggr_new(2),  theta_nxt(n_quad), M_vec(n_quad, n_I), diff
    next_state_new = next_state 
    diff = 1.0
    outer_iter =0 
    do while (diff > conv_crit .or. outer_iter < 20)
            outer_iter = 1 + outer_iter
            ! interpolate future values 
            ! use MKL Library for multiplication, same as
            ! nxt_mat        = matmul(polyn_points, smol_coeffs)
            polyn_points = Smolyak_Polynomial2(next_state, n_active_dims, &
                                      n_quad, n_states, max_smol_level, smol_elem_ani)

            CALL DGEMM('N','N', n_quad, n_nxt, n_states, 1.0_dp, polyn_points, n_quad, & 
                               smol_coeffs, n_states, 0.0_dp, nxt_mat, n_quad)   

             ! non interpolated values
             nxt_mat_2(:, 1) = k_next
             nxt_mat_2(:, 2) = next_zf
             nxt_mat_2(:, 3) = next_omg
            
             nom_i_vec(1) = ih_fixed
             infl_vec(1) = infl_h
                
            ! calculate current state equilibrium 
            call calc_equilibrium_ifixed( nxt_mat, n_nxt, nxt_mat_2, sss, outer_iter, &
                        c_spending_vec, s, l_aggr_vec, q_current, infl_vec, &
                        nom_i_vec, share_vec, share_vec_new, bF_share_vec, inv_h, inv_f, & 
                        bond_supply_shock, bond_supply_shock_nxt, &
                        q_new, s_new, k_next_new, c_spending_vec_new, & 
                        infl_new, v_new, mc_new, l_aggr_new, w_choice_new, theta_nxt, M_vec, excess_b)


                counter = 0
                
                if (vector_mus_dimensions(1) > 0) then
                    counter = counter + 1
                    next_state_new(:, counter)   = (k_next_new/exp(dz_vec_adj) -k_grid_mean )/k_grid_dev
                endif
    
                if (vector_mus_dimensions(2) > 0) then
                    counter = counter + 1
                    next_state_new(:, counter) = (theta_nxt - tht_h_grid_mean)/tht_h_grid_dev
                endif
    
    
                if (vector_mus_dimensions(3) > 0) then
                counter = counter + 1
                endif
           
    
                if (vector_mus_dimensions(4) > 0) then
                    counter = counter + 1
                    next_state_new(:, counter)  = (w_choice_new(1)/exp(dz_vec_adj)/ & 
                                                                  (varsigma_vec(1)  + & 
                                                                  (1.0-varsigma_vec(1))*s**(sigma-1.0))**(1.0/(sigma-1.0)) & 
                                                                  -wh_grid_mean) /wh_grid_dev
                endif
    
                if (vector_mus_dimensions(5) > 0) then
                    counter = counter + 1
                    next_state_new(:, counter)  = (w_choice_new(2)/exp(dz_vec_adj)*s/ & 
                                                                  (varsigma_vec(2)*s**(1.0_dp-sigma) + (1.0-varsigma_vec(2)))**(1.0/(sigma-1.0)) & 
                                                                  -wf_grid_mean) /wf_grid_dev
                endif
    
                if (vector_mus_dimensions(6) > 0) then
                counter = counter + 1
                endif
    
                if (vector_mus_dimensions(7) > 0) then
                    counter = counter + 1
                    next_state_new(:, counter)  = ((nom_i_vec(1)-1.0) - ih_grid_mean) /ih_grid_dev 
                endif
    
                if (vector_mus_dimensions(8) > 0) then
                    counter = counter + 1
                    next_state_new(:, counter)  = ((nom_i_vec(1) + nom_i_vec(2)-1.0) - if_grid_mean) /if_grid_dev 
                endif

                if (vector_mus_dimensions(9) > 0) then
                counter = counter + 1
                endif
    

        ! constrain transition matrix
        where(next_state_new > 1.0_dp) 
            next_state_new   = 1.0_dp 
        endwhere
        where(next_state_new < -1.0_dp)
            next_state_new   = -1.0_dp
        endwhere
        

        diff =  maxval([ maxval(abs((c_spending_vec_new-c_spending_vec))), &
                         maxval(abs(log(k_next_new/k_next))), &
                         maxval(abs(log(l_aggr_new/l_aggr_vec))), &
                         (abs(log(q_new/q_current))), &
                         (abs(log(infl_new(2)/infl_vec(2)))), &
                         maxval(abs(next_state_new - next_state)) ])

        ! update value function
        next_state = next_state +  0.2_dp*(next_state_new-next_state)
         
        if (abs((q_new-q_current)) > 0.01) then
           q_current = q_current + 0.002_dp*sign(1.0_dp, q_new - q_current)
        else
           q_current = q_current + 0.1_dp*(q_new-q_current)
       endif
    
        if (abs((s-s_new)) > 0.05) then
           s = s + 0.025_dp*sign(1.0_dp, s_new - s)
        else
           s = s + 0.5_dp*(s_new-s)
       endif
    
        where (abs(l_aggr_vec-l_aggr_new) > 0.01)
           l_aggr_vec = l_aggr_vec + 0.002_dp*sign(1.0_dp, l_aggr_new - l_aggr_vec)
        elsewhere
           l_aggr_vec = l_aggr_vec + 0.2_dp*(l_aggr_new-l_aggr_vec)
        endwhere
    
        c_spending_vec = c_spending_vec  + 0.50*(c_spending_vec_new - c_spending_vec)
        infl_vec(2)    = infl_vec(2)     + 0.20*(infl_new(2) - infl_vec(2))
        k_next         = k_next          + 0.20*(k_next_new - k_next)

        enddo
        share_vec      = share_vec_new 

end subroutine 


end module mod_calc

