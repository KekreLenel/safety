! -------------------------------------------------------------------------
! mod_results.f90: results module, public subroutines:
! - create_results: given model solution, creates and stores simulations and IRFs 
! -------------------------------------------------------------------------
! authors:         Rohan Kekre and Moritz Lenel
! for updates see: https://github.com/KekreLenel/safety
! -------------------------------------------------------------------------
module mod_results

use omp_lib
use base_lib, only: dp, sqrt_eps

implicit none
private

public :: create_results

contains

subroutine create_results()

use mod_smolyak, only : Smolyak_Polynomial
use mod_param, only: n_states, idx_k, idx_thh, idx_zf, n_shocks, & 
                     idx_zf, idx_wh, idx_wf, idx_dis, idx_omg, shock_grid, & 
                     state_grid, n_quad, aalpha, smolyak_d, vector_mus_dimensions, &
                     next_zf_mat, smol_polynom, results_folder, &  
                     max_smol_level, smol_elem_ani, k_grid_mean, k_grid_dev, n_uni_quad_zf, &
                     tht_h_grid_mean, tht_h_grid_dev, wh_grid_mean, wh_grid_dev, & 
                     wf_grid_mean, wf_grid_dev, zf_grid_mean, zf_grid_dev, n_interp, &
                     sidx_zf, wealth_share_grid, smol_grid, &
                     varsigma_vec, sigma, n_active_dims, n_I, no_shock_idx, &
                     dis_grid_dev, dis_grid_mean, quad_weight_vec, &
                     idx_ih, idx_if, ih_grid_mean, ih_grid_dev, if_grid_mean, if_grid_dev, irf_indices, &
                     sidx_z, sidx_dis, sidx_p, sidx_omg, irf_shock_sizes, irf_ll_length, &
                     irf_indices_ll, sidx_bg, omg_grid_mean, omg_grid_dev, run_bg, run_samp, & 
                     n_uni_quad_p, n_uni_quad_omg,  n_uni_quad_z, n_sample, n_sims, n_irf_periods, n_sim_periods, n_burn, n_sample_periods


real(dp) :: v_mat(n_states, n_I), k_next_mat(n_quad, n_states),  shock_mat_tmp(n_shocks, n_sim_periods + n_burn), & 
            tht_next_mat(n_quad, n_states), q_mat(n_states),  sample_mat_tmp(n_shocks, n_sample_periods), &
            kappa_mat(n_states,n_I), l_aggr_mat(2, n_states), &
            next_state_mat(n_quad, n_active_dims, n_states), & 
            next_state_irf_mat(n_quad, n_active_dims, n_states), & 
            next_state_irftransition_mat(1, n_active_dims, n_states), & 
            s_mat(n_states)

real(dp) :: start_vec(n_active_dims), diff
integer  :: qqq, n_periods, ttt, aaa
real(dp) :: interp_input_mat(n_states, n_interp), interp_input_irf_mat(n_states, n_interp), & 
            state_variables(smolyak_d), other_variables(n_interp), valuation_input_mat(14,n_states), & 
            interp_input_irf_f_mat(n_states, n_interp), valuation_irf_mat(14,n_states)

character :: command_input*2,  filename*100, n_shock_str*20
integer :: narg, command_len, sss, fff, sample_idx, sample_counter, ppp, ddd, counter, iii, n_active_shocks

INTEGER, parameter :: GENID = 3, SUBID = 1, LSEED = 1
INTEGER            :: SEED(LSEED), STATE(633), LSTATE, IFAIL
integer, parameter :: n_spread = 20
real(dp) :: state_spread_vec(n_spread) 
real(dp), allocatable :: state_series(:,:,:), other_vars_series(:,:,:), shock_series(:,:,:), simshock_series(:,:), &
                              state_spread_series(:,:,:), value_series(:,:,:)

integer  ::  nrhs, nrhs2, nrhs3, nrhs4, lda, ldaf, ldb, ldx
integer  ::  ipiv(n_states), iwork(n_states)
integer  ::  info

real(dp) ::  a(n_states,n_states), af(n_states,n_states), r(n_states), & 
             c(n_states),  work(4*n_states), polyn_points(1,n_states)

real(dp), allocatable :: b(:,:), interp_coeffs(:,:), interp_irf_f_coeffs(:,:), interp_irf_coeffs(:,:), & 
                         state_coeffs(:,:), smol_coeffs(:,:,:), smol_irf_coeffs(:,:,:), value_coeffs(:,:), & 
                         smol_irf_f_coeffs(:,:,:), smol_irftransition_f_coeffs(:,:,:), value_irf_coeffs(:,:), & 
                         smol_transition_eta_h_coeffs(:,:,:), smol_transition_eta_f_coeffs(:,:,:), & 
                         smol_irftransition_coeffs(:,:,:), ferr(:), berr(:)
real(dp), allocatable :: b2(:,:), ferr2(:), berr2(:), b3(:,:), ferr3(:), berr3(:), b4(:,:), ferr4(:), berr4(:)

real(dp) ::  rcond
character(1) :: equed
character :: shock_char*1, shock_char2*2, iter_char*4
 

real(dp) :: current_state(1,n_active_dims), next_state(1,n_active_dims), stochsteady_state(1,n_active_dims), checker , &
            p_dis, big_weight_vec(n_quad), next_state_lo(1,n_active_dims), next_state_hi(1,n_active_dims), state_tmp(smolyak_d,1)
integer  :: q_nxt, sample_length
integer :: lwork
real(dp) :: next_state_tmp(n_quad-3,n_active_dims) 
integer, allocatable :: starting_locs(:), sample_vector(:)
real(dp), allocatable :: starting_states(:,:), state_collection(:,:)
real(dp) , allocatable :: rq(:), shock_grid_temp(:,:), shock_tmp(:), placeholder(:,:), sample_mat(:,:)
real(dp), allocatable :: state_series_N(:,:,:), other_vars_series_N(:,:,:), shock_series_N(:,:,:), value_series_N(:,:,:)
real(dp), allocatable :: state_series_N_helper(:,:,:), other_vars_series_N_helper(:,:,:), value_series_N_helper(:,:,:), shock_series_N_helper(:,:,:), irf_factor2(:), irf_factor3(:), irf_factor4(:)
real(dp), allocatable :: irf_state_series(:,:), irf_other_vars_series(:,:), irf_shock_series(:,:), irf_value_series(:,:)
integer, allocatable ::   iq(:)

        narg=command_argument_count()
        if (narg>0) then
            call get_command_argument(1,command_input,command_len)
        else
            write(*,*) 'ERROR: No run index specified.'
            stop
        endif

    lda  = n_states; ldaf = n_states 
    ldb  = n_states; ldx  = n_states


        ! store grid sizes
        open (unit = 10, file = trim(results_folder) // 'grid.csv', ACTION="write",  &
                & FORM="formatted", ACCESS="sequential")
        write(10,'(7i6)')  n_I, n_states, n_active_dims, n_interp, n_shocks, n_spread, smolyak_d
        close(10)

        ! store grid sizes
        open (unit = 10, file = trim(results_folder) // 'grid_locs.csv', ACTION="write",  &
                & FORM="formatted", ACCESS="sequential")
        write(10,'(14f10.4)')  k_grid_mean, k_grid_dev, tht_h_grid_mean, tht_h_grid_dev, &
                              wh_grid_mean, wh_grid_dev, wf_grid_mean, wf_grid_dev, & 
                              ih_grid_mean, ih_grid_dev, if_grid_mean, if_grid_dev, omg_grid_mean, omg_grid_dev
        close(10)
        write(*,*)
        write(*,*) 'READING RESULTS'

       ! Get all solution files
        open (unit = 10, file = trim(results_folder) // 'next_state_mat.dat', ACTION="read",  &
                & FORM="unformatted", ACCESS="STREAM")
        read(10) next_state_mat; close(10)

        ! print all results
        open (unit = 10, file = trim(results_folder) // 'results_mat.dat', ACTION="read", &
                & FORM="unformatted", ACCESS="STREAM")
        read(10) interp_input_mat; close(10)

        open (unit = 10, file = trim(results_folder) // 'valuation_mat.dat', ACTION="read", &
                & FORM="unformatted", ACCESS="STREAM")
        read(10) valuation_input_mat; close(10)

        write(*,*)
        write(*,*) 'FIND STOCHASTIC STEADY STATE'

        nrhs = n_active_dims
        allocate(b(n_states,nrhs), smol_irftransition_coeffs(1,n_states,nrhs), &
            smol_irftransition_f_coeffs(n_quad,n_states,nrhs), smol_irf_f_coeffs(n_quad,n_states,nrhs), &
            smol_irf_coeffs(n_quad,n_states,nrhs), smol_coeffs(n_quad,n_states,nrhs), ferr(nrhs), &  
             berr(nrhs))
        
        ! polynomial coefficients for transition matrix
        do qqq = 1,n_quad

        b = transpose(next_state_mat(qqq,:,:))

        call F07ABF( 'E','N',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,      &
                     equed, r, c, b, ldb, smol_coeffs(qqq,:,:), & 
                     ldx, rcond, ferr, berr, work, iwork, info)
        
        enddo
        

        ! polynomial coefficients for non-state variables 
        nrhs4 = 14
        nrhs2 = n_interp
        allocate(b2(n_states,nrhs2), interp_coeffs(n_states,nrhs2), & 
                 interp_irf_f_coeffs(n_states,nrhs2), interp_irf_coeffs(n_states,nrhs2), ferr2(nrhs2), &  
                 berr2(nrhs2), value_irf_coeffs(n_states, nrhs4))

        b2 = (interp_input_mat)
         
        call F07ABF( 'E','N',n_states, nrhs2, smol_polynom, lda,af,ldaf,ipiv,      &
                     equed, r, c, b2, ldb, interp_coeffs, & 
                     ldx, rcond, ferr2, berr2, work, iwork, info )
        

        nrhs3 = smolyak_d
        allocate(b3(n_states,nrhs3), state_coeffs(n_states,nrhs3), ferr3(nrhs3), &  
             berr3(nrhs3))
        b3 = transpose(state_grid)
         
        call F07ABF( 'E','N',n_states, nrhs3, smol_polynom, lda,af,ldaf,ipiv,      &
                     equed, r, c, b3, ldb, state_coeffs, & 
                     ldx, rcond, ferr3, berr3, work, iwork, info )

        allocate(b4(n_states,nrhs4), value_coeffs(n_states,nrhs4), ferr4(nrhs4), &  
             berr4(nrhs4))
        b4 = transpose(valuation_input_mat)
         
        call F07ABF( 'E','N',n_states, nrhs4, smol_polynom, lda,af,ldaf,ipiv,      &
                     equed, r, c, b4, ldb, value_coeffs, & 
                     ldx, rcond, ferr4, berr4, work, iwork, info )
        
        ! start from some interior grid point
        current_state(1,:) = smol_grid(1,:)

        diff = 1.0_dp 
        do while (diff > sqrt_eps) 
            
            polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                max_smol_level, smol_elem_ani) 
            
            CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                               smol_coeffs(no_shock_idx,:,:), n_states, 0.0_dp, & 
                               next_state, 1)   
            
            diff = maxval(abs(current_state - next_state))
            current_state = next_state
        enddo
        
        stochsteady_state = current_state
        
        CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                               interp_coeffs, n_states, 0.0_dp, & 
                               other_variables, 1)   

        CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                               state_coeffs, n_states, 0.0_dp, & 
                               state_variables, 1)   


        write(*,*) 'Stochastic steady state -------------------------------------' 
        write(*,*) '-------------------------------------------------------------' 

        write(*,"(A12,F10.4)") 'Capital   :', state_variables(idx_k)
        write(*,"(A12,F10.4)") 'Wealth h  :', state_variables(idx_thh) 
        write(*,"(A12,F10.4)") 'wh        :', state_variables(idx_wh) 
        write(*,"(A12,F10.4)") 'wf        :', state_variables(idx_wf) 
        write(*,"(A12,F10.4)") 'Disaster p:', exp(state_variables(idx_dis))
        write(*,"(A12,F10.4)") 'Zeta      :', state_variables(idx_zf)
        write(*,"(A12,F10.4)") 'ih        :', state_variables(idx_ih)
        write(*,"(A12,F10.4)") 'if        :', state_variables(idx_if)


        write(*,*)
        write(*,*) 'Simulations  ------------------------------------------------' 
        write(*,*) '-------------------------------------------------------------' 
        ! seed random number generator
        LSTATE  = 633
        SEED    = 712
        ifail   = 0
        call G05KFF( GENID, SUBID, SEED, LSEED, STATE, LSTATE, IFAIL)
        
        ! --------------------------------------------------------------------------
        ! First create burn-in sample 
        n_periods = 2*n_burn
        allocate( simshock_series(n_periods,1), state_collection(n_active_dims, n_burn ))

        ! generate random numbers for burn-in
        ifail = 0
        call g05saf( n_periods, STATE, simshock_series, IFAIL)

        current_state = stochsteady_state 
        
        do ttt = 1,n_periods
            
            polyn_points = smolyak_polynomial(current_state, n_active_dims, &
                                max_smol_level, smol_elem_ani) 

            
            call dgemm('n','n', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                               state_coeffs, n_states, 0.0_dp, & 
                               state_tmp(:,1), 1)   
       

            p_dis = exp(state_tmp(idx_dis,1))

            ! current weight vec
            big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
            big_weight_vec(n_quad)     = p_dis
            
            checker = 0.0_dp
            q_nxt   = 0
            do while (checker <= simshock_series(ttt,1) .and. q_nxt < n_quad) 
                q_nxt   = q_nxt + 1
                checker = checker + big_weight_vec(q_nxt) 
            enddo

            if (ttt > n_burn) then 
                state_collection(:,ttt-n_burn) = current_state(1,:)
            endif 
            
            call dgemm('n','n', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                               smol_coeffs(q_nxt,:,:), n_states, 0.0_dp, & 
                               next_state, 1)   
            where (next_state > 1.0_dp)
                next_state = 1.0_dp 
            elsewhere (next_state < -1.0_dp)
                next_state = -1.0_dp 
            endwhere

            current_state = next_state
            
        enddo
        
        ! --------------------------------------------------------------------------

        write(*,*)
        write(*,*) 'Sampling Simulation Starting Points -------------------------' 
        write(*,*) '-------------------------------------------------------------' 

        ! Generate random integers to sample starting points
        ifail = 0
        allocate(starting_locs(n_sims), starting_states(n_active_dims,n_sims))

        sample_length = n_burn
        allocate(sample_vector(sample_length))
        sample_vector = (/(ttt, ttt=1,sample_length, 1)/) 
        call g05ndf(sample_vector, sample_length, starting_locs, n_sims, STATE, IFAIL)
        ! Obtain starting states
        starting_states = state_collection(:,starting_locs)

        n_periods = n_sim_periods
        deallocate(simshock_series)
        allocate( simshock_series(n_periods,n_sims), shock_series(n_shocks, n_periods, n_sims), value_series(14,n_periods,n_sims), &
        state_series(smolyak_d, n_periods, n_sims), other_vars_series(n_interp, n_periods, n_sims) )

        ! generate random numbers:
        ifail = 0
        call g05saf( n_periods*n_sims, STATE, simshock_series, IFAIL)

        
        sample_counter = 0

        write(*,*)
        write(*,*) 'Business cycle moments  -------------------------------------' 
        write(*,*) '-------------------------------------------------------------' 
         
        ! loop through
        do sss = 1,n_sims

        
        current_state(1,:) = starting_states(:,sss)

        shock_series(:,1,sss) = shock_grid(no_shock_idx(1),:)

        do ttt = 1,n_periods
            
            polyn_points = smolyak_polynomial(current_state, n_active_dims, &
                                max_smol_level, smol_elem_ani) 

            
            call dgemm('n','n', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                               state_coeffs, n_states, 0.0_dp, & 
                               state_series(:,ttt,sss), 1)   

            call dgemm('n','n', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                               interp_coeffs, n_states, 0.0_dp, & 
                               other_vars_series(:,ttt,sss), 1)   
            
            call dgemm('n','n', 1, nrhs4, n_states, 1.0_dp, polyn_points, 1, & 
                               value_coeffs, n_states, 0.0_dp, & 
                               value_series(:,ttt,sss), 1)   

            checker = 0.0_dp
            q_nxt   = 0
            do while (checker <= simshock_series(ttt,sss) .and. q_nxt < n_quad-1) 
                q_nxt   = q_nxt + 1
                checker = checker + quad_weight_vec(q_nxt)
            enddo
            
            call dgemm('n','n', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                               smol_coeffs(q_nxt,:,:), n_states, 0.0_dp, & 
                               next_state, 1)   
            where (next_state > 1.0_dp)
                next_state = 1.0_dp 
            elsewhere (next_state < -1.0_dp)
                next_state = -1.0_dp 
            endwhere
                
           
            if (ttt<n_periods) then 
                shock_series(:,ttt+1,sss) = shock_grid(q_nxt,:)
            endif

            current_state = next_state
            
        enddo
        enddo

        open (unit = 10, file = trim(results_folder) // 'sim_state_series.txt', ACTION="write", STATUS="replace", &
                & FORM="formatted", ACCESS="STREAM")
        write(10,*) state_series; close(10)
       !   
        open (unit = 10, file = trim(results_folder) // 'sim_vars_series.txt', ACTION="write", STATUS="replace", &
                & FORM="formatted", ACCESS="STREAM")
        write(10,*) other_vars_series; close(10)
        
        open (unit = 10, file = trim(results_folder) // 'sim_value_series.txt', ACTION="write", STATUS="replace", &
                & FORM="formatted", ACCESS="STREAM")
        write(10,*) value_series; close(10)

        open (unit = 10, file = trim(results_folder) // 'sim_shock_series.txt', ACTION="write", STATUS="replace", &
        & FORM="formatted", ACCESS="STREAM")
        write(10,*) shock_series; close(10)

        write(*,*)
        write(*,*) 'Business cycle moments - including disaster -----------------' 
        write(*,*) '-------------------------------------------------------------' 

        ! loop through
        do sss = 1,n_sims
        
        current_state(1,:) = starting_states(:, sss)

        shock_series(:,1,sss) = shock_grid(no_shock_idx(1),:)

        do ttt = 1,n_periods
            
            polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                max_smol_level, smol_elem_ani) 

            
            CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                               state_coeffs, n_states, 0.0_dp, & 
                               state_series(:,ttt,sss), 1)   
            
            CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                               interp_coeffs, n_states, 0.0_dp, & 
                               other_vars_series(:,ttt,sss), 1)   
            
            CALL DGEMM('N','N', 1, nrhs4, n_states, 1.0_dp, polyn_points, 1, & 
                               value_coeffs, n_states, 0.0_dp, & 
                               value_series(:,ttt,sss), 1)   

            p_dis = exp(state_series(idx_dis, ttt, sss) )

            ! current weight vec
            big_weight_vec(1:n_quad-1) = quad_weight_vec*(1.0-p_dis)
            big_weight_vec(n_quad)     = p_dis
            
            checker = 0.0_dp
            q_nxt   = 0
            do while (checker <= simshock_series(ttt,sss) .and. q_nxt < n_quad) 
                q_nxt   = q_nxt + 1
                checker = checker + big_weight_vec(q_nxt) 
            enddo
            
            CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                               smol_coeffs(q_nxt,:,:), n_states, 0.0_dp, & 
                               next_state, 1)   
                
           
            if (ttt<n_periods) then 
                shock_series(:,ttt+1,sss) = shock_grid(q_nxt,:)
            endif
            
            current_state        = next_state
            
        enddo

        enddo


        open (unit = 10, file = trim(results_folder) // 'sim_dis_state_series.txt', ACTION="write", STATUS="replace", &
                & FORM="formatted", ACCESS="STREAM")
        write(10,*) state_series; close(10)
       !   
        open (unit = 10, file = trim(results_folder) // 'sim_dis_vars_series.txt', ACTION="write", STATUS="replace", &
                & FORM="formatted", ACCESS="STREAM")
        write(10,*) other_vars_series; close(10)

        open (unit = 10, file = trim(results_folder) // 'sim_dis_shock_series.txt', ACTION="write", STATUS="replace", &
        & FORM="formatted", ACCESS="STREAM")
        write(10,*) shock_series; close(10)

        open (unit = 10, file = trim(results_folder) // 'sim_dis_value_series.txt', ACTION="write", STATUS="replace", &
        & FORM="formatted", ACCESS="STREAM")
        write(10,*) value_series; close(10)
        

    ! ------------------------------------------------------------ !
    ! Sample n_sample IRF starting points
    ! ------------------------------------------------------------ !
    write(*,*)
    write(*,*) 'Sampling IRF Starting Points -------------------------' 
    write(*,*) '-------------------------------------------------------------' 
    sample_length = n_burn

    ! Generate random integers to sample starting points
    ifail = 0
    deallocate(starting_locs, starting_states, sample_vector)
    allocate(starting_locs(n_sample), starting_states(n_active_dims,n_sample), sample_vector(sample_length))

    sample_vector = (/(ttt, ttt=1,sample_length, 1)/) 
    call g05ndf(sample_vector, sample_length, starting_locs, n_sample, STATE, IFAIL)
    ! Obtain starting states
    if (n_sample > 1) then  
        starting_states = state_collection(:,starting_locs)
    endif    
    ! always make 1st average
    starting_states(:,1) = sum(state_collection,2)/max(1,size(state_collection,2))

    write(*,*)
    write(*,*) 'Create Impulse Responses (MIT shocks)  ----------------------' 
    write(*,*) '-------------------------------------------------------------' 



        ! impulse responses no shock
        
        n_periods = n_irf_periods
        deallocate(value_series, state_series, other_vars_series, shock_series )
        allocate(value_series(14,n_periods,n_sample),state_series(smolyak_d,n_periods,n_sample), other_vars_series(n_interp,n_periods,n_sample), shock_series(n_shocks,n_periods, n_sample) )
        allocate(irf_value_series(14,n_periods),irf_state_series(smolyak_d,n_periods), irf_other_vars_series(n_interp,n_periods), irf_shock_series(n_shocks,n_periods) )
        allocate(value_series_N(14,n_periods,n_sample),state_series_N(smolyak_d,n_periods,n_sample), other_vars_series_N(n_interp,n_periods,n_sample), shock_series_N(n_shocks,n_periods, n_sample) )
        allocate(value_series_N_helper(14,n_periods,n_sample),state_series_N_helper(smolyak_d,n_periods,n_sample), other_vars_series_N_helper(n_interp,n_periods,n_sample), shock_series_N_helper(n_shocks,n_periods, n_sample) )

          ! now loop over other IRFs
          do fff = 1,size(irf_indices,1)
          write(*,*) fff 
                
                write(shock_char,'(i1)') irf_indices(fff)
       ! Get all solution files
        open (unit = 10, file = trim(results_folder) // 'next_state_mat.dat', ACTION="read",  &
                & FORM="unformatted", ACCESS="STREAM")
        read(10) next_state_irf_mat; close(10)

        ! print all results
        open (unit = 10, file = trim(results_folder) // 'results_mat.dat', ACTION="read", &
                & FORM="unformatted", ACCESS="STREAM")
        read(10) interp_input_irf_mat; close(10)

        open (unit = 10, file = trim(results_folder) // 'valuation_mat.dat', ACTION="read", &
                & FORM="unformatted", ACCESS="STREAM")
        read(10) valuation_irf_mat; close(10)
  
  
               ! Get all solution files
                open (unit = 10, file = trim(results_folder) // 'next_state_shocktrans_mat_' // shock_char // '.dat', ACTION="read", &
                        & FORM="unformatted", ACCESS="STREAM")
                read(10) next_state_irftransition_mat; close(10)
  
  
              ! polynomial coefficients for transition matrix
              do qqq = 1,n_quad
  
                  b = transpose(next_state_irf_mat(qqq,:,:))
  
                  call F07ABF( 'E','N',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,      &
                               equed, r, c, b, ldb, smol_irf_coeffs(qqq,:,:), & 
                               ldx, rcond, ferr, berr, work, iwork, info)

  
              enddo
          
                  b = transpose(next_state_irftransition_mat(1,:,:))
  
                  call F07ABF( 'E','N',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,      &
                               equed, r, c, b, ldb, smol_irftransition_coeffs(1,:,:), & 
                               ldx, rcond, ferr, berr, work, iwork, info)
  
           
              b2 = interp_input_irf_mat
               
              call F07ABF( 'E','N',n_states, nrhs2, smol_polynom, lda,af,ldaf,ipiv,      &
                           equed, r, c, b2, ldb, interp_irf_coeffs, & 
                           ldx, rcond, ferr2, berr2, work, iwork, info )
                
               b4 = transpose(valuation_irf_mat)
                
               call F07ABF( 'E','N',n_states, nrhs4, smol_polynom, lda,af,ldaf,ipiv,      &
                            equed, r, c, b4, ldb, value_irf_coeffs, & 
                            ldx, rcond, ferr4, berr4, work, iwork, info )
  

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED( & 
        !$OMP starting_states, n_active_dims, smol_elem_ani, nrhs2, interp_coeffs, n_states, other_vars_series, max_smol_level, &
        !$OMP nrhs3, nrhs4, state_coeffs, value_coeffs, state_series, value_series, nrhs, smol_irftransition_coeffs, shock_grid, &
        !$OMP interp_irf_coeffs, smol_irf_coeffs, value_irf_coeffs, smol_coeffs, no_shock_idx, shock_series, irf_shock_sizes, irf_indices, fff, n_quad, n_periods) &
        !$OMP PRIVATE( &
        !$OMP current_state, polyn_points, next_state, ttt )
        !$OMP DO SCHEDULE(static)
        do iii = 1,n_sample
        

            ! Starting State
            current_state(1,:) = starting_states(:,iii) 
              polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                      max_smol_level, smol_elem_ani) 
  
              ! state and variables in stochastic steady state
              CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                                     interp_coeffs, n_states, 0.0_dp, & 
                                     other_vars_series(:,1,iii), 1)   
  
              CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                                     state_coeffs, n_states, 0.0_dp, & 
                                     state_series(:,1,iii), 1)   

              CALL DGEMM('N','N', 1, nrhs4, n_states, 1.0_dp, polyn_points, 1, & 
                                     value_coeffs, n_states, 0.0_dp, & 
                                     value_series(:,1,iii), 1)   
  
              ! stay one period in stochastic steady state
              other_vars_series(:,2,iii) = other_vars_series(:,1,iii)
              state_series(:,2,iii)      = state_series(:,1,iii)
              value_series(:,2,iii)      = value_series(:,1,iii)
              
              ! transition to shocked state
              CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                                     smol_irftransition_coeffs(1,:,:), n_states, 0.0_dp, & 
                                     next_state, 1)   
              
              ! in that state find state and variables
              current_state = next_state 
              polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                      max_smol_level, smol_elem_ani) 
              
              CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                                     interp_irf_coeffs, n_states, 0.0_dp, & 
                                     other_vars_series(:,3,iii), 1)   

              CALL DGEMM('N','N', 1, nrhs4, n_states, 1.0_dp, polyn_points, 1, & 
                                     value_irf_coeffs, n_states, 0.0_dp, & 
                                     value_series(:,3,iii), 1)   
  
              CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                                     state_coeffs, n_states, 0.0_dp, & 
                                     state_series(:,3,iii), 1)   
  
              ! transition out of this state
              CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                                     smol_irf_coeffs(no_shock_idx,:,:), n_states, 0.0_dp, & 
                                     next_state, 1)   
              
              shock_series(:,1,iii) = shock_grid(no_shock_idx(1),:)
              shock_series(:,2,iii) = shock_grid(no_shock_idx(1),:)
              shock_series(:,3,iii) = shock_grid(no_shock_idx(1),:)
  
              if (irf_indices(fff) == sidx_z) then 
                  shock_series(sidx_z,3,iii) = irf_shock_sizes(fff)
              elseif (irf_indices(fff) == sidx_dis) then 
                  shock_series(:,3,iii) = shock_grid(n_quad,:)
              endif
  
              current_state = next_state
  
              ! loop through
              do ttt = 4,n_periods
                   
                  shock_series(:,ttt,iii) = shock_grid(no_shock_idx(1),:)
  
                  polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                      max_smol_level, smol_elem_ani) 
                  
                  CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                                     smol_coeffs(no_shock_idx,:,:), n_states, 0.0_dp, & 
                                     next_state, 1)   
                      
                  CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                                     interp_coeffs, n_states, 0.0_dp, & 
                                     other_vars_series(:,ttt,iii), 1)   

                  CALL DGEMM('N','N', 1, nrhs4, n_states, 1.0_dp, polyn_points, 1, & 
                                     value_coeffs, n_states, 0.0_dp, & 
                                     value_series(:,ttt,iii), 1)   
                  
                  CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                                     state_coeffs, n_states, 0.0_dp, & 
                                     state_series(:,ttt,iii), 1)   
                  
                  current_state        = next_state
                  
              enddo

 
          enddo
        !$OMP END DO
        !$OMP END PARALLEL

        if (irf_indices(fff) == 0) then 
            state_series_N      = state_series
            other_vars_series_N = other_vars_series
            shock_series_N      = shock_series
            value_series_N      = value_series
            
            do ttt = 1,n_periods 
                state_series_N_helper(:,ttt,:)      = state_series_N(:,1,:)      
                other_vars_series_N_helper(:,ttt,:) = other_vars_series_N(:,1,:) 
                shock_series_N_helper(:,ttt,:)      = shock_series_N(:,1,:)      
                value_series_N_helper(:,ttt,:)      = value_series_N(:,1,:)      
            enddo

            irf_state_series      = sum(state_series_N     ,3)/n_sample 
            irf_other_vars_series = sum(other_vars_series_N,3)/n_sample 
            irf_shock_series      = sum(shock_series_N     ,3)/n_sample 
            irf_value_series      = sum(value_series_N     ,3)/n_sample 
        else 
          ! Take Averages across sample points, Shocked
          irf_state_series      = sum(state_series      - state_series_N       + state_series_N_helper      ,3)/n_sample 
          irf_other_vars_series = sum(other_vars_series - other_vars_series_N  + other_vars_series_N_helper ,3)/n_sample 
          irf_shock_series      = sum(shock_series      - shock_series_N       + shock_series_N_helper      ,3)/n_sample 
          irf_value_series      = sum(value_series      - value_series_N       + value_series_N_helper      ,3)/n_sample 
        endif

              
              if (irf_indices(fff) == sidx_z) then 
                  open (unit = 10, file = trim(results_folder) // 'g_irf_state_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_state_series; close(10)
                    
                  open (unit = 10, file = trim(results_folder) // 'g_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_other_vars_series; close(10)
  
                  open (unit = 10, file = trim(results_folder) // 'g_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_shock_series; close(10)
  
                  open (unit = 10, file = trim(results_folder) // 'g_irf_value_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_value_series; close(10)
              elseif (irf_indices(fff) == sidx_p) then 
                  open (unit = 10, file = trim(results_folder) // 'p_irf_state_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_state_series; close(10)
                    
                  open (unit = 10, file = trim(results_folder) // 'p_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_other_vars_series; close(10)
  
                  open (unit = 10, file = trim(results_folder) // 'p_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_shock_series; close(10)
  
                  open (unit = 10, file = trim(results_folder) // 'p_irf_value_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_value_series; close(10)
              elseif (irf_indices(fff) == sidx_omg) then 
                  open (unit = 10, file = trim(results_folder) // 'omg_irf_state_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_state_series; close(10)
                    
                  open (unit = 10, file = trim(results_folder) // 'omg_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_other_vars_series; close(10)
  
                  open (unit = 10, file = trim(results_folder) // 'omg_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_shock_series; close(10)

                  open (unit = 10, file = trim(results_folder) // 'omg_irf_value_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_value_series; close(10)
              elseif (irf_indices(fff) == sidx_zf) then 
                  open (unit = 10, file = trim(results_folder) // 'zf_irf_state_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_state_series; close(10)
                    
                  open (unit = 10, file = trim(results_folder) // 'zf_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_other_vars_series; close(10)
  
                  open (unit = 10, file = trim(results_folder) // 'zf_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_shock_series; close(10)

                  open (unit = 10, file = trim(results_folder) // 'zf_irf_value_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_value_series; close(10)
              elseif (irf_indices(fff) == 0) then 
                  open (unit = 10, file = trim(results_folder) // 'none_irf_state_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_state_series; close(10)
                    
                  open (unit = 10, file = trim(results_folder) // 'none_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_other_vars_series; close(10)
  
                  open (unit = 10, file = trim(results_folder) // 'none_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_shock_series; close(10)

                  open (unit = 10, file = trim(results_folder) // 'none_irf_value_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_value_series; close(10)
              elseif (irf_indices(fff) == sidx_dis) then 
                  open (unit = 10, file = trim(results_folder) // 'dis1_irf_state_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_state_series; close(10)
                    
                  open (unit = 10, file = trim(results_folder) // 'dis1_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_other_vars_series; close(10)
  
                  open (unit = 10, file = trim(results_folder) // 'dis1_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_shock_series; close(10)

                  open (unit = 10, file = trim(results_folder) // 'dis1_irf_value_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_value_series; close(10)
              endif
          
          enddo
         
  
          if (run_bg == 1) then
          write(*,*)
          write(*,*) 'Create long-lived IRFs                 ----------------------' 
          write(*,*) '-------------------------------------------------------------' 

          ! now loop over other IRFs
          do fff = 1,size(irf_indices_ll,1)
              write(*,*) fff   
              write(shock_char2,'(i2)') 10+fff
              write(iter_char, '(i4)') 1001
  
              ! Get all solution files
              open (unit = 10, file = trim(results_folder) // 'next_state_shocktrans_mat_' // shock_char2  // '_'  // iter_char  // '.dat', ACTION="read", &
                        & FORM="unformatted", ACCESS="STREAM")
              read(10) next_state_irftransition_mat; close(10)
  
              b = transpose(next_state_irftransition_mat(1,:,:))
  
              call F07ABF( 'E','N',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,      &
                            equed, r, c, b, ldb, smol_irftransition_coeffs(1,:,:), & 
                            ldx, rcond, ferr, berr, work, iwork, info)
  
           do iii = 1,1 ! no sampling, start from average 
        

            ! Starting State
            current_state(1,:) = starting_states(:,iii)
              
              polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                      max_smol_level, smol_elem_ani) 
  
              ! state and variables in stochastic steady state
              CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                                     interp_coeffs, n_states, 0.0_dp, & 
                                     other_vars_series(:,1,iii), 1)   
  
              CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                                     state_coeffs, n_states, 0.0_dp, & 
                                     state_series(:,1,iii), 1)   
  
              ! stay one period in stochastic steady state
              other_vars_series(:,2,iii) = other_vars_series(:,1,iii)
              state_series(:,2,iii)      = state_series(:,1,iii)
              
              ! transition to shocked state
              CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                                     smol_irftransition_coeffs(1,:,:), n_states, 0.0_dp, & 
                                     next_state, 1)   
  
  
              shock_series(:,1,iii) = shock_grid(no_shock_idx(1),:)
              shock_series(:,2,iii) = shock_grid(no_shock_idx(1),:)
              shock_series(:,3,iii) = shock_grid(no_shock_idx(1),:)
  
              current_state = next_state 
              do ttt = 1,irf_ll_length 
               write(iter_char, '(i4)') 1000+ttt
  
               ! Get all solution files
                open (unit = 10, file = trim(results_folder) // 'next_state_shock_mat_' // shock_char2 // '_'  // iter_char  // '.dat', ACTION="read",  &
                        & FORM="unformatted", ACCESS="STREAM")
                read(10) next_state_irf_mat; close(10)
  
                ! print all results
                open (unit = 10, file = trim(results_folder) // 'results_shock_mat_' // shock_char2  // '_'  // iter_char  // '.dat', ACTION="read",  &
                        & FORM="unformatted", ACCESS="STREAM")
                read(10) interp_input_irf_mat; close(10)
              
  
              ! polynomial coefficients for transition matrix
              do qqq = 1,n_quad
  
                  b = transpose(next_state_irf_mat(qqq,:,:))
  
                  call F07ABF( 'E','N',n_states, nrhs, smol_polynom, lda,af,ldaf,ipiv,      &
                               equed, r, c, b, ldb, smol_irf_coeffs(qqq,:,:), & 
                               ldx, rcond, ferr, berr, work, iwork, info)
  
  
              enddo
          
              b2 = interp_input_irf_mat
               
              call F07ABF( 'E','N',n_states, nrhs2, smol_polynom, lda,af,ldaf,ipiv,      &
                           equed, r, c, b2, ldb, interp_irf_coeffs, & 
                           ldx, rcond, ferr2, berr2, work, iwork, info )
                
              
              ! in that state find state and variables
              polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                      max_smol_level, smol_elem_ani) 
              
              CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                                     interp_irf_coeffs, n_states, 0.0_dp, & 
                                     other_vars_series(:,2+ttt,iii), 1)   
              
              CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                                     state_coeffs, n_states, 0.0_dp, & 
                                     state_series(:,2+ttt,iii), 1)   
  
              ! transition out of this state
              CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                                     smol_irf_coeffs(no_shock_idx,:,:), n_states, 0.0_dp, & 
                                     next_state, 1)   
              
              shock_series(:,ttt+1,iii) = shock_grid(no_shock_idx(1),:)
  
  
              current_state = next_state
              
              enddo
  
              ! loop through
              do ttt = 3+irf_ll_length,n_periods
                   
                  shock_series(:,ttt,iii) = shock_grid(no_shock_idx(1),:)
  
                  polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                      max_smol_level, smol_elem_ani) 
                  
                  CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                                     smol_coeffs(no_shock_idx,:,:), n_states, 0.0_dp, & 
                                     next_state, 1)   
                      
                  CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                                     interp_coeffs, n_states, 0.0_dp, & 
                                     other_vars_series(:,ttt,iii), 1)   
                  
                  CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                                     state_coeffs, n_states, 0.0_dp, & 
                                     state_series(:,ttt,iii), 1)   
                  
                  current_state        = next_state
                  
              enddo
            
          enddo

          ! Take Averages across sample points, Shocked
          irf_state_series      = state_series(:,:,1)      - state_series_N(:,:,1)       + state_series_N_helper(:,:,1) 
          irf_other_vars_series = other_vars_series(:,:,1) - other_vars_series_N(:,:,1)  + other_vars_series_N_helper(:,:,1) 
          irf_shock_series      = shock_series(:,:,1)      - shock_series_N(:,:,1)       + shock_series_N_helper(:,:,1)     
              
        !       if (irf_indices_ll(fff) == sidx_bg) then 
                if (fff == 1) then 
                    open (unit = 10, file = trim(results_folder) // 'bg1_irf_state_series.txt', ACTION="write", STATUS="replace", &
                            & FORM="formatted", ACCESS="STREAM")
                    write(10,*) irf_state_series; close(10)
                      
                    open (unit = 10, file = trim(results_folder) // 'bg1_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                            & FORM="formatted", ACCESS="STREAM")
                    write(10,*) irf_other_vars_series; close(10)
  
                    open (unit = 10, file = trim(results_folder) // 'bg1_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                            & FORM="formatted", ACCESS="STREAM")
                    write(10,*) irf_shock_series; close(10)
                elseif (fff == 2) then 
                  open (unit = 10, file = trim(results_folder) // 'fix_irf_state_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_state_series; close(10)
                    
                  open (unit = 10, file = trim(results_folder) // 'fix_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_other_vars_series; close(10)
  
                  open (unit = 10, file = trim(results_folder) // 'fix_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_shock_series; close(10)
                elseif (fff == 3) then 
                  open (unit = 10, file = trim(results_folder) // 'fx0_irf_state_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_state_series; close(10)
                    
                  open (unit = 10, file = trim(results_folder) // 'fx0_irf_vars_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_other_vars_series; close(10)
  
                  open (unit = 10, file = trim(results_folder) // 'fx0_irf_shock_series.txt', ACTION="write", STATUS="replace", &
                          & FORM="formatted", ACCESS="STREAM")
                  write(10,*) irf_shock_series; close(10)
                endif
          
          enddo
          endif



    if (run_samp == 1) then

    write(*,*)
    write(*,*) 'Create sample path (interpolated)  ----------------------' 
    write(*,*) '-------------------------------------------------------------' 
  
      ! establish how many shock dimensions are activated 
      counter = 0
      if (n_uni_quad_z > 1) then 
          counter = counter + 1
      endif 
      if (n_uni_quad_zf > 1) then 
          counter = counter + 1
      endif 
      if (n_uni_quad_p > 1) then 
          counter = counter + 1
      endif 
      if (n_uni_quad_omg > 1) then 
          counter = counter + 1
      endif 
      
      allocate(shock_grid_temp(counter, n_quad-2), shock_tmp(counter), sample_mat(counter,n_sample_periods))

      counter = 0 
      if (n_uni_quad_z > 1) then 
        counter = counter + 1
        shock_grid_temp(counter,1:n_quad-2) = shock_grid(1:n_quad-2,sidx_z)    
      endif 
      if (n_uni_quad_zf > 1) then 
        counter = counter + 1
        shock_grid_temp(counter,1:n_quad-2) = shock_grid(1:n_quad-2,sidx_zf)    
      endif 
      if (n_uni_quad_p > 1) then 
          counter = counter + 1
      shock_grid_temp(counter,1:n_quad-2) = shock_grid(1:n_quad-2,sidx_p)    
      endif 
      if (n_uni_quad_omg > 1) then 
          counter = counter + 1
      shock_grid_temp(counter,1:n_quad-2) = shock_grid(1:n_quad-2,sidx_omg)    
      endif 
      
      allocate(iq(2*(n_quad-2) + 1), rq((counter+1)*(counter+2)/2*(n_quad-2) + 2*counter + 1  ), placeholder(n_quad-2,1))
      n_active_shocks = counter
    
    n_periods = n_sample_periods 
    deallocate( irf_state_series, irf_other_vars_series, irf_shock_series )
    allocate( irf_state_series(smolyak_d, n_periods), irf_other_vars_series(n_interp, n_periods), irf_shock_series(n_shocks, n_periods) )
    deallocate(value_series, state_series, other_vars_series, shock_series )
    allocate(value_series(14,n_periods,n_sample),state_series(smolyak_d,n_periods,n_sample), other_vars_series(n_interp,n_periods,n_sample), shock_series(n_shocks,n_periods, n_sample) )
    deallocate(state_series_N, other_vars_series_N, state_series_N_helper, other_vars_series_N_helper)
    allocate(state_series_N(smolyak_d,n_periods,n_sample), other_vars_series_N(n_interp,n_periods,n_sample))
    allocate(state_series_N_helper(smolyak_d,n_periods,n_sample), other_vars_series_N_helper(n_interp,n_periods,n_sample))


    do fff = 1,3 
        write(*,*) fff 
        write(filename,"(A20,I0.1,A8)") '../src/params/sample', fff, '_mat.txt'
        open (unit = 10, file = trim(filename), ACTION="read",  &
                & FORM="formatted", ACCESS="sequential")
        write(n_shock_str,'(i6)') n_shocks*(n_sample_periods)
        read(10,'(' // trim(n_shock_str) // 'f100.0)') sample_mat_tmp
        close(10)
        
        counter = 0 
        if (n_uni_quad_z > 1) then 
          counter = counter + 1
          sample_mat(counter,:) = sample_mat_tmp(sidx_z,:)
        endif 
        if (n_uni_quad_zf > 1) then 
            counter = counter + 1
          sample_mat(counter,:) = sample_mat_tmp(sidx_zf,:)
        endif 
        if (n_uni_quad_p > 1) then 
            counter = counter + 1
          where (sample_mat_tmp(sidx_p,2:n_sample_periods) > maxval(shock_grid(:,sidx_p)))
          sample_mat_tmp(sidx_p,2:n_sample_periods) = maxval(shock_grid(:,sidx_p))
          elsewhere (sample_mat_tmp(sidx_p,2:n_sample_periods) < minval(shock_grid(:,sidx_p)))
          sample_mat_tmp(sidx_p,2:n_sample_periods) = minval(shock_grid(:,sidx_p))
          endwhere
          sample_mat(counter,:) = sample_mat_tmp(sidx_p,:)
        endif 
        if (n_uni_quad_omg > 1) then 
            counter = counter + 1
          where (sample_mat_tmp(sidx_omg,2:n_sample_periods) > maxval(shock_grid(:,sidx_omg)))
          sample_mat_tmp(sidx_omg,2:n_sample_periods) = maxval(shock_grid(:,sidx_omg))
          elsewhere (sample_mat_tmp(sidx_omg,2:n_sample_periods) < minval(shock_grid(:,sidx_omg)))
          sample_mat_tmp(sidx_omg,2:n_sample_periods) = minval(shock_grid(:,sidx_omg))
          endwhere
          sample_mat(counter,:) = sample_mat_tmp(sidx_omg,:)
        endif 


    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED( & 
    !$OMP starting_states, n_active_dims, smol_elem_ani, nrhs2, interp_coeffs, n_states, other_vars_series, sample_mat_tmp, max_smol_level, &
    !$OMP vector_mus_dimensions, n_active_shocks, dis_grid_mean, dis_grid_dev, omg_grid_dev, omg_grid_mean, shock_grid_temp,  qqq,  sample_mat,  &
    !$OMP nrhs3, nrhs4, state_coeffs, value_coeffs, state_series, value_series, nrhs, smol_irftransition_coeffs, shock_grid, &
    !$OMP interp_irf_coeffs, smol_irf_coeffs, value_irf_coeffs, smol_coeffs, no_shock_idx, shock_series, irf_shock_sizes, irf_indices, fff, n_quad, n_periods) &
    !$OMP PRIVATE( &
    !$OMP ddd, counter, current_state, next_state_tmp, shock_tmp, polyn_points, next_state, ttt , placeholder, ifail, iq, rq )
    !$OMP DO SCHEDULE(static)
    do iii = 1,n_sample !  sampling prohibitively slow 

        current_state(1,:) = starting_states(:,iii) 
        
        ! assign starting states for p and omega
        counter = 0
        do ddd = 1,smolyak_d 
        if (vector_mus_dimensions(ddd) > 0) then 
            counter = counter + 1 
            if (ddd == idx_dis) then  
            current_state(1,counter) = (sample_mat_tmp(sidx_p,1)- dis_grid_mean)/dis_grid_dev
            elseif (ddd == idx_omg) then  
            current_state(1,counter) = (sample_mat_tmp(sidx_omg,1)- omg_grid_mean)/omg_grid_dev
            endif
        endif
        enddo

        where (current_state > 1 - sqrt_eps) 
            current_state = 1 - sqrt_eps 
        elsewhere  (current_state < -1 + sqrt_eps)
            current_state = -1 + sqrt_eps 
        endwhere 

        ! loop through
        do ttt = 1,n_periods
            
            ! polynomial coeffs given current state 
            polyn_points = Smolyak_Polynomial(current_state, n_active_dims, &
                                max_smol_level, smol_elem_ani) 

            CALL DGEMM('N','N', 1, nrhs2, n_states, 1.0_dp, polyn_points, 1, & 
                               interp_coeffs, n_states, 0.0_dp, & 
                               other_vars_series(:,ttt,iii), 1)   
            
            CALL DGEMM('N','N', 1, nrhs3, n_states, 1.0_dp, polyn_points, 1, & 
                               state_coeffs, n_states, 0.0_dp, & 
                               state_series(:,ttt,iii), 1)   
            if (ttt > 1) then 
                shock_series(:,ttt, iii) = sample_mat_tmp(:,ttt)
            else 
                shock_series(:,ttt, iii) = shock_grid(no_shock_idx(1),:)
            endif
            
           if (ttt < n_periods) then 
            ! calculate next state for every quadrature shock given current state
            do qqq = 1, n_quad - 2
            
            CALL DGEMM('N','N', 1, nrhs, n_states, 1.0_dp, polyn_points, 1, & 
                               smol_coeffs(qqq,:,:), n_states, 0.0_dp, & 
                               next_state_tmp(qqq, :) , 1)   
            enddo      
            ! interpolate for desired shock for every state
            shock_tmp = sample_mat(:,ttt+1)
            
            do ddd = 1,n_active_dims
            ifail = -1
                call e01zmf (n_active_shocks, n_quad-2, shock_grid_temp, next_state_tmp(:,ddd), -1, -1, iq, rq, ifail)
                ifail = -1
                call e01znf (n_active_shocks, n_quad-2, shock_grid_temp, next_state_tmp(:,ddd), iq, rq, 1, shock_tmp, current_state(1,ddd), placeholder, ifail)   
            enddo 
            endif 

        where (current_state > 1 - sqrt_eps) 
            current_state = 1 - sqrt_eps 
        elsewhere  (current_state < -1 + sqrt_eps)
            current_state = -1 + sqrt_eps 
        endwhere 
        


        enddo 
        
        enddo 
        !$OMP END DO
        !$OMP END PARALLEL

        irf_state_series      = sum(state_series     ,3)/n_sample 
        irf_other_vars_series = sum(other_vars_series,3)/n_sample 
        irf_shock_series      = sum(shock_series,3)/n_sample 
        
        write(filename,"(A4,I0.1,A17)") 'samp', fff, '_state_series.txt'
        open (unit = 10, file = trim(results_folder) // trim(filename), ACTION="write", STATUS="replace", &
                & FORM="formatted", ACCESS="STREAM")
        write(10,*) irf_state_series; close(10)
          
        write(filename,"(A4,I0.1,A16)") 'samp', fff, '_vars_series.txt'
        open (unit = 10, file = trim(results_folder) // trim(filename), ACTION="write", STATUS="replace", &
                & FORM="formatted", ACCESS="STREAM")
        write(10,*) irf_other_vars_series; close(10)

        write(filename,"(A4,I0.1,A17)") 'samp', fff, '_shock_series.txt'
        open (unit = 10, file = trim(results_folder) // trim(filename), ACTION="write", STATUS="replace", &
        & FORM="formatted", ACCESS="STREAM")
        write(10,*) irf_shock_series; close(10)

        enddo

        endif

end subroutine create_results


end module mod_results

