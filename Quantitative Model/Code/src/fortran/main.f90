! -------------------------------------------------------------------------
! main.f90: main program 
! -------------------------------------------------------------------------
! authors:         Rohan Kekre and Moritz Lenel
! for updates see: https://github.com/KekreLenel/safety
! -------------------------------------------------------------------------
program main

    use nag_library, only: x06xaf, x06aaf, x06acf
    use mkl_service
    use omp_lib
    use base_lib,    only: dp, Fill_linspace_dp
    use mod_param,   only: init_setup, grid_setup
    use mod_calc,    only: calc_steady, calc_sol
    use mod_results, only: create_results 

    implicit none

    ! timing variables
    real(dp)  :: timer_start, timer_end, time_tot 
    character :: time_string*22, date*8, time*10

    write(*,*) "---------------------------------------------------------------"
    write(*,*) "---------------------------------------------------------------"
    write(*,*) "---------------------------------------------------------------"
    write(*,*) ""
    write(*,*) ""

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    ! write date and time
    call date_and_time(DATE=date,TIME=time)
    time_string = date(5:6)//'-'//date(7:8)//'-'//date(1:4)//' at '//time(1:2)//':'//time(3:4)//':'//time(5:6)
    write(*,*) 'Started on ', time_string, "."
    write(*,*) ''

    ! initialize OMP timer
    timer_start = omp_get_wtime()

    write(*,*) 'Available number of threads:', OMP_GET_NUM_PROCS()
    write(*,*) 'Requested number of threads:', OMP_GET_MAX_THREADS()
    write(*,*) 'Available threads for MKL  :', mkl_get_max_threads()   
    write(*,*) 
    
    ! set parameters
    call init_setup()
    
    ! calc steady state 
    call calc_steady()
    
    ! setup Smolyak grid
    call grid_setup()
    
    ! solve model and impulse responses 
    call calc_sol()
        
    ! simulate and create result files 
    call create_results()

    ! Stopping timer
    timer_end   = omp_get_wtime()

    ! Overall computation time:
    time_tot    = timer_end - timer_start

    write(*,*) ''
    write(*,"(A24, F12.6, A9)") ' Total computation time: ', time_tot/60.0_dp,' minutes.'
    write(*,*) '----------------------------------------------------------------'


end program main




