! This module implements the MM method of Gospodinov and Lkhagvasuren, 2013

module mod_markov
use nag_library, only: nag_wp
use base_lib, only: Fill_linspace_dp

implicit none
private
public :: calc_rouwen, make_markov_2d, make_markov_3d

integer,    parameter :: dp     = nag_wp
real(dp),   parameter :: eps    = epsilon(eps)

contains

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

subroutine make_markov_3d( A , Omg , n_grid , n_tune, stoch_nodes, stoch_probs)

    integer,        parameter  ::   n_stoch = 3

    ! Inputs
    integer,        intent(in) ::   n_grid, n_tune
    real(kind=dp),  intent(in) ::   A(n_stoch, n_stoch), Omg(n_stoch, n_stoch)

    ! Output
    real(kind=dp),  intent(out)::   stoch_nodes(n_grid**n_stoch, n_stoch), &
                                  & stoch_probs(n_grid**n_stoch, n_grid**n_stoch)

    integer       :: h, i, j, k, n_a, n_b, dummy_exceed, ix, ixx, n_ax, n_bx, &
                   & dummy_exceedx, iz, B_int(3, n_tune), n_az, n_bz, ix1, ix2, &
                   & j1, j2, j3, i1, i2, i3
    real(kind=dp) :: probtemp(n_grid,n_grid), prob1(n_grid,n_grid), z(n_grid), y(n_grid), &
                   & pmat(n_stoch, n_grid, n_grid, n_grid, n_grid), p, px, v1, v1x, &
                   & eps, mu, vact, vactx,  r, pz, vactz, stoch_nodes_temp(n_grid**n_stoch, n_stoch), &
                   & B(3, n_tune), bvectemp(n_tune), rz, probz(n_grid,n_grid)

    real(kind=dp), dimension(n_stoch, n_stoch) ::   A_new, V_y_new, V_y_old, V_e_new

    eps = epsilon(eps)

    ! Get rouwen for iid
    call calc_rouwen(0.0_dp,0.0_dp,1.0_dp, n_grid, probtemp, y)

    ! normalize the initial var so unconditional variances are 1.
    call norm_var(A, Omg, n_stoch, A_new, V_y_new, V_y_old, V_e_new)

    pmat = 0.0_dp
    px   = 0.0_dp

    do h = 1,n_grid
    do i = 1,n_grid
        do j = 1,n_grid
            do k = 1,n_stoch

                ! Conditional expected value
                mu = A_new(k,1)*y(h)+ A_new(k,2)*y(i)+ A_new(k,3)*y(j)

                ! Pick variance of current variable
                vact = V_e_new(k,k);

                ! Calculate AR coeff
                r = sqrt(1.0_dp - vact);

                ! Univariate rouwen
                call calc_rouwen(r, 0.0_dp, 1.0_dp, n_grid, prob1, z)

                ! Conditional variance
                call calc_cond_var(mu, vact, n_grid, z,  v1, p, n_a, n_b, dummy_exceed);

                ! If no tuning
                if (n_tune < 2) then

                    if ( n_a == n_b ) then ! if mu is outside of the grids
                        pmat(k,h,i,j,:) = prob1(:,n_a);
                    else ! more relevant case
                        pmat(k,h,i,j,:) = p*prob1(:,n_a)+(1-p)*prob1(:,n_b)
                    end if

                else ! tune

                    if (n_a==n_b) then  ! if mu is outside of the grids
                        pmat(k,h,i,j,:)=prob1(:,n_a);
                    else  ! begining of the more relevant

                        B   = huge(dp)
                        ixx = 0

                        do ix = 1, n_tune

                            vactx= maxval( [ sqrt(eps), vact*(1.0_dp-real((ix-1),dp)/real((n_tune-1),dp ) ) ] )

                            call calc_cond_var(mu, vactx, n_grid, z, v1x, px, n_ax,n_bx, dummy_exceedx)

                            if ( abs(dummy_exceedx)< 0.5 ) then
                                ixx = ixx+1
                                B(:,ixx)     = [v1x, px, vactx]
                                B_int(:,ixx) = [n_ax, n_bx, dummy_exceedx]
                            end if

                        end do


                        if (ixx<1) then
                            pmat(k,h,i,j,:)=p*prob1(:,n_a)+(1-p)*prob1(:,n_b)
                        else

                            bvectemp = abs(B(1,:)-vact)
                            iz = maxval(minloc(bvectemp))


                            pz      = B(2,iz)
                            vactz   = B(3, iz)
                            n_az    = B_int(1, iz)
                            n_bz    = B_int(2, iz)

                            rz      =sqrt(1-vactz)

                            call calc_rouwen(rz, 0.0_dp, 1.0_dp, n_grid, probz, z)

                            pmat(k,h,i,j,:) = pz*probz(:,n_az)+(1-pz)*probz(:,n_bz)

                        end if ! end of the more relevane
                    end if
                end if


            end do
        end do
    end do
    end do

    ! convert the transition probabilities into a conventional form
    stoch_probs = 0.0_dp

    ! combine probabilities. ! i1, i2, i3 denote current state
                             ! j1, j2, j3 denote tomorrow's state
    ix2=0
    do i1=1,n_grid
        do i2=1,n_grid
            do i3 = 1,n_grid
                ix2=ix2+1;
                do j1 = 1,n_grid
                    do j2 = 1,n_grid
                        do j3 = 1,n_grid
                        ix1 = (j1-1)*n_grid*n_grid + (j2-1)*n_grid+j3
                        stoch_probs(ix1,ix2)=pmat(1, i1,i2,i3, j1)*pmat(2, i1,i2,i3, j2) &
                                            & *pmat(2, i1,i2,i3, j3)
                        end do
                    end do
                end do
            end do
        end do
    end do

    do i = 1,n_grid**n_stoch
        stoch_probs(:,i) = stoch_probs(:,i) / sum(stoch_probs(:,i));
    end do



    ix=0
    stoch_nodes_temp = 0.0_dp
    do j1=1,n_grid
        do j2=1,n_grid
            do j3 = 1,n_grid
                ix=ix+1
                stoch_nodes_temp(ix,:) = [y(j1), y(j2), y(j3)]
            end do
        end do
    end do

    stoch_nodes(:,1) = stoch_nodes_temp(:,1)*sqrt(V_y_old(1,1))
    stoch_nodes(:,2) = stoch_nodes_temp(:,2)*sqrt(V_y_old(2,2))
    stoch_nodes(:,3) = stoch_nodes_temp(:,3)*sqrt(V_y_old(3,3))




end subroutine make_markov_3d


subroutine make_markov_2d( A , Omg , n_stoch, n_grid , n_tune, stoch_nodes, stoch_probs)

    ! Inputs
    integer,        intent(in) ::   n_grid, n_tune, n_stoch
    real(kind=dp),  intent(in) ::   A(n_stoch, n_stoch), Omg(n_stoch, n_stoch)

    ! Output
    real(kind=dp),  intent(out)::   stoch_nodes(n_grid**n_stoch, n_stoch), &
                                  & stoch_probs(n_grid**n_stoch, n_grid**n_stoch)

    integer       :: i, j, k, n_a, n_b, dummy_exceed, ix, ixx, n_ax, n_bx, &
                   & dummy_exceedx, iz, B_int(3, n_tune), n_az, n_bz, ix1, ix2, i1, i2, i3, i4
    real(kind=dp) :: probtemp(n_grid,n_grid), prob1(n_grid,n_grid), z(n_grid), y(n_grid), &
                   & pmat(n_stoch, n_grid, n_grid, n_grid), p, px, v1, v1x, &
                   & eps, mu, vact, vactx,  r, pz, vactz, stoch_nodes_temp(n_grid**n_stoch, n_stoch), &
                   & B(3, n_tune), bvectemp(n_tune), rz, probz(n_grid,n_grid)

    real(kind=dp), dimension(n_stoch, n_stoch) ::   A_new, V_y_new, V_y_old, V_e_new

    eps = epsilon(eps)

    ! Get rouwen for iid
    call calc_rouwen(0.0_dp,0.0_dp,1.0_dp, n_grid, probtemp, y)

    ! normalize the initial var so unconditional variances are 1.
    call norm_var(A, Omg, n_stoch, A_new, V_y_new, V_y_old, V_e_new)

    pmat = 0.0_dp
    px   = 0.0_dp

    do i = 1,n_grid
        do j = 1,n_grid
            do k = 1,n_stoch

                ! Conditional expected value
                mu = A_new(k,1)*y(i)+ A_new(k,2)*y(j)

                ! Pick variance of current variable
                vact = V_e_new(k,k);

                ! Calculate AR coeff
                r = sqrt(1.0_dp - vact);

                ! Univariate rouwen
                call calc_rouwen(r, 0.0_dp, 1.0_dp, n_grid, prob1, z)

                ! Conditional variance
                call calc_cond_var(mu, vact, n_grid, z,  v1, p, n_a, n_b, dummy_exceed);

                ! If no tuning
                if (n_tune < 2) then

                    if ( n_a == n_b ) then ! if mu is outside of the grids
                        pmat(k,i,j,:) = prob1(:,n_a);
                    else ! more relevant case
                        pmat(k,i,j,:) = p*prob1(:,n_a)+(1-p)*prob1(:,n_b)
                    end if

                else ! tune

                    if (n_a==n_b) then  ! if mu is outside of the grids
                        pmat(k,i,j,:)=prob1(:,n_a);
                    else  ! begining of the more relevant

                        B   = huge(dp)
                        ixx = 0

                        do ix = 1, n_tune

                            vactx= maxval( [ sqrt(eps), vact*(1.0_dp-real((ix-1),dp)/real((n_tune-1),dp ) ) ] )

                            call calc_cond_var(mu, vactx, n_grid, z, v1x, px, n_ax,n_bx, dummy_exceedx)

                            if ( abs(dummy_exceedx)< 0.5 ) then
                                ixx = ixx+1
                                B(:,ixx)     = [v1x, px, vactx]
                                B_int(:,ixx) = [n_ax, n_bx, dummy_exceedx]
                            end if

                        end do


                        if (ixx<1) then
                            pmat(k,i,j,:)=p*prob1(:,n_a)+(1-p)*prob1(:,n_b)
                        else

                            bvectemp = abs(B(1,:)-vact)
                            iz = maxval(minloc(bvectemp))


                            pz      = B(2,iz)
                            vactz   = B(3, iz)
                            n_az    = B_int(1, iz)
                            n_bz    = B_int(2, iz)

                            rz      =sqrt(1-vactz)

                            call calc_rouwen(rz, 0.0_dp, 1.0_dp, n_grid, probz, z)

                            pmat(k,i,j,:) = pz*probz(:,n_az)+(1-pz)*probz(:,n_bz)

                        end if ! end of the more relevane
                    end if
                end if


            end do
        end do
    end do


    ! convert the transition probabilities into a conventional form
    stoch_probs = 0.0_dp

    ix2=0
    do i1=1,n_grid
        do i2=1,n_grid
            ix2=ix2+1;
            do i3=1,n_grid
                do i4=1,n_grid
                    ix1 = (i3-1)*n_grid+i4
                    stoch_probs(ix1,ix2)=pmat(1,i1,i2,i3)*pmat(2,i1,i2,i4)
                end do
            end do
        end do
    end do

    do i = 1,n_grid*n_grid
        stoch_probs(:,i) = stoch_probs(:,i) / sum(stoch_probs(:,i));
    end do



    ix=0
    stoch_nodes_temp = 0.0_dp
    do i=1,n_grid
        do j=1,n_grid
            ix=ix+1
            stoch_nodes_temp(ix,:) = [y(i), y(j)]
        end do
    end do

    stoch_nodes(:,1) = stoch_nodes_temp(:,1)*sqrt(V_y_old(1,1))
    stoch_nodes(:,2) = stoch_nodes_temp(:,2)*sqrt(V_y_old(2,2))




end subroutine make_markov_2d



! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

subroutine calc_rouwen(rho_Rouw, mu_uncond, sig_uncond, n_R, P_Rouw, z_Rouw)

    ! Inputs
    integer,        intent(in) ::   n_R
    real(kind=dp),  intent(in) ::   rho_Rouw, mu_uncond, sig_uncond

    ! Outputs
    real(kind=dp),  intent(out)::   P_Rouw(n_R, n_R), z_Rouw(n_R)

    real(kind=dp) :: step_R, p
    real(kind=dp), dimension(n_R,n_R) :: a1R, a2R, a3R, a4R
    integer       :: iii


    if ( abs(rho_Rouw)>1.0_dp ) then
        write(*,*) ''
        write(*,*) 'WARNING:'
        write(*,*) 'The persistence parameter, rho, must be less than one in absolute value.'
        write(*,*) ''
        stop
    end if

    ! Create grids
    step_R = sig_uncond*sqrt(real(n_R,dp) - 1.0_dp)


    call Fill_linspace_dp(-1.0_dp, 1.0_dp, z_Rouw)

    z_Rouw = mu_uncond+step_R*z_Rouw

    ! CONSTRUCTION OF THE TRANSITION PROBABILITY MATRIX
    p = (rho_Rouw + 1.0_dp)/2.0_dp

    P_Rouw(1,1) = p
    P_Rouw(1,2) = 1-p
    P_Rouw(2,1) = 1-p
    P_Rouw(2,2) = p

    a1R = 0.0_dp
    a2R = 0.0_dp
    a3R = 0.0_dp
    a4R = 0.0_dp

    do iii = 2,n_R-1

        a1R(1:iii,1:iii)    = P_Rouw(1:iii,1:iii) ! of interest iii +1
        a2R(1:iii,2:iii+1)  = P_Rouw(1:iii,1:iii)
        a3R(2:iii+1,1:iii)  = P_Rouw(1:iii,1:iii)
        a4R(2:iii+1,2:iii+1)= P_Rouw(1:iii,1:iii)

        P_Rouw = p*a1R+(1-p)*a2R+(1-p)*a3R+p*a4R
        P_Rouw(2:iii, :) = P_Rouw(2:iii, :)/2
    end do

    P_Rouw = transpose(P_Rouw)

    do iii = 1,n_R
        P_Rouw(:,iii) = P_Rouw(:,iii)/sum(P_Rouw(:,iii))
    end do


end subroutine calc_rouwen


subroutine norm_var(A, Omg, n_stoch, A_new, V_y_new, V_y_old, V_e_new)

    ! Inputs
    integer,        intent(in)  :: n_stoch
    real(kind=dp),  intent(in)  :: A(n_stoch, n_stoch), Omg(n_stoch, n_stoch)

    ! Outputs
    real(kind=dp), intent(out)  :: A_new(n_stoch, n_stoch), V_y_new(n_stoch, n_stoch), &
                                 & V_y_old(n_stoch, n_stoch), V_e_new(n_stoch, n_stoch)

    ! variable declarations
    integer         :: i, j
    real(kind=dp)   :: dif, V0(n_stoch, n_stoch), V(n_stoch, n_stoch), eps

    eps = epsilon(eps)

    dif = 100_dp
    V0  = 0.0_dp

    ! Find unconditional
    do while ( dif > sqrt(eps) )
        V = matmul( matmul(A,V0) , transpose(A) ) + Omg
        dif=  maxval( V-V0 );
        V0 = V
    end do

    V_y_old = V0

    V_e_new = 0.0_dp
    A_new   = 0.0_dp

    do i=1, n_stoch
        V_e_new(i,i)=Omg(i,i)/V_y_old(i,i);
        do j = 1, n_stoch
            A_new(i,j)=A(i,j)*sqrt(V_y_old(j,j))/sqrt(V_y_old(i,i));
        end do
    end do

    V_y_new = 0.0_dp

    do i=1, n_stoch
        do j=1, n_stoch
        V_y_new(i,j)=V_y_old(i,j)/(sqrt(V_y_old(i,i))*sqrt(V_y_old(j,j)) );
        end do
    end do



end subroutine norm_var



! calc_cond_var
! This function calculates the conditional variance of the mixture
! distribution given the conditional mean mu and the conditional variance
! v0 of the mass distributions on the n grids given by z.

! For details, see Nikolay Gospodinov and Damba Lkhagvasuren, 2013

subroutine calc_cond_var( mu, v0, n_grid, z,  v1, p, n_a, n_b, dummy_exceed)

    ! Input
    integer, intent(in)       :: n_grid
    real(kind=dp), intent(in) :: mu, v0, z(n_grid)

    ! Output
    integer, intent(out)      :: n_a, n_b, dummy_exceed
    real(kind=dp), intent(out):: v1, p

    ! Declare vars
    real(kind=dp) :: r, zm(n_grid), zax, zbx

    r = sqrt(1.0_dp-v0);

    zm=z*r;
    if (mu+eps>=zm(n_grid)) then
        dummy_exceed=1;
        n_a = n_grid
        n_b = n_grid
        p   = 0.0_dp
        v1  = v0
    elseif (mu-eps<=zm(1)) then
        dummy_exceed = -1
        n_a = 1
        n_b = 1
        p   = 1.0_dp
        v1  = v0
    else
        dummy_exceed=0

        n_a = 1 + floor( (mu-zm(1))/(zm(2)-zm(1)) )

        n_b = n_a+1
        zax = zm(n_a)
        zbx = zm(n_b)

        p = (zbx-mu)/(zbx-zax)

        v1 = v0+p*(1-p)*((zbx-zax)**2)

    end if

end subroutine calc_cond_var


end module mod_markov