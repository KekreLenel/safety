module mod_normal
! S15ABF = Normal CDF, F08FAF = eigenvalues and vectors, F08JEF eigen tridiagonal
use nag_library, only: S15ABF, F08FAF
use base_lib, only : dp, pi, m_choose_r, choldc
! use csv_file
implicit none
private
public ::     trunc_normal_moment_ub, normal_pdf, get_quadrature_points, trunc_normal_moment_lb, &
            & trunc_normal_moment_lb_ub

INTERFACE get_quadrature_points
    MODULE PROCEDURE    get_quadrature_points_ub, get_quadrature_points_lb, &
                    &   get_quadrature_points_lb_ub, get_quadrature_points_nb
END INTERFACE


contains

subroutine get_quadrature_points_nb(mu, sigma, n, points, weights)

real(kind=dp), intent(IN) :: mu, sigma
integer, intent(IN) :: n
real(kind=dp), intent(OUT) :: points(n), weights(n)
integer :: i, info
real(kind=dp) ::   workspace(4*n)
real(kind=dp) :: diag(n), off_diag(n-1), J_mat(n,n)

    diag(:) = 0.0_dp;


    do i = 1,n-1
        off_diag(i) = sqrt(1.0_dp*i)
    end do
    call F08JEF('I', n, diag, off_diag,  J_mat, n,  workspace, INFO)
    !call F08JEF('V', 'U', n, J_mat, n, eigenvalues, WORK, -1, INFO)
    ! allocate(dummy(int(WORK(1))))
    !call F08JEF('V', 'U', n, J_mat, n, eigenvalues, dummy, int(WORK(1)), INFO)
    points = diag
    points = mu+sigma*points
    weights = J_mat(1,1:n)**2


end subroutine get_quadrature_points_nb



subroutine get_quadrature_points_ub(mu, sigma, lb, ub, n, points, weights)

real(kind=dp), intent(IN) :: mu, sigma, ub
integer, intent(IN) :: n, lb
real(kind=dp), intent(OUT) :: points(n), weights(n)
integer :: i
real(kind=dp) ::    moments_vec(2*n+1)

    if (lb /= -99) then
        write(*,*) 'ERROR: lb is badly specified'
        stop
    end if

    do i = 0,2*n
        call trunc_normal_moment_ub(mu, sigma, ub, i, moments_vec(i+1))
    end do

    call make_points_and_weights(n, moments_vec, points, weights)

end subroutine get_quadrature_points_ub


subroutine get_quadrature_points_lb(mu, sigma, lb, ub, n, points, weights)

real(kind=dp), intent(IN) :: mu, sigma, lb
integer, intent(IN) :: n, ub
real(kind=dp), intent(OUT) :: points(n), weights(n)
integer :: i
real(kind=dp) ::    moments_vec(2*n+1)

    if (ub /= -99) then
        write(*,*) 'ERROR: ub is badly specified'
        stop
    end if

    do i = 0,2*n
        call trunc_normal_moment_lb(mu, sigma, lb, i, moments_vec(i+1))
    end do

    call make_points_and_weights(n, moments_vec, points, weights)

end subroutine get_quadrature_points_lb


subroutine get_quadrature_points_lb_ub(mu, sigma, lb, ub, n, points, weights)

real(kind=dp), intent(IN) :: mu, sigma, lb, ub
integer, intent(IN) :: n
real(kind=dp), intent(OUT) :: points(n), weights(n)
integer :: i
real(kind=dp) ::    moments_vec(2*n+1)

    do i = 0,2*n
        call trunc_normal_moment_lb_ub(mu, sigma, lb, ub, i, moments_vec(i+1))
    end do

    call make_points_and_weights(n, moments_vec, points, weights)

end subroutine get_quadrature_points_lb_ub



subroutine make_points_and_weights(n, moments_vec, points, weights)
integer, intent(IN) :: n
real(kind=dp), intent(IN) :: moments_vec(2*n+1)
real(kind=dp), intent(OUT) :: points(n), weights(n)
integer :: info, i, ifail
real(kind=dp) ::    moments_mat(n+1,n+1), alpha_vec(n), beta_vec(n-1), &
                &   J_mat(n,n), eigenvalues(n), work(1), diag(n+1)
real(kind=dp), allocatable :: dummy(:)

    call matrix_maker(n,moments_vec,moments_mat)

    ifail = 0
    call choldc(moments_mat,n+1,diag,ifail)
    if (ifail == 1) then
        write(*,*) 'WARNING: CHOLDC FAILED'
    end if

    alpha_vec(1) = moments_mat(2,1) / diag(1)
    do i = 2,n
        alpha_vec(i) = moments_mat(i+1,i) / diag(i) - moments_mat(i,i-1) / diag(i-1)
    end do

    do i = 1,n-1
        beta_vec(i) = diag(i+1) / diag(i)
    end do

    J_mat(:,:) = 0.0_dp;
    do i = 1,n
        J_mat(i,i) = alpha_vec(i)
    end do
    do i = 1,n-1
        J_mat(i,i+1) = beta_vec(i);
        J_mat(i+1,i) = beta_vec(i);
    end do

    call F08FAF('V', 'U', n, J_mat, n, eigenvalues, WORK, -1, INFO)
    allocate(dummy(int(WORK(1))))
    call F08FAF('V', 'U', n, J_mat, n, eigenvalues, dummy, int(WORK(1)), INFO)
    points = eigenvalues
    weights = J_mat(1,1:n)**2

end subroutine make_points_and_weights



! Based on notes by Phoebus J. Dhrymes
subroutine trunc_normal_moment_ub(mu, sigma, b, m, moment)
real(kind=dp), intent(in)   :: mu, sigma, b
integer, intent(in)         :: m
real(kind=dp), intent(out)  :: moment
real(kind=dp)               :: h
real(kind=dp), allocatable :: I_vec(:)
integer :: ifail, r, i, combo
    allocate(I_vec(m+1))
    h = (b-mu)/sigma
    moment = 0.0_dp
    do  i = 1,m+1
        r=i-1
        combo = m_choose_r(m,r)
        if (i==1) then
            I_vec(i) = 1.0_dp
        else
                if (i==2) then
                    ifail = 0
                    I_vec(i) = -normal_pdf(h)/S15ABF(h, ifail)
                else
                        I_vec(i) = -(h**(r-1))*normal_pdf(h)/S15ABF(h, ifail)+(r-1)*I_vec(i-2)
                end if
        end if
        if (m /= r) then
            moment = moment+combo*(mu**(m-r))*(sigma**r)*I_vec(i)
        else
            moment = moment+combo*1*(sigma**r)*I_vec(i)
        end if
    end do

end subroutine trunc_normal_moment_ub


subroutine trunc_normal_moment_lb(mu, sigma, a, m, moment)
real(kind=dp), intent(in)   :: mu, sigma, a
integer, intent(in)         :: m
real(kind=dp), intent(out)  :: moment
real(kind=dp)               :: h
real(kind=dp), allocatable :: I_vec(:)
integer :: ifail, r, i, combo
    allocate(I_vec(m+1))
    h = (a-mu)/sigma
    moment = 0.0_dp
    do  i = 1,m+1
        r=i-1
        combo = m_choose_r(m,r)
        if (i==1) then
            I_vec(i) = 1.0_dp
        else
                if (i==2) then
                    ifail = 0
                    I_vec(i) = normal_pdf(h)/(1-S15ABF(h, ifail))
                else
                        I_vec(i) = (h**(r-1))*normal_pdf(h)/(1-S15ABF(h, ifail))+(r-1)*I_vec(i-2)
                end if
        end if
        if (m /= r) then
            moment = moment+combo*(mu**(m-r))*(sigma**r)*I_vec(i)
        else
            moment = moment+combo*1*(sigma**r)*I_vec(i)
        end if
    end do

end subroutine trunc_normal_moment_lb


subroutine trunc_normal_moment_lb_ub(mu, sigma, a, b, m, moment)

use nag_library, only: S15ABF, F08FAF

real(kind=dp), intent(in)   :: mu, sigma, a, b
integer, intent(in)         :: m
real(kind=dp), intent(out)  :: moment
real(kind=dp)               :: h1, h2
real(kind=dp), allocatable :: I_vec(:)
integer :: ifail, r, i, combo
    allocate(I_vec(m+1))
    h1 = (a-mu)/sigma
    h2 = (b-mu)/sigma
    moment = 0.0_dp
    do  i = 1,m+1
        r=i-1
        combo = m_choose_r(m,r)
        if (i==1) then
            I_vec(i) = 1.0_dp
        else
                if (i==2) then
                    ifail = 0
                    I_vec(i) = -(normal_pdf(h2)-normal_pdf(h1))/(S15ABF(h2, ifail)-S15ABF(h1, ifail))
                else
                        I_vec(i) = -((h2**(r-1))*normal_pdf(h2)-(h1**(r-1))*normal_pdf(h1))/ &
                                    & (S15ABF(h2, ifail)-S15ABF(h1, ifail))+(r-1)*I_vec(i-2)
                end if
        end if
        if (m /= r) then
            moment = moment+combo*(mu**(m-r))*(sigma**r)*I_vec(i)
        else
            moment = moment+combo*1*(sigma**r)*I_vec(i)
        end if
    end do

end subroutine trunc_normal_moment_lb_ub



function normal_pdf(x)
    real(kind=dp), intent(IN)  :: x
    real(kind=dp)              :: normal_pdf
        normal_pdf = (1/(sqrt(2*pi)))*exp(-(x**2)/2)
end function normal_pdf

subroutine matrix_maker(m,vec,gram_mat)
  integer, intent(IN) :: m
  real(kind=dp), intent(IN) :: vec(2*m+1)
  real(kind=dp), intent(OUT) :: gram_mat(m+1,m+1)
  integer :: i, j

  do i = 1,m+1
     do j = 1,m+1
        gram_mat(i,j) = vec((i+j)-1)
     end do
  end do

end subroutine matrix_maker


end module mod_normal
