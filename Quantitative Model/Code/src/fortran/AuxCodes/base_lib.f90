module base_lib
use omp_lib
use nag_library, only: nag_wp
implicit none
private
public :: Fill_linspace_dp, m_choose_r, pi, dp, eps, sqrt_eps, choldc, updateBrent, & 
          QsortC, calc_var, calc_cov, remove_dups_2d, remove_dups_1d_real, initBrent

integer, parameter   :: dp       = nag_wp

real(dp), parameter  :: pi       = 4*atan(1.0_dp), &
                     &  eps      = epsilon(eps),   &
                     &  sqrt_eps = sqrt(eps)

contains

   ! Stolen from https://github.com/astrofrog/fortranlib/blob/master/src/lib_array.f90
   subroutine Fill_linspace_dp(xmin,xmax,x)
      implicit none
      real(dp),intent(in) :: xmin,xmax
      real(dp),intent(inout) :: x(:)
      integer :: i,n
      n = size(x)
      if (n == 1) then
         if(xmin /= xmax) then
            write(0,'("ERROR: Cannot call linspace with n=1 and xmin /= xmax")')
            stop
         else
            x = xmin
         end if
      else
         x(1) = xmin
         do i=2,n-1
             x(i) = (xmax-xmin) * real(i-1,dp) / real(n-1,dp) + xmin
         end do
         x(n) = xmax
      end if
   end subroutine Fill_linspace_dp


    RECURSIVE FUNCTION m_choose_r(m,r_temp) RESULT(res)
    INTEGER :: m, r_temp, r, res
    IF (m < r_temp) THEN
         write(*,*) 'Warning M_CHOOSE_R'
         stop
    END IF

    r = min(m-r_temp, r_temp)
    IF (m == r) THEN
         res = 1
    ELSE IF (r == 0) THEN
         res = 1
    ELSE IF (r == 1) THEN
         res = m
    ELSE
         res = m*m_choose_r_second(m-1,r-1)/r
    END IF
   END FUNCTION m_choose_r

    RECURSIVE FUNCTION m_choose_r_second(m,r) RESULT(res)
    INTEGER :: m, r, res

    IF (r == 1) THEN
         res = m
    ELSE
         res = m*m_choose_r_second(m-1,r-1)/r
    END IF
    END FUNCTION m_choose_r_second


   SUBROUTINE choldc(a,n,p,ifail)
   INTEGER, intent(IN) :: n
   REAL(kind=dp), intent(INOUT) :: a(n,n)
   REAL(kind=dp), intent(OUT) :: p(n)
   INTEGER, intent(INOUT) :: ifail
   INTEGER :: i,j,k
   REAL(kind=dp) :: sum
   do i=1,n
      do j=i,n
         sum=a(i,j)
         do k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
         enddo
         if(i.eq.j)then
            if(sum.le.0.) then
               ifail = 1
               sum = sqrt(sum**2)
            end if
            p(i)=sqrt(sum)
         else
            a(j,i)=sum/p(i)
         endif
      enddo
   enddo
   return
   END

recursive subroutine QsortC(A)
  real(dp), intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real(dp), intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real(dp) :: temp
  real(dp) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition


function calc_mean(vec, n_vec) result(mean)
    integer, intent(in)   :: n_vec
    real(dp), intent(in)  :: vec(n_vec)
    real(dp) :: mean

    mean = sum(vec)/real(n_vec,dp)

end function calc_mean

function calc_var(vec, n_vec) result(var)
    integer, intent(in)   :: n_vec
    real(dp), intent(in)  :: vec(n_vec)
    real(dp) :: var

    real(dp) :: mean

    mean = calc_mean(vec, n_vec)
    var  = sum((vec-mean)**2)/real(n_vec-1,dp)

end function calc_var

function calc_cov(vec1, vec2, n_vec) result(cov)
    integer, intent(in)   :: n_vec
    real(dp), intent(in)  :: vec1(n_vec), vec2(n_vec)
    real(dp) :: cov

    real(dp) :: mean1, mean2

    mean1 = calc_mean(vec1, n_vec)
    mean2 = calc_mean(vec2, n_vec)
    cov  = sum((vec1-mean1)*(vec2-mean2))/real(n_vec-1,dp)

end function calc_cov

! adapted from https://rosettacode.org for 2d matrices along rows
function remove_dups_2d(input_mat) result(output_mat)

    implicit none
    integer, intent(in)  :: input_mat(:,:)  ! The input
    integer, allocatable :: output_mat(:,:) ! The output
    integer :: k                            ! The number of unique elements
    integer :: i, j
 
    k = 1
    output_mat = input_mat
    output_mat = -9999
    output_mat(1,:) = input_mat(1,:)
    outer: do i=2,size(input_mat,1)
        do j=1,k
            if (all(output_mat(j,:).eq.input_mat(i,:))) then
                ! Found a match so start looking again
                cycle outer
            end if
        end do
        ! No match found so add it to the output
        k = k + 1
        output_mat(k,:) = input_mat(i,:)
    end do outer
  
  output_mat = output_mat(1:k,:)

end function remove_dups_2d


! adapted from https://rosettacode.org for 2d matrices along rows
function remove_dups_1d_real(input_mat) result(output_mat)

    implicit none
    real(dp), intent(in)  :: input_mat(:)  ! The input
    real(dp), allocatable :: output_mat(:) ! The output
    integer :: k                          ! The number of unique elements
    integer :: i, j
 
    k = 1
    output_mat = input_mat
    output_mat = -9999
    output_mat(1) = input_mat(1)
    outer: do i=2,size(input_mat)
        do j=1,k
            if ( output_mat(j) == input_mat(i) ) then
                ! Found a match so start looking again
                cycle outer
            end if
        end do
        ! No match found so add it to the output
        k = k + 1
        output_mat(k) = input_mat(i)
    end do outer
  
    output_mat = output_mat(1:k)

end function remove_dups_1d_real



subroutine initBrent(xbracket, ybracket, xnew, ynew, xvec, yvec, mflag, sflag, iter)
    
    implicit none 
    real(dp), intent(inout) :: xbracket(2), ybracket(2)
    real(dp), intent(out)   :: xvec(4), yvec(3), xnew, ynew
    logical, intent(out)    :: mflag
    integer, intent(out)    :: iter, sflag

    ! check that root is bracketed
    if ((ybracket(1)*ybracket(2)) > 0) then 
        write(*,*) "ERROR: Initial bracket does not contain root."
        stop
    endif

    if (abs(ybracket(1)) < abs(ybracket(2))) then
        ! Swap lower and upper bounds.
        xbracket = xbracket(2:1:-1) 
        ybracket = ybracket(2:1:-1) 
    endif

    ! prepare brent storage
    xvec(1:2) = xbracket 
    xvec(3:4) = xbracket(1)
    yvec(1:2) = ybracket
    yvec(3)   = ybracket(1)

    ! prepare flags
    xnew    = 0.0
    ynew    = 0.0
    mflag   = .true.
    sflag   = 0
    iter    = 0

end subroutine



subroutine updateBrent(xnew, ynew, x, y, mflag, flag, iter, xtol, ytol, maxiter)

    implicit none 
    real(dp), intent(inout)   :: x(4), y(3), xnew, ynew
    logical, intent(inout)    :: mflag
    integer, intent(inout)    :: iter, flag
    real(dp), intent(in)      :: xtol, ytol
    integer, intent(in)       :: maxiter 
    real(dp)                  :: delta, min_1, min_2, min_3
        iter = iter + 1
       
        if (iter > 1) then  ! if already updated once

            ! update history
            x(4) = x(3)
            x(3) = x(2)
            y(3) = y(2)
    
            if (y(1)*ynew < 0) then 
                x(2) = xnew
                y(2) = ynew
            else
                x(1) = xnew
                y(1) = ynew
            endif
    
            if (abs(y(1)) < abs(y(2))) then 
                ! Swap lower and upper bounds.
                x(1:2) = x(2:1:-1) 
                y(1:2) = y(2:1:-1) 
            endif
    
            if (iter > maxiter) then 
                flag = 3
                x = x(2)
                y = y(2)
            elseif (abs(x(2)-x(1)) < xtol) then 
                flag = 2 
                x = x(2)
                y = y(2)
            elseif (abs(ynew) < ytol) then 
                flag = 1
                x = x(2)
                y = y(2)
            else 
                flag = 0
            endif
        else 
            flag = 0
        endif 

        if (flag == 0) then ! if no solution yet 

            ! inverse quadratic interpolation if f(x(1))!=f(x(2))!=f(x(3))
            ! and linear interpolation (secant method) otherwise.
            if (abs(y(1)-y(3)) > ytol .and. abs(y(2)-y(3)) > ytol) then 
                xnew = x(1)*y(2)*y(3)/((y(1)-y(2))*(y(1)-y(3))) + & 
                       x(2)*y(1)*y(3)/((y(2)-y(1))*(y(2)-y(3))) + & 
                       x(3)*y(1)*y(2)/((y(3)-y(1))*(y(3)-y(2)))
            else
                xnew = x(2) - y(2) * (x(2)-x(1))/(y(2)-y(1))
            endif

            delta   = abs(2*eps*abs(x(2)))
            min_1    = abs(xnew-x(2))
            min_2    = abs(x(2)-x(3))
            min_3    = abs(x(3)-x(4))

            ! use bisection method if satisfies one of the conditions.
            if (( ((xnew > (3.0*x(1)+x(2))/4.0) .and. (xnew < x(2))    )     == .false.) .or. & 
                ( ((mflag == .true.           ) .and. (min_1 >= min_2/2.0) ) == .true. ) .or. &
                ( ((mflag == .false.          ) .and. (min_1 >= min_3/2.0) ) == .true. ) .or. &
                ( ((mflag == .true.           ) .and. (min_2 < delta)   )    == .true. ) .or. &
                ( ((mflag == .false.          ) .and. (min_3 < delta)   )    == .true. ) ) then 
                xnew = (x(1)+x(2))/2
                mflag = .true.
            else
                mflag = .false.
            endif

        endif

end subroutine


!!! #### TEST BRENT 
!!! function testBrent()
!!! 
!!!     f(x) = (x+3)*(x-1)^2
!!!     xbracket = [-4.0, 4/3]
!!!     ybracket = funfun.(xbracket)
!!!     
!!!     #initialize brent
!!!     xnew, ynew, xvec, yvec, mflag, sflag, iter = initBrent(xbracket, ybracket)
!!! 
!!!     while sflag == 0 
!!!         
!!!         # get next brent
!!!         xnew, mflag, sflag, iter = updateBrent!(xnew, ynew, xvec, yvec, mflag, iter)
!!! 
!!!         # update function value
!!!         ynew = funfun(xnew)
!!! 
!!! 
!!!     end
!!!     
!!!     println([xvec[1],yvec[1]])
!!! 
!!! end



end module base_lib
