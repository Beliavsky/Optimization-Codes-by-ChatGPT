! Prompt givent to ChatGPT-4:
! "Write a subroutine that uses gradient descent with the arguments
! 
!     subroutine random_min(funk, xh, xx, min_value)
! 
! where xx(:) is an intent(in out) argument that is a guess for the function argument on input and is the best function argument on output and where xh is a parameter that determines when the value of xx has converged."
! The prompts given to implement the other methods were shorter but followed this prompt
! in the same conversation
module gradient_descent_min_mod
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    private
    public :: gradient_descent_min, dp

    contains

    subroutine gradient_descent_min(funk, xh, step_size, iter_max, xx, min_value, &
                                    grad, niter, nfunc_eval)
        implicit none
        interface
            function funk(x)
                import :: dp
                real(kind=dp), dimension(:), intent(in) :: x
                real(kind=dp) :: funk
            end function funk
        end interface
        real(kind=dp), intent(in) :: xh
        real(kind=dp), intent(in) :: step_size
        integer      , intent(in) :: iter_max
        real(kind=dp), dimension(:), intent(inout) :: xx
        real(kind=dp), intent(out) :: min_value
        real(kind=dp), intent(out) :: grad(:)
        integer      , intent(out) :: niter, nfunc_eval
        real(kind=dp) :: norm_grad, xtry(size(xx)), xchange(size(xx))
        integer :: n
        integer, parameter :: nfunc_eval_max = 10**6
        niter = 0
        nfunc_eval = 0
        n = size(xx)

        main_loop: do
            ! Calculate gradient
            if (niter >= iter_max) exit
            call gradient(funk, xx, grad)
            niter = niter + 1
            ! print*,"xx =",xx ! debug
            ! Check for convergence
            norm_grad = sqrt(sum(grad**2))
            if (norm_grad < xh) exit

            ! Update xx using gradient descent
            xchange = step_size * grad
            change_size_loop: do
               nfunc_eval = nfunc_eval + 1
               if (nfunc_eval >= nfunc_eval_max) exit main_loop
               xtry = xx - xchange
               if (funk(xtry) < funk(xx)) then
                  xx = xtry
                  exit change_size_loop
               else
                  xchange = xchange/2
               end if
            end do change_size_loop
        end do main_loop

        ! Calculate the minimum value
        min_value = funk(xx)
    end subroutine gradient_descent_min

    subroutine gradient(funk, x, grad)
        interface
            function funk(x)
                import :: dp
                real(kind=dp), dimension(:), intent(in) :: x
                real(kind=dp) :: funk
            end function funk
        end interface
        real(kind=dp), dimension(:), intent(inout) :: x
        real(kind=dp), intent(out) :: grad(size(x))
        real(kind=dp) :: h, f1, f2
        integer :: i
        h = 1.0e-8_dp

        do i = 1, size(x)
            x(i) = x(i) + h
            f1 = funk(x)
            x(i) = x(i) - 2.0_dp * h
            f2 = funk(x)
            x(i) = x(i) + h
            grad(i) = (f1 - f2) / (2.0_dp * h)
        end do
    end subroutine gradient

end module gradient_descent_min_mod
