! ChatGPT-4 prompt: "Implement the Hooke-Jeeves direct search method."
module hooke_jeeves_mod
use optim_util_mod, only: dp
implicit none
public :: hooke_jeeves
contains
subroutine hooke_jeeves(funk, x0, step, epsilon, max_iter, xmin, fmin)
    ! Hooke-Jeeves direct search method for unconstrained optimization
    interface
        function funk(x)
            import :: dp
            real(kind=dp), dimension(:), intent(in) :: x
            real(kind=dp) :: funk
        end function funk
    end interface
    real(kind=dp), dimension(:), intent(in)  :: x0
    real(kind=dp), intent(in out) :: step
    real(kind=dp), intent(in)  :: epsilon
    integer, intent(in) :: max_iter
    real(kind=dp), dimension(:), intent(out) :: xmin
    real(kind=dp), intent(out) :: fmin
    real(kind=dp), dimension(:), allocatable :: x, x_old, x_new
    integer :: iter, i, n
    n = size(x0)
    allocate(x(n), x_old(n), x_new(n))
    x = x0
    x_new = x
    fmin = funk(x)
    do iter = 1, max_iter
        x_old = x
        ! Exploratory move
        do i = 1, n
            x_new(i) = x(i) + step
            if (funk(x_new) < fmin) then
                x(i) = x_new(i)
                fmin = funk(x_new)
            else
                x_new(i) = x(i) - step
                if (funk(x_new) < fmin) then
                    x(i) = x_new(i)
                    fmin = funk(x_new)
                else
                    x_new(i) = x(i)
                end if
            end if
        end do
        ! Pattern move
        x_new = x + (x - x_old)
        if (funk(x_new) < fmin) then
            x = x_new
            fmin = funk(x_new)
        endif
        ! Check for convergence
        if (maxval(abs(x - x_old)) < epsilon) exit
        ! Reduce step size
        if (funk(x_new) >= funk(x_old)) step = step / 2.0_dp
    end do
    xmin = x
end subroutine hooke_jeeves
end module hooke_jeeves_mod