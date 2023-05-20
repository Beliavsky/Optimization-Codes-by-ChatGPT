module simulated_annealing_mod
use optim_util_mod, only: dp
implicit none
private
public :: simulated_annealing
contains
subroutine simulated_annealing(funk, x0, T0, T_min, alpha, max_iter, xmin, fmin)
    ! ChatGPT-4 prompt: "Implement simulated annealing."
    ! Simulated Annealing algorithm for unconstrained optimization
    interface
        function funk(x)
            import :: dp
            real(kind=dp), dimension(:), intent(in) :: x
            real(kind=dp) :: funk
        end function funk
    end interface
real(kind=dp), dimension(:), intent(in)  :: x0    ! Initial guess for the minimum point
real(kind=dp), intent(in)  :: T0                  ! Initial temperature, should be set considering the problem's scale. Reasonable default could be 1.0_dp
real(kind=dp), intent(in)  :: T_min               ! Minimum temperature, annealing will stop when temperature is below this. Reasonable default is 1e-3_dp
real(kind=dp), intent(in)  :: alpha               ! Cooling rate, a value between (0, 1), typically closer to 1, maybe 0.9_dp
integer, intent(in) :: max_iter                   ! Maximum number of iterations
real(kind=dp), dimension(:), intent(out) :: xmin  ! Minimum point found
real(kind=dp), intent(out) :: fmin                ! Function value at the minimum point
    real(kind=dp), dimension(:), allocatable :: x, x_new
    real(kind=dp) :: T, f, f_new, p, rand_num
    integer :: i, n
    n = size(x0)
    allocate(x(n), x_new(n))
    x = x0
    f = funk(x)
    T = T0
    do i = 1, max_iter
        call random_number(rand_num)
        x_new = x + T * (2.0_dp * rand_num - 1.0_dp)
        f_new = funk(x_new)

        if (f_new < f) then
            x = x_new
            f = f_new
        else
            call random_number(rand_num)
            p = exp((f - f_new) / T)
            if (rand_num < p) then
                x = x_new
                f = f_new
            end if
        end if
        T = alpha * T
        if (T < T_min) exit
    end do
    xmin = x
    fmin = funk(xmin)
end subroutine simulated_annealing
end module simulated_annealing_mod