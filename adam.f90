! ChatGPT-4 told to "implement the Adam optimization method."
module adam_mod
    use optim_util_mod, only: dp, numerical_gradient
    implicit none
    private
    public :: adam_min
    
contains

    subroutine adam_min(funk, x0, alpha, beta1, beta2, epsilon, max_iter, xmin, fmin)
        ! Adam optimization method for unconstrained optimization

        interface
            function funk(x)
                import :: dp
                real(kind=dp), dimension(:), intent(in) :: x
                real(kind=dp) :: funk
            end function funk
        end interface

        real(kind=dp), dimension(:), intent(in)  :: x0     ! Initial guess for the solution
        real(kind=dp), intent(in)  :: alpha                ! Step size (learning rate), default: 0.001
        real(kind=dp), intent(in)  :: beta1                ! Decay rate for the first moment estimate, default: 0.9
        real(kind=dp), intent(in)  :: beta2                ! Decay rate for the second moment estimate, default: 0.999
        real(kind=dp), intent(in)  :: epsilon              ! Small constant for numerical stability, default: 1.0e-8
        integer, intent(in) :: max_iter                    ! Maximum number of iterations
        real(kind=dp), dimension(:), intent(out) :: xmin   ! Optimal solution
        real(kind=dp), intent(out) :: fmin                 ! Minimum function value

        real(kind=dp), dimension(:), allocatable :: m, v, x, g, m_hat, v_hat
        real(kind=dp) :: t
        integer :: iter, n
        n = size(x0)
        allocate(x(n), m(n), v(n), g(n), m_hat(n), v_hat(n))
        x = x0
        m = 0.0_dp
        v = 0.0_dp
        t = 0.0_dp
        do iter = 1, max_iter
            g = numerical_gradient(funk, x)
            t = t + 1.0_dp
            m = beta1 * m + (1.0_dp - beta1) * g
            v = beta2 * v + (1.0_dp - beta2) * g**2
            m_hat = m / (1.0_dp - beta1**t)
            v_hat = v / (1.0_dp - beta2**t)
            x = x - alpha * m_hat / (sqrt(v_hat) + epsilon)
        end do
        xmin = x
        fmin = funk(xmin)
    end subroutine adam_min
end module adam_mod
