! ChatGPT-4 told to "implement the Adam optimization method."
module adam_mod
    use optim_util_mod, only: dp, numerical_gradient
    implicit none
    private
    public :: adam_min, momentum_optimizer, rmsprop, adagrad, adadelta
    
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
            v = beta2 * v + (1.0_dp - beta2) * sum(g**2)
            m_hat = m / (1.0_dp - beta1**t)
            v_hat = v / (1.0_dp - beta2**t)
            x = x - alpha * m_hat / (sqrt(v_hat) + epsilon)
        end do
        xmin = x
        fmin = funk(xmin)
    end subroutine adam_min

subroutine momentum_optimizer(funk, x0, alpha, beta, max_iter, xmin, fmin)
    ! ChatGPT-4 prompt: 'Implement the "Momentum Optimizer" optimization algorithm'.
    ! Momentum optimization method for unconstrained optimization
    interface
        function funk(x)
            import :: dp
            real(kind=dp), dimension(:), intent(in) :: x
            real(kind=dp) :: funk
        end function funk
    end interface
    real(kind=dp), dimension(:), intent(in)  :: x0
    real(kind=dp), intent(in)  :: alpha
    real(kind=dp), intent(in)  :: beta
    integer, intent(in) :: max_iter
    real(kind=dp), dimension(:), intent(out) :: xmin
    real(kind=dp), intent(out) :: fmin
    real(kind=dp), dimension(:), allocatable :: v, x, g
    real(kind=dp) :: f_old, f_new, x_new(size(x0))
    real(kind=dp), parameter :: epsilon = 1.0d-6
    integer :: iter, n

    n = size(x0)
    allocate(x(n), v(n), g(n))

    x = x0
    v = 0.0_dp
    f_old = funk(x)
    
    do iter = 1, max_iter
        g = numerical_gradient(funk, x)
        v = beta * v - alpha * g
        x_new = x + v
        f_new = funk(x_new)

        if (f_new < f_old) then
            x = x_new
            f_old = f_new
        else
            v = v / 2.0_dp
            if (norm2(v) < epsilon) exit
        end if
    end do

    xmin = x
    fmin = funk(xmin)
end subroutine momentum_optimizer

subroutine rmsprop(funk, x0, alpha, beta, epsilon, max_iter, xmin, fmin)
    ! ChatGPT-4 prompt: "Implement the Root Mean Square Propagation (RMSprop) algorithm."
    ! RMSProp optimization method for unconstrained optimization
    interface
        function funk(x)
            import :: dp
            real(kind=dp), dimension(:), intent(in) :: x
            real(kind=dp) :: funk
        end function funk
    end interface

    real(kind=dp), dimension(:), intent(in)  :: x0       ! initial guess for the parameters
    real(kind=dp), intent(in)  :: alpha     ! learning rate (a typical default value is 0.001)
    real(kind=dp), intent(in)  :: beta      ! decay rate for the moving average (a typical default value is 0.9)
    real(kind=dp), intent(in)  :: epsilon   ! small constant for numerical stability (a typical default value is 1e-8)
    integer, intent(in) :: max_iter           ! maximum number of iterations
    real(kind=dp), dimension(:), intent(out) :: xmin    ! output: minimizer
    real(kind=dp), intent(out) :: fmin        ! output: minimum function value

    real(kind=dp), dimension(:), allocatable :: v, x, g
    integer :: iter, n

    n = size(x0)
    allocate(x(n), v(n), g(n))

    x = x0
    v = 0.0_dp

    do iter = 1, max_iter
        g = numerical_gradient(funk, x)
        v = beta * v + (1.0_dp - beta) * sum(g**2)
        x = x - alpha * g / (sqrt(v) + epsilon)
    end do

    xmin = x
    fmin = funk(xmin)
end subroutine rmsprop

subroutine adagrad(funk, x0, alpha, epsilon, max_iter, xmin, fmin)
    ! Adagrad optimization method for unconstrained optimization
    interface
        function funk(x)
            import :: dp
            real(kind=dp), dimension(:), intent(in) :: x
            real(kind=dp) :: funk
        end function funk
    end interface
    real(kind=dp), dimension(:), intent(in)  :: x0       ! initial guess for the parameters
    real(kind=dp), intent(in)  :: alpha     ! learning rate (a typical default value is 0.01)
    real(kind=dp), intent(in)  :: epsilon   ! small constant for numerical stability (a typical default value is 1e-8)
    integer, intent(in) :: max_iter           ! maximum number of iterations
    real(kind=dp), dimension(:), intent(out) :: xmin    ! output: minimizer
    real(kind=dp), intent(out) :: fmin        ! output: minimum function value
    real(kind=dp), dimension(:), allocatable :: g, x, s
    integer :: iter, n
    n = size(x0)
    allocate(x(n), g(n), s(n))
    x = x0
    s = 0.0_dp
    do iter = 1, max_iter
        g = numerical_gradient(funk, x)
        s = s + sum(g**2)
        x = x - alpha * g / (sqrt(s) + epsilon)
    end do
    xmin = x
    fmin = funk(xmin)
end subroutine adagrad

subroutine adadelta(funk, x0, rho, epsilon, max_iter, xmin, fmin)
    ! ChatGPT-4 prompt: "Implement the Adadelta algorithm."
    ! Adadelta optimization method for unconstrained optimization
    interface
        function funk(x)
            import :: dp
            real(kind=dp), dimension(:), intent(in) :: x
            real(kind=dp) :: funk
        end function funk
    end interface

    real(kind=dp), dimension(:), intent(in)  :: x0      ! initial guess for the parameters
    real(kind=dp), intent(in)  :: rho       ! decay rate (a typical default value is 0.95)
    real(kind=dp), intent(in)  :: epsilon   ! small constant for numerical stability (a typical default value is 1e-6)
    integer, intent(in) :: max_iter            ! maximum number of iterations
    real(kind=dp), dimension(:), intent(out) :: xmin  ! output: minimizer
    real(kind=dp), intent(out) :: fmin         ! output: minimum function value

    real(kind=dp), dimension(:), allocatable :: g, x, s, delta_x, delta_x_prev
    integer :: iter, n

    n = size(x0)
    allocate(x(n), g(n), s(n), delta_x(n), delta_x_prev(n))

    x = x0
    s = 0.0_dp
    delta_x = 0.0_dp

    do iter = 1, max_iter
        g = numerical_gradient(funk, x)
        s = rho * s + (1.0_dp - rho) * sum(g**2)
        delta_x_prev = delta_x
        delta_x = sqrt((delta_x_prev + epsilon) / (s + epsilon)) * g
        x = x - delta_x
    end do

    xmin = x
    fmin = funk(xmin)
end subroutine adadelta

end module adam_mod
