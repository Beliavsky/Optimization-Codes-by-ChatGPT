module conjugate_gradient_min_mod
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    private
    public :: cg_min

    contains

    subroutine cg_min(funk, xx, max_iter, tol, min_point, min_value)
        implicit none
        interface
            function funk(x)
                import :: dp
                real(kind=dp), dimension(:), intent(in) :: x
                real(kind=dp) :: funk
            end function funk
        end interface
        real(kind=dp), dimension(:), intent(in) :: xx
        integer, intent(in) :: max_iter
        real(kind=dp), intent(in) :: tol
        real(kind=dp), dimension(:), intent(out) :: min_point
        real(kind=dp), intent(out) :: min_value
        real(kind=dp), dimension(:), allocatable :: x_curr, x_new, grad_curr, grad_new, direction
        real(kind=dp) :: alpha, beta
        integer :: n, i

        ! Initialize x_curr, grad_curr, and direction
        n = size(xx)
        allocate(x_curr(n), x_new(n), grad_curr(n), grad_new(n), direction(n))

        x_curr = xx
        grad_curr = numerical_gradient(funk, x_curr)
        direction = -grad_curr

        ! Conjugate gradient optimization loop
        do i = 1, max_iter
            ! Compute step size
            alpha = line_search(funk, x_curr, grad_curr, direction)

            ! Update x_new and grad_new
            x_new = x_curr + alpha * direction
            grad_new = numerical_gradient(funk, x_new)

            ! Check convergence
            if (maxval(abs(grad_new)) < tol) exit

            ! Update beta and direction
            beta = dot_product(grad_new, grad_new - grad_curr) / dot_product(grad_curr, grad_curr)
            direction = -grad_new + beta * direction

            ! Update x_curr and grad_curr
            x_curr = x_new
            grad_curr = grad_new
        end do

        ! Set output values
        min_point = x_curr
        min_value = funk(x_curr)

        ! Deallocate memory
        deallocate(x_curr, x_new, grad_curr, grad_new, direction)
    end subroutine cg_min

    function numerical_gradient(funk, x) result(grad)
        interface
            function funk(x)
                import :: dp
                real(kind=dp), dimension(:), intent(in) :: x
                real(kind=dp) :: funk
            end function funk
        end interface
        real(kind=dp), dimension(:), intent(in) :: x
        real(kind=dp), dimension(size(x)) :: grad
        real(kind=dp), dimension(size(x)) :: x_plus, x_minus
        real(kind=dp) :: delta, eps
        integer :: i

        eps = 1.0e-8
        x_plus = x
        x_minus = x

        do i = 1, size(x)
            delta = eps * (1.0_dp + abs(x(i)))
            x_plus(i) = x(i) + delta
            x_minus(i) = x(i) - delta
            grad(i) = (funk(x_plus) - funk(x_minus)) / (2.0_dp * delta)
            x_plus(i) = x(i)
            x_minus(i) = x(i)
        end do
    end function numerical_gradient

    function line_search(funk, x_curr, grad_curr, x_new) result(step_size)
        interface
            function funk(x)
                import :: dp
                real(kind=dp), dimension(:), intent(in) :: x
                real(kind=dp) :: funk
            end function funk
        end interface
        real(kind=dp), dimension(:), intent(in) :: x_curr, grad_curr, x_new
        real(kind=dp) :: step_size
        real(kind=dp) :: alpha, beta, f_curr, f_new
        integer :: i

        alpha = 1.0_dp
        beta = 0.5_dp
        f_curr = funk(x_curr)

        do i = 1, 100
            f_new = funk(x_curr + alpha * (x_new - x_curr))
            if (f_new < f_curr) exit
            alpha = alpha * beta
        end do

        step_size = alpha
    end function line_search

end module conjugate_gradient_min_mod
