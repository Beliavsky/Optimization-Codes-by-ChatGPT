! ChatGPT-4 told to "implement the BFGS algorithm."
module bfgs_min_mod
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    private
    public :: bfgs_min

    contains

    subroutine bfgs_min(funk, xx, max_iter, tol, min_point, min_value)
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
        real(kind=dp), dimension(:), allocatable :: x_curr, x_new
        real(kind=dp), dimension(:), allocatable :: grad_curr, grad_new
        real(kind=dp), dimension(:,:), allocatable :: hessian_inv
        real(kind=dp) :: f_curr, f_new, step_size, delta_f
        integer :: n, i

        ! Initialize x_curr, f_curr, grad_curr, and hessian_inv
        n = size(xx)
        allocate(x_curr(n), x_new(n), grad_curr(n), grad_new(n))
        allocate(hessian_inv(n, n))

        x_curr = xx
        f_curr = funk(x_curr)
        grad_curr = numerical_gradient(funk, x_curr)
        hessian_inv = eye(n)

        ! BFGS optimization loop
        do i = 1, max_iter
            ! Compute the search direction and step size
            x_new = x_curr - matmul(hessian_inv, grad_curr)
            step_size = line_search(funk, x_curr, grad_curr, x_new)

            ! Update x_new and f_new
            x_new = x_curr + step_size * (x_new - x_curr)
            f_new = funk(x_new)

            ! Check convergence
            delta_f = abs(f_new - f_curr)
            if (delta_f < tol) exit

            ! Update grad_new
            grad_new = numerical_gradient(funk, x_new)

            ! Update hessian_inv
            call update_hessian_inv(hessian_inv, x_curr, x_new, grad_curr, grad_new)

            ! Update x_curr, f_curr, and grad_curr
            x_curr = x_new
            f_curr = f_new
            grad_curr = grad_new
        end do

        ! Set output values
        min_point = x_curr
        min_value = f_curr

        ! Deallocate memory
        deallocate(x_curr, x_new, grad_curr, grad_new, hessian_inv)
    end subroutine bfgs_min

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

subroutine update_hessian_inv(hessian_inv, x_curr, x_new, grad_curr, grad_new)
    real(kind=dp), dimension(:,:), intent(inout) :: hessian_inv
    real(kind=dp), dimension(:), intent(in) :: x_curr, x_new, grad_curr, grad_new
    real(kind=dp), dimension(size(x_curr)) :: s, y
    real(kind=dp) :: rho

    s = x_new - x_curr
    y = grad_new - grad_curr
    rho = 1.0_dp / dot_product(y, s)

    hessian_inv = hessian_inv - rho * outer_product(s, s) &
                  + (1.0_dp + rho * dot_product(y, matmul(hessian_inv, y))) * outer_product(s, s) * rho &
                  - rho * (outer_product(s, matmul(transpose(hessian_inv), y)) + outer_product(matmul(y, hessian_inv), s))
end subroutine update_hessian_inv

function eye(n) result(identity_matrix)
    integer, intent(in) :: n
    real(kind=dp), dimension(n, n) :: identity_matrix
    integer :: i

    identity_matrix = 0.0_dp
    do i = 1, n
        identity_matrix(i, i) = 1.0_dp
    end do
end function eye

function outer_product(a, b) result(outer_prod)
    integer, parameter :: dp = kind(1.0d0)
    real(kind=dp), dimension(:), intent(in) :: a, b
    real(kind=dp), dimension(size(a), size(b)) :: outer_prod
    integer :: i, j

    do i = 1, size(a)
        do j = 1, size(b)
            outer_prod(i, j) = a(i) * b(j)
        end do
    end do
end function outer_product

end module bfgs_min_mod