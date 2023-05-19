! function created by ChatGPT-4 needed to implement steepest descent
module optim_util_mod
implicit none
private
public :: dp, numerical_gradient
integer, parameter :: dp = kind(1.0d0)
contains
    function numerical_gradient(funk, x, eps) result(grad)
        interface
            function funk(x)
                import :: dp
                real(kind=dp), dimension(:), intent(in) :: x
                real(kind=dp) :: funk
            end function funk
        end interface
        real(kind=dp), dimension(:), intent(in) :: x
        real(kind=dp), intent(in), optional :: eps ! width used in finite difference
        real(kind=dp), dimension(size(x)) :: grad
        real(kind=dp), dimension(size(x)) :: x_plus, x_minus
        real(kind=dp) :: delta, eps_
        integer :: i
        if (present(eps)) then 
           eps_ = eps
        else
           eps_ = 1.0e-8
        end if
        x_plus = x
        x_minus = x

        do i = 1, size(x)
            delta = eps_ * (1.0_dp + abs(x(i)))
            x_plus(i) = x(i) + delta
            x_minus(i) = x(i) - delta
            grad(i) = (funk(x_plus) - funk(x_minus)) / (2.0_dp * delta)
            x_plus(i) = x(i)
            x_minus(i) = x(i)
        end do
    end function numerical_gradient
end module optim_util_mod