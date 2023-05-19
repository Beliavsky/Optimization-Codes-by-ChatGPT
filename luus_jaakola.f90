! ChatGPT-4 told to "write a subroutine that implements the Luus-Jaakola 
! direct search optimization method."
module luus_jaakola_min_mod
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    private
    public :: luus_jaakola_min, dp

    contains

subroutine luus_jaakola_min(funk, xx, r, n_iter, alpha, min_point, min_value)
    interface
        function funk(x)
            import :: dp
            real(kind=dp), dimension(:), intent(in) :: x
            real(kind=dp) :: funk
        end function funk
    end interface
    real(kind=dp), dimension(:), intent(in) :: xx
    real(kind=dp), intent(in out) :: r
    integer, intent(in) :: n_iter
    real(kind=dp), intent(in) :: alpha
    real(kind=dp), dimension(:), intent(out) :: min_point
    real(kind=dp), intent(out) :: min_value
    real(kind=dp), dimension(:), allocatable :: x_new
    real(kind=dp) :: f_new, r_min_
    integer :: n, i
    r_min_ = 1.0d-6
    ! Initialize min_point and min_value
    min_point = xx
    min_value = funk(min_point)

    ! Allocate memory for x_new
    n = size(xx)
    allocate(x_new(n))

    ! Luus-Jaakola optimization loop
    do i = 1, n_iter
        ! Generate a new candidate point
        call random_number(x_new)
        x_new = min_point + r * (2 * x_new - 1)

        ! Evaluate the function at the new point
        f_new = funk(x_new)
      
        ! Update min_point and min_value if the new point is better
        if (f_new < min_value) then
            min_point = x_new
            min_value = f_new
            ! Decrease the search radius
            r = r * alpha
            if (r < r_min_) return
            ! print*,r, min_value, min_point ! debug
        end if
    end do
end subroutine luus_jaakola_min

end module luus_jaakola_min_mod
