! ChatGPT-4 told to "Implement the Nelder-Mead Simplex method."
module nelder_mead_mod
    use iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: nelder_mead

contains

    subroutine nelder_mead(funk, x0, tol, max_iter, xmin, fmin)
        ! Nelder-Mead Simplex method for optimization

        interface
            real(kind=dp) function funk(x)
                import :: dp
                real(kind=dp), dimension(:), intent(in) :: x
            end function funk
        end interface

        real(kind=dp), dimension(:), intent(in)  :: x0
        real(kind=dp), intent(in)  :: tol
        integer, intent(in) :: max_iter
        real(kind=dp), dimension(:), intent(out) :: xmin
        real(kind=dp), intent(out) :: fmin

        real(kind=dp), dimension(size(x0)+1,size(x0)) :: simplex
        real(kind=dp), dimension(size(x0)+1) :: fvalues
        real(kind=dp), dimension(size(x0)) :: xnew, xcentroid
        real(kind=dp) :: fr, fe, fc, distance, max_distance
        integer :: i, iter, r_index, e_index, c_index, n
        n = size(x0)

        ! Initialize the simplex
        simplex(1,:) = x0
        fvalues(1) = funk(x0)
        do i = 1, n
            simplex(i+1,:) = x0
            simplex(i+1,i) = simplex(i+1,i) + 1.0_dp
            fvalues(i+1) = funk(simplex(i+1,:))
        end do

        iter = 0
        do while (iter < max_iter)
            ! Sort the simplex and function values
            call sort_simplex(simplex, fvalues)

            ! Check for convergence
            max_distance = 0.0_dp
            do i = 2, n+1
               distance = sqrt(sum((simplex(i,:) - simplex(1,:))**2))
               max_distance = max(max_distance, distance)
            end do
            if (max_distance <= tol) exit

            ! Compute centroid of the n best points
            xcentroid = sum(simplex(1:n,:)) / n

            ! Reflection
            xnew = xcentroid + 1.0_dp * (xcentroid - simplex(n+1,:))
            fr = funk(xnew)
            r_index = n+1

            if (fr < fvalues(n)) then
                ! Expansion
                xnew = xcentroid + 2.0_dp * (xcentroid - simplex(n+1,:))
                fe = funk(xnew)
                e_index = n+1

                if (fe < fvalues(n)) then
                    ! Accept the expansion point
                    simplex(e_index,:) = xnew
                    fvalues(e_index) = fe
                else
                    ! Accept the reflection point
                    simplex(r_index,:) = xnew
                    fvalues(r_index) = fr
                end if
            else
                ! Contraction
                xnew = xcentroid + 0.5_dp * (xcentroid - simplex(n+1,:))
                fc = funk(xnew)
                c_index = n+1

                if (fc < fvalues(n+1)) then
                    ! Accept the contraction point
                    simplex(c_index,:) = xnew
                    fvalues(c_index) = fc
                else

                    ! Shrink the simplex
                    do i = 2, n+1
                        simplex(i,:) = simplex(1,:) + 0.5_dp * (simplex(i,:) - simplex(1,:))
                        fvalues(i) = funk(simplex(i,:))
                    end do
                end if
            end if
            iter = iter + 1
        end do

        ! Output the best point and its function value
        xmin = simplex(1,:)
        fmin = fvalues(1)
    end subroutine nelder_mead

    subroutine sort_simplex(simplex, fvalues)
        real(kind=dp), dimension(:,:) :: simplex
        real(kind=dp), dimension(:) :: fvalues
        integer :: i, j
        do i = 1, size(fvalues)-1
            do j = i+1, size(fvalues)
                if (fvalues(j) < fvalues(i)) then
                    call swap_rows(simplex, i, j)
                    call swap(fvalues(i), fvalues(j))
                end if
            end do
        end do
    end subroutine sort_simplex

    subroutine swap_rows(mat, row1, row2)
        real(kind=dp), dimension(:,:) :: mat
        integer :: row1, row2
        real(kind=dp), dimension(size(mat,2)) :: tmp
        tmp = mat(row1,:)
        mat(row1,:) = mat(row2,:)
        mat(row2,:) = tmp
    end subroutine swap_rows

    subroutine swap(a, b)
        real(kind=dp) :: a, b, tmp
        tmp = a
        a = b
        b = tmp
    end subroutine swap

end module nelder_mead_mod

