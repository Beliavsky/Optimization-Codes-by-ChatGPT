program xminimize
    ! try various optimization methods on a battery of test functions
    use adam_mod, only: adam_min, momentum_optimizer, rmsprop, adagrad, adadelta
    use nelder_mead_mod, only: nelder_mead
    use conjugate_gradient_min_mod, only: cg_min
    use bfgs_min_mod, only: bfgs_min
    use gradient_descent_min_mod, only : gradient_descent_min, dp
    use luus_jaakola_min_mod, only: luus_jaakola_min
    use simulated_annealing_mod, only: simulated_annealing
    use hooke_jeeves_mod, only: hooke_jeeves
    implicit none
    integer, parameter :: iter_max = 10**6, ndim=10, ipow_step_max=3, &
       ifunk_try(*) = [1, 2, 3, 4] ! determines which functions to optimize
    real(kind=dp), dimension(ndim) :: xx, grad, min_point, true_min_point
    real(kind=dp) :: min_value, step_size, r, alpha, tinit, tfinal, t1, t2, tol_bfgs, & 
                     tol_cg, tol_nm, true_min_value
    real(kind=dp), parameter :: xh = 1.0d-6, funk_scale = 1.0d0, pi = 3.141592653589793238462643_dp
    integer :: ipow_step, niter, nfunc_eval, niter_lj, ifunk, jfunk, nfunk_try
!   set constants that determine which optimization methods will be tried
    logical, parameter :: do_lj = .true., do_bfgs = .true., do_cg = .true., &
                          do_nm = .true., do_adam = .true., do_mo = .true., &
                          do_rmsprop = .true., do_adagrad = .true., do_hj = .true., &
                          do_adadelta = .true., do_simann = .true., print_true = .true.
    character (len=20), parameter :: func_names(*) = [character(len=20) :: "quadratic", "Rosenbrock", & ! 1:2
   "Styblinski-Tang", "Ackley", "Rastrigin", "Sphere", "Booth", "Beale", "Three-hump camel", &
   "Easom", "Goldstein-Price", "Levy"]
    call cpu_time(tinit)
    nfunk_try = size(ifunk_try)
do_ifunk: do jfunk=1, nfunk_try ! loop over functions to optimize
    ifunk = ifunk_try(jfunk)
    print "(/,a)", "function: " // trim(func_names(ifunk))
    print "('#dimensions = ',i0, /)", ndim
    if (do_lj) then
       xx = 0.0_dp
       r = 1.0_dp
       niter_lj = 10**7
       alpha = 0.98_dp
       call cpu_time(t1)
       call luus_jaakola_min(funk, xx, r, niter_lj, alpha, min_point, min_value)
       call cpu_time(t2)
       print "(/,a)", "Luus-Jaakola method"
       print *, "r: "            , r
       print *, "Minimum point: ", min_point
       print *, "Minimum value: ", min_value
       print *, "time elapsed: ", t2-t1
    end if
    if (do_hj) then
       xx = 0.0_dp
       block
       real(kind=dp), parameter :: epsilon = 1.0d-8
       real(kind=dp) :: step
       call cpu_time(t1)
       step = 0.0001_dp
       call hooke_jeeves(funk, xx, step, epsilon, iter_max, min_point, min_value)
       call cpu_time(t2)
       print "(/,a)", "Hooke-Jeeves method"
       print *, "Minimum point: ", min_point
       print *, "Minimum value: ", min_value
       print *, "time elapsed: ", t2-t1          
       end block
    end if
    if (do_simann) then
       xx = 0.0_dp
       block
       real(kind=dp), parameter :: alpha = 0.9_dp, T0 = 1.0_dp, T_min = 1.0d-3
       call cpu_time(t1)
       call simulated_annealing(funk, xx, T0, T_min, alpha, iter_max, min_point, min_value)
       call cpu_time(t2)
       print "(/,a)", "simulated annealing"
       print *, "Minimum point: ", min_point
       print *, "Minimum value: ", min_value
       print *, "time elapsed: ", t2-t1          
       end block
    end if
    if (do_adam) then
       xx = 0.0_dp
       block
       real(kind=dp), parameter :: alpha = 0.001_dp, beta1 = 0.9_dp, &
          beta2 = 0.999_dp, epsilon = 1.0d-8
       call cpu_time(t1)
       call adam_min(funk, xx, alpha, beta1, beta2, epsilon, iter_max, min_point, min_value)
       call cpu_time(t2)
       print "(/,a)", "Adam method"
       print *, "Minimum point: ", min_point
       print *, "Minimum value: ", min_value
       print *, "time elapsed: ", t2-t1          
       end block
    end if
    if (do_mo) then
       xx = 0.0_dp
       block
       real(kind=dp), parameter :: alpha = 0.01_dp, beta = 0.9_dp
       print "(/,a)", "Momentum Optimization method"
       call cpu_time(t1)
       call momentum_optimizer(funk, xx, alpha, beta, iter_max, min_point, min_value)
       call cpu_time(t2)
       print *, "Minimum point: ", min_point
       print *, "Minimum value: ", min_value
       print *, "time elapsed: ", t2-t1          
       end block
    end if
    if (do_rmsprop) then
       xx = 0.0_dp
       block
       real(kind=dp), parameter :: alpha = 0.001_dp, beta = 0.9_dp, epsilon = 1.0d-8
       print "(/,a)", "RMSProp method"
       call cpu_time(t1)
       call rmsprop(funk, xx, alpha, beta, epsilon, iter_max, min_point, min_value)
       call cpu_time(t2)
       print *, "Minimum point: ", min_point
       print *, "Minimum value: ", min_value
       print *, "time elapsed: ", t2-t1          
       end block
    end if
    if (do_adagrad) then
       xx = 0.0_dp
       block
       real(kind=dp), parameter :: alpha = 0.01_dp, epsilon = 1.0d-8
       print "(/,a)", "Adagrad method"
       call cpu_time(t1)
       call adagrad(funk, xx, alpha, epsilon, iter_max, min_point, min_value)
       call cpu_time(t2)
       print *, "Minimum point: ", min_point
       print *, "Minimum value: ", min_value
       print *, "time elapsed: ", t2-t1          
       end block
    end if
    if (do_adadelta) then
       xx = 0.0_dp
       block
       real(kind=dp), parameter :: rho = 0.95_dp, epsilon = 1.0d-6
       print "(/,a)", "Adadelta method"
       call cpu_time(t1)
       call adadelta(funk, xx, rho, epsilon, iter_max, min_point, min_value)
       call cpu_time(t2)
       print *, "Minimum point: ", min_point
       print *, "Minimum value: ", min_value
       print *, "time elapsed: ", t2-t1          
       end block
    end if
    if (do_bfgs) then
       xx = 0.0_dp
       tol_bfgs = 1.0d-10
       call cpu_time(t1)
       call bfgs_min(funk, xx, iter_max, tol_bfgs, min_point, min_value)
       call cpu_time(t2)
       print "(/,a)", "BFGS method"
       print *, "Minimum point: ", min_point
       print *, "Minimum value: ", min_value
       print *, "time elapsed: ", t2-t1
    end if
    if (do_cg) then
       xx = 0.0_dp
       tol_cg = 1.0d-10
       call cpu_time(t1)
       call cg_min(funk, xx, iter_max, tol_cg, min_point, min_value)
       call cpu_time(t2)
       print "(/,a)", "conjugate gradient method"
       print *, "Minimum point: ", min_point
       print *, "Minimum value: ", min_value
       print *, "time elapsed: ", t2-t1
    end if
    if (do_nm) then
       xx = 0.0_dp
       tol_nm = 1.0d-10
       call cpu_time(t1)
       call nelder_mead(funk, xx, tol_nm, iter_max, min_point, min_value)
       call cpu_time(t2)
       print "(/,a)", "Nelder-Mead"
       print *, "Minimum point: ", min_point
       print *, "Minimum value: ", min_value
       print *, "time elapsed: ", t2-t1
    end if
    print "(/,a)", "steepest descent"
    do ipow_step=0, ipow_step_max
       ! Set initial guess and convergence parameter
       step_size = 10.0**(-ipow_step)
       xx = 0.0_dp
       ! Call gradient_descent_min subroutine
       call cpu_time(t1)
       call gradient_descent_min(funk, xh, step_size, iter_max, xx, min_value, &
            grad, niter, nfunc_eval)
       call cpu_time(t2)
       print *
       print *, "Minimum point: ", xx
       print *, "Minimum value: ", min_value
       print *, "L2 norm of gradient: ",sqrt(sum(grad**2))
       print *, "step size: ",step_size
       print *, "#iterations: ", niter
       print *, "#func eval: ", nfunc_eval
       print *, "time elapsed: ", t2-t1
    end do
    if (print_true) then
       call funk_min(true_min_point, true_min_value)
       print "(/,a)", "True values"
       print *, "Minimum point: ", true_min_point
       print *, "Minimum value: ", true_min_value
    end if
    call cpu_time(tfinal)
end do do_ifunk
print "(/,a,f10.4)", "time elapsed (s): ", tfinal-tinit

contains

! Prompt given to ChatGPT-4:
! Fill in the following subroutine with the known locations and values of function minima 
! for the functions defined in funk.
! 
!     subroutine funk_min(ifunk, xmin, value_min)
!     integer, intent(in) :: ifunk ! function in funk(x) for which minimum found
!     real(kind=dp), intent(out) :: xmin(:)   ! location of minimum
!     real(kind=dp), intent(out) :: value_min ! function value at minimum
!     real(kind=dp), parameter   :: bad_value = -999.0_dp ! value to set xmin and value_min if they are unknown
!     end subroutine funk_min

subroutine funk_min(xmin, value_min)
    real(kind=dp), intent(out) :: xmin(:)   ! location of minimum
    real(kind=dp), intent(out) :: value_min ! function value at minimum
    real(kind=dp), parameter   :: bad_value = -999.0_dp ! value to set xmin and value_min if they are unknown
    integer :: i, n
    n = size(xmin)
    select case (ifunk)
        case (1)
            value_min = 0.0_dp
            xmin = [(real(i, kind=dp), i = 1, n)]
        case (2) ! Rosenbrock function
            value_min = 0.0_dp
            xmin = [(1.0_dp, i = 1, n)]
        case (3) ! Styblinski-Tang Function
            value_min = -39.16599_dp * n
            xmin = [(-2.903534_dp, i = 1, n)]
        case (4) ! Ackley function
            value_min = 0.0_dp
            xmin = [(0.0_dp, i = 1, n)]
        case (5) ! Rastrigin function
            value_min = 0.0_dp
            xmin = [(0.0_dp, i = 1, n)]
        case (6) ! Sphere function
            value_min = 0.0_dp
            xmin = [(0.0_dp, i = 1, n)]
        case (7) ! Booth function
            value_min = 0.0_dp
            xmin(1:2) = [1.0_dp, 3.0_dp]
        case (8) ! Beale function
            value_min = 0.0_dp
            xmin(1:2) = [3.0_dp, 0.5_dp]
        case (9) ! Three-hump camel function
            value_min = 0.0_dp
            xmin(1:2) = [0.0_dp, 0.0_dp]
        case (10) ! Easom function
            value_min = -1.0_dp
            xmin(1:2) = [pi, pi]
        case (11) ! Goldstein-Price function
            value_min = 3.0_dp
            xmin(1:2) = [0.0_dp, -1.0_dp]
        case (12) ! Levy function
            value_min = 0.0_dp
            xmin(1:2) = [(1.0_dp, i = 1, n)]
        case default
            value_min = bad_value
            xmin = [(bad_value, i = 1, n)]
    end select
end subroutine funk_min

    function funk(x)
        ! Prompt to ChatGPT-4 was "Add more test functions to funk, especially those that have 
        ! been studied in the literature on numerical optimization and which have known minima."
        ! test functions for unconstrained optimization
        ! https://en.wikipedia.org/wiki/Test_functions_for_optimization has been consulted
        real(kind=dp), dimension(:), intent(in) :: x
        real(kind=dp) :: funk
        integer :: i
        select case (ifunk)
           case (1)   
              funk = funk_scale * sum([((x(i)-i)**2, i=1,size(x))])
           case (2) ! Rosenbrock function
              funk = 0.0_dp
              do i=1,size(x)-1
                 funk = funk + 100*(x(i+1)-x(i))**2 + (1-x(i))**2
              end do
           case (3) ! Styblinski-Tang Function
              funk = 0.0_dp
              do i=1,size(x)
                 funk = funk + x(i)**4 - 16*x(i)**2 + 5*x(i)
              end do
              funk = funk/2
           case (4) ! Ackley function
              funk = -20*exp(-0.2*sqrt(sum(x**2)/size(x))) - exp(sum(cos(2*pi*x))/size(x)) + 20 + exp(1.0_dp)
           case (5) ! Rastrigin function
               funk = 10*size(x) + sum(x**2 - 10*cos(2*pi*x))
           case (6) ! Sphere function
              funk = sum((x-1)**2)
        end select
    end function funk

end program xminimize
