# Optimization Codes by ChatGPT
The prompts given to ChatGPT-4 are shown in the codes. Sometimes cleanup was required.

Compile with 

`gfortran optim_util.f90 adam.f90 nelder_mead.f90 conjugate_gradient_min.f90 bfgs_min.f90 luus_jaakola.f90 gradient_descent.f90 xminimize.f90` 

Works with gfortran 12.0.1 20220213 and ifort version 2021.6.0.

When ChatGPT-4 generates code, it often forgets the `import` statement that is needed in the interface below.
```
    interface
        function funk(x)
            import :: dp
            real(kind=dp), dimension(:), intent(in) :: x
            real(kind=dp) :: funk
        end function funk
    end interface`
```
