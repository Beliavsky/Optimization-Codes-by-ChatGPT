# Optimization Codes by ChatGPT
The prompts given to ChatGPT-4 are shown in the codes. Sometimes cleanup was required. This is an experiment in code generation by ChatGPT -- the codes should not be relied upon.

Compile with 

`gfortran optim_util.f90 adam.f90 nelder_mead.f90 conjugate_gradient_min.f90 bfgs_min.f90 luus_jaakola.f90 gradient_descent.f90 xminimize.f90` 

Works with gfortran 12.0.1 20220213 and ifort version 2021.6.0.

## Some ChatGPT-4 quirks

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
It can produce code with `implicit none` and undeclared variables that does not compile. It tends to declare a variable `i` as integer even in a procedure that does not use the variable `i`.

When ChatGPT-4 translated some algorithms from the literature with `g*g'` or equivalently `g*transpose(g)`, it erroneously translated this to `g**2` in Fortran instead of `sum(g**2)`.

ChatGPT-4 sometimes forgets that the `random_number()` intrinsic is a subroutine, not a function, so that you cannot write 
```
x = random_number()
```
ChatGPT-4 sometimes declares as `intent(in)` an argument that needs to be `intent(in out)`, since its initial value in changed in the procedure.
