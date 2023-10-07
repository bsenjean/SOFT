subroutine secant(f,a,b,x,tol,maxiter,iter,ier)

integer, parameter :: dp=kind(1.d0)
real(dp),external :: f
real(dp),intent(inout) :: a
real(dp),intent(inout) :: b
real(dp),intent(inout) :: x
real(dp),intent(in) :: tol
integer,intent(in) :: maxiter
integer,intent(out) :: iter
integer,intent(out) :: ier

real(dp) :: fa,fb,xzero,fxzero,deltax

fa = f(a)
fb = f(b)

if (fa*fb .gt. 0.0_dp) then
 print*,"ERROR IN SECANT METHOD !!! BAD INTERVAL !!!"
end if

xzero = b - fb*(b - a)/(fb - fa)
fxzero = f(xzero)
iter = 1
deltax = 1.0_dp

do while (deltax.gt.tol.and.iter.lt.maxiter)
 x = xzero - fxzero*(xzero - b)/(fxzero - fb)
 deltax = abs(xzero - x)
 b = xzero
 xzero = x
 fb = fxzero
 fxzero = f(x)
 iter = iter + 1
end do

 if (iter.gt.maxiter) then
  ier = 1
  print*,"ERROR IN SECANT METHOD !!! TOO MANY ITERATIONS !!!"
 else
  ier = 0
 end if

end subroutine
