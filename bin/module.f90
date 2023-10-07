module beta_module

integer,parameter :: dp=kind(1.d0)
real, parameter :: Pi=3.14159265358979323846264338327950288419716939937510
real(dp) :: t_tmp,U_tmp

contains
function integrand(x) result(res)

real(dp),intent(in) :: x
real(dp) :: res
real(dp) :: num,den

num = bessel_j0(x)*bessel_j1(x)
den = x*(1.0_dp + exp(x*U_tmp/2.0_dp/t_tmp))

res = num/den

end function

function f(x) result(res)

real(dp),intent(in) :: x
real(dp) :: res
real(dp) :: a,b

integer,parameter :: limit=1000
integer,parameter :: lenw=limit*4
real(dp) :: abserr
real(dp),parameter :: epsabs=0.0_dp
real(dp),parameter :: epsrel=0.00001_dp
integer :: ier
integer :: iwork(limit)
integer,parameter :: inf=1
integer :: last
integer :: neval
real(dp) :: work(lenw)

! double, quadrature de Gauss pour integrer de 0 a inf. (inf = 1 
! equivaut a infini
call dqagi(integrand,0.0_dp,inf,epsabs,epsrel,a,abserr,neval, &
           ier,limit,lenw,last,iwork,work)

b = -sin(Pi/x)*2.0_dp*t_tmp*x/Pi

res = b + 4.0_dp*t_tmp*a

end function

end module beta_module
