! This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.

program beta_and_derivatives

!---------------!
  !
  ! ==== Global Data ====
  !
  USE beta_module
  !
  ! ==== Local Data ====
  !

real(dp) :: a
real(dp) :: b
real(dp) :: U,t
real(dp) :: x,beta,dbetadt,dbetadU
real(dp),parameter :: tol=1.0e-9_dp
integer,parameter :: maxiter=300
integer :: iter
integer :: ier

read(*,*) U,t

a = 1.0_dp
b = 2.0_dp

U_tmp=U
t_tmp=t
call secant(f,a,b,x,tol,maxiter,iter,ier)
beta = x

a = 1.0_dp
b = 2.0_dp
U_tmp=U+0.000001_dp
t_tmp=t
call secant(f,a,b,x,tol,maxiter,iter,ier)
dbetadU=(x-beta)/0.000001_dp
U_tmp=U

open (97, file='beta_dbetadU.dat',access='sequential')
 write(97,15) beta,dbetadU
close(97)

    15 format(f15.10,f15.10,f15.10)
end program
