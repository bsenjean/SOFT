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

program SOET_DMET
!---------------!
  !
  ! ==== Local Data ====
  !
 implicit none
 integer,parameter :: dp=kind(1.d0), OUT=99
 integer :: i,j,L,Nelec,nocc,imp
 real(dp) :: t,U
 DOUBLE PRECISION, allocatable :: D_SOFT(:,:),D_SOFT_diag(:,:), &
                                  D(:,:), D_F(:,:),D_E(:,:),D_unocc_F(:,:), D_occ_E(:,:),D_unocc_E(:,:), &
                                  U_SVD(:,:),V_SVD(:,:),S(:),D_occ_F(:,:), &
                                  sqrt_n0(:,:),n0(:,:), sqrt_un_n0(:,:), inv_sqrt_un_n0(:,:), &
                                  V_F(:,:),V_E(:,:),Socc(:,:),tildeC_F(:,:),C_B(:,:), &
                                  D_provisoir(:,:),eig(:),P(:,:),Pbar(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         Allocation and Initialization           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 open(10,file='mat.dat',status='old')
 read(10,*)L, Nelec, U, t, imp ! #sites, #electrons, t, U, #impurities
 allocate(D_F(imp,L))
 allocate(D_E(L-imp,L))
 allocate(D_SOFT(L,L))
 allocate(D_SOFT_diag(L,L))
 allocate(D(L,L))
 allocate(eig(L))

 open ( UNIT = OUT, FILE = 'Schmidt_decomposition.out', ACCESS = 'SEQUENTIAL' )

 write(OUT,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
 write(OUT,*)'% This program performs the Schmidt Decomposition %'
 write(OUT,*)'% on the Slater determinant obtained by KS-SOFT.  %'
 write(OUT,*)'% It returns the embedded problem to be solved    %'
 write(OUT,*)'% either by DMRG for 2 or more than 2 impurities, %'
 write(OUT,*)'% or analytically for a single impurity site.     %'
 write(OUT,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
 write(OUT,*)
 write(OUT,*)'Number sites                :', L
 write(OUT,*)'Number of impurity sites    :', imp
 write(OUT,*)'Number of electrons         :', Nelec
 write(OUT,*)'On-site Coulomb repulsion U :', U
 write(OUT,*)'Hopping parameter t         :', t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             READ MATRIX ELEMENTS                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 write(OUT,*)
 write(OUT,*)'Step 1. Read the matrix representation of :'
 write(OUT,*)'        - the 1RDM of the KS-SOFT system, D_soft (size L x L)'
!!! SOFT density matrix
 do i=1,L
    read(10,*)(D_SOFT(i,j),j=1,L)
 enddo
 close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Diagonalize the non-interacting Hamiltonian   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(OUT,*)
 write(OUT,*)'Step 2. Diagonalize D_soft. Determine the transformation matrix D.'
 write(OUT,*)'        D_diag = D^T D_soft D'

! Initialize the transformation D to be set in diasym.
! The minus sign to sort from the biggest to the lowest one.
 D=-D_SOFT
 call diasym(D,eig,L)
 D_SOFT_diag=matmul(matmul(transpose(D),D_SOFT),D)

! Determine nocc
 nocc = 0
 do i=1,L
    if (D_SOFT_diag(i,i) > 0.99999) then
       nocc = nocc + 1 
    endif
 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Split transformation matrix into         !
!         fragment (F) + environment (E)          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(OUT,*)
 write(OUT,*)'Step 3. Split the transformation matrix D (L x L) into:'
 write(OUT,*)'         -  a fragment part D_F (n_imp x L),'
 write(OUT,*)'         -  an environment part D_E (L-n_imp x L).'
 do i=1,imp
  do j=1,L
    D_F(i,j)=D(i,j)
  enddo
 enddo

 do i=1,L-imp
  do j=1,L
     D_E(i,j)=D(i+imp,j)
  enddo
 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Split transformation matrix into         !
!            occupied and unoccupied              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(OUT,*)
 write(OUT,*)'Step 4. Split D_F and D_E into occupied and unoccupied parts:'
 write(OUT,*)'        -  D_occ_F (n_imp x n_occ),'
 write(OUT,*)'        -  D_unocc_F (n_imp x L-n_occ),' 
 write(OUT,*)'        -  D_occ_E(L-n_imp x n_occ),'
 write(OUT,*)'        -  D_unocc_E(L-n_imp x L-n_occ).'
 allocate(D_occ_F(imp,nocc))
 allocate(D_unocc_F(imp,L-nocc))
 allocate(D_occ_E(L-imp,nocc))
 allocate(D_unocc_E(L-imp,L-nocc))

 do i=1,imp
  do j=1,nocc
     D_occ_F(i,j)=D_F(i,j)
  enddo
 enddo

 do i=1,imp
  do j=1,L-nocc
     D_unocc_F(i,j)=D_F(i,j+nocc)
  enddo
 enddo

 do i=1,L-imp
  do j=1,nocc
     D_occ_E(i,j)=D_E(i,j)
  enddo
 enddo

 do i=1,L-imp
  do j=1,L-nocc
     D_unocc_E(i,j)=D_E(i,j+nocc)
  enddo
 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          Singular Value Decomposition           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(OUT,*)
 write(OUT,*)'Step 5. Perform a Singular Value Decomposition of D_occ_F.'
 write(OUT,*)'        D_occ_F = U_SVD sqrt_n0 V_SVD with size:'
 write(OUT,*)'        - U_SVD (n_imp x n_imp),'
 write(OUT,*)'        - sqrt_n0 (n_imp x n_imp),'
 write(OUT,*)'        - V_SVD (n_occ x n_occ).'
 allocate(D_provisoir(imp,nocc))
 allocate(S(nocc)) ! temporary matrix for having sqrt_n0 in SVD calculation.
 allocate(sqrt_n0(imp,imp))
 D_provisoir=D_occ_F ! because the input is destroyed in DGSEDD, and I don't want to destroy D_occ_F
 allocate(U_SVD(imp,imp))
 allocate(V_SVD(nocc,nocc))

 call SVD(D_provisoir,U_SVD,S,sqrt_n0,V_SVD,imp,nocc)
 
 deallocate(D_provisoir)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Decomposition of V_SVD into          !
!           fragment (F) + environment (E)        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(OUT,*)
 write(OUT,*)'Step 6. Decompose V_SVD into:'
 write(OUT,*)'        - a fragment part V_F (n_occ x n_imp),'
 write(OUT,*)'        - an environment part V_E (n_occ x n_occ-n_imp).'
 allocate(V_F(nocc,imp))
 allocate(V_E(nocc,nocc-imp))

 do i=1,nocc
  do j=1,imp
     V_F(i,j)=V_SVD(i,j)
  enddo
  do j=1,nocc-imp
     V_E(i,j)=V_SVD(i,j+imp)
  enddo
 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Compute the overlap matrix Socc          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(OUT,*)
 write(OUT,*)'Step 7. Compute the overlap matrix, S_occ (n_occ x n_occ)'
 write(OUT,*)'        S_occ = V_F n0 V_F^T, '
 write(OUT,*)'        where n0 = sqrt_n0^T sqrt_n0.'

 allocate(Socc(nocc,nocc))
 allocate(n0(imp,imp))
 n0=matmul(transpose(sqrt_n0),sqrt_n0)
 Socc=matmul(matmul(V_F,n0),transpose(V_F))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Compute the Projector                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(OUT,*)
 write(OUT,*)'Step 8. Compute the Projector P (L x 2n_imp):'
 write(OUT,*)'        P = Pbar x Block_matrix[[U_SVD^T,0],[0,U_SVD^T]],'
 write(OUT,*)'        where Pbar is the Projector in the basis of natural orbital:'
 write(OUT,*)'        Pbar = Block_matrix[[U_SVD,0],[0,C_B]],'
 write(OUT,*)'        C_B = D_occ_E V_F (sqrt(1 - n0))^-1.'
 allocate(tildeC_F(imp,imp))
 allocate(sqrt_un_n0(imp,imp))
 allocate(inv_sqrt_un_n0(imp,imp))
 allocate(C_B(L-imp,imp))
 allocate(Pbar(L,2*imp))
 allocate(P(L,2*imp))
 allocate(D_provisoir(2*imp,2*imp))

 tildeC_F=matmul(D_occ_F,V_F)

 sqrt_un_n0=0
 do i=1,imp
    sqrt_un_n0(i,i)=sqrt(1-n0(i,i))
 enddo

 call inverse(sqrt_un_n0,inv_sqrt_un_n0,imp)

 C_B = matmul(matmul(D_occ_E,V_F),inv_sqrt_un_n0)

 Pbar=0
 do i=1,imp
  do j=1,imp
    Pbar(i,j)=U_SVD(i,j)
  enddo
 enddo

 do i=imp+1,L
  do j=imp+1,2*imp
    Pbar(i,j)=C_B(i-imp,j-imp)
  enddo
 enddo

 D_provisoir=0
 do i=1,imp
  do j=1,imp
     D_provisoir(i,j)=U_SVD(j,i)
     D_provisoir(i+imp,j+imp)=U_SVD(j,i)
  enddo
 enddo

 P=matmul(Pbar,D_provisoir)

 write(OUT,*)
 write(OUT,*)' Projector P :'
 do i=1,L
   write(OUT,10)(P(i,j),j=1,2*imp)
 enddo
 write(OUT,*)
 write(OUT,*)'End of file.'
 open (98,FILE='projector.dat',ACCESS ='SEQUENTIAL')
 do i=1,L
   write(98,10)(P(i,j),j=1,2*imp)
 enddo
 close(98)


close(OUT)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      Format                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    10 format(500f25.15)

 end program SOET_DMET
