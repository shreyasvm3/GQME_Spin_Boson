! MPI data division
subroutine para_range(n1,n2,nprocs,irank,ista,iend)
implicit none
integer, intent(in) :: n1,n2,nprocs,irank
integer, intent(out) :: ista,iend
integer :: iwork1,iwork2
iwork1 = (n2 - n1 + 1)/nprocs
iwork2 = mod(n2 - n1 + 1, nprocs)
ista = irank*iwork1 + n1 + min(irank, iwork2)
iend = ista + iwork1 -1
if (iwork2 .gt. irank) then
iend = iend + 1
endif
return
end subroutine para_range

! initiate seed
subroutine init_random_seed(myrank)
integer :: myrank
integer :: i, n, clock
integer, allocatable :: seed(:)

call random_seed(size = n)
allocate(seed(n))

call system_clock(count=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /) + myrank*clock

call random_seed(put = seed)

deallocate(seed)
end subroutine init_random_seed
!#########################################
!      Left over bits of Awig for F kernel
!#########################################
subroutine AhatF(x,p,AF)
use parameters, only : Nel, Iu
implicit none

real*8, intent(in)      :: x(Nel), p(Nel)
complex*16, intent(out) :: AF(Nel,Nel)
integer                 :: i, j
real*8                  :: krdelta

AF = 0.d0

do i = 1, Nel
  do j = 1, Nel 
     AF(i,j) = x(i)*x(j) + p(i)*p(j) + Iu*(x(i)*p(j)-p(i)*x(j)) &
             - 0.5d0*krdelta(i,j)
  end do
end do

AF = 2.d0*AF        

end subroutine AhatF
!#########################################
!      Left over bits of Awig for Fdot kernel
!#########################################
subroutine AhatFdot(x,p,RR,PP,AFdot)
use parameters, only : Nel, Iu, Nnuc, OmegaN, CouplingN, Eps, GammaC, kN, Beta
implicit none

real*8, intent(in)      :: x(Nel), p(Nel), RR(Nnuc), PP(Nnuc)
complex*16, intent(out) :: AFdot(Nel,Nel)
integer                 :: i, j, l
complex*16              :: HRhoWig(Nel,Nel), RhoHWig(Nel,Nel)
real*8                  :: krdelta

AFdot = 0.d0

!Define the wigner transforms of HRho and RhoH with the sampling function
!already divided out
!Specific for spin-boson model. Will change based on H
HRhoWig      = 0.d0
RhoHWig      = 0.d0

HRhoWig(1,2) =  GammaC
HRhoWig(2,1) =  GammaC
RhoHWig(1,2) =  GammaC
RhoHWig(2,1) =  GammaC

!{First caculate terms that are common to all
do i = 1, Nnuc
   HRhoWig(1,1) = HRhoWig(1,1) + ((PP(i)**2 + kN(i)*RR(i)**2)/2.d0  &
                + OmegaN(i)/4.d0*dsinh(Beta*OmegaN(i)))/dcosh(Beta*OmegaN(i)/2.d0)**2
end do
HRhoWig(2,2) = HRhoWig(1,1)
RhoHWig(1,1) = HRhoWig(1,1)
RhoHWig(2,2) = HRhoWig(1,1)
!}

!Now add in the rest of the terms
do i = 1, Nnuc 
   HRhoWig(1,1) = HRhoWig(1,1) + Iu*CouplingN(i)*PP(i)/OmegaN(i)*dtanh(Beta*OmegaN(i)/2.d0) &
                - CouplingN(i)*RR(i)
   HRhoWig(2,2) = HRhoWig(2,2) - Iu*CouplingN(i)*PP(i)/OmegaN(i)*dtanh(Beta*OmegaN(i)/2.d0) &
                + CouplingN(i)*RR(i)
   RhoHWig(1,1) = RhoHWig(1,1) - Iu*CouplingN(i)*PP(i)/OmegaN(i)*dtanh(Beta*OmegaN(i)/2.d0) &
                - CouplingN(i)*RR(i)
   RhoHWig(2,2) = RhoHWig(2,2) + Iu*CouplingN(i)*PP(i)/OmegaN(i)*dtanh(Beta*OmegaN(i)/2.d0) &
                + CouplingN(i)*RR(i)
end do 

HRhoWig(1,1) = HRhoWig(1,1) + Eps
HRhoWig(2,2) = HRhoWig(2,2) - Eps
RhoHWig(1,1) = RhoHWig(1,1) + Eps
RhoHWig(2,2) = RhoHWig(2,2) - Eps
!}

!{Calculate A_wig with sampling function left out.
do i = 1, Nel
   do j = 1, Nel
      do l = 1, Nel
         AFdot(i,j) = AFdot(i,j) &
                    + (x(l)*x(j)+p(l)*p(j)+Iu*(x(l)*p(j)-p(l)*x(j)) - 0.5d0*krdelta(l,j))*HRhoWig(l,i) &
                    - (x(i)*x(l)+p(i)*p(l)+Iu*(x(i)*p(l)-p(i)*x(l)) - 0.5d0*krdelta(i,l))*RhoHWig(j,l) 
      end do
   end do
end do
          
AFdot = 2.d0*AFdot          

end subroutine AhatFdot
!#########################################
!      Bwig for F and Fdot kernels 
!#########################################
subroutine Bhat(x,p,RR,PP,Bwig)
use parameters, only : Nel, Iu, Nnuc
implicit none

real*8, intent(in)      :: x(Nel), p(Nel), RR(Nnuc), PP(Nnuc)
complex*16, intent(out) :: Bwig(Nel,Nel)
integer                 :: j,k,l
real*8                  :: krdelta, H(Nel,Nel), id(Nel,Nel)
real*8                  :: U, dU(Nnuc), V(Nel,Nel), dV(Nel,Nel,Nnuc)

Bwig = 0.d0

call potbits(RR,U,dU,V,dV)

!Construct the Hamiltonian matrix
id = 0.d0
do j = 1, Nel
   id(j,j) = 1.d0
end do
H = (U + dot_product(PP,PP)/2.d0)*id + V

do j = 1, Nel
  do k = 1, Nel 
      do l = 1, Nel
         Bwig(j,k) = Bwig(j,k) &
                   + (x(k)*x(l)+p(k)*p(l)+Iu*(x(k)*p(l)-p(k)*x(l)) - 0.5d0*krdelta(k,l))*H(j,l) &
                   - (x(l)*x(j)+p(l)*p(j)+Iu*(x(l)*p(j)-p(l)*x(j)) - 0.5d0*krdelta(l,j))*H(l,k)
      end do
  end do
end do

Bwig = (2.d0**(Nel+1))*dexp(-dot_product(x,x)-dot_product(p,p))*Bwig        

end subroutine Bhat

!##############################################
!     Kronecker delta function
!##############################################
double precision function krdelta(i,j)

implicit none
integer, intent(in)     :: i, j 

if (i.eq.j) then
   krdelta = 1.d0
elseif (i.ne.j) then
   krdelta = 0.d0
endif

end function krdelta



