! ######################### !
! ### Global parameters ### !
! ######################### !
! last edited 10/21/2022 S. Malpthak
module parameters
implicit none

! pi and the imaginary unit
real*8, parameter     :: pi = dacos(-1.d0)
complex*16, parameter :: Iu = (0.d0,1.d0)

! propagator parameters
integer :: Ntime, MaxIt              ! number of timesteps, max. # of iterations
real*8  :: TimeStep, ConvTol         ! Timestep, Convergrence tolerance for memory kernel

integer :: Nnuc, Nel                 !Sub-system dimensionalities

real*8                  :: Eps, GammaC                    ! parameters for electronic sub-system
real*8,allocatable      :: OmegaN(:), CouplingN(:)        ! Freqs. and couplings strenghts 
                                                          ! and spring constant
                                                          ! = Mass*Omega**2 for bath modes
real*8                  :: Beta                           ! Beta, inverse temp.

real*8                  :: OmegaC, Eta, Xi                ! Ohmic bath parameters. Eta = pi*hbar*xi/2

complex*16, allocatable :: Fk(:,:,:), Fdotk(:,:,:)        !Auxiliary kernels 
real*8,allocatable      :: L0(:,:)                        !bath-averaged Louivillian
end module

! ########################### !
! ### Read the input file ### !
! ########################### !
subroutine input
use parameters
implicit none

integer                 :: i, t
real*8                  :: temp
real*8,allocatable      :: rInp(:,:,:), iInp(:,:,:)
character*75            :: infostr

open(555,file='input',status='old')
read(555,'(a75)') infostr
read(555,*) Nnuc, Nel

!Ndof = Nnuc + Nel

read(555,'(a75)') infostr
read(555,*) Eps, GammaC

read(555,'(a75)') infostr
read(555,*) Beta

read(555,'(a75)') infostr
read(555,*) OmegaC, Xi
Eta = pi*Xi/2.d0

read(555,'(a75)') infostr
read(555,*) TimeStep, Ntime, ConvTol, MaxIt

close(555)

!Timestep = Timestep*4.d0
!Ntime    = int(dble(Ntime)/4.d0)

allocate(OmegaN(Nnuc),CouplingN(Nnuc))
allocate(Fk(Nel**2,Nel**2,0:Ntime),Fdotk(Nel**2,Nel**2,0:Ntime))
allocate(L0(Nel**2,Nel**2))
allocate(rInp(Nel**2,Nel**2,0:Ntime),iInp(Nel**2,Nel**2,0:Ntime))

OmegaN          = 0.d0
CouplingN       = 0.d0
L0              = 0.d0
Fk              = 0.d0
Fdotk           = 0.d0
rInp            = 0.d0
iInp            = 0.d0

!Generate freqs and coupling of the bath
do i = 1, Nnuc
   OmegaN(i)    =-OmegaC*log((dfloat(i)-0.5d0)/Nnuc)
   CouplingN(i) = OmegaN(i)*dsqrt(2.d0*Eta*1.d0*OmegaC/pi/Nnuc) ! m = 1.d0
enddo

! Read in elements of F kernal
open(1005,file='fort.1005',status='old')
open(2005,file='fort.2005',status='old')

do t = 0, Ntime
   read(1005,*) temp, rInp(1,1,t), rInp(1,2,t), rInp(1,3,t), rInp(1,4,t), &
                            rInp(2,1,t), rInp(2,2,t), rInp(2,3,t), rInp(2,4,t), &
                            rInp(3,1,t), rInp(3,2,t), rInp(3,3,t), rInp(3,4,t), &
                            rInp(4,1,t), rInp(4,2,t), rInp(4,3,t), rInp(4,4,t)
 !  if (t.ne.Ntime) then
 !     read(1005,*)  
 !     read(1005,*) 
 !     read(1005,*) 
!   end if
   read(2005,*) temp, iInp(1,1,t), iInp(1,2,t), iInp(1,3,t), iInp(1,4,t), &
                            iInp(2,1,t), iInp(2,2,t), iInp(2,3,t), iInp(2,4,t), &
                            iInp(3,1,t), iInp(3,2,t), iInp(3,3,t), iInp(3,4,t), &
                            iInp(4,1,t), iInp(4,2,t), iInp(4,3,t), iInp(4,4,t)
!   if (t.ne.Ntime) then
!     read(2005,*) 
!     read(2005,*) 
!     read(2005,*) 
!   end if
end do
close(1005)
close(2005)
! Construct F kernal from real and imaginary parts
Fk = rInp + Iu*iInp


rInp = 0.d0
iInp = 0.d0
! Read in elements of Fdot kernal
open(3005,file='fort.3005',status='old')
open(4005,file='fort.4005',status='old')

do t = 0, Ntime
   read(3005,*) temp, rInp(1,1,t), rInp(1,2,t), rInp(1,3,t), rInp(1,4,t), &
                            rInp(2,1,t), rInp(2,2,t), rInp(2,3,t), rInp(2,4,t), &
                            rInp(3,1,t), rInp(3,2,t), rInp(3,3,t), rInp(3,4,t), &
                            rInp(4,1,t), rInp(4,2,t), rInp(4,3,t), rInp(4,4,t)

!   if (t.ne.Ntime) then
!     read(3005,*) 
!     read(3005,*) 
!     read(3005,*) 
!   end if
   read(4005,*) temp, iInp(1,1,t), iInp(1,2,t), iInp(1,3,t), iInp(1,4,t), &
                            iInp(2,1,t), iInp(2,2,t), iInp(2,3,t), iInp(2,4,t), &
                            iInp(3,1,t), iInp(3,2,t), iInp(3,3,t), iInp(3,4,t), &
                            iInp(4,1,t), iInp(4,2,t), iInp(4,3,t), iInp(4,4,t)
!   if (t.ne.Ntime) then
!     read(4005,*)  
!     read(4005,*) 
!     read(4005,*) 
!   end if

end do
close(3005)
close(4005)
! Construct Fdot kernal from real and imaginary parts
! Fdot is actually i*Fdot
Fdotk = rInp + Iu*iInp

do t = 0, Ntime
  Fk(:,:,t) = dsin(5.d0*dble(t)*TimeStep)*dexp(-dble(t)*Timestep/4.d0)
end do

!Generate L0 operator
call GenerateL0()



end subroutine input
!################################################################
subroutine GenerateL0()
use parameters, only : L0, Nel, Eps, GammaC, OmegaN, Beta, Nnuc

integer :: i, j 
real*8  :: H0(Nel,Nel)

H0 = 0.d0
!{ Define the average hamiltonian
H0(1,2) = GammaC
H0(2,1) = GammaC

do i = 1, Nnuc
   H0(1,1) = H0(1,1) + OmegaN(i)/2.d0/dtanh(Beta*OmegaN(i)/2.d0) 
end do
H0(2,2) = H0(1,1)

H0(1,1) = H0(1,1) + Eps
H0(2,2) = H0(2,2) - Eps
!}
!{ Define the average Louivillian
L0(1,1) =  0.d0
L0(1,2) = -H0(2,1)
L0(1,3) =  H0(1,2)
L0(1,4) =  0.d0
L0(2,1) =  L0(1,2)
L0(2,2) =  H0(1,1)-H0(2,2)
L0(2,3) =  0.d0
L0(2,4) =  H0(1,2)
L0(3,1) =  L0(1,3)
L0(3,2) =  L0(2,3)
L0(3,3) = -L0(2,2)
L0(3,4) = -H0(2,1)
L0(4,1) =  L0(1,4)
L0(4,2) =  L0(2,4)
L0(4,3) =  L0(3,4)
L0(4,4) =  0.d0
!}
end subroutine GenerateL0


