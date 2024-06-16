! ######################### !
! ### Global parameters ### !
! ######################### !
! last edited 10/21/2022 S. Malpthak
module parameters
implicit none

! pi and the imaginary unit
real*8, parameter     :: pi = dacos(-1.d0)
complex*16, parameter :: Iu = (0.d0,1.d0)

!
integer :: NumberMCsteps ! number of trajectory pairs
integer :: Ncheck        ! check convergence after using 'NumberMCsteps'
                         ! trajectories, repeat 'Ncheck' times

! propagator parameters
integer :: Ntime                     ! number of timesteps
real*8  :: TimeStep, EnergyTolerance ! timestep, energy tolerance

integer :: Ndof ! full system dimensionality
integer :: Nnuc, Nel

real*8                  :: Eps, GammaC                    ! parameters for electronic sub-system
real*8,allocatable      :: OmegaN(:), CouplingN(:), kN(:) ! Freqs. and couplings strenghts 
                                                          ! and spring constant
                                                          ! = Mass*Omega**2 for bath modes
real*8                  :: Beta                           ! Beta, inverse temp.

real*8                  :: OmegaC, Eta, Xi                ! Ohmic bath parameters. Eta = pi*hbar*xi/2

end module

! ########################### !
! ### Read the input file ### !
! ########################### !
subroutine input
use parameters
implicit none

integer         :: i, j, k
character*75    :: infostr

open(555,file='input',status='old')
read(555,'(a75)') infostr
read(555,*) Nnuc, Nel

Ndof = Nnuc + Nel

allocate(OmegaN(Nnuc),CouplingN(Nnuc),kN(Nnuc))

read(555,'(a75)') infostr
read(555,*) Eps, GammaC

read(555,'(a75)') infostr
read(555,*) Beta

read(555,'(a75)') infostr
read(555,*) OmegaC, Xi

Eta = pi*Xi/2.d0

read(555,'(a75)') infostr
read(555,*) TimeStep, Ntime, EnergyTolerance

read(555,'(a75)') infostr
read(555,*) NumberMCsteps

read(555,'(a75)') infostr
read(555,*) Ncheck

close(555)

do i = 1, Nnuc
   OmegaN(i)    =-OmegaC*log((dfloat(i)-0.5d0)/Nnuc)
   CouplingN(i) = OmegaN(i)*dsqrt(2.d0*Eta*1.d0*OmegaC/pi/Nnuc) ! m = 1.d0
   kN(i)        = OmegaN(i)**2 ! Spring constant for each bath HO
enddo

open(111,file="Frequency.out",status='unknown')
write(111,*) "Discretized bath frequencies and couplings"
write(111,*)
do i = 1, Nnuc
  write(111,*) OmegaN(i), CouplingN(i)
enddo 
close(111)

end subroutine input
