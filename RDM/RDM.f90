! ################################################### !
! ### this program computes the memory kernal K   ### !
! ### from the auxiliary kernals F and Fdot       ### !
! ### by solving a Volerra equation               ### ! 
! ### of the second kind                          ### !
! ################################################### !
! last edited 12/15/2022 S.V.Malpathak
program RDM
use parameters
implicit none

!include 'mpif.h'

integer                 :: i, j, t, flag, itcount, TotTS
real*8                  :: dt, tr
complex*16, allocatable :: rho(:,:), K1(:), K2(:), K3(:)
complex*16, allocatable :: integral(:,:), K4(:)
complex*16, allocatable :: rhotemp(:)
call input

TotTS = int(dble(Ntime)/2.d0) ! # of timesteps for rho
dt    = Timestep*2.d0         ! timestep for rho
 
!Timestep for rho is twice the timestep that K is calculated with 
allocate(rho(Nel**2,0:TotTS),K1(Nel**2),K2(Nel**2))
allocate(K3(Nel**2),K4(Nel**2))
allocate(integral(Nel**2,0:TotTS),rhotemp(Nel**2))

rho     = 0.d0
rhotemp = 0.d0
K1      = 0.d0
K2      = 0.d0
K3      = 0.d0
K4      = 0.d0
integral= 0.d0
flag    = 1
itcount = 0

!Initialize
rho(1,0) = 1.d0
 
do t = 0, TotTS-1
  call Calc_Int(t,rho(:,0:t),integral(:,t))
 
  !Get K1
  K1 = -Iu*matmul(L0,rho(:,t)) - integral(:,t)
  
  ! Get K2
  rhotemp = rho(:,t) + Timestep*K1
  K2 = -Iu*matmul(L0,rhotemp) - integral(:,t) &
       -Timestep*matmul(Kk(:,:,2*t+1),rho(:,0))
  
  !Get K3
  rhotemp = rho(:,t) + Timestep*K2
  K3 = -Iu*matmul(L0,rhotemp) - integral(:,t) &
       -Timestep*matmul(KK(:,:,2*t+1),rho(:,0))
  
  !Get K4
  rhotemp = rho(:,t) + dt*K3
  K4 = -Iu*matmul(L0,rhotemp) - integral(:,t) &
       -dt*matmul(KK(:,:,2*(t+1)),rho(:,0))
  
  rho(:,t+1) = rho(:,t) + (K1+2.d0*K2+2.d0*K3+K4)*dt/6.d0
  
  tr = 0.d0
  !Specific for Nel = 1
  tr = dble(rho(1,t+1) + rho(4,t+1))
   
  if(dabs(tr-1.d0).gt.ConvTol) then
     print*, "Trace deviating from tolerance", "Time:", t*dt, "Trace:", tr
  end if
end do

!print memory kernal
open(1111,file="RhoReal.out",status='unknown')
open(2222,file="RhoImag.out",status='unknown')

do i = 0, TotTS
   write(1111,'(E15.6)',advance='no'), dble(i*dt)
   write(2222,'(E15.6)',advance='no'), dble(i*dt)
   do j = 1, (Nel**2)-1
         write(1111,'(E15.6)',advance='no'), dble(rho(j,i))
         write(2222,'(E15.6)',advance='no'), dimag(rho(j,i))
   end do 
   write(1111,'(E15.6)'), dble(rho(Nel**2,i))
   write(2222,'(E15.6)'), dimag(rho(Nel**2,i))
end do

close(1111)
close(2222)


end program RDM

!#######################################################
!       Calculate 1D integral with trapezoid rule      !
!       from xmin to xmax
!#######################################################
subroutine TrapRule(xmin,xmax,fx,dx,intgrl)
use parameters, only : Nel
implicit none

integer, intent(in)     :: xmin, xmax
real*8                  :: dx
complex*16,intent(in)   :: fx(Nel**2,xmin:xmax)
integer                 :: i
complex*16,intent(out)  :: intgrl(Nel**2)

intgrl = 0.d0

if (xmin.eq.xmax) then
   intgrl = 0.d0 
else
   do i = xmin+1, xmax-1
      intgrl = intgrl + fx(:,i)
   end do
   intgrl = intgrl + 0.5d0*(fx(:,xmin)+fx(:,xmax))
   intgrl = intgrl*dx
end if

end subroutine TrapRule
!######################################################!
!               Calculate rhodot given rho             !       
!######################################################!
subroutine Calc_Int(t,dens,integral)
use parameters, only : Timestep, Kk, Nel
implicit none

integer, intent(in)   :: t
complex*16,intent(in) :: dens(Nel**2,0:t)
integer               :: i
complex*16            :: intg(Nel**2,0:t)
complex*16,intent(out):: integral(Nel*2)

intg = 0.d0
do i = 0, t
   intg(:,i) = matmul(Kk(:,:,2*i),dens(:,t-i))
end do
!Get intergal
call TrapRule(0,t,intg(:,0:t),2.d0*Timestep,integral) 

end subroutine Calc_Int

