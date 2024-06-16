! ################################################### !
! ### this program computes the memory kernal K   ### !
! ### from the auxiliary kernals F and Fdot       ### !
! ### by solving a Volerra equation               ### ! 
! ### of the second kind                          ### !
! ################################################### !
! last edited 12/15/2022 S.V.Malpathak
program MemoryKernel
use parameters
implicit none

!include 'mpif.h'

integer             :: i, j, k, l, flag, itcount

complex*16, allocatable :: Kk(:,:,:), Gk(:,:,:), Kold(:,:,:)
complex*16, allocatable :: intg(:,:,:)
call input

allocate(Kk(Nel**2,Nel**2,0:Ntime),Gk(Nel**2,Nel**2,0:Ntime))
allocate(Kold(Nel**2,Nel**2,0:Ntime),intg(Nel**2,Nel**2,0:Ntime))

Kk      = 0.d0
Gk      = 0.d0
Kold    = 0.d0
intg    = 0.d0
flag    = 1
itcount = 0

!Initialize 

do i = 0, Ntime
  Gk(:,:,i) = Fdotk(:,:,i) - matmul(Fk(:,:,i),L0)
end do
Kold = Gk
!print*, "G0", Gk(:,:,0)
!print*, "G1000", Gk(:,:,10)

do while (flag.eq.1)
  
  do i = 0, Ntime
    intg = 0.d0
    do j = 0, i
      ! Get the integrand
      intg(:,:,j) = Iu*matmul(Fk(:,:,i-j),Kold(:,:,j))
    end do
!      if(i.eq.0.or.i.eq.10) then
!        print*, "i", i, "int0", intg(:,:,0), "inti", intg(:,:,i)
!      end if
    !integrate
    call TrapRule(0,i,intg(:,:,0:i),Timestep,Kk(:,:,i)) 
!    if(i.eq.0.or.i.eq.10) then
!      print*, "after int", "i", i, "K12", Kk(1,2,i), "K22", KK(2,2,i)
!    end if 
  end do
  !Add in the g-term
  Kk = Kk + Gk

  !Check convergence
  if (maxval(cdabs(Kk-Kold)).gt.ConvTol) goto 101
  ! If all satisfy tolerance, leave loop
  flag = 0
  101 continue
  itcount = itcount + 1 
  if (itcount.ge.MaxIt) then
    flag = 0
    print*, "Max. iterations reached", itcount
    print*, "Max.diff. found in this iteration", maxval(cdabs(Kk-Kold))
  end if
  ! Guess for next iteration
  Kold  = KK
  ! Count iteration
end do

!print memory kernal
open(1111,file="KReal.out",status='unknown')
open(2222,file="KImag.out",status='unknown')

do i = 0, Ntime
   write(1111,'(E15.6)',advance='no'), dble(i*TimeStep)
   write(2222,'(E15.6)',advance='no'), dble(i*TimeStep)
   do j = 1, Nel**2
      do k = 1, Nel**2
          if (j+k.eq.2*Nel**2) goto 501 
               write(1111,'(E15.6)',advance='no'), dble(Kk(j,k,i))
               write(2222,'(E15.6)',advance='no'), dimag(Kk(j,k,i))
      end do
   end do 
   501 continue 
   write(1111,'(E15.6)'), dble(Kk(Nel**2,Nel**2,i))
   write(2222,'(E15.6)'), dimag(Kk(Nel**2,Nel**2,i))
end do

close(1111)
close(2222)

print*, "# of iterations = ", itcount

end program MemoryKernel


!#######################################################
!       Calculate 1D integral with trapezoid rule      !
!       from xmin to xmax
!#######################################################
subroutine TrapRule(xmin,xmax,fx,dx,intgrl)
use parameters, only : Nel
implicit none

integer, intent(in)     :: xmin, xmax
real*8                  :: dx
complex*16,intent(in)   :: fx(Nel**2,Nel**2,xmin:xmax)
integer                 :: i
complex*16,intent(out)  :: intgrl(Nel**2,Nel**2)

intgrl = 0.d0

if (xmin.eq.xmax) then
   intgrl = 0.d0 
else
   do i = xmin+1, xmax-1
      intgrl = intgrl + fx(:,:,i)
   end do
   intgrl = intgrl + 0.5d0*(fx(:,:,xmin)+fx(:,:,xmax))
   intgrl = intgrl*dx
end if

end subroutine TrapRule



