! ##################################################
! ### sample initial conditions with a Gaussian ###
! ### random number generator.                  ###
! ##################################################
subroutine LSCsample(PScoord)
use parameters, only : Ndof, Nnuc, Nel, OmegaN, Beta
implicit none

real*8, intent(out)     ::      PScoord(2,Ndof)
real*8                  ::      gauss, sigmaP, sigmaQ
integer                 ::      i 

! nuclear dof
do i = 1, Nnuc
   sigmaP = dsqrt(OmegaN(i)/2.d0/dtanh(Beta*OmegaN(i)/2.d0))
   sigmaQ = dsqrt(0.5d0/OmegaN(i)/dtanh(Beta*OmegaN(i)/2.d0))
   PScoord(1,i) = gauss(sigmaQ) 
   PScoord(2,i) = gauss(sigmaP) 
end do
! initially excited oscillator
  PScoord(1,Nnuc+1) = gauss(dsqrt(0.5d0))
  PScoord(2,Nnuc+1) = gauss(dsqrt(0.5d0)) 
! initial ground state oscillators
do i = Nnuc + 2, Ndof
  PScoord(1,i) = gauss(dsqrt(0.5d0))    
  PScoord(2,i) = gauss(dsqrt(0.5d0)) 
end do
end subroutine
!##############################################
!     Gaussian random number generator
!##############################################
double precision function gauss(sigma)

implicit none
real*8, intent(in)      :: sigma
real*8                  :: w, v1, v2, l
real*8                  :: s1, s2

w = 2.d0

do
        call random_number(s1)
        call random_number(s2)

        v1 = 2.d0*s1 - 1.d0
        v2 = 2.d0*s2 - 1.d0
        w = v1*v1 + v2*v2

        if (w.lt.1.d0) exit
end do

l     = v1*sqrt(-2.d0*log(w)/(w))
l     = sigma*l
gauss = l

end function gauss
