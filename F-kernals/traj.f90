! ############################################################
! ############################################################
!       this subroutine stores a classical trajectory
!       in the 'coord' array, the corresponding
!       monodromy matrix in the 'Monodromy' array, and
!       the action integral in the 'Sf' array, given
!       initial phase space points 'Initialq' and 'Initialp'.
! #############################################################
! #############################################################
! parts of original code supplied by T.J.H.Hele

subroutine PropagateFwd(Initialq,Initialp,coord,Sf,flagE)
use parameters, only: TimeStep, Ntime, EnergyTolerance, Ndof, Nnuc, Nel
implicit none

real*8, intent(in)     :: Initialq(Ndof), Initialp(Ndof)
real*8, intent(out)    :: coord(2,Ndof,0:Ntime), Sf(0:Ntime)
integer, intent(inout) :: flagE
integer                :: i, j, k
real*8                 :: InitialEnergy, dtp(3), Energy
real*8                 :: U, dU(Nnuc), V(Nel,Nel), dV(Nel,Nel,Nnuc)
real*8                 :: hdtm, qdot(Nel), H2
real*8, allocatable    :: q(:), p(:)

  allocate(q(Ndof),p(Ndof))
  flagE  = 0
  ! initial conditions
  q = Initialq
  p = Initialp

  ! total initial energy
  call Hamiltonian2(q(1:Nnuc),q(Nnuc+1:Ndof),p(Nnuc+1:Ndof),H2)
  InitialEnergy = 0.5d0*dot_product(p(1:Nnuc),p(1:Nnuc)) + H2

  ! initial conditions
  Sf(0)        = 0.d0
  coord(1,:,0) = Initialq
  coord(2,:,0) = Initialp
  hdtm         = 0.5d0*timestep ! mass_j = 1

  do i = 1, Ntime

     ! evolve nuclear position and classical action
     ! under H1 for time dt/2
     q(1:Nnuc)  = q(1:Nnuc) + p(1:Nnuc)*hdtm
     Sf(i)      = Sf(i-1) + 0.5d0*dot_product(p(1:Nnuc),p(1:Nnuc))*hdtm

     ! get V(R) matrix elements with new nuclear position
     call potbits(q(1:Nnuc),U,dU,V,dV)

     ! evolve nuclear momentum and mapping variables under H2 for time dt
     call H2evolve(q(1:Nnuc),p(1:Nnuc),q(Nnuc+1:Ndof),p(Nnuc+1:Ndof),U,dU,V,dV)

     qdot    = matmul(V,p(Nnuc+1:Ndof))
     call Hamiltonian2(q(1:Nnuc),q(Nnuc+1:Ndof),p(Nnuc+1:Ndof),H2)

     ! evolve classical action under H2 for time dt
     Sf(i)  = Sf(i) + (dot_product(p(Nnuc+1:Ndof),qdot) - H2)*timestep

     ! evolve nuclear position under H1 for time dt/2
     q(1:Nnuc) = q(1:Nnuc) + p(1:Nnuc)*hdtm
     ! evolve classical action under H1 for time dt/2
     Sf(i) = Sf(i) + 0.5d0*dot_product(p(1:Nnuc),p(1:Nnuc))*hdtm

     ! energy conservation check
     call Hamiltonian2(q(1:Nnuc),q(Nnuc+1:Ndof),p(Nnuc+1:Ndof),H2)
     Energy = 0.5d0*dot_product(p(1:Nnuc),p(1:Nnuc)) + H2

     ! store trajectory 
      coord(1,:,i)       = q
      coord(2,:,i)       = p

     ! check energy conservation
     if (dabs(1.d0-Energy/InitialEnergy).ge.EnergyTolerance) then
        flagE = 1
        goto 112
     endif
     ! proceed to drop trajectory of not conserved

  end do

112 continue

end subroutine PropagateFwd

! ############################################################
! ############################################################
! ### this subroutine evolves nuclear momentum 'PP' under
! ### H2 fort time 'timestep'
! #############################################################
! #############################################################
subroutine H2evolve(RR,PP,qs,ps,U,dU,V,dV)
use parameters, only : timestep, kN, Nnuc, Nel
implicit none

integer                         :: i
real*8, intent(inout)           :: PP(Nnuc), qs(Nel), ps(Nel)
real*8, intent(in)              :: U, dU(Nnuc), V(Nel,Nel), dV(Nel,Nel,Nnuc), RR(Nnuc)
real*8, dimension(Nel,Nel)      :: Smat, Cmat, Dmat
real*8, dimension(Nel,Nel,Nnuc) :: Emat, Fmat, Gammat, Ximat, Wmat
real*8, dimension(2,2)          :: dGammat, dXimat
real*8, dimension(Nel)          :: lam, qc, pc, co, si
real*8                          :: dt, sinterm, costerm, delp(Nnuc), lam12

  dt = timestep

  ! get eigenvalues 'lam' and eigenvectors 'Smat' of V(R), 
  ! as well as their derivatives 'dlam' and 'dSmat'
  call eigensys(V,dV,lam,Smat)

  do i = 1, Nel
     co(i) = dcos( lam(i)*dt)
     si(i) = dsin(-lam(i)*dt) 
  enddo

  ! copy the original phase space variables
  qc = qs
  pc = ps

  ! matrix multiplication subroutine
  call Threm(Smat,co,transpose(Smat),Cmat)
  call Threm(Smat,si,transpose(Smat),Dmat)

  ! propagate mapping variables
  qs = matmul(Cmat,qc) - matmul(Dmat,pc)
  ps = matmul(Cmat,pc) + matmul(Dmat,qc)


  ! rotate V'(R) in adiabatic basis
  do i = 1, Nnuc
    Wmat(:,:,i)  = matmul(transpose(Smat),matmul(dV(:,:,i),Smat))
  end do
  lam12        = lam(1) - lam(2)
  Ximat(:,:,:) = 0.d0
  Gammat       = Wmat*dt

  if (lam12.ne.0.d0) then
     sinterm       =  sin(lam12*dt)
     Gammat(1,2,:) =  sinterm*Wmat(1,2,:)/lam12
     Gammat(2,1,:) =  Gammat(1,2,:)
     costerm       =  1.d0-cos(lam12*dt)
     Ximat(1,2,:)  =  Wmat(1,2,:)*costerm/lam12
     Ximat(2,1,:)  = -Ximat(1,2,:) ! skew-symmetric
  endif

  do i = 1, Nnuc
     Emat(:,:,i) = matmul(matmul(Smat,Gammat(:,:,i)),transpose(Smat))
     Fmat(:,:,i) = matmul(matmul(Smat,Ximat(:,:,i)),transpose(Smat))
  end do

  ! evolve nuclear momentum
  do i = 1, Nnuc
    delp(i) = dot_product(qc,matmul(Emat(:,:,i),qc)) + dot_product(pc,matmul(Emat(:,:,i),pc)) & 
            - 2.d0*dot_product(qc,matmul(Fmat(:,:,i),pc))
    PP(i) = PP(i) - dU(i)*dt - 0.5d0*delp(i) + 0.5d0*(dV(1,1,i) + dV(2,2,i))*dt
  end do
 
  return
end subroutine H2evolve

! ######################################## !
! ### compute elements of V as well    ### !
! ### the derivatives w.r.t. nuclear R ### ! 
! ######################################## !
subroutine potbits(R,U,dU,V,dV)
use parameters, only   : Nnuc, Nel, kN, CouplingN, Eps, GammaC
implicit none

real*8, intent(in)  :: R(Nnuc)
real*8, intent(out) :: U, dU(Nnuc), V(Nel,Nel), dV(Nel,Nel,Nnuc)
integer             :: i

U       = 0.d0
dU      = 0.d0
V(1,1)  = Eps
V(1,2)  = GammaC
V(2,1)  = GammaC
V(2,2)  =-Eps
dV      = 0.d0

do i = 1, Nnuc
   U         = U + 0.5d0*kN(i)*R(i)**2
   dU(i)     = kN(i)*R(i) 
   V(1,1)    = V(1,1) - CouplingN(i)*R(i) 
   V(2,2)    = V(2,2) + CouplingN(i)*R(i)
   dV(1,1,i) = -CouplingN(i)
   dV(2,2,i) =  CouplingN(i)
end do

return
end subroutine potbits
! ######################################## !
! ### compute Hamiltonian              ### !
! ###   H2                             ### !
! ######################################## !
subroutine Hamiltonian2(Rnuc,x,p,H2)
use parameters, only : Nnuc, Nel
implicit none

real*8, intent(in)  :: Rnuc(Nnuc), x(Nel), p(Nel)
real*8              :: U, dU(Nnuc), V(Nel,Nel), dV(Nel,Nel,Nnuc)
real*8, intent(out) :: H2
integer             :: i

call potbits(Rnuc,U,dU,V,dV)

H2 = U + 0.5d0*(dot_product(p,matmul(V,p)) + dot_product(x,matmul(V,x))) 

do i = 1, Nel
   H2 = H2 - 0.5d0*V(i,i)
end do

end subroutine Hamiltonian2
! ######################################## !
! ### compute the eigenvalues/vectors  ### ! 
! ### of V as well as the derivatives  ### !
! ######################################## !
subroutine eigensys(V,dV,lam,Smat)
use parameters, only : pi, Nel, Nnuc
implicit none

real*8, intent(in)   :: V(Nel,Nel), dV(Nel,Nel,Nnuc)
real*8               :: hp, hm, hpr(Nnuc), hmr(Nnuc), theta, disc
real*8, intent(out)  :: lam(Nel), Smat(Nel,Nel)!, dlam(Nel,Nnuc)


  hp   = 0.5d0*(V(1,1)+V(2,2))
  hm   = 0.5d0*(V(1,1)-V(2,2))
  disc = dsqrt(hm**2 + V(1,2)**2)

  lam(1) = hp + disc ! eigenvalue 1
  lam(2) = hp - disc ! eigenvalue 2

!  hpr = 0.5d0*(dV(1,1,:)+dV(2,2,:))
!  hmr = 0.5d0*(dV(1,1,:)-dV(2,2,:))

!  dlam(1,:) = hpr(:) + (hmr(:)*hm + V(1,2)*dV(1,2,:))/disc ! derivative e value 1
!  dlam(2,:) = hpr(:) - (hmr(:)*hm + V(1,2)*dV(1,2,:))/disc ! derivative e value 2

  if (V(1,2).eq.0.d0) then
     if (V(1,1).ge.V(2,2)) then
        theta = 0.d0
     else 
        theta = pi*0.5d0
     endif
  elseif (V(1,1).eq.V(2,2)) then
     if (V(1,2).gt.0) then
        theta = pi*0.25d0
     else
        theta = pi*0.75d0
     endif
  else

     theta = datan(2*V(1,2)/(V(1,1)-V(2,2)))
     if (theta.le.0) theta = theta + pi
     theta = theta*0.5d0
     if (V(1,2).le.0) theta = theta + pi*0.5d0

  end if

  ! matrix of column eigenvectors
  Smat(1,1)  =  dcos(theta)
  Smat(2,2)  =  Smat(1,1)
  Smat(1,2)  = -dsin(theta)
  Smat(2,1)  = -Smat(1,2)

  return
end subroutine eigensys

! matrix multiplier
subroutine threm(Amat,lam,Bmat,Tmat)
implicit none

integer             :: j, k
real*8, intent(in)  :: Amat(2,2), lam(2), Bmat(2,2)
real*8, intent(out) :: Tmat(2,2)

  do j = 1, 2
     do k = 1, 2
        Tmat(j,k) = Amat(j,1)*lam(1)*Bmat(1,k) + Amat(j,2)*lam(2)*Bmat(2,k)
     enddo
  enddo

  return
end subroutine threm
