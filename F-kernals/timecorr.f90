! ################################################### !
! ### this program computes the electronic state  ### !
! ### populations as well as the final nuclear    ### !
! ### momentum distribution of a 1D particle in   ### ! 
! ### a 2-level system with the DF implementation ### !
! ### of MQC-IVR.                                 ### !
! ################################################### !
! last edited 10/21/2022 S.V.Malpathak
program CF
use parameters
implicit none

include 'mpif.h'

integer             :: i, j, k, l, m, u, v
integer             :: flagE
integer             :: countE, brokenE
integer*8           :: brokenEF, trajused
integer             :: ierr, myrank, nprocs, istart, iend, imc
 
real*8, allocatable :: PScoord(:,:), qpfwd(:,:,:)
real*8, allocatable :: Initialq(:), Initialp(:)

complex*16, allocatable :: probs(:), AwigF(:,:), AwigFdot(:,:), Bwig(:,:,:)
complex*16, allocatable :: F(:,:,:,:,:), Fdot(:,:,:,:,:), FT(:,:,:,:,:), FdotT(:,:,:,:,:)
complex*16, allocatable :: Fker(:,:,:,:,:,:), Fdotker(:,:,:,:,:,:)
call input

allocate(PScoord(2,Ndof),qpfwd(2,Ndof,0:Ntime))
allocate(Initialq(Ndof),Initialp(Ndof))
allocate(AwigF(Nel,Nel),AwigFdot(Nel,Nel),Bwig(Nel,Nel,0:Ntime))
allocate(F(Nel,Nel,Nel,Nel,0:Ntime),Fdot(Nel,Nel,Nel,Nel,0:Ntime))
allocate(FT(Nel,Nel,Nel,Nel,0:Ntime),FdotT(Nel,Nel,Nel,Nel,0:Ntime))
allocate(Fker(Ncheck,Nel,Nel,Nel,Nel,0:Ntime),Fdotker(Ncheck,Nel,Nel,Nel,Nel,0:Ntime))

! normalization constant for LSCsample
!normC = 1.d0/(pi)**Nel

F       = 0.d0
Fdot    = 0.d0
FT      = 0.d0
FdotT   = 0.d0
Fker    = 0.d0
Fdotker = 0.d0
AwigF   = 0.d0
AwigFdot= 0.d0
Bwig    = 0.d0

! accumulate broken trajectory info
brokenEF = 0    ! broken energy
trajused = 0    ! total trajectories used
! set up parallelization
call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nprocs,ierr)
call mpi_comm_rank(mpi_comm_world,myrank,ierr)

if (myrank.eq.0) then
   !write percentage of broken trajectories
   open(111,file='Traj_Info.out',status='unknown')
end if

do m = 1, Ncheck ! write correlation function at each m

 countE  = 0     ! count broken energy
 F       = 0.d0
 Fdot    = 0.d0
 call para_range(1,NumberMCsteps,nprocs,myrank,istart,iend)
 call mpi_barrier(mpi_comm_world,ierr)

 call init_random_seed(myrank)

 do imc = istart, iend  ! loop over trajectories
   ! sample scheme for LSC
   call LSCsample(PScoord)
   flagE  = 0

   ! Initial conditions with LSCsample sampling
   Initialq = PScoord(1,:)
   Initialp = PScoord(2,:)

   ! propagate first trajectory
   call PropagateFwd(Initialq,Initialp,qpfwd,flagE)
   ! check for conservation of ...
   if (flagE.eq.1) then                 ! energy
      countE = countE + 1
      goto 222
   endif

   flagE = 0 
  
   ! The leftover parts of rho_wig/LSCsample function
   call AhatF(PSCoord(1,Nnuc+1:Ndof),PSCoord(2,Nnuc+1:Ndof),AwigF)
   call AhatFdot(PScoord(1,Nnuc+1:Ndof),PScoord(2,Nnuc+1:Ndof),PScoord(1,1:Nnuc),PScoord(2,1:Nnuc),AwigFdot)

   do i = 0, Ntime
       call Bhat(qpfwd(1,Nnuc+1:Ndof,i),qpfwd(2,Nnuc+1:Ndof,i),qpfwd(1,1:Nnuc,i),qpfwd(2,1:Nnuc,i),Bwig(:,:,i))
       do j = 1, Nel
          do k = 1, Nel
             do u = 1, Nel
                do v = 1, Nel
                   F(j,k,u,v,i)     = F(j,k,u,v,i)    + AwigF(u,v)*Bwig(j,k,i)
                   Fdot(j,k,u,v,i)  = Fdot(j,k,u,v,i) + AwigFdot(u,v)*Bwig(j,k,i)
                   !Fdot is actually i*Fdot
                end do
             end do
          end do
       end do
   end do
        
  222 continue

 enddo
 call mpi_barrier(mpi_comm_world,ierr)
! Gather F and Fdot kernals
 call mpi_reduce(F,FT,(Ntime+1)*(Nel**4),mpi_double_complex,mpi_sum,0,mpi_comm_world,ierr)
 call mpi_reduce(Fdot,FdotT,(Ntime+1)*(Nel**4),mpi_double_complex,mpi_sum,0,mpi_comm_world,ierr)
! gather conservation counters
 call mpi_reduce(countE,brokenE,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)

 if (myrank.eq.0) then
   Fker(m,:,:,:,:,:)    = FT(:,:,:,:,:)
   Fdotker(m,:,:,:,:,:) = FdotT(:,:,:,:,:)

   brokenEF   = brokenEF + brokenE
   trajused   = trajused + (NumberMCsteps - brokenE)

   FT    = 0.d0
   FdotT = 0.d0
   do j = 1, Ncheck
     FT    = FT + Fker(j,:,:,:,:,:)
     FdotT = FdotT + Fdotker(j,:,:,:,:,:)
   end do
   FT    = FT/dble(trajused)
   FdotT = FdotT/dble(trajused)
  
   do i = 0, Ntime
      write(1000+m,'(E15.6)',advance='no'), dble(i*TimeStep)
      do j = 1, Nel
         do k = 1, Nel
            do u = 1, Nel
               do v = 1, Nel
                  if (j+k+u+v.eq.4*Nel) goto 501 
                  write(1000+m,'(E15.6)',advance='no'), dble(FT(j,k,u,v,i))
               end do
            end do
         end do
      end do 
      501 continue 
      write(1000+m,'(E15.6)'), dble(FT(2,2,2,2,i))
 
      write(2000+m,'(E15.6)',advance='no'), dble(i*TimeStep)
      do j = 1, Nel
         do k = 1, Nel
            do u = 1, Nel
               do v = 1, Nel
                  if (j+k+u+v.eq.4*Nel) goto 502 
                  write(2000+m,'(E15.6)',advance='no'), dimag(FT(j,k,u,v,i))
               end do
            end do
         end do
      end do 
      502 continue 
      write(2000+m,'(E15.6)'), dimag(FT(2,2,2,2,i))

      write(3000+m,'(E15.6)',advance='no'), dble(i*TimeStep)
      do j = 1, Nel
         do k = 1, Nel
            do u = 1, Nel
               do v = 1, Nel
                  if (j+k+u+v.eq.4*Nel) goto 503
                  write(3000+m,'(E15.6)',advance='no'), dble(FdotT(j,k,u,v,i))
               end do
            end do
         end do
      end do 
      503 continue 
      write(3000+m,'(E15.6)'), dble(FdotT(2,2,2,2,i))
 
      write(4000+m,'(E15.6)',advance='no'), dble(i*TimeStep)
      do j = 1, Nel
         do k = 1, Nel
            do u = 1, Nel
               do v = 1, Nel
                  if (j+k+u+v.eq.4*Nel) goto 504 
                  write(4000+m,'(E15.6)',advance='no'), dimag(FdotT(j,k,u,v,i))
               end do
            end do
         end do
      end do 
      504 continue 
      write(4000+m,'(E15.6)'), dimag(FdotT(2,2,2,2,i))

    enddo
   
   ! write broken trajectory information
   write(111,*) m
   write(111,*) 'Broken energy     :', brokenEF
   write(111,*) 'Total traj used   :', trajused
   write(111,*) ''
 endif

enddo

if (myrank.eq.0) then
  close(111)
end if

call mpi_finalize(ierr)

stop

end program CF
