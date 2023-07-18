#include "dprec.fh"
!
!********************************************************************
!
!  Subroutines to work in Liouville space
!
!  Converted to fortran 90 and interfaced with sqm(AmberTools) by
!  Kirill Velizhanin (kirill@lanl.gov)
!
!********************************************************************
!
!--------------------------------------------------------------------
!
!  Acting by TDHF tetradic operator L on an arbitrary single-particle
!  matrix xi.
!       
!  Notation closely resembles one in S. Tretiak and S. Mukamel, Chem. Rev.,
!  102, 3171-3212 (2002)
!
!--------------------------------------------------------------------
! ! seems that it takes advanthe of 0 block strucutre of guess fector in MO basis
  !
   subroutine site2mo(qm2ds,zz,xi,v1)
   use qm2_davidson_module
   implicit none
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds
   _REAL_ f0,f1
   parameter (f0 = 0.0)
   parameter (f1 = 1.0)
   _REAL_ xi(qm2ds%Nb,qm2ds%Nb),v1(qm2ds%Nrpa)
   _REAL_ zz(qm2ds%Nb,qm2ds%Nb)
   
   !write (6,*) 'BEFORE SITE2MO'
   !dgemm C = alpha*A*B + beta*C

  ! call print_1d(v1, 'v1')
   call dgemm ('T','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,f1,xi,qm2ds%Nb, &
      qm2ds%vhf,qm2ds%Nb,f0,zz,qm2ds%Nb)

  ! write (6,*) 'after 1 DGEMM'
 !  call print_2d(zz, 'zz')

   
   call dgemm ('T','N',qm2ds%Nh,qm2ds%Np,qm2ds%Nb,f1, &
      qm2ds%vhf(:,qm2ds%Np+1:qm2ds%Nb),qm2ds%Nb,zz,qm2ds%Nb,f0,v1,qm2ds%Nh)

 !  write (6,*) 'after 2 DGEMM'
  ! call print_1d(v1, 'v1')

   call dgemm ('T','N',qm2ds%Nh,qm2ds%Np,qm2ds%Nb,f1, &
      zz(1,qm2ds%Np+1),qm2ds%Nb,qm2ds%vhf,qm2ds%Nb,f0, &
      v1(qm2ds%Ncis+1),qm2ds%Nh)      

  ! write (6,*) 'after 3 DGEMM'
 !  call print_1d(v1, 'v1')
   return
   end subroutine
!


   subroutine mo2site (qm2ds,v1,xi,zz)

   use qm2_davidson_module
   implicit none
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds
   _REAL_ f0,f1
   integer i,j
   parameter (f0 = 0.0) ! do not add beta*C, f0=eta=0
   parameter (f1 = 1.0) ! alpha, no scaling 

   _REAL_ xi(qm2ds%Nb,qm2ds%Nb),v1(qm2ds%Nrpa)
   _REAL_ zz(qm2ds%Nb,qm2ds%Nb)

   real(8) :: v1_2d(qm2ds%Nh,qm2ds%Np)

   ! v1 = u1 = eta from davidson
   ! xi = qm2ds%xi
   ! zz = qm2ds%eta
   ! call mo2site(qm2ds,u1,qm2ds%xi,qm2ds%eta)


! write (6,*) 'BEFORE MO2SITE'

!  call print_2d(qm2ds%vhf, 'qm2ds%vhf')
!  call print_2d( qm2ds%vhf(:,qm2ds%Np+1:qm2ds%Nb), ' qm2ds%vhf(:,qm2ds%Np+1:qm2ds%Nb)')
!     write (6,*) 'v1', v1
!     write (6,*) 'v1 size', size(v1)
!     write (6,*) 'qm2ds%xi', qm2ds%xi
!       write (6,*) 'qm2ds%eta or zz', qm2ds%eta
   
   
!    ! call print_2d(xi, 'xi ')
!    call print_2d(zz, 'zz ')
! call print_2d(reshape(v1, [qm2ds%Nh, qm2ds%Np]), ' v1 reshape')



   !dgemm C = alpha*A*B + beta*C
   call dgemm ('T', 'T', &! TRANS A B
               qm2ds%Np, qm2ds%Nb, qm2ds%Nh,& ! M N K
               f1, & ! alpha
 !              reshape(v1, [qm2ds%Nh, qm2ds%Np]), &   ! A
               v1(1), &   ! A    ! seems like v1(1) or v1 are the same, probably the index at which start array -> 2d transformation
 !              v1_2d, &   ! A
               qm2ds%Nh, &! LDA - dimension
               qm2ds%vhf(:,qm2ds%Np+1:qm2ds%Nb), & ! B
               qm2ds%Nb, &! LDB - dimension
               f0,       &! beta
               zz,       &! C, output
               qm2ds%Nb)  !LDC - dimension

   ! write (6,*) 'after 1 DGEMM'

   ! ! call print_2d(qm2ds%vhf, 'qm2ds%vhf')
   ! ! call print_2d(xi, 'xi ')
   ! call print_2d(zz, 'zz ')

   
 !  call print_2d(v1, 'v1 ')
 !---------------------------------------------------------------------
   ! write (6,*) 'v1(qm2ds%Ncis+1)',v1(qm2ds%Ncis+1)
   ! write (6,*) 'shape v1(qm2ds%Ncis+1)', shape(v1(qm2ds%Ncis+1))
   call dgemm ('N','T', &
               qm2ds%Nh,qm2ds%Nb,qm2ds%Np, &
               f1, &
               v1(qm2ds%Ncis+1), & ! A
               
               qm2ds%Nh, &
               qm2ds%vhf, &        ! B
               qm2ds%Nb, &
               f0, &
               zz(qm2ds%Np+1,1), &
               qm2ds%Nb)

!    write (6,*) 'after 2 DGEMM'            
! !  !  call print_2d(v1, 'v1 ')
!     call print_2d(zz, 'zz ')
   
!     write (6,*) 'after 3 DGEMM'

!    call print_2d(xi, 'xi ')
!    call print_2d(zz, 'zz ')

               
   call dgemm ('N','N', &
               qm2ds%Nb,qm2ds%Nb,qm2ds%Nb, &
               f1,qm2ds%vhf, &
                qm2ds%Nb, &
                zz,qm2ds%Nb, &
                f0, &
                xi, &
                qm2ds%Nb)

   ! call print_2d(qm2ds%vhf, 'qm2ds%vhf after mo2site')
!  call print_2d(xi, 'xi after mo2site')
   ! call print_2d(zz, 'zz after mo2site')


   return
   end subroutine
!
!********************************************************************
!
   subroutine site2moph (Nb,Np,Nh,M4,uu,xi,v1,zz)
   implicit none
   integer Nb,Np,Nh,M4
   _REAL_ uu(Nb*Nb),xi(Nb,Nb),zz(Nb,Nb),v1(2*M4)
   _REAL_ f0,f1
   parameter (f0 = 0.0)
   parameter (f1 = 1.0)
   call dgemm ('T','N',Nb,Nb,Nb,f1,xi,Nb,uu,Nb,f0,zz,Nb)
   call dgemm ('T','N',Nh,Np,Nb,f1,uu(Nb*Np+1),Nb,zz,Nb,f0,v1,Nh)
   call dgemm ('T','N',Nh,Np,Nb,f1,zz(1,Np+1),Nb,uu,Nb,f0,v1(M4+1),Nh)      
   return
   end subroutine
        
   subroutine mo2siteph (Nb,Np,Nh,M4,uu,v1,xi,zz)
   implicit none
   integer Nb,Np,Nh,M4
   _REAL_ uu(Nb*Nb),xi(Nb,Nb),zz(Nb,Nb),v1(2*M4)
   _REAL_ f0,f1
   call dgemm ('T','T',Np,Nb,Nh,f1,v1(1),Nh,uu(Nb*Np+1),Nb,f0,zz,Nb)
   call dgemm ('N','T',Nh,Nb,Np,f1,v1(M4+1),Nh,uu,Nb,f0,zz(Np+1,1),Nb)
   call dgemm ('N','N',Nb,Nb,Nb,f1,uu,Nb,zz,Nb,f0,xi,Nb)
   return
   end subroutine
!
!********************************************************************
!
   subroutine site2mof (Nb,uu,xi,eta,zz)
   implicit none
   integer Nb
   _REAL_ uu(Nb,Nb),xi(Nb,Nb),zz(Nb,Nb),eta(Nb,Nb)
   _REAL_ f0,f1
   parameter (f0 = 0.0)
   parameter (f1 = 1.0)
   call dgemm ('N','N',Nb,Nb,Nb,f1,xi,Nb,uu,Nb,f0,zz,Nb)
   call dgemm ('T','N',Nb,Nb,Nb,f1,uu,Nb,zz,Nb,f0,eta,Nb)
   return
   entry mo2sitef (Nb,uu,xi,eta,zz)
   call dgemm ('N','T',Nb,Nb,Nb,f1,xi,Nb,uu,Nb,f0,zz,Nb)
   call dgemm ('N','N',Nb,Nb,Nb,f1,uu,Nb,zz,Nb,f0,eta,Nb)
   return
   end
!
!********************************************************************
!
!********************************************************************
!
   subroutine Vxi(qm2_params,qmmm_mpi,qm2_struct,qm2ds,qmmm_struct, &
                                                            xi,eta)
   use qm2_davidson_module
   use qmmm_struct_module, only : qmmm_struct_type
   use qmmm_module,only:qm2_structure, qmmm_mpi_structure
   use qm2_params_module,  only : qm2_params_type

   implicit none
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qm2_structure),intent(inout) :: qm2_struct
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds
   _REAL_ xi(qm2ds%Nb,qm2ds%Nb),eta(qm2ds%Nb,qm2ds%Nb)
   _REAL_ coef
   integer i,j,l
   qm2ds%xis(:)=0.d0
   qm2ds%etas(:)=0.d0
   coef = 1.000000;

   write (6,*) '============== Vxi ====================='

   ! call print_2d(xi, 'xi ')
   ! call print_2d(eta, 'eta ')

! symmetric part:
   l=0
   do i=1,qm2ds%Nb
      ! write (6,*) 'i', i
    do j=1,i
      ! write (6,*) 'j', j
      !  write (6,*) 'l', l
       l=l+1
       qm2ds%xis(l)=coef*0.5*(xi(i,j)+xi(j,i))
    end do
   end do

   ! call print_1d(qm2ds%xis, 'qm2ds%xis before Vxi pack')
   ! call print_1d(qm2ds%etas, 'qm2ds%etas before Vxi pack')

   

   call Vxi_pack(qm2_params,qmmm_mpi,qm2_struct,qm2ds,qmmm_struct,qm2ds%xis,qm2ds%etas)
  ! call print_1d(qm2ds%xis, 'qm2ds%xis after Vxi pack')
  ! call print_1d(qm2ds%etas, 'qm2ds%etas after Vxi pack')

    l=0 ! probably triangular -> 2d array
    do i=1,qm2ds%Nb
      do j=1,i-1
       !  write (6,*) 'j', j
         l = l + 1
      !   write (6,*) 'l', l
         eta(i,j)=eta(i,j)+qm2ds%etas(l)
         eta(j,i)=eta(j,i)+qm2ds%etas(l)
      end do
      ! write (6,*) 'l2', l
      l=l+1
      eta(i,i)=qm2ds%etas(l)
   end do
   call print_2d(eta, 'G SYM')

!  antisymmetric part:
   l=0
   do i=1,qm2ds%Nb
    do j = 1,i
      !  write (6,*) 'j', j
      l = l + 1
      qm2ds%xis(l)=coef*0.5*(xi(i,j)-xi(j,i))
    end do
   end do


  ! call print_1d(qm2ds%etas, ' antisym qm2ds%etas before Vxi packA')

   call Vxi_packA(qm2_params,qmmm_mpi,qm2_struct,qm2ds,qmmm_struct,qm2ds%xis,qm2ds%etas)
 !  call print_1d(qm2ds%etas, ' antisym qm2ds%etas AFTER Vxi packA')
 !  call print_2d(eta, 'eta symmetric after Vxi')
   l=0
   do i=1,qm2ds%Nb
      do j = 1,i-1
         l = l + 1
         eta(i,j) = eta(i,j) + qm2ds%etas(l)
         eta(j,i) = eta(j,i) - qm2ds%etas(l)
      end do
      l = l + 1
   end do
   call print_2d(eta, 'G SYM + ANTISYM')
   !call print_2d(eta, 'eta ANTIsymmetric after 1 Vxi')

   return





   entry Vxi_symm (xi,eta)   ! eta = Vxi
!  pack:

   l = 0
   do i = 1,qm2ds%Nb
    do j = 1,i
       l = l + 1
       qm2ds%xis(l) = xi(j,i)
    end do
   end do

!  multiply:
   call Vxi_pack(qm2_params,qmmm_mpi,qm2_struct,qm2ds,qmmm_struct,qm2ds%xis,qm2ds%etas)
!  unpack:
   l = 0
   do i = 1,qm2ds%Nb
    do j = 1,i
       l = l + 1
       eta(i,j) = qm2ds%etas(l)
       eta(j,i) = qm2ds%etas(l)
    end do
   end do
   end
!
!********************************************************************
!
!********************************************************************
!
   subroutine Vxi_pack(qm2_params,qmmm_mpi,qm2_struct,qm2ds,qmmm_struct,xi,eta)
   use qmmm_module,only: qm2_structure, qmmm_mpi_structure
   use qm2_davidson_module
   use qmmm_struct_module, only : qmmm_struct_type
   use qm2_params_module,  only : qm2_params_type

   implicit none
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qm2_structure),intent(inout) :: qm2_struct
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
   _REAL_ xi(qm2ds%Lt),eta(qm2ds%Lt)

! --- if this is the first time in this routine, load coulomb matrix
   write(6,*) 'qm2ds%iderivfl', qm2ds%iderivfl  
   if (qm2_params%vxi_first) then
      write(6,*) 'Vxi: first time in this routine'

      if (index(qm2_params%keywr,'INDO').NE.0.AND.index(qm2_params%keywr,'MINDO').EQ.0) then
         write(6,*) qm2_params%keywr,' hamiltonian requested'
         write(6,*)  'Use *Z program for ', qm2_params%keywr
         stop
      endif 
         qm2_params%vxi_first=.false.
         write(6,*) 'm2_params%vxi_first', qm2_params%vxi_first
   endif
   eta(:)=0.0
      
   !write(6,*) 'calling fock2'
  ! call print_1d(qm2ds%W, 'qm2ds%W before fock2')

  ! call print_1d(eta, 'eta before fock2')
   call qm2_fock2(qmmm_mpi,qm2_params,qm2_struct, qmmm_struct, & 
                                                      eta, &     ! F
                                                      xi, &      ! Ptot
                                                      qm2ds%W, & ! W 2c2e
                                                      qm2_params%orb_loc)

  ! call print_1d(eta, 'eta AFTER fock2')
                                        


   if (qm2ds%iderivfl.eq.0) then ! We are not in analytic derivatives   
    !  write(6,*) 'calling fock1'
     ! call print_1d(eta, 'eta BEFORE fock1')
      call qm2_fock1(qmmm_mpi,qm2_params,qmmm_struct, eta,xi) 
    !  call print_1d(eta, 'eta AFTER fock1')

   endif

  ! write(6,*) 'before eta(:)=eta(:)*2.0'
   eta(:)=eta(:)*2.0 !Why *2.0? Is it for the commutator in L(xi)?

   return
   end
!
!********************************************************************
!
!********************************************************************
!
   subroutine Vxi_packA(qm2_params,qmmm_mpi,qm2_struct,qm2ds,qmmm_struct,xi,eta)
   use qmmm_module,only: qm2_structure, qmmm_mpi_structure
   use qm2_davidson_module
   use qmmm_struct_module, only : qmmm_struct_type
   use qm2_params_module,  only : qm2_params_type

   implicit none
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qm2_structure),intent(inout) :: qm2_struct
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   _REAL_ xi(qm2ds%Lt),eta(qm2ds%Lt)

! --- if this is the first time in this routine, load coulomb matrix
   if (qm2_params%vxia_first) then
      if (index(qm2_params%keywr,'INDO').NE.0.AND.index(qm2_params%keywr,'MINDO').EQ.0) then
       write(6,*)  qm2_params%keywr,' hamiltonian requested'
       write(6,*)  'Use *Z program for ', qm2_params%keywr
       stop
      endif 
         qm2_params%vxia_first=.false.
   endif
      eta(:)=0.0
   !   call print_1d(xi, 'xi before fock2 in Vxi_packA')
   call qm2_fock2(qmmm_mpi,qm2_params,qm2_struct,qmmm_struct, eta,xi,qm2ds%W,qm2_params%orb_loc)
  ! call print_1d(eta, 'eta after fock2 in Vxi_packA')
   if (qm2ds%iderivfl.eq.0) then ! We are not in analytic derivatives
      call qm2_fock1_skew(qm2_params,qmmm_mpi,qmmm_struct, eta,xi)
   endif
      eta(:)=eta(:)*2.0
   return
   end
! MODULE ArrayModule
!    contains


! END MODULE ArrayModule
      
