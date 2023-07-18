#include "dprec.fh"
!
!********************************************************************
!
!  Davidson diagonalization
!
!  Converted to fortran 90 (from fortran 77)
!  and interfaced with sqm(AmberTools) by
!  Kirill A. Velizhanin (kirill@lanl.gov)
!
!********************************************************************	
!
!********************************************************************
!   
   subroutine davidson(qm2_params,qmmm_nml,qmmm_mpi,cosmo_c_struct,qm2_struct,qm2ds, qmmm_struct)
   use qmmm_module,only:qm2_structure, qmmm_mpi_structure !cml-test
   use qm2_davidson_module
   use qmmm_struct_module, only : qmmm_struct_type  
   use cosmo_C, only : cosmo_C_structure            !! solvents ??? 
   use qm2_params_module,  only : qm2_params_type   !! SQM constant 
   use qmmm_nml_module   , only : qmmm_nml_type     !! QM/MM interface???

   implicit none
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi 
   type(cosmo_C_structure), intent (inout) :: cosmo_c_struct
   type(qm2_structure),intent(inout) :: qm2_struct
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds

   _REAL_ random,rranset,rranf
   _REAL_ ferr1,ferr2,fn
   _REAL_ f,f1,f2,ddot,f0,f11

   integer i,j,lprint
   integer iseed,ip,ih,j0,Mx0,imin,one,iloops
   integer kflag,nd,nd1,nd1_old,j1

   integer Lt,Mb
!
!--------------------------------------------------------------------
!
!  Initialize Davidson
!
!--------------------------------------------------------------------
   write(6,*) 'INIT'
   ! Nb - basis size, total number of atomic orbitals
   Lt=qm2ds%Nb*(qm2ds%Nb+1)/2 !! Lt - number of non-zero element for tridiagonal matrix , Nb here is matrix dim R(M,M), not used
   Mb=qm2ds%Nb**2             !!! ??? Not used? Basis squared?? 
   one=1 
   lprint=qm2ds%verbosity     !! verbosity flag
   nd=qm2ds%nd                ! Dimension of Krylov expansion in Davidson

   iloops=0

   if(qm2ds%dav_guess==0) then  ! do not use previoud Davison as guess
      qm2ds%istore=0 ! overwriting qm2ds%istore to not use guess
   end if

! Split Mx into batches of j1 size
   j1=nd/qm2ds%idavfrac   !! idavfrac = Fraction of the davidson space for 1 batch
   j1=min(j1,qm2ds%Mx)    !! Mx = number of excited states to be computed 
   nd1=min(j1+2,nd/2)     !! 
   if (qm2ds%mdflag.ge.0) j1=qm2ds%Mx !! batch size = number of excited states to be computed??
   
   ! No shift is allowed for MD points
                                      !!?? What is 
      

   if (qm2ds%irflag.gt.1) then ! if > 1
!     Calculations will start from irflag state for ceo
      j0=qm2ds%irflag-1
      !write(6,*)"begin at ", j0,qm2ds%Mx,qm2ds%irflag
      if (j0.gt.qm2ds%Mx) then ! Mx = number of excited states
         write(6,*) 'Looks like irflag= ', qm2ds%irflag
         write(6,*) 'Show that states Mx= ',qm2ds%Mx
         write(6,*) 'are already calculated, exiting'

         stop 
      end if
!----------------- Read units -----------------------------
!     Read irflag-1 states 
      open (qm2ds%modes_unit,file=trim(qm2ds%modes_b),form='unformatted',status='old')
      open (qm2ds%ee_unit,file=trim(qm2ds%ee_b),status='old')

      do j=1,j0
         read (qm2ds%modes_unit) (qm2ds%v0(i,j),i=1,qm2ds%Nrpa)
         read (qm2ds%ee_unit,*) qm2ds%e0(j)
      enddo
   
      close(qm2ds%modes_unit)
      close(qm2ds%ee_unit)

      if (j0.eq.qm2ds%Mx) goto 70

   else    ! Assume start from the beginning 
      j0=0
   end if
!-----------------------------------------------------------------

   qm2ds%Mj=j0 !! what the hell is Mj??

   if (lprint.gt.1) then
      write(6,*) 'Davidson parameters'
      write(6,*) 'mdflag=',qm2ds%mdflag
      write(6,*) 'irflag=',qm2ds%irflag
      write(6,*) 'Ncis=',qm2ds%Ncis       ! size of CIS matrix
      write(6,*) 'Mx=',qm2ds%Mx         ! number of excited states to be computed 
      write(6,*) 'Mj=',qm2ds%Mj         ! ???
      write(6,*) 'nd=',nd               !
      write(6,*) 'nd1=',nd1
      write(6,*) 'j1=',j1
      write(6,*)'qm2ds%istore=',qm2ds%istore
      write(6,*)'qm2ds%istore_M=',qm2ds%istore_M
   end if


   if(qm2ds%istore.gt.0) then ! MD point only!!!!	 
!     recover excited state vectors from AO representation in v2
!     recovered state vectors are put to v0
      do j=1,qm2ds%istore
         call site2mo(qm2ds,qm2ds%rrwork,qm2ds%v2(1,j),qm2ds%v0(1,j))
      end do
   endif

! kav: the lines below are commented out since they are not present in original davidson
! 
! CML TEST 7/15/12
!
!--------------------------------------------------------------------
!  Begin big loop
!--------------------------------------------------------------------
! find many vectors:
10 continue !! 10 statement label, probably for goto or loop termination
   if (j0.lt.qm2ds%Mx) then ! less than Mx, N excited states to be computed
   write(6,*) '====== BIG LOOP 10 ======, ITER ',iloops
   iloops=iloops+1
   j1=min(j1,(qm2ds%Mx-j0)) !! found vectors - requested vectors
   nd1=min(j1+2,nd/2)

   if (lprint.gt.3) then
      write(6,*)
      write(6,*) 'Davidson batch ', iloops
      write(6,*) 'So far found ', j0, ' states'
      write(6,*) 'out of requested ',qm2ds%Mx, ' states'
      write(6,*) 'This batch will seek',j1,' vectors'

      write(6,*) 'qm2ds%fs', qm2ds%fs
      write(6,*) 'qm2ds%e0', qm2ds%e0(j0)
      if(j0.gt.0) write(6,*) 'Shift is',qm2ds%fs+qm2ds%e0(j0),' eV'
   end if 

! order quasidiagonal:
   ! _REAL_,pointer::ehf(:) !!
   ! qm2ds%ehf=>qmmm_scratch%mat_diag_workspace(1:qm2_struct%norbs,1)
   write(6,*) 'D qm2ds%ehf', qm2ds%ehf !!
   !write(6,*) 'D qm2ds%ehf.shape', shape(ehf) !!
   i=0
   do ip=1,qm2ds%Np
    do ih=1,qm2ds%Nh
       i=i+1
      !  write(6,*) 'i in order quasidiagonal',i
      !  write(6,*) 'qm2ds%ehf(ih+qm2ds%Np)', qm2ds%ehf(ih+qm2ds%Np)
      !  write(6,*) 'qm2ds%ehf(ip)', qm2ds%ehf(ip) 
       qm2ds%rrwork(i)=qm2ds%ehf(ih+qm2ds%Np)-qm2ds%ehf(ip) ! Lancos vectors(i) = 
       write(6,*) 'qm2ds%rrwork(i)', qm2ds%rrwork(i)
    end do
   end do

   write(6,*) 'qm2ds%rrwork before sort' !!
   write(6,*)  qm2ds%rrwork !!
   write(6,*)'qm2ds%rrwork before sort - SHAPE', shape(qm2ds%rrwork) !!

   call rrdpsort(qm2ds%rrwork,qm2ds%Ncis,qm2ds%ix,1) !!  return the permutation vector resulting from
   write(6,*)'ix ', qm2ds%ix !!
   ! write(6,*) 'ix', qm2ds%ix
   ! write(6,*) 'ix shape', shape(qm2ds%ix)
   !!!!  sorting dx in increasing order and do not sort dx.
!  Account for found vectors
   do j=1,j0
    do i=1,qm2ds%Ncis
       qm2ds%rrwork(i)=qm2ds%rrwork(i)+(qm2ds%fs+qm2ds%e0(j0)) & !! ?? Wilkinson shift maybe??
          *abs(qm2ds%v0(i,j)**2-qm2ds%v0(qm2ds%Ncis+i,j)**2)     !! ??
    end do
   end do

   write(6,*) 'qm2ds%rrwork AFTER sort' !!
   write(6,*)  qm2ds%rrwork !!


! try to find vector new vectors in the batch:
   call davidson0(qm2_params,qmmm_nml,qmmm_mpi, cosmo_c_struct, qm2_struct,qm2ds,&
      qmmm_struct,qm2ds%Ncis,lprint,qm2ds%ftol0,qm2ds%ftol1,qm2ds%ferr, &
      qm2ds%Np,qm2ds%Nh,j0,j1, &
      qm2ds%e0,qm2ds%v0,kflag,qm2ds%ix, &
      qm2ds%rrwork(1),qm2ds%rrwork(2*qm2ds%Ncis+1), &
      qm2ds%rrwork(4*qm2ds%Ncis+1), &
      nd,nd1,qm2ds%vexp1,qm2ds%vexp,qm2ds%ray,qm2ds%rayv, &
      qm2ds%rayvL,qm2ds%rayvR,qm2ds%raye,qm2ds%raye1, &
      qm2ds%ray1,qm2ds%ray1a,qm2ds%ray2,qm2ds%idav,qm2ds%istore)

! Printing out found eigenvalues, error and tolerance
   if(lprint>0) then
      write(6,*)' i, e0(i), ferr(i), ftol0'
      do i=1,j0
         write(6,111) i,' +++ ',qm2ds%e0(i),qm2ds%ferr(i),qm2ds%ftol0
      end do
      write(6,*)'-------------------------------------------------'
      write(6,*) 
   end if

   if (qm2ds%mdflag.le.-3) qm2ds%Mj=j0
111   format (i3,a,g24.16,2(' ',e8.2))
      call flush(6)

! Write vectors only for BIG sizes in the case of crash/restart       	  
   
      if(qm2ds%mdflag.lt.0.and.qm2ds%Nb.gt.100) then
         open (qm2ds%modes_unit,file=trim(qm2ds%modes_b),form='unformatted')
         open (qm2ds%ee_unit,file=trim(qm2ds%ee_b))

         do j=1,j0
            write (qm2ds%modes_unit) (qm2ds%v0(i,j),i=1,qm2ds%Nrpa)
            write (qm2ds%ee_unit,*)  qm2ds%e0(j)
         end do

         close(qm2ds%modes_unit)
         close(qm2ds%ee_unit)
      end if

   goto 10
   end if
       
!
!--------------------------------------------------------------------
!
! End big loop
!
!--------------------------------------------------------------------
 write(6,*)'D I am HERE 70'
 70   continue

! end find many vectors
!--------------------------------------------------------------------
!
   if(qm2ds%mdflag.ne.0) then   ! initial MD point, store vector
!  normalize vectors
      do j=1,qm2ds%Mx
         fn=0
! CIS
         if(qm2ds%idav==1) then
            do i=1,qm2ds%Ncis
               qm2ds%v0(i+qm2ds%Ncis,j)=0.0
!               write(6,*)'D m2ds%v0(i+qm2ds%Ncis,j)', qm2ds%v0(i+qm2ds%Ncis,j) !!
               
            end do
         end if
! end CIS

         fn=ddot(qm2ds%Ncis,qm2ds%v0(1,j),one,qm2ds%v0(1,j),one) &
            -ddot(qm2ds%Ncis,qm2ds%v0(qm2ds%Ncis+1,j),one, &
            qm2ds%v0(qm2ds%Ncis+1,j),one)

         write(6,*)'D fn', fn !!
         f=1/sqrt(abs(fn))
         call dscal(qm2ds%Nrpa,f,qm2ds%v0(1,j),one)
!         write(6,*)'D qm2ds%v0(1,j) after dscal', qm2ds%v0(1,j) !!
      end do
      write(6,*)'D qm2ds%v0', qm2ds%v0 !!
!     write to the hard disk
      write(6,*)'D qm2ds%e0 before sorting', qm2ds%e0 !!
      call rrdpsort(qm2ds%e0,qm2ds%Mx,qm2ds%kx,2)
      write(6,*)'D qm2ds%e0 AFTER sorting', qm2ds%e0 !!

      if (qm2ds%mdflag.lt.0) then
         open (qm2ds%modes_unit,file=trim(qm2ds%modes_b),form='unformatted')
         open (qm2ds%ee_unit,file=trim(qm2ds%ee_b))
!
         do j=1,qm2ds%Mx
            write (qm2ds%modes_unit) (qm2ds%v0(i,qm2ds%kx(j)),i=1,qm2ds%Nrpa)
            write (qm2ds%ee_unit,*)  qm2ds%e0(j)
         end do
 
         close(qm2ds%modes_unit)
         close(qm2ds%ee_unit)
      end if
   end if

   if (qm2ds%mdflag.ge.0) then ! MD point only!!!!
!     store excited state vectors in AO representation in v2
!     qm2ds%istore=min(Mx,Mx_ss)
      if (qm2ds%istore.le.qm2ds%Mx) qm2ds%istore=qm2ds%Mx
      if (qm2ds%istore_M.le.qm2ds%Mx) qm2ds%istore_M=qm2ds%Mx
      if (qm2ds%istore.le.qm2ds%istore_M) qm2ds%istore=qm2ds%istore_M

      do j=1,qm2ds%Mx
         call mo2site(qm2ds,qm2ds%v0(1,j),qm2ds%v2(1,j),qm2ds%rrwork)
      end do
   end if

   return     
   end

!
!********************************************************************
!
!********************************************************************
!
   subroutine davidson0(qm2_params,qmmm_nml,qmmm_mpi,cosmo_c_struct,qm2_struct,qm2ds,&
      qmmm_struct,Ncis,lprint,ftol0,ftol1,ferr, &
      Np,Nh,j0,j1,&
      e0,v0,kflag,ix,&
      ee2,eta,&
      xi, &
      nd,nd1,vexp1,vexp,ray,rayv,&
      rayvL,rayvR,raye,raye1, &
      ray1,ray1a,ray2,idav,istore)

   !  ix is actually ix, args are positional !!!!!!!!!!!!!!!!!!!!!

   use qm2_davidson_module   
   use cosmo_C,only:cosmo_C_structure
   use qmmm_struct_module, only : qmmm_struct_type
   use qmmm_module,only:qm2_structure, qmmm_mpi_structure
   use qm2_params_module,  only : qm2_params_type
    use qmmm_nml_module   , only : qmmm_nml_type

   implicit none
     type(qmmm_nml_type),intent(inout) :: qmmm_nml
     type(cosmo_C_structure), intent (inout) :: cosmo_c_struct
     type(qmmm_struct_type), intent(inout) :: qmmm_struct
     type(qm2_davidson_structure_type), intent(inout) :: qm2ds
     type(qm2_structure),intent(inout) :: qm2_struct
     type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
     type(qm2_params_type),intent(inout) :: qm2_params

   !logical check_symmetry; !!JAB Testing

   integer Np,Nh,Ncis,lprint,kflag,nd,nd1,nd1_old,j0,j1
   integer one,i,j,k,n,m,icount,idav,istore,iloop,itarget
   integer info,ix(qm2ds%Ncis)
   integer test(qm2ds%Ncis)
   integer aa, bb ! for debugging loops

   !integer l,u,c
   _REAL_ ftol0,ftol1,ferr(j1+j0),ddot,fn,f2m
   _REAL_ f0,f1,f2,f3,f4,tresh2
   _REAL_ e0(Ncis),v0(qm2ds%Nrpa,qm2ds%Mx),ee2(qm2ds%Ncis)
   _REAL_ eta(qm2ds%Nrpa),xi(qm2ds%Nrpa)
!  Davidson expansion vectors:
   _REAL_ vexp1(qm2ds%Ncis,nd),vexp(qm2ds%Nrpa,nd)

! Davidson Rayleigh matrices:
   _REAL_ ray1(nd,nd),ray2(nd,nd),ray1a(nd,nd)
   _REAL_ ray(nd,nd),raye(nd),raye1(nd)
   _REAL_ rayv(nd,nd),rayvL(nd,nd),rayvR(nd,nd)
   _REAL_, allocatable :: dtmp(:)
   !_REAL_ f11;

   if (lprint.gt.4) write(6,*)' Entering davidson0'

   call clearing(nd*nd,ray1(1,1))
   call clearing(nd*nd,ray1a(1,1))
   call clearing(nd*nd,ray2(1,1))
   call clearing(nd*nd,ray(1,1))
   call clearing(nd*nd,rayv(1,1))
   call clearing(nd*nd,rayvL(1,1))
   call clearing(nd*nd,rayvR(1,1))
   call clearing(nd,raye)
   call clearing(nd,raye1)

   one=1

   icount=0
   nd1_old=0
   m=0
   n=0
   iloop=0

   if (istore.gt.0) then    ! Use vectors from the previous step
      j1=min(j1,istore)
      goto 70 
   end if 

! **** Assign Davidson trial vectors

85 continue ! used for goto 85
!---------------------------DEBUG VEXP1----------------------------------------
write(6,*)'!---------------------------DEBUG VEXP1---------------------------------------- '
write(6,*)'D vexp1 shape', shape(vexp1) !!
write(6,*)'D line ~388 checking vexp1', vexp1 ! seems to be 0s only   

do i = 1, size(vexp1, 1)
   do j = 1, size(vexp1, 2)
      write(6, ' (F8.6, 3x)', advance='no')   vexp1(j,i) ! Fortran - column major
   end do
   write(*, *)
end do

!---------------------------DEBUG VEXP1----------------------------------------


tresh2=qm2ds%tresh*0.9
do j=1,nd   
   call clearing (2*Ncis,vexp(1,j))
   call clearing (Ncis,vexp1(1,j))
end do 

! vexp and vexp1 are now zero !!JAB
! assign trial vectors based on the MO
itarget=0

write(6,*)'D nd1 line ~412', nd1
do j=1,nd1
   write(6,*)'j', j
!80    continue ! NEVER used as goto, probably just a label like end do

   write(6,*)'D HERE 80'

   itarget=itarget+1
   write(6,*)'itarget', itarget
!      if (itarget.ge.Np*Nh) goto 85 ! Restart vectors
   ! if (itarget.ge.Np*Nh) then 
   !    write(6,*)'Restart vectors: size > CIS  matrix' 
   !    goto 85 ! Restart vectors if size > CIS  matrix 
   ! end if
   if (itarget.ge.Np*Nh) goto 85 ! Restart vectors
   ! if (itarget.ge.Np*Nh) then 
   !    write(6,*)'Restart vectors: size > CIS  matrix' 
   !    goto 85 ! Restart vectors if size > CIS  matrix 
   ! end if


   f1=0.0
   ! write(6,*)'D ix', ix !!
   ! write(6,*)'D ix shape', shape(ix) !!
   ! write(6,*)'D Ncis', qm2ds%Ncis !!
   write(6,*)'DDD j0', j0 ! taken from irflag, which is 0 be default? 

   do i=1,j0
      write(6,*)'DDD i', i !!
      f1=f1+v0(ix(itarget),i)**2+v0(ix(itarget)+Ncis,i)**2
      write(6,*)'D f1', f1 !!
   end do

   if (f1.ge.tresh2) cycle  ! MO pair is not accepted !! taken from aimc_openshell
   if (lprint.gt.2) write(6,*)'++ START ',itarget,ix(itarget)
   write(6,*)'D ix', ix !!
   vexp1(ix(itarget),j) = vexp1(ix(itarget),j)+1.0



end do



! !---------------------------DEBUG VEXP1----------------------------------------
! write(6,*)'!---------------------------DEBUG VEXP1---------------------------------------- '
! write(6,*)'D vexp1 shape', shape(vexp1) !!
! write(6,*)'D line ~450 checking vexp1', vexp1 ! seems to be 0s only   

! do i = 1, size(vexp1, 1)
!    do j = 1, size(vexp1, 2)
!       write(6, ' (F8.6, 3x)', advance='no')   vexp1(i,j) ! Fortran - column major
!    end do
!    write(*, *)
! end do

! !---------------------------DEBUG VEXP1----------------------------------------



!  Orthogonolize trial vectors (vexp1) to found eigenvectors, some strange
!  function XI(V0) !!JAB
 
  !!ORTHONORMALIZE WITH XI=X+Y??
   write(6,*)'D HERE ORTHONORMALIZE WITH XI=X+Y' 
   do i=1,j0            !  does not inter loop while none found
      write(6,*)'D i', i
      !!MAKE XI
      do k=1,Ncis
 !        write(6,*)'D Ncis', Ncis
 !        write(6,*)'D test1'
         xi(k)=v0(k,i)+v0(k+Ncis,i) !xi=X+Y??
         write(6,*)'k xi', k 
      end do
      !!NORMALIZE
      f2=ddot(Ncis,xi,one,xi,one) !xi.xi
      write(6,*)'f2 ', f2
      f2=1.0D0/sqrt(abs(f2)) !normalization constant

      call dscal(Ncis,f2,xi,one) !normalize xi
     
      do j=1,nd1
         !!ORTHOGONALIZE
         f1=-ddot(Ncis,xi,one,vexp1(1,j),one) !xi.vexp1
         call daxpy(Ncis,f1,xi,one,vexp1(1,j),one) !vexp1=vexp1+xi.vexp1*xi
         !    daxpy(N, DA, DX, INCX, DY, INCY)
      end do

      f2=ddot(Ncis,vexp(1,j),one,vexp1(1,j),one) !vexp.vexp1 but vexp=0 if started from 85??
      f2=1.0D0/sqrt(abs(f2)) !normalization constant
      call dscal(Ncis,f2,vexp1(1,j),one) !normalize

      !!SAME AS ABOVE BUT WITH XI=X-Y
      do k=1,Ncis
         xi(k)=v0(k,i)-v0(k+Ncis,i) !xi=X-Y??
      end do

      f2=ddot(Ncis,xi,one,xi,one) !xi.xi
      f2=1.0D0/sqrt(abs(f2)) 
      call dscal(Ncis,f2,xi,one) !normalize xi

      do j=1,nd1   
         f1=-ddot(Ncis,xi,one,vexp1(1,j),one) 
         call daxpy(Ncis,f1,xi,one,vexp1(1,j),one)
      end do

      f2=ddot(Ncis,vexp(1,j),one,vexp1(1,j),one)
      f2=1.0D0/sqrt(abs(f2))
      call dscal(Ncis,f2,vexp1(1,j),one)
   end do

!  ORTHOGONALIZE AND NORMALIZE TRIAL VECTORS
90     continue
   write(6,*)'D HERE 90'
   do j=1,nd1

      !!NORMALIZE
      write(6,*)'vexp1(1,j)', vexp1(1,j)
      f2=ddot(Ncis,vexp1(1,j),one,vexp1(1,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(Ncis,f2,vexp1(1,j),one)

      do i=1,j-1
         
         !!ORTHOGONALIZE
         f1=ddot(Ncis,vexp1(1,i),one,vexp1(1,j),one)
         if(abs(f1).gt.qm2ds%tresh1) then !!REMOVE VERY NONORTHOGONAL VECTORS
            call dcopy(Ncis,vexp1(1,nd1),one,vexp1(1,j),one) 
            nd1=nd1-1 !!NUMBER OF VECTORS
            if(lprint.gt.2) write(6,*) 'Trial removed #',j,i,f1,nd1
            goto 90
         end if
         
         call daxpy(Ncis,-f1,vexp1(1,i),one,vexp1(1,j),one) !!ORTHOGONALIZE
        !  daxpy(N, DA, DX, INCX, DY, INCY) : DY = DY + DA * DX

      end do

      f2=ddot(Ncis,vexp1(1,j),one,vexp1(1,j),one) !!NORMALIZE AGAIN
      f2=1/sqrt(abs(f2))
      call dscal(Ncis,f2,vexp1(1,j),one)
   end do 

!  CHECK ORTHOGONALITY OF EXPANSIONS
   if(lprint.gt.0) write(6,*) 'Check expansion to expansion'
   do j=1,nd1
      do i=1,nd1
         f1= ddot(Ncis,vexp1(1,i),one,vexp1(1,j),one)
         if(f1.gt.1.0E-15.and.i.ne.j) goto 90
      end do
   end do

   if (nd1.le.j1) j1=nd1
! **** End assign Davidson trial vectors

10 continue !! go to very beggining?

! **** Write some things, check the number of iterations, etc.
   icount=icount+1 !iteration counter
   write(6,*) '====='
   write(6,*) 'icount', icount
   write(6,*) '====='
   if(lprint.gt.0) write(6,*) 'COUNT=',icount,'Exp=',nd1
   if(icount.gt.qm2ds%icount_M) then
         write(6,*) "Number of davidson iterations exceeded, exiting"
         stop
   endif
   if(lprint.gt.4) write(6,*) 'nd1,nd1_old,j0,j1',nd1,nd1_old,j0,j1

! **** Davidson Restart
   if(nd1.gt.nd) then !if the number of Krylov subspace expansions is  
      iloop=iloop+1
      if(lprint.gt.4) write(6,*)
      write(6,*) 'Davidson not converged, expansion=',nd
   end if

70     continue

   if(nd1.gt.nd.or.istore.gt.0) then  ! Use vectors from the previous step
      istore=0
      if (lprint.gt.4) write(6,*) 'Restart Davidson with guesses from the previous step'
      if(lprint.gt.4) write(6,*) 'Restart loop = ',iloop
      icount=0
      nd1_old=0
      m=0
      n=0

      if(idav.eq.2) nd1=min((nd-2),2*j1)  ! RPA
      if(idav.eq.1) nd1=min((nd-2),j1)  ! CIS
      if(lprint.gt.1) write(6,*) 'Currently have ',j1,' states'
      if(lprint.gt.1) write(6,*) 'With ',nd1,' initial guesses'

      !NEW EXPANSION VECTORS
      if(idav.eq.2) then ! RPA 
         do j=1,j1
          do i=1,Ncis ! Ncis = Ncis
             if(j.le.nd1) vexp1(i,j)=v0(i,j+j0)+v0(i+Ncis,j+j0)
             if((j+j1).le.nd1) vexp1(i,j+j1)=v0(i,j+j0)-v0(i+Ncis,j+j0)
          end do
         end do
      end if

      if(idav.eq.1) then ! CIS 
         do j=1,j1
          do i=1,Ncis
             if(j.le.nd1) vexp1(i,j)=v0(i,j+j0) 
          end do
         end do
      end if


!  Orthogonalize and normalize new expansion vectors
12    continue
      do j=1,nd1
         f2=ddot(Ncis,vexp1(1,j),one,vexp1(1,j),one)
         f2=1/sqrt(abs(f2))
         call dscal(Ncis,f2,vexp1(1,j),one)

         do i=1,j-1
            f1=ddot(Ncis,vexp1(1,i),one,vexp1(1,j),one)
            if (abs(f1).gt.qm2ds%tresh1) then
               call dcopy(Ncis,vexp1(1,nd1),one,vexp1(1,j),one) 
               nd1=nd1-1
               if(lprint.gt.2) write(6,*) 'New expansion removed #',j,i,f1,nd1
               goto 12
            end if

            call daxpy(Ncis,-f1,vexp1(1,i),one,vexp1(1,j),one)
         end do

         f2=ddot(Ncis,vexp1(1,j),one,vexp1(1,j),one)
         f2=1/sqrt(abs(f2))
         call dscal(Ncis,f2,vexp1(1,j),one)
      end do


!  Check
      do j=1,nd1
       do i=1,nd1
          f1=ddot(Ncis,vexp1(1,i),one,vexp1(1,j),one)
          if (f1.gt.1.0E-15.and.i.ne.j.and.lprint.gt.0) write(6,*) i,j,f1
          if(f1.gt.1.0E-15.and.i.ne.j) goto 12
       end do
      end do

      goto 10
   end if

! **** End Davidson Restart	
   
! Next is calculated vexp=L(vexp1), so vexp becomes the Lxi of the guess xi !!JAB
! Vexp1 was the orthonormalized guess expansion vectors from the previous
! section. Throughout comments using A and B are related to  L(xi)=L([X;Y]) similarly 
! to [A,B;B,A][X;Y] from RPA eigenvalue problem, but definately not the band 
! ABBA as one might mistakenly assume !!JAB

! !---------------------------DEBUG VEXP1 before LXi----------------------------------------
!    write(6,*)'!---------------------------DEBUG VEXP1---------------------------------------- '
!    write(6,*)'D vexp1 shape', shape(vexp1) !!
!    write(6,*)'D line ~670 checking vexp1', vexp1 ! seems to be 0s only   
   
!    do i = 1, size(vexp1, 1)
!       do j = 1, size(vexp1, 2)
!          write(6, ' (F8.6, 3x)', advance='no')   vexp1(i,j) ! Fortran - column major
!       end do
!       write(*, *)
!    end do
!    write(*, *)
!    write(6, *) 'vexp1(1,1)', vexp1(1,1)
!    !---------------------------DEBUG VEXP1----------------------------------------



   do i=nd1_old+1,nd1

      write(6, *) 'i=nd1_old+1,nd1', i ! 1-4 iterations
      call clearing (2*Ncis,eta)


   !   call print_2d(vexp, 'vexp before Lxi_testing')
      
      ! write(6, *) 'eta before dcopy', eta
      call dcopy(Ncis,vexp1(1,i),one,eta,one) ! eta is copy of vexp1(1,i)


      call print_2d(vexp1, 'vexp1 before Lxi_testing')
     ! call print_2d(vexp, 'vexp before Lxi_testing')
     ! write(6, *) 'vexp1', vexp1
     ! write(6, *) 'eta', eta
      call Lxi_testing(qm2_params,qmmm_nml,qmmm_mpi,cosmo_c_struct,qm2_struct, &
                       qm2ds,qmmm_struct,&
                       eta,vexp(1,i),cosmo_c_struct%solvent_model)
      write(6, *) '=============== calling Lxi i ============', i                      
      ! call print_2d(vexp, 'vexp after Lxi_testing')

      ! MY GUESS: vexp is returned snd filled only here

      ! write(6, *) 'vexp1(1,1)', vexp1(1,1)
      ! write(6, *) 'calling Lxi i', i           

      ! call print_2d(vexp, 'vexp BEFORE X-Y')
      
       
    !  write(6, *) 'vexp after Lxi', vexp
! CIS - set Y=0
      !call print_2d(vexp, 'vexp BEFORE CIS')
      if(idav.eq.1) then
         write(6, *) 'CIS'
         do j=1,Ncis
            vexp(Ncis+j,i)=0.0
         enddo
      endif
      !call print_2d(vexp, 'vexp after CIS')

! form vexp(1,i) having Ab and Bb to (A+B)b and (A-B)b	  
! thus vexp was L([X;Y], now it is [LX+LY;LX-LY] !!JAB

   
      do j=1,Ncis
         f1=vexp(j,i)
         vexp(j,i)=vexp(j,i)+vexp(Ncis+j,i) ! A + B
         vexp(Ncis+j,i)=f1-vexp(Ncis+j,i)
      end do
   end do
   call print_2d(vexp, 'vexp ALL AFTER CIS') 
   ! SOMEHOW, AFTER CIS LOWER part of columns is not 0 !!! MAYBE BUG?

!       stop
! **** Operations in Krylov space
! form ray1=b(A-B)b and ray2=b(A+B)b	
!That is, these are formed from vexp=L(xi) and the new guess expansions Vexp1
   write(6, *) 'FORMING RAY1 AND RAY2'
   write(6, *) 'nd1', nd1
   call print_2d(vexp, 'vexp before ray1')
   call print_2d(vexp1, 'vexp1 before ray1')
   do i=1,nd1
      do j=nd1_old+1,nd1
         ray1(i,j)=ddot(Ncis,vexp1(1,i),one,vexp(1,j),one)
         ray2(i,j)=ddot(Ncis,vexp1(1,i),one,vexp(Ncis+1,j),one)
      end do
   end do

   call print_2d(ray1, 'ray1')
   call print_2d(ray2, 'ray2')
  
   if(nd1_old.ne.0) then
      write(6, *) 'nd1_old NOT EQUAL 0'
      do i=nd1_old+1,nd1
       do j=1,nd1
          ray1(i,j)=ddot(Ncis,vexp1(1,i),one,vexp(1,j),one)
          ray2(i,j)=ddot(Ncis,vexp1(1,i),one,vexp(Ncis+1,j),one)
       end do
      end do
   end if

   nd1_old=nd1
   write(6,*)'ind1_old',nd1_old


! form ray1a=sqrt(b(A-B)b)
   f1=1.0
   f0=0.0

   call dcopy(nd*nd,ray1,one,rayvR,one)


   write(6,*)'info00',info 

   allocate(dtmp(4*nd1))
   call dsyev ('v','u',nd1,rayvR,nd,raye,dtmp,4*nd1,info) !! dsyev - eigensolver
        !Eigenvalues of ray1 in raye and eigenvectors in rayvR
   write(6,*)'info000',info
   write(6,*)'info Nrpa,nd,nd1',qm2ds%Nrpa,nd,nd1
   write(6,*)'info shapes',shape(rayvR),shape(raye),shape(xi)

   call print_1d(raye, 'raye')
   call print_2d(rayvR, 'rayvR')
 
   do j=1,nd1
      raye1(j)=Sign(Sqrt(Abs(raye(j))),raye(j)) !make raye1 sqrt(eigenvalues) with same sign as eigenvalues
   end do

   ! @@@@@@@@@@@@@@@@@@@@@@@@@@@
   write(6,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@'

   call print_1d(raye1, 'raye1')

   ! PRECONDITIONER??? 
   do j=1,nd1
    do i=1,nd1
       rayv(i,j)=rayvR(j,i)*raye1(i)
    end do
   end do
   call print_2d(rayv, 'rayv')

   call dgemm('N','N',nd1,nd1,nd1,f1,rayvR,nd,rayv,nd,f0,ray1a,nd) !
   call print_2d(ray1a, 'ray1a')

   call dgemm('N','N',nd1,nd1,nd1,f1,ray2,nd,ray1a,nd,f0,rayv,nd)
   call print_2d(rayv, 'rayv')

   call dgemm('N','N',nd1,nd1,nd1,f1,ray1a,nd,rayv,nd,f0,ray,nd)

   call print_2d(ray, 'ray')
   call dcopy(nd*nd,ray,one,rayv,one)
   call symmetr(nd,rayv) ! symmetrize the matrix
   call print_2d(rayv, 'rayv AFTER SYMMETRIZATION')

! find eigenvalues and eigenvectors of ray
   write(6,*)'info0',info
   call dsyev ('v','u',nd1,rayv,nd,raye,dtmp,4*nd1,info)
   do j=1,nd1
      raye(j) = Sign(Sqrt(Abs(raye(j))),raye(j))
   enddo
   deallocate(dtmp)
   call print_1d(raye, 'raye AFTER find eigenvalues and eigenvectors of ray')



! Solve for Right EigenVector  
! rayvR = |X+Y>=ray1a*rayv
   call dgemm('N','N',nd1,nd1,nd1,f1,ray1a,nd,rayv,nd,f0,rayvR,nd)
   call print_2d(rayvR, 'rayvR')

! Solve for Left EigenVector  
! rayvL = |X-Y> = 1/E b(A+B)b|X+Y> =1/raye * ray2 * rayvR
   call dgemm('N','N',nd1,nd1,nd1,f1,ray2,nd,rayvR,nd,f0,rayvL,nd)
   call print_2d(rayvL, 'rayvL')



   do j=1,nd1
    do i=1,nd1
       rayvL(i,j)=rayvL(i,j)/raye(j)
    end do
   end do
   call print_2d(rayvL, 'rayvL supposedly after preconditioner')

   if(lprint.gt.2) then 
      write(6,*)'info',info
      write(6,910) (raye(i),i=1,j1)
   end if 

   if(raye(1).le.0.1) goto 100           ! RPA bad behavior
! **** End operations in Krylov space

!  Form approximate eigenvectors
   do i=j0+1,j0+j1
      call clearing(2*Ncis,v0(1,i))
   end do

   do k=j0+1,j0+j1
      write(6,*)'k',k

 write(6,*)' === DAXPY LOOP ===',rayvR
 !new eigenvectors from converged onward (j0 is the number of converged vectors)
 !vexp1 is calculated at the beginning of the routine. The previous Krylov space
 !calculations result in f1 and f2 below which are used to mix vexp1 with
 !unconverged vectors in v0
 ! daxpy: Y := A * X + Y
      do i=1,nd1
        ! write(6,*)'i', i
         f1=rayvR(i,(k-j0))
         !write(6,*)'f1', f1
         f2=rayvL(i,(k-j0))
         ! write(6,*)'f2', f2
         call daxpy(Ncis,f1,vexp1(1,i),one,v0(1,k),one) !V0=V0+f1*Vexp1
         call daxpy(Ncis,f2,vexp1(1,i),one,v0(1+Ncis,k),one)
      end do

      do i=1,Ncis !Whatever is going on here, it looks like some X-Y, X+Y stuff
         f3=v0(i,k)
        ! write(6,*)'f3', f3
         v0(i,k)=-v0(i,k)-v0(Ncis+i,k)
         v0(Ncis+i,k)=v0(Ncis+i,k)-f3
      end do
      !call print_2d(v0, 'v0')
! CIS : Y=0
      if(idav.eq.1) then
         do i=1,Ncis
            v0(i+Ncis,k)=0.0
         end do
      end if
   end do

   call print_2d(v0, 'v0')

   write(6,*)' === END DAXPY LOOP ==='

!  Orthogonalize and normalize approximate eigenvectors
   do j=j0+1,j0+j1
      do i=j0+1,j-1
         f1=-ddot(Ncis,v0(1,i),one,v0(1,j),one) &
            +ddot(Ncis,v0(1+Ncis,i),one,v0(1+Ncis,j),one)
         call daxpy(2*Ncis,f1,v0(1,i),one,v0(1,j),one)
      end do

      f2=ddot(Ncis,v0(1,j),one,v0(1,j),one) &
         -ddot(Ncis,v0(1+Ncis,j),one,v0(1+Ncis,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(2*Ncis,f2,v0(1,j),one)


   end do
   write(6,*)'=== AFTER ORTHOGONALIZATION ==='
   call print_2d(v0, 'v0')
 
! FIND eigenvalue, residual vectors and residual norm:
! print results of this iteration and check for converged vectors

   if (lprint.gt.1) write(6,*)'eigenvalues and residual norm'
   n=0
   m=0

   write(6,*) '===== SECOND Lxi TESTING ======='
   ! call print_1d(eta, 'eta')
   do j=j0+1,j0+j1

!      ! write(6,*), 'v0(1,j)', v0(1,j)
! !      write(6,*), 'eta', eta
!       write(6,*) '==============================='
      ! "LXI FROM ABOVE FOR COMPARISON
      ! call Lxi_testing(qm2_params,qmmm_nml,qmmm_mpi,cosmo_c_struct,qm2_struct, &
      !                  qm2ds,qmmm_struct,&
      !                  eta,vexp(1,i),cosmo_c_struct%solvent_model)


      call Lxi_testing(qm2_params,qmmm_nml,qmmm_mpi,cosmo_c_struct, qm2_struct,qm2ds,qmmm_struct, &
                       v0(1,j),eta,cosmo_c_struct%solvent_model) !L(xi) output in eta !!JAB

      write(6,*) '===== AFTER SECOND Lxi TESTING ======='
      ! call print_1d(eta, 'eta')
      ! write(6,*) '=====after Lxi testing ======='
      ! write(6,*), 'v0(1,j)', v0(1,j)
      ! write(6,*), 'eta', eta
      ! write(6,*) '==============================='     

      f1=ddot(Ncis,v0(1,j),one,eta(1),one) &
         -ddot(Ncis,v0(1+Ncis,j),one,eta(1+Ncis),one)

      ! write(6,*), '!!! !!! !!! f1', f1       


! CIS
      if (idav.eq.1) f1=ddot(Ncis,v0(1,j),one,eta(1),one) !E=<xi|L(xi)> now in f1 !!JAB
      call dcopy(2*Ncis,eta,one,xi,one)
      call daxpy(2*Ncis,-f1,v0(1,j),one,xi,one)
      call print_1d(xi, 'xi')
      f2=ddot(Ncis,xi(1),one,xi(1),one)
      ! write (6,*) 'f2', f2
      f3=ddot(Ncis,xi(1+Ncis),one,xi(1+Ncis),one)
      ! write (6,*) 'f3', f3
! CIS
      if(idav.eq.1) then !Clear Y after Liouville equation
         f3=0.0 
         call clearing(Ncis,xi(1+Ncis))
         call clearing(Ncis,eta(1+Ncis))
      endif
      write (6,*) 'ftol0', ftol0
      write (6,*) 'e0', e0
      write (6,*) 'ferr', ferr
      write (6,*) 'ee2', ee2

      f2=f2+f3 !f(x) + f(y)
      write (6,*) 'f2 +f3', f2+f3

      if(f2.le.ftol0) then ! Converged vector
         write (6,*) '== HERE ftol0'
         n=n+1
         call dcopy(2*Ncis,v0(1,j0+n),one,v0(1,j),one) !Move converged vector to beginning
         write (6,*) 'f1', f1
         e0(j0+n)= f1 !move converged energy to beginning
         ferr(j0+n)=abs(f2)+abs(f3) !

         if(lprint.gt.3) then
            write(6,"(3i5,1x,2f14.9,2g10.3,2x,A)") &
               j,n,m,raye(j-j0),f1,f2,f3,' Converged!'
      end if
! CIS
      else if(idav.eq.1) then       

 ! FORM perturbed residual vectors for new initial guesses
         if((nd1+m).eq.nd) goto 45
         m=m+1
        
         !Not obvious how this is done, but f1/(f1-ee2(i) is some sort of energy
         !weighting and eta is L(xi) v0
         do i=1,Ncis
            vexp1(i,nd1+m)=(eta(i)-f1*v0(i,j))/(f1-ee2(i))
         end do

         f4=ddot(Ncis,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
         f4=1/sqrt(abs(f4))
         call dscal(Ncis,f4,vexp1(1,nd1+m),one)

         if(lprint.gt.3) then 
                write(6,"(3i5,1x,2f14.9,2g10.3)") &
               j,n,m,raye(j-j0),f1,f2,f3
         end if
      else
         write (6,*) '== HERE else'
         if((nd1+m).eq.nd) goto 45
         m=m+1

         do i=1,Ncis
            vexp1(i,nd1+m)=(eta(i)-eta(i+Ncis)-f1*(v0(i,j)-v0(i+Ncis,j))) &
               /(f1-ee2(i))
         end do

         f4=ddot(Ncis,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
         f4=1/sqrt(abs(f4))
         call dscal(Ncis,f4,vexp1(1,nd1+m),one)

         ! if((nd1+m).eq.nd) goto 45
         ! m=m+1

         ! do i=1,Ncis
         !    vexp1(i,nd1+m)=(eta(i)+eta(i+Ncis)-f1*(v0(i,j)+v0(i+Ncis,j))) &
         !       /(f1-ee2(i))
         ! end do

         ! f4=ddot(Ncis,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
         ! f4=1/sqrt(abs(f4))
         ! call dscal(Ncis,f4,vexp1(1,nd1+m),one)
         if(lprint.gt.1) write(6,"(3i5,1x,2f14.9,2g10.3)") &
            j,n,m,raye(j-j0),f1,f2,f3
      end if          
   end do
     
   write(6,*) "EIGENVALUES"
   call print_1d(e0, 'e0')

45 continue
        
   if((j1-n).eq.0) then
      if(lprint.gt.4) write(6,*) 'All vectors found after loop' &
         ,iloop, ', Expansion ', nd1
      goto 100
   end if

   if(m.eq.0) then ! Restart Davidson
      nd1=nd+1
      goto 10
   end if

   if(lprint.gt.4) write(6,*) 'New perturbed m=',m
  
15 continue
   
   
   do j=1+nd1,nd1+m
      write(6,*) '== HERE after 15'
      f2=ddot(Ncis,vexp1(1,j),one,vexp1(1,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(Ncis,f2,vexp1(1,j),one)

      do i=1,j-1
         f1=ddot(Ncis,vexp1(1,i),one,vexp1(1,j),one)
         if(abs(f1).gt.qm2ds%tresh) then
            call dcopy(Ncis,vexp1(1,nd1+m),one,vexp1(1,j),one) 
            m=m-1     
            if(lprint.gt.2) write(6,*) 'exp removed #',j,i,f1,nd1+m
            goto 15
         end if

         call daxpy(Ncis,-f1,vexp1(1,i),one,vexp1(1,j),one)
      end do

      f2=ddot(Ncis,vexp1(1,j),one,vexp1(1,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(Ncis,f2,vexp1(1,j),one)
   end do
!  Check
   do j=1+nd1,nd1+m
    do i=1,j-1
       f1=ddot(Ncis,vexp1(1,i),one,vexp1(1,j),one)
       if(f1.gt.1.0E-15.and.i.ne.j) goto 15          
    end do
   end do


   if(lprint.gt.1) write(6,*) 'Perturbed left m=',m

   if(m.eq.0) then      
      write(6,*) 'Run out of expansion vectors'
      goto 100 
   end if

   if(iloop.gt.qm2ds%iloop_M) then
      write(6,*)'Davidson not converged after ',iloop,' loops'
      goto 100 
   end if

   if((nd1+m).gt.(nd-1).and.nd1.ne.nd) then
      nd1=nd
   else
      nd1=nd1+m
   endif

   goto 10

100   continue
   j0=j0+n
   if (n.eq.0) then
      write(6,*)'Could not go further ',j0,' vectors'
      kflag=1
   end if
    
   if(lprint.gt.0) then
      write(6,*) '@@@@ Davidson subroutine Found vectors',j0
   end if

910   format(' ', 100g11.3)

   return
   end
!
