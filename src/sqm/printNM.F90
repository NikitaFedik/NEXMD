!!********************************************************************
!
!  Subroutine to print normal modes after Davidson diagonlization
!  that is, print the array qm2ds%v2, the normal modes in AO rep
!
!********************************************************************
!
!--------------------------------------------------------------------
!
!  Josiah A. Bjorgaard NDSU 2013
!
!--------------------------------------------------------------------

   subroutine printNM_AO(qm2ds,filenumber)
   use qm2_davidson_module
   implicit none
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds
   integer i,j,k,filenumber,n
   n=qm2ds%Nb
   write(filenumber,*)'Atomic orbital basis'
   do i=1,qm2ds%Mx
   	write(filenumber,*) 'Printing normal mode',i,'with',n,'by',n,'values'
   	do j=0,n-1
   		write(filenumber,"(1000F7.4)") (qm2ds%v2(j*n+k,i),k=1,n)
   	end do
        write(6,*)'v2,test2,',qm2ds%v2(:,i)
   end do
   flush(filenumber)
   end

   subroutine printNM_MO(qm2ds,filenumber)
   use qm2_davidson_module
   implicit none
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds
   integer i,n,m,filenumber
   write(filenumber,*)'Molecular orbital basis'
   write(filenumber,*) 'Nocc:',qm2ds%Np,'Nvirt:',qm2ds%Nh
   do i=1,qm2ds%Mx
   	write(filenumber,*) 'Printing p-h for normal mode',i
   	do n=0,qm2ds%Np-1
   		write(filenumber,"(1000F7.4)") (qm2ds%v0(n*qm2ds%Nh+m,i),m=1,qm2ds%Nh)
   	end do   
        write(filenumber,*) 'Printing h-p for normal mode',i
        do n=0,qm2ds%Np-1
                write(filenumber,"(1000F7.4)") (qm2ds%v0(n*qm2ds%Nh+m+qm2ds%Ncis,i),m=1,qm2ds%Nh)
        end do     
   end do
   flush(filenumber)
   end subroutine
    

!!*********************************************************************
!
!  Subroutine to calculate the point charges associated with the
!  Multipolar expansions of the one-center two-orbital charge
!  Densities
!  See Dewar, M.; Thiel, W., Theoret. Chim. Acta 1977 46, 89-104 for an
!  explanation of various multipolar moments used in this expansion
!
!**********************************************************************
!
! Output:
!
!       
!-----------------------------------------------------------------------
!
!   Josiah A. Bjorgaard NDSU 2013
!
!----------------------------------------------------------------------


subroutine printCfitNM(qm2_params,qm2ds,qmmm_struct,filenumber)
 use qmmm_module
! use constants
 use qm2_davidson_module
 use qmmm_struct_module, only : qmmm_struct_type
   use qm2_params_module,  only : qm2_params_type

      implicit none
      type(qm2_params_type),intent(inout) :: qm2_params
      type(qm2_davidson_structure_type), intent(inout) :: qm2ds
      type(qmmm_struct_type), intent(in) :: qmmm_struct
      real :: dipd(3),dipod(3,2),coords(3),coords0(3),charge,D1,D2,tcharge,trace,&
       BohrtoA
      integer :: filenumber,orb_beg,orb_end,norb1,norb2,nquant,i,nmode,norbs_atom
      character*2 :: orbt(4)

      ! Check if there are any atoms with d orbitals
      do i = 1,qmmm_struct%nquant_nlink
         if ( qm2_params%natomic_orbs(i) > 5 ) return
      end do
1000 format (F15.10,2X,F15.10,2X,F15.10,2X,E17.10,2X,A2,2X,A2,2X,A2,1X,I4)

BohrtoA=0.529177249! A/Bohr

 orbt(1)='s'
 orbt(2)='px'
 orbt(3)='py'
 orbt(4)='pz'

 !loop over the number of normal modes to represent as charges
   do nmode=1,qm2ds%Mx
   
   call mo2site(qm2ds,qm2ds%v0(:,nmode), qm2ds%xi_scratch, qm2ds%eta_scratch)

   write(filenumber,*) 'Mode',qm2ds%Mx!!Test
   tcharge=0.0
   trace=0.0
   dipd=0.d0
   dipod=0.d0
      !Print header
      write(filenumber,*)
      write(filenumber,*)
      write(filenumber,*) 'NDDO Style Multipole Charges for Normal Mode ',nmode
      write(filenumber,*) '    X           Y          Z            Charge    Pole Orb1 Orb2 Atom' 

     ! loop over number of atoms (qmmm_struct%nquant_nlink
      do nquant=1,qmmm_struct%nquant_nlink
            !atomic coordinates for atom nquant
            coords0=qmmm_struct%qm_coords(:,nquant)

            !charge distance coefficients D1,D2 for atom (see. desc. below)
            D1=qm2_params%multip_2c_elec_params(1,nquant)*BohrtoA !Bohr to Angstrom
            D2=qm2_params%multip_2c_elec_params(2,nquant)*BohrtoA !Bohr to Angstrom
            !beginning and final indices of an atom's orbitals in the ground state
            !eigenvectors or, sort of, in the normal modes
            orb_beg=qm2_params%orb_loc(1,nquant)
            orb_end=qm2_params%orb_loc(2,nquant)
            !write(6,*) 'orb beginning,orbending',orb_beg,orb_end,size(qm2ds%xi_scratch)
            ! Determe of the values difference between the orbitals
             norbs_atom=orb_end-orb_beg+1

          ! Loop over the number of orbitals in this atom twice to select the blocks
        ! of the transition density associated with only this atom, since NDDO
          ! approximation. This captures all interatomic basis function pairs

          do norb1 = 1 , norbs_atom
             do norb2 = 1 , norbs_atom
               charge=qm2ds%xi_scratch((orb_beg+norb1-2)*qm2ds%Nb+orb_beg+norb2-1)! coef from dens which

               !There are four cases of multipolar moments
               ! 1.) ss = monopole
               ! 2.) sp = dipole
               ! 3,) pp = monopole + linear quadropole
               ! 4.) pp = monopole + square quadropole
               ! Distance from charges to atomic center are related to D1 for
               ! case (2) and D2 for case (3,4) where definitions can be found in
               ! the ref. above, and in the MOPAC manual defined as DD and QQ.
               ! D1 and D2 parameters are contained in the variable
               ! qm2_params%multip_2c_elec_params
               !
               ! The orbitals are numbered as s,px,py,pz or just s in the list of
               ! coefficients from the normal mode or eigenvector

               ! M=monopole, D=dipole, LQ=linear quadropole, SQ=square
               ! quadropole
                
               ! Case 1:ss-M
               if (norb1 == 1 .AND. norb2 == 1) then
               coords=coords0 !xyz of M chg
               dipd=dipd-charge*coords/BohrtoA
               write(filenumber,1000) coords,charge,'M',orbt(norb1),orbt(norb2),nquant
               trace=trace+charge
                           
               ! Case 2:sp-D
               elseif ((norb1 > 1 .AND. norb2 == 1) .OR. &
                      (norb1 ==1 .AND. norb2>1)) then
               coords=coords0 !central coordinates
               if (norb2==1) then
               dipod(norb1-1,1)=dipod(norb1-1,1)-charge*D1/BohrtoA
               coords(norb1-1)=coords0(norb1-1)+D1 !xyz of first charge of dipole +
               write(filenumber,1000)coords,charge*0.5,'D',orbt(norb1),orbt(norb2),nquant
               coords(norb1-1)=coords0(norb1-1)-D1 !xyz of second dip charge -
               write(filenumber,1000)coords,charge*(-0.5),'D',orbt(norb1),orbt(norb2),nquant
               elseif (norb1==1) then
               dipod(norb2-1,2)=dipod(norb2-1,2)-charge*D1/BohrtoA
               coords(norb2-1)=coords0(norb2-1)+D1 !xyz of first charge of dipole +
               write(filenumber,1000) coords,charge*0.5,'D',orbt(norb1),orbt(norb2),nquant
               coords(norb2-1)=coords0(norb2-1)-D1 !xyz of second dip charge -
               write(filenumber,1000) coords,charge*(-0.5),'D',orbt(norb1),orbt(norb2),nquant
               endif

               ! Case 2:pp-M,LQ
               elseif (norb1 > 1 .AND. norb2 > 1 .AND. norb1 == norb2) then
               coords=coords0 !xyz of M and first LQ chg
               dipd=dipd-charge*coords/BohrtoA
               write(filenumber,1000) coords,charge , 'M',orbt(norb1),orbt(norb2),nquant
               write(filenumber,1000) coords,charge*(-0.5),'LQ',orbt(norb1),orbt(norb2),nquant
               coords(norb1-1)=coords0(norb1-1)+2*D2 !xyz of second LQ chg +
               write(filenumber,1000) coords,charge*0.25,'LQ',orbt(norb1),orbt(norb2),nquant
               coords(norb1-1)=coords0(norb1-1)-2*D2 !xyz of third LQ chg -
               write(filenumber,1000) coords,charge*0.25,'LQ',orbt(norb1),orbt(norb2),nquant
               trace=trace+charge              

               ! Case 4 :pp-SQ
               elseif (norb1 > 1 .AND. norb2 > 1 .AND. norb1 .NE. norb2) then
               coords=coords0 !central coordinates
               coords(norb1-1)=coords0(norb1-1)+D2
               coords(norb2-1)=coords0(norb2-1)+D2 !xyz of first SQ chg ++
               write(filenumber,1000) coords,charge*0.25,'SQ',orbt(norb1),orbt(norb2),nquant
               coords(norb1-1)=coords0(norb1-1)-D2 !xyz of second SQ chg +-
               write(filenumber,1000) coords,charge*(-0.25) ,'SQ',orbt(norb1),orbt(norb2),nquant
               coords(norb2-1)=coords0(norb2-1)-D2 !xyz of third SQ chg --
               write(filenumber,1000) coords,charge*0.25 ,'SQ',orbt(norb1),orbt(norb2),nquant
               coords(norb1-1)=coords0(norb1-1)+D2 !xyz of fourth SQ chg -+
               write(filenumber,1000) coords,charge*(-0.25) ,'SQ',orbt(norb1),orbt(norb2),nquant
               
               else
               write(filenumber,*) 'ERROR'
              
               end if

             enddo !over orbitals
          enddo !over orbitals

        enddo !over atoms
    enddo !over nms
flush(filenumber)
return
end subroutine

!!********************************************************************
!
!  Subroutine to print normal modes in binary (MO representation)
!  (and energies in ee.b)
!  
!********************************************************************
!
!--------------------------------------------------------------------
!
!  Josiah A. Bjorgaard LANL 2015
!
!--------------------------------------------------------------------

   subroutine printNM_binary(qm2_struct,qm2ds,fn1,fn2,fn3)
   use qm2_davidson_module
   use qmmm_module
   implicit none
   type(qm2_structure),intent(inout) :: qm2_struct
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds
   integer i,j,fn1,fn2,fn3
   open (fn1,file='modes.b',form='unformatted')
   open (fn2,file='ee.b')
   do j=1,qm2ds%Mx
      write (fn1) (qm2ds%v0(i,j),i=1,qm2ds%Nrpa)
      write (fn2,*)  qm2ds%e0(j)
   enddo
   open (fn3,file='rhogr.b',form='unformatted')
   do j=1,qm2_struct%matsize
      write (fn3) qm2_struct%den_matrix(j)
   enddo
   close(fn1)
   close(fn2)
   close(fn3)
   end subroutine

