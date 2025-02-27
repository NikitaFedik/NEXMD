! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "def_time.h"
subroutine qm2_energy(qm2_params,qmmm_nml,qmewald,qm2_rij_eqns, qm_gb, qmmm_mpi, qmmm_opnq, qmmm_scratch, &
	xlbomd_struct, cosmo_c_struct,qm2_struct, qm2ds,qmmm_struct, escf,scf_mchg,&
	natom,born_radii, one_born_radii, coords, scaled_mm_charges)

   ! qm2_energy calculates the energy of the QMMM system in KCal/mol and places
   ! the answer in escf.

   !     Variables for qm-mm:
   !
   !     qmmm_struct%nquant          - Number of REAL qm atoms.
   !     coords: natom array of cartesian coordinates (Amber x array - only used occasionally,
   !                                                   most routines want qm_xcrd which is the coordinates
   !                                                   of the pair list.)
   !
   use qmmm_module, only : qm2_structure, qm_ewald_structure, &
         qm2_rij_eqns_structure, qm_gb_structure, qmmm_mpi_structure, &
         qmmm_opnq_structure, qmmm_scratch_structure
   use qm2_pm6_hof_module
   use dh_correction_module, only : dh_correction
   use constants, only : EV_TO_KCAL, zero
   use cosmo_C, only : cosmo_C_structure !cosmo_c_struct%solvent_model, cosmo_c_struct%potential_type
   use xlbomd_module, only : xlbomd_structure
   use qmmm_struct_module, only : qmmm_struct_type
   use qm2_davidson_module, only : qm2_davidson_structure_type
   use qm2_params_module,  only : qm2_params_type
   use qmmm_nml_module   , only : qmmm_nml_type
   implicit none

#ifdef MPI
   include 'mpif.h'
   integer :: ier
   _REAL_ :: tmp_send, tmp_recv
#endif

   !Passed in
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qm2_rij_eqns_structure),intent(inout) :: qm2_rij_eqns
   type(qm_gb_structure),intent(inout) :: qm_gb
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
   type(qmmm_opnq_structure),intent(inout) :: qmmm_opnq
   type(qmmm_scratch_structure),intent(inout) :: qmmm_scratch
   type(qm_ewald_structure),intent(inout) :: qmewald

   type(xlbomd_structure),intent(inout) :: xlbomd_struct
   type(qm2_structure),intent(inout) :: qm2_struct
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   type(cosmo_C_structure),intent(inout) :: cosmo_c_struct
   integer, intent(in) :: natom
   _REAL_, intent(out) :: escf
   _REAL_, intent(inout) :: scf_mchg(qmmm_struct%nquant_nlink)
   _REAL_, intent(in) :: born_radii(natom), one_born_radii(natom)
   _REAL_, intent(in) :: coords(3*natom)
   _REAL_, intent(in) :: scaled_mm_charges(natom);
   !Local
   integer i
   _REAL_ ANGLE, HTYPE
   _REAL_ ,allocatable, dimension(:,:) :: density_matrix_unpacked;

   !Initialisations required on every call...
   qm2_struct%hmatrix = zero  !Note this zeros the entire array
   qmmm_struct%enuclr_qmqm = zero
   qmmm_struct%enuclr_qmmm = zero
   !End initialisations required on every call

#ifndef SQM
   !==============================================
   !  Start Ewald Potential Calculation (non-scf)
   !==============================================
   !If qm_ewald is on calculate the potential at each QM atom due to the MM atoms
   if (qmmm_nml%qm_ewald>0) then
      call timer_start(TIME_QMMMEWALDENERGY)
      !Parallel
      call qm_ewald_mm_pot(qmmm_struct%qm_xcrd, qmmm_struct%qm_mm_pairs, qmmm_struct%qm_coords, natom, &
            qmmm_nml%qmmmrij_incore, qm2_rij_eqns%qmmmrijdata, &
            scaled_mm_charges, qmewald%kvec)
#ifdef MPI
      if (qmmm_nml%qmtheory%DFTB) then
        !At present only the master thread does DFTB calculations so if we are doing DFTB
        !Then we need to reduce the mmpot array from it's current distributed form. Hopefully
        !this can go away once the DFTB is properly parallelized.
# ifdef USE_MPI_IN_PLACE
        if (qmmm_mpi%commqmmm_master) then
          call mpi_reduce(MPI_IN_PLACE,qmewald%mmpot,qmmm_struct%nquant_nlink, &
                        MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
        else
          call mpi_reduce(qmewald%mmpot,0,qmmm_struct%nquant_nlink, &
                        MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
        end if
# else
           call mpi_reduce(qmewald%mmpot,qmmm_scratch%matsize_red_scratch,qmmm_struct%nquant_nlink, &
                             MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
           if (qmmm_mpi%commqmmm_master) &
             qmewald%mmpot(1:qmmm_struct%nquant_nlink) = qmmm_scratch%matsize_red_scratch(1:qmmm_struct%nquant_nlink)
# endif
      end if
#endif

      !if we are keeping the QM image charges fixed (qm_ewald==2) we need to calculate the
      !potential at each QM atom due to the images here. If we are not keeping the charges
      !fixed (qm_ewald==1) we do this inside the scf routine at each scf step.
      !note on the very first step of MD scf_mchg will be zero here and we will also do
      !the calculation during the SCF and so use that result. For all other steps
      !scf_mchg at this point contains the charges from converged density matrix
      !calculated from the previous MD structure.
      if (qmmm_nml%qm_ewald==2) then
         ! The following calculation of Mulliken charges is needed only if NOT using
         ! DFTB. If we are using DFTB, the charges will already be calculated anyway, and
         ! will already be stored into scf_mchg.
         if (.not. qmmm_nml%qmtheory%DFTB) then
            !if we are doing QM_Ewald but are keeping the images fixed we don't calculate the Mulliken charges during
            !the SCF but keep them fixed at the previous step's values. Hence we need to calculate them here using
            !what should be the density matrix from the previous step.
            !All threads do this to save the time involved in reducing it.
 	  
   	    allocate(density_matrix_unpacked(qm2_struct%norbs,qm2_struct%norbs));
	    
            call unpacking(qm2_struct%norbs,qm2_struct%den_matrix,density_matrix_unpacked,'s');

            do i=1,qmmm_struct%nquant_nlink
               call qm2_calc_mulliken(qm2_params,qm2_struct,i,scf_mchg(i),density_matrix_unpacked);
            end do

            deallocate(density_matrix_unpacked);

         endif
         !Parallel
         call qm_ewald_qm_pot(qmmm_struct%nquant, qmmm_struct%nlink, scf_mchg, &
               qmmm_struct%qm_coords,qmewald%kvec)
#ifdef MPI
         if (qmmm_nml%qmtheory%DFTB) then
           !At present only the master thread does DFTB calculations so if we are doing DFTB
           !Then we need to reduce the qmpot array from it's current distributed form. Hopefully
           !this can go away once the DFTB is properly parallelized.
# ifdef USE_MPI_IN_PLACE
           if (qmmm_mpi%commqmmm_master) then
             call mpi_reduce(MPI_IN_PLACE,qmewald%qmpot,qmmm_struct%nquant_nlink, &
                             MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
           else
             call mpi_reduce(qmewald%qmpot,0,qmmm_struct%nquant_nlink, &
                             MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
           end if
# else
           call mpi_reduce(qmewald%qmpot,qmmm_scratch%matsize_red_scratch,qmmm_struct%nquant_nlink, &
                             MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
           if (qmmm_mpi%commqmmm_master) &
             qmewald%qmpot(1:qmmm_struct%nquant_nlink) = qmmm_scratch%matsize_red_scratch(1:qmmm_struct%nquant_nlink)
# endif
         end if
#endif
      end if
      call timer_stop(TIME_QMMMEWALDENERGY)
   end if
   !==============================================
   !  END Ewald Potential Calculation (non-scf)
   !==============================================
#endif  /* not SQM */

   !================================================
   !  Start Generalised Born Calculation (non-scf)
   !================================================
   if (qmmm_nml%qmgb == 2) then
      call timer_start(TIME_QMMMGBENERGY)
      !Step 1 - calculate the GB potential at each QM atom due to all the MM atoms
      !Needs 2 real scratch arrays of at least natom long.
      !and 1 int scratch array of at least natom long.
      !also needs FULL natom long coordinate array (coords) NOT qm_xcrd.
      !Parallel
      qm_gb%gb_mmpot(1:qmmm_struct%nquant_nlink)=zero

#ifndef SQM
      call qmgb_calc_mm_pot(qm_gb, qmmm_mpi, qmmm_nml, & 
           & qmmm_struct, natom, qm_gb%gb_mmpot, qmmm_struct%atom_mask, &
           & scaled_mm_charges,qmmm_scratch%qm_real_scratch(1), &
           & qmmm_scratch%qm_real_scratch(natom+1), qmmm_scratch%qm_int_scratch(1), &
           & qmmm_struct%qm_coords, &
           & coords, born_radii, one_born_radii, qmmm_struct%iqmatoms)
#endif

      !Step 2 - calculate 1/fij values for QM-QM pairs. Needed in SCF but only dependent
      !         on Rij.
      !Parallel
      call qmgb_calc_qmqm_onefij(qm_gb, qmmm_mpi,qmmm_struct%nquant_nlink, qm_gb%qmqm_onefij, qmmm_struct%iqmatoms, born_radii, &
            one_born_radii, &
            qmmm_struct%qm_coords)

#ifdef MPI
      if (qmmm_nml%qmtheory%DFTB) then
        !At present only the master thread does DFTB calculations so if we are doing DFTB
        !Then we need to reduce the gb_mmpot array from the current distributed form.
        !Hopefully this can go away once the DFTB is properly parallelized. Note, the qmqm_onefij
        !array can stay distributed since dftb calculated gb_qmpot by calling the regular
        !qm_gb.f routine for this which knows about running in parallel.
# ifdef USE_MPI_IN_PLACE
        if (qmmm_mpi%commqmmm_master) then
          call mpi_reduce(MPI_IN_PLACE,qm_gb%gb_mmpot,qmmm_struct%nquant_nlink, &
                          MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
        else
          call mpi_reduce(qm_gb%gb_mmpot,0,qmmm_struct%nquant_nlink, &
                          MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
        end if
# else
        call mpi_reduce(qm_gb%gb_mmpot,qmmm_scratch%matsize_red_scratch,qmmm_struct%nquant_nlink, &
                          MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
        if (qmmm_mpi%commqmmm_master) &
          qm_gb%gb_mmpot(1:qmmm_struct%nquant_nlink) = qmmm_scratch%matsize_red_scratch(1:qmmm_struct%nquant_nlink)
# endif

      end if
#endif

      call timer_stop(TIME_QMMMGBENERGY)
   end if
   !================================================
   !   End Generalised Born Calculation (non-scf)
   !================================================
   if (qmmm_nml%qmtheory%DFTB) then  !Else we do regular semiempirical QMMM below.
      call qm2_dftb_energy(qm2_params, qmmm_nml,qmmm_scratch,qmmm_mpi,qm_gb,qmewald,qm2_struct,qmmm_struct, escf,scf_mchg)
   else
      !ALL THREADS CURRENTLY NEED TO DO HCORE_QMQM AND HCORE_QMMM to get erepul

      !       Get the QM-QM interactions:
      !       qm2_struct%hmatrix = 1-electron matrix
      !       qm2_struct%qm_qm_2e_repul = Electron repulsion integrals (2-electrons terms)
      !       ENUCL = core-core repulsions
      call timer_start(TIME_QMMMENERGYHCORE)
      
      call timer_start(TIME_QMMMENERGYHCOREQM)
      !Parallel
      call qm2_hcore_qmqm(qm2_params,qmmm_nml, qmmm_mpi, qm2_rij_eqns, &
            cosmo_c_struct,qm2_struct,qmmm_struct, &
	    qmmm_struct%qm_coords,qm2_struct%hmatrix,qm2_struct%qm_qm_2e_repul, &
            qmmm_struct%enuclr_qmqm)

      call timer_stop_start(TIME_QMMMENERGYHCOREQM,TIME_QMMMENERGYHCOREQMMM)
      !       Add on the QM-MM interactions to qm2_struct%hmatrix and qmmm_enuclr:
      !
      !        - Updated hmatrix includes interactions between QM electrons and
      !          partial charges on MM atoms (1-electron terms)
      !
      !        - Updated ENUCLR includes interactions between core QM charges
      !          and partial atomic MM charges (QM-MM core-core interactions)
      !          qmmm_int = 0 Skip QM-MM interaction
      !                   = 1 (def) Calculate FULL Multipole interaction.
      !                   = 2 Calculate Full Multipole and include AM1/PM3 Extra
      !                       core-core Gaussian terms.
      !                   = 5 Mechanical embedding: no QM-MM interaction
      if (qmmm_nml%qmmm_int>0 .and. (qmmm_nml%qmmm_int /= 5) ) then
         !Parallel
         call qm2_hcore_qmmm(qm2_params,qmmm_nml,qmmm_mpi,qm2_rij_eqns, qm2_struct,qmmm_struct, qm2_struct%hmatrix, &
               qmmm_struct%enuclr_qmmm,qmmm_struct%qm_xcrd)

         if (qmmm_nml%qmmm_switch) then
            call qm2_hcore_add_switched(qm2_params,qmmm_mpi, qm2_struct%hmatrix, qmmm_struct%switched_mmpot)
         end if
      end if
      call timer_stop(TIME_QMMMENERGYHCOREQMMM)

      call timer_stop_start(TIME_QMMMENERGYHCORE,TIME_QMMMENERGYSCF)


      !========================================================
      !   Setup Fock Matrix Prediction Prior to Calling SCF
      !========================================================
      !By default if we do not have a Fock matrix prediction routine
      !defined (qmmm_nml%fock_predict==0) then we don't do anything
      !and the density matrix from the previous SCF step will be used
      !to build the first Fock matrix. Otherwise we build a Fock matrix
      !here and we skip the first building of the Fock matrix.
      if (qmmm_nml%fock_predict == 1) then
        call timer_start(TIME_QMMMENERGYSCFFOCKPRED)
        call qm2_fock_predict(qmmm_nml,qmmm_struct%num_qmmm_calls, qm2_struct%hmatrix, qm2_struct%matsize, &
                              qm2_struct%fock_matrix,qm2_struct%fock_mat_final1, &
                              qm2_struct%fock_mat_final2, qm2_struct%fock_mat_final3, &
                              qm2_struct%fock_mat_final4)
        call timer_stop(TIME_QMMMENERGYSCFFOCKPRED)
      end if
      !=============================================================
      !   End Setup Fock Matrix Prediction Prior to Calling SCF
      !=============================================================

      !==========================
      !   Calculate SCF Energy
      !==========================
      !Parallel
      call qm2_scf(qmmm_opnq, qm2_params, qmewald, qmmm_nml, qm_gb, qmmm_mpi, qmmm_scratch, &  
            xlbomd_struct,cosmo_c_struct,qm2_struct, qm2ds, qmmm_struct, &
            qm2_struct%fock_matrix, qm2_struct%hmatrix,qm2_struct%qm_qm_2e_repul,escf, &
            qm2_struct%den_matrix,scf_mchg,qmmm_struct%num_qmmm_calls)
      !==========================
      ! End calculate SCF Energy
      !==========================
#ifndef SQM
      if (qmmm_nml%qm_ewald>0) then

         !We need to calculate the periodic interaction of QM core charges with MM Ewald potential
         call timer_start(TIME_QMMMENERGYHCOREQMMM)
         !In parallel all threads do 1->nquant in ewald core but have partial mmpot and qmpot arrays
         !so return partial ewald_core energies. Gets reduced later in main sander.
         call qm_ewald_core(qmewald%ewald_core,qm2_params%core_chg, qmewald%mmpot, qmewald%qmpot, scf_mchg, qmewald%coulpot)
         !Need to adjust the SCF and enuclr_qmmm energies by the value qm_ewald_core returns.
         qmmm_struct%enuclr_qmmm = qmmm_struct%enuclr_qmmm+qmewald%ewald_core
         call timer_stop(TIME_QMMMENERGYHCOREQMMM)
      end if
#endif
      !Add the nuclear-nuclear energy into the scf energy
      escf = escf + (qmmm_struct%enuclr_qmqm+qmmm_struct%enuclr_qmmm)*EV_TO_KCAL
      if (qmmm_opnq%useOPNQ) then
         escf=escf+(qmmm_opnq%vdWCorrection+qmmm_opnq%OPNQCorrection)*EV_TO_KCAL
      end if
      !Add on the heat of formation.
      if (qmmm_mpi%commqmmm_master) then
         escf = escf + qm2_params%tot_heat_form
         ! Add PM6 corrections to Heat of Formation
         if (qmmm_nml%qmtheory%PM6) then
            escf = escf + hofCorrection(qmmm_struct)
         end if
         ! Add on dispersion correction
         if (qmmm_nml%qmtheory%DISPERSION .or. qmmm_nml%qmtheory%DISPERSION_HYDROGENPLUS) then
            call dh_correction(qm2_params,qm2_struct, qmmm_struct%nquant_nlink, qmmm_struct%qm_coords, &
                 qmmm_struct%iqm_atomic_numbers,qmmm_nml%qmtheory, &
                 qmmm_struct%dCorrection, qmmm_struct%hCorrection)
            escf = escf + qmmm_struct%dCorrection + qmmm_struct%hCorrection
         endif
      end if
      call timer_stop(TIME_QMMMENERGYSCF)
   end if
   RETURN
end subroutine qm2_energy


