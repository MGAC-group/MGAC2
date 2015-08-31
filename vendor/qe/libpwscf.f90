! This code is adapted from the QE groups
! COUPLE and pw.x codes. The copyright of the
! QE group is included to show 


!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE pwstart(lib_comm,nim,npt,npl,nta,nbn,ndg,retval,infile) BIND(C)
  !----------------------------------------------------------------------------
  !
  ! ... Library interface to the Plane Wave Self-Consistent Field code
  !
  USE ISO_C_BINDING
  USE environment,       ONLY : environment_start
  USE mp_global,         ONLY : mp_startup
  USE read_input,        ONLY : read_input_file
  USE command_line_options, ONLY: set_command_line
  USE parallel_include
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE cell_base,        ONLY : fix_volume, fix_area
  USE control_flags,    ONLY : conv_elec, gamma_only, lscf, twfcollect
  USE control_flags,    ONLY : conv_ions, istep, nstep, restart, lmd, lbfgs
  USE force_mod,        ONLY : lforce, lstres, sigma, force
  USE check_stop,       ONLY : check_stop_init, check_stop_now
  USE mp_images,        ONLY : intra_image_comm
  USE qmmm,             ONLY : qmmm_initialization, qmmm_shutdown, &
                               qmmm_update_positions, qmmm_update_forces
  
  !
  IMPLICIT NONE
  !INTEGER, INTENT(IN)    :: lib_comm, nim, npt, npl, nta, nbn, ndg
  !INTEGER, INTENT(INOUT) :: retval
  !CHARACTER(LEN=80)      :: infile_

  INTEGER (kind=C_INT), VALUE :: lib_comm, nim, npt, npl, nta, nbn, ndg
  INTEGER (kind=C_INT), INTENT(OUT) :: retval
  CHARACTER (kind=C_CHAR), INTENT(IN) :: infile(*)
  INTEGER  :: i, lib_comm_, nim_, npt_, npl_, nta_, nbn_, ndg_, retval_
  CHARACTER(LEN=1024)  :: infile_



  INTEGER :: exit_status

  lib_comm_ = lib_comm
  nim_ = nim
  npt_ = npt
  npl_ = npl
  nta_ = nta
  nbn_ = nbn
  ndg_ = ndg
  retval = 0
  infile_ = ' '
  !
  ! ... Copying a string from C to Fortran is a bit ugly.
  DO i=1,1024
      IF (infile(i) == C_NULL_CHAR) EXIT
      infile_ = TRIM(infile_) // infile(i)
  END DO


  !
  !
  CALL set_command_line( nimage=nim_, npot=npt_, npool=npl_, ntg=nta_, &
      nband=nbn_, ndiag=ndg_ )
  CALL mp_startup ( my_world_comm=lib_comm_ )
  CALL environment_start ( 'PWSCF' )
  !
  CALL read_input_file ('PW', infile_ )

  !stripped from runpwscf, up to DO loop

  exit_status = 0
  IF ( ionode ) WRITE( unit = stdout, FMT = 9010 ) ntypx, npk, lmaxx
  !
  IF (ionode) CALL plugin_arguments()
  CALL plugin_arguments_bcast( ionode_id, intra_image_comm )
  !
  ! ... needs to come before iosys() so some input flags can be
  !     overridden without needing to write PWscf specific code.
  ! 
  CALL qmmm_initialization()
  !
  ! ... convert to internal variables
  !
  CALL iosys()
  !
  IF ( gamma_only ) WRITE( UNIT = stdout, &
     & FMT = '(/,5X,"gamma-point specific algorithms are used")' )
  !
  ! call to void routine for user defined / plugin patches initializations
  !
  CALL plugin_initialization()
  !
  CALL check_stop_init()
  !
  CALL setup ()
  !
  CALL qmmm_update_positions()
  !
  CALL init_run()
  !
  ! ... dry run: code will stop here if called with exit file present
  ! ... useful for a quick and automated way to check input data
  !
  IF ( check_stop_now() ) THEN
     CALL punch( 'config' )
     exit_status = 255
  !   RETURN
  ENDIF
  retval_=exit_status
  retval=retval_
  !
  ! ... Perform actual calculation
  !
  !CALL run_pwscf  ( retval )
  !
  !CALL stop_run( retval )
  !
  RETURN
  !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !
END SUBROUTINE pwstart
!
!
!----------------------------------------------------------------------------
SUBROUTINE pwstep(retval) BIND(C)
  !----------------------------------------------------------------------------
  !
  ! ... Library interface to the Plane Wave Self-Consistent Field code
  !
  USE ISO_C_BINDING
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE cell_base,        ONLY : fix_volume, fix_area
  USE control_flags,    ONLY : conv_elec, gamma_only, lscf, twfcollect
  USE control_flags,    ONLY : conv_ions, istep, nstep, restart, lmd, lbfgs
  USE force_mod,        ONLY : lforce, lstres, sigma, force
  USE check_stop,       ONLY : check_stop_init, check_stop_now
  USE mp_images,        ONLY : intra_image_comm
  USE qmmm,             ONLY : qmmm_initialization, qmmm_shutdown, &
                               qmmm_update_positions, qmmm_update_forces
  
  !
  IMPLICIT NONE
  !INTEGER, INTENT(IN)    :: lib_comm, nim, npt, npl, nta, nbn, ndg
  !INTEGER, INTENT(INOUT) :: retval
  !CHARACTER(LEN=80)      :: infile_

  INTEGER (kind=C_INT), INTENT(OUT) :: retval
  INTEGER :: exit_status
 
  exit_status=0
 
!   main_loop: DO
     !
     ! ... electronic self-consistency or band structure calculation
     !
     IF ( .NOT. lscf) THEN
        CALL non_scf ()
     ELSE
        CALL electrons()
     END IF
     !
     ! ... code stopped by user or not converged
     !
     IF ( check_stop_now() .OR. .NOT. conv_elec ) THEN
        IF ( check_stop_now() ) exit_status = 255
        IF ( .NOT. conv_elec )  exit_status =  2
        ! workaround for the case of a single k-point
        twfcollect = .FALSE.
        CALL punch( 'config' )
        retval=exit_status
        RETURN
     ENDIF
     !
     ! ... ionic section starts here
     !
     CALL start_clock( 'ions' )
     conv_ions = .TRUE.
     !
     ! ... recover from a previous run, if appropriate
     !
     !IF ( restart .AND. lscf ) CALL restart_in_ions()
     !
     ! ... file in CASINO format written here if required
     !
     IF ( lmd ) CALL pw2casino()
     !
     ! ... force calculation
     !
     IF ( lforce ) CALL forces()
     !
     ! ... stress calculation
     !
     IF ( lstres ) CALL stress ( sigma )
     !
     ! ... send out forces to MM code in QM/MM run
     !
     CALL qmmm_update_forces(force)
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        if (fix_volume) CALL impose_deviatoric_stress(sigma)
        !
        if (fix_area)  CALL  impose_deviatoric_stress_2d(sigma)
        !
        ! ... ionic step (for molecular dynamics or optimization)
        !
        CALL move_ions()
        !
        ! ... then we save restart information for the new configuration
        !
        IF ( istep < nstep .AND. .NOT. conv_ions ) &
           CALL punch( 'config' )
        !
     END IF
     !
     CALL stop_clock( 'ions' )
     !
     ! ... exit condition (ionic convergence) is checked here
     !
     IF ( conv_ions ) THEN
        exit_status=128
        retval=exit_status
        RETURN
     END IF
     !
     ! ... receive new positions from MM code in QM/MM run
     !
     CALL qmmm_update_positions()
     !
     ! ... terms of the hamiltonian depending upon nuclear positions
     ! ... are reinitialized here
     !
     IF ( lmd .OR. lbfgs ) CALL hinit1()
     !
!  END DO main_loop
   retval=exit_status


END SUBROUTINE pwstart
!
!
!
!
!----------------------------------------------------------------------------
SUBROUTINE pwend(retval) BIND(C)
  !----------------------------------------------------------------------------
  !
  ! ... Library interface to the Plane Wave Self-Consistent Field code
  !
  USE ISO_C_BINDING
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE cell_base,        ONLY : fix_volume, fix_area
  USE control_flags,    ONLY : conv_elec, gamma_only, lscf, twfcollect
  USE control_flags,    ONLY : conv_ions, istep, nstep, restart, lmd, lbfgs
  USE force_mod,        ONLY : lforce, lstres, sigma, force
  USE check_stop,       ONLY : check_stop_init, check_stop_now
  USE mp_images,        ONLY : intra_image_comm
  USE qmmm,             ONLY : qmmm_initialization, qmmm_shutdown, &
                               qmmm_update_positions, qmmm_update_forces
  
  !
  IMPLICIT NONE
  !INTEGER, INTENT(IN)    :: lib_comm, nim, npt, npl, nta, nbn, ndg
  !INTEGER, INTENT(INOUT) :: retval
  !CHARACTER(LEN=80)      :: infile_

  INTEGER (kind=C_INT), INTENT(OUT) :: retval
  INTEGER :: exit_status
 
  exit_status=0

  IF ( .not. lmd) CALL pw2casino()
  CALL punch('all')
  !
  CALL qmmm_shutdown()
  !
  IF ( .NOT. conv_ions )  exit_status =  3
  retval=exit_status

END SUBROUTINE pwend











