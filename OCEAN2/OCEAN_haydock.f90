! Copyright (C) 2015 - 2025 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_haydock
  use AI_kinds
  use OCEAN_timekeeper
  use iso_c_binding
  implicit none
  private
  save  


!  REAL(DP), ALLOCATABLE :: a( : )
!  REAL(DP), ALLOCATABLE :: b( : )
  REAL(DP), ALLOCATABLE :: real_a( : )
  REAL(DP), ALLOCATABLE :: imag_a( : )
  REAL(DP), ALLOCATABLE :: real_b( : )
  REAL(DP), ALLOCATABLE :: imag_b( : )
  REAL(DP), ALLOCATABLE :: real_c( : )
  REAL(DP), ALLOCATABLE :: imag_c( : )
 

  REAL(DP) :: el, eh, gam0, eps, nval,  ebase, gamGauss
  REAL(DP) :: gres, gprc, ffff, ener
  REAL(DP) :: e_start, e_stop, e_step
  REAL(DP), ALLOCATABLE :: e_list( : )

  real(DP) :: eps1Conv( 3 )

  
  INTEGER  :: haydock_niter = 0
  INTEGER  :: ne
  INTEGER  :: nloop
  INTEGER  :: inv_loop
  

  CHARACTER(LEN=3) :: calc_type
  LOGICAL  :: echamp
  LOGICAL  :: project_absspct
  LOGICAL  :: is_first = .true.

  LOGICAL  :: val_loud = .true.
  LOGICAL  :: complex_haydock = .true.

  public :: OCEAN_haydock_setup, OCEAN_haydock_do

  contains

  subroutine OCEAN_haydock_pseudoHerm( sys, hay_vec, ierr )
    use AI_kinds, only : DP
    use OCEAN_energies
    use OCEAN_system, only : o_system 
    use OCEAN_psi 
    use OCEAN_action, only : OCEAN_xact
    use OCEAN_mpi, only : myid, root, comm
    use OCEAN_filenames, only : OCEAN_filenames_spectrum
    use OCEAN_constants, only : Hartree2eV 
    
    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    !JTV need to figure out a work-around. Right now hay_vec is inout because of
    ! a depndency tracing back to calling copy and possibly copy_min, and
    ! possibly needing to go min->full, copy full, full->min
    type( ocean_vector ), intent( inout ) :: hay_vec
      
    complex(DP) :: ctmp, psqrtc, psqrtc2, ctmp2
    real(DP) :: aitmp, atmp, ibtmp, rbtmp, ictmp, rctmp, rb0, ib0
    integer :: iter
    type( ocean_vector ) :: psi_s, psi_t, psi_r
    type( ocean_vector ) :: psi_tmp
    
    
    real(DP), allocatable :: ReKrylovOverlaps( : ), ImKrylovOverlaps( : )
    complex(DP),allocatable :: overlaps(:)
    character( LEN = 40 ) :: abs_filename

    call OCEAN_psi_new( psi_s, ierr, hay_vec )
    if( ierr .ne. 0 ) return

    call OCEAN_energies_sfact_copy( sys, hay_vec, psi_s, ierr )
    if( ierr .ne. 0 ) return
    
    call OCEAN_psi_new( psi_r, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_zero_min( psi_r, ierr )
    if( ierr .ne. 0 ) return
    
    call OCEAN_psi_new( psi_t, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_zero_min( psi_t, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_new( psi_tmp, ierr )
    if( ierr .ne. 0 ) return

    allocate( ReKrylovOverlaps( 0:haydock_niter ), ImKrylovOverlaps( 0:haydock_niter ) )
    ReKrylovOverlaps( : ) = 0.0_DP
    ImKrylovOverlaps( : ) = 0.0_DP


    call OCEAN_xact( sys, sys%interactionScale, psi_s, psi_tmp, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_prep_min2full( psi_tmp, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_start_min2full( psi_tmp, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_finish_min2full( psi_tmp, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_copy_min( psi_t, psi_tmp, ierr )
    call OCEAN_energies_allow( sys, psi_t, ierr, sfact=.true.)


    call OCEAN_psi_dot( psi_s, psi_t, rb0, ierr, ib0, involution=.true. )
    if( ierr .ne. 0 ) return
    write(6,*) rb0, ib0

    ctmp = sqrt(cmplx(rb0,ib0,DP))
    rb0 = real(ctmp,DP)
    ib0 = aimag(ctmp)
    if( myid .eq. 0 ) write(6,*) 'b0:', rb0, ib0

    
    call OCEAN_psi_divide( psi_s, ierr, rb0, ib0 )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_divide( psi_t, ierr, rb0, ib0 )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_dot( psi_s, hay_vec, ReKrylovOverlaps( 0 ), ierr, ImKrylovOverlaps( 0 ) )
    if( ierr .ne. 0 ) return

    real_b(0) = 0.0_DP
    imag_b(0) = 0.0_DP

    do iter = 1, haydock_niter
      if( sys%cur_run%have_val ) then
        if( myid .eq. root ) write(6,*)   " iter. no.", iter-1
      endif

      call OCEAN_psi_copy_min( psi_tmp, psi_t, ierr )
      call OCEAN_energies_allow( sys, psi_tmp, ierr, sfact=.true. )
      call OCEAN_psi_dot( psi_t, psi_tmp, real_a(iter-1), ierr, imag_a(iter-1) )
      if( ierr .ne. 0 ) return

      call OCEAN_psi_axpy( -real_a(iter-1), psi_s, psi_t, ierr, -imag_a(iter-1))
      if( ierr .ne. 0 ) return

      call OCEAN_psi_axpy( -real_b(iter-1), psi_r, psi_t, ierr, -imag_b(iter-1))
      if( ierr .ne. 0 ) return

      call OCEAN_psi_copy_min( psi_r, psi_s, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_copy_min( psi_s, psi_t, ierr )
      if( ierr .ne. 0 ) return

      call OCEAN_psi_prep_min2full( psi_s, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_start_min2full( psi_s, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_finish_min2full( psi_s, ierr )
      if( ierr .ne. 0 ) return 
 
      call OCEAN_xact( sys, sys%interactionScale, psi_s, psi_tmp, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_prep_min2full( psi_tmp, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_start_min2full( psi_tmp, ierr ) 
      if( ierr .ne. 0 ) return
      call OCEAN_psi_finish_min2full( psi_tmp, ierr )
      if( ierr .ne. 0 ) return

      call OCEAN_psi_copy_min( psi_t, psi_tmp, ierr )
      call OCEAN_energies_allow( sys, psi_t, ierr, sfact=.true.)
      call OCEAN_psi_dot( psi_s, psi_t, rbtmp, ierr, ibtmp, involution=.true. )

      ctmp = sqrt( cmplx( rbtmp, ibtmp, DP ) )
      real_b(iter) = real(ctmp,DP)
      imag_b(iter) = aimag(ctmp)
      real_c(iter) = real_b(iter)
      imag_c(iter) = 0.0_DP

      
      call OCEAN_psi_divide( psi_s, ierr, real_b(iter), imag_b(iter) )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_divide( psi_t, ierr, real_b(iter), imag_b(iter) )
      if( ierr .ne. 0 ) return


      call OCEAN_psi_dot( psi_s, hay_vec, ReKrylovOverlaps( iter ), ierr, ImKrylovOverlaps( iter ) )
      if( ierr .ne. 0 ) return

      if( myid .eq. 0 ) then
        write ( 6, '(1x,6(f20.13,2x),i6)' ) real_a(iter-1)*Hartree2eV, imag_a(iter-1) * Hartree2eV, &
                                                      real_b(iter) * Hartree2eV, imag_b(iter) * Hartree2eV, &
                                                      real_c(iter) * Hartree2eV, imag_c(iter) * Hartree2eV, iter
        write(6,*) ReKrylovOverlaps( iter-1 ), ImKrylovOverlaps( iter-1)
      endif
    enddo

    if( myid .eq. 0 ) then
      call write_lanczos( haydock_niter, sys, hay_vec%kpref, ierr )
      if( ierr .ne. 0 ) return

      ! later parallelize over energy points

      call OCEAN_filenames_spectrum( sys, abs_filename, ierr )
      if( ierr .ne. 0 ) return
      open( unit=99, file=abs_filename, form='formatted', status='unknown' )
      rewind 99
      allocate( overlaps( haydock_niter ) )
      overlaps( 1:haydock_niter) = rb0*cmplx( ReKrylovOverlaps( 0:haydock_niter-1), & 
                                              ImKrylovOverlaps(0:haydock_niter-1),DP )
      call write_val_tri( 99, haydock_niter, hay_vec%kpref, sys%celvol, sys%valence_ham_spin, &
                          sys%cur_run%semiTDA, sys%cur_run%backf, overlaps, ierr )
      close( 99 )
      deallocate( overlaps )
      if( ierr .ne. 0 ) return

    endif
    
    
  end subroutine OCEAN_haydock_pseudoHerm



  subroutine OCEAN_haydock_do( sys, hay_vec, restartBSE, newEps, ierr )
    use OCEAN_system, only : o_system
    use OCEAN_psi, only : ocean_vector

    integer, intent( inout ) :: ierr
    logical, intent( inout ) :: restartBSE
    real(DP), intent( inout ) :: newEps
    type( o_system ), intent( in ) :: sys
    !JTV need to figure out a work-around. Right now hay_vec is inout because of
    ! a depndency tracing back to calling copy and possibly copy_min, and
    ! possibly needing to go min->full, copy full, full->min
    type( ocean_vector ), intent( inout ) :: hay_vec


    if( sys%bwflg ) then
      write(6,*) 'pseudo'
      call OCEAN_haydock_pseudoHerm( sys, hay_vec, ierr )
    else
      call OCEAN_haydock_Herm_do( sys, hay_vec,  restartBSE, newEps, ierr )
    endif

  end subroutine OCEAN_haydock_do


  subroutine OCEAN_haydock_Herm_do( sys, hay_vec, restartBSE, newEps, ierr )
    use AI_kinds
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range
    use OCEAN_action, only : OCEAN_xact


    implicit none
    integer, intent( inout ) :: ierr
    logical, intent( inout ) :: restartBSE
    real(DP), intent( inout ) :: newEps
    type( o_system ), intent( in ) :: sys
    !JTV need to figure out a work-around. Right now hay_vec is inout because of
    ! a depndency tracing back to calling copy and possibly copy_min, and
    ! possibly needing to go min->full, copy full, full->min
    type( ocean_vector ), intent( inout ) :: hay_vec

    real(DP) :: imag_a, maxDiff, relArea 
    integer :: iter, haydock_niter_actual, iter2

!    character( LEN=21 ) :: lanc_filename

    type( ocean_vector ) :: psi, old_psi, new_psi

    logical :: prevConv
    
    prevConv = .false.

!    if( myid .eq. root ) write(6,*) 'entering haydock'

    call OCEAN_psi_new( psi, ierr, hay_vec )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi'
!    call OCEAN_energies_allow( sys, psi, ierr, .true. )
!    if( ierr .ne. 0 ) return

    call OCEAN_psi_new( new_psi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'new_psi'

    call OCEAN_psi_new( old_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_zero_min( old_psi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'old_psi'

    if( myid .eq. root ) then 
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', hay_vec%kpref
      write(6,*) sys%interactionScale, haydock_niter
    endif
    call MPI_BARRIER( comm, ierr )



!    call OCEAN_tk_start( tk_psisum )

    do iter = 1, haydock_niter
      if( sys%cur_run%have_val ) then
        if( myid .eq. root ) write(6,*)   " iter. no.", iter-1
      endif
!        call OCEAN_energies_allow( sys, psi, ierr )
!        if( ierr .ne. 0 ) return
!      endif

      call OCEAN_xact( sys, sys%interactionScale, psi, new_psi, ierr )
      if( ierr .ne. 0 ) return
!      if( myid .eq. root ) write(6,*) 'Done with ACT'


!      call OCEAN_energies_allow( sys, new_psi, ierr, .true. )
!      if( ierr .ne. 0 ) return

      ! This should be hoisted back up here
      call ocean_hay_ab_twoterm( sys, psi, new_psi, old_psi, iter, restartBSE, newEps, ierr )
      if( ierr .ne. 0 ) return
      if( restartBSE ) goto 11

      ! test to check convergence
      if( sys%earlyExit .and. iter .gt. sys%haydockConvergeSpacing ) then
        if( myid .eq. root ) then
          iter2 = iter - sys%haydockConvergeSpacing
          call check_convergence( iter, iter2, sys, hay_vec%kpref, maxDiff, relArea )
        endif
        call MPI_BCAST( relArea, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
        
        if( relArea .lt. sys%haydockConvergeThreshold ) then
          if( prevConv ) then
            if( myid .eq. root ) write(6,*) 'Convergence: ', iter, maxDiff, relArea
            haydock_niter_actual = iter
            goto 11
          else
            prevConv = .true.
          endif
        else
          if( myid .eq. root ) write(6,*) 'Not converged', iter, maxDiff, relArea
          prevConv = .false.
        endif
      endif
      haydock_niter_actual = iter
          

    enddo

11  continue

    call OCEAN_tk_stop( tk_psisum )
    call MPI_BARRIER( comm, ierr )

    if( myid .eq. 0 .and. .not. restartBSE ) then
!      write(lanc_filename, '(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'lanceig_', sys%cur_run%elname, &
!        '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
      call haydump( haydock_niter_actual, sys, hay_vec%kpref, ierr )
      call redtrid(  haydock_niter_actual, sys, hay_vec%kpref, ierr )
    endif

    call OCEAN_psi_kill( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( new_psi, ierr )
    if( ierr .ne. 0 ) return
    
    call OCEAN_psi_kill( old_psi, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BARRIER( comm, ierr )

  end subroutine OCEAN_haydock_Herm_do




  subroutine OCEAN_hay_ab_twoterm( sys, psi, hpsi, old_psi, iter, restartBSE, newEps, ierr )
#ifdef __HAVE_F03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif
    use OCEAN_system, only : O_system
    use OCEAN_psi
    use OCEAN_mpi, only : root, myid, comm, &
                          MPI_DOUBLE_PRECISION, MPI_LOGICAL, MPI_INTEGER, MPI_STATUS_IGNORE
    use OCEAN_constants, only : Hartree2eV

    implicit none
    integer, intent(inout) :: ierr                  
    integer, intent(in) :: iter                     
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi 
    logical, intent(inout) :: restartBSE
    real(DP), intent(inout) :: newEps

    real(dp) :: btmp, atmp, aitmp
    integer :: ialpha, ikpt, arequest, airequest, brequest 
   
    ! hpsi -= b(i-1) * psi^{i-1}
    btmp = -real_b(iter-1)
    ! y:= a*x + y
    call OCEAN_psi_axpy( btmp, old_psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    ! calc ctmp = < hpsi | psi > and begin Iallreduce
    call OCEAN_psi_dot( hpsi, psi, atmp, ierr, ival=aitmp, rrequest=arequest, irequest=airequest )
    if( ierr .ne. 0 ) return

    ! finish allreduce to get atmp
    ! we want iatmp too (for output/diagnostics), but that can wait
    call MPI_WAIT( arequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    real_a(iter-1) = atmp

    ! hpsi -= a(i) * psi^{i}
    atmp = -atmp
    call OCEAN_psi_axpy( atmp, psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    !
    call OCEAN_psi_dot( hpsi, hpsi, btmp, ierr, rrequest=brequest )
    if( ierr .ne. 0 ) return

    ! copies psi onto old_psi
    call OCEAN_psi_copy_min( old_psi, psi, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( brequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    real_b(iter) = sqrt( btmp )
    real_c(iter) = real_b(iter)
    btmp = 1.0_dp / real_b( iter )
    call OCEAN_psi_scal( btmp, hpsi, ierr )
    if( ierr .ne. 0 ) return

    ! copies hpsi onto psi
    call OCEAN_psi_copy_min( psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_prep_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_start_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( airequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    imag_a(iter-1) = aitmp

    if( myid .eq. 0 ) then
      write ( 6, '(1x,6(f20.13,2x),i6)' ) real_a(iter-1)*Hartree2eV, imag_a(iter-1) * Hartree2eV, &
                                                    real_b(iter) * Hartree2eV, imag_b(iter) * Hartree2eV, &
                                                    real_c(iter) * Hartree2eV, imag_c(iter) * Hartree2eV, iter

      if( mod( iter, 10 ) .eq. 0 ) then 
        call haydump( iter, sys, psi%kpref, ierr )
        if( ierr .ne. 0 ) return
        call write_lanczos( iter, sys, psi%kpref, ierr )
        if( ierr .ne. 0 ) return

        ! need to sync first and maybe need to write out old_psi above where it is (maybe?) 
        ! already distributed?
!        call OCEAN_psi_write( sys, psi, 'psi_', .false., ierr )
      endif
#ifdef __HAVE_F03
      if( ieee_is_nan( real_a(iter-1) ) ) then
#else
      if( real_a(iter-1) .ne. real_a(iter-1) ) then
#endif
        write(6,*) 'NaN detected'
        ierr = -1
        return
      endif

      if( sys%convEps .and. sys%cur_run%calc_type .eq. 'VAL' ) then
        call testConvergeEps( iter, sys, psi%kpref, sys%celvol, sys%valence_ham_spin, restartBSE, newEps )
      endif
!      call haydump( iter, sys, ierr )
    endif

    if( sys%convEps ) then
      call MPI_BCAST( restartBSE, 1, MPI_LOGICAL, root, comm, ierr )
      call MPI_BCAST( newEps, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    endif

    ! Might be moved up & out?
    call OCEAN_psi_finish_min2full( psi, ierr )
    if( ierr .ne. 0 ) return


  end subroutine  OCEAN_hay_ab_twoterm





  subroutine check_convergence( iter1, iter2, sys, kpref, maxDiff, relArea )
    use OCEAN_system, only : o_system
    implicit none
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: iter1, iter2
    real(DP), intent( in ) :: kpref
    real(DP), intent( out ) :: maxDiff, relArea

    real(DP), allocatable :: sp1(:,:), sp2(:,:)
    real(DP) :: area1, area2, diff
    integer :: i

    allocate( sp1(3,ne), sp2(3,ne) )

    select case( sys%cur_run%calc_type)
      case( 'XES', 'XAS' )
        call calc_spect_core( sp1, iter1, kpref )
        call calc_spect_core( sp2, iter2, kpref )
      case( 'VAL', 'RXS' )
        call calc_spect_val( sp1, iter1, kpref, sys%celvol, sys%valence_ham_spin, sys%cur_run%backf )
        call calc_spect_val( sp2, iter2, kpref, sys%celvol, sys%valence_ham_spin, sys%cur_run%backf )

      case default
        call calc_spect_core( sp1, iter1, kpref )
        call calc_spect_core( sp2, iter2, kpref )

    end select
    
    maxDiff = 0.0_DP
    relArea = 0.0_DP
    area1 = 0.0_DP
    area2 = 0.0_DP

    do i = 1, ne
      diff = abs(sp1(3,i) - sp2(3,i) )
      if( diff .gt. maxDiff ) maxDiff = diff
      relArea = relArea + diff
      area1 = area1 + abs(sp1(3,i))
      area2 = area2 + abs(sp2(3,i))
    enddo

    write(6,*) relArea, area1, area2
    relArea = 2.0_DP * relArea / ( area1 + area2 )


!    open(unit=99,file='check.txt',form='formatted', status='unknown')
!    do i = 1, ne
!      write(99,*) sp1(1,i), sp1(3,i), sp2(3,i)
!    enddo
!    close(99)
  
    deallocate( sp1, sp2 )

  end subroutine check_convergence



  subroutine haydump( iter, sys, kpref, ierr )
    use OCEAN_system, only : o_system
    use OCEAN_constants, only : Hartree2eV
    use OCEAN_filenames, only : OCEAN_filenames_spectrum
    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: iter
    real(DP), intent( in ) :: kpref

    integer :: ie, jdamp, jj
    real(DP), external :: gamfcn
    real(DP) :: e, gam, dr, di, ener, spct( 0 : 1 ), spkk, pi
    complex(DP) :: rm1, ctmp, disc, delta

    character( LEN=40 ) :: abs_filename
    
    call OCEAN_filenames_spectrum( sys, abs_filename, ierr )
    if( ierr .ne. 0 ) return
    
!    rm1 = -1; rm1 = sqrt( rm1 ); pi = 4.0d0 * atan( 1.0d0 )
!    open( unit=99, file='absspct', form='formatted', status='unknown' )
    open( unit=99, file=abs_filename, form='formatted', status='unknown' )
    rewind 99

    select case ( sys%cur_run%calc_type)
      case( 'XES', 'XAS' )
        call write_core( 99, iter, kpref, sys%oldXASbroaden, sys%celvol )
      case( 'VAL', 'RXS' )
        call write_val( 99, iter, kpref, sys%celvol, sys%valence_ham_spin, sys%cur_run%semiTDA, sys%cur_run%backf )

      case default
        call write_core( 99, iter, kpref, sys%oldXASbroaden, sys%celvol )
    
    end select

    close(unit=99)
    !
    return
  end subroutine haydump

  subroutine write_val( fh, iter, kpref , ucvol, val_ham_spin, semiTDA, backf )
    use OCEAN_constants, only : Hartree2eV, bohr, alphainv
    implicit none
    integer, intent( in ) :: fh, iter, val_ham_spin
    real(DP), intent( in ) :: kpref, ucvol
    logical, intent( in ) :: semiTDA, backf
    !
    integer :: ie, i
    real(DP) :: ere, reeps, imeps, lossf, fact, mu, reflct
    complex(DP) :: ctmp, arg, rp, rm, rrr, al, be, eps, refrac

    fact = kpref * real( 2 / val_ham_spin, DP ) * ucvol

    write(fh,"(a)") "#   omega (eV)      epsilon_1       epsilon_2       n"// &
      "               kappa           mu (cm^(-1))    R"//  &
      "               epsinv"

!p%kpref = 4.0d0 * pi * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 )
    do ie = 1, 2 * ne, 2
      ere = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )
!      if( backf ) ere = ere**2
#if(1)
      ctmp = cmplx( ere, gam0, DP )
      if( backf ) ctmp = ctmp**2
      if( backf ) ere = ere**2

      arg = ( ere - real_a( iter - 1 ) ) ** 2 - 4.0_dp * real_b( iter ) ** 2
      arg = sqrt( arg )

      rp = 0.5_dp * ( ere - real_a( iter - 1 ) + arg )
      rm = 0.5_dp * ( ere - real_a( iter - 1 ) - arg )
      if( aimag( rp ) .lt. 0.0_dp ) then
        rrr = rp
      else
        rrr = rm
      endif

      al =  ctmp - real_a( iter - 1 ) - rrr
      be = -ctmp - real_a( iter - 1 ) - rrr
#else
      ctmp = cmplx( ere - real_a( iter - 1 ), gam0 - imag_a( iter -1 ), DP )
      arg = sqrt( ctmp ** 2 - 4.0_dp * cmplx( real_b( iter ), imag_b( iter ), DP ) &
                                     * cmplx( real_c( iter ), -imag_c( iter ), DP ) )
      rp = 0.5_dp * ( ctmp + arg )
      if( aimag( rp ) .lt. 0.0_dp ) then
        rrr = 0.5_dp * ( ctmp + arg )
      else
        rrr = 0.5_dp * ( ctmp - arg )
      endif

      ctmp = cmplx( ere, gam0, DP )
  
      al = ctmp - cmplx( real_a( iter-1 ), imag_a( iter-1 ), DP ) - rrr
      be = -ctmp - cmplx( real_a( iter-1 ), imag_a( iter-1 ), DP ) - rrr

#endif

      do i = iter-1, 0, -1
!        al =  ctmp - real_a( i ) - real_b( i + 1 ) ** 2 / al
!        be = -ctmp - real_a( i ) - real_b( i + 1 ) ** 2 / be
        al = ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
           - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / al
        be = -ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
           - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / be
!        al =  ctmp - cmplx( real_a( i ), imag_a( i ), DP ) - real_b( i + 1 ) **2 / al
!        be = -ctmp - cmplx( real_a( i ), imag_a( i ), DP ) - real_b( i + 1 ) **2 / be
      enddo

      if( semiTDA ) then
        eps = 1.0_dp - fact / al - fact / be
      else
        eps = 1.0_dp - fact / al 
      endif

      reeps = dble( eps )
      imeps = aimag( eps )
!      rad = sqrt( reeps ** 2 + imeps ** 2 )
!      theta = acos( reeps / rad ) / 2
!      indref = sqrt( rad ) * cos( theta )
!      indabs = sqrt( rad ) * sin( theta )
!      ref = ( ( indref - 1 ) ** 2 + indabs ** 2 ) /
!   &        ( ( indref + 1 ) ** 2 + indabs ** 2 )
      lossf = imeps / ( reeps ** 2 + imeps ** 2 )

      if( backf )ere = sqrt(ere)
      refrac = sqrt(eps)
      reflct = abs((refrac-1.0d0)/(refrac+1.0d0))**2
      mu = 2.0d0 * ere * Hartree2eV * aimag(refrac) / ( bohr * alphainv * 1000 )

      write(fh,'(8(1E24.16,1X))') ere*Hartree2eV, reeps, imeps, refrac-1.0d0, mu, reflct, lossf

    enddo

  end subroutine write_val

  subroutine write_core( fh, iter, kpref, oldXASbroaden, ucVol )
    use OCEAN_constants, only : Hartree2eV, alpha, eV2Hartree, bohr, PI_DP, au2sec!, bohr2cm
    implicit none
    integer, intent( in ) :: fh, iter
    real(DP), intent( in ) :: kpref
    real(DP), intent( in ) :: ucVol
    logical, intent( in ) :: oldXASbroaden
    !
    integer :: ie, jdamp, jj, jdampStop
    real(DP), external :: gamfcn
    real(DP) :: e, gam, dr, di, ener, spct( 0 : 1 ), spkk( 0 : 1 )
    real(DP) :: alphaByVolumeHartree, sigma, tpa, tpacm4
    complex(DP) :: ctmp, disc, delta, rm1

    real(DP), parameter :: bohrSq2CmSq = bohr * bohr * 1.d-16
    !
    if( gamGauss .gt. (gam0/10.0_DP ) ) then
      call write_core_Voigt( fh, iter, kpref )
      return
    endif

    ! prefactor for comverting to cross-section. The energy will be in eV, so add that factor here
    alphaByVolumeHartree = eV2Hartree * alpha * ucVol 

    tpa = (alpha * alpha * 2.0_DP * PI_DP )* ucVol
    tpacm4 = tpa * bohrSq2CmSq * bohrSq2CmSq * au2sec
!    tpacm4 = tpa * ( bohr2cm**4 ) * au2sec

    jdampStop = 0
    if( oldXASbroaden ) then
      jdampStop = 1
    endif
    write( fh, '(A,1i5,A,1e15.8,A,1e15.8,A,1e15.8)' ) '#   iter=', iter, '   gam=', gam0, '   kpref=', kpref, 'vol=', ucVol
    write( fh, '(9(A15,1x))' ) '#   Energy', 'Spect', 'Spect(0)', 'SPKK', 'SPKK(0)', 'sigma (a.u.^2)', 'sigma( cm^2 )', &
               'sig.^2 (a.u.^5)', 'sig.^2 cm^4 sec'
    rm1 = -1; rm1 = sqrt( rm1 )
    do ie = 0, 2 * ne, 2
       e = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )
       do jdamp = 0, jdampStop
          gam= gam0 + gamfcn( e, nval, eps ) * dble( jdamp )
!          ctmp = e - a( iter - 1 ) + rm1 * gam

!          if( .true. ) then
#if(1)
          ctmp = cmplx( e - real_a( iter - 1 ), gam + imag_a( iter - 1 ), DP )  
          disc = sqrt( ctmp ** 2 - 4 * cmplx( real_b( iter ), imag_b( iter ) ) & 
                                     * cmplx( real_c( iter ), imag_c( iter ) ) )
          if( aimag( disc ) .gt. 0.0d0 ) then
            delta = (ctmp + disc ) / 2.0_dp
          else
            delta = (ctmp - disc ) / 2.0_dp
          endif

#else
            ctmp = e - real_a( iter - 1 ) + rm1 * gam
            disc = sqrt( ctmp ** 2 - 4 * real_b( iter ) ** 2 )
            di= -rm1 * disc
            if ( di .gt. 0.0d0 ) then
               delta = ( ctmp + disc ) / 2
            else
               delta = ( ctmp - disc ) / 2
            end if
#endif

          do jj = iter - 1, 0, -1
!             delta = e - a( jj ) + rm1 * gam - b( jj + 1 ) ** 2 / delta
!           if( .false. ) then
!             delta = e - real_a( jj ) + rm1 * gam - real_b( jj + 1 ) ** 2 / delta
!           else
            ctmp = cmplx( real_b( jj+1 ), imag_b( jj+1 ) ) * cmplx( real_c( jj+1 ), imag_c( jj+1 ) )
            delta = cmplx( e - real_a( jj ), gam + imag_a( jj ) ) - ctmp / delta
!           endif
          end do
          dr = delta
!          di = -rm1 * delta
!          di = abs( di )
          di = abs(aimag( delta ) )
          spct( jdamp ) = kpref * di / ( dr ** 2 + di ** 2 )
          spkk( jdamp ) = kpref * dr / ( dr ** 2 + di ** 2 )
       end do
       ! Making the broadening no longer defaults, but keeping the file format
       ! unchanged to minimize changes to the user
       if( .not. oldXASbroaden ) then
          spct(1) = spct(0)
          spkk(1) = spkk(0)
       endif
!       write ( fh, '(4(1e15.8,1x),1i5,1x,2(1e15.8,1x),1i5)' ) ener, spct( 1 ), spct( 0 ), spkk, iter, gam, kpref, ne
       
       ener = ebase + Hartree2eV * e
       sigma = alphaByVolumeHartree * ener * spct( 0 ) 
       write ( fh, '(9(1e15.8,1x))') ener, spct( 1 ), spct( 0 ), spkk( 1 ), spkk( 0 ), sigma, sigma*bohrSq2CmSq, &
              spct(0) * tpa, spct(0) * tpacm4
    end do
  end subroutine write_core

  subroutine write_core_Voigt( fh, iter, kpref ) 
    use OCEAN_constants, only : Hartree2eV
    implicit none
    integer, intent( in ) :: fh, iter
    real(DP), intent( in ) :: kpref
    !
    integer :: ie, jdamp, jj, ii, nstep, nratio, istart
    real(DP), external :: gamfcn
    real(DP) :: e, gam, dr, di, ener, spct( 0 : 1 ), spkk( 0 : 1 )
    complex(DP) :: ctmp, disc, delta, rm1
    !
    real(DP) :: estep, prefactor, inv2sigsquared, ee, pref, de
    real(DP), allocatable :: spct_array( :, :), spkk_array( :, : )

    de = ( eh - el ) / dble(ne)
    estep = min( gam0, gamGauss ) / 2.0_DP
    nratio = ( ( de - estep ) / estep ) + 1
    estep = de / real(nratio,DP)

    nstep = ( ( 2.99999_DP * gamGauss ) / estep ) + 1
    prefactor = estep / ( gamGauss * 2.506628274631_dp )
    inv2sigsquared = 1.0_DP / ( 2.0_DP * gamGauss * gamGauss )

    write(6,*) gam0, gamGauss, nstep, nratio
    allocate( spct_array( -nstep : nstep, 0:1 ), spkk_array( -nstep : nstep, 0:1 ) )

    write( fh, '(A,1i5,A,1e15.8,A,1e15.8)' ) '#   iter=', iter, '   gam=', gam0, '   kpref=', kpref
    write( fh, '(5(A15,1x))' ) '#   Energy', 'Spect', 'Spect(0)', 'SPKK', 'SPKK(0)'
    rm1 = -1; rm1 = sqrt( rm1 ) 
    do ie = 0, 2 * ne, 2
      spct(:) = 0.0_DP
      spkk(:) = 0.0_DP

      ee = el + ( eh - el ) * dble( ie ) / dble( 2 * ne ) 


      ! For each energy step of the output spectra (ie), we are going to integrate the over 
      ! neighboring energy points to include the Gaussian broadening using a finer grid with
      ! step size 'estep' instead of step size 'de'. Some (most) of this grid can be reused
      ! as we move from ie to ie + 1, except of course the first time.
      ! Here we set the bounds of the fine grid that need to be recalculated, and shift the 
      ! part that can be reused.
      if( ie .eq. 0 ) then
        istart = -nstep
      else
        istart = nstep - nratio
        do ii = -nstep, nstep - nratio
          spct_array( ii, 0 ) = spct_array( ii+nratio, 0 )
          spct_array( ii, 1 ) = spct_array( ii+nratio, 1 )
          spkk_array( ii, 0 ) = spkk_array( ii+nratio, 0 )
          spkk_array( ii, 1 ) = spkk_array( ii+nratio, 1 )
        enddo
      endif

      do ii = istart, nstep
        e = ee + real(ii,DP) * estep
      
        do jdamp = 0, 1
          gam= gam0 + gamfcn( e, nval, eps ) * dble( jdamp )
          ctmp = cmplx( e - real_a( iter - 1 ), gam + imag_a( iter - 1 ), DP )
          disc = sqrt( ctmp ** 2 - 4 * cmplx( real_b( iter ), imag_b( iter ) ) & 
                                     * cmplx( real_c( iter ), imag_c( iter ) ) )
          if( aimag( disc ) .gt. 0.0d0 ) then
            delta = (ctmp + disc ) / 2.0_dp
          else
            delta = (ctmp - disc ) / 2.0_dp
          endif
      
    
          do jj = iter - 1, 0, -1
            ctmp = cmplx( real_b( jj+1 ), imag_b( jj+1 ) ) * cmplx( real_c( jj+1 ), imag_c( jj+1 ) )
            delta = cmplx( e - real_a( jj ), gam + imag_a( jj ) ) - ctmp / delta
          end do
          dr = delta
          di = abs(aimag( delta ) )
          spct_array( ii, jdamp ) = prefactor * kpref * di / ( dr ** 2 + di ** 2 )
          spkk_array( ii, jdamp ) = prefactor * kpref * dr / ( dr ** 2 + di ** 2 )
        end do
      end do

      do ii = -nstep, nstep
        pref = exp( -real(ii**2,DP)* estep**2 * inv2sigsquared )
        spct( 0 ) = spct( 0 ) + pref * spct_array( ii, 0 )
        spct( 1 ) = spct( 1 ) + pref * spct_array( ii, 1 )
        spkk( 0 ) = spkk( 0 ) + pref * spkk_array( ii, 0 )
        spkk( 1 ) = spkk( 1 ) + pref * spkk_array( ii, 1 )
      enddo

      ener = ebase + Hartree2eV * e
      write ( fh, '(5(1e15.8,1x))') ener, spct( 1 ), spct( 0 ), spkk( 1 ), spkk( 0 )
    end do 

    deallocate( spct_array, spkk_array )
  end subroutine write_core_Voigt

  subroutine calc_spect_core( sp, iter, kpref )
    use OCEAN_constants, only : Hartree2eV, eV2Hartree
    implicit none
    real(DP), intent( out ) :: sp(:,:)
    integer, intent( in ) :: iter
    real(DP), intent( in ) :: kpref

    integer :: ie, jj
    real(DP) :: e, dr, di
    complex(DP) :: ctmp, disc, delta

    do ie = 1, ne
      e = el + ( eh - el ) * real( 2*(ie-1)+1, DP ) / real( 2 * ne, DP )
      ctmp = cmplx( e - real_a( iter - 1 ), gam0 + imag_a( iter - 1 ), DP )
      disc = sqrt( ctmp ** 2 - 4 * cmplx( real_b( iter ), imag_b( iter ) ) &
                                 * cmplx( real_c( iter ), -imag_c( iter ) ) )
      if( aimag( disc ) .gt. 0.0d0 ) then
        delta = (ctmp + disc ) / 2.0_dp
      else
        delta = (ctmp - disc ) / 2.0_dp
      endif

      do jj = iter -1, 0, -1
        ctmp = cmplx( real_b( jj+1 ), imag_b( jj+1 ) ) * cmplx( real_c( jj+1 ), -imag_c( jj+1 ) )
        delta = cmplx( e - real_a( jj ), gam0 + imag_a( jj ) ) - ctmp / delta
      enddo

      dr = delta
      di = abs( aimag( delta ) )
      sp(1,ie) = ebase + Hartree2eV * e
      sp(2,ie) = kpref * dr / ( dr ** 2 + di ** 2 )
      sp(3,ie) = kpref * di / ( dr ** 2 + di ** 2 )
    enddo

  

  end subroutine calc_spect_core

  subroutine calc_spect_val( sp, iter, kpref, celvol, nspin, backf )
    use OCEAN_constants, only : Hartree2eV, eV2Hartree
    implicit none
    real(DP), intent( out ) :: sp(:,:)
    integer, intent( in ) :: iter
    real(DP), intent( in ) :: kpref
    real(DP), intent( in ) :: celvol
    integer, intent( in ) :: nspin
    logical, intent( in ) :: backf

    integer :: ie, jj, i
    real(DP) :: e, dr, di, fact
    complex(DP) :: ctmp, disc, delta, arg, rp, rm , rrr, al, be, eps

    fact = kpref * real( 2 / nspin, DP ) * celvol

    do ie = 1, ne
      e = el + ( eh - el ) * real( 2*(ie-1)+1, DP ) / real( 2 * ne, DP )

      ctmp = cmplx( e, gam0, DP )
      if( backf ) then
        e = e**2
        ctmp = ctmp**2
      endif

      arg =  ( e - real_a( iter - 1 ) )**2 - 4.0_dp * real_b( iter ) ** 2 
      arg = sqrt(arg)

      rp = 0.5_dp * ( e - real_a( iter - 1 ) + arg )
      rm = 0.5_dp * ( e - real_a( iter - 1 ) - arg )
      if( aimag( rp ) .lt. 0.0_dp ) then
        rrr = rp
      else
        rrr = rm
      endif 
      al =  ctmp - real_a( iter - 1 ) - rrr
      be = -ctmp - real_a( iter - 1 ) - rrr

      do i = iter-1, 0, -1
        al = ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
           - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / al
        be = -ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
           - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / be
      enddo

      eps = 1.0_dp - fact / al - fact / be
      dr = real( eps, DP )
      di = aimag( eps ) 
      sp(1,ie) = ebase + Hartree2eV * e
      sp(2,ie) = dr
      sp(3,ie) = di
!      sp(2,ie) = fact * dr / ( dr ** 2 + di ** 2 )
!      sp(3,ie) = fact * di / ( dr ** 2 + di ** 2 )
    enddo

    

  end subroutine calc_spect_val

  subroutine OCEAN_haydock_setup( sys, ierr )
    use OCEAN_mpi
    use OCEAN_constants, only : Hartree2eV, eV2Hartree
    use OCEAN_system
    implicit none

    type(o_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer :: dumi, iter, ierr_
    character(len=4) :: inv_style
    real( DP ) :: dumf
    real( DP ), parameter :: default_gam0 = 0.1_DP

    gamGauss = sys%gaussBroaden * eV2Hartree

    if( .not. is_first ) goto 10
    is_first = .false.

    if( myid .eq. root ) then

      open(unit=99,file='bse.in',form='formatted',status='old')
      rewind(99)
      read(99,*) dumi
      read(99,*) dumf
      read(99,*) calc_type
!      select case ( calc_type )
!        case('hay')
          read(99,*) haydock_niter, ne, el, eh, gam0, ebase
          call checkBroadening( sys, gam0, default_gam0 )

!          el = el / 27.2114d0
!          eh = eh / 27.2114d0
!          gam0 = gam0 / 27.2114d0
          el = el * eV2Hartree
          eh = eh * eV2Hartree
          gam0 = gam0 * eV2Hartree
!          inv_loop = 1
!          allocate( e_list( inv_loop ) )
      close(99)


      open(unit=99,file='epsilon',form='formatted',status='old')
      rewind 99
      read(99,*) eps
      close(99)

      open( unit=99, file='nval.h', form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) nval
      close( unit=99 )
    endif

#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( haydock_niter, 1, MPI_INTEGER, root, comm, ierr )

!    call MPI_BCAST( nloop, 1, MPI_INTEGER, root, comm, ierr )
!    call MPI_BCAST( gres, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( gprc, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( ffff, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( ener, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( eps, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( nval, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )


!    call MPI_BCAST( e_start, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( e_stop, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( e_step, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( inv_loop, 1, MPI_INTEGER, root, comm, ierr )
!    if( myid .ne. root ) allocate( e_list( inv_loop ) )
!    call MPI_BCAST( e_list, inv_loop, MPI_DOUBLE_PRECISION, root, comm, ierr )


    call MPI_BCAST( echamp, 1, MPI_LOGICAL, root, comm, ierr )
    call MPI_BCAST( project_absspct, 1, MPI_LOGICAL, root, comm, ierr )
#endif

10 continue

!    if( allocated( a ) ) deallocate( a )
!    if( allocated( b ) ) deallocate( b )
    if( allocated( real_a ) ) deallocate( real_a )
    if( allocated( imag_a ) ) deallocate( imag_a )
    if( allocated( real_b ) ) deallocate( real_b )
    if( allocated( imag_b ) ) deallocate( imag_b )
    if( allocated( real_c ) ) deallocate( real_c )
    if( allocated( imag_c ) ) deallocate( imag_c )
    if( haydock_niter .gt. 0 ) then
!      allocate( a( 0 : haydock_niter ) )
!      allocate( b( 0 : haydock_niter ) )
      allocate( real_a( 0 : haydock_niter ) )
      allocate( imag_a( 0 : haydock_niter ) )
      allocate( real_b( 0 : haydock_niter ) )
      allocate( imag_b( 0 : haydock_niter ) )
      allocate( real_c( 0 : haydock_niter ) )
      allocate( imag_c( 0 : haydock_niter ) )
!      a(:) = 0.0_DP
!      b(:) = 0.0_DP
      real_a(:) = 0.0_DP
      imag_a(:) = 0.0_DP
      real_b(:) = 0.0_DP
      imag_b(:) = 0.0_DP
      real_c(:) = 0.0_DP
      imag_c(:) = 0.0_DP
!    else
!      allocate( a(1), b(1) )
    endif

  end subroutine OCEAN_haydock_setup


  subroutine checkBroadening( sys, broaden, default_broaden )
    use OCEAN_corewidths, only : returnLifetime
    use OCEAN_system
    use OCEAN_mpi, only : myid, root
    implicit none
    type( o_system ), intent( in ) :: sys
    real(DP), intent( inout ) :: broaden
    real(DP), intent( in ) :: default_broaden
    
    if( broaden .gt. 0.0_dp ) return


    select case ( sys%cur_run%calc_type )
      case( 'VAL' )
        broaden = default_broaden 
      case( 'XAS' , 'XES', 'RXS' ) 
        call returnLifetime( sys%cur_run%ZNL(1), sys%cur_run%ZNL(2), sys%cur_run%ZNL(3), broaden )
        if( broaden .le. 0 ) broaden = default_broaden
        if( myid .eq. root ) write(6,*) 'Default broadening used: ', broaden
      case default
        broaden = default_broaden
    end select
    write(6,*) 'Default requested for broadening: ', broaden

    end subroutine checkBroadening

  subroutine redtrid(n,sys, kpref, ierr)
    use OCEAN_system, only : o_system
    use OCEAN_filenames, only : OCEAN_filenames_lanc
    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) ::  n
    real(DP), intent( in ) :: kpref

    double precision, allocatable :: ar(:,:),ai(:,:)
    double precision, allocatable :: w(:),zr(:,:),zi(:,:)
    double precision, allocatable :: fv1(:),fv2(:),fm1(:)
    integer :: matz,nm,i,j,nn


    nn=n+1
    nm=nn+10
    allocate(ar(nm,nm),ai(nm,nm),w(nm),zr(nm,nm),zi(nm,nm))
    allocate(fv1(nm),fv2(nm),fm1(2*nm))
    do i=1,n+1
      do j=1,n+1
        ar(j,i)=0.d0
        ai(j,i)=0.d0
      end do
    end do
    do i=1,n+1
      ar(i,i)=real_a(i-1)
      ai(i,i)=imag_a(i-1)
      if (i.le.n) then
        ar(i+1,i)=real_c(i)
        ar(i,i+1)=real_b(i)
        ai(i+1,i)=imag_c(i)
        ai(i,i+1)=imag_b(i)
      end if
    end do
    matz=0
    call elsch(nm,nn,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
    if( n .lt. 0 ) return

    call write_lanczos( n, sys, kpref, ierr, w )

    deallocate(ar,ai,w,zr,zi,fv1,fv2,fm1)
    return
  end subroutine redtrid

  subroutine write_lanczos( n, sys, kpref, ierr, w )
    use OCEAN_system, only : o_system
    use OCEAN_filenames, only : OCEAN_filenames_lanc
    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) ::  n
    real(DP), intent( in ) :: kpref
    real(DP), intent( in ), optional :: w(:)

    integer :: i
    character( LEN=40 ) :: lanc_filename

    call OCEAN_filenames_lanc( sys, lanc_filename, ierr )
    if( ierr .ne. 0 ) return

    open(unit=99,file=lanc_filename,form='formatted',status='unknown')
    rewind 99
    write ( 99, '(1i8,1x,1ES24.17)' ) n, kpref
    if( complex_haydock ) then
      write ( 99, '(2(2x,ES24.17))' ) real_a( 0 ), imag_a( 0 )
      do i = 1, n
        write ( 99, '(2x,6ES24.17)' ) real_a( i ), imag_a( i ), real_b( i ), imag_b( i ), &
                                     real_c( i ), imag_c( i )
      enddo
    else
      do i = 0, n
        if ( i .eq. 0 ) then
          write ( 99, '(2x,ES24.17)' ) real_a( i )
        else
          write ( 99, '(2(2x,ES24.17))' ) real_a( i ), real_b( i )
        end if
      end do
    endif

  
    if( present( w ) ) then
      write (99,'(2x,2i5,1f20.10)') (i,n+1,w(i),i=1,n+1)
    endif
    close(unit=99)
  end subroutine write_lanczos



  subroutine testConvergeEps( iter, sys, kpref, ucvol, val_ham_spin, restartBSE, newEps )
    use OCEAN_system, only : o_system
    use OCEAN_constants, only : PI_DP
    implicit none
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: iter, val_ham_spin
    real(DP), intent( in ) :: kpref, ucvol
    logical, intent( inout ) :: restartBSE
    real(DP), intent( out ) :: newEps

    complex(DP) :: ctmp, arg, rrr, rp, rm, al, be, eps
    real(DP) :: tcEps, fact, oldEps, epsErr
    integer :: i

    fact = kpref * real( 2 / val_ham_spin, DP ) * ucvol
    ctmp = cmplx( 0, gam0, DP )

    arg = real_a( iter - 1 )** 2 - 4.0_dp * real_b( iter ) ** 2
    arg = sqrt( arg )
  
    rp = 0.5_dp * ( - real_a( iter - 1 ) + arg )
    rm = 0.5_dp * ( - real_a( iter - 1 ) - arg )
    if( aimag( rp ) .lt. 0.0_dp ) then
      rrr = rp
    else
      rrr = rm
    endif

    al = ctmp - real_a( iter - 1 ) - rrr
    be = -ctmp - real_a( iter - 1 ) - rrr

    do i = iter-1, 0, -1
      al = ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
         - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / al
      be = -ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
         - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / be
    enddo

    eps = 1.0_dp - fact / al - fact / be

    tcEps = abs( dble( eps ) )

    if( iter .lt. 3 ) then
      eps1Conv( iter ) = tcEps
!      write(6,*) 'Estimated eps1(0): ', tcEps
    else
      eps1Conv( 1 ) = eps1Conv( 2 )
      eps1Conv( 2 ) = eps1Conv( 3 )
      eps1Conv( 3 ) = tcEps

      tcEps = sum(eps1Conv(:)) / 3.0_DP
      if( max( sys%epsilon0, eps1Conv( 3 ) ) .lt. 100.0d0 ) then
        write(6,'(3(A,F9.4,1X))') 'Est. eps1(0): ', eps1Conv( 3 ), ';  Avg: ', tcEps, &
                                  ';  Current: ', sys%epsilon0
      else
        write(6,'(3(A,E24.12,1X))') 'Est. eps1(0): ', eps1Conv( 3 ), ';  Avg: ', tcEps, &
                                    ';  Current: ', sys%epsilon0
      endif

      ! change to percentage
      if( ( maxval(eps1Conv(:)) - minval(eps1Conv(:)) )/tcEps .gt. 0.05_dp ) return
      if( abs( sys%epsilon0 - tcEps ) / ( sys%epsilon0 + tcEps - 2.0_dp ) &
                  .lt. 10.0_DP * sys%epsConvergeThreshold ) then
        if( ( maxval(eps1Conv(:)) - minval(eps1Conv(:)) )/tcEps .gt. 0.5_dp * sys%epsConvergeThreshold ) then
!            write(6,*) 'C', ( maxval(eps1Conv(:)) - minval(eps1Conv(:)) )/tcEps
            return
        endif
      endif
          

      if( abs( sys%epsilon0 - tcEps ) .gt. ( maxval(eps1Conv(:)) - minval(eps1Conv(:)) ) .and. &
          abs( sys%epsilon0 - tcEps ) / ( sys%epsilon0 + tcEps - 2.0_dp ) & 
                  .gt. 0.5_DP * sys%epsConvergeThreshold ) then  
        newEps = ( 2.0_DP * eps1Conv( 3 ) + eps1Conv(2) ) / 3.0_DP
#if 0
        if( abs( newEps - sys%epsilon0 )/( newEps + sys%epsilon0 ) .gt. 0.15_dp ) then
          newEps = 0.95_dp * newEps + 0.05_DP * sys%epsilon0
        elseif( abs( newEps - sys%epsilon0 )/( newEps + sys%epsilon0 ) .gt. 0.02_dp ) then
          newEps = 0.98_dp * newEps + 0.02_DP * sys%epsilon0
        endif
#else
        epsErr = 20.0_dp * abs( newEps - sys%epsilon0 )/( newEps + sys%epsilon0 - 2.0_dp )
        epsErr = ( 0.2_DP / PI_DP ) * atan( epsErr )
        write(6,*) newEps, sys%epsilon0, epsErr
        newEps = (1.0_DP-epsErr)*newEps + epsErr * sys%epsilon0
#endif
        
        
        restartBSE = .true.
        write(6,*) 'Restart eps1(0): ', newEps, sys%epsilon0
        return
      endif
    endif

  end subroutine testConvergeEps

  subroutine write_val_tri( fh, iter, kpref , ucvol, val_ham_spin, semiTDA, backf, overlaps, ierr )
    use OCEAN_constants, only : Hartree2eV, bohr, alphainv
    implicit none
    integer, intent( in ) :: fh, iter, val_ham_spin
    real(DP), intent( in ) :: kpref, ucvol
    logical, intent( in ) :: semiTDA, backf
    complex(DP), intent( in ) :: overlaps( iter )
    integer, intent( inout ) :: ierr
    !
    integer :: ie, i
    real(DP) :: ere, reeps, imeps, lossf, fact, mu, reflct
    complex(DP) :: ctmp, arg, rp, rm, rrr, al, be, eps, refrac
    complex(DP), allocatable :: tmp_d(:), tmp_du(:), tmp_dl(:), tmp_b(:)


    write(6,*) iter
    allocate( tmp_d( iter), tmp_b(iter), tmp_du(iter-1), tmp_dl(iter-1) )

    fact = kpref * real( 2 / val_ham_spin, DP ) * ucvol

    write(fh,"(a)") "#   omega (eV)      epsilon_1       epsilon_2       n"// &
      "               kappa           mu (cm^(-1))    R"//  &
      "               epsinv"

    do ie = 1, 2 * ne, 2
      ere = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )

      tmp_b(:) = 0.0_DP
      tmp_b(1) = 1.0_DP
      do i = 1, iter
        tmp_d( i ) = cmplx( ere, gam0, DP ) - cmplx(real_a(i-1), imag_a(i-1), DP )
      enddo
      do i = 1, iter -1
        tmp_du( i ) = -cmplx( real_b(i), imag_b(i) )
        tmp_dl( i ) = -cmplx( real_c(i), imag_c(i) )
      enddo

      call ZGTSV( iter, 1, tmp_dl, tmp_d, tmp_du, tmp_b, iter, ierr )
      if( ierr .ne. 0 ) then
        write(6,*) 'Failed at ie, ere:', ie, ere, ierr
        return
      endif
    
      ctmp = fact * dot_product( overlaps, tmp_b )

      if( semiTDA ) then
        tmp_b(:) = 0.0_DP
        tmp_b(1) = 1.0_DP
        do i = 1, iter
          tmp_d( i ) = -cmplx( ere, gam0, DP ) - cmplx(real_a(i-1), imag_a(i-1), DP )
        enddo
        do i = 1, iter -1
          tmp_du( i ) = -cmplx( real_b(i), imag_b(i) )
          tmp_dl( i ) = -cmplx( real_c(i), imag_c(i) )
        enddo

        call ZGTSV( iter, 1, tmp_dl, tmp_d, tmp_du, tmp_b, iter, ierr )
        if( ierr .ne. 0 ) then
          write(6,*) 'Failed at ie, ere:', ie, ere, ierr
          return
        endif
        
        ctmp = ctmp + fact * dot_product( overlaps, tmp_b )
      endif
      eps = 1.0_DP - ctmp
    
      reeps = dble( eps )
      imeps = aimag( eps )
      lossf = imeps / ( reeps ** 2 + imeps ** 2 )
      refrac = sqrt(eps)
      reflct = abs((refrac-1.0d0)/(refrac+1.0d0))**2
      mu = 2.0d0 * ere * Hartree2eV * aimag(refrac) / ( bohr * alphainv * 1000 )

      write(fh,'(8(1E24.16,1X))') ere*Hartree2eV, reeps, imeps, refrac-1.0d0, mu, reflct, lossf

    enddo

    deallocate( tmp_d, tmp_du, tmp_dl, tmp_b )
  end subroutine write_val_tri

end module OCEAN_haydock
