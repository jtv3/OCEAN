! Copyright (C) 2018 - 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
!
! This is beginning as a straight up copy of the QE54 file.
module ocean_qe62_files
  use ai_kinds, only : DP
#ifdef MPI_F08
  use mpi_f08, only : MPI_COMM, MPI_REQUEST
#endif


  implicit none
  private
  save


  logical :: is_init = .false.
  logical :: is_gamma
  logical, parameter :: GammaFullStorage = .false.
  character( len=128 ) :: prefix

  integer :: bands(2)
  integer :: brange(4)
  integer :: kpts(3)
  integer :: nspin 
  integer :: nfiles


#ifdef MPI_F08
  type( MPI_COMM ) :: inter_comm
  type( MPI_COMM ) :: pool_comm
#else
  integer :: inter_comm
  integer :: pool_comm
#endif

  integer :: inter_myid
  integer :: inter_nproc
  integer, parameter :: inter_root = 0

  integer :: npool
  integer :: mypool
  integer :: pool_myid
  integer :: pool_nproc
  integer, parameter :: pool_root = 0
  
  integer :: pool_nbands
  integer :: pool_val_nbands
  integer :: pool_con_nbands

  public :: qe62_read_init, qe62_read_at_kpt, qe62_clean, qe62_get_ngvecs_at_kpt, &
            qe62_read_energies_single
  public :: qe62_kpts_and_spins, qe62_return_my_bands, qe62_return_my_val_bands, qe62_return_my_con_bands, &
            qe62_is_my_kpt
  public :: qe62_nprocPerPool, qe62_getPoolIndex, qe62_returnGlobalID
  public :: qe62_getAllBandsForPoolID, qe62_getValenceBandsForPoolID, qe62_getConductionBandsForPoolID
  public :: qe62_npool, qe62_universal2KptAndSpin, qe62_poolComm, qe62_poolID

  contains

  pure function qe62_poolComm()
#ifdef MPI_F08
    type( MPI_COMM ) :: qe62_poolComm
#else
    integer :: qe62_poolComm
#endif

    qe62_poolComm = pool_comm
  end function qe62_poolComm

  pure function qe62_npool() result( npool_ )
    integer :: npool_
    npool_ = npool
  end function qe62_npool

  pure function qe62_nprocPerPool() result( nproc )
    integer :: nproc

    nproc = pool_nproc
  end function qe62_nprocPerPool

  pure function qe62_poolID() result( pid )
    integer :: pid

    pid = pool_myid
  end function qe62_poolID

  pure function qe62_getPoolIndex( ispin, ikpt ) result( poolIndex )
    integer, intent( in ) :: ispin, ikpt
    integer :: poolIndex
    integer :: kptCounter
    !
    kptCounter = ikpt + ( ispin - 1 ) * product(kpts(:))
    poolIndex = mod( kptCounter, npool )
  end function qe62_getPoolIndex

  pure function qe62_getAllBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

    bands_remain = bands(2)-bands(1)+1

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function qe62_getAllBandsForPoolID

  pure function qe62_getValenceBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

!    bands_remain = brange(4)-brange(3)+brange(2)-brange(1)+2
    bands_remain = brange(2) - brange(1) + 1

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function qe62_getValenceBandsForPoolID

  pure function qe62_getConductionBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i
    
!    bands_remain = brange(4)-brange(3)+brange(2)-brange(1)+2
    bands_remain = brange(4) - brange(3) + 1
  
    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function qe62_getConductionBandsForPoolID

  pure function qe62_returnGlobalID( poolIndex, poolID ) result( globalID )
    integer, intent( in ) :: poolIndex, poolID
    integer :: globalID

    globalID = poolIndex * pool_nproc + poolID
  end function qe62_returnGlobalID


  subroutine qe62_universal2KptAndSpin( uni, ikpt, ispin )
    integer, intent( in ) :: uni
    integer, intent( out ) :: ispin, ikpt
    !
    integer :: i, ierr
    logical :: is_kpt

    i = 0
    do ispin = 1, nspin
      do ikpt = 1, product(kpts(:))
        call qe62_is_my_kpt( ikpt, ispin, is_kpt, ierr )
        if( is_kpt ) i = i + 1
        if( uni .eq. i ) return
      enddo
    enddo

    ikpt = 0
    ispin = 0
  end subroutine qe62_universal2KptAndSpin

  subroutine qe62_is_my_kpt( ikpt, ispin, is_kpt, ierr )
    integer, intent( in ) :: ikpt, ispin
    logical, intent( out ) :: is_kpt
    integer, intent( inout ) :: ierr
    !
    if( .not. is_init ) then
      ierr = 2
      return
    endif
    if( ispin .lt. 1 .or. ispin .gt. nspin ) then
      ierr = 3 
      return
    endif
    if( ikpt .lt. 1 .or. ikpt .gt. product(kpts(:)) ) then
      ierr = 4
      return
    endif

!    if( mod( ikpt + (ispin-1)*product(kpts(:)), npool ) .eq. mypool ) then
    if( qe62_getPoolIndex( ispin, ikpt ) .eq. mypool ) then
      is_kpt = .true.
    else
      is_kpt = .false.
    endif

  end subroutine qe62_is_my_kpt

  subroutine qe62_return_my_bands( nbands, ierr )
    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr
    
    if( .not. is_init ) ierr = 1
    nbands = pool_nbands
  end subroutine qe62_return_my_bands

  subroutine qe62_return_my_val_bands( nbands, ierr )
    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr

    if( .not. is_init ) ierr = 1
    nbands = pool_val_nbands
  end subroutine qe62_return_my_val_bands

  subroutine qe62_return_my_con_bands( nbands, ierr )
    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr

    if( .not. is_init ) ierr = 1
    nbands = pool_con_nbands
  end subroutine qe62_return_my_con_bands

  integer function qe62_kpts_and_spins()
    qe62_kpts_and_spins = product(kpts(:)) * nspin / npool
    return
  end function qe62_kpts_and_spins

  subroutine qe62_read_energies_single( myid, root, comm, energies, ierr )
#ifdef MPI
    use ocean_mpi, only : MPI_DOUBLE_PRECISION
#endif
    integer, intent( in ) :: myid, root
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    real( DP ), intent( out ) :: energies( :, :, : )
    integer, intent( inout ) :: ierr
    !
    integer :: ispn, ikpt, i, bstop2 !bstop1, bstart2, bstop2
    real(DP) :: dumf
    character(len=128) :: lineBurn

    if( is_init .eqv. .false. ) then
      write(6,*) 'qe62_read_energies_single called but was not initialized'
      ierr = 4
      return
    endif 

    if( size( energies, 1 ) .lt. ( bands(2)-bands(1)+1 ) ) then
      write(6,*) 'Error! Band dimension of energies too small in qe62_read_energies_single'
      ierr = 1
      return
    endif
    if( size( energies, 2 ) .lt. product(kpts(:)) ) then
      write(6,*) 'Error! K-point dimension of energies too small in qe62_read_energies_single'
      ierr = 2
      return
    endif
    if( size( energies, 3 ) .lt. nspin ) then
      write(6,*) 'Error! Spin dimension of energies too small in qe62_read_energies_single'
      ierr = 3
      return
    endif

    if( myid .eq. root ) then

      open( unit=99, file='enkfile', form='formatted', status='old' )

      do ispn = 1, nspin
        do ikpt = 1, product(kpts(:))
          read(99,*) energies( :, ikpt, ispn )
        enddo
      enddo

      close( 99 )
      energies(:,:,:) = energies(:,:,:) * 2.0_DP
    endif

#ifdef MPI
    call MPI_BCAST( energies, size(energies), MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

  end subroutine qe62_read_energies_single

  pure function qe62_wfcFile( ikpt, ispin ) result( fileName )
    integer, intent ( in ) :: ikpt, ispin

    integer :: totKpt
    character(len=512) :: fileName

    if( ispin .gt. 1 ) then
      totKpt = ikpt + product( kpts(:) )
    else
      totKpt = ikpt
    endif

    write( fileName, '(a,a3,i0,a4)' ) trim( prefix ), 'wfc', totKpt, '.dat'

  end function qe62_wfcFile

#if 0
  pure function qe62_eigFile( ikpt, ispin ) result( fileName )
    integer, intent ( in ) :: ikpt, ispin
    character(len=512) :: fileName

    if( nspin .eq. 1 ) then
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'eigenval.xml' 
    elseif( ispin .eq. 1 ) then
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'eigenval1.xml'
    else
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'eigenval2.xml'
    endif

  end function qe62_eigFile

  pure function qe62_gkvFile( ikpt, ispin ) result( fileName )
    integer, intent ( in ) :: ikpt, ispin
    character(len=512) :: fileName

    write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'gkvectors.dat'

  end function qe62_gkvFile

  pure function qe62_evcFile( ikpt, ispin ) result( fileName )
    integer, intent ( in ) :: ikpt, ispin
    character(len=512) :: fileName

    if( nspin .eq. 1 ) then
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'evc.dat'
    elseif( ispin .eq. 1 ) then
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'evc1.dat'
    else
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'evc2.dat'
    endif

  end function qe62_evcFile
#endif

#if 0
  subroutine qe62_read_energies( myid, root, comm, nbv, nbc, nkpts, nspns, val_energies, con_energies, ierr )
#ifdef MPI
    use ocean_mpi, only : MPI_DOUBLE_PRECISION
#endif
    integer, intent( in ) :: myid, root, nbv, nbc, nkpts, nspns
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    real( DP ), intent( out ) :: val_energies( nbv, nkpts, nspns ), con_energies( nbv, nkpts, nspns )
    integer, intent( inout ) :: ierr
    !
    integer :: ispn, ikpt

    if( myid .eq. root ) then
      open( unit=99, file='enkfile', form='formatted', status='old' )
      do ispn = 1, nspns
        do ikpt = 1, nkpts
          read( 99, * ) val_energies( :, ikpt, ispn )
          read( 99, * ) con_energies( :, ikpt, ispn )
        enddo
      enddo
      close( 99 )
    endif

#ifdef MPI
    call MPI_BCAST( val_energies, nbv*nkpts*nspns, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( con_energies, nbc*nkpts*nspns, MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

  end subroutine qe62_read_energies
#endif


  subroutine qe62_clean( ierr )
    integer, intent( inout ) :: ierr
    !
    nfiles = 0

#ifdef MPI
    call MPI_COMM_FREE( pool_comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_FREE( inter_comm, ierr )
    if( ierr .ne. 0 ) return
#endif

  end subroutine qe62_clean

  ! Read the universal little files
  subroutine qe62_read_init( comm, isGamma, isFullStorage, ierr )
    use ocean_mpi, only : MPI_INTEGER, MPI_CHARACTER, MPI_LOGICAL
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    logical, intent( out ) :: isGamma, isFullStorage
    integer, intent( inout ) :: ierr
    !
    integer :: i
    logical :: ex
    character(len=128) :: tmp
    !
    ! Set the comms for file handling
    call MPI_COMM_DUP( comm, inter_comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_SIZE( inter_comm, inter_nproc, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_RANK( inter_comm, inter_myid, ierr )
    if( ierr .ne. 0 ) return


    if( inter_myid .eq. inter_root ) then
      open( unit=99, file='prefix', form='formatted', status='old' )
      read(99,*)  tmp
      close( 99 )
      write( prefix, '(a,a,a)' ) 'Out/', trim(tmp), '.save/'

      open(unit=99,file='brange.ipt', form='formatted', status='old' )
      read(99,*) brange(:)
      close(99)

      inquire( file='bands.ipt', exist=ex )
      if( ex ) then
        open(unit=99,file='bands.ipt', form='formatted', status='old' )
        read(99,*) bands(:)
        close(99)
      else
!        open(unit=99,file='brange.ipt', form='formatted', status='old' )
!        read(99,*) brange(:)
!        close(99)
        bands(2) = brange(4)
        bands(1) = brange(1)
      endif
      open(unit=99,file='kmesh.ipt', form='formatted', status='old' )
      read(99,*) kpts(:)
      close(99)
      open(unit=99,file='nspin', form='formatted', status='old' )
      read(99,*) nspin
      close(99)

      inquire( file='gamma', exist=ex )
      if( ex ) then
        open(unit=99,file='gamma', form='formatted', status='old' )
        read(99,*) is_gamma
        close(99)
      else
        is_gamma = .false.
      endif

      nfiles = nspin * product(kpts(:) )
    endif

#ifdef MPI
    call MPI_BCAST( nfiles, 1, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( bands, 2, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( kpts, 3, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( nspin, 1, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( prefix, len(prefix), MPI_CHARACTER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( is_gamma, 1, MPI_LOGICAL, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return
#endif

    call set_pools( ierr )
    if( ierr .ne. 0 ) return

    call writeDiagnostics( )

!    write(6,*) 'qe62_read_init was successful'
    is_init = .true.
    isGamma = is_gamma

    if( isGamma ) then
      isFullStorage = gammaFullStorage
    else
      isFullStorage = .true.
    endif
  
  end subroutine  qe62_read_init 

  subroutine writeDiagnostics( )
    if( inter_myid .eq. inter_root ) then
      write( 6, '(A)' ) "    #############################"
      write( 6, '(A,I8)' ) " Npools:      ", npool
      write( 6, '(A,I8)' ) " Nprocs/pool: ", pool_nproc
    endif

    write(1000+inter_myid, '(A)' ) "    #############################"
    write(1000+inter_myid, '(A,I8)' ) " Npools:      ", npool
    write(1000+inter_myid, '(A,I8)' ) " Nprocs/pool: ", pool_nproc
    write(1000+inter_myid, '(A,I8)' ) " My pool:     ", mypool
    write(1000+inter_myid, '(A,I8)' ) " My pool id:  ", pool_myid
    write(1000+inter_myid, '(A,I8)' ) " My bands:    ", pool_nbands


    write(1000+inter_myid, '(A)' ) "    #############################"
    flush(1000+inter_myid)
  end subroutine 

  subroutine set_pools( ierr )
    integer, intent( inout ) :: ierr
    !
    integer :: i

    if( nfiles .ge. inter_nproc ) then
      mypool = inter_myid
      npool = inter_nproc

    else
      do i = 2, inter_nproc
        if( mod( inter_nproc, i ) .eq. 0 ) then
          if( inter_myid .eq. 0 ) write(6,*) i, inter_nproc
          write(1000+inter_myid,*)  i, inter_nproc, nfiles
          if( nfiles .ge. (inter_nproc/i) ) then
            npool = inter_nproc/i
            mypool = inter_myid/ i 
#ifdef DEBUG
            write(6,*) '*', inter_myid, npool, mypool
#endif
            goto 11
          endif
        endif
      enddo
11    continue
    endif
    
    call MPI_COMM_SPLIT( inter_comm, mypool, inter_myid, pool_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_COMM_RANK( pool_comm, pool_myid, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_SIZE( pool_comm, pool_nproc, ierr )
    if( ierr .ne. 0 ) return


    pool_nbands = qe62_getAllBandsForPoolID( pool_myid )
!    nbands_left = brange(4)-brange(3)+brange(2)-brange(1)+2
!    do i = 0, pool_nproc-1
!      nbands = nbands_left / ( pool_nproc - i )
!      if( i .eq. pool_myid ) pool_nbands = nbands
!      nbands_left = nbands_left - nbands
!    enddo

  end subroutine set_pools


  subroutine qe62_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
    use OCEAN_mpi
    integer, intent( in ) :: ikpt, ispin
    integer, intent( out ) :: gvecs
    integer, intent( inout ) :: ierr
    !
    integer :: i, crap
    logical :: is_kpt

    call qe62_is_my_kpt( ikpt, ispin, is_kpt, ierr )
    if( ierr .ne. 0 ) return
    if( is_kpt .eqv. .false. ) then 
      gvecs = 0
      return
    endif

    if( pool_myid .eq. pool_root ) then
      write(myid+1000,*) 'Opening file', trim(qe62_wfcFile( ikpt, ispin ))
      flush(myid+1000)
      open( unit=99,file=trim(qe62_wfcFile( ikpt, ispin )), form='unformatted', status='old' )
      read( 99 )
      read(99) crap, gvecs
      close( 99 )
      ! Do we want to expand the coeffs inside this module?
      if( gammaFullStorage .and. is_gamma ) gvecs = 2 * gvecs - 1
    endif

!    write(6,*) 'gvecs', pool_root, pool_comm
111 continue
    
#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, i )
    if( ierr .ne. 0 ) return
    if( i .ne. 0 ) then
      ierr = i
      return
    endif

    call MPI_BCAST( gvecs, 1, MPI_INTEGER, pool_root, pool_comm, ierr )
    if( ierr .ne. 0 ) return
#endif
!    write(6,*) 'gvecs', pool_root, pool_comm

  end subroutine qe62_get_ngvecs_at_kpt

  subroutine qe62_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
#ifdef MPI
    use OCEAN_mpi, only : MPI_INTEGER, MPI_DOUBLE_COMPLEX, MPI_STATUSES_IGNORE, myid, MPI_REQUEST_NULL, &
                          MPI_STATUS_IGNORE, MPI_UNDEFINED
!    use OCEAN_mpi
#endif
    use SCREEN_timekeeper, only : SCREEN_tk_start, SCREEN_tk_stop
    integer, intent( in ) :: ikpt, ispin, ngvecs, my_bands
    integer, intent( out ) :: gvecs( 3, ngvecs )
    complex( DP ), intent( out ) :: wfns( ngvecs, my_bands )
    integer, intent( inout ) :: ierr
    !
    complex( DP ), allocatable, dimension( :, :, : ) :: cmplx_wvfn
    real( DP ) :: bvecs(3,3), kpt(3), scal
    integer :: test_gvec, itarg, nbands_to_send, nr, ierr_, nbands, id, start_band, crap, i, j, k, bufferSize, maxBands
    integer :: ikpt_, ispin_, dumi, npolar, my_bands_
    logical :: is_gamma_
#ifdef MPI_F08
    type( MPI_REQUEST ), allocatable :: requests( : )
    type( MPI_DATATYPE ) :: newType
#else
    integer, allocatable :: requests( : )
    integer :: newType
#endif

    if( qe62_getPoolIndex( ispin, ikpt ) .ne. mypool ) return
    !
    if( qe62_getAllBandsForPoolID( pool_myid ) .ne. my_bands ) then
      ierr = 1
      return
    endif

    nbands = bands(2)-bands(1)+1

    if( pool_myid .eq. pool_root ) then
    
      ! get gvecs first
      open( unit=99, file=qe62_wfcFile( ikpt, ispin), form='unformatted', status='old' )
      read( 99 ) ikpt_, kpt, ispin_, is_gamma_, scal
      read( 99 ) dumi, test_gvec, npolar, my_bands_
      ! Are we expanding the wave functions here?
      if( gammaFullStorage .and. is_gamma ) then
        if( ( 2 * test_gvec - 1 ) .ne. ngvecs ) then
          ierr = -2
          write(6,*) (2*test_gvec-1), ngvecs
          return
        endif
      else
        if( test_gvec .ne. ngvecs ) then
          ierr = -2
          write(6,*) test_gvec, ngvecs
          return
        endif
      endif

      write(6,*) dumi, test_gvec, npolar, my_bands_

      ! Error synch. Also ensures that the other procs have posted their recvs
      call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, ierr_ )
      if( ierr .ne. 0 ) return
      if( ierr_ .ne. 0 ) then
        ierr = ierr_
        return
      endif

      if( npolar .ne. 1 ) then
        write(6,*) 'Polarization ne 1 not yet allowed'
        ierr = 10
        return
      endif

      ! Using test_gvec for gamma support
      read( 99 ) bvecs( :, : )
      read( 99 ) gvecs( :, 1:test_gvec )

      maxBands = qe62_getAllBandsForPoolID( 0 )
      do id = 1, pool_nproc - 1
        maxBands = max( maxBands, qe62_getAllBandsForPoolID( id ) )
      enddo
      !
      ! no more than 8GB @ 16Byte/complex
!      bufferSize = floor( 36870912.0_DP /  ( real( test_gvec, DP ) * real( maxBands, DP ) ) )
      bufferSize = floor( 536870912.0_DP /  ( real( test_gvec, DP ) * real( maxBands, DP ) ) )
!      bufferSize = floor( 1073741824.0_DP /  ( real( test_gvec, DP ) * real( maxBands, DP ) ) )
      bufferSize = max( bufferSize, 1 )
      bufferSize = min( bufferSize, pool_nproc - 1 )
      allocate( cmplx_wvfn( test_gvec, maxBands, bufferSize ) )

!      nr = pool_nproc 
      nr = bufferSize
!      nr = 2 * pool_nproc
      allocate( requests( 0:nr ) )
      requests(:) = MPI_REQUEST_NULL
#ifdef MPI
      call MPI_IBCAST( gvecs, 3*test_gvec, MPI_INTEGER, pool_root, pool_comm, requests( 0 ), ierr )
#endif

      write(1000+myid,*) '***Reading k-point: ', ikpt, ispin
      write(1000+myid,*) '   Ngvecs: ', ngvecs
      write(1000+myid,*) '   Nbands: ', nbands
      write(1000+myid,*) '   Nbuffer:', bufferSize

      ! if we are skipping bands, do that here
!      do i = 1, bands(1)-1
!        read( 99 ) crap, cmplx_wvfn( 1, 1 )
!      enddo



      start_band = 1

#ifdef MPI
      ! loop over each proc in this pool to send wavefunctions
      do id = 0, pool_nproc - 1
        nbands_to_send = qe62_getAllBandsForPoolID( id )

        ! If mine, read directly to wvfn
        if( id .eq. pool_myid ) then
          call SCREEN_tk_start("dft-read")
          do i = 1, nbands_to_send
            read( 99 ) wfns( 1:test_gvec * npolar, i )
          enddo
          call SCREEN_tk_stop("dft-read")
        
        ! IF not mine, find open buffer, read to buffer, send buffer
        else
          j = 0
          do i = 1, bufferSize
            if( requests( i ) .eq. MPI_REQUEST_NULL ) then
              j = i
              exit
            endif
          enddo
          if( j .eq. 0 ) then
            call MPI_WAITANY( bufferSize, requests(1:bufferSize), j, MPI_STATUS_IGNORE, ierr )
            write(1000+myid, * ) 'Waited for buffer:', j
            if( j .eq. MPI_UNDEFINED ) j = 1
          endif
          

          call SCREEN_tk_start("dft-read")
          do i = 1, nbands_to_send
            read( 99 ) cmplx_wvfn( 1:test_gvec, i, j )
          enddo
          call SCREEN_tk_stop("dft-read")

          write(1000+myid,'(A,3(1X,I8))') '   Sending ...', id, start_band, nbands_to_send
          call MPI_IRSEND( cmplx_wvfn( 1, 1, j ), nbands_to_send*test_gvec, MPI_DOUBLE_COMPLEX, &
                         id, 1, pool_comm, requests( j ), ierr )
          if( ierr .ne. 0 ) return
      
        endif

#if 0
        !TODO
        ! Use MPI_TYPE or VECTOR to only send the first test_gvec 
        ! For now, no worries
        ! don't send if I am me
        if( id .ne. pool_myid ) then
          write(1000+myid,'(A,3(1X,I8))') '   Sending ...', id, start_band, nbands_to_send
          call MPI_IRSEND( cmplx_wvfn( 1, start_band ), nbands_to_send*test_gvec, MPI_DOUBLE_COMPLEX, &
                         id, 1, pool_comm, requests( id ), ierr )
          ! this might not sync up
          if( ierr .ne. 0 ) return
        else
          write(1000+myid,'(A,3(1X,I8))') "   Don't Send: ", start_band, nbands_to_send, my_bands
          wfns( 1:test_gvec, : ) = cmplx_wvfn( 1:test_gvec, start_band : nbands_to_send + start_band - 1 )
        endif
#endif

        start_band = start_band + nbands_to_send
      enddo
#endif

      close( 99 )

      nr=nr+1
    else
      nr = 2
      allocate( requests( nr ), cmplx_wvfn(1,1,1) )
      requests(:) = MPI_REQUEST_NULL
      write(1000+myid,*) '***Receiving k-point: ', ikpt, ispin
      write(1000+myid,*) '   Ngvecs: ', ngvecs
      write(1000+myid,*) '   Nbands: ', nbands

      if( gammaFullStorage .and. is_gamma ) then
        test_gvec = ( ngvecs + 1 ) / 2
      else
        test_gvec = ngvecs 
      endif

      write(1000+myid,*) '   Ngvecs: ', test_gvec
      write(1000+myid,*) '   Gamma : ', is_gamma

#ifdef MPI
      call MPI_TYPE_VECTOR( my_bands, test_gvec, ngvecs, MPI_DOUBLE_COMPLEX, newType, ierr )
      if( ierr .ne. 0 ) return
      call MPI_TYPE_COMMIT( newType, ierr )
      if( ierr .ne. 0 ) return


      call MPI_IRECV( wfns, 1, newType, pool_root, 1, pool_comm, requests( 1 ), ierr )
      if( ierr .ne. 0 ) return
      call MPI_TYPE_FREE( newType, ierr )
      if( ierr .ne. 0 ) return

      call MPI_IBCAST( gvecs, 3*test_gvec, MPI_INTEGER, pool_root, pool_comm, requests( 2 ), ierr )
      call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, ierr_ )
      if( ierr .ne. 0 .or. ierr_ .ne. 0 ) then
        call MPI_CANCEL( requests( 1 ) , ierr )
        call MPI_CANCEL( requests( 2 ) , ierr )
        ierr = 5
        return
      endif

#endif


    endif


    call MPI_WAITALL( nr, requests, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    deallocate( cmplx_wvfn, requests )

    if( gammaFullStorage .and. is_gamma ) then
      j = test_gvec
      do i = 1, test_gvec
        if( gvecs(1,i) .eq. 0 .and. gvecs(2,i) .eq. 0 .and. gvecs(3,i) .eq. 0 ) then
          k = i + 1
          exit
        endif
        j = j + 1
        gvecs(:,j) = -gvecs(:,i)
        wfns(j,:) = conjg( wfns(i,:) )
      enddo

      do i = k, test_gvec
        j = j + 1
        gvecs(:,j) = -gvecs(:,i)
        wfns(j,:) = conjg( wfns(i,:) )
      enddo
    endif

    write(1000+myid,*) '***Finishing k-point: ', ikpt, ispin
    call MPI_BARRIER( pool_comm, ierr )

  end subroutine qe62_read_at_kpt



end module ocean_qe62_files