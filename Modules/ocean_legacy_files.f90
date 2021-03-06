! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
module ocean_legacy_files
  use ai_kinds, only : DP
#ifdef MPI_F08
  use mpi_f08, only : MPI_COMM, MPI_REQUEST
#endif


  implicit none
  private
  save


  logical :: is_init = .false.
  integer :: nfiles
  integer, allocatable :: file_indx( : )
  character( len=12 ), allocatable :: file_names( : )

  integer :: bands(2)
  integer :: brange(4)
  integer :: kpts(3)
  integer :: nspin 


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

  public :: olf_read_init, olf_read_at_kpt, olf_clean, olf_read_energies, olf_get_ngvecs_at_kpt, &
            olf_read_energies_single
  public :: olf_read_at_kpt_split
  public :: olf_kpts_and_spins, olf_is_my_kpt
  public :: olf_return_my_bands, olf_return_my_val_bands, olf_return_my_con_bands
  public :: olf_nprocPerPool, olf_getPoolIndex, olf_returnGlobalID
  public :: olf_getAllBandsForPoolID, olf_getValenceBandsForPoolID, olf_getConductionBandsForPoolID
  public :: olf_npool, olf_universal2KptAndSpin, olf_poolID, olf_poolComm

  contains

  subroutine olf_return_my_val_bands( nb, ierr)
    integer, intent( out ) :: nb
    integer, intent( inout ) :: ierr
    if( .not. is_init ) then
      ierr = 1
      return
    endif
    nb = olf_getValenceBandsForPoolID( pool_myid )
  end subroutine olf_return_my_val_bands

  subroutine olf_return_my_con_bands( nb, ierr)
    integer, intent( out ) :: nb
    integer, intent( inout ) :: ierr
    if( .not. is_init ) then
      ierr = 1
      return
    endif
    nb = olf_getConductionBandsForPoolID( pool_myid )
  end subroutine olf_return_my_con_bands

  pure function olf_poolComm()
#ifdef MPI_F08
    type( MPI_COMM ) :: olf_poolComm
#else
    integer :: olf_poolComm
#endif

    olf_poolComm = pool_comm
  end function olf_poolComm

  pure function olf_npool() result( npool_ )
    integer :: npool_

    npool_ = npool
  end function olf_npool

  pure function olf_nprocPerPool() result( nproc )
    integer :: nproc
    
    nproc = pool_nproc
  end function olf_nprocPerPool

  pure function olf_poolID() result( pid )
    integer :: pid

    pid = pool_myid
  end function olf_poolID


  pure function olf_getPoolIndex( ispin, ikpt ) result( poolIndex )
    integer, intent( in ) :: ispin, ikpt
    integer :: poolIndex
    integer :: kptCounter
    !
    kptCounter = ikpt + ( ispin - 1 ) * product(kpts(:))
    poolIndex = mod( kptCounter, npool )
  end function olf_getPoolIndex

  pure function olf_getAllBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

!    bands_remain = brange(4)-brange(3)+brange(2)-brange(1)+2
    bands_remain = bands(2) - bands(1) + 1

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function olf_getAllBandsForPoolID

  pure function olf_getValenceBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

!    bands_remain = brange(4)-brange(3)+brange(2)-brange(1)+2
    bands_remain = brange(2) - brange(1) + 1

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function olf_getValenceBandsForPoolID

  pure function olf_getConductionBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

!    bands_remain = brange(4)-brange(3)+brange(2)-brange(1)+2
    bands_remain = brange(4) - brange(3) + 1

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function olf_getConductionBandsForPoolID

  pure function olf_returnGlobalID( poolIndex, poolID ) result( globalID )
    integer, intent( in ) :: poolIndex, poolID
    integer :: globalID

    globalID = poolIndex * pool_nproc + poolID
  end function olf_returnGlobalID

  subroutine olf_universal2KptAndSpin( uni, ikpt, ispin )
    integer, intent( in ) :: uni
    integer, intent( out ) :: ispin, ikpt
    !
    integer :: i, ierr
    logical :: is_kpt

    i = 0
    do ispin = 1, nspin
      do ikpt = 1, product(kpts(:))
        call olf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
        if( is_kpt ) i = i +1
        if( uni .eq. i ) return
      enddo
    enddo

    ikpt = 0
    ispin = 0
  end subroutine olf_universal2KptAndSpin

  subroutine olf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
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
    if( olf_getPoolIndex( ispin, ikpt ) .eq. mypool ) then
      is_kpt = .true.
    else
      is_kpt = .false.
    endif

  end subroutine olf_is_my_kpt

  subroutine olf_return_my_bands( nbands, ierr )
    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr
    
    if( .not. is_init ) ierr = 1
    nbands = pool_nbands
  end subroutine olf_return_my_bands

  integer function olf_kpts_and_spins()
    olf_kpts_and_spins = product(kpts(:)) * nspin / npool
    return
  end function olf_kpts_and_spins

  subroutine olf_read_energies_single( myid, root, comm, energies, ierr )
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
    integer :: ispn, ikpt, bstop1, bstart2, bstop2

    if( is_init .eqv. .false. ) then
      write(6,*) 'olf_read_energies_single called but was not initialized'
      ierr = 4
      return
    endif 

    if( size( energies, 1 ) .lt. ( bands(2) - bands(1)+1 ) ) then
      write(6,*) 'Error! Band dimension of energies too small in olf_read_energies_single'
      ierr = 1
      return
    endif
    if( size( energies, 2 ) .lt. product(kpts(:)) ) then
      write(6,*) 'Error! K-point dimension of energies too small in olf_read_energies_single'
      ierr = 2
      return
    endif
    if( size( energies, 3 ) .lt. nspin ) then
      write(6,*) 'Error! Spin dimension of energies too small in olf_read_energies_single'
      ierr = 3
      return
    endif


    if( myid .eq. root ) then
      bstop1  = brange(2) - brange(1) + 1
      bstart2 = brange(3) - brange(1) + 1
      bstop2  = brange(4) - brange(1) + 1

      open( unit=99, file='enkfile', form='formatted', status='old' )
      do ispn = 1, nspin
        do ikpt = 1, product(kpts(:))
          read( 99, * ) energies( 1 : bstop1, ikpt, ispn )
          read( 99, * ) energies( bstart2 : bstop2, ikpt, ispn )
        enddo
      enddo
      close( 99 )
      ! We are using Hartree
      energies( :, :, : ) = 0.5_dp * energies( :, :, : )
    endif

#ifdef MPI
    call MPI_BCAST( energies, size(energies), MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

  end subroutine olf_read_energies_single



  subroutine olf_read_energies( myid, root, comm, nbv, nbc, nkpts, nspns, val_energies, con_energies, ierr )
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
      ! We are using Hartree
      val_energies( :, :, : ) = 0.5_dp * val_energies( :, :, : )
      con_energies( :, :, : ) = 0.5_dp * con_energies( :, :, : )
    endif

#ifdef MPI
    call MPI_BCAST( val_energies, nbv*nkpts*nspns, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( con_energies, nbc*nkpts*nspns, MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

  end subroutine olf_read_energies



  subroutine olf_clean( ierr )
    integer, intent( inout ) :: ierr
    !
    nfiles = 0
    if( allocated( file_indx ) ) deallocate( file_indx )
    if( allocated( file_names ) ) deallocate( file_names )

#ifdef MPI
    call MPI_COMM_FREE( pool_comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_FREE( inter_comm, ierr )
    if( ierr .ne. 0 ) return
#endif

  end subroutine olf_clean

  ! Read the universal little files
  ! 
  subroutine olf_read_init( comm, ierr, bands_ )
    use ocean_mpi, only : MPI_INTEGER, MPI_CHARACTER
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    integer, intent( inout ) :: ierr
    integer, intent( in ), optional :: bands_( 2 )
    !
    integer :: i
    !
    ! Set the comms for file handling
    call MPI_COMM_DUP( comm, inter_comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_SIZE( inter_comm, inter_nproc, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_RANK( inter_comm, inter_myid, ierr )
    if( ierr .ne. 0 ) return


    if( inter_myid .eq. inter_root ) then
      open( unit=99, file='masterwfile', form='formatted', status='old' )
      read( 99, * ) nfiles
      close( 99 )
      allocate( file_names( nfiles ), file_indx( nfiles ) )
      open( unit=99, file='listwfile', form='formatted', status='old' )
      do i = 1, nfiles
        read( 99, * ) file_indx( i ), file_names( i )
      enddo
      close( 99 )

      open(unit=99,file='brange.ipt', form='formatted', status='old' )
      read(99,*) brange(:)
      close(99)
      open(unit=99,file='kmesh.ipt', form='formatted', status='old' )
      read(99,*) kpts(:)
      close(99)
      open(unit=99,file='nspin', form='formatted', status='old' )
      read(99,*) nspin
      close(99)

      if( present( bands_ ) ) then
        bands(:) = bands_
        if( bands( 1 ) .lt. brange( 1 ) .or. bands( 2 ) .gt. brange( 4 ) ) then
          ierr = 1
          return
        endif
      else
        bands( 1 ) = brange( 1 )
        bands( 2 ) = brange( 4 )
      endif
    endif

#ifdef MPI
    call MPI_BCAST( nfiles, 1, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    if( inter_myid .ne. inter_root ) then
      allocate( file_names( nfiles ), file_indx( nfiles ) )
    endif

    call MPI_BCAST( file_indx, nfiles,  MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( file_names, 12 * nfiles, MPI_CHARACTER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( brange, 4, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( kpts, 3, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( nspin, 1, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( bands, 2, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return
#endif

    call set_pools( ierr )
    if( ierr .ne. 0 ) return

    call writeDiagnostics( )

!    write(6,*) 'olf_read_init was successful'
    is_init = .true.
  
  end subroutine  olf_read_init 

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

    npool = 1
    mypool = 0
    
    call MPI_COMM_SPLIT( inter_comm, mypool, inter_myid, pool_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_COMM_RANK( pool_comm, pool_myid, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_SIZE( pool_comm, pool_nproc, ierr )
    if( ierr .ne. 0 ) return


    pool_nbands = olf_getAllBandsForPoolID( pool_myid )
!    nbands_left = brange(4)-brange(3)+brange(2)-brange(1)+2
!    do i = 0, pool_nproc-1
!      nbands = nbands_left / ( pool_nproc - i )
!      if( i .eq. pool_myid ) pool_nbands = nbands
!      nbands_left = nbands_left - nbands
!    enddo

  end subroutine set_pools


  subroutine olf_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
    use OCEAN_mpi
    integer, intent( in ) :: ikpt, ispin
    integer, intent( out ) :: gvecs
    integer, intent( inout ) :: ierr
    !
    integer :: i, itarg
    logical :: is_kpt

    call olf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
    if( ierr .ne. 0 ) return
    if( is_kpt .eqv. .false. ) then 
      gvecs = 0
      return
    endif

    if( pool_myid .eq. pool_root ) then
      ! determine which file we want (usually we have 1 kpt per file)

!      itarg = 0
!      if( file_indx( nfiles ) .gt. ikpt ) then
!        do i = min( ikpt, nfiles-1), 1, -1
!          if( file_indx( i ) .le. ikpt .and. file_indx( i + 1 ) .gt. ikpt ) then
!            itarg = i
!            exit
!          endif
!        enddo
!      else
!        itarg = nfiles
!      endif
      itarg = wvfn_file_indx( ikpt, ispin )

      if( itarg .eq. 0 ) then
        write( 6, * ) "Couldn't figure out which file to open", ikpt
        ierr = 12
        goto 111
      endif
      open( unit=99,file=file_names(itarg), form='unformatted', status='old' )
      read( 99 ) gvecs
      close( 99 )
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

  end subroutine olf_get_ngvecs_at_kpt

  subroutine olf_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
#ifdef MPI
    use OCEAN_mpi, only : MPI_INTEGER,MPI_DOUBLE_PRECISION, MPI_STATUSES_IGNORE, myid, MPI_REQUEST_NULL
!    use OCEAN_mpi
#endif
    integer, intent( in ) :: ikpt, ispin, ngvecs, my_bands
    integer, intent( out ) :: gvecs( 3, ngvecs )
    complex( DP ), intent( out ) :: wfns( ngvecs, my_bands )
    integer, intent( inout ) :: ierr
    !
    real( DP ), allocatable, dimension( :, : ) :: re_wvfn, im_wvfn
    integer, allocatable, dimension( :, : ) :: trans_gvecs
    integer :: test_gvec, itarg, nbands_to_send, nr, ierr_, nbands_to_read, id, start_band, stop_band, band_overlap, i
#ifdef MPI_F08
    type( MPI_REQUEST ), allocatable :: requests( : )
#else
    integer, allocatable :: requests( : )
#endif

    if( olf_getPoolIndex( ispin, ikpt ) .ne. mypool ) return
    !
    if( olf_getAllBandsForPoolID( pool_myid ) .ne. my_bands ) then
      ierr = 1
      return
    endif

    nbands_to_read = brange(4)-brange(3)+brange(2)-brange(1)+2

    if( pool_myid .eq. pool_root ) then
      itarg = wvfn_file_indx( ikpt, ispin )
      if( itarg .lt. 1 ) then
        ierr = -1
        goto 10
      endif

      open( unit=99, file=file_names( itarg ), form='unformatted', status='old' )
      read( 99 ) test_gvec
      if( test_gvec .ne. ngvecs ) then
        ierr = -2
      endif

      ! Error synch. Also ensures that the other procs have posted their recvs
10    continue
      call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, ierr_ )
      if( ierr .ne. 0 ) return
      if( ierr_ .ne. 0 ) then
        ierr = ierr_
        return
      endif

      allocate( trans_gvecs( ngvecs, 3 ) )
      read( 99 ) trans_gvecs
      gvecs = transpose( trans_gvecs )
      deallocate( trans_gvecs )

      nr = 2 * ( pool_nproc - 1 ) + 1 
!      nr = 2 * pool_nproc
      allocate( requests( 0:nr ) )
      requests(:) = MPI_REQUEST_NULL
#ifdef MPI
      call MPI_IBCAST( gvecs, 3*ngvecs, MPI_INTEGER, pool_root, pool_comm, requests( nr ), ierr )
#endif
      allocate( re_wvfn( ngvecs, nbands_to_read ), im_wvfn( ngvecs, nbands_to_read ) )

      write(1000+myid,*) '***Reading k-point: ', ikpt, ispin
      write(1000+myid,*) '   Ngvecs: ', ngvecs
      write(1000+myid,*) '   Nbands: ', nbands_to_read


      read( 99 ) re_wvfn

      ! slow way, but easy to program
      stop_band = 1
      if( brange( 2 ) .ge. brange( 3 ) ) then
        start_band = brange(2) - brange(1) + 1
        band_overlap = brange(2) - brange(3) + 1
        stop_band = brange(1) + band_overlap
        do i = start_band, stop_band, - 1
          write(1000+myid,*) 'Sending', i-band_overlap, 'to' , i
          re_wvfn( :, i ) = re_wvfn( :, i - band_overlap )
        enddo
      endif

      start_band = stop_band
#ifdef MPI
      ! loop over each proc in this pool to send wavefunctions
      do id = 0, pool_nproc - 1
        nbands_to_send = olf_getAllBandsForPoolID( id )

        ! don't send if I am me
        if( id .ne. pool_myid ) then
          write(1000+myid,'(A,3(1X,I8))') '   Sending ...', id, start_band, nbands_to_send
          call MPI_IRSEND( re_wvfn( 1, start_band ), nbands_to_send*ngvecs, MPI_DOUBLE_PRECISION, &
                         id, 1, pool_comm, requests( id ), ierr )
          ! this might not sync up
          if( ierr .ne. 0 ) return
        else
          write(1000+myid,'(A,3(1X,I8))') "   Don't Send: ", start_band, nbands_to_send, my_bands
        endif

        if( mod( id+1, 16 ) .eq. 0 ) then
          call MPI_WAITALL( 16, requests( id-15 ), MPI_STATUSES_IGNORE, ierr )
          ! this might not sync up
          if( ierr .ne. 0 ) return
          write(1000+myid, '(A,1X,I8)' ) '   triggered intermediate send', id
        endif


        start_band = start_band + nbands_to_send
      enddo
#endif

      read( 99 ) im_wvfn
      stop_band = 1
      if( brange( 2 ) .ge. brange( 3 ) ) then
        start_band = brange(2) - brange(1) + 1
        band_overlap = brange(2) - brange(3) + 1
        stop_band = brange(1) + band_overlap
        do i = start_band, stop_band, - 1
          im_wvfn( :, i ) = im_wvfn( :, i - band_overlap )
        enddo
      endif
      start_band = stop_band

#ifdef MPI
     ! loop over each proc in this pool to send imag wavefunctions
      do id = 0, pool_nproc - 1
        nbands_to_send = olf_getAllBandsForPoolID( id )

        ! don't send if I am me
        if( id .ne. pool_myid ) then
          call MPI_IRSEND( im_wvfn( 1, start_band ), nbands_to_send*ngvecs, MPI_DOUBLE_PRECISION, &
                         id, 2, pool_comm, requests( id + pool_nproc - 1), ierr )
          ! this might not sync up
          if( ierr .ne. 0 ) return
        endif

        if( mod( id+1, 16 ) .eq. 0 ) then
          call MPI_WAITALL( 16, requests( id + pool_nproc - 16 ), MPI_STATUSES_IGNORE, ierr )
          ! this might not sync up
          if( ierr .ne. 0 ) return
          write(1000+myid, '(A,1X,I8)' ) '   triggered intermediate send', id
        endif

        start_band = start_band + nbands_to_send
      enddo
#endif

      close( 99 )

    else
      nr = 3
      allocate( requests( nr ) )
      requests(:) = MPI_REQUEST_NULL
      write(1000+myid,*) '***Receiving k-point: ', ikpt, ispin
      write(1000+myid,*) '   Ngvecs: ', ngvecs
      write(1000+myid,*) '   Nbands: ', bands(2)-bands(1)+1


      allocate( re_wvfn( ngvecs, my_bands ), im_wvfn( ngvecs, my_bands ) )
#ifdef MPI
      call MPI_IRECV( re_wvfn, ngvecs*my_bands, MPI_DOUBLE_PRECISION, pool_root, 1, pool_comm, & 
                      requests( 1 ), ierr )
      call MPI_IRECV( im_wvfn, ngvecs*my_bands, MPI_DOUBLE_PRECISION, pool_root, 2, pool_comm, &
                      requests( 2 ), ierr )

      call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, ierr_ )
      if( ierr .ne. 0 .or. ierr_ .ne. 0 ) then
        call MPI_CANCEL( requests( 1 ) , ierr )
        call MPI_CANCEL( requests( 2 ) , ierr )
        ierr = 5
        return
      endif

      call MPI_IBCAST( gvecs, 3*ngvecs, MPI_INTEGER, pool_root, pool_comm, requests( 3 ), ierr )
#endif


    endif


    call MPI_WAITALL( nr, requests, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    if( pool_myid .eq. pool_root ) then
      wfns( :, : ) = cmplx( re_wvfn( :, stop_band:my_bands+stop_band-1 ), &
                            im_wvfn( :, stop_band:my_bands+stop_band-1 ), DP )
    else
      wfns( :, : ) = cmplx( re_wvfn( :, 1:my_bands ), im_wvfn( :, 1:my_bands ), DP )
    endif
    deallocate( re_wvfn, im_wvfn, requests )

    write(1000+myid,*) '***Finishing k-point: ', ikpt, ispin
    call MPI_BARRIER( pool_comm, ierr )

  end subroutine olf_read_at_kpt



  subroutine olf_read_at_kpt_split( ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands, &
                                    valGvecs, conGvecs, valUofG, conUofG, ierr )

#ifdef MPI
    use OCEAN_mpi, only : MPI_INTEGER,MPI_DOUBLE_PRECISION, MPI_STATUSES_IGNORE, myid, & 
                          MPI_REQUEST_NULL, MPI_STATUS_IGNORE
#endif
    integer, intent( in ) :: ikpt, ispin, valNGvecs, conNGvecs, valBAnds, conBands
    integer, intent( out ) :: valGvecs( 3, valNgvecs ), conGvecs( 3, conNgvecs )
    complex( DP ), intent( out ) :: valUofG( valNGvecs, valBands ), conUofG( conNGvecs, conBands )
    integer, intent( inout ) :: ierr
    !
    real( DP ), allocatable, dimension( :, : ) :: re_wvfn, im_wvfn
    integer, allocatable, dimension( :, : ) :: trans_gvecs
    integer :: test_gvec, itarg, nbands_to_send, nr, ierr_, nbands_to_read, id, start_band, &
               stop_band, i
#ifdef MPI_F08
    type( MPI_REQUEST ), allocatable :: requests( :, : )
    type( MPI_REQUEST ) :: gvecReq
#else
    integer, allocatable :: requests( :, : )
    integer :: gvecReq
#endif

    ! Test to make sure I'm in the right place
    if( olf_getPoolIndex( ispin, ikpt ) .ne. mypool ) return
    !
    ! Check input consistency
    if( olf_getValenceBandsForPoolID( pool_myid ) .ne. valBands ) then
      write(1000+myid,*) 'Valence-band mismatch', olf_getValenceBandsForPoolID( pool_myid ), valBands
      ierr = 1
      return
    endif
    if( olf_getConductionBandsForPoolID( pool_myid ) .ne. conBands ) then
      ierr = 2
      return
    endif
    if( valNGvecs .ne. conNGvecs ) then
      ierr = 3
      return
    endif

    ! includes any overlap
    nbands_to_read = brange(4)-brange(3)+brange(2)-brange(1)+2

    if( pool_myid .eq. pool_root ) then
      itarg = wvfn_file_indx( ikpt, ispin )
      if( itarg .lt. 1 ) then
        ierr = -1
        goto 10
      endif

      open( unit=99, file=file_names( itarg ), form='unformatted', status='old' )
      read( 99 ) test_gvec
      if( test_gvec .ne. valNGvecs ) then
        ierr = -2
      endif

      ! Error synch. Also ensures that the other procs have posted their recvs
10    continue
      call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, ierr_ )
      if( ierr .ne. 0 ) return
      if( ierr_ .ne. 0 ) then
        ierr = ierr_
        return
      endif

      allocate( trans_gvecs( valNGvecs, 3 ) )
      read( 99 ) trans_gvecs
      valgvecs = transpose( trans_gvecs )
      deallocate( trans_gvecs )

!      nr = 2 * ( pool_nproc - 1 ) + 1
      nr = 2 * pool_nproc 
      allocate( requests( 0:nr-1, 2 ) )
      requests(:,:) = MPI_REQUEST_NULL
#ifdef MPI
      call MPI_IBCAST( valGvecs, 3*valNGvecs, MPI_INTEGER, pool_root, pool_comm, gvecReq, ierr )
#endif
      allocate( re_wvfn( valNGvecs, nbands_to_read ), im_wvfn( valNGvecs, nbands_to_read ) )

      write(1000+myid,*) '***Reading k-point: ', ikpt, ispin
      write(1000+myid,*) '   Ngvecs: ', valNGvecs
      write(1000+myid,*) '   Nbands: ', nbands_to_read


      read( 99 ) re_wvfn

      ! don't delete any overlap for split
!      ! slow way, but easy to program
!      stop_band = 1
!      if( brange( 2 ) .ge. brange( 3 ) ) then
!        start_band = brange(2) - brange(1) + 1
!        band_overlap = brange(2) - brange(3) + 1
!        stop_band = brange(1) + band_overlap
!        do i = start_band, stop_band, - 1
!          write(1000+myid,*) 'Sending', i-band_overlap, 'to' , i
!          re_wvfn( :, i ) = re_wvfn( :, i - band_overlap )
!        enddo
!      endif

      start_band = 1
#ifdef MPI
      ! loop over each proc in this pool to send wavefunctions

      ! first do the valence bands
      do i = 1, 2
        do id = 0, pool_nproc - 1
  !        nbands_to_send = olf_getAllBandsForPoolID( id )
          if( i .eq. 1 ) then
            nbands_to_send = olf_getValenceBandsForPoolID( id )
          else
            nbands_to_send = olf_getConductionBandsForPoolID( id )
          endif

          ! don't send if I am me
          if( id .ne. pool_myid ) then
            write(1000+myid,'(A,3(1X,I8))') '   Sending ...', id, start_band, nbands_to_send
            call MPI_IRSEND( re_wvfn( 1, start_band ), nbands_to_send*valngvecs, MPI_DOUBLE_PRECISION, &
                           id, i, pool_comm, requests( id, i ), ierr )
            ! this might not sync up
            if( ierr .ne. 0 ) return
          else
            if( i .eq. 1 ) then
              write(1000+myid,'(A,3(1X,I8))') "   Don't Send: ", start_band, nbands_to_send, valBands
            else
              write(1000+myid,'(A,3(1X,I8))') "   Don't Send: ", start_band, nbands_to_send, conBands
            endif
          endif

          if( mod( id+1, 16 ) .eq. 0 ) then
            call MPI_WAITALL( 16, requests( id-15, i ), MPI_STATUSES_IGNORE, ierr )
            ! this might not sync up
            if( ierr .ne. 0 ) return
            write(1000+myid, '(A,1X,I8)' ) '   triggered intermediate send', id
          endif


          start_band = start_band + nbands_to_send
        enddo
      enddo
#endif

      read( 99 ) im_wvfn
!      stop_band = 1
!      if( brange( 2 ) .ge. brange( 3 ) ) then
!        start_band = brange(2) - brange(1) + 1
!        band_overlap = brange(2) - brange(3) + 1
!        stop_band = brange(1) + band_overlap
!        do i = start_band, stop_band, - 1
!          im_wvfn( :, i ) = im_wvfn( :, i - band_overlap )
!        enddo
!      endif
!      start_band = stop_band

#ifdef MPI
     ! loop over each proc in this pool to send imag wavefunctions
      start_band = 1
      do i = 1, 2
        do id = 0, pool_nproc - 1
!          nbands_to_send = olf_getAllBandsForPoolID( id )
          if( i .eq. 1 ) then
            nbands_to_send = olf_getValenceBandsForPoolID( id )
          else
            nbands_to_send = olf_getConductionBandsForPoolID( id )
          endif

          ! don't send if I am me
          if( id .ne. pool_myid ) then
            call MPI_IRSEND( im_wvfn( 1, start_band ), nbands_to_send*valngvecs, MPI_DOUBLE_PRECISION, &
                           id, 2+i, pool_comm, requests( id + pool_nproc - 1, i), ierr )
            ! this might not sync up
            if( ierr .ne. 0 ) return
          endif

          if( mod( id+1, 16 ) .eq. 0 ) then
            call MPI_WAITALL( 16, requests( id + pool_nproc - 16, i ), MPI_STATUSES_IGNORE, ierr )
            ! this might not sync up
            if( ierr .ne. 0 ) return
            write(1000+myid, '(A,1X,I8)' ) '   triggered intermediate send', id
          endif

          start_band = start_band + nbands_to_send
        enddo
      enddo
#endif

      close( 99 )

      else
      nr = 2
      allocate( requests( nr,2 ) )
      requests(:,:) = MPI_REQUEST_NULL
      write(1000+myid,*) '***Receiving k-point: ', ikpt, ispin
      write(1000+myid,*) '   Ngvecs: ', valNGvecs
      write(1000+myid,*) '   Nbands: ', valBands, conBands


      allocate( re_wvfn( valngvecs, valBands + conBands ), im_wvfn( valngvecs, valBands + conBands ) )
#ifdef MPI
      call MPI_IRECV( re_wvfn, valngvecs*valBands, MPI_DOUBLE_PRECISION, pool_root, 1, pool_comm, &
                      requests( 1, 1 ), ierr )
      call MPI_IRECV( re_wvfn(1,valBands+1), valngvecs*conBands, MPI_DOUBLE_PRECISION, pool_root, 2, pool_comm, &
                      requests( 1, 2 ), ierr )
      call MPI_IRECV( im_wvfn, valngvecs*valBands, MPI_DOUBLE_PRECISION, pool_root, 3, pool_comm, &
                      requests( 2, 1 ), ierr )
      call MPI_IRECV( im_wvfn(1,valBands+1), valngvecs*conBands, MPI_DOUBLE_PRECISION, pool_root, 4, pool_comm, &
                      requests( 2, 2 ), ierr )

      call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, ierr_ )
      if( ierr .ne. 0 .or. ierr_ .ne. 0 ) then
        call MPI_CANCEL( requests( 1, 1 ) , ierr )
        call MPI_CANCEL( requests( 1, 2 ) , ierr )
        call MPI_CANCEL( requests( 2, 1 ) , ierr )
        call MPI_CANCEL( requests( 2, 2 ) , ierr )
        ierr = 5
        return
      endif

      call MPI_IBCAST( valgvecs, 3*valngvecs, MPI_INTEGER, pool_root, pool_comm, gvecReq, ierr )
#endif


    endif


    call MPI_WAITALL( nr*2, requests, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    ! legacy version has uniform G-vector lists, copy val to con
    call MPI_WAIT( gvecReq, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    conGvecs(:,:) = valGvecs(:,:)

    if( pool_myid .eq. pool_root ) then
!      wfns( :, : ) = cmplx( re_wvfn( :, stop_band:my_bands+stop_band-1 ), &
!                            im_wvfn( :, stop_band:my_bands+stop_band-1 ), DP )
      valUofG( :, : ) = cmplx( re_wvfn( :, 1 : valBands ), im_wvfn( :, 1 : valBands ), DP )
      start_band = brange(2) - brange(1) + 2
      stop_band = start_band + conBands - 1
      conUofG( :, : ) = cmplx( re_wvfn( :, start_band:stop_band ), &
                               im_wvfn( :, start_band:stop_band ), DP )
    else
!      wfns( :, : ) = cmplx( re_wvfn( :, 1:my_bands ), im_wvfn( :, 1:my_bands ), DP )
      valUofG( :, : ) = cmplx( re_wvfn( :, 1 : valBands ), im_wvfn( :, 1 : valBands ), DP )
      start_band = valBands + 1
      stop_band = start_band + conBands - 1
      conUofG( :, : ) = cmplx( re_wvfn( :, start_band:stop_band ), &
                               im_wvfn( :, start_band:stop_band ), DP )
    endif
    deallocate( re_wvfn, im_wvfn, requests )

    write(1000+myid,*) '***Finishing k-point: ', ikpt, ispin
    call MPI_BARRIER( pool_comm, ierr )

  end subroutine olf_read_at_kpt_split


  integer function wvfn_file_indx( ikpt_, ispin_ )
    integer, intent( in ) :: ikpt_, ispin_
    !
    integer :: i, ikpt
    !
    wvfn_file_indx = 0

    ikpt = ikpt_ + ( ispin_ - 1 ) * product( kpts(:) )
    wvfn_file_indx = 0
    if( file_indx( nfiles ) .gt. ikpt ) then
      do i = min( ikpt, nfiles-1), 1, -1
        if( file_indx( i ) .le. ikpt .and. file_indx( i + 1 ) .gt. ikpt ) then
          wvfn_file_indx = i
          return
        endif
      enddo
    else
      wvfn_file_indx = nfiles
    endif

    return
  end function

end module ocean_legacy_files
