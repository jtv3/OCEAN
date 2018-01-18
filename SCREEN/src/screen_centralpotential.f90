module screen_centralPotential
  use ai_kinds, only : DP

  implicit none
  private
  save

  type potential
    real(DP), allocatable :: pot(:)
    real(DP), allocatable :: rad(:)
    integer :: z, n, l
!    real(DP) :: deltaR  ! this is the spacing on the written out induced potential
  end type potential



  type( potential ), allocatable, target :: znlPotentials( : )

  public :: potential
  public :: screen_centralPotential_newScreenShell, screen_centralPotential_loadAll
  public :: screen_centralPotential_prepAll, screen_centralPotential_freePot
  public :: screen_centralPotential_loadInternal, screen_centralPotential_freeInternal
  public :: screen_centralPotential_findNextByZ, screen_centralPotential_countByZ
!  public :: screen_centralPotential_NOutGrid, screen_centralPotential_MakeOutGrid

  contains

#if 0
  pure function screen_centralPotential_NOutGrid( Pot, rmax ) result( n )
    type( potential ), intent( in ) :: Pot
    real( DP ), intent( in ) :: rmax
    integer :: n

    n = max( floor( rmax / Pot%deltaR ) + 1, 1 )
  end function screen_centralPotential_NOutGrid

  subroutine screen_centralPotential_MakeOutGrid( Pot, rmax, rlist, ierr )
    type( potential ), intent( in ) :: Pot
    real( DP ), intent( in ) :: rmax
    real( DP ), intent( out ) :: rlist( : )
    integer, intent( inout ) :: ierr

    integer :: i, n  

    n = screen_centralPotential_NOutGrid( Pot, rmax )

    if( size( rlist, 1 ) .lt. n ) then
      ierr = 1
      return
    endif

    do i = 1, n
      rlist( i ) = real( i - 1, DP ) * Pot%deltaR
    enddo

  end subroutine screen_centralPotential_MakeOutGrid
#endif

  pure function screen_centralPotential_countByZ( Z ) result( n )
    integer, intent( in ) :: Z
    integer :: n
    !
    integer :: PotSize, i
    !
    n = 0
    PotSize = size( znlPotentials )
    do i = 1, PotSize
      if( znlPotentials( i )%z .eq. Z ) n = n + 1
    enddo
  end function screen_centralPotential_countByZ

  subroutine screen_centralPotential_findNextByZ( Z, i, PotPointers, ierr )
    integer, intent( in ) :: Z
    integer, intent( inout ) :: i
    type( potential ), pointer, intent( out ) :: PotPointers
    integer, intent( inout ) :: ierr
    !
    integer :: PotSize, j

!    if( size( PotPointers ) .lt. screen_centralPotential_countByZ( Z ) ) then
!      write(6,*) 'PotPointers too small!'
!      ierr = 1
!      return
!    endif

    PotSize = size( znlPotentials )
    
    do j = i+1, PotSize
      if( znlPotentials( j )%z .eq. Z ) then
        PotPointers => znlPotentials( j )
        return
      endif
    enddo

    ierr = 1

  end subroutine screen_centralPotential_findNextByZ
    

  subroutine screen_centralPotential_loadInternal( ierr )
    use ocean_mpi, only : myid, root
    integer, intent( inout ) :: ierr
    !
    integer, allocatable :: tmp_znl( :, : )
    integer :: znlLength, i
    !
    ! load up potentials to be screened later
    ! 2000 should be more than enough for the entire periodic table at the moment
    if( myid .eq. root ) then
      allocate( tmp_znl( 3, 2000 ), stat=ierr )
    else
      allocate( tmp_znl( 1, 1 ), stat=ierr )
    endif
    call screen_centralPotential_prepAll( tmp_znl, znlLength, ierr )
    if( ierr .ne. 0 ) return
    if( myid .eq. root ) then
      do i = 1, znlLength
        write(6,*) tmp_znl( :, i )
      enddo
    endif
    allocate( znlPotentials( znlLength ), stat=ierr )
    if( ierr .ne. 0 ) return
    call screen_centralPotential_loadAll( znlPotentials, tmp_znl, ierr )
    if( ierr .ne. 0 ) return
    deallocate( tmp_znl )

  end subroutine screen_centralPotential_loadInternal

  subroutine screen_centralPotential_freeInternal()
    integer :: i

    do i = 1, size( znlPotentials ) 
      call screen_centralPotential_freePot( znlPotentials( i ) )
    enddo

    if( allocated( znlPotentials ) ) deallocate( znlPotentials )
  end subroutine screen_centralPotential_freeInternal

  subroutine screen_centralPotential_freePot( pot )
    type( potential ), intent( inout ) :: pot


    if( allocated( pot%pot ) ) deallocate( pot%pot )
    if( allocated( pot%rad ) ) deallocate( pot%rad )

  end subroutine screen_centralPotential_freePot

  subroutine screen_centralPotential_loadAll( znlPot, znl, ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER
    type( potential ), intent( out ) :: znlPot( : )
    integer, intent( in ) :: znl(:,:)
    integer, intent( inout ) :: ierr

    integer :: i, Nznl, ierr_

    Nznl = size( znlPot, 1 )
    ierr_ = 0

    do i = 1, Nznl
      if( myid .eq. root ) then
        call screen_centralPotential_load( znl(1,i), znl(2,i), znl(3,i), znlPot(i), ierr )
      endif
#ifdef MPI
      call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
#endif
      if( ierr .ne. 0 ) return
      if( ierr_ .ne. 0 ) then
        ierr = ierr_
        return
      endif
      call screen_centralPotential_share( znlPot(i), myid, root, comm, ierr )
      if( ierr .ne. 0 ) return
    enddo
  end subroutine screen_centralPotential_loadAll

  subroutine screen_centralPotential_share( pot, myid_, root_, comm_, ierr )
    use ocean_mpi, only : MPI_INTEGER, MPI_DOUBLE_PRECISION
    type( potential ), intent( inout ) :: pot
    integer, intent( in ) :: myid_, root_, comm_
    integer, intent( inout ) :: ierr
    !
    integer :: intArray( 4 )

    if( myid_ .eq. root_ ) then
      intArray(4) = size( pot%pot, 1 )
      intArray(1) = pot%z
      intArray(2) = pot%n
      intArray(3) = pot%l
    endif

    call MPI_BCAST( intArray, 4, MPI_INTEGER, root_, comm_, ierr )
    if( ierr .ne. 0 ) return
    if( myid_ .ne. root_ ) then
      allocate( pot%pot( intArray(4) ), pot%rad( intArray(4) ) )
    endif
    call MPI_BCAST( pot%pot, intArray(4), MPI_DOUBLE_PRECISION, root_, comm_, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( pot%rad, intArray(4), MPI_DOUBLE_PRECISION, root_, comm_, ierr )
    if( ierr .ne. 0 ) return

  end subroutine screen_centralPotential_share

      

  subroutine screen_centralPotential_prepAll( znl, n, ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER
    integer, intent( out ) :: znl(:,:)
    integer, intent( out ) :: n
    integer, intent( inout ) :: ierr

    integer :: i, maxSize, j

    if( myid .ne. root ) then
#ifdef MPI
      call MPI_BCAST( n, 1, MPI_INTEGER, root, comm, ierr )
#endif
      return
    endif

    maxSize = size( znl, 2 )
    open( unit=99, file='edgelist', form='formatted', status='old' )

    do i = 1, maxSize
      read(99,*,iostat=j) znl(1:3, i )
      if( j .lt. 0 ) then
        n = i - 1
        close( 99 )
#ifdef MPI
        call MPI_BCAST( n, 1, MPI_INTEGER, root, comm, ierr )
#endif
        return
      elseif( j .gt. 0 ) then
        ierr = j
        return
      endif
    enddo
    close( 99 )

    ! if this isn't here the other mpi will hang above ...
#ifdef MPI
        call MPI_BCAST( n, 1, MPI_INTEGER, root, comm, ierr )
#endif

    write(6,'(A,I8)') 'Too many unique edges! Greater than: ', maxSize
    ierr = 1
  end subroutine 

  subroutine screen_centralPotential_newScreenShell( pot, newPot, rad, ierr )
    type( potential ), intent( in ) :: pot
    type( potential ), intent( out ) :: newPot
    real(DP), intent( in ) :: rad
    integer, intent( inout ) :: ierr

    integer :: i
    real(DP) :: invRad
    character(len=40) :: filename

    if( ( .not. allocated( pot%pot ) ) .or. ( .not. allocated( pot%rad ) ) ) then
      ierr = 5
      return
    endif

    allocate( newPot%pot( size( pot%pot ) ), newPot%rad( size( pot%rad ) ) )
    newPot%rad(:) = pot%rad(:)

    
    invRad = 1.0_DP / rad
    do i = 1, size( pot%pot ) 
      if( newPot%rad( i ) .lt. rad ) then
        newPot%pot( i ) = pot%pot( i ) + invRad
      else
        newPot%pot( i : size( newPot%pot ) ) = 0.0_DP
        exit
      endif
    enddo

    write(6,*) 'newScreenShell', rad, newPot%rad( size( newPot%rad ) )
    write(filename,'(A,F4.2)') 'vpert.', rad
    open(unit=99,file=filename, form='formatted' )
    do i = 1, size( pot%pot ) 
      write(99,* ) newPot%rad( i ), newPot%pot( i ), pot%pot( i )
    enddo
    close( 99 )

    newPot%z = pot%z
    newPot%n = pot%n
    newPot%l = pot%l
!    newPot%deltaR = pot%deltaR

  end subroutine screen_centralPotential_newScreenShell

  subroutine screen_centralPotential_load( z, n, l, pot, ierr )
    integer, intent( in ) :: z, n, l
    type( potential ), intent( out ) :: pot
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: tmpPot(:), tmpRad(:), tmp(:)
    integer :: maxLength, curLength, fh

    character(len=26) :: fileName

    pot%z = z
    pot%n = n
    pot%l = l
!    pot%deltaR = 0.02_DP
    
    write(fileName,'(A17,I3.3,A1,I2.2,A1,I2.2)') 'zpawinfo/vc_barez', z, 'n', n, 'l', l
    fh = 99
    open( unit=fh, file=fileName, form='formatted', status='old' )
    rewind( 99 )

    maxLength = 4000
    curLength = 1

    allocate( tmpPot( maxLength ), tmpRad( maxLength ), STAT=ierr )
    if( ierr .ne. 0 ) return

    do while( maxLength .le. 512000 )

      call doRead( fh, curLength, MaxLength, tmpRad, tmpPot, ierr )

      if( curLength .lt. 0 ) goto 11

      allocate( tmp( maxLength ), stat=ierr )
      if( ierr .ne. 0 ) return
      tmp(:) = tmpPot(:)
      deallocate( tmpPot )
      allocate( tmpPot( maxLength*2 ), stat=ierr )
      if( ierr .ne. 0 ) return
      tmpPot(1:maxLength) = tmp(:)
      tmp(:) = tmpRad(:)
      deallocate( tmpRad )
      allocate( tmpRad( maxLength*2 ), stat=ierr )
      tmpRad(1:maxLength) = tmp(:)

      curLength = maxLength + 1
      maxLength = 2*maxLength
      deallocate( tmp )
    enddo
    !    If we exit the loop then we've given up not reached the end of the file
    ierr = -1
    return

11  continue
    curLength = abs( curLength )
    allocate( pot%rad( curLength ), pot%pot( curLength ), stat=ierr )
    if( ierr .ne. 0 ) return
    pot%rad(:) = tmpRad(1:curLength)
    pot%pot(:) = tmpPot(1:curLength)
    deallocate( tmpPot, tmpRad )

    close( fh )

    write(6,* ) z, n, l
    write(6,*) curLength
    write(6,*) pot%rad(1), pot%pot(1)
    write(6,*) pot%rad(curLength), pot%pot(curLength)

  end subroutine screen_centralPotential_load

  subroutine doRead( fh, curLength, maxLength, rad, pot, ierr )
    integer, intent( in ) :: fh, maxLength
    integer, intent( inout ) :: curLength
    real(DP), intent( inout ) :: rad(:), pot(:)
    integer, intent( inout ) :: ierr
    !
    integer :: i, j, start

    start = curLength

    do i = start, maxLength
      read(fh,*,iostat=j) rad(i), pot(i)
      if( j .lt. 0 ) then
        curLength = -( i - 1 )
        return
      elseif( j .gt. 0 ) then
        ierr = j
        return
      endif
      ! if j .eq. 0 then keep going!
    enddo

  end subroutine doRead

end module screen_centralPotential