! Copyright (C) 2015, 2018 SCREEN collaboration
!
! This file is part of the SCREEN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module SCREEN_timekeeper
  implicit none
  save
  private
  integer(8) :: max_count, count_rate

!  integer, private, parameter :: ndivs = 6
!  integer, public, parameter :: tk_lr = 1, &
!                                tk_mult = 2, &
!                                tk_e0 = 3, &
!                                tk_psisum = 4, &
!                                tk_inv = 5, &
!                                tk_buffer2min = 6
!
!  character(LEN=10), private, parameter :: tk_label( ndivs ) = (/ 'long range', 'multiplet ', &
!    'energies  ', 'sum vector', 'inversion ', 'buffer2min' /)
!
!  integer(8), private :: total( ndivs )
!  integer(8), private :: prev( ndivs )



  integer :: ndivs
  integer :: maxdivs
  character(LEN=32), allocatable :: tk_label( : )
  integer(8), allocatable :: total( : ), prev( : )

  INTERFACE SCREEN_tk_start
    module procedure SCREEN_tk_start_str
    module procedure SCREEN_tk_start_id
  END INTERFACE

  INTERFACE SCREEN_tk_stop
    module procedure SCREEN_tk_stop_str
    module procedure SCREEN_tk_stop_id
  END INTERFACE SCREEN_tk_stop


  public :: SCREEN_tk_start, SCREEN_tk_stop, SCREEN_tk_init, SCREEN_tk_printtimes

  contains

  subroutine string2id( str, id, canAdd )
    character(LEN=*), intent( in ) :: str
    integer, intent( out ) :: id
    logical, intent( in ) :: canAdd

    integer :: i

    do i = 1, ndivs
      if( str .eq. tk_label( i ) ) then
        id = i 
        return
      endif
    enddo

    if( canAdd ) then

      ndivs = ndivs + 1
      if( ndivs .gt. maxdivs ) then
        ! lazy
        id = -1
        return
      endif

      tk_label( ndivs ) = str
      id = ndivs

    else
      id = -1
    endif
    return

  end subroutine string2id

  subroutine id2string( str, id )
    character(LEN=*), intent( out ) :: str
    integer, intent( in ) :: id

    if( id .gt. ndivs .or. id .lt. 1 ) then
      str = ''
      return
    endif
    str = tk_label( id ) 

  end subroutine id2string

  subroutine SCREEN_tk_init
    implicit none
    integer(8) :: cl
    call SYSTEM_CLOCK( cl, count_rate, max_count )

    maxdivs = 128
    ndivs = 0
    allocate( tk_label( maxdivs ), total( maxdivs ), prev( maxdivs ) )

    total( : ) = 0
  end subroutine SCREEN_tk_init

  subroutine SCREEN_tk_start_str( str )
    implicit none
    character(LEN=*), intent( in ) :: str
    integer :: id

    call string2id( str, id, .true. )
    if( id .lt. 1 ) then
      write(6,*) 'Programming bug, too many time keepers'
      return
    endif
    call SCREEN_tk_start_id( id )
  end subroutine SCREEN_tk_start_str

  subroutine SCREEN_tk_start_id( id )
    implicit none
    integer, intent(in) :: id
    !
    call SYSTEM_CLOCK( prev( id ) )
    !
  end subroutine SCREEN_tk_start_id

  subroutine SCREEN_tk_stop_str( str )
    implicit none
    character(LEN=*), intent( in ) :: str
    integer :: id

    call string2id( str, id, .false. )
    if( id .lt. 1 ) then
      write(6,*) 'Programming bug. String stopped, but not started', str
      return
    endif
    call SCREEN_tk_stop_id( id )
  end subroutine SCREEN_tk_stop_str

  subroutine SCREEN_tk_stop_id( id )
    implicit none
    integer, intent(in) :: id
    integer(8) :: cl
    !
    call SYSTEM_CLOCK( cl )
    if( prev( id ) .gt. cl ) then
      total( id ) = total( id ) + ( cl - ( prev( id ) - max_count ) )
    else
      total( id ) = total( id ) + ( cl - prev( id ) )
    endif
    !
  end subroutine SCREEN_tk_stop_id


  subroutine SCREEN_tk_printtimes( myid, comm )
#ifdef MPI
#ifndef __OLD_MPI

#ifdef MPI_F08
    use mpi_f08
#else
    use mpi
#endif

#endif
#endif
    use ai_kinds, only :DP
    implicit none
#ifdef MPI
#ifdef __OLD_MPI
    include 'mpif.h'
#endif
#endif
    integer, intent( in ) :: myid
#ifdef MPI_F08
    type( MPI_COMM ), optional, intent(in) :: comm
#else
    integer, optional, intent(in) :: comm
#endif

    character(LEN=16) :: filnam
    integer :: iter, nproc, ierr, i, j, nReport, maxReport
    integer(8) :: iavg, imax, imin
    real(DP) :: avg, stddev
    integer, allocatable :: rankCount(:)
    integer(8), allocatable :: allTotal(:,:), rankTotal(:)
    character(LEN=32), allocatable :: allLabel(:,:), reportLabel(:)
#ifdef MPI_F08
    TYPE(MPI_Status) :: stat
#else
    integer :: stat(MPI_STATUS_SIZE)
#endif

!   This can be changed if for some reason we run on all the computers
    if( myid .gt. 999999 ) return 

    ! New behavior
    nproc = 1
    if( present( comm ) ) then
      call MPI_COMM_SIZE( comm, nproc, ierr )
      if( ierr .ne. 0) return
    endif

    if( nproc .gt. 1 ) then

      if( myid .eq. 0 ) then
        allocate( rankCount(0:nproc-1) )
        allocate( allTotal( maxdivs, 0:nproc-1 ), allLabel( maxdivs, 0:nproc-1 ) )
        rankCount(:) = 0
        allTotal(:,:) = 0
        allLabel(:,:) = ''
        rankCount(0) = ndivs
        allTotal(1:ndivs,0) = total(1:ndivs)
        allLabel(1:ndivs,0) = tk_label(1:ndivs)
        do iter = 1, nproc-1
          call MPI_RECV( rankCount(iter), 1, MPI_INTEGER, iter, 1, comm, stat, ierr )
          if( rankCount(iter) .gt. 0 ) then
            call MPI_RECV( allLabel(1:rankCount(iter),iter), 32*rankCount(iter), &
                           MPI_CHARACTER, iter, 2, comm, stat, ierr )
            call MPI_RECV( allTotal(1:rankCount(iter),iter), rankCount(iter), &
                           MPI_INTEGER8, iter, 3, comm, stat, ierr )
          endif
        enddo
      else
        call MPI_SEND( ndivs, 1, MPI_INTEGER, 0, 1, comm, ierr )
        if( ndivs .gt. 0 ) then
          call MPI_SEND( tk_label(1:ndivs), 32*ndivs, MPI_CHARACTER, 0, 2, comm, ierr )
          call MPI_SEND( total(1:ndivs), ndivs, MPI_INTEGER8, 0, 3, comm, ierr )
        endif
      endif

      if( myid .eq. 0 ) then
        maxReport = sum( rankCount(:) )
        allocate( reportLabel( maxReport ), rankTotal( 0:nproc-1 ) )
        nReport = 0
        do i = 1, rankCount(0)
          nReport = nReport + 1
          reportLabel( nReport ) = allLabel(i,0)
        enddo
        do iter = 1, nproc-1
          do i = 1, rankCount(iter)
            do j = 1, nReport
              if( allLabel(i,iter) .eq. reportLabel(j) ) exit
            enddo
            if( j .gt. nReport ) then
              nReport = nReport + 1
              reportLabel( nReport ) = allLabel(i,iter)
            endif
          enddo
        enddo

        write(filnam,'(A6,I6.6,A4)' ) 'timing', myid, '.txt'
        open(unit=99,file=filnam,form='formatted',status='unknown')
        rewind(99)

        ! 32 + 1 + 20 + 1 + 5
        write(99,'(A1,A69,A,A22,A17,A17,A17,A17)') '#', '', 'Root', 'Avg', 'StdDev', 'Min', 'Max'
        do iter = 1, nReport
          ! Missing labels contribute zero ticks for that rank.
          rankTotal(:) = 0
          do i = 0, nproc-1
            do j = 1, rankCount(i)
              if( allLabel(j,i) .eq. reportLabel(iter) ) then
                rankTotal(i) = allTotal(j,i)
                exit
              endif
            enddo
          enddo
          iavg = sum( rankTotal(:) )
          imax = maxval( rankTotal(:) )
          imin = minval( rankTotal(:) )
          avg = dble( iavg ) / dble( nproc )
          stddev = 0.0_DP
          do i = 0, nproc-1
            stddev = stddev + ( avg - dble( rankTotal(i) ) ) ** 2
          enddo
          stddev = sqrt( stddev/dble(nproc) )
          ! 10 digits before the decimal should get us 10 years in seconds
          write( 99, '(A,1X,I20,1X,A5,F16.6,1X,A5,1X,F16.6,1X,F16.6,1X,F16.6,1X,F16.6)' ) &
                  reportLabel( iter ), rankTotal(0), ' tics', &
                  (dble( rankTotal(0) )/dble(count_rate)), 'secs', &
                  avg/dble(count_rate), stddev/dble(count_rate), &
                  dble(imin)/dble(count_rate), dble(imax)/dble(count_rate)
        enddo
        close(99)

        deallocate( rankCount, allTotal, allLabel, reportLabel, rankTotal )
      endif

    else
      write(filnam,'(A6,I6.6,A4)' ) 'timing', myid, '.txt'
      open(unit=100+myid,file=filnam,form='formatted',status='unknown')
      rewind(100+myid)

      do iter = 1, ndivs
        write(100+myid,'(A,1X,I20,1X,A,F24.6,1X,A)') tk_label( iter ), total( iter ), ' tics', &
                  (dble( total( iter ) )/dble(count_rate)), 'secs'
      enddo

      close(100+myid)
    endif

  end subroutine SCREEN_tk_printtimes


end module SCREEN_timekeeper
