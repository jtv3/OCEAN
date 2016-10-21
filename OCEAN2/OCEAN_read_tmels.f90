subroutine OCEAN_read_tmels( sys, p, file_selector, ierr )
  use OCEAN_mpi
  use OCEAN_system
  use OCEAN_psi

  implicit none

  type(O_system), intent( in ) :: sys
  type(OCEAN_vector), intent( inout ) :: p
  integer, intent( in ) :: file_selector
  integer, intent( inout ) :: ierr


  integer :: nbc(2), nbv, nk, ik, fh, elements, ic, i, j

  real(dp) :: inv_qlength, qinb(3), max_psi, su
  complex(dp), allocatable :: psi_in(:,:)
  real(dp), allocatable :: psi_transpose( :, : ), re_wgt(:,:,:), im_wgt(:,:,:)

#ifdef MPI
  integer(MPI_OFFSET_KIND) :: offset
#endif

  max_psi = 0.0_dp

! qlength
  qinb(:) = sys%qinunitsofbvectors(:)
    
  inv_qlength = (qinb(1) * sys%bvec(1,1) + qinb(2) * sys%bvec(1,2) + qinb(3) * sys%bvec(1,3) ) ** 2 &
              + (qinb(1) * sys%bvec(2,1) + qinb(2) * sys%bvec(2,2) + qinb(3) * sys%bvec(2,3) ) ** 2 &
              + (qinb(1) * sys%bvec(3,1) + qinb(2) * sys%bvec(3,2) + qinb(3) * sys%bvec(3,3) ) ** 2 
  inv_qlength = 1.0_dp / dsqrt( inv_qlength )
!  inv_qlength = dsqrt( inv_qlength )


  if( sys%cur_run%num_bands .ne. ( sys%brange(4)-sys%brange(3)+1 ) ) then
    if(myid .eq. root ) write(6,*) 'Conduction band mismatch!', sys%cur_run%num_bands, ( sys%brange(4)-sys%brange(3)+1 ) 
    ierr = -1
    return
  endif
  if( sys%cur_run%val_bands .ne. ( sys%brange(2)-sys%brange(1)+1 ) ) then
    if(myid .eq. root ) write(6,*) 'Valence band mismatch!', sys%cur_run%val_bands, ( sys%brange(2)-sys%brange(1)+1 )
    ierr = -1
    return
  endif


  select case (file_selector )
  
  case( 1 )

    if( myid .eq. root ) then
      write(6,*) 'Inverse Q-length:', inv_qlength
      open(unit=99,file='tmels.info',form='formatted',status='old')
      read(99,*) nbv, nbc(1), nbc(2), nk
      close(99)
      if( nk .ne. sys%nkpts ) then
        write(6,*) 'tmels.info mismatch: nkpts'
        ierr = -1
        return
      endif
    endif

#ifdef MPI
    call MPI_BCAST( nbc, 2, MPI_INTEGER, 0, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( nbv, 1, MPI_INTEGER, 0, comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_FILE_OPEN( comm, 'ptmels.dat', MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    offset = 0
    call MPI_FILE_SET_VIEW( fh, offset, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, & 
                            'native', MPI_INFO_NULL, ierr)
    if( ierr .ne. MPI_SUCCESS ) return
#else
    open( unit=99,file='ptmels.dat',form='binary',status='old')
#endif
    allocate( psi_in( nbv, nbc(1):nbc(2) ), psi_transpose( sys%cur_run%val_bands, sys%cur_run%num_bands ) )


    do ik = 1, sys%nkpts
      ! Read in the tmels
#ifdef MPI
      offset = nbv * ( nbc(2) - nbc(1) + 1 ) * ( ik - 1 )
      elements = nbv * ( nbc(2) - nbc(1) + 1 )
      call MPI_FILE_READ_AT_ALL( FH, offset, psi_in, elements, &
                                 MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr )
      if( ierr .ne. MPI_SUCCESS) return
#else
      read(99) psi_in
#endif

!      max_psi = max( max_psi, maxval( real(psi_in(:,:) ) ) )

      psi_transpose( :, : ) = inv_qlength * real( psi_in( sys%brange(1):sys%brange(2), sys%brange(3):sys%brange(4) ), DP )
      p%valr(1:sys%cur_run%num_bands,1:sys%cur_run%val_bands,ik,1) = transpose( psi_transpose )

      psi_transpose( :, : ) = inv_qlength * real( aimag( psi_in( sys%brange(1):sys%brange(2), sys%brange(3):sys%brange(4) ) ), DP )
      p%vali(1:sys%cur_run%num_bands,1:sys%cur_run%val_bands,ik,1) = transpose( psi_transpose )


      max_psi = max( max_psi, maxval( p%valr(1:sys%cur_run%num_bands,1:sys%cur_run%val_bands,ik,1) ) )
    enddo

#ifdef MPI
    call MPI_FILE_CLOSE( fh, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
#else
    close(99)
#endif
    deallocate( psi_in, psi_transpose )
  case( 0 )

    if( myid .eq. root ) then
      write(6,*) 'Inverse Q-length:', inv_qlength
      open(unit=99,file='tmels.info',form='formatted',status='old')
      read(99,*) nbv, nbc(1), nbc(2), nk
      close(99)
      if( nk .ne. sys%nkpts ) then
        write(6,*) 'tmels.info mismatch: nkpts'
        ierr = -1
        return
      endif


      allocate( re_wgt( 0:3, nbv, nbc(2)-nbc(1)+1 ), im_wgt( 0:3, nbv, nbc(2)-nbc(1)+1 ) )
      su = 0.0_DP
      fh = 99
      open( unit=fh, file='tmels', status='unknown' )
      rewind fh
      do ik = 1, sys%nkpts
!        do ic = 0, 3
          do j = 1, nbc(2) - nbc(1) + 1
            do i = 1, nbv
!              read ( fh, '(2(1x,1e22.15))' ) re_wgt( 0, i, j ), im_wgt( 0, i, j )
              read ( fh, * ) re_wgt( 0, i, j ), im_wgt( 0, i, j )
            end do
          end do
!        end do
        do i = 1, sys%cur_run%val_bands
          do j = 1, sys%cur_run%num_bands
             p%valr( j, i, ik, 1 ) = re_wgt( 0, i, j )
             p%vali( j, i, ik, 1 ) = im_wgt( 0, i, j )
             su = su + re_wgt( 0, i, j ) * re_wgt( 0, i, j ) + im_wgt( 0, i, j ) * im_wgt( 0, i, j )
          end do
        end do
      end do
      close( unit=fh )
      !
      deallocate( re_wgt, im_wgt )
      su = sqrt( su )
      write(6,*) 'pnorm=', su, sys%cur_run%val_bands, sys%cur_run%num_bands

      su = 0.0_DP
      do ik = 1, sys%nkpts
        do i = 1, sys%cur_run%val_bands
          do j = 1, sys%cur_run%num_bands
            su = su + p%valr( j, i, ik, 1 ) * p%valr( j, i, ik, 1 ) + p%vali( j, i, ik, 1 ) * p%vali( j, i, ik, 1 )
          enddo
        enddo
      enddo
      su = sqrt( su )
      write(6,*) 'pnorm=', su


    endif
    
#ifdef MPI
    ! BCAST psi for parallel
    if( nproc .gt. 1 ) ierr = -1000
#endif

  
  case default
    if( myid .eq. root ) write(6,*) "Unsupported tmels format requested: ", file_selector
    ierr = -1
    return
  end select

  if( myid .eq. root ) write(6,*) 'Max Psi:', max_psi

end subroutine OCEAN_read_tmels