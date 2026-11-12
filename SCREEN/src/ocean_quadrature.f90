module ocean_quadrature
  use ai_kinds, only : DP 

  implicit none
  private
  save

!  type quadHolder
!    real(WP), allocatable, dimension( : ) :: abscissae
!    real(WP), allocatable, dimension( : ) :: weights
!  end type quadHolder

!  type( quadHolder ), allocatable, protected :: LegendreQuad( : )

!  public :: quadHolder
  public :: ocean_quadrature_loadLegendre

  contains

  subroutine ocean_quadrature_loadLegendre( order, abscissae, weights, ierr )
    integer, intent( in ) :: order
    real(DP), intent( out ) :: abscissae( : )
    real(DP), intent( out ) :: weights( : )
    integer, intent( inout ) :: ierr

    integer :: i, j

    ierr = 0
    if( mod( order, 2 ) .ne. 0 ) then
      ierr = 1
      return
    endif
    if( order .lt. 2 .or. size( abscissae ) .lt. order .or. size( weights ) .lt. order ) then
      ierr = 2
      return
    endif

    open( unit=99, file='EvenQuadHalf.txt', form='formatted', status='old', iostat=ierr )
    if( ierr .ne. 0 ) return
    do
      read(99,*,iostat=ierr) i
      if( ierr .ne. 0 ) return
      if( i .lt. 2 .or. mod( i, 2 ) .ne. 0 ) then
        ierr = 1
        return
      endif
      if( i .eq. order ) then
        do j = 1, i/2
          read(99,*,iostat=ierr) abscissae(j), weights(j)
          if( ierr .ne. 0 ) return
          abscissae(i-j+1) = -abscissae(j)
          weights(i-j+1) = weights(j)
        enddo
        exit
      else
        do j = 1, i/2
          read(99,*,iostat=ierr)
          if( ierr .ne. 0 ) return
        enddo
      endif
    enddo

    close( 99, iostat=ierr )

  end subroutine ocean_quadrature_loadLegendre

end module
