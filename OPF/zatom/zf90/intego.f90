! Copyright (C) 2010,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine intego( e, l, kap, n, nn, is, ief, x0, phi, z, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead, der0, der1 )
  implicit none
  !
  integer :: l, n, nn, is, ief, nr
  real( kind = kind( 1.0d0 ) ) :: e, kap, x0, z, dl, rel, plead, der0, der1
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: phi, v, xm1, xm2, r, dr, r2
  !
  integer :: i, nnideal
  real( kind = kind( 1.0d0 ) ) :: p0, p1, p2
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: xm, extf, midf
  logical :: cont, ready
  !
  ! we shall set ief to -1 if energy is too low, +1 if too high.
  ief=0
  !
  ! prepare certain parameters outside loop
  call setxkdk( e, l, kap, rel, dl, nr, r, r2, v, xm1, xm2, xm, extf, midf )
  !
  ! see Desclaux for origin of equations. here, get two points...
  p0 = r( 1 ) ** ( plead - 0.5d0 ) * ( 1.0d0 - r( 1 ) * z / dble( l + 1 ) ) / sqrt( xm( 1 ) )
  p1 = r( 2 ) ** ( plead - 0.5d0 ) * ( 1.0d0 - r( 2 ) * z / dble( l + 1 ) ) / sqrt( xm( 2 ) )
  phi( 1 ) = p0 * sqrt( xm( 1 ) * r( 1 ) )
  phi( 2 ) = p1 * sqrt( xm( 2 ) * r( 2 ) )
  !
  ! initialize number of nodes, and determine the ideal number.
  nn = 0
  nnideal = n - l - 1
  !
  ! integ out.  count nodes, stop along way if there are too many.
  cont = .true.
  ready = .false.
  i = 3 
  do while ( cont .and. ( ief .eq. 0 ) )
     if ( extf( i ) .lt. 0.0d0 ) ief = -1
     p2 = ( midf( i - 1 ) * p1 - extf( i - 2 ) * p0 ) / extf( i )
     phi( i ) = p2 * sqrt( xm( i ) * r( i ) )
     if ( abs( p2 ) .gt. 1.0d10 ) call temper( 1, i, phi( 1 ), p0, p1, p2 )
     if ( p2 * p1 .lt. 0.0d0 ) nn = nn + 1
     if ( ( nn .gt. nnideal ) .and. ( ief .eq. 0 ) ) ief = 1
     if ( ( nn .eq. nnideal ) .and. ( p2 * p1 .gt. 0.0d0 ) .and. ( p2 / p1 .lt. 1.0d0 ) ) ready = .true.
     if ( ( is .eq. 0 ) .and. ( ready ) ) then
        if ( e .lt. v( i ) + dble( l * ( l + 1 ) ) / ( 2.0d0 * r2( i ) ) ) then
           is = i - 2
           cont = .false.
        end if
     end if
     if ( ( is .ne. 0 ) .and. ( i .eq. is + 2 ) ) cont = .false.
     p0 = p1
     p1 = p2
     i = i + 1
     if ( ( i .gt. nr ) .and. ( ief .eq. 0 ) ) ief = -1
  end do
  !
  if ( ief .eq. 0 ) then
     call getxval( phi( is - 2 ), dl, r( is ), x0 )
     der0 = phi( is )
     der1 = ( 8.0d0 * ( phi( is + 1 ) - phi( is - 1 ) ) - ( phi( is + 2 ) - phi( is - 2 ) ) ) / ( 12.0d0 * dl * r( is ) )
  end if
  !
  return
end subroutine intego
