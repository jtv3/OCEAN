! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine confmel( nsphpt, xsph, ysph, zsph, wsph, prefs )
  implicit none
  !
  integer :: nsphpt
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 )
  real( kind = kind( 1.0d0 ) ), dimension( nsphpt ) :: xsph, ysph, zsph, wsph
  !
  integer :: l, m, iproj, mc, npmax, nqproj, atno, nc, lc, lmin, lmax, ip, powmax
  real( kind = kind( 1.0d0 ) ) :: dqproj, qhat( 3 ), ehat( 3 ), ephotev, lam, pi, q, dummy
! real( kind = kind( 1.0d0 ) ) :: qtr, mul, div
  integer, allocatable :: nproj( : )
  character * 4 :: add04
  character * 15 :: spcttype
  character * 17 :: f17
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, : ) :: ifcn
  complex( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, :, : ) :: nbsemel
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  ! get core level info
  open( unit=99, file='ZNL', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) atno, nc, lc
  close( unit=99 )
  write ( add04, '(1a1,1i3.3)' ) 'z', atno
  !
  ! input information about the projectors 
  call nbseprjprep( lmin, lmax, npmax, dqproj, nqproj, add04 )
  allocate( nproj( lmin : lmax ) )
  call nbseprjnproj( lmin, lmax, nproj, add04 )
  allocate( nbsemel( npmax, -lmax : lmax, lmin : lmax, -lc : lc ) )
  !
  ! input information about radial part of matrix elements
  write( f17, '(1a7,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'melfile', 'z', atno, 'n', nc, 'l', lc
  open( unit=99, file=f17, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) powmax   
  allocate( ifcn( npmax, 0 : powmax, lmin : lmax ) ) 
  ifcn( :, :, : ) = 0.0d0 
  do l = lmin, lmax
     do ip = 0, powmax
        read ( 99, * ) ifcn( 1 : nproj( l ), ip, l )
     end do
  end do 
  close( unit=99 )
  !
  ! input information about the probe
  open( unit=99, file='spectfile', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) spcttype
  call fancyvector( ehat, dummy, 99 )
  call fancyvector( qhat, q, 99 )
! read ( 99, * ) ehat( : )
! read ( 99, * ) qhat( : )
  read ( 99, * ) ephotev
! read ( 99, * ) qtr, mul, div
  close( unit=99 )
  ehat( : ) = ehat( : ) / sqrt( sum( ehat( : ) ** 2 ) )
  qhat( : ) = qhat( : ) / sqrt( sum( qhat( : ) ** 2 ) )
  write ( 6, '(3(1f10.5))' ) ehat( : ), qhat( : )
  select case( spcttype )
  case( 'transverse' )
     lam = 10000.0d0 / ( ephotev / 0.806554445d0 )
     q = 0.529177d0 * ( 2.0d0 * pi / lam )
  case( 'longitudinal' )
 !   q = qtr * mul / div
  end select
  write ( 6, '(1a4,1f10.5)' ) 'q = ', q
  !
  ! calculate full matrix elements
  call elssub( lmin, lmax, nproj, npmax, lc, nbsemel, powmax, ifcn, spcttype, ehat, qhat, q, &
       nsphpt, xsph, ysph, zsph, wsph, prefs )
  !
  ! output matrix elements in desired order
  open( unit=99, file='mels', form='formatted', status='unknown' )
  rewind 99
  do mc = -lc, lc
     do l = lmin, lmax
        do m = -l, l
           do iproj = 1, nproj( l )
              write ( 99, '(2(1x,1e15.8),5i4)' ) nbsemel( iproj, m, l, mc ), nproj( l ), iproj, l, m, mc
           end do
        end do
     end do
  end do
  close( unit=99 )
  !
  return
end subroutine confmel
