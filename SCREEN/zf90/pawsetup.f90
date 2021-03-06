! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
!  program pawsetup
!
!  creates the site list based upon the inputs
!
!TODO This should really be moved to a subroutine 
program pawsetup
      use ocean_phys, only : ophys_fixCoords
      use periodic
      implicit none
!
      integer :: ntypat, natom, nedges, counter, counter2, uniquepsp, ierr
      integer, allocatable :: znucl(:), typat(:), edges(:,:), sites(:), &
     &    atomcount(:), sitenum(:), pspused(:),inppopts(:), inppfill(:)
      real(kind=kind(1.d0)) :: avecs(3,3)
      real(kind=kind(1.d0)), allocatable :: xred(:,:), xinput(:,:), xcoord(:,:)
      character(len=50), allocatable :: pplist(:),ppopts(:),ppfill(:)
      character(len=50) :: atemp1, atemp2
      character(len=4) :: smode
      character(len=8) :: coordFormat
      logical :: ex, val
!
      character(len=9), parameter :: f9='formatted'
!
      open(unit=99,file='ntype',form=f9,status='old')
      read(99,*) ntypat
      close(99)
!
      allocate(znucl(ntypat))
      open(unit=99,file='znucl',form=f9,status='old')
      read(99,*) znucl(:)
      close(99)
!
      open(unit=99,file='natoms',form=f9,status='old')
      read(99,*) natom
      close(99)
!
      allocate(typat(natom))
      open(unit=99,file='typat',form=f9,status='old')
      read(99,*) typat(:)
      close(99)
!
      allocate(xinput(3,natom), xcoord(3, natom), xred(3,natom))
      open(unit=99,file='taulist',form=f9,status='old')
      read(99,*) xinput(:,:)
      close(99)
!
      open(unit=99,file='coord',form=f9,status='old')
      read(99,*) coordFormat
      close(99)
!
      open(unit=99,file='avecsinbohr.ipt',form=f9,status='old')
      read(99,*) avecs(:,:)
      close(99)
!
      call ophys_fixCoords( avecs, coordFormat, xinput, xred, xcoord, ierr )
!
      open(unit=99,file='nedges',form=f9,status='old')
      read(99,*) nedges
      close(99)
!
      allocate(edges(3,nedges))
      open(unit=99,file='edges',form=f9,status='old')
      read(99,*) edges(:,:)
      close(99)
!
! For each element in the file this will write out the element symbol and the
!  reduced coordinate location
      open(unit=99,file='xyz.wyck',form=f9,status='unknown')
      write(99,*) natom
      do counter=1,natom
        write(99,'(A3,3(1X,F14.10))') elements(znucl(typat(counter))),   &
     &                               xred(:,counter)
      enddo
      close(99)

      val = .false.
      inquire(file='screen.mode',exist=ex)
      if( ex ) then
        open(unit=99, file='screen.mode', form=f9, status='old' )
        read(99,*) smode
        close( 99 )
        if( smode .eq. 'grid' ) val = .true.
      endif
      if( val ) goto 111
!
      allocate(sites(natom),sitenum(natom),pspused(ntypat))
      sites(:) = 0
      sitenum(:) = 0
      pspused(:) = 0
! For each edge listed, flags the site that it corresponds to
      do counter=1,nedges
        sites(edges(1,counter)) = 1 !sites(edges(1,counter)) + 1
        !sites(edges(1,counter)) = sites(edges(1,counter)) + 1
        pspused(typat(edges(1,counter))) = 1
      enddo
!      
!
! Gets the psp names for all the elements      
      allocate(pplist(ntypat))
      open(unit=99,file='pplist',form=f9,status='old')
      read(99,*) pplist(:)
      close(99)
!
! Gets the psp opts
      uniquepsp = sum(pspused(:))
! apr20
!      allocate(ppopts(uniquepsp), ppfill(uniquepsp), inppopts(uniquepsp)&
!     &        ,inppfill(uniquepsp))
      allocate(ppopts(ntypat), ppfill(ntypat), inppopts(ntypat)&
     &        ,inppfill(ntypat))
      inppopts(:) = 0
      inppfill(:) = 0
      ppopts(:) = ''
      ppfill(:) = ''
      open(unit=99,file='opf.opts',form=f9,status='old')
      do counter=1, ntypat
        read(99,*,end=100,err=100) inppopts(counter), ppopts(counter)
      enddo
  100 continue
      close(99)
!
! Gets the psp fillinpaw details
      open(unit=99,file='opf.fill',form=f9,status='old')
      do counter=1, ntypat
        read(99,*,end=101,err=101) inppfill(counter), ppfill(counter)
      enddo
  101 continue
      close(99)
!
! Write out a list of the psp for which cores will be run 
      open(unit=99,file='pspoptions',form=f9,status='unknown')
      do counter=1,ntypat
        if (pspused(counter) .ne. 0 ) then
          do counter2=1, ntypat
            if (inppopts(counter2) .eq. znucl(counter) )                   &
     &            atemp1 = ppopts(counter2)
            if (inppfill(counter2) .eq. znucl(counter) )                   &
     &            atemp2 = ppfill(counter2)
          enddo
          write(99,'(I3,1X,A2,3(1X,A50))') znucl(counter),                &
     &      elements(znucl(counter)), pplist(counter), atemp1, atemp2
        endif
      enddo
      close(99)
! Slightly convoluted, but ...
! Gives the total number of unique sites to consider
! Then lists each site, giving the symbol and its rank in the xyz.wyck file
      open(unit=99,file='sitelist',form=f9,status='unknown')
      write(99,*) sum(sites)
!      allocate(atomcount(ntypat))
      ! Needs to be more than possible Z
      allocate(atomcount(160))
      atomcount(:) = 0
      do counter=1,natom
        atomcount(znucl(typat(counter))) = atomcount(znucl(typat(counter))) + 1
        if (sites(counter) .ne. 0 ) then
          write(99,*) elements(znucl(typat(counter))), znucl(typat(counter)),  &
     &               atomcount(znucl(typat(counter))), pplist(typat(counter))
        endif
        sitenum(counter) = atomcount(znucl(typat(counter)))
      enddo
      close(99)
!
!
! Writes out each psp name, Z, n, l, Symbol, rank (in xyz.wyck)
!  for each edge to be considered.
      open(unit=99,file='hfinlist',form=f9,status='unknown')
      do counter=1,nedges
        write(99,'(A50,3(1X,I3),1X,A2,1X,I4)')pplist(typat(edges(1,counter))),&
     &            znucl(typat(edges(1,counter))),  edges(2:3,counter),  &
     &  elements(znucl(typat(edges(1,counter)))), sitenum(edges(1,counter))
      enddo
      close(99)
!
      deallocate(znucl,typat,xred,xinput,xcoord,edges,sites,atomcount, pplist,sitenum)
111   continue
!
end program pawsetup 
!!!!!!!!
