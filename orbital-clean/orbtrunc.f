!------------------------------------------------------------------------
      program orbtrunc
      implicit none
      integer norbs,nbasis,ntrunc,nverts,natoms

      call get_params(norbs,natoms,nbasis,ntrunc,nverts)
      write(5,9000)norbs,nbasis,natoms,ntrunc,nverts     

      call trunc_mov_invert(norbs,nbasis,ntrunc,nverts,natoms)

 9000 format(/5x,'total of',i5,' orbitals created from',i5,' basis fns',
     $     ' centered on',i5,' atoms'/
     $    20x,'           truncating:',i5,' orbitals'/
     $    20x,'copying and inverting:',i5,' orbitals'/)
      end program orbtrunc
!------------------------------------------------------------------------
      subroutine get_params(norbs,natoms,nbasis,ntrunc,nverts)
      implicit none
      integer norbs,nbasis,ntrunc,nverts,natoms
! input file trunc_params is as follows:
! first line: #orbs #atoms #basis_fns #trucated_orbs #move_and_inverted_orbs
! next #trunc lines: orbital_id iatom1 iatom2
! next #nvert lines: iorb1 iorb2 iatom1(2)      

      open(1,status='unknown',file='trunc_params')
      read(1,*)norbs,natoms,nbasis,ntrunc,nverts
      close(1)
      return
      end subroutine get_params
!------------------------------------------------------------------------
      subroutine trunc_mov_invert(norbs,nbasis,ntrunc,nverts,natoms)
      implicit none
!passed in variables
      integer norbs,nbasis,ntrunc,nverts,natoms
!local storage/vaiables      
      character*4 symbol
      integer iorb,ibasis,print_upto,junk1,junk2,iatom,ratom,itemp,
     $     molvec,iatom1,iatom2,itrunc,ibasis1,ibasis2,iorb1,iorb2,
     $     ivert,ic,max,min,modj,modic,j,i,ijunk
      integer orb_id(natoms,2),orb_trunc(ntrunc+nverts,3)
      double precision vects(norbs,nbasis+5),temp_vect(nbasis+5)
      double precision v0,v1,v2,v3,v4

! read vectors into memory
      OPEN(1,STATUS='UNKNOWN',FILE='orbital_vects')      
      do iorb =1,norbs
          do ibasis = 1,nbasis,5
             read(1,9000)junk1,junk2,v0,v1,v2,v3,v4
             vects(iorb,ibasis)   =v0
             vects(iorb,ibasis+1) =v1
             vects(iorb,ibasis+2) =v2
             vects(iorb,ibasis+3) =v3
             vects(iorb,ibasis+4) =v4
          enddo
      enddo
      close(1)
 9000 format(i2,i3,5e15.8)

! read/determine information on which basis fns belong to which atoms
      OPEN(1,STATUS='UNKNOWN',FILE='orbital_ids')      
      do iatom =1,natoms
         itemp = iatom/100
         itemp = iatom - itemp*100
         do ibasis=1,nbasis
            read(1,9001)ijunk,symbol,ratom
            if(ratom.eq.itemp)then
               orb_id(iatom,1)=ijunk
               goto 10
            endif
         enddo
 10      continue 
      enddo
      close(1)
! assign last basis fn for each atoms
      do iatom = 1,natoms-1
         orb_id(iatom,2)=orb_id(iatom+1,1)-1
      enddo
      orb_id(natoms,2)=nbasis
 9001 format(i5,2x,a2,i2)

! read in which orbitals and atoms to keep basis fns for
      open(1,status='unknown',file='trunc_params')
      rewind(1)
      read(1,*)
      do itrunc=1,ntrunc
         read(1,*)molvec,iatom1,iatom2
         orb_trunc(itrunc,1)=molvec
         orb_trunc(itrunc,2)=iatom1
         orb_trunc(itrunc,3)=iatom2
      enddo
! read in which orbtials to copy from and atoms to invert on the copy
      do itrunc=ntrunc+1,ntrunc+nverts
         read(1,*)iorb1,iorb2,iatom1
         orb_trunc(itrunc,1)=iorb1
         orb_trunc(itrunc,2)=iorb2
         orb_trunc(itrunc,3)=iatom1
      enddo
      close(1)

! do actual truncation
      do itrunc=1,ntrunc
         iorb = orb_trunc(itrunc,1)
         do ibasis =1,nbasis
            temp_vect(ibasis)=vects(iorb,ibasis)
            vects(iorb,ibasis)=0.0d+00
         enddo
         do iatom = 2,3
            iatom1= orb_trunc(itrunc,iatom)
            ibasis1=orb_id(iatom1,1)
            ibasis2=orb_id(iatom1,2)
            do ibasis = ibasis1,ibasis2
               vects(iorb,ibasis)=temp_vect(ibasis)
            enddo
         enddo
      enddo

! do orbital movement
      do ivert = 1,nverts
         iorb1 = orb_trunc(ntrunc+ivert,1)
         iorb2 = orb_trunc(ntrunc+ivert,2)
         do ibasis = 1,nbasis
            vects(iorb2,ibasis) = vects(iorb1,ibasis)
         enddo
! switch signs for atom1 (creating antibonding orbital)
           iatom1 =orb_trunc(ntrunc+ivert,3)
           ibasis1=orb_id(iatom1,1)
           ibasis2=orb_id(iatom1,2)
         do ibasis = ibasis1,ibasis2
            vects(iorb2,ibasis) =(-1)*vects(iorb2,ibasis)
         enddo
      enddo

! write out final orbitals
      open(1,status='unknown',file='modified_orbitals')
      write(1,'("MODIFIED ORBITALS"/x,"$VEC")')      
      do j = 1,norbs
         ic = 0
         max = 0
 100     continue
         min = max+1
         max = max+5
         ic = ic+1
         if (max .gt. nbasis) max = nbasis
         modj  = mod(j ,100 )
         modic = mod(ic,1000)
         write (1,9002) modj,modic,(vects(j,i),i = min,max)
         if (max.lt.nbasis) go to 100
      enddo
      write(1,'(x,"$END")')
      close(1)
 9002 format(i2,i3,1p,5e15.8)
      return
      end subroutine trunc_mov_invert
!--------------------------------------------------------
