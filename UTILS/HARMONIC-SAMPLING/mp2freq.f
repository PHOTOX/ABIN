!  this program will take a portion of the molpro output file
!  and convert it to a Frequency.dat file
!  input format:
!  on the first line put the number of atoms, the number of modes,
!  and the number of the first mode you want to include in the 
!  Frequency.dat file (1 to include all modes)
!  Then, without skipping a line append the list of frequencies and 
!  modes from the molpro output, starting with the line 
!  "Normal Modes"
      program mp2freq
      integer MaxNumAtoms
      parameter (MaxNumAtoms=100)
      integer NumAtoms, NumModes, i, j, k, istartpt
      real*8 freq(MaxNumAtoms*3-6)
      real*8 modes(MaxNumAtoms*3-6,MaxNumAtoms*3)
      character*16 junk, junk2
      character*256 junk3

 1003 format(i10,i10)
 1004 format(4f20.10)

      read(5,*)NumAtoms,NumModes
      istartpt=1 
10    read(5,*) junk
      if (junk.eq.'Normal         ') goto 20
      goto 10
 
20    continue 
      write(*,1003) NumAtoms, NumModes-istartpt+1

      do i=1, (floor(dble(NumModes)/5.0d0))
         read(5,*) junk3
         read(5,*) junk, junk2,(freq(k), k=5*i-4,5*i)
         read(5,*) junk3
         read(5,*) junk3
         do j=1, (NumAtoms*3)
            read(5,*) junk, (modes(k,j), k=5*i-4,5*i)
         enddo
      enddo


      if (ceiling(dble(NumModes)/5.0d0) .gt. (i-1)) then
         i=ceiling(dble(NumModes)/5.0d0)
         read(5,*) junk3
         read(5,*) junk, junk2, (freq(k), k=5*i-4,NumModes)
         read(5,*) junk3
         read(5,*) junk3
         do j=1, (NumAtoms*3)
            read(5,*) junk, (modes(k,j), k=5*i-4,NumModes)
         enddo
      endif

      NumModes=NumModes-istartpt +1
      do i=1, NumModes
         freq(i)=freq(i+istartpt-1)
         do j=1, NumAtoms*3
            modes(i,j)=modes(i+istartpt-1,j)
         enddo
      enddo

      do i=1, (int(floor(dble(NumModes)/4.0d0)))
         write(*,1004) (freq(k), k=4*i-3,4*(i))
      enddo
      i=int(ceiling(dble(NumModes)/4.0d0))
      if (4*i.ne.NumModes) write(*,1004) (freq(k), k=4*i-3,NumModes)

      do k=1, NumModes
         do i=1, (floor(dble(NumAtoms*3)/4.0d0))
            write(*,1004) (modes(k,j), j=4*i-3, 4*i)
         enddo
         i=int(ceiling(dble(NumAtoms*3)/4.0d0))
         if (4*i.ne.NumAtoms*3) write(*,1004) (
     $        modes(k,j), j=4*i-3, NumAtoms*3)
      enddo

      end program
