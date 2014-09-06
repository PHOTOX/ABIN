!---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
!
!      set of random number routines extracted from toolbox
!                                          m. lewerenz may/2001
!-----------------------------------------------------------------------
module mod_random
!---------------------------- ranlfg.in! -------------------------------
!---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
!
!      parameters for lagged fibonacci generators and  generator state
!
!-----------------------------------------------------------------------
!
!      possible (np,nq) values, (np,np-nq) is also valid:
!              (17,5), (250,103), (521,158), (1279,418),
!              (2281,715), (4423,1393), (1279,1063)
!      ref.: Bhanot et al., phys. rev b 33, 7841 (1986);
!            Zierler, inf. control 15, 67 (1961)
!
!      mersenne prime primitive trinomials:
!      (heringa et al. int.j.mod.phys.c 3, 561 (1992))
!
!      (89,38)
!      (127,1), (127,7), (127,15), (127,30), (127,63)
!      (521,32), (521,48), (521,158), (521,168)
!      (607,105), (607,147), (607, 273)
!      (1279,216), (1279,418)
!      (2281,715), (2281,915), (2281,1029)
!      (3217,67), (3217,576)
!      (4423,271), (4423,369), (4423,370), (4423,649), (4424,1393),
!                  (4423,1419), (4423,2098)
!      (9689,84), (9689,471), (9689,1836), (9689,2444), (9689,4187)
!      (19937,881), (19937,7083), (19937,9842)
!      (23209,1530), (23209,6619), (23209,9739)
!      (44497,8575), (44497,21034)
!      (110503,25230), (110503,53719)
!      (132049,7000), (132049,33912), (132049,41469), (132049,52549),
!                     (132049,54454)
!
!      another pair from brent92 who recommends q=0.618p : (258,175) 
!      brent's ranu4 uses (132049,79500)

!-----------------------------------------------------------------------
!     parameter (np=250,nq=103)
!     parameter (np=2281,nq=715)
!     parameter (np=4423,nq=1393)
      use mod_const, only: DP
      private
      public    :: gautrg,vranf,rsavef
      integer,parameter :: np=1279,nq=418
      real(DP)  :: x(np)
      integer   :: last,init
      !gautrg variables determining its state
      real(DP)  :: gsave
      integer   :: isave=-1
      save
      contains
!---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
!----------------------------- last line -------------------------------
      subroutine gautrg(gran,nran,iseed,iout)
         !DH WARNING: initialiazation from gautrg and vranf are different
         !gautrg gran is shifted by 1
!
!      vectorized portable unit standard deviation gaussian random
!      number generator using the box-mueller transformation method.
!      this method is faster than the central limit method, when the
!      uniform random numbers are comparatively expensive.
!      version working on a square with sine and cosine functions. the
!      same sequence is produced independent of the sequence of calls
!      if no calls to other vranf based prn generators are made.
!
!      gran  : vector of length nran returning the random numbers
!      nran  : number of desired random numbers, with nran=0 and iseed
!              not 0 only generator initialization is done. no action
!              when both are zero.
!      iseed : if not 0, integer to start generator. use 0 to continue
!              a previously used/initialized random sequence, unchanged
!              on output.
!      iout  : unit number for messages.
!
!      times for the generation of 10**6 prn`s:
!      rs6000/350 41 mhz    2.31 s 
!      r8000      75 mhz    1.43 s
!      axp21064  200 mhz    1.20 s
!      axp21164  250 mhz    0.572 s
!      t3e-900   450 mhz    0.52 s
!      cray-t90  450 mhz    0.0728 s   1.28*10**7 prn/s
!
!      subroutines called: vranf, r1mach
!      m. lewerenz 6/may/90, modified 17/jun/91, mar/95
!
      implicit real(DP) (a-h,o-z)
      real(DP),parameter  :: two=2.d0,twom=-two,one=1.d0
      integer,parameter :: npoly=11
      real(DP) :: gran(nran),scf(npoly),xran(1),ccf(npoly)
      real(DP) :: cx,sx,y
!     save isave,gsave,tiny,twopi,pi4,scf,ccf
      save tiny,twopi,pi4,scf,ccf
!     data isave/-1/

!      POLYNOMIAL FROM CHEBYSHEV APPROXIMATION ON [ 0.000, 0.790]
!      FOR COS(X) WITH ABSOLUTE ERROR LESS THAN 0.2220E-14

!      real(DP),parameter,dimension(npoly) ::  ccf = &
       DATA ccf/ 0.9999999999999986D+00, 0.6612476846390664D-13, &
               -0.4999999999989523D+00,-0.5434658088910759D-10, &
                0.4166666737609693D-01,-0.4648977428396692D-08, &
               -0.1388871052129944D-02,-0.4228394738587799D-07, &
                0.2486361945804866D-04,-0.5317743184071916D-07, &
               -0.2539224676809412D-06/

!      POLYNOMIAL FROM CHEBYSHEV APPROXIMATION ON [ 0.000, 0.790]
!      FOR SIN(X) WITH ABSOLUTE ERROR LESS THAN 0.2220E-14

      DATA scf/-0.9443414574112906D-15, 0.1000000000000244D+01, &
               -0.1224236196202217D-10,-0.1666666664242968D+00, &
               -0.2471495821870120D-08, 0.8333348067492644D-02, &
               -0.5482536616811601D-07,-0.1982815612039858D-03, &
               -0.2017619095413939D-06, 0.2948964053761139D-05, &
               -0.1051448397925916D-06/

      if(isave.lt.0) then
        isave=0
        tiny=r1mach(0)**2
        pi4=atan(one)
        twopi=pi4*8
      end if
      if(iseed.ne.0) call vranf(gran,0,iseed,iout)

      if(nran.gt.0) then
        newran=nran-isave
        if(isave.eq.1) gran(1)=gsave
        call vranf(gran(isave+1),newran,0,iout)
        do 100 i=1,newran-1,2
        fac=sqrt(twom*log(gran(isave+i)+tiny))
        !DHmod: rename x to y,so that in does not conflict with x(np) from module
        y=pi4*gran(isave+i+1)
        cx=(((((((((ccf(11)*y+ccf(10))*y+ccf(9))*y+ccf(8))*y &
                            +ccf(7))*y+ccf(6))*y+ccf(5))*y  &
                            +ccf(4))*y+ccf(3))*y+ccf(2))*y+ccf(1)
        sx=(((((((((scf(11)*y+scf(10))*y+scf(9))*y+scf(8))*y &
                            +scf(7))*y+scf(6))*y+scf(5))*y  &
                            +scf(4))*y+scf(3))*y+scf(2))*y+scf(1)

        sxi=(two*cx)*sx
        cxi=(two*cx)*cx-one
        sxi=(two*cxi)*sxi
        cxi=(two*cxi)*cxi-one
        sxi=(two*cxi)*sxi
        cxi=(two*cxi)*cxi-one
        gran(isave+i)=fac*sxi
        gran(isave+i+1)=fac*cxi
  100   continue

        if(mod(newran,2).eq.1) then
          call vranf(xran,1,0,iout)
          fac=sqrt(twom*log(gran(nran)+tiny))
          trig=twopi*xran(1)
          gran(nran)=fac*sin(trig)
          gsave=fac*cos(trig)
          isave=1
        else
          isave=0
        end if
      end if
      return
      end subroutine gautrg
!-----------------------------------------------------------------------
!-------------------- ranf/vranf uniform random package ----------------
!-----------------------------------------------------------------------

      subroutine vranf(ranv,nran,iseed,iout)

!      machine independent portable uniform random number generator for
!      the interval [0,1) based on a floating point subtractive lagged
!      fibonacci method similar to the feedback shift register method
!      proposed by kirkpatrick/stoll. subtractive (or additive) version
!      gives much better prn quality than the original xor operation.
!      an additive variant is given in v r a n f a.

!      v r a n f  and r a n f use the same method and can be used to
!      work together on the same random sequence. either of them can be
!      used for initialization. the state of the generator can be saved
!      or retrieved for restart with routine r s a v e f.
! ----------------------------------------------------------------------

!      ranv  : vector of length nran for random numbers; output
!      nran  : number of desired random numbers. nran=0 and iseed.ne.0
!              -> only generator initialization. no action when both 
!              are zero.; input
!      iseed : if not 0, integer to start generator. use 0 to continue
!              a previously used/initialized random sequence, unchanged
!              on output; input
!      iout  : unit number for error messages. silent for iout.le.0
!
! ----------------------------------------------------------------------
!      method:
!      x(k+np)=x(k)-x(k+np-nq), initial x array generated by xuinit.
!      ieee standard requires real(DP) to have at least 48 mantissa bits.
!      with nbit=48 this generator is entirely machine independent
!      and will always give the same random sequence. you can change
!      the period of the generator by setting a different nbit value or
!      changing np and nq appropriately.

! ----------------------------------------------------------------------
!      this is a floating implementation of a generator described in:
!      M.H. KALOS, P.A. WHITLOCK, MONTE CARLO METHODS, APPENDIX,
!                  WILEY 1986
!      D.W. HEERMANN, COMPUTER SIMULATION METHODS, 2ND ED.,SPRINGER 199
!                  APPENDIX A1
!      d. stauffer, f.w. hehl, v. winkelmann, j.g. zabolitzky,
!                  computer simulation & computer algebra, section 2.2

!      original references:
!      s. kirkpatrick, e.p. stoll, j. comput. phys. 40, 517 (1981)
!      r.c. tausworthe, random numbers generated by linear recurrence
!                  modulo 2, math. comp. 19, 201 (1965)
!      t.g. lewis, w.h. payne, generalized feedback shift register
!                  pseudorandom number algorithm, j. acm 20, 456 (1973)

! ----------------------------------------------------------------------
!      other (np,nq) values: (17,5), (250,103), (521,158), (1279,418),
!                            (2281,715), (4423,1393), (1279,1063)
!          ref.: Bhanot et al., phys. rev b 33, 7841 (1986);
!                Zierler, inf. control 15, 67 (1961)
! ----------------------------------------------------------------------
!      alternative additive formulation bypassing if statements:
!      temp=x(k)+x(k+np-nq)
!      x(k)=temp-float(int(temp))

!      alternative subtractive formulation bypassing if statements:
!      temp=x(k)-x(k+np-nq)
!      x(k)=(temp+one)-float(int(temp+one))
! ----------------------------------------------------------------------

!      timing in s for 10**7 random numbers:

!      machine           mhz   unroll(4)     no unrolling
!      ibm rs6000/350     41      2.67            3.9
!      dec axp3800       200      1.17
!      dec axp3600       150      1.35
!      dec alpha 21164   250      0.45
!      t3e alpha 21164   450      0.40
!      sun ultra1        ???      1.38
!      ibm 3090vf         58      ---        2.47(vec), 4.57(sc)
!      cray t90          450      0.31      ranf() takes 0.13 s
!      
!      unrolling to depth 6 gives a slight speed increase.
!----------------------------------------------------------------------
!      subroutines called: xuinit,errprt    m. lewerenz may/91 & nov/93

      implicit real(DP) (a-h,o-z)
      parameter (nratio=np/nq,nexec=4,mroll=4,zero=0.d0,one=1.d0)
      dimension ranv(nran)
      common /doctrl/ nroll

!      table initialization by xuinit

      if(iseed.ne.0) then
        call xuinit(x,np,nq,0,nexec,iseed,init,last,iout)
      end if

!      fibonacci generator updates elements of x in a cyclic fashion
!      and copies them into ranv in blocks of max. length np.
!      loop split into chunks of max. length nq to avoid recurrence.
!      unrolling improves performance on superscalar machines.

      if(nran.gt.0) then
        if(init.ne.0) then
          j=0
          left=nran
   10     continue
          if(nroll.gt.1) then
            loop=mod((min(nq,left+last)-last),mroll)
          else
            loop=min(nq,left+last)-last
          end if

          do 500 i=last+1,last+loop
          x1=x(i)-x(i+np-nq)
          if(x1.lt.zero) x1=x1+one
          x(i)=x1
          j=j+1
          ranv(j)=x1
  500     continue
          if(nroll.gt.1) then
            do 501 i=last+loop+1,min(nq,left+last),mroll
            x1=x(i)-x(i+np-nq)
            x2=x(i+1)-x(i+1+np-nq)
            x3=x(i+2)-x(i+2+np-nq)
            x4=x(i+3)-x(i+3+np-nq)
            if(x1.lt.zero) x1=x1+one
            if(x2.lt.zero) x2=x2+one
            if(x3.lt.zero) x3=x3+one
            if(x4.lt.zero) x4=x4+one
            x(i)=x1
            x(i+1)=x2
            x(i+2)=x3
            x(i+3)=x4
            ranv(j+1)=x1
            ranv(j+2)=x2
            ranv(j+3)=x3
            ranv(j+4)=x4
            j=j+4
  501       continue
          end if

          if(last.lt.nratio*nq) then
            do 650 k=1,nratio-1
            limit=min((k+1)*nq,left+last)
            if(nroll.gt.1) then
              loop=mod((limit-max(k*nq,last)),mroll)
            else
              loop=limit-max(k*nq,last)
            end if

            do 600 i=max(k*nq,last)+1,max(k*nq,last)+loop
            x1=x(i)-x(i-nq)
            if(x1.lt.zero) x1=x1+one
            x(i)=x1
            j=j+1
            ranv(j)=x1
  600       continue
            if(nroll.gt.1) then
              do 601 i=max(k*nq,last)+loop+1,limit,mroll
              x1=x(i)-x(i-nq)
              x2=x(i+1)-x(i+1-nq)
              x3=x(i+2)-x(i+2-nq)
              x4=x(i+3)-x(i+3-nq)
              if(x1.lt.zero) x1=x1+one
              if(x2.lt.zero) x2=x2+one
              if(x3.lt.zero) x3=x3+one
              if(x4.lt.zero) x4=x4+one
              x(i)=x1
              x(i+1)=x2
              x(i+2)=x3
              x(i+3)=x4
              ranv(j+1)=x1
              ranv(j+2)=x2
              ranv(j+3)=x3
              ranv(j+4)=x4
              j=j+4
  601         continue
            end if
  650       continue
          end if

          limit=min(np,left+last)
          if(nroll.gt.1) then
            loop=mod((limit-max(nratio*nq,last)),mroll)
          else
            loop=limit-max(nratio*nq,last)
          end if

          do 700 i=max(nratio*nq,last)+1,max(nratio*nq,last)+loop
          x1=x(i)-x(i-nq)
          if(x1.lt.zero) x1=x1+one
          x(i)=x1
          j=j+1
          ranv(j)=x1
  700     continue
          if(nroll.gt.1) then
            do 701 i=max(nratio*nq,last)+loop+1,limit,mroll
            x1=x(i)-x(i-nq)
            x2=x(i+1)-x(i+1-nq)
            x3=x(i+2)-x(i+2-nq)
            x4=x(i+3)-x(i+3-nq)
            if(x1.lt.zero) x1=x1+one
            if(x2.lt.zero) x2=x2+one
            if(x3.lt.zero) x3=x3+one
            if(x4.lt.zero) x4=x4+one
            x(i)=x1
            x(i+1)=x2
            x(i+2)=x3
            x(i+3)=x4
            ranv(j+1)=x1
            ranv(j+2)=x2
            ranv(j+3)=x3
            ranv(j+4)=x4
            j=j+4
  701       continue
          end if

          last=mod(limit,np)
          left=nran-j
          if(left.gt.0) goto 10
        else
          call errprt(iout,'vranf','incorrect initialization',-1)
        end if
      end if
      return
      end subroutine vranf

!-----------------------------------------------------------------------

      subroutine xuinit(x,np,nq,mode,nexec,iseed,init,last,iout)

!      initializes a (np,nq) lagged fibonacci generator table
!      with random bits generated by a congruential generator using
!      l'ecuyers decomposition. ref.: bratley et al. p. 214

!      x     : vector of length np for initial random number table;
!              output
!      np,nq : parameters p and q of feedback shift register generator;
!              input
!      mode  : operation for lfg:
!              mode=<0 -> subtractive generator, mode=1 additive
!      nexec : number of warm up cycles for the table. nexec*nbit*np
!              random numbers are generated and discarded; input
!      iseed : integer seed for congruential generator generating
!              the bits of the initial table entries; input
!      init  : returns updated seed of congruential generator.
!              0 if table was not initialized, > 0 if ok; output
!      last  : pointer to the last used number in the table; output
!      iout  : unit number for messages, silent for iout.le.0; input
!      subroutines called : none             m. lewerenz mar/93, mar/98

      implicit real(DP) (a-h,o-z)
      parameter (ia=40692,ib=52774,ic=3791,ip=2147483399)
      parameter (zero=0.d0,half=0.5d0,iphalf=ip/2,nbit=48)
      logical high
      dimension x(np)

      if(nq.ge.np.or.iseed.eq.0) then
        call errprt(iout,'xuinit','illegal seed parameter(s)',-1)
      else

!      set table to zero and exercise the bit generator a little

        ix=iabs(iseed)
        if(ix.ne.0) then
          do i=1,np
            x(i)=zero
            k1=ix/ib
            ix=ia*(ix-k1*ib)-k1*ic
            if(ix.lt.0) ix=ix+ip
          end do

!      assemble np numbers with mantissa length nbit from random bits
!      'high' toggle compensates for bias from odd ip

          high=.true.
          do i=1,np
            add=half
            do j=1,nbit
              k1=ix/ib
              ix=ia*(ix-k1*ib)-k1*ic
              if(ix.lt.0) ix=ix+ip
              if(high) then
                if(ix.ge.iphalf) x(i)=x(i)+add
                high=.false.
              else
                if(ix.gt.iphalf) x(i)=x(i)+add
                high=.true.
              end if
              add=add*half
            end do
          end do
          if(nexec.gt.0) call xuwarm(x,np,nq,mode,nbit*nexec,iout)
        end if
        init=ix
        last=0
      end if
      return
      end subroutine xuinit

!-----------------------------------------------------------------------

      subroutine xuwarm(x,np,nq,mode,nexec,iout)

!      warms up (p,q) lagged fibonacci generators (lfg) by nexec rounds

!      x     : vector of length np for initial random number table;
!              output
!      np,nq : parameters p and q of feedback shift register generator;
!              input
!      mode  : operation for lfg:
!              mode=<0 -> subtractive generator, mode=1 additive
!      nexec : number of warm up cycles for the table. nexec*nbit*np
!              random numbers are generated and discarded; input
!      iout  : unit number for messages, silent for iout.le.0; input
!      subroutines called : errprt                   m. lewerenz mar/98

      implicit real(DP) (a-h,o-z)
      parameter (zero=0.d0,one=1.d0)
      dimension x(np)

      if(nq.ge.np.or.np.eq.0.or.nq.eq.0) then
        call errprt(iout,'xuwarm','illegal table parameter(s)',-1)
      else

!      exercise the generator for nexec rounds of np prn's
!      separate sections for subtractive or additive version

          if(mode.le.0) then
            do k=1,nexec
              do i=1,nq
                x(i)=x(i)-x(i+np-nq)
                if(x(i).lt.zero) x(i)=x(i)+one
              end do
              do i=nq+1,np
                x(i)=x(i)-x(i-nq)
                if(x(i).lt.zero) x(i)=x(i)+one
              end do
            end do
          else
            do k=1,nexec
              do i=1,nq
                x(i)=x(i)+x(i+np-nq)
                if(x(i).ge.one) x(i)=x(i)-one
              end do
              do i=nq+1,np
                x(i)=x(i)+x(i-nq)
                if(x(i).ge.one) x(i)=x(i)-one
              end do
            end do
          end if
      end if
      return
      end subroutine xuwarm

!---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
!----------------------------------------------------------------------
!---+----|----+----|----+----|----+----|----+----|----+----|----+----|--

      FUNCTION R1MACH (IDUM)
      IMPLICIT real(DP) (A-H,O-Z)
      PARAMETER (ONE=1.D0,TWO=2.D0,HALF=0.5D0)
      SAVE ICALL,EPS
      DATA ICALL,EPS/0,ONE/
! ---------------------------------------------------------------------
! THIS ROUTINE COMPUTES THE UNIT ROUNDOFF OF THE MACHINE.
! THIS IS DEFINED AS THE SMALLEST POSITIVE MACHINE NUMBER
! U SUCH THAT 1.0 + U .NE. 1.0E0
! ---------------------------------------------------------------------
      IF(ICALL.EQ.0) THEN
        ICALL=1
        U = ONE
  10    U = U*HALF
        COMP = ONE + U
        IF(COMP .NE. ONE) GO TO 10
        EPS = U*TWO
      END IF
      R1MACH = EPS
      RETURN
      END FUNCTION R1MACH

!-----------------------------------------------------------------------
!-------------------------- error handling -----------------------------
!-----------------------------------------------------------------------

      subroutine errprt(iout,pgname,text,icode)

!      prints error messages from library subroutines

!      iout   : unit number for message output, 0-> no output; input
!      pgname : name of the subroutine calling errprt; input
!      text   : message text; input
!      icode  : severity code: 0 -> warning, < 0 -> fatal error with
!               abort, else -> error but execution continues
!      subroutines called : strlen                   m. lewerenz dec/93

      character pgname*(*),text*(*),header*20,tail*40
      save nerror,nwarn,icall
      common /errcnt/ maxerr,maxwrn
      data icall/0/

      if(icall.eq.0) then
        icall=1
        nerror=0
        nwarn=0
      end if

      if(icode.lt.0) then
        header='  *** fatal error,'
        tail=', execution aborted ***'
      else if(icode.eq.0) then
        header='  *** warning,'
        tail=' ***'
        nwarn=nwarn+1
      else
        header='  *** error,'
        tail=', return without action ***'
        nerror=nerror+1
      end if

!      write the message on unit iout

      if(iout.gt.0) then
        call getstr(pgname,lname,iout)
        call getstr(text,ltext,iout)
        call getstr(header,lhead,iout)
        call getstr(tail,ltail,iout)
        write(iout,'(/6a/)') header(1:lhead),' ',text(1:ltext),' in ', &
     &                      pgname(1:lname),tail(1:ltail)
        call flush(iout)
      end if

      jcode=icode
      if(maxerr.gt.0.and.nerror.ge.maxerr) then
        if(iout.gt.0) write(iout,'(/a)')  &
     &  '  *** maximum number of errors exceeded, program stopped *** '
        jcode=-1
      end if
      if(maxwrn.gt.0.and.nwarn.ge.maxwrn) then
        if(iout.gt.0) write(iout,'(/a)') &
     &  '  *** maximum number of warnings exceeded, program stopped ***'
        jcode=-1
      end if
      if(iout.gt.0) call flush(iout)
!
      if(jcode.lt.0) stop
      return
!
! ---------------------------------------------------------------------
!      error report, returns current number of errors and warning 
!
      entry errnum(nerr,nwrn)
      nerr=nerror
      nwrn=nwarn
      return
      end subroutine errprt
!
!---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
!----------------------------------------------------------------------
!---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
!
      subroutine getstr(string,ls,iout)
!
!      eliminates leading and trailing blanks from string and returns
!      length of remaining string
!
!      string : character string; input & output
!               on output leading and trailing blanks of the original
!               string are missing
!      ls     : length of non blank portion of string; output
!      iout   : fortran unit for messages, silent if iout.le.0; input
!      subroutines called: errprt                    m. lewerenz jul/00
!
      character string*(*)
!
      ls=len(string)
      if(ls.gt.0) then
        i=0
    5   i=i+1
        if(i.le.ls) then
          if(string(i:i).eq.' ') then
            goto 5
          else
            ls=ls+1-i
            if(i.gt.1) then
              do j=1,ls
                string(j:j)=string(j+i-1:j+i-1)
              end do
            end if
   10       continue
            if(string(ls:ls).eq.' ') then
              ls=ls-1
              goto 10
            end if
          end if
        else
          ls=0
          call errprt(iout,'getstr','empty string',0)
        end if
      end if
      return
      end subroutine getstr

      subroutine rsavef(iout,lread)
         use mod_utils, only: abinerror
         integer,intent(in) :: iout  !where do we write the state
!         integer,intent(in) :: isave !0=read,1=save
         logical,intent(out),optional :: lread
         character(len=*),parameter   :: chprng='PRNG STATE (OPTIONAL)'
         character(len=50)  :: readstring
         integer            :: iost,i
         logical            :: lexist,lopened

         if(present(lread)) lread=.false.

         inquire(unit=iout,exist=lexist,opened=lopened)
         if(lexist.and..not.lopened)then
            write(*,*)'Unit for PRNG state must be opened!Exiting...'
            call abinerror('rsavef')
         end if

!         if (isave.eq.1)then
!         if (isave.eq.1)then
!            write(iout,'(A)')chprng
!            write(iout,*)init,last
!            do i=1,np
!               write(iout,*)x(i)
!            end do
!         end if

         if (present(lread))then
            read(iout,'(A)',iostat=iost)readstring
            if(iost /= 0)then
               write(*,*)'PRNG state not present in restart.xyz. Ignoring...'
               write(*,*)'Random seed from input parameter file will be used.'
            else if (chprng.ne.trim(readstring)) then
               write(*,*)'ERROR: PRNG STATE IN RESTART FILE SEEMS TO BE BROKEN.'
               write(*,*)'Expected:',chprng
               write(*,*)'Got',readstring
               call abinerror('rsavef')
            else
               read(iout,*)init,last
               read(iout,*)isave,gsave
               do i=1,np
                  read(iout,*)x(i)
               end do
               lread=.true.
            end if
         else
            write(iout,'(A)')chprng
            write(iout,*)init,last
            write(iout,*)isave,gsave
            do i=1,np
               write(iout,*)x(i)
            end do
         end if
      end subroutine rsavef
           
!---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
!---------------------------- last line --------------------------------
end module mod_random
