c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
c      set of random number routines extracted from toolbox
c                                          m. lewerenz may/2001
c
c-----------------------------------------------------------------------
c
      subroutine gautrg(gran,nran,iseed,iout)
c
c      vectorized portable unit standard deviation gaussian random
c      number generator using the box-mueller transformation method.
c      this method is faster than the central limit method, when the
c      uniform random numbers are comparatively expensive.
c      version working on a square with sine and cosine functions. the
c      same sequence is produced independent of the sequence of calls
c      if no calls to other vranf based prn generators are made.
c
c      gran  : vector of length nran returning the random numbers
c      nran  : number of desired random numbers, with nran=0 and iseed
c              not 0 only generator initialization is done. no action
c              when both are zero.
c      iseed : if not 0, integer to start generator. use 0 to continue
c              a previously used/initialized random sequence, unchanged
c              on output.
c      iout  : unit number for messages.
c
c      times for the generation of 10**6 prn`s:
c      rs6000/350 41 mhz    2.31 s 
c      r8000      75 mhz    1.43 s
c      axp21064  200 mhz    1.20 s
c      axp21164  250 mhz    0.572 s
c      t3e-900   450 mhz    0.52 s
c      cray-t90  450 mhz    0.0728 s   1.28*10**7 prn/s
c
c      subroutines called: vranf, r1mach
c      m. lewerenz 6/may/90, modified 17/jun/91, mar/95
c
      implicit real*8 (a-h,o-z)
      parameter (two=2.d0,twom=-two,one=1.d0,npoly=11)
      dimension gran(nran),ccf(npoly),scf(npoly),xran(1)
      save isave,gsave,tiny,twopi,pi4,ccf,scf
      data isave/-1/
C
C      POLYNOMIAL FROM CHEBYSHEV APPROXIMATION ON [ 0.000, 0.790]
C      FOR COS(X) WITH ABSOLUTE ERROR LESS THAN 0.2220E-14
C
      DATA ccf/ 0.9999999999999986D+00, 0.6612476846390664D-13,
     #         -0.4999999999989523D+00,-0.5434658088910759D-10,
     #          0.4166666737609693D-01,-0.4648977428396692D-08,
     #         -0.1388871052129944D-02,-0.4228394738587799D-07,
     #          0.2486361945804866D-04,-0.5317743184071916D-07,
     #         -0.2539224676809412D-06/
C
C      POLYNOMIAL FROM CHEBYSHEV APPROXIMATION ON [ 0.000, 0.790]
C      FOR SIN(X) WITH ABSOLUTE ERROR LESS THAN 0.2220E-14
C
      DATA scf/-0.9443414574112906D-15, 0.1000000000000244D+01,
     #         -0.1224236196202217D-10,-0.1666666664242968D+00,
     #         -0.2471495821870120D-08, 0.8333348067492644D-02,
     #         -0.5482536616811601D-07,-0.1982815612039858D-03,
     #         -0.2017619095413939D-06, 0.2948964053761139D-05,
     #         -0.1051448397925916D-06/
c
      if(isave.lt.0) then
        isave=0
        tiny=r1mach(0)**2
        pi4=atan(one)
        twopi=pi4*8
      end if
      if(iseed.ne.0) call vranf(gran,0,iseed,iout)
c
      if(nran.gt.0) then
        newran=nran-isave
        if(isave.eq.1) gran(1)=gsave
        call vranf(gran(isave+1),newran,0,iout)
        do 100 i=1,newran-1,2
        fac=sqrt(twom*log(gran(isave+i)+tiny))
        x=pi4*gran(isave+i+1)
        cx=(((((((((ccf(11)*x+ccf(10))*x+ccf(9))*x+ccf(8))*x
     #                       +ccf(7))*x+ccf(6))*x+ccf(5))*x
     #                       +ccf(4))*x+ccf(3))*x+ccf(2))*x+ccf(1)
        sx=(((((((((scf(11)*x+scf(10))*x+scf(9))*x+scf(8))*x
     #                       +scf(7))*x+scf(6))*x+scf(5))*x
     #                       +scf(4))*x+scf(3))*x+scf(2))*x+scf(1)
        sxi=(two*cx)*sx
        cxi=(two*cx)*cx-one
        sxi=(two*cxi)*sxi
        cxi=(two*cxi)*cxi-one
        sxi=(two*cxi)*sxi
        cxi=(two*cxi)*cxi-one
        gran(isave+i)=fac*sxi
        gran(isave+i+1)=fac*cxi
  100   continue
c
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
      end
c-----------------------------------------------------------------------
c-------------------- ranf/vranf uniform random package ----------------
c-----------------------------------------------------------------------
c
      subroutine vranf(ranv,nran,iseed,iout)
c
c      machine independent portable uniform random number generator for
c      the interval [0,1) based on a floating point subtractive lagged
c      fibonacci method similar to the feedback shift register method
c      proposed by kirkpatrick/stoll. subtractive (or additive) version
c      gives much better prn quality than the original xor operation.
c      an additive variant is given in v r a n f a.
c
c      v r a n f  and r a n f use the same method and can be used to
c      work together on the same random sequence. either of them can be
c      used for initialization. the state of the generator can be saved
c      or retrieved for restart with routine r s a v e f.
c ----------------------------------------------------------------------
c
c      ranv  : vector of length nran for random numbers; output
c      nran  : number of desired random numbers. nran=0 and iseed.ne.0
c              -> only generator initialization. no action when both 
c              are zero.; input
c      iseed : if not 0, integer to start generator. use 0 to continue
c              a previously used/initialized random sequence, unchanged
c              on output; input
c      iout  : unit number for error messages. silent for iout.le.0
c
c ----------------------------------------------------------------------
c      method:
c      x(k+np)=x(k)-x(k+np-nq), initial x array generated by xuinit.
c      ieee standard requires real*8 to have at least 48 mantissa bits.
c      with nbit=48 this generator is entirely machine independent
c      and will always give the same random sequence. you can change
c      the period of the generator by setting a different nbit value or
c      changing np and nq appropriately.
c
c ----------------------------------------------------------------------
c      this is a floating implementation of a generator described in:
C      M.H. KALOS, P.A. WHITLOCK, MONTE CARLO METHODS, APPENDIX,
C                  WILEY 1986
C      D.W. HEERMANN, COMPUTER SIMULATION METHODS, 2ND ED.,SPRINGER 199
C                  APPENDIX A1
c      d. stauffer, f.w. hehl, v. winkelmann, j.g. zabolitzky,
c                  computer simulation & computer algebra, section 2.2
c
c      original references:
c      s. kirkpatrick, e.p. stoll, j. comput. phys. 40, 517 (1981)
c      r.c. tausworthe, random numbers generated by linear recurrence
c                  modulo 2, math. comp. 19, 201 (1965)
c      t.g. lewis, w.h. payne, generalized feedback shift register
c                  pseudorandom number algorithm, j. acm 20, 456 (1973)
c
c ----------------------------------------------------------------------
c      other (np,nq) values: (17,5), (250,103), (521,158), (1279,418),
c                            (2281,715), (4423,1393), (1279,1063)
c          ref.: Bhanot et al., phys. rev b 33, 7841 (1986);
c                Zierler, inf. control 15, 67 (1961)
c ----------------------------------------------------------------------
c      alternative additive formulation bypassing if statements:
c      temp=x(k)+x(k+np-nq)
c      x(k)=temp-float(int(temp))
c
c      alternative subtractive formulation bypassing if statements:
c      temp=x(k)-x(k+np-nq)
c      x(k)=(temp+one)-float(int(temp+one))
c ----------------------------------------------------------------------
c
c      timing in s for 10**7 random numbers:
c
c      machine           mhz   unroll(4)     no unrolling
c      ibm rs6000/350     41      2.67            3.9
c      dec axp3800       200      1.17
c      dec axp3600       150      1.35
c      dec alpha 21164   250      0.45
c      t3e alpha 21164   450      0.40
c      sun ultra1        ???      1.38
c      ibm 3090vf         58      ---        2.47(vec), 4.57(sc)
c      cray t90          450      0.31      ranf() takes 0.13 s
c      
c      unrolling to depth 6 gives a slight speed increase.
c ----------------------------------------------------------------------
c      subroutines called: xuinit,errprt    m. lewerenz may/91 & nov/93
c
      implicit real*8 (a-h,o-z)
      include 'ranlfg.inc'
      parameter (nratio=np/nq,nexec=4,mroll=4,zero=0.d0,one=1.d0)
      dimension ranv(nran)
      common /doctrl/ nroll
c
c      table initialization by xuinit
c
      if(iseed.ne.0) then
        call xuinit(x,np,nq,0,nexec,iseed,init,last,iout)
      end if
c
c      fibonacci generator updates elements of x in a cyclic fashion
c      and copies them into ranv in blocks of max. length np.
c      loop split into chunks of max. length nq to avoid recurrence.
c      unrolling improves performance on superscalar machines.
c
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
convex, cray, and ibm directives
c$dir no_recurrence
cdir$ ivdep
cibmdir ignore recrdeps
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
c
          if(last.lt.nratio*nq) then
            do 650 k=1,nratio-1
            limit=min((k+1)*nq,left+last)
            if(nroll.gt.1) then
              loop=mod((limit-max(k*nq,last)),mroll)
            else
              loop=limit-max(k*nq,last)
            end if
convex, cray, and ibm directives
c$dir no_recurrence
cdir$ ivdep
cibmdir ignore recrdeps
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
c
          limit=min(np,left+last)
          if(nroll.gt.1) then
            loop=mod((limit-max(nratio*nq,last)),mroll)
          else
            loop=limit-max(nratio*nq,last)
          end if
convex, cray, and ibm directives
c$dir no_recurrence
cdir$ ivdep
cibmdir ignore recrdeps
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
c
          last=mod(limit,np)
          left=nran-j
          if(left.gt.0) goto 10
        else
          call errprt(iout,'vranf','incorrect initialization',-1)
        end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine xuinit(x,np,nq,mode,nexec,iseed,init,last,iout)
c
c      initializes a (np,nq) lagged fibonacci generator table
c      with random bits generated by a congruential generator using
c      l'ecuyers decomposition. ref.: bratley et al. p. 214
c
c      x     : vector of length np for initial random number table;
c              output
c      np,nq : parameters p and q of feedback shift register generator;
c              input
c      mode  : operation for lfg:
c              mode=<0 -> subtractive generator, mode=1 additive
c      nexec : number of warm up cycles for the table. nexec*nbit*np
c              random numbers are generated and discarded; input
c      iseed : integer seed for congruential generator generating
c              the bits of the initial table entries; input
c      init  : returns updated seed of congruential generator.
c              0 if table was not initialized, > 0 if ok; output
c      last  : pointer to the last used number in the table; output
c      iout  : unit number for messages, silent for iout.le.0; input
c      subroutines called : none             m. lewerenz mar/93, mar/98
c
      implicit real*8 (a-h,o-z)
      parameter (ia=40692,ib=52774,ic=3791,ip=2147483399)
      parameter (zero=0.d0,half=0.5d0,iphalf=ip/2,nbit=48)
      logical high
      dimension x(np)
c
      if(nq.ge.np.or.iseed.eq.0) then
        call errprt(iout,'xuinit','illegal seed parameter(s)',-1)
      else
c
c      set table to zero and exercise the bit generator a little
c
        ix=iabs(iseed)
        if(ix.ne.0) then
          do i=1,np
            x(i)=zero
            k1=ix/ib
            ix=ia*(ix-k1*ib)-k1*ic
            if(ix.lt.0) ix=ix+ip
          end do
c
c      assemble np numbers with mantissa length nbit from random bits
c      'high' toggle compensates for bias from odd ip
c
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
      end
c
c-----------------------------------------------------------------------
c
      subroutine xuwarm(x,np,nq,mode,nexec,iout)
c
c      warms up (p,q) lagged fibonacci generators (lfg) by nexec rounds
c
c      x     : vector of length np for initial random number table;
c              output
c      np,nq : parameters p and q of feedback shift register generator;
c              input
c      mode  : operation for lfg:
c              mode=<0 -> subtractive generator, mode=1 additive
c      nexec : number of warm up cycles for the table. nexec*nbit*np
c              random numbers are generated and discarded; input
c      iout  : unit number for messages, silent for iout.le.0; input
c      subroutines called : errprt                   m. lewerenz mar/98
c
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0)
      dimension x(np)
c
      if(nq.ge.np.or.np.eq.0.or.nq.eq.0) then
        call errprt(iout,'xuwarm','illegal table parameter(s)',-1)
      else
c
c      exercise the generator for nexec rounds of np prn's
c      separate sections for subtractive or additive version
c
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
      end
C
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
C----------------------------------------------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
C
      FUNCTION R1MACH (IDUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ONE=1.D0,TWO=2.D0,HALF=0.5D0)
      SAVE ICALL,EPS
      DATA ICALL,EPS/0,ONE/
C ---------------------------------------------------------------------
C THIS ROUTINE COMPUTES THE UNIT ROUNDOFF OF THE MACHINE.
C THIS IS DEFINED AS THE SMALLEST POSITIVE MACHINE NUMBER
C U SUCH THAT 1.0 + U .NE. 1.0E0
C ---------------------------------------------------------------------
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
      END
c
c-----------------------------------------------------------------------
c-------------------------- error handling -----------------------------
c-----------------------------------------------------------------------
c
      subroutine errprt(iout,pgname,text,icode)
c
c      prints error messages from library subroutines
c
c      iout   : unit number for message output, 0-> no output; input
c      pgname : name of the subroutine calling errprt; input
c      text   : message text; input
c      icode  : severity code: 0 -> warning, < 0 -> fatal error with
c               abort, else -> error but execution continues
c      subroutines called : strlen                   m. lewerenz dec/93
c
      character pgname*(*),text*(*),header*20,tail*40
      save nerror,nwarn,icall
      common /errcnt/ maxerr,maxwrn
      data icall/0/
c
      if(icall.eq.0) then
        icall=1
        nerror=0
        nwarn=0
      end if
c
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
c
c      write the message on unit iout
c
      if(iout.gt.0) then
        call getstr(pgname,lname,iout)
        call getstr(text,ltext,iout)
        call getstr(header,lhead,iout)
        call getstr(tail,ltail,iout)
        write(iout,'(/6a/)') header(1:lhead),' ',text(1:ltext),' in ',
     #                      pgname(1:lname),tail(1:ltail)
        call flush(iout)
      end if
c
      jcode=icode
      if(maxerr.gt.0.and.nerror.ge.maxerr) then
        if(iout.gt.0) write(iout,'(/a)')
     #  '  *** maximum number of errors exceeded, program stopped *** '
        jcode=-1
      end if
      if(maxwrn.gt.0.and.nwarn.ge.maxwrn) then
        if(iout.gt.0) write(iout,'(/a)')
     #  '  *** maximum number of warnings exceeded, program stopped ***'
        jcode=-1
      end if
      if(iout.gt.0) call flush(iout)
c
      if(jcode.lt.0) stop
      return
c
c ---------------------------------------------------------------------
c      error report, returns current number of errors and warning 
c
      entry errnum(nerr,nwrn)
      nerr=nerror
      nwrn=nwarn
      return
      end
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c----------------------------------------------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine getstr(string,ls,iout)
c
c      eliminates leading and trailing blanks from string and returns
c      length of remaining string
c
c      string : character string; input & output
c               on output leading and trailing blanks of the original
c               string are missing
c      ls     : length of non blank portion of string; output
c      iout   : fortran unit for messages, silent if iout.le.0; input
c      subroutines called: errprt                    m. lewerenz jul/00
c
      character string*(*)
c
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
      end
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c---------------------------- last line --------------------------------
