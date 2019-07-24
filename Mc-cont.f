      PROGRAM SAMC_SZ
      use DynamicResize
      use ifport
      implicit none
      
      include 'mpif.h'
      include 'var.common'

      integer nmax
      parameter (nmax=512)

      real*8 X(3,nmax),AM(nmax),ssmean,SSMEANC,S2,S,R2END,K1,K2
      real*8 SS(nmax),OTEV(3),V(3,nmax),G(nmax),rg(nmax),rand0
      real*8 Phi(nmax),CopPhi(2,nmax)
      real*8 WLFact
      real*8 rndctrl(100),avrnd,sumrndc,sdrnd
      real*8 aECCB,aEB,aNECCB,aNEB,aLoc,aNLoc,aRep,aNRep,aCB,aNCB
      real*8 tlocp,teccbp,tcbp,tebp,trepp,timestart,timefinish
      real*8 SSA,K1A,K2A,OTA(3),SSB,K1B,K2B,OTB(3),ssma,ssmb

      integer(kind=8) MAXST,NSTEP,tryECCBf,tryCBf,cntECCBf,cntCBf
      integer(kind=8) cntEBf,cntRepf,cntLocf

      integer neighb(nmax,nmax),lastneig(nmax),iseed,xkey,nxc
      integer ntot,I,J,NG,NB,II,k,N,NST,nconffile,nres,NHist,EKey
      integer CopType(nmax),StifKey(nmax),NofCop
      integer E1,E1Ctrl,Ecur,E2,E2Ctrl,ErrE1,ErrE2,nDir
      integer MKey,WLRes,MPKey,MPerr,EofKey,MPrnkTrue,DirRnkmn,DirRnkmx
      integer cnfe1(100),cnfe2(100),cnfkey(100),cnfdel(100),y
      integer XShift(6),NSTD,EAA,EAB,EBB
      integer icur,iprev,i2cur,i2prev
      


      character WLResC
      character(Len=256) :: confname,vrmlname,logfname,mtrxfname
      character(Len=256) :: mc_wl,samcf,xconfname,prrunf,parallf
      character(Len=256) :: genm,henm,genmtmp,henmtmp,rep
      character(Len=1024) :: fmt,fmt2,path

      data AM/nmax*1.D0/

      include 'ranmar.common'

      call MPI_Init(MPerr)
      call MPI_Comm_Size(MPI_COMM_WORLD,MPsize,MPerr)
      call MPI_Comm_Rank(MPI_COMM_WORLD,MPrnkTrue,MPerr)
      if (MPrnkTrue.eq.0) write(*,*)'Running ',MPsize,' processes'

      DirRnkmn=0
      DirRnkmx=0
      open(12,FILE='inmc-dir.dat',IOStat=EofKey,Status='old')
       if(EofKey.ne.0) goto 15
       read(12,'(I4)')nDir
       read(12,*)
       do i=1,nDir
        if(i.gt.1) DirRnkmn=DirRnkmx+1
        read(12,'(I9,A28)')DirRnkmx,path

        if(MprnkTrue.le.DirRnkmx) exit
       enddo
      close(12)

      path=adjustl(path)
      EofKey= chdir(path)
C      call getcwd(fmt)
C      write(*,*)Mprnk,EofKey,fmt
C      stop''
15    Mprnk=MPrnkTrue-DirRnkmn
      open(12,FILE='inmc-input.dat')
      READ(12,'(I15/I15/I15/I15/I15/F15.1/I15/F15.6/I15/I15/I15)')
     *KEY,N,NC,NST,MAXST,AX,iseed,rand0,NG,NofCop,LenA
      READ(12,'(I15/I15/F15.5/F15.3/F15.3/F15.3/I15/I15)')LenB,
     *SMT0,SMgam0,AngStifA,DevStifA,AngStifB,EWLbound,Wkey
      read(12,'(I15/I15)')EPSST,xkey
      close(12,STATUS='KEEP')

      write(logfname,'(A7,I0.2,A4)')'mc-log_',MPrnk,'.dat'
      write(xconfname,'(A5,I0.2,A4)')'x_cnf',MPrnk,'.dat'
      write(prrunf,'(A11,I0.2,A4)')'mc-prodrun_',MPrnk,'.dat'
      write(parallf,'(A15)')'inmc-parall.dat'
      write(genm,'(A7,I0.2,A4)')'mc-gnm_',MPrnk,'.dat'
      write(henm,'(A7,I0.2,A4)')'mc-hnm_',MPrnk,'.dat'
      write(genmtmp,'(A10,I0.2,A4)')'mc-gnmtmp_',MPrnk,'.dat'
      write(henmtmp,'(A10,I0.2,A4)')'mc-hnmtmp_',MPrnk,'.dat'

      distr=0.D0
      NB=N/NC
      SIG=1.D0
      EPS=1
      BOND=1.D0
      RCUT=BOND+0.5
      RMAX=2.5D0
      dr=0.1D0
      BndMax=1.25*BOND
      BndMin=0.8*BOND
      BndMax3=BndMax**3
      BndMin3=BndMin**3
      Mxctet=(2.D0*BndMax**2-SIG**2)/(2.D0*BndMax**2)
C      write(*,*)Mxctet
      RTop=RCUT+dr*sqrt(3.D0)
      RTop=Rtop*Rtop
      Rbot=RCUT-dr*sqrt(3.D0)
      Rbot=Rbot*Rbot
      Rbot2=SIG+dr*sqrt(3.D0)
      Rbot2=Rbot2*Rbot2
      NFOR=10
      AY=AX
      AZ=AX
      SAMCPST=0
      WLFact=2.D0
      Sumgam=0.D0
      rndctrl=0.D0
      CDeepKey=0
      cntLoc=0
      cntECCB=0
      cntEB=0
      cntRep=0
      cntCB=0
      cntNLoc=0
      cntNECCB=0
      cntNEB=0
      cntNRep=0
      ECCBf=0
      cntNCB=0
      tryECCB=0
      tryCB=0
      tryECCBf=0
      tryCBf=0
      cntECCBf=0
      cntCBf=0
      pcb=1.d0
      peccb=1.d0
      cntEBf=0
      cntRepf=0
      cntLocf=0
      tlocp=0.d0
      teccbp=0.d0
      tcbp=0.d0
      tebp=0.d0
      trepp=0.d0
      XShift=0
      ngwinbnd=0

      if(xkey.eq.1) then
       call OutPrintCurCnf(N,X,E1,E2,CopType,0)
      endif

      if(KEY.GE.4) then
       nconffile=1
       nres=1
       write(confname,'(A5,I0.2,A1,I0.3,A4)')'xcnf_',MPrnk,'_'
     *  ,nconffile,'.dat'
       open (34,FILE=confname)
      endif

      iseed=iseed+MPrnk*10
C      Initialize random number generator
      call rmarin(iseed)


C      if(MPsize.ne.1) then
      MPKey=0
      call IniParall(N,NB,X,rand0,CopType,StifKey,CopPhi,
     *      E1,E2,WLFact,neighb,lastneig,MPKey,parallf)
      if(MPKey.eq.1) goto 605
C      endif
C      if(MPsize.eq.1) then
C       call BEGIN(N,NB,X,rand0,CopType,Phi,StifKey,CopPhi,MPrnk)
C       WinKey=2
C       nWin=1
C      endif
C      goto 610


      call E1Full(N,X,E1,1,CopType)
      call E2Full(N,X,E2,Phi,CopPhi,CopType,StifKey,1)

      Ecur=E1+E2

      if(KEY.eq.0) then
       nWLHf=1
       call SAMC(E1,0,E2,0,-1,0.D0,WLRes,0.D0,0,0,0,N)
      else
       call SAMC(E1,0,E2,0,-1,0.D0,WLRes,WLFact,0,0,0,N)
C       if(KEY.eq.1)
       if(KEY.eq.2) then
        SumGam=0.d0
        nWLHf=0
        WLH=int(Mask)
        SAMCPST=0
       endif
       if(KEY.eq.3) then
        stop 'Warning! Do not use KEY=3 mode!!!'
       endif
       if(KEY.eq.4) then
        allocate(PrWLH(Nbin,N2bin))
        allocate(oWLH(Nbin,N2bin))
        WLH=int(mask)
        PrWLH=0
        oWLH=0
        do i=1,Nbin
         EAA=i-Nbin
         if(EAA.gt.-700) then
          if(MOD(i,50).ne.0)cycle
         endif
         if(EAA.gt.-1100.and.EAA.lt.-700) then
          if(MOD(i,20).ne.0)cycle
         endif
         if(EAA.gt.-1300.and.EAA.lt.-1100) then
          if(MOD(i,10).ne.0)cycle
         endif
         do k=1,N2bin
          if(k.gt.20) then
           if(MOD(k,10).ne.0) cycle
          endif
          oWLH(i,k)=5
         enddo
        enddo
        SAMCPST=0
       endif
      endif

      call OutIni(MPrnkTrue,N,CopType,StifKey,E1,E2,logfname,prrunf,
     *            genmtmp,henmtmp)
C     MAIN CYCLE
C     ==================================================================
      NSTEP=1
      MKey=0
      E1cmax=E1
      E1cmin=E1
      E2cmax=E2
      E2cmin=E2
      Efcmax=E1cmax+E2cmax
      Efcmin=E1cmin+E2cmin
      cntLoc=0
      cntECCB=0
      cntEB=0
      cntRep=0
      cntCB=0
      cntNLoc=0
      cntNECCB=0
      cntNEB=0
      cntNRep=0
      cntECCBf=0
      cntNCB=0
      tryECCB=0
      tryCB=0
      NHist=10
      cnfdel=0
      icur=0
      iprev=0
      i2cur=0
      i2prev=0
      tloc=0.D0
      teb=0.D0
      teccb=0.D0
      trept=0.D0
      call cpu_time(timestart)
      if(KEY.ge.4) open(36,FILE=prrunf)

      do while (MKey.EQ.0)
       call Step(N,X,V,NSTEP,neighb,lastneig,CopType,E1,WLFact,Phi,
     *          CopPhi,StifKey,E2,rndctrl)
       Ecur=E1+E2

       if(xkey.eq.1)then
        call OutPrintCurCnf(N,X,E1,E2,CopType,1)
       endif

       if(MOD(NSTEP,1000000).eq.0) then
        call CheckShift(N,X,V,XShift)
C        call WHLP(N,X,V,ENV,EStif,Phi,CopPhi,CopType,StifKey,neighb,
C     *            lastneig,NSTEP,WLFact)
C        call Exchange(N,X,V,ENV,EStif,Phi,CopPhi,CopType,StifKey,neighb,
C     *            lastneig,NSTEP,WLFact)
       endif


       if(MOD(NSTEP,10*NST).EQ.0) then
        call SAMC(0,0,0,0,-2,0.d0,WLRes,0.d0,0,0,0,N)
        call OutPrintGHTMP(henmtmp,genmtmp,WLRes,NSTEP,WLFact)
       endif

       if(CDeepKey.eq.2.and.E1.le.EDeepBnd-5)then
        CDeepKey=0
        EofKey=0
        write(rep,'(A6,I0.2,A4)')'x_cnfi',MPrnk+1,'.dat'
        open(50,FILE=rep,Status='old',IOStat=EofKey)
        if(EofKey.ne.0) then
         close(50)
         open(50,FILE=rep)
         do i=1,N
          write(50,'(3E18.10)')X(1,i),X(2,i),X(3,i)
         enddo
         close(50)
        endif
       endif


       if(MOD(NSTEP,NST).EQ.0)then
        WLResC='n'
        SSMEAN=0.
        do I=1,NC
         J=1+NB*(I-1)
         SS(I)=S2(NB,X(1,J),AM(J))
         SSMEAN=SSMEAN+SS(I)
        enddo
        SSMEAN=SSMEAN/NC
        call SAMC(0,0,0,0,-2,0.d0,WLRes,0.d0,0,0,0,N)
        if(WLRes.eq.2) WLResC='y'
       if(KEY.eq.3.and.MOD(NSTEP,SMT0).eq.0)call Umbrella

        call E1Full(N,X,E1Ctrl,1,CopType)
        call E2Full(N,X,E2Ctrl,Phi,CopPhi,CopType,StifKey,1)

        ErrE1=E1-E1Ctrl
        ErrE2=E2-E2ctrl
        if(ErrE1.ne.0)E1=E1Ctrl
        aLoc=dble(cntLoc)/dble(NST*N)*100.D0
        aECCB=dble(cntECCB)/dble(tryECCB)*100.D0
        aCB=dble(cntCB)/dble(tryCB)*100.D0
        aEB=dble(cntEB)/dble(2*NST)*100.D0
C        if(MPrnk.eq.0)write(*,*)cntEB,NST,aEB,cntREP,aRep
        aRep=dble(cntRep)/dble(2*NST)*100.D0
        tryECCBf=tryECCBf+tryECCB
        tryCBf=tryCBf+tryCB
        cntECCBf=cntECCBf+cntECCB
        cntCBf=cntCBf+cntCB
        cntEBf=cntEBf+cntEB
        cntRepf=cntRepf+cntRep
        cntLocf=cntLocf+cntLoc
        cntEB=0
        cntRep=0
        cntLoc=0
        cntCB=0
        cntECCB=0
        tryECCB=0
        tryCB=0
        if(teccb-teccbp.gt.0.00001) then
         peccb=(aECCB/(teccb-teccbp))/(aLoc/(tloc-tlocp))
         if(peccb.gt.1.d0)peccb=1.d0
         if(peccb.lt.0.001)peccb=0.001
        else
         peccb=1.d0
        endif
        if(tcb-tcbp.gt.0.00001) then
         pcb=(aCB/(tcb-tcbp))/(aLoc/(tloc-tlocp))
         if(pcb.gt.1.d0)pcb=1.d0
         if(pcb.lt.0.001)pcb=0.001
        else
         pcb=1.d0
        endif
        tlocp=tloc
        teccbp=teccb
        tcbp=tcb
        tebp=teb
        trepp=trept

       write(*,'(I3,I15,I4,2I6,I4,4I6,A4,10F8.2,2F6.3)')MPrnkTrue,NSTEP,
     *   ErrE1,E1cmax,E1cmin,ErrE2,E2cmax,E2cmin,Efcmax,Efcmin,
     *   WLResC,aLoc,tloc/60.D0,aECCB,teccb/60.D0,aCB,tcb/60.D0,aEB,
     *   teb/60.D0,aRep,trept/60.D0,peccb,pcb
     
        open(31,FILE=logfname,position="append")
        write(31,'(I15,E10.3,I4,2I8,I4,4I8,A4,10F8.2,2F6.3)')NSTEP,
     *   SSMEAN,ErrE1,E1cmax,E1cmin,ErrE2,E2cmax,E2cmin,Efcmax,
     *   Efcmin,WLResC,aLoc,tloc/60.D0,aECCB,teccb/60.D0,aCB,
     *   tcb/60.D0,aEB,teb/60.D0,aRep,trept/60.D0,peccb,pcb
        close(31)
       endif

       if(KEY.GE.4.and.MOD(NSTEP,1000).eq.0)then
        call RTEND(X, N, NC, NB, R2END)
        call tensor(N,NC,X,CopType,SSMEANC,K1,K2,OTEV,SSA,K1A,K2A,OTA,
     *    SSB,K1B,K2B,OTB)
        call EnergyContr(N,V,CopType,neighb,lastneig,EAA,EAB,EBB)
        SSMEAN=0.
        do I=1,NC
         J=1+NB*(I-1)
         SS(I)=S2(NB,X(1,J),AM(J))
         SSMEAN=SSMEAN+SS(I)
        enddo
        SSMEAN=SSMEAN/NC

        write(36,'(I15,3I5,20E13.5,3I5)')NSTEP,Ecur,E1,E2,SSMEAN,
     *    R2END,SSMEANC,K1,K2,OTEV(1),OTEV(2),OTEV(3),SSA,K1A,K2A,
     *    OTA(1),OTA(2),OTA(3),SSB,K1B,K2B,OTB(1),OTB(2),OTB(3),
     *    EAA,EAB,EBB
C       write(*,'(2I5,6E13.5)')ENV,EStif,SSMean,SSMEANC,ssa,ssma,ssb,ssmb
C           call RDF(N,NG,AX,AY,AZ,X,G,RG,0)
C        if(MOD(NSTEP,5000000).eq.0)then
         icur=E1-Emin+1
         i2cur=E2-E2min+1
         if(PrWLH(icur,i2cur).lt.oWLH(icur,i2cur)) then
          if(icur.ne.iprev.or.i2cur.ne.i2prev) then
          nres=nres+1
          if (MOD(nres,1000).eq.0) then
           close(34)
           nconffile=nconffile+1
           write(confname,'(A5,I0.2,A1,I0.3,A4)')'xcnf_',MPrnk,'_',
     *     nconffile,'.dat'
           open(34,FILE=confname)
          endif
         
          Write(34,*)'STEP: ',NSTEP,E1,E2
          do i=1,N
           Write(34,'(3E14.6)')X(1,I),X(2,I),X(3,I)
          enddo
          endif
         endif
         PrWLH(icur,i2cur)=PrWLH(icur,i2cur)+1
         iprev=icur
         i2prev=i2cur
       endif

       NSTEP=NSTEP+1
       if (NSTEP.GT.MAXST) Mkey=1
30    end do
      call cpu_time(timefinish)

C     END OF MAIN CYCLE
C===================================================================
C      close(580)



      NSTEP=NSTEP-1

      WLResC='n'
      call SAMC(0,0,0,0,-2,0.d0,WLRes,0.d0,0,0,0,N)
      if(WLRes.eq.2) WLResC='y'

C      if(KEY.GE.2) then
C       call RDF(N,NG,AX,AY,AZ,X,G,RG,1)
C      endif

      open(10,FILE=xconfname)
C     Saving of current values of parameters and coordinates
600   do I=1,N
         write(10,'(3E18.10)')X(1,I),X(2,I),X(3,I)
      enddo
      close(10)

      open(31,FILE=logfname,position="append")

      sumrndc=sum(rndctrl)
      rndctrl=rndctrl/sumrndc
      avrnd=sum(rndctrl)/100.D0
      sdrnd=dsqrt(sum((rndctrl-avrnd)**2))/100.D0
      if(sdrnd.ge.0.00001)then
       write(vrmlname,'(A4,I0.2,A4)')'rnd_',MPrnk,'.dat'
       open(11,FILE=vrmlname)
       do i=1,100
        write(11,'(F10.7)')rndctrl(i)
       enddo
       close(11)
      endif
      write(31,'(A20,F8.5,A8,F8.5)')'Random generator mu ',avrnd,
     * ' st.dev ',sdrnd


      write(31,'(A6,2(A15,A7),A10)')'','Accepted steps','%','E1<>E2',
     * '%','Time,m'
      aLoc=dble(cntLocf)/dble(NSTEP)/dble(N)*100.D0
      aNLoc=dble(cntNLoc)/dble(NSTEP)/dble(N)*100.D0
      aECCB=dble(cntECCBf)/dble(NSTEP)*100.D0
      aEB=dble(cntEBf)/dble(NSTEP)*100.D0
      aCB=dble(cntCBf)/dble(NSTEP)*100.D0
      aNECCB=dble(cntNECCB)/dble(NSTEP)*100.D0
      aNEB=dble(cntNEB)/dble(NSTEP)*100.D0
      aNCB=dble(cntNCB)/dble(NSTEP)*100.D0
      aRep=dble(cntRepf)/dble(NSTEP)*100.D0
      aNRep=dble(cntNRep)/dble(NSTEP)*100.D0
      Write(31,'(A6,2(I15,F7.2),F10.2)')'Local',cntLocf,aLoc,cntNLoc,
     *aNLoc,tloc/60.d0
      Write(31,'(A6,2(I15,F7.2),F10.2)')'ECCB',cntECCBf,aECCB,cntNECCB,
     *aNECCB,teccb/60.d0
      Write(31,'(A6,2(I15,F7.2),F10.2)')'EB',cntEBf,aEB,cntNEB,aNEB,
     *teb/60.d0
      Write(31,'(A6,2(I15,F7.2),F10.2)')'Rept',cntRepf,aRep,cntNRep,
     * aNRep,trept/60.d0
      Write(31,'(A6,2(I15,F7.2),F10.2)')'CB',cntCBf,aCB,cntNCB,aNCB,
     * tcb/60.d0

      Write(*,*)'Process ',MPrnkTrue,' CPU Time: ',timefinish-timestart
      Write(31,*)'CPU Time: ',timefinish-timestart
C      Write(*,*)'CPU Time - Stepe: ',tstep
C      Write(*,*)'CPU Time - Stepi: ',tall
C      Write(*,*)'CPU Time - Table: ',ttab
C      Write(*,*)'CPU Time - ENV  : ',tenv
C      Write(*,*)'CPU Time - Ebond: ',teb
C      Write(*,*)'CPU Time - Estif: ',test
C      Write(*,*)'CPU Time - SAMC : ',tSMC
C      Write(*,'(A8,I2,A11,F8.2)')' Process ',MPrnk,' Sum gamma: ',SumGam
      
C      Write(31,*)'CPU Time - Step : ',tstep
C      Write(31,*)'CPU Time - Stepi: ',tall
C      Write(31,*)'CPU Time - Table: ',ttab
C      Write(31,*)'CPU Time - ENV  : ',tenv
C      Write(31,*)'CPU Time - Ebond: ',teb
C      Write(31,*)'CPU Time - Estif: ',test
C      Write(31,*)'CPU Time - SAMC : ',tSMC
      Write(31,*)'Sum of gamma SAMC: ',SumGam
      if(ECCBf.ne.0)Write(*,*)'ECCB fails',ECCBf
      if(ECCBf.ne.0)Write(31,*)'ECCB fails',ECCBf
      Write(31,*)'Coordinate shifts made'
      Write(31,'(6I10)')(XShift(i),i=1,6)

      if(KEY.eq.4)then
       write(genm,'(A9,I0.2,A4)')'mc-prodg_',MPrnk,'.dat'
       write(henm,'(A9,I0.2,A4)')'mc-prodh_',MPrnk,'.dat'
      endif
      
      open(40,FILE=genm)
      open(41,FILE=henm)
      write(40,'(I15,A18,A4)')NSTEP+SAMCPST,' steps. G(Env,Est)',WLResC
      write(41,'(I15,A18,A4)')NSTEP+SAMCPST,' steps. H(Env,Est)',WLResC
      write(40,'(6I8,E15.6,I8,E15.8)')Nbin,Emin,Emax,N2bin,E2min,E2max,
     * WLFact,nWLHf,SumGam
      write(41,'(6I8,E15.6,I8,E15.8)')Nbin,Emin,Emax,N2bin,E2min,E2max,
     * WLFact,nWLHf,SumGam

      write(fmt,'(A4,I0,A4)')'(A8,',N2bin,'I18)'
      write(40,fmt)' ENV\Est',(E2min+k-1,k=1,N2bin)
      write(41,fmt)' ENV\Est',(E2min+k-1,k=1,N2bin)
      write(fmt2,'(A4,I0,A4)')'(I8,',N2bin,'I18)'
      write(fmt,'(A4,I0,A6)')'(I8,',N2bin,'E18.8)'

      do i=1,Nbin
       write(40,fmt)Emin+i-1,(WLG(i,k),k=1,N2bin)
       write(41,fmt2)Emin+i-1,(WLH(i,k),k=1,N2bin)
      enddo
      close(40)
      close(41)
      close(31)
      close(34)
      if (KEY.ge.4) then
       close(36)
       if(allocated(PrWLH))deallocate(PrWLH)
       if(allocated(oWLH))deallocate(oWLH)
      endif
      
      write(vrmlname,'(A11)')'mc-VRML.wrl'
      call ToVRML(N,NB,X,SIG,CopType,vrmlname)
      if(allocated(WLH)) deallocate(WLH)
      if(allocated(WLG)) deallocate(WLG)
      if(allocated(PW)) deallocate(PW)
C      if(allocated(stan)) then
C       call StepAnalysisOut
C      endif
C      call Distr_BONDS(N,X,int(Estif+ENV),1)
C      call Distr_NBHD(N,X,int(Estif+ENV),1)

605   if(MPsize.gt.1) then
       call MPI_Barrier(MPI_COMM_WORLD,MPerr)
C       call FinParall
      endif
      call MPI_Finalize(MPerr)
C610   call fullstatZ(N,NC,X,AX)
C      close(120)
      if(MPrnkTrue.eq.0) stop 'Successfull end '
      if(MPrnkTrue.ne.0) stop

505   end
Cedum
