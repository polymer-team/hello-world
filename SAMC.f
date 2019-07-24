      subroutine SAMC(Eold,Enew,E2old,E2new,WLKey,rnd,WLResult,SMgam,
     * AccSt,AccNst,NSTEP,N)
      use DynamicResize
      implicit none
      include 'var.common'

      integer(kind=8)HSUM,HMIN,HMAX,AVGH,VARH,avh,NSTEP,AccSt,AccNSt
      integer iold,inew,i,WLKey,j1,j2,WLResult,k,i2new,i2old,m,j
      integer Enew,Eold,E2old,E2new,N,NgKey
      real rnd
      real*8 VARHr,WLPROB,SMest,hold,SMgam2,Nstgap,Maxsh
      real*8 SMgam,SMshift,SMshift1,SMshift2,delt,GMIN

      character(Len=256) :: genm,henm,newstate
      character(Len=1024) :: fmt,fmt2

      Data Maxsh/0.D0/
      WLResult=0

      if (WLKey.Eq.0.and.KEY.lt.3) then
       iold=Eold-Emin+1
       inew=Enew-Emin+1
       i2old=E2old-E2min+1
       i2new=E2new-E2min+1
       if(i2new.lt.1.or.i2new.gt.N2bin.or.inew.lt.1.or.inew.gt.Nbin)then
        if(Enew.lt.Eminbnd) then
         if(WinKey.eq.0.or.WinKey.eq.1) then
          WLResult=0
          return
         endif
        endif
        if(Enew.gt.Emaxbnd) then
         if(WinKey.eq.0.or.WinKey.eq.-1) then
          WLResult=0
          return
         endif
        endif
        if(E2new.lt.E2minbnd) then
C         write(*,*)'E2min'
         if(WinKey.ne.2) then
          WLResult=0
          return
         endif
        endif
        if(E2new.gt.E2maxbnd) then
C         write(*,*)'E2max',E2new,E2maxbnd
         if(WinKey.ne.2) then
          WLResult=0
          return
         endif
        endif

        Nbinnew=Nbin
        N2binnew=N2bin
        j1=1
        j2=1
        if(i2new.lt.1) then
         N2binnew = N2bin-i2new+1
         j2=N2binnew-N2bin+1
         i2new=1
         i2old=i2old+j2-1
         E2min=E2new
        endif
        if(i2new.gt.N2bin) then
         N2binnew = i2new
         E2max=E2new
        endif
        if(inew.lt.1) then
         Nbinnew = Nbin-inew+1
         j1=Nbinnew-Nbin+1
         inew=1
         iold=iold+j1-1
         Emin=Enew
        endif
        if(inew.gt.Nbin) then
         Nbinnew = inew
         Emax=Enew
        endif
        allocate(WLGtmp(Nbin,N2bin))
        allocate(WLHtmp(Nbin,N2bin))
        allocate(PWtmp(Nbin,N2bin))
        allocate(Masktmp(Nbin,N2bin))
        
C        allocate(gemodtmp(6,Nbin,N2bin))
C        gemodtmp=gemod
C        deallocate(gemod)
C        allocate(gemod(6,Nbinnew,N2binnew))
C        gemod=0.D0
C     gemod(1:6,j1:Nbin+j1-1,j2:N2bin+j2-1)=gemodtmp(1:6,1:Nbin,1:N2bin)
C        deallocate(gemodtmp)


        WLGtmp=WLG
        WLHtmp=WLH
        PWtmp=PW
        MaskTmp=Mask
        deallocate(WLH)
        deallocate(WLG)
        deallocate(PW)
        deallocate(Mask)
        allocate(WLH(Nbinnew,N2binnew))
        allocate(WLG(Nbinnew,N2binnew))
        allocate(PW(Nbinnew,N2binnew))
        allocate(Mask(Nbinnew,N2binnew))

        WLH=0
        WLH(j1:Nbin+j1-1,j2:N2bin+j2-1)=WLHtmp(1:Nbin,1:N2bin)
        Mask=0.d0
        Mask(j1:Nbin+j1-1,j2:N2bin+j2-1)=masktmp(1:Nbin,1:N2bin)

        WLG=0.D0
        WLG(j1:Nbin+j1-1,j2:N2bin+j2-1)=WLGtmp(1:Nbin,1:N2bin)

        PW=0.D0
        PW(j1:Nbin+j1-1,j2:N2bin+j2-1)=PWtmp(1:Nbin,1:N2bin)

        N2bin=N2binnew
        Nbin=Nbinnew

        deallocate(WLGtmp)
        deallocate(WLHtmp)
        deallocate(PWtmp)
        deallocate(Masktmp)
       endif

       if(Mask(inew,i2new).lt.0.5) then
        allocate(PWtmp(Nbin,n2bin))
        PWtmp=PW
        Mask(inew,i2new)=1.D0
        call Weights(Nbin,N2bin)
        deallocate(PWtmp)
        if(Smgam.lt.0.00001) then
         if(samcKey.ne.6) then
          hsum=sum(WLH)
          varhr=sum(Mask)
          AVGH=HSUM/int(varhr)
          SMshift=50.D0
          if(WLG(iold,i2old).lt.50.D0) then
           SMshift=WLG(iold,i2old)
          endif
          varhr=dble(WLH(iold,i2old))/dble(AVGH)
          if(varhr.gt.1.D0)varhr=1.D0
          SMshift=SMshift*varhr
          WLG(inew,i2new)=WLG(iold,i2old)-SMshift
C          write(*,*)inew-iold,i2new-i2old
         else
          SMshift1=0.D0
          Smshift2=0.D0
          j2=0
          j1=0
          hsum=sum(WLH(inew,1:N2bin))
          if(hsum.ne.0) then
           NgKey=0
           i=i2new+1
           do while(NgKey.eq.0)
            if(i.ge.N2bin+1) exit
            if(WLH(inew,i).ne.0) then
             NgKey=1
             SMShift2=50.D0
             j2=i
             if(i+1.le.N2bin) then
              delt=WLG(inew,i+1)-WLG(inew,i)
              SMshift2=delt*(i-i2new)
             endif
            endif
            i=i+1
           enddo
          endif
          hsum=sum(WLH(1:Nbin,i2new))
          if(hsum.ne.0) then
           NgKey=0
           i=inew-1
           do while(NgKey.eq.0)
            if(i.le.0) exit
            if(WLH(inew,i).ne.0) then
             NgKey=1
             SMShift1=50.D0
             j1=i
             if(i-1.ge.1) then
              delt=WLG(inew,i-1)-WLG(inew,i)
              SMshift1=delt*(inew-i)
             endif
            endif
            i=i-1
           enddo
          endif
          if(j1*j2.eq.0) then
           if(j1.eq.0) then
            WLG(inew,i2new)=WLG(inew,j2)-SMshift2
           endif
           if(j2.eq.0) then
            WLG(inew,i2new)=WLG(j1,i2new)-SMshift1
           endif
          else
           if(inew-j1.lt.j2-inew) then
            WLG(inew,i2new)=WLG(j1,i2new)-SMshift1
           else
            WLG(inew,i2new)=WLG(inew,j2)-SMshift2
           endif
          endif
         endif
        endif
        if(IniKey.eq.0) then
         write(newstate,'(A7,I0.2,A4)')'mc-nst_',MPrnk,'.dat'
         open(70,FILE=newstate,position="append")
         write(70,'(4I6,I3,I15,3E18.8)')Enew,E2new,Eold,E2old,samcKey,
     *   NSTEP,SMgam,WLG(inew,i2new),WLG(iold,i2old)
         close(70)
        endif
        
        if (E1cmax.lt.Enew) E1cmax=Enew
        if (E1cmin.gt.Enew) E1cmin=Enew
        if (E2cmax.lt.E2new) E2cmax=E2new
        if (E2cmin.gt.E2new) E2cmin=E2new
        if (Efcmax.lt.Enew+E2new) Efcmax=Enew+E2new
        if (Efcmin.gt.Enew+E2new) Efcmin=Enew+E2new
       endif

C       if(iold.NE.inew.or.i2old.ne.i2new) NewBinRate=NewBinRate+1.D0

       if(WLG(iold,i2old).gt.WLG(inew,i2new)) then
        WLPROB=1.D0
       else
        WLPROB=dexp(WLG(iold,i2old)-WLG(inew,i2new))
       endif
       SumGam=SumGam+SMgam
C       gemod(SAMCKey,1:Nbin,1:N2bin)=gemod(SAMCKey,1:Nbin,1:N2bin)-WLG
C       WLG=WLG-SMgam*PW
       if (WLPROB.LE.rnd) then
C         no transition
        WLH(iold,i2old)=WLH(iold,i2old)+1
        WLG(iold,i2old)=WLG(iold,i2old)+SMgam
        WLResult=0
       else
C         transition to new BIN
        WLH(inew,i2new)=WLH(inew,i2new)+1
        WLG(inew,i2new)=WLG(inew,i2new)+SMgam
        WLResult=1
C        AcRate=AcRate+1.D0
        AccSt=AccSt+1
        if(iold.NE.inew.or.i2old.ne.i2new) AccNst=AccNst+1
       endif
C       gemod(SAMCKey,1:Nbin,1:N2bin)=gemod(SAMCKey,1:Nbin,1:N2bin)+WLG
       return
      endif


      if (WLKey.GT.0.and.KEY.lt.3) then
       iold=Eold-Emin+1
       i2old=E2old-E2min+1
       SumGam=SumGam+SMgam
C       gemod(SAMCKey,1:Nbin,1:N2bin)=gemod(SAMCKey,1:Nbin,1:N2bin)-WLG
C       WLG=WLG-SMgam*PW
       WLH(iold,i2old)=WLH(iold,i2old)+1
       WLG(iold,i2old)=WLG(iold,i2old)+SMgam
       WLResult=0
C       gemod(SAMCKey,1:Nbin,1:N2bin)=gemod(SAMCKey,1:Nbin,1:N2bin)+WLG
       return
      endif

C checking for flatness of H(E)
      if (WLKey.eq.-2) then
C        Write(*,*)'MAXSHIFT  ',Maxsh
        WLResult=0
C        HSUM=0
        hsum=sum(WLH)
        varhr=sum(Mask)
        HMIN=WLH(1,1)
        GMIN=10000.D0
C        HMAX=WLH(1,1)
        hmax=maxval(WLH)
        do 20 i=1,Nbin
        do 20 k=1,N2bin
C           HSUM=HSUM+WLH(i,1)
         if(WLH(i,k).ne.0) then
          if(HMIN.eq.0)HMIN=WLH(i,k)
          if(WLH(i,k).LT.HMIN) HMIN=WLH(i,k)
C          if(GMIN.lt.1.D-6)GMIN=WLG(i,k)
          if(WLG(i,k).LT.GMIN) GMIN=WLG(i,k)
         endif
20      continue
        if(GMIN.gt.500.D0) then
         WLG=WLG-(GMIN+100.D0)*Mask
        endif
C        AVGH=INT(HSUM/Nbin)
        AVGH=HSUM/int(varhr)
        VARH=HMAX-AVGH
        VARHr=REAL(VARH)/REAL(AVGH)
        if (VARHr.GT.0.2) then
         return
        endif
        VARH=AVGH-HMIN
        VARHr=REAL(VARH)/REAL(AVGH)
        if (VARHr.GT.0.2) then
         return
        endif
        WLResult=2
        return
      endif

      if (WLkey.eq.0.and.KEY.ge.3) then
       iold=Eold-Emin+1
       inew=Enew-Emin+1
       i2old=E2old-E2min+1
       i2new=E2new-E2min+1

       if (inew.lt.1.or.inew.gt.Nbin) then
        WLResult=0
        return
       endif
       if (i2new.lt.1.or.i2new.gt.N2bin) then
        WLResult=0
        return
       endif
       if (Mask(inew,i2new).lt.0.5) then
        WLResult=0
        return
       endif

C      Mask! Esli popal kuda ne popadali!!! Return!
C       if(iold.NE.inew.or.i2old.NE.i2new) NewBinRate=NewBinRate+1.D0


       if(WLG(iold,i2old).gt.WLG(inew,i2new)) then
        WLPROB=1.D0
       else
        WLPROB=dexp(WLG(iold,i2old)-WLG(inew,i2new))

       endif

       if (WLPROB.LE.rnd) then
C         no transition
        WLH(iold,i2old)=WLH(iold,i2old)+1
        WLResult=0
       else
C         transition to new BIN
        WLH(inew,i2new)=WLH(inew,i2new)+1
        WLResult=1
C        AcRate=AcRate+1.D0
        if(iold.NE.inew.or.i2old.NE.i2new) AccNst=AccNst+1
        AccSt=AccSt+1
       endif
C       write(*,'(E15.8,F10.6I2,2I4)')WLPROB,rnd,WLResult,int(Eold),
C     *  int(Enew)
       return
      endif

      if (WLKey.GT.0.and.KEY.ge.3) then
       iold=Eold-Emin+1
       i2old=E2old-E2min+1

       WLH(iold,i2old)=WLH(iold,i2old)+1
       WLResult=0
       return
      endif


      if (WLKey.EQ.-1) then
       if(KEY.eq.0) then
       select case(Winkey)
        case(0)
         Nbin=2
         if(Eold.lt.Emaxbnd) then
          Emin=Eold
          Emax=Eold+1.D0
         else
          Emax=Eold
          Emin=Eold-1.D0
         endif
        case(1)
         Nbin=2
         if(Eold.gt.Eminbnd) then
          Emax=Eold
          Emin=Eold-1.D0
         else
          Emin=Eold
          Emax=Eold+1.D0
         endif
C         Nbin=Eold-Eminbnd+1
C         Emin=Eminbnd
C         Emax=Eold
        case(-1)
         Nbin=2
         if(Eold.lt.Emaxbnd) then
          Emin=Eold
          Emax=Eold+1.D0
         else
          Emax=Eold
          Emin=Eold-1.D0
         endif
C         Nbin=Emaxbnd-Eold+1
C         Emin=Eold
C         Emax=Emaxbnd
        case(2)
        if (Eold.le.0) then
         Nbin=2
         Emax=Eold
         Emin=Eold-1.D0
        else
         Nbin=2
         Emax=Eold+1.D0
         Emin=Eold
        endif

        case default
        write(*,*)'SAMC: Error, wrong ident. of window'
       end select

       N2bin=1
       E2max=E2old
       E2min=E2old
       Nbinnew=Nbin
       allocate(WLH(Nbin,N2bin))
       allocate(WLG(Nbin,N2bin))
       allocate(PW(Nbin,N2bin))
       allocate(mask(Nbin,N2bin))
       
C       allocate(gemod(6,Nbin,N2bin))
C       gemod=0.d0

C       iold=Eold-Emin+1
C       i2old=E2old-E2min+1

       mask=1.D0
       WLG=0.D0
       WLH=0
       SumGam=0.D0

       call Weights(Nbin,N2bin)


       Return
       else
        write(genm,'(A7,I0.2,A4)')'mc-gnm_',MPrnk,'.dat'
        write(henm,'(A7,I0.2,A4)')'mc-hnm_',MPrnk,'.dat'
        open(40,FILE=genm)
        open(41,FILE=henm)
        read(40,'(I15)')SAMCPST
        read(41,*)
        read(41,*)
        read(41,*)

        read(40,'(6I8,E15.6,I8,E15.8)')Nbin,Emin,Emax,N2bin,E2min,E2max,
     *       SMgam,nWLHf,SumGam
        read(40,*)

        write(fmt,'(A4,I0,A4)')'(I8,',N2bin,'I18)'
        write(fmt2,'(A4,I0,A6)')'(I8,',N2bin,'E18.8)'

        allocate(WLH(Nbin,N2bin))
        allocate(WLG(Nbin,N2bin))
        allocate(PW(Nbin,N2bin))
        allocate(Mask(Nbin,N2bin))

C       allocate(gemod(Nbin,N2bin))
C       gemod=0.D0

        Mask=0.D0

        do 47 i=1,Nbin
         read(40,fmt2)j,(WLG(i,k),k=1,N2bin)
         read(41,fmt) j,(WLH(i,k),k=1,N2bin)
         do 47 k=1,N2bin
          if(WLH(i,k).ne.0) Mask(i,k)=1.d0
47      continue
        call Weights(Nbin,N2bin)

        close(40)
        close(41)
       endif
      endif
      end
      
      subroutine Weights(L1,L2)
      use DynamicResize
      implicit none
      include 'var.common'
      integer K,i,j,L1,L2
      real*8 Denom,Sumc,Acc,PWL,Den2(L2)

      Denom=0.D0
      Sumc=0.D0
      Den2=0.d0
C      K=L
      K=0
C      do 50 i=1,L1
C      do 50 j=1,L2
C       if(WLH(i,j).ne.0)K=K+1
C50    continue
C      if(incKey.eq.1)K=K+1
C      if(incKey.eq.0)K=L1

C      if(i1new.eq.L1.or.i2new.eq.L2)K=K+1
C      if(i1new.lt.1.or.i2new.lt.1)K=K+1
C       write(*,'(20F5.2)')(Mask(i,1),i=1,Nbin)

       select case(Wkey)
       case(0)
C       const weights
C        Denom=dble(K)
         Denom=sum(Mask)
C        do 10 i=1,K
        PWL=1.D0/Denom
C        if(incKey.eq.0) then
C         do 11 i=1,L1
C          PW(i,1)=PWL
C          Sum=Sum+PWL
C11       continue
C        else
         do 10 i=1,L1
         do 10 j=1,L2
C          if(WLH(i,j).ne.0.or.(i.eq.i1new.and.j.eq.i2new)) then
           PW(i,j)=PWL*Mask(i,j)
           Sumc=Sumc+PWL*Mask(i,j)
C          else
C           PW(i,j)=0.D0
C          endif
10       continue
C        endif
       case(1)
C      linear weights
C        stop'Net obrabotki vesa 1'
        do 20 i=1,L2
C        do 20 k=1,L1
         Denom=Denom+dble(i)
         do 20 k=1,L1
         Den2(i)=Den2(i)+Mask(k,i)
C         Denom=Denom+dble(i)/L1*Mask(k,i)
20      continue

C        do 20 i=1,K
C         Denom=Denom+dble(i)
C20      continue
        do 21 k=1,L2
        PWL=(L2-k+1)/Denom
        do 21 i=1,L1
C         PWL=dble(K+1-i)/Denom

          PW(i,k)=PWL/Den2(k)*Mask(i,k)
         Sumc=Sumc+PWL/Den2(k)*Mask(i,k)
21      continue

C        do 21 i=1,K
C         PWL=dble(K+1-i)/Denom
C          PW(i,1)=PWL
C         Sumc=Sumc+PWL
C21      continue
       case(2)
C      quadratic weights
C        stop'Net obrabotki vesa 2'
        do 30 i=1,L2
         Denom=Denom+dble(i*i)
        do 30 k=1,L1
         Den2(i)=Den2(i)+Mask(k,i)
30      continue
        do 31 i=1,L2
         PWL=dble((L2+1-i)*(L2+1-i))/Denom
        do 31 k=1,L1
          PW(k,i)=PWL/Den2(i)*Mask(k,i)
         Sumc=Sumc+PWL/Den2(i)*Mask(k,i)
31      continue

C        do 30 i=1,K
C         Denom=Denom+dble(i*i)
C30      continue
C        do 31 i=1,K
C         PWL=dble((K+1-i)*(K+1-i))/Denom
C          PW(i,1)=PWL
C         Sumc=Sumc+PWL
C31      continue
       case(3)
        do 40 i=1,L2
         Denom=Denom+dexp(dble(i))
        do 40 k=1,L1
         Den2(i)=Den2(i)+Mask(k,i)
40      continue
        do 41 i=1,L2
         PWL=dexp(dble(L2+1-i))/Denom
        do 41 k=1,L1
          PW(k,i)=PWL/Den2(i)*Mask(k,i)
         Sumc=Sumc+PWL/Den2(i)*Mask(k,i)
41      continue

       case(4)
        do 50 i=1,L2
         Denom=Denom+dexp(dble(i))
        do 50 k=1,L1
         if(Mask(k,i).ge.0.5) Den2(i)=Den2(i)+dble((L1+1-k)*(L1+1-k))
50      continue
        do 51 i=1,L2
         PWL=dexp(dble(L2+1-i))/Denom
        do 51 k=1,L1
         if(Mask(k,i).ge.0.5) then
          PW(k,i)=PWL/Den2(i)*dble((L1+1-k)*(L1+1-k))
          Sumc=Sumc+PWL/Den2(i)*dble((L1+1-k)*(L1+1-k))
         else
          PW(k,i)=0.D0
         endif
51      continue

       case(5)
C      double quadratic weights
        do 60 i=1,L2
         Denom=Denom+dble(i*i)
        do 60 k=1,L1
         if(Mask(k,i).ge.0.5) Den2(i)=Den2(i)+dble((L1+1-k)*(L1+1-k))
60      continue
        do 61 i=1,L2
         PWL=dble((L2+1-i)*(L2+1-i))/Denom
        do 61 k=1,L1
         if(Mask(k,i).ge.0.5) then
          PW(k,i)=PWL/Den2(i)*dble((L1+1-k)*(L1+1-k))
          Sumc=Sumc+PWL/Den2(i)*dble((L1+1-k)*(L1+1-k))
         else
          PW(k,i)=0.D0
         endif
61      continue

       case(6)
C      double linaer weights
        do 70 i=1,L2
         Denom=Denom+dble(i)
        do 70 k=1,L1
         if(Mask(k,i).ge.0.5) Den2(i)=Den2(i)+dble(L1+1-k)
70      continue
        do 71 i=1,L2
         PWL=dble(L2+1-i)/Denom
        do 71 k=1,L1
         if(Mask(k,i).ge.0.5) then
          PW(k,i)=PWL/Den2(i)*dble(L1+1-k)
          Sumc=Sumc+PWL/Den2(i)*dble(L1+1-k)
         else
          PW(k,i)=0.D0
         endif
71      continue


       case default
       write(*,*)'Durdulet'
       end select

C      open(15,FILE='PW.dat')
C      open(20,FILE='Mask.dat')

C      do 80 i=1,L1
C       write(15,'(100F10.7)') (PW(i,k),k=1,L2)
C       write(20,'(100F10.7)') (Mask(i,k),k=1,L2)
C80    continue
C      close(15)
C      close(20)
C      stop''

      Acc=1.D-6
      if (Sumc.lt.1.D0-Acc.or.Sumc.gt.1.D0+Acc) then
       write(*,*)Sumc
       stop 'incorrect sum of weights'
      endif
      end


      subroutine Umbrella
      use DynamicResize
      implicit none
      include 'var.common'

      WLG=WLG+Dlog(dble(WLH))
      WLH=int(Mask)
      return
      end
