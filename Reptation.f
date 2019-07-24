      Subroutine Reptation(N,X,V,neighb,lastneig,CopType,Phi,CopPhi,
     *                     ENV,EST,WLFact,NSTEP)
      implicit none
      include 'var.common'

      real pi
      parameter (pi=3.1415926)

      integer (kind=8)NSTEP
      integer i,k,N,J,ECKey,BrdCop(128),BrdN,CopType(N)
      integer lastneig(N),lntmp(N),neigtmp(N,N),Neighb(N,N),tmp(N)
      integer ENV,EST,DENV,DEST,ENVn,ESTn,Efull,IniRep,WLResult
      integer ngh,oe,nch,e1,e2,ReptKey
      integer stENV,stEST

      real rand(5),rnd1,rnd2,rnd3,rnd4,rnd5
      real*8 xs,ys,zs,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,l,f1,f2,f3
      real*8 X(3,N),V(3,N),XI,YI,ZI,Phi(N),CopPhi(2,N),COS1,atet,aphi
      real*8 DX,DY,DZ,RX,RY,RZ,R2,WLFact,XTMP(3,N),PhiTmp(N),VXI,VYI,VZI
      real*8 ro,atetmax,costet,csin,ccos,rndctrl(100)

      Data BrdN,IniRep /0,0/
      Data BrdCop /128*0/


      call CheckReptation(ReptKey)
      if(ReptKey.ne.0)return

      if(IniRep.eq.0) then
C     Identifying monomers on borders of blocks
C     mojet menyats9 ot modeli (Est ili ENV )
       BrdCop=0
       BrdN=0
       IniRep=1
       do 30 i=1,N-1
        if(CopType(i).ne.CopType(i+1)) then
         BrdN=BrdN+1
         BrdCop(BrdN)=i
        endif
30     continue
C      dumat s goto nachalo
C       return
      endif

      ECKey=0
      call ranmar(rand,5)
      rnd1=rand(1)
      rnd2=rand(2)
      rnd3=rand(3)
      rnd4=rand(4)
      rnd5=rand(5)
C      rnd4=0.0
      if(rnd4.ge.0.5) then
C     delete Nth, add 1st
       ngh=1
       oe=N
       nch=1
       e1=1
       e2=0
C       ps1=1
      else
C     delete 1st, add Nth
       ngh=N
       oe=1
       nch=-1
       e1=0
       e2=1
C       ps1=N
      endif

C       if(NSTEP.eq.1) then
C        open(80,FILE='RepErr1.dat')
C        do i=1,N
C         write(80,'(66I4)')i,(neighb(i,k),k=1,10)
C        enddo
C        close(80)
C        write(*,*)J,oe
C       endif

      DX=X(1,ngh)-X(1,ngh+nch)
      DY=X(2,ngh)-X(2,ngh+nch)
      DZ=X(3,ngh)-X(3,ngh+nch)
      R2=DX**2+DY**2+DZ**2
      R2=dSQRT(R2)
      c1=DX/R2
      c2=DY/R2
      c3=DZ/R2
C     building normal to O'z' = A
      d1=0.D0
      d2=0.D0
      d3=0.D0
      if(c1.gt.0.8)then
       d2=1.D0
      else
       d1=1.D0
      endif
      f1=c2*d3-c3*d2
      f2=c3*d1-c1*d3
      f3=c1*d2-c2*d1

      l=dsqrt((f2*c1-f1*c2)**2+(f3*c2-f2*c3)**2+(f1*c3-f3*c1)**2)
      a1=(f3*c2-f2*c3)/l
      a2=(f1*c3-f3*c1)/l
      a3=(f2*c1-f1*c2)/l
C      building normal to O'z' and O'x' = B
      l=dsqrt((a2*c1-a1*c2)**2+(a1*c3-a3*c1)**2+(a3*c2-a2*c3)**2)
      b1=(a3*c2-a2*c3)/l
      b2=(a1*c3-a3*c1)/l
      b3=(a2*c1-a1*c2)/l

C      len2=0.512+rand(ind+3)*1.441125
C      len2=len2**(1.D0/3.D0)
C      ccos2=(rand(ind+4)-0.5)*2.D0
C      csin2=sqrt(1.D0-ccos2**2)
C      aphi2=rand(ind+5)*360.D0/180.D0*pi


      ro=BndMin3+rnd1*(BndMax3-BndMin3)
      ro=ro**(1.D0/3.D0)
C      ro=1.D0
C      costet=-(SIG**2-R2**2-ro**2)/(2.d0*R2*ro)
C      ccos=-1.D0+rnd2*(1.D0+costet+0.2)
C      if(ccos.gt.1.D0) ccos=1.D0
      ccos=-1.D0+rnd2*(1.D0+Mxctet)
C      ccos=costet
C      write(*,*)costet,ccos

C      if(MPrnk.eq.0) ccos=(rnd2-0.5)*2.D0
C      i=int((ccos+1.D0)/2.D0*100.D0)+1
C      rndctrl(i)=rndctrl(i)+1.D0

      csin=dsqrt(1.D0-ccos**2)
      aphi=2.d0*PI*rnd3

C      aphi=rand(ind+2)*360.D0/180.D0*pi

C      costet=-(SIG**2-R2**2-ro**2)/(2.d0*R2*ro)
C      atetmax=acos(costet)
C      atet=(PI-atetmax)*dble(rnd2)
C      atet=PI*rnd2


      xs=ro*csin*DCOS(aphi)
      ys=ro*csin*DSIN(aphi)
      zs=-ro*ccos

C      XI=X(1,ngh)+BOND*DSIN(atet)*DCOS(aphi)
C      YI=X(2,ngh)+BOND*DSIN(atet)*DSIN(aphi)
C      ZI=X(3,ngh)+BOND*DCOS(atet)

      XI=X(1,ngh)+xs*a1+ys*b1+zs*c1
      YI=X(2,ngh)+xs*a2+ys*b2+zs*c2
      ZI=X(3,ngh)+xs*a3+ys*b3+zs*c3

C      call angles(XI,YI,ZI,X(1,ngh),X(2,ngh),X(3,ngh),X(1,ngh+nch),
C     * X(2,ngh+nch),X(3,ngh+nch),COS1)

C      write(*,'(20F10.3)')atet,atetmax,XI,YI,ZI,X(1,ngh),X(2,ngh),
C     * X(3,ngh),ro

      VXI=XI-AX*(DINT(2.D0*XI/AX)-DINT(XI/AX))
      VYI=YI-AY*(DINT(2.D0*YI/AY)-DINT(YI/AY))
      VZI=ZI-AZ*(DINT(2.D0*ZI/AZ)-DINT(ZI/AZ))

      DENV=0
      DEST=0
      neigtmp=0
      lntmp=0

      do 100 i=1,N-1
       neigtmp(i+e1,1:lastneig(i+e2))=neighb(i+e2,1:lastneig(i+e2))+nch
       lntmp(i+e1)=lastneig(i+e2)
100   continue
       lntmp(ngh)=0

C     Old NV Energy
C      write(*,*)'Opa!'
C      lastneig(oe)=0
C       if(NSTEP.eq.152647) then
C        open(80,FILE='RepErr.dat')
C        do 15 i=1,N
C         write(80,'(66I4)')i,(neighb(i,j),j=1,lastneig(i))
C15      continue
C        close(80)
C       endif

C       if(NSTEP.eq.1) then
C        open(80,FILE='RepErr2.dat')
C        do i=1,N
C         write(80,'(66I4)')i,(neigtmp(i,k),k=1,10)
C        enddo
C        close(80)
C        write(*,*)J,oe
C       endif


      do 10 i=1,lastneig(oe)
C       write(*,*)'Opachki!'
C       stop'stope'
       J=neighb(oe,i)
       if(oe.eq.N)then
        neigtmp(j+1,lntmp(j+1))=0
        lntmp(j+1)=lntmp(j+1)-1
       else
        if(lntmp(j-1).gt.1)then
         if(j-1.gt.N.or.j-1.lt.1) write(*,*)'j-1=',j-1,lntmp(j-1)
         if(lntmp(j-1).gt.N) write(*,*)'lntmp=',lntmp(j-1),J,oe,NSTEP
         neigtmp(j-1,1:lntmp(j-1)-1)=neigtmp(j-1,2:lntmp(j-1))
         neigtmp(j-1,lntmp(j-1))=0
         lntmp(j-1)=lntmp(j-1)-1

        else
         neigtmp(j-1,1)=0
         lntmp(j-1)=0
        endif
       endif
       DX=V(1,J)-V(1,oe)
       DY=V(2,J)-V(2,oe)
       DZ=V(3,J)-V(3,oe)
       RX=DX-AX*DINT(2.D0*DX/AX)
       RY=DY-AY*DINT(2.D0*DY/AY)
       RZ=DZ-AZ*DINT(2.D0*DZ/AZ)
       R2=RX*RX+RY*RY+RZ*RZ
       R2=dsqrt(R2)
       if(R2.LE.RCUT) then
          if (CopType(oe).eq.1) then
             if (CopType(J).eq.1) then
               DENV=DENV+EPS
             else  
               DENV=DENV+2*EPS
             endif
          else 
             if (CopType(J).eq.1) then
               DENV=DENV+2*EPS
             else 
               DENV=DENV+4*EPS
             endif
          endif
      endif
10    continue

C     New NV energy and new neighb list
C      do 20 i=2,N-1
      i=2
      do 20 while (i.le.N-1)
       DX=VXI-V(1,i)
       DY=VYI-V(2,i)
       DZ=VZI-V(3,i)
       RX=DX-AX*DINT(2.D0*DX/AX)
       RY=DY-AY*DINT(2.D0*DY/AY)
       RZ=DZ-AZ*DINT(2.D0*DZ/AZ)
       R2=RX*RX+RY*RY+RZ*RZ
       R2=dsqrt(R2)
       if (R2.LE.RMAX) then
        lntmp(ngh)=lntmp(ngh)+1
        k=lntmp(i+nch)
        if(oe.eq.N) then
         tmp(1:k)=neigtmp(i+e1,1:k)
         neigtmp(i+e1,2:k+1)=tmp(1:k)
         neigtmp(i+e1,1)=1
        else
         neigtmp(i+nch,k+1)=N
        endif
        neigtmp(ngh,lntmp(ngh))=i+nch
        lntmp(i+nch)=lntmp(i+nch)+1
        if(R2.GT.SIG.and.R2.LE.RCUT) then
          if (CopType(i).eq.1) then
             if (CopType(ngh).eq.1) then
               DENV=DENV-EPS
             else  
               DENV=DENV-2*EPS
             endif
          else 
             if (CopType(ngh).eq.1) then
               DENV=DENV-2*EPS
             else 
               DENV=DENV-4*EPS
             endif
          endif
        endif
        if(R2.LE.SIG) then
         ECKey=1
         goto 60
        endif
       else
        i=i+(R2-RMAX)/BndMax
       endif
       i=i+1
20    continue

C     old Stiffnes
      if(Phi(oe-nch).gt.CopPhi(1,oe-nch).or.Phi(oe-nch).lt.
     * CopPhi(2,oe-nch)) then
       DEST=DEST-EPSST
      endif
C     new Stiffnes
      call angles(XI,YI,ZI,X(1,ngh),X(2,ngh),X(3,ngh),X(1,ngh+nch),
     * X(2,ngh+nch),X(3,ngh+nch),COS1)
      if(COS1.gt.CopPhi(1,ngh+nch).or.COS1.lt.CopPhi(2,ngh+nch))
     * DEST=DEST+EPSST

C     change of EStif due to switching monomer types
      do 50 i=1,BrdN
C     old Stiffnes
       k=BrdCop(i)
      if(Phi(k+e2).gt.CopPhi(1,k+e2).or.Phi(k+e2).lt.CopPhi(2,k+e2))then
        DEST=DEST-EPSST
       endif
      if(Phi(k+e2).gt.CopPhi(1,k+e1).or.Phi(k+e2).lt.CopPhi(2,k+e1))then
         DEST=DEST+EPSST
       endif
C       endif
50    continue
C       if(NSTEP.eq.23865)write(*,*)Mprnk,NSTEP,cntEB,'rept pered SAMC'
C       if(NSTEP.eq.1) then
C        open(80,FILE='RepErr3.dat')
C        do i=1,N
C         write(80,'(66I4)')i,(neigtmp(i,k),k=1,10)
C        enddo
C        close(80)
C        write(*,*)J,oe
C       endif

      ENVn=ENV+DENV
      ESTn=EST+DEST
60    call SAMC(ENV,ENVn,EST,ESTn,ECKey,rnd5,WLResult,WLFact,cntRep,
     *          cntNRep,NSTEP+SAMCPST,N)
C       if(NSTEP.eq.23865)write(*,*)Mprnk,NSTEP,cntEB,'rept posle SAMC'
C       stENV=40-(Emaxbnd-ENV)
C       stEst=Est+1

      if(WLResult.EQ.1) then
C     ***********Reptation accepted*************************************
C       if(NSTEP.eq.23865)write(*,*)Mprnk,NSTEP,cntEB,'rept posle sRAR'
C       stRAR(3,1,stENV,stEst)=stRAR(3,1,stENV,stEst)+1
C       stRAR(3,2,stENV,stEst)=stRAR(3,2,stENV,stEst)+ENVn-ENV
C       stRAR(3,3,stENV,stEst)=stRAR(3,3,stENV,stEst)+Estn-Est

       XTMP=X
       X(1:3,1+e1:N-e2)=XTMP(1:3,1+e2:N-e1)
       XTMP=V
       V(1:3,1+e1:N-e2)=XTMP(1:3,1+e2:N-e1)
       X(1,ngh)=XI
       X(2,ngh)=YI
       X(3,ngh)=ZI
       V(1,ngh)=VXI
       V(2,ngh)=VYI
       V(3,ngh)=VZI

       PhiTmp(2+e1:N-1-e2)=Phi(2+e2:N-1-e1)
       PhiTmp(ngh+nch)=COS1
       PhiTmp(1)=Phi(1)
       PhiTmp(N)=Phi(N)
       Phi=PhiTmp
C       if(ABS(ccos-cos1).gt.1.D-5)write(*,*)ccos,cos1
C      i=int((ccos+1.D0)/2.D0*100.D0)+1
C      rndctrl(i)=rndctrl(i)+1.D0

C       R2=sqrt((X(1,1)-X(1,3))**2+(X(2,1)-X(2,3))**2+(X(3,1)-X(3,3))**2)
C       if(Abs(R2-SIG).gt.1.D-6)write(*,*)R2,'Durulet'
C       write(*,*)ccos,cos1,ccos-cos1


C       call TABLE(N,NC,RMAX,AX,AY,AZ,V,neighb,lastneig)
C       write(*,*)''
C       do 15 i=1,N
C        if(lastneig(i).ne.lntmp(i)) then
C         write(*,*)'Err: lastneig',i
C         write(*,'(20I4)')(lastneig(k),k=1,N)
C         write(*,'(20I4)')(lntmp(k),k=1,N)
C        endif
C        do 15 k=1,lastneig(i)
C        if(neighb(i,k).ne.neigtmp(i,k)) then
C         write(*,*)'Err: neig',i,k
C         write(*,'(20I4)')i,(neighb(i,j),j=1,lastneig(i))
C         write(*,'(20I4)')i,(neigtmp(i,j),j=1,lntmp(i))
C        endif
C15     continue


       neighb=neigtmp
       lastneig=lntmp

       ENV=ENVn
       EST=ESTn
       Efull=ENV+EST
C       if(NSTEP.eq.1) then
C        open(80,FILE='RepErr4.dat')
C        do i=1,N
C         write(80,'(66I4)')i,(neighb(i,k),k=1,10)
C        enddo
C        close(80)
C        write(*,*)J,oe
C       endif

       return
      endif
C       if(NSTEP.eq.23865)write(*,*)Mprnk,NSTEP,cntEB,'rept end'
C       stRAR(3,4,stENV,stEst)=stRAR(3,4,stENV,stEst)+1

      return
      

C     result of monomer type change (tolko esli vzaimod. otlichaetsya)
C     result of monomer type change ENV(tolko esli vzaimod. otlichaets)
C     Dodelat esli pomenyaets9 model
C      do 40 i=1,BrdN
C       do 41 j=1,lntmp(i)
C        k=neigtmp(i,j)
C        DX=V(1,i)-V(1,k)
C        DY=V(2,i)-V(2,k)
C        DZ=V(3,i)-V(3,k)
C        RX=DX-AX*DINT(2.D0*DX/AX)
C        RY=DY-AY*DINT(2.D0*DY/AY)
C        RZ=DZ-AZ*DINT(2.D0*DZ/AZ)
C        R2=RX*RX+RY*RY+RZ*RZ
C        R2=dsqrt(R2)
C        if(R2.GT.SIG.and.R2.LE.RCUT) then
C     old energy
C     Dodelat !!! Smotret 4tobi k-ii ne izmenilsya toje
C         if(CopType(i).ne.CopType(k)) DENV=DENV-3*EPS
C         if(CopType(i).eq.CopType(k).and.CopType(k).eq.1) then
C          DENV=DENV+4*EPS
C         endif
C         DENV=DENV+EPS
C        endif
C41     continue
C40    continue

      end
