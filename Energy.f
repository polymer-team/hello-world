C  This version of energy file is for system with 2 energy contributions
C  Energy 1 is for Non-valent, energy 2 is for stiffness.
C  There are 3 additional moves: reptation, CCB, Endbridging
      subroutine Energy1(mono,N,X,V,neig,ln,VXI,VYI,VZI,EKey,DENV,
     * CopType)
      implicit none
      include 'var.common'

      integer N,mono,neig(N,N),ln(N),EKey,i,j,DENV,k,CopType(N)
      real*8 V(3,N),VXI,VYI,VZI,DXn,DYn,DZn,DXo,DYo,DZo,RX,RY,RZ,R2
      real*8 SIG2,RCUT2,X(3,N)
      
      SIG2=SIG*SIG
      RCUT2=RCUT*RCUT

      DENV=0
      do 10 I=1,ln(mono)
       J=neig(mono,I)

       DXo=V(1,J)-V(1,mono)
       DYo=V(2,J)-V(2,mono)
       DZo=V(3,J)-V(3,mono)
C     Old Energy
       call CalcDV(RX,RY,RZ,DXo,DYo,DZo)
       R2=RX*RX+RY*RY+RZ*RZ
       if(R2.GT.RTop.or.(R2.GT.Rbot2.and.R2.LT.Rbot)) cycle
       if(R2.GT.SIG2.and.R2.LE.RCUT2) then
          if (CopType(mono).eq.1) then
             if (CopType(J).eq.1) then
               DENV=DENV+EPS
             else
	       DENV=DENV+2*EPS
             endif
          else 
             if  (CopType(J).eq.1) then
               DENV=DENV+2*EPS
             else
               DENV=DENV+4*EPS
             endif
          endif
       endif
       if (R2.LE.SIG2) then
        call OutErrorEnergyNV(mono,J,N,R2,V,neig,ln)
       endif


       DXn=V(1,J)-VXI
       DYn=V(2,J)-VYI
       DZn=V(3,J)-VZI
C     New Energy
       call CalcDV(RX,RY,RZ,DXn,DYn,DZn)
       R2=RX*RX+RY*RY+RZ*RZ
       if(R2.GT.SIG2.and.R2.LE.RCUT2) then
          if (CopType(mono).eq.1) then
             if (CopType(J).eq.1) then
               DENV=DENV-EPS
             else  
               DENV=DENV-2*EPS
             endif
          else 
             if  (CopType(J).eq.1) then
               DENV=DENV-2*EPS
             else  
               DENV=DENV-4*EPS
             endif
          endif
       endif
       if (R2.LE.SIG2) then
        EKey=Ekey+1
        return
       endif
10    continue
      end

      subroutine Energy2(mono,N,X,XI,YI,ZI,DEstif,StifKey,Phi,
     * CopPhi,COS1,COS2,COS3)
      implicit none
      include 'var.common'

      integer mono,N,StifKey(N),m,DEstif,Estifctrl
      real*8 X(3,N),XI,YI,ZI,COS1,COS2,COS3,Phi(N),CopPhi(2,N)
      m=mono
      DEstif=0
C     new angles
      select case(StifKey(mono))
       case(-2)
        call angles(XI,YI,ZI,X(1,m+1),X(2,m+1),X(3,m+1),
     *  X(1,m+2),X(2,m+2),X(3,m+2),COS1)
C        if (Phi(m+1).gt.CopPhi(m+1)) DEstif=DEstif-1
C        if (COS1.gt.CopPhi(m+1)) DEstif=DEstif+1

        if(Phi(m+1).le.CopPhi(1,m+1).and.Phi(m+1).ge.CopPhi(2,m+1)) then
         DEStif=DEStif+EPSST
        endif
        if(COS1.le.CopPhi(1,m+1).and.COS1.ge.CopPhi(2,m+1)) then
         DEStif=DEStif-EPSST
        endif

       case(-1)
        call angles(XI,YI,ZI,X(1,m+1),X(2,m+1),X(3,m+1),
     *  X(1,m+2),X(2,m+2),X(3,m+2),COS1)
        call angles(X(1,m-1),X(2,m-1),X(3,m-1),XI,YI,ZI,
     *  X(1,m+1),X(2,m+1),X(3,m+1),COS2)
C        if (Phi(m+1).gt.CopPhi(m+1)) DEstif=DEstif-1
C        if (COS1.gt.CopPhi(m+1)) DEstif=DEstif+1
        if(Phi(m+1).le.CopPhi(1,m+1).and.Phi(m+1).ge.CopPhi(2,m+1)) then
         DEStif=DEStif+EPSST
        endif
        if(COS1.le.CopPhi(1,m+1).and.COS1.ge.CopPhi(2,m+1)) then
         DEStif=DEStif-EPSST
        endif
        if(Phi(m).le.CopPhi(1,m).and.Phi(m).ge.CopPhi(2,m)) then
         DEStif=DEStif+EPSST
        endif
        if(COS2.le.CopPhi(1,m).and.COS2.ge.CopPhi(2,m)) then
         DEStif=DEStif-EPSST
        endif

       case(0)
        call angles(XI,YI,ZI,X(1,m+1),X(2,m+1),X(3,m+1),
     *  X(1,m+2),X(2,m+2),X(3,m+2),COS1)
        call angles(X(1,m-1),X(2,m-1),X(3,m-1),XI,YI,ZI,
     *  X(1,m+1),X(2,m+1),X(3,m+1),COS2)
        call angles(X(1,m-2),X(2,m-2),X(3,m-2),X(1,m-1),X(2,m-1),
     *  X(3,m-1),XI,YI,ZI,COS3)

        if(Phi(m+1).le.CopPhi(1,m+1).and.Phi(m+1).ge.CopPhi(2,m+1)) then
         DEStif=DEStif+EPSST
        endif
        if(COS1.le.CopPhi(1,m+1).and.COS1.ge.CopPhi(2,m+1)) then
         DEStif=DEStif-EPSST
        endif
        if(Phi(m).le.CopPhi(1,m).and.Phi(m).ge.CopPhi(2,m)) then
         DEStif=DEStif+EPSST
        endif
        if(COS2.le.CopPhi(1,m).and.COS2.ge.CopPhi(2,m)) then
         DEStif=DEStif-EPSST
        endif
        if(Phi(m-1).le.CopPhi(1,m-1).and.Phi(m-1).ge.CopPhi(2,m-1)) then
         DEStif=DEStif+EPSST
        endif
        if(COS3.le.CopPhi(1,m-1).and.COS3.ge.CopPhi(2,m-1)) then
         DEStif=DEStif-EPSST
        endif

       case(1)
        call angles(X(1,m-1),X(2,m-1),X(3,m-1),XI,YI,ZI,
     *  X(1,m+1),X(2,m+1),X(3,m+1),COS2)
        call angles(X(1,m-2),X(2,m-2),X(3,m-2),X(1,m-1),X(2,m-1),
     *  X(3,m-1),XI,YI,ZI,COS3)
        if(Phi(m).le.CopPhi(1,m).and.Phi(m).ge.CopPhi(2,m)) then
         DEStif=DEStif+EPSST
        endif
        if(COS2.le.CopPhi(1,m).and.COS2.ge.CopPhi(2,m)) then
         DEStif=DEStif-EPSST
        endif
        if(Phi(m-1).le.CopPhi(1,m-1).and.Phi(m-1).ge.CopPhi(2,m-1)) then
         DEStif=DEStif+EPSST
        endif
        if(COS3.le.CopPhi(1,m-1).and.COS3.ge.CopPhi(2,m-1)) then
         DEStif=DEStif-EPSST
        endif

       case(2)
        call angles(X(1,m-2),X(2,m-2),X(3,m-2),X(1,m-1),X(2,m-1),
     *  X(3,m-1),XI,YI,ZI,COS3)
        if(Phi(m-1).le.CopPhi(1,m-1).and.Phi(m-1).ge.CopPhi(2,m-1)) then
         DEStif=DEStif+EPSST
        endif
        if(COS3.le.CopPhi(1,m-1).and.COS3.ge.CopPhi(2,m-1)) then
         DEStif=DEStif-EPSST
        endif

       case(3)
        call angles(X(1,m-1),X(2,m-1),X(3,m-1),XI,YI,ZI,
     *  X(1,m+1),X(2,m+1),X(3,m+1),COS2)
        if(Phi(m).le.CopPhi(1,m).and.Phi(m).ge.CopPhi(2,m)) then
         DEStif=DEStif+EPSST
        endif
        if(COS2.le.CopPhi(1,m).and.COS2.ge.CopPhi(2,m)) then
         DEStif=DEStif-EPSST
        endif
       case default
      end select
      end



      SubRoutine E1Full(N,X,ENV,Ekey,CopType)
      implicit none
      include 'var.common'
      integer i,j,N,NB,ENV,Ekey,CopType(N)
      real*8 X(3,N),V(3,N),DX,DY,DZ,RX,RY,RZ,R2,SIG2,RCUT2

      ENV=0
      NB=N/NC
      SIG2=SIG*SIG
      RCUT2=RCUT*RCUT
      
      call CalcV(N,X,V)

      do 10 i=1,N-1
       do 11 j=i+1,N
C     not from same chain!!!!
        if(((j-i).eq.1).and.((i-1)/NB.Eq.(j-1)/NB)) cycle
        DX=V(1,J)-V(1,I)
        DY=V(2,J)-V(2,I)
        DZ=V(3,J)-V(3,I)

        call CalcDV(RX,RY,RZ,DX,DY,DZ)
        R2=RX*RX+RY*RY+RZ*RZ

      if (R2.GT.SIG2.and.R2.LE.RCUT2) then
          if (CopType(i).eq.1) then
             if (CopType(j).eq.1) then
               ENV=ENV-EPS
             else  
               ENV=ENV-2*EPS
             endif
          else 
             if  (CopType(j).eq.1) then
               ENV=ENV-2*EPS
             else  
               ENV=ENV-4*EPS
             endif
          endif
      endif
        if (R2.LE.SIG2) then
         Write(*,*)'E1Full - Distance Error: Mono',i,j,MPrnk
        endif
11     continue
10    continue
      end



      SubRoutine E2Full(N,X,Estif,Phi,CopPhi,CopType,StifKey,Ekey)
      implicit none
      include 'var.common'
      integer i,j,N,NB,StifAr(N),CopType(N),Estif,StifKey(N),SK,Ekey
      real*8 X(3,N),V(3,N),DX,DY,DZ,RX,RY,RZ,R2
      real*8 Ksi(N),Phi(N),CopPhi(2,N)

      Estif=0
      call AllAngles(N,NB,Ksi,X)
      do 12 i=1,N
       SK=StifKey(i)
       if(SK.eq.-2.or.SK.eq.2.or.SK.ge.4)cycle
       if(Ksi(i).gt.CopPhi(1,i).or.Ksi(i).lt.CopPhi(2,i)) then
        Estif=Estif+EPSST
       endif
12    continue
      Phi=Ksi
      end

      
      subroutine EnergyCCB(N,V,ln,neig,DENV,DEST,POS1,POS2,Phi,CopPhi,
     * CopType)
      implicit none
      include 'var.common'
      
      integer N,DENV,DEST,POS1,POS2,ln(N),neig(N,N),CopType(N)
      integer I,k,m,ch,j
      real*8 V(3,N),RCUT2,DX,DY,DZ,RX,RY,RZ,R,Phi(N),CopPhi(2,N)
      
      RCUT2=RCUT*RCUT
      DENV=0
      DEST=0
      ch=1
      if(POS2.ne.N) ch=-1

      DO 20 I=POS1,POS2,ch
       do 10 k=1,ln(i)
        m=neig(i,k)
        if(m.ge.pos1.and.m.le.i+1.and.pos2.eq.N) cycle
        if(m.le.pos1.and.m.ge.i-1.and.pos2.eq.1) cycle
        DX=V(1,i)-V(1,m)
        DY=V(2,i)-V(2,m)
        DZ=V(3,i)-V(3,m)
        call CalcDV(RX,RY,RZ,DX,DY,DZ)
        R=RX*RX+RY*RY+RZ*RZ
        if (R.LE.RCUT2) then
          if (CopType(i).eq.1) then
             if (CopType(m).eq.1) then
               DENV=DENV-EPS
             else  
               DENV=DENV-2*EPS
             endif
          else 
             if (CopType(m).eq.1) then
               DENV=DENV-2*EPS
             else  
               DENV=DENV-4*EPS
             endif
          endif
       endif
        
10     continue
       J=I-ch
       if(Phi(J).gt.CopPhi(1,J).or.Phi(J).lt.CopPhi(2,J)) then
        DEST=DEST+EPSST
       endif
20    continue
      end
      
      subroutine EnergyEB(N,X,DENV,DEST,Phi,CopPhi,j,oe,ch,COS1,COS2,
     * CopType)
      implicit none
      include 'var.common'
      integer N,DENV,DEST,i,j,k,oe,ch,CopType(N),xtemp
      real*8 X(3,N),Phi(N),CopPhi(2,N),COS1,COS2,XI,YI,ZI,R,R1
      
      do 101 i=oe+ch,J-2*ch,ch
      if(Phi(i).gt.CopPhi(1,i).or.Phi(i).lt.CopPhi(2,i)) DEST=DEST-EPSST
      k=J-ch+oe-i
      if(Phi(k).gt.CopPhi(1,i).or.Phi(k).lt.CopPhi(2,i)) DEST=DEST+EPSST
101   continue

      if(Phi(J-ch).gt.CopPhi(1,J-ch).or.Phi(J-ch).lt.CopPhi(2,J-ch))then
       DEST=DEST-EPSST
      endif
      if(Phi(J).gt.CopPhi(1,J).or.Phi(J).lt.CopPhi(2,J)) DEST=DEST-EPSST

      call ANGLES(X(1,oe+ch),X(2,oe+ch),X(3,oe+ch),X(1,oe),X(2,oe),
     *  X(3,oe),X(1,J),X(2,J),X(3,J),COS1)
      if(COS1.gt.CopPhi(1,J-ch).or.COS1.lt.CopPhi(2,J-ch)) then
       DEST=DEST+EPSST
      endif
      call ANGLES(X(1,oe),X(2,oe),X(3,oe),X(1,J),X(2,J),X(3,J),
     *  X(1,J+ch),X(2,J+ch),X(3,J+ch),COS2)
      if(COS2.gt.CopPhi(1,J).or.COS2.lt.CopPhi(2,J)) DEST=DEST+EPSST

      XI=X(1,oe)-X(1,J)
      YI=X(2,oe)-X(2,J)
      ZI=X(3,oe)-X(3,J)
      R=dsqrt(XI*XI+YI*YI+ZI*ZI)/BOND
      if(R.lt.RCUT)then
          if (CopType(oe).eq.1) then
             if (CopType(J).eq.1) then
               DENV=DENV+EPS
             else  
               DENV=DENV+2*EPS
             endif
          else 
               if  (CopType(J).eq.1) then
                DENV=DENV+2*EPS
               else  
                DENV=DENV+4*EPS
               endif
          endif
      endif
      XI=X(1,J-ch)-X(1,J)
      YI=X(2,J-ch)-X(2,J)
      ZI=X(3,J-ch)-X(3,J)
      R1=dsqrt(XI*XI+YI*YI+ZI*ZI)/BOND
      xtemp=J-ch
      if(R1.lt.RCUT) then
          if (CopType(xtemp).eq.1) then
             if (CopType(J).eq.1) then
               DENV=DENV-EPS
             else  
               DENV=DENV-2*EPS
             endif
          else 
             if (CopType(J).eq.1) then
               DENV=DENV-2*EPS
               else  
               DENV=DENV-4*EPS
             endif
          endif
      endif
      end


      SubRoutine EnergyContr(N,V,CopType,neig,ln,EAA,EAB,EBB)
      implicit none
      include 'var.common'
      integer i,j,k,N,NB,CopType(N),EAA,EAB,EBB
      integer neig(n,n),ln(n)
      real*8 V(3,N),DX,DY,DZ,RX,RY,RZ,R2

      EAA=0
      EAB=0
      EBB=0
      do 10 i=1,N-2
       do 11 k=1,ln(i)
        j=neig(i,k)
        if(j.lt.i) cycle
        DX=V(1,J)-V(1,I)
        DY=V(2,J)-V(2,I)
        DZ=V(3,J)-V(3,I)
        call CalcDV(RX,RY,RZ,DX,DY,DZ)
        R2=RX*RX+RY*RY+RZ*RZ
        R2=dsqrt(R2)
         if(R2.LE.RCUT) then
          if(CopType(i).eq.1) then
           if(CopType(j).eq.1) then
            EAA=EAA+EPS
           else
            EAB=EAB+2*EPS
           endif
          else
           if(CopType(j).eq.1) then
            EAB=EAB+2*EPS
           else
            EBB=EBB+4*EPS
           endif
          endif
         endif
11     continue
10    continue
      return
      end


      subroutine StiffKeyArray(N,StifKey)
      implicit none
      include 'var.common'

      integer N,NB,i,nchain,mnmn,mnmx,StifKey(N)
      NB=N/NC
      if (NB.Gt.3) then
      do 100 i=1,N
       nchain=(i-1)/NB
       mnmn=nchain*NB+1
       mnmx=(nchain+1)*NB
       StifKey(i)=0
       if (i.eq.mnmn) StifKey(i)=-2
       if (i.eq.mnmn+1) StifKey(i)=-1
       if (i.eq.mnmx-1) StifKey(i)=1
       if (i.eq.mnmx) StifKey(i)=2
100    continue
      endif
      if (NB.eq.3) then
      do 101 i=1,N
       if (mod(i,3).eq.0) StifKey(i)=2
       if (mod(i,3).eq.1) StifKey(i)=-2
       if (mod(i,3).eq.2) StifKey(i)=3
101    continue
      endif
      if (NB.lt.3) then
       do 102 i=1,N
        StifKey(i)=4
102    continue
      endif
      end

