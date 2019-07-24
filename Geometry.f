      subroutine CalcV(N,X,V)
      implicit none
      include 'var.common'

      integer i,N
      real*8 X(3,N),V(3,N)
      
      do 1000 I=1,N
       V(1,I)=X(1,I)-AX*(DINT(2.D0*X(1,I)/AX)-DINT(X(1,I)/AX))
       V(2,I)=X(2,I)-AY*(DINT(2.D0*X(2,I)/AY)-DINT(X(2,I)/AY))
       V(3,I)=X(3,I)-AZ*(DINT(2.D0*X(3,I)/AZ)-DINT(X(3,I)/AZ))
1000  continue
      end
      
      
      subroutine CalcVI(XI,YI,ZI,VXI,VYI,VZI)
      implicit none
      include 'var.common'

      real*8 XI,YI,ZI,VXI,VYI,VZI

      VXI=XI-AX*(DINT(2.D0*XI/AX)-DINT(XI/AX))
      VYI=YI-AY*(DINT(2.D0*YI/AY)-DINT(YI/AY))
      VZI=ZI-AZ*(DINT(2.D0*ZI/AZ)-DINT(ZI/AZ))
      end
      
      subroutine CalcDV(RX,RY,RZ,DX,DY,DZ)
      implicit none
      include 'var.common'

      real*8 RX,RY,RZ,DX,DY,DZ

       RX=DX-AX*DINT(2.D0*DX/AX)
       RY=DY-AY*DINT(2.D0*DY/AY)
       RZ=DZ-AZ*DINT(2.D0*DZ/AZ)
      end

      subroutine CheckGeometry(XI,YI,ZI,CheckKey)
      implicit none
      integer CheckKey
      real*8 XI,YI,ZI

      CheckKey=0
      return
      end
      
      subroutine ChooseMono(NB,mono,j,rnd,nj)
      integer NB,mono,j,nj
      real rnd
      
      mono = j + int(real(nb)*rnd)

      if(mono.lt.j.or.mono.gt.nj)then
       write(*,*)'Wring mono ',mono
       stop 'wrong mono !!!'
      endif
      end
      
      subroutine FirstMono(X1,Y1,Z1)
      implicit none
      real*8 X1,Y1,Z1
      
      X1=0.D0
      Y1=0.D0
      Z1=0.D0
      end
      
      subroutine CheckReptation(ReptKey)
      implicit none
      integer ReptKey
      
C     if ReptKey is not 0, there will be no reptation moves
      ReptKey=0
      end
      
      subroutine POS2CCB(N,rnd2,POS2)
      implicit none
      real rnd2
      integer N,POS2
      
      if(rnd2.lt.0.5) then
       pos2=1
      else
       pos2=N
      endif
      end
      
      subroutine ChooseEndEB(N,X,ln,neig,BridgeCNT,BridgeNBH,m)
      implicit none
      include 'var.common'
      integer N,ln(N),neig(N,N),BridgeCNT,BridgeNBH(N),i,j,m
      real*8 X(3,N),R,XI,YI,ZI,R1

      if(ln(1).gt.0)then
       do 100 i=1,ln(1)
        J=neig(1,i)
        if(J.eq.3) cycle
        if(J.eq.N) exit
        XI=X(1,1)-X(1,J)
        YI=X(2,1)-X(2,J)
        ZI=X(3,1)-X(3,J)
        R=dsqrt(XI*XI+YI*YI+ZI*ZI)/BOND
        XI=X(1,J-1)-X(1,J)
        YI=X(2,J-1)-X(2,J)
        ZI=X(3,J-1)-X(3,J)
        R1=dsqrt(XI*XI+YI*YI+ZI*ZI)/BOND
        if(R.lt.BndMax.and.R1.gt.SIG) then
         BridgeCNT=BridgeCNT+1
         BridgeNBH(BridgeCNT)=J
        endif
100    continue
      endif
300   m=BridgeCNT
C      goto 200
      if(ln(N).gt.0)then
       do 110 i=1,ln(N)
        J=neig(N,i)
        if(J.eq.1) cycle
        if(J.eq.N-2) exit
        XI=X(1,N)-X(1,J)
        YI=X(2,N)-X(2,J)
        ZI=X(3,N)-X(3,J)
        R=dsqrt(XI*XI+YI*YI+ZI*ZI)/BOND
        XI=X(1,J+1)-X(1,J)
        YI=X(2,J+1)-X(2,J)
        ZI=X(3,J+1)-X(3,J)
        R1=dsqrt(XI*XI+YI*YI+ZI*ZI)/BOND
        if(R.lt.BndMax.and.R1.gt.SIG) then
         BridgeCNT=BridgeCNT+1
         BridgeNBH(BridgeCNT)=J
        endif
110    continue
      endif

      end
      
      subroutine BondLength(mono,N,X,Ekey,XI,YI,ZI)
      implicit none

      integer N,mono,NB,Ekey
      real*8 X(3,N),RX,RY,RZ,R2,XI,YI,ZI
      integer cnt

      DATA cnt/0/

      include 'var.common'

      cnt=cnt+1
      NB=N/NC
C      Write(*,*)'EB Ekey = ',Ekey
C       if(cnt.gt.3470*N) then
C      if(mono.eq.128)Write(*,*)'128:'
C      endif
      if(MOD((mono-1),NB).NE.0) then
C     checking left neighbour
        RX=XI-X(1,mono-1)
        RY=YI-X(2,mono-1)
        RZ=ZI-X(3,mono-1)
        R2=RX*RX+RY*RY+RZ*RZ
        R2=dsqrt(R2)
        R2=R2/BOND

        if(R2.LT.BndMin.or.R2.GT.BndMax) then
          Ekey=Ekey+1
        endif
      endif
      if(MOD(mono,NB).NE.0) then
C     Checking right neighbour
        RX=XI-X(1,mono+1)
        RY=YI-X(2,mono+1)
        RZ=ZI-X(3,mono+1)
        R2=RX*RX+RY*RY+RZ*RZ
        R2=dsqrt(R2)
        R2=R2/BOND
        if(R2.LT.BndMin.or.R2.GT.BndMax) then
          Ekey=Ekey+1
        endif
      endif
      end

