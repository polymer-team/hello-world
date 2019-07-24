      SUBROUTINE BEGIN(N,NB,X,rand0,CopType,Phi,StifKey,CopPhi)
      Implicit none
      include 'var.common'

      real pi
      parameter (pi=3.1415926)

      integer N,NB,I,J,k,CopType(N),tmpAB
      
      Real*8 X(3,N),rand0
      Real*8 Phi(N),StifKey(N),CopPhi(2,N),CosStifA1,CosStifA2,CosStifB
      character(Len=256) :: xcfname,xconfname
      
      write(xconfname,'(A5,I0.2,A4)')'x_cnf',MPrnk,'.dat'

      if (KEY.eq.0) then
       call BUILDCONF(N,NB,X,rand0)
C      X(1,1)=0.D0
C      X(2,1)=0.D0
C      X(3,1)=0.D0
C      X(1,2)=1.D0
C      X(2,2)=0.D0
C      X(3,2)=0.D0
C      X(1,3)=2.D0
C      X(2,3)=0.D0
C      X(3,3)=0.D0
C       write(*,*)'Begin: THIS IS LINEAR CONFORMATION!!!!!!'
       open(10,FILE=xconfname)
       do 41 I=1,N
        write(10,'(3E18.10)')X(1,I),X(2,I),X(3,I)
41     continue
       close(10)

       i=0
       do 110 k=1,NC
        do 120 j=1,NB
         i=i+1
         tmpAB=mod(j,(LenA+LenB))
         if ((tmpAB.LE.LenA).AND.(tmpAB.NE.0)) then
          CopType(i)=1
         else
          CopType(i)=2
         endif
120     continue
110    continue
       if(MPrnk.eq.0) then
        open(20,FILE='mc-copol.dat')
        do 130 k=1,N
         write(20,'(I1)') CopType(k)
130      continue
        close(20)
       endif
      endif
      if (KEY.ne.0)then
C       open(20,FILE='mc-copol.dat')
       open(10,FILE=xconfname)
       do 502 i=1,N
C        read(20,*)CopType(i)
        read(10,*)X(1,I),X(2,I),X(3,I)
        tmpAB=mod(i,(LenA+LenB))
        if ((tmpAB.LE.LenA).AND.(tmpAB.NE.0)) then
         CopType(i)=1
        else
         CopType(i)=2
        endif
502    continue
C       close(20)
       close(10)
C       call ToVRML(N,NB,X,SIG,CopType,'beg.wrl')
      endif

      CosStifA1=DCOS((AngStifA-DevStifA)*pi/180.D0)
      CosStifA2=DCOS((AngStifA+DevStifA)*pi/180.D0)
      if (AngStifA+DevStifA.gt.180.d0) CosStifA2=-1.5D0
      if (DevStifA.eq.0.d0) CosStifA2=-1.5D0
      CosStifB=DCOS(AngStifB*pi/180.D0)
C      write(*,*)CosStifA1,CosStifA2,CosStifB

      do 101 i=1,N
       if(CopType(i).eq.1) then
        CopPhi(1,i)=CosStifA1
        CopPhi(2,i)=CosStifA2
       endif
       if(CopType(i).eq.2) then
        CopPhi(1,i)=CosStifB
        CopPhi(2,i)=-1.5D0
       endif

C       if(CopType(i).eq.1) CopPhi(i)=CosStifA
C       if(CopType(i).eq.2) CopPhi(i)=CosStifB
101   continue
C      call ToVRML(N,NB,X,SIG,CopType,'xconf.wrl')
C      stop''

C      call AllAngles(N,NB,Phi,X)
      call StiffKeyArray(N,StifKey)

      RETURN
      END

      SUBROUTINE BUILDCONF(N,NB,X,rand0)
      implicit none
      include 'var.common'

      integer N,NB,I,J,k,CheckKey

      Real*8 X(3,N),rand0,ran,R1,R2,R3,R,XI,YI,ZI,VXI,VYI,VZI
      Real*8 V(3,N),DX,DY,DZ,RNDM

      ran=abs(rand0)-int(abs(rand0))
      CALL RNDM1(RAN)

      call FirstMono(XI,YI,ZI)
      X(1,1)=XI
      X(2,1)=YI
      X(3,1)=ZI

      DO 20 I=2,N
10     R1=1.D0-2.D0*RNDM(-1)
       R2=1.D0-2.D0*RNDM(-1)
       R3=1.D0-2.D0*RNDM(-1)
       R=R1*R1+R2*R2+R3*R3
C      write(*,*) 'BEGIN: I = ', I, '   R = ',R
       IF(R.GT.1.D0)GO TO 10
       R=DSQRT(R)
       J=I-1
       XI=X(1,J)+BOND*R1/R
       IF(DABS(XI).GE.AX/2.D0)GO TO 10
       YI=X(2,J)+BOND*R2/R
       IF(DABS(YI).GE.AY/2.D0)GO TO 10
       ZI=X(3,J)+BOND*R3/R
       IF(DABS(ZI).GE.AZ/2.D0)GO TO 10
       X(1,I)=XI
       X(2,I)=YI
       X(3,I)=ZI

       call CheckGeometry(XI,YI,ZI,CheckKey)
       if(CheckKey.ne.0) goto 10

       call CalcVI(XI,YI,ZI,VXI,VYI,VZI)

       V(1,I)=VXI
       V(2,I)=VYI
       V(3,I)=VZI

       do 100 k=1,I-1
        if ((i-k.EQ.1).and.((i-1)/NB.EQ.(k-1)/NB)) go to 100
        DX=V(1,k)-V(1,I)
        DY=V(2,k)-V(2,I)
        DZ=V(3,k)-V(3,I)
        
        call CalcDV(R1,R2,R3,DX,DY,DZ)

        R=R1*R1+R2*R2+R3*R3
        R=DSQRT(R)
        if (R.LE.SIG) go to 10
        DX=X(1,k)-X(1,I)
        DY=X(2,k)-X(2,I)
        DZ=X(3,k)-X(3,I)
        R=DX*DX+DY*DY+DZ*DZ
        R=DSQRT(R)
        if (R.LE.SIG) go to 10
100    continue
20    CONTINUE
      end

      SUBROUTINE RNDM1(RAND0)
      implicit none
      integer NRAND
      real*8 RAND0,RAND
      real*8 RNDM
      NRAND=1+INT(783637.D0*RAND0)
      RAND=RNDM(NRAND)
      RETURN
      END

      REAL*8 FUNCTION RNDM(KEY)
      implicit none
      integer KEY,NRAND
      data nrand/0/
      IF(KEY.GE.0)GO TO 10
      NRAND=125*NRAND
      NRAND=NRAND-(NRAND/2796203)*2796203
      RNDM=(NRAND)/2796202.D0
      RETURN
   10 NRAND=KEY
      RNDM=(NRAND)/2796202.D0
      RETURN
      END



