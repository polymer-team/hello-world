      subroutine ANGLES(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,COSPHI)
      implicit none
      
      real*8 X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,COSPHI
      real*8 L1,M1,N1,L2,M2,N2

      L1=X1-X2
      M1=Y1-Y2
      N1=Z1-Z2
C      L2=X2-X3
C      M2=Y2-Y3
C      N2=Z2-Z3
      L2=X3-X2
      M2=Y3-Y2
      N2=Z3-Z2

      COSPHI=(L1*L2+M1*M2+N1*N2)/(DSQRT(L1*L1+M1*M1+N1*N1)*DSQRT(L2*L2+
     *M2*M2+N2*N2))
      return
      end
      
      subroutine MAXANGL(N,CT,CL,CO,CPL,CPR,NC)
      implicit none
      
      integer N,NC,k,NB
      integer*2 CT(N)
      real*8 CL(N),CO(N),CPL(N-1),CPR(N-1),DISCR
      
      NB=N/NC
      do 10 k=1,N
      if ((mod(k,NB).EQ.1).OR.(mod(k,NB).EQ.0).OR.(k.EQ.N)) then
      CPL(k)=0.D0
      CPR(k)=0.D0
      else
      DISCR=1-(CL(k-1)/CL(k)*DSIN(CO(k-1))**2)
      if(DISCR.lt.0) then
      CPL(k)=0.D0
      else
      CPL(k)=CL(k-1)/CL(k)*(DSIN(CO(k-1)))**2-
     *DCOS(CO(k-1))*DSQRT(DISCR)
      endif
      DISCR=1-(CL(k+1)/CL(k)*DSIN(CO(k+1))**2)
      if(DISCR.lt.0) then
      CPR(k)=0.D0
      else
      CPR(k)=CL(k+1)/CL(k)*(DSIN(CO(k+1)))**2-
     *DCOS(CO(k+1))*DSQRT(1-(CL(k+1)/CL(k)*DSIN(CO(k+1))**2))
      endif
      endif
      write(*,'(I8,I8,I8,F8.2,F8.2)') k,CT(k),mod(k,NB),CPL(k),CPR(k)
10    continue
      return
      end
      
      subroutine AllAngles(N,NB,Phi,X)
      implicit none
      integer N,NB,i
      real*8 Phi(N),X(3,N)
      
      Phi(1)=-2.D0
      Phi(N)=-2.D0
      do 10 i=2,N-1
       if (mod(i,NB).eq.0.or.mod(i,NB).eq.1) then
        Phi(i)=-2.D0
       else
        call angles(X(1,i-1),X(2,i-1),X(3,i-1),X(1,i),X(2,i),X(3,i),
     * X(1,i+1),X(2,i+1),X(3,i+1),Phi(i))
       endif
10    continue

      end
