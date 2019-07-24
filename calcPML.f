      Program locq6
      implicit none
      include 'var.common'
      integer N,i,k,j,m,l,Nmax
      parameter (Nmax=1000)
      real*8 X(3,Nmax),DX,DY,DZ,R,RQ6C,q6(Nmax),q4(Nmax)
      real*8 costet,sintet,cosphi,phi,yr6(13),yi6(13),syr6(13),syi6(13)
      real*8 nn,pi,q6mr(Nmax,13),q6mi(Nmax,13),d6(Nmax),nn6(Nmax)
      real*8 nnc6(Nmax),D6pair,tanphi,d6trh
      integer StateKey(Nmax)
      
      integer Ncryst
      
      PI=4.D0*DATAN(1.D0)

C      write(*,*)'RCUT'
C      read(*,*)RQ6C
C      write(*,*)'D treshold'
C      read(*,*)d6trh

C      X(1,1)=0.714178
C      X(2,1)=0.714178
C      X(3,1)=0.D0
C      X(1,2)=0.714178
C      X(2,2)=0.D0
C      X(3,2)=0.714178
C      X(1,3)=0.D0
C      X(2,3)=0.714178
C      X(3,3)=0.714178
C      X(1,4)=0.D0
C      X(2,4)=0.D0
C      X(3,4)=0.D0
C      write(*,'(F10.6)')-0.5
      X(2,2)=1.D0
      X(2,3)=0.5D0
      X(1,3)=sqrt(3.D0)/2.D0
      X(2,4)=-0.5D0
      X(1,4)=sqrt(3.D0)/2.D0
      X(2,5)=-1.D0
      X(2,6)=-0.5D0
      X(1,6)=-sqrt(3.D0)/2.D0
      X(2,7)=0.5D0
      X(1,7)=-sqrt(3.D0)/2.D0
      X(3,8)=sqrt(6.D0)/3.D0
      X(1,8)=-sqrt(3.d0)/3.D0
      X(3,9)=sqrt(6.D0)/3.D0
      X(1,9)=sqrt(3.d0)/6.D0
      X(2,9)=0.5
      X(3,10)=sqrt(6.D0)/3.D0
      X(1,10)=sqrt(3.d0)/6.D0
      X(2,10)=-0.5
      X(3,11)=-sqrt(6.D0)/3.D0
      X(1,11)=-sqrt(3.d0)/3.D0
      X(3,12)=-sqrt(6.D0)/3.D0
      X(1,12)=sqrt(3.d0)/6.D0
      X(2,12)=0.5
      X(3,13)=-sqrt(6.D0)/3.D0
      X(1,13)=sqrt(3.d0)/6.D0
      X(2,13)=-0.5
C      do i=1,N
C       read(15,'(3F10.6)')X(1,i),X(2,i),X(3,i)
C       write(*,'(I3,3F10.6)')i,X(1,i),X(2,i),X(3,i)
C      enddo
      write(*,*)'Number of monomers'
      read(*,*)N

      open(15,FILE='inp.dat')
      do i=1,N
       read(15,*)X(1,I),X(2,I),X(3,I)
C       write(*,'(I3,3F10.6)')i,X(1,i),X(2,i),X(3,i)
      enddo
      close(15)

C      RQ6C=1.2
      d6trh=0.5

C      do 120 m=1,10
C      RQ6C=1.2
C      do 130 l=1,10
      RQ6C=1.4
      write(*,*)'RCUT'
      read(*,*)RQ6C

      syr6=0.D0
      syi6=0.D0
      q6=0.D0
      q6mr(N,13)=0.D0
      q6mi(N,13)=0.D0
      d6=0.D0
      nn6=0.D0
      nnc6=0.D0
      Ncryst=0
      StateKey=0

      do 11 i=1,N
C      nn=0.D0
       syr6=0.D0
       syi6=0.D0

       do 12 k=1,N
        if(k.eq.i) cycle
        DX=X(1,k)-X(1,i)
        DY=X(2,k)-X(2,i)
        DZ=X(3,k)-X(3,i)
        R=sqrt(DX**2+DY**2+DZ**2)
        if(R.lt.RQ6C) then
C         if(i.eq.1)write(*,*)i,R
         costet=DZ/R
C         sintet=sqrt(1.D0-costet**2)
         if(Abs(costet).eq.1.D0)then
          phi=pi/2.D0
          phi=0.D0
C         cosphi=DX/(R*sintet)
         else
          if(DX.ne.0.D0) then
           tanphi=DY/DX
           phi=atan(tanphi)
           if(DX.lt.0.d0) then
            phi=phi+pi
           endif
          else
           if(DY.gt.0.D0)phi=pi/2.D0
           if(DY.lt.0.D0)phi=-pi/2.D0
          endif
         endif

C         if(cosphi.gt.1.D0)cosphi=1.D0
C         if(cosphi.lt.-1.D0)cosphi=-1.D0
C         phi=acos(cosphi)
         call CalcYLM(6,costet,phi,yr6,yi6)
C         if(i.ge.10) then
C          write(*,'(3F10.5)')DX,R,sintet
C          write(*,'(I2,F20.15,3F10.5)')i,costet,cosphi,yr6(1),yi6(1)
C         endif
         syr6=syr6+yr6
         syi6=syi6+yi6
C         if(i.eq.2) then
C          if(k.eq.6.or.k.eq.8)write(*,*)costet,phi
C          write(*,'(I3,13F8.4)')k,(yr6(j),j=1,13)
C         endif
C        nn=nn+1.D0
         nn6(i)=nn6(i)+1.D0
        endif
12     continue
       if(nn6(i).gt.0.5) then
        syr6=syr6/nn6(i)
        syi6=syi6/nn6(i)
C        if(i.eq.2) then
C         write(*,'(I3,13F8.4)')I,(syr6(j),j=1,13)
C         write(*,'(I3,13F8.4)')I,(syi6(j),j=1,13)
C        endif

        do 13 k=1,13
         q6(i)=q6(i)+syr6(k)**2+syi6(k)**2
13      continue
C       q6(i)=q6(i)*4.D0*pi/(13.D0)
        q6(i)=sqrt(q6(i))
        q6mr(i,1:13)=syr6/q6(i)
        q6mi(i,1:13)=syi6/q6(i)
       endif
C       write(*,'(I4,F5.1,2F18.8)')i,nn6(i),q6(i)
11    continue

      open(100,FILE='q6.dat')
      do 20 i=1,N
      do 21 k=1,N
       if(k.eq.i) cycle
       DX=X(1,i)-X(1,k)
       DY=X(2,i)-X(2,k)
       DZ=X(3,i)-X(3,k)
       R=sqrt(DX**2+DY**2+DZ**2)
C       if(i.eq.1)write(*,*)i,R
       if(R.lt.RQ6C) then
        D6pair=0.D0
        do j=1,13
         D6pair=D6pair+q6mr(i,j)*q6mr(k,j)+q6mi(i,j)*q6mi(k,j)
C         d6(i)=d6(i)+q6mr(i,j)*q6mr(k,j)+q6mi(i,j)*q6mi(k,j)
        enddo
        d6(i)=d6(i)+D6pair
        if(D6pair.gt.D6trh)nnc6(i)=nnc6(i)+1.D0
       endif
21    continue
      d6(i)=d6(i)/nn6(i)
C      write(*,*)i,nn6(i),d6(i)
      if(nn6(i).ge.5.and.nnc6(i).ge.nn6(i)-1) then
C     crystal
      StateKey(i)=1
      else
       if(nn6(i).le.4) then
        StateKey(i)=2
       endif
      endif
      write(100,'(3F10.5,I4)')q6(i),d6(i),nn6(i),StateKey(i)
      if(d6(i).gt.d6trh.and.nn6(i).ge.12) Ncryst=Ncryst+1
20    continue
      Write(*,'(2F10.3,I5)')RQ6C,d6trh,Ncryst
      close(100)
C      RQ6C=RQ6C+0.05
C130   continue
C      d6trh=d6trh+0.05
C120   continue

      end

      subroutine CalcYLM(l,cosx,phi,ylmr,ylmi)
      implicit none
C      real pi
C      parameter (pi=3.1415926)

      integer l,m,i
      integer:: factor
      real*8 pln,cosx,ylmr(l*2+1),ylmi(l*2+1),ylm4r(9),ylm4i(9),phi,sign
      real*8 pi

      PI=4.D0*DATAN(1.D0)
      
      ylmr=0.D0
      ylmi=0.D0
      ylm4r=0.D0
      ylm4i=0.D0

C      l=6
      do 11 m=0,l
       call plgndr(l,m,cosx,pln)
C       write(*,*)factor(m)
       ylmr(l+m+1)=pln*dcos(dble(m)*phi)*dsqrt(dble(2*l+1)/(4.D0*pi)*
     *              dble(factor(l-m))/dble(factor(l+m)))
       ylmi(l+m+1)=pln*dsin(dble(m)*phi)*dsqrt(dble(2*l+1)/(4.D0*pi)*
     *              dble(factor(l-m))/dble(factor(l+m)))
       if(m.ne.0) then
        sign=1.D0
        if(Mod(m,2).ne.0)sign=-1.D0
        ylmr(l-m+1)=ylmr(l+m+1)*sign
        ylmi(l-m+1)=-ylmi(l+m+1)*sign
       endif
11    continue
C      do 12 i=1,13
C       write(*,'(I4,3F18.8)')i-7,ylm6r(i),ylm6i(i),dsin(dble(i-7)*phi)
C12    continue
      end
      
      ReCursive FUNCTION Factor(n)  RESULT(Fact)

      IMPLICIT NONE
      INTEGER :: Fact
      INTEGER, INTENT(IN) :: n

      IF (n.eq.0) THEN
       Fact = 1
      ELSE
       Fact = n * Factor(n-1)
      END IF
      END FUNCTION Factor

      subroutine plgndr(l,m,x,pln)
      
      implicit none
      integer l,m
      real*8 x,pln
      integer i,ll
      real*8 fact,pll,pmm,pmmp1,somx2
      
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.D0) stop'bad arg'
      pmm=1.D0
      if(m.gt.0) then
       somx2=dsqrt((1.D0-x)*(1.D0+x))
       fact=1.D0
       do i=1,m
        pmm=-pmm*fact*somx2
        fact=fact+2.D0
       enddo
      endif
      if(l.eq.m) then
C       plgndr=pmm
       pln=pmm
      else
       pmmp1=x*(2*m+1)*pmm
       if(l.eq.m+1) then
C        plgndr=pmmp1
        pln=pmmp1
       else
        do ll=m+2,l
         pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
         pmm=pmmp1
         pmmp1=pll
        enddo
C        plgndr=pll
        pln=pll
       endif
      endif
      return
      end
      

