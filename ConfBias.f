      subroutine CCB(N,X,V,neighb,lastneig,Phi,CopPhi,E1,E2,
     *                WLFact,NSTEP,CopType)
      implicit none
      include 'ranmar.common'
      include 'var.common'

      integer N,NB,i,k,ECKey,E1,E2,E1n,E2n,WLResult,CopType(N)
      integer pos1,pos2,tmp,DE1,DE2
      integer neighb(N,N),lastneig(N),neigtmp(N,N),lntmp(N)
      integer (kind=8) NSTEP

      real*8 X(3,N),V(3,N),Phi(N),CopPhi(2,N),XTMP(3,N),VTMP(3,N)
      real*8 WLFact,Phitmp(N)

      real rand(3),rnd1,rnd2,rnd3

      character(Len=256) :: vrmlname

      ECKey=0
      NB=N/NC
      call ranmar(rand,3)
      rnd1=rand(1)
      rnd2=rand(2)
      rnd3=rand(3)
      pos1=int(real(NB-4)*rnd1)+1+2
      call POS2CCB(N,rnd2,POS2)
C      if(rnd2.lt.0.5) then
C       pos2=1
C      else
C       pos2=N
C      endif

      if(pos1.ge.N-1.or.pos1.le.2) then
C       write(*,*)'Err: CB, random gen: pos1 is not within interval',pos1
       return
      endif

      XTMP=X
      VTMP=V
      Phitmp=Phi
      lntmp=lastneig
      neigtmp=neighb

      call EnergyCCB(N,V,lastneig,neighb,DE1,DE2,POS1,POS2,Phi,CopPhi,
     * CopType)

      E1n=E1-DE1
      E2n=E2-DE2

      call BUILDCCB(N,XTMP,VTMP,POS1,POS2,ECKey,neigtmp,lntmp,
     *              Phitmp,CopPhi,NSTEP)

      call EnergyCCB(N,VTMP,lntmp,neigtmp,DE1,DE2,POS1,POS2,Phitmp,
     *               CopPhi,CopType)

      E1n=E1n+DE1
      E2n=E2n+DE2
      
      call SAMC(E1,E1n,E2,E2n,ECKey,rnd3,WLResult,WLFact,cntCB,
     *          cntNCB,NSTEP+SAMCPST,N)


      if(WLResult.EQ.1) then

       X=XTMP
       V=VTMP
       Phi=Phitmp
       E1=E1n
       E2=E2n
       neighb=neigtmp
       lastneig=lntmp
       return
      endif

      return
      end


      SUBROUTINE BUILDCCB(N,X,V,POS1,POS2,res,neig,ln,
     *                    Phi,CopPhi,NSTEP)
      implicit none
      include 'var.common'
      include 'ranmar.common'

      real pi
      parameter (pi=3.1415926)

      integer(kind=8)NSTEP
      integer N,I,J,k,POS1,POS2,ind,res,m,ch,Ekey
      integer ln(N),neig(N,N),lntmp(N),neigtmp(N,N),tmp(N),lput(N)
      integer FINDNBH
      
      Real*8 X(3,N),V(3,N),R1,R2,R3,R,XI,YI,ZI,DX,DY,DZ,Phi(N)
      Real*8 xs,ys,zs,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,l,f1,f2,f3
      Real*8 CopPhi(2,N),VXI,VYI,VZI
      Real*8 atet,aphi,ro,COS1,costet,atetmax,ccos,csin
      real rand(N*3),rnd1,rnd2,rnd3,rnd4,rnd5

      res=1
      call ranmar(rand,3*N)
      ind=0
      lntmp=ln
      neigtmp=neig
C      lput=-1

      if(pos2.eq.N)then
       ch=1
       do 40 i=1,pos1-1
        if(ln(i).ne.0) then
         tmp(1:ln(i))=neig(i,1:ln(i))
         k=FINDNBH(ln(i),tmp,pos1)
         neigtmp(i,k:ln(i))=0
         lntmp(i)=k-1
C        else
C         lput(i)=1
        endif
C        if(ln(i).eq.0.and.lntmp(i).eq.1)write(*,*)'Durdulet'
40     continue
      else
       ch=-1
       do 41 i=N,pos1+1,-1
        if(ln(i).ne.0) then
         tmp(1:ln(i))=neig(i,1:ln(i))
         k=FINDNBH(ln(i),tmp,pos1)
         if(neig(i,k).ne.pos1)k=k-1
         if(ln(i)-k.gt.0) then
          tmp(1:ln(i)-k)=neig(i,k+1:ln(i))
          neigtmp(i,1:ln(i)-k)=tmp(1:ln(i)-k)
          neigtmp(i,ln(i)-k+1:ln(i))=0
          lntmp(i)=ln(i)-k
         else
          lntmp(i)=0
          neig(i,1:ln(i))=0
         endif
C        new>>>>
         neigtmp(i,lntmp(i)+1:N)=0
C        <<<<< new
        endif
41     continue
      endif

      DO 20 I=POS1,POS2,ch
C      old energy
       J=I-ch
       if(ln(i).ne.0) neigtmp(i,1:ln(i))=0
       lntmp(i)=0

       DX=X(1,J)-X(1,J-ch)
       DY=X(2,J)-X(2,J-ch)
       DZ=X(3,J)-X(3,J-ch)
       R=DX**2+DY**2+DZ**2
       R=dSQRT(R)
       c1=DX/R
       c2=DY/R
       c3=DZ/R
C      building normal to O'z' = A
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

C       aphi=2.d0*PI*rand(ind+1)
C       ro=BndMin+rand(ind+2)*(BndMax-BndMin)
C       costet=-(SIG**2-R**2-ro**2)/(2.d0*R*ro)
C       atetmax=acos(costet)
C       atet=(PI-atetmax)*rand(ind+3)
C       atet=PI*rand(ind+1)
       ro=BndMin3+rand(ind+1)*(BndMax3-BndMin3)
       ro=ro**(1.D0/3.D0)
       ccos=-1.D0+rand(ind+2)*(1.D0+Mxctet)
       csin=dsqrt(1.D0-ccos**2)
       aphi=2.d0*PI*rand(ind+3)

C       xs=ro*SIN(atet)*COS(aphi)
C       ys=ro*SIN(atet)*SIN(aphi)
C       zs=ro*COS(atet)
       xs=ro*csin*DCOS(aphi)
       ys=ro*csin*DSIN(aphi)
       zs=-ro*ccos

       XI=X(1,J)+xs*a1+ys*b1+zs*c1
       YI=X(2,J)+xs*a2+ys*b2+zs*c2
       ZI=X(3,J)+xs*a3+ys*b3+zs*c3
       X(1,I)=XI
       X(2,I)=YI
       X(3,I)=ZI

       call CheckGeometry(XI,YI,ZI,EKey)
       if(Ekey.ne.0) return

       call CalcVI(XI,YI,ZI,VXI,VYI,VZI)

       V(1,I)=VXI
       V(2,I)=VYI
       V(3,I)=VZI

       if(pos2.eq.N) then
       k=1
       do 100 while (k.le.J-1)
        DX=V(1,k)-V(1,I)
        DY=V(2,k)-V(2,I)
        DZ=V(3,k)-V(3,I)
        call CalcDV(R1,R2,R3,DX,DY,DZ)
        R=R1*R1+R2*R2+R3*R3
        R=DSQRT(R)

        if(R.le.RMAX)then
C         if(R.LE.RCUT) then
C          DENV=DENV-EPS
          if (R.LE.SIG) then
C          cntCBer=cntCBer+1
C          if(k.eq.J-1)cntCBer2=cntCBer2+1
          return
          endif
C         endif
         lntmp(i)=lntmp(i)+1
         neigtmp(i,lntmp(i))=k
         lntmp(k)=lntmp(k)+1
         neigtmp(k,lntmp(k))=i
        else
         k=k+(R-RMAX)/BndMax
        endif
        k=k+1
100    continue
       else
       k=N
       do 101 while (k.ge.J+1)
C       if(NSTEP.eq.3825398.and.MPrnk.eq.2)write(*,*)'brk4',i,k,lntmp(i),
C     *lntmp(k)
        DX=V(1,k)-V(1,I)
        DY=V(2,k)-V(2,I)
        DZ=V(3,k)-V(3,I)
        call CalcDV(R1,R2,R3,DX,DY,DZ)
        R=R1*R1+R2*R2+R3*R3
        R=DSQRT(R)
        if(R.le.RMAX)then
C         if(R.LE.RCUT) then
C          DENV=DENV-EPS
          if (R.LE.SIG) then
C          cntCBer=cntCBer+1
C          if(k.eq.J+1)cntCBer2=cntCBer2+1
          return
          endif
C         endif
         if(lntmp(i).ne.0) then
          tmp(1:lntmp(i))=neigtmp(i,1:lntmp(i))
          neigtmp(i,2:lntmp(i)+1)=tmp(1:lntmp(i))
          lntmp(i)=lntmp(i)+1
          neigtmp(i,1)=k
         else
          lntmp(i)=1
          neigtmp(i,1)=k
         endif
C         lntmp(i)=lntmp(i)+1
C         neigtmp(i,lntmp(i))=k
C       if(lntmp(k).eq.0) write(*,*)'lntmp(k) nyjna proverka!'
         if(lntmp(k).ne.0) then
          tmp(1:lntmp(k))=neigtmp(k,1:lntmp(k))
          neigtmp(k,2:lntmp(k)+1)=tmp(1:lntmp(k))
          lntmp(k)=lntmp(k)+1
          neigtmp(k,1)=i
         else
          lntmp(k)=1
          neigtmp(k,1)=i
         endif
        else
         k=k-int((R-RMAX)/BndMax)
        endif
        k=k-1
101    continue
C      if(NSTEP.eq.3825398.and.MPrnk.eq.2)write(*,*)'brk5'

       endif
       call ANGLES(X(1,J-ch),X(2,J-ch),X(3,J-ch),X(1,J),X(2,J),X(3,J),
     *   X(1,I),X(2,I),X(3,I),COS1)
C        if(COS1.gt.CopPhi(1,J).or.COS1.lt.CopPhi(2,J)) then
C         DEST=DEST+EPSST
C        endif
        Phi(J)=COS1
       ind=ind+3
20    CONTINUE
C      do 71 i=1,N
C       write(*,'(20I4)')i,lntmp(i),(neigtmp(i,k),k=1,lntmp(i))
C71    continue
      neig=neigtmp
      ln=lntmp
      res=0
      end       
