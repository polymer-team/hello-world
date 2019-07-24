      subroutine EndBridging(N,X,V,neighb,lastneig,Phi,CopPhi,ENV,EST,
     *                       WLFact,NSTEP,CopType)
      implicit none

      include 'var.common'
      include 'ranmar.common'

      integer (kind=8)NSTEP
      integer N,WLResult,NB,K,I,m,J,oe,e1,e2,ch,se,lp,rp,CopType(N)
      integer ENV,DENV,EST,DEST,Efull,DEfull,ENVn,ESTn
      integer neighb(N,N),lastneig(N),neigtmp(N,N),lntmp(N),tmp(N)
      integer BridgeNBH(N),BridgeCNT
      integer stENV,stEST
      integer FINDNBH


      real*8 X(3,N),XI,YI,ZI,V(3,N),XTMP(3,N),VTMP(3,N),R,R1
      real*8 Phi(N),CopPhi(2,N),COS1,COS2,COS3,PhiTmp(N)
      real*8 WLFact

      real rand(2),rnd4,rnd5

      character(Len=256) :: vrmlname,rep


      BridgeCNT=0
      BridgeNBH=0
      DEST=0
      DENV=0
      XTMP=0.D0
      call ranmar(rand,2)
      rnd4=rand(1)-1.0D-6
      rnd5=rand(2)

      call ChooseEndEB(N,X,lastneig,neighb,BridgeCNT,BridgeNBH,m)

200   if(BridgeCNT.eq.0) return

      i=int(real(BridgeCNT)*rnd4)+1
      if(i.gt.BridgeCNT) write(*,*)MPrnk,'Error EB! i>CNT'

      J=BridgeNBH(i)
      if(i.le.m)then
C     left end
       oe=1
       e1=1
       e2=0
       ch=1
       se=N
       lp=1
       rp=J
      else
C     right end
       oe=N
       e1=0
       e2=1
       ch=-1
       se=1
       lp=J
       rp=N
      endif

      call EnergyEB(N,X,DENV,DEST,Phi,CopPhi,j,oe,ch,COS1,COS2,CopType)

      ESTn=EST+DEST
      ENVn=ENV+DENV
      call SAMC(ENV,ENVn,EST,ESTn,0,rnd5,WLResult,WLFact,cntEB,cntNEB,
     *          NSTEP+SAMCPST,N)

C      stENV=40-(Emaxbnd-ENV)
C      stEst=EST+1

      if(WLResult.EQ.1) then
C       stRAR(2,1,stENV,stEst)=stRAR(2,1,stENV,stEst)+1
C       stRAR(2,3,stENV,stEst)=stRAR(2,3,stENV,stEst)+DEST

       neigtmp=neighb
       lntmp=lastneig
       do 102 i=lp+e2,rp-e1
        XTMP(1:3,i)=X(1:3,rp-i+lp+e2-e1)
        VTMP(1:3,i)=V(1:3,rp-i+lp+e2-e1)
        PhiTMP(i)=Phi(rp-i+lp+e2-e1)
        lntmp(i)=lastneig(rp-i+lp+e2-e1)
C        neigtmp(i,1:lntmp(i))=neighb(rp-i+lp+e2-e1,1:lntmp(i))
        neigtmp(i,1:N)=neighb(rp-i+lp+e2-e1,1:N)
102    continue

       X(1:3,lp+e2:rp-e1)=XTMP(1:3,lp+e2:rp-e1)
       V(1:3,lp+e2:rp-e1)=VTMP(1:3,lp+e2:rp-e1)
       PhiTMP(oe)=Phi(oe)
       PhiTMP(J-ch)=COS1
       PhiTMP(J)=COS2
       Phi(lp:rp)=PhiTMP(lp:rp)
       neighb(lp:rp,1:N)=neigtmp(lp:rp,1:N)
       lastneig(lp:rp)=lntmp(lp:rp)

C       if(NSTEP.eq.39) then
C        open(80,FILE='EBErr2.dat')
C        do i=1,N
C         write(80,'(66I4)')i,lastneig(i),(neighb(i,k),k=1,20)
C        enddo
C        close(80)
C        write(*,*)J,oe
C       endif

C       delete 1st/Nth from neighbour list of J
       if(oe.eq.1) neighb(J,1:lastneig(J)-1)=neighb(J,2:lastneig(J))
       neighb(J,lastneig(J))=0
       lastneig(J)=lastneig(J)-1

C       delete Jst from neighbour list of 1st/Nth (now (J-1)st/(J+1)st)
       do i=1,lastneig(J-ch)
        if(neighb(J-ch,i).eq.J) exit
       enddo
       neighb(J-ch,i:lastneig(J-ch)-1)=neighb(J-ch,i+1:lastneig(J-ch))
       neighb(J-ch,lastneig(J-ch))=0
       lastneig(J-ch)=lastneig(J-ch)-1

C       if(NSTEP.eq.39) then
C        open(80,FILE='EBErr3.dat')
C        do i=1,N
C         write(80,'(66I4)')i,(neighb(i,k),k=1,20)
C        enddo
C        close(80)
C        write(*,*)J,oe
C       endif

C       renumerate elements
       if(oe.eq.1) then
        do 103 i=1,N
         lntmp=0
         do 104 k=1,lastneig(i)
          if(neighb(i,k).ge.J) exit
          lntmp(k)=J-neighb(i,k)
104      continue
         do 105 m=1,k-1
          neighb(i,m)=lntmp(k-m)
105      continue
103     continue
       else
        do 106 i=1,N
         lntmp=0
         do 107 k=lastneig(i),1,-1
          if(neighb(i,k).le.J) exit
          lntmp(k)=J+1+N-neighb(i,k)
107      continue
         do 108 m=lastneig(i),k+1,-1
          neighb(i,m)=lntmp(k+1+(lastneig(i)-m))
108      continue
106     continue
       endif

C       if(NSTEP.eq.75435) then
C        open(80,FILE='EBErr4.dat')
C        do i=1,N
C         write(80,'(66I4)')i,(neighb(i,k),k=1,lastneig(i))
C        enddo
C        close(80)
C        write(*,*)J,oe
C       endif

C       add J to neigbour list of new 1st/Nth
       do i=1,lastneig(oe)
        if(neighb(oe,i).gt.J) exit
       enddo
       lntmp(1:lastneig(oe)-i+1)=neighb(oe,i:lastneig(oe))
       neighb(oe,i+1:lastneig(oe)+1)=lntmp(1:lastneig(oe)-i+1)
       neighb(oe,i)=J
       lastneig(oe)=lastneig(oe)+1

       if(oe.eq.1)then
        lntmp(1:lastneig(J))=neighb(J,1:lastneig(J))
        neighb(J,2:lastneig(J)+1)=lntmp(1:lastneig(J))
        neighb(J,1)=1
       else
        neighb(J,lastneig(J)+1)=N
       endif
       lastneig(J)=lastneig(J)+1
       EST=ESTn
       ENV=ENVn

C       if(NSTEP.eq.22861) then
C        open(80,FILE='EBErr5.dat')
C        do i=1,N
C         write(80,'(66I4)')i,(neighb(i,k),k=1,20)
C        enddo
C        close(80)
C        write(*,*)J,oe
C       endif

C       check if number of neighbours is greater than 1
       XI=X(1,J+1)-X(1,J-1)
       YI=X(2,J+1)-X(2,J-1)
       ZI=X(3,J+1)-X(3,J-1)
       R=dsqrt(XI*XI+YI*YI+ZI*ZI)/BOND
       if(R.le.RMAX)then
C      podumat, mojet nado -ch
        if(lastneig(J+ch).gt.0) then
        tmp(1:lastneig(J+ch))=neighb(J+ch,1:lastneig(J+ch))
        k=FINDNBH(lastneig(J+ch),tmp,J-ch)
        if(neighb(J+ch,k).ne.J-ch) then
         tmp(1:lastneig(J+ch)-k+1)=neighb(J+ch,k:lastneig(J+ch))
         neighb(J+ch,k)=J-ch
         neighb(J+ch,k+1:lastneig(J+ch)+1)=tmp(1:lastneig(J+ch)-k+1)
         lastneig(J+ch)=lastneig(J+ch)+1
         if(lastneig(J-ch).gt.0) then
          tmp(1:lastneig(J-ch))=neighb(J-ch,1:lastneig(J-ch))
          k=FINDNBH(lastneig(J-ch),tmp,J+ch)
          tmp(1:lastneig(J-ch)-k+1)=neighb(J-ch,k:lastneig(J-ch))
          neighb(J-ch,k)=J+ch
          neighb(J-ch,k+1:lastneig(J-ch)+1)=tmp(1:lastneig(J-ch)-k+1)
          lastneig(J-ch)=lastneig(J-ch)+1
         else
          lastneig(J-ch)=1
          neighb(J-ch,1)=J+ch
         endif
         endif
C         write(*,'(20I3)')J+ch,(neighb(J+ch,k),k=1,lastneig(J+ch))
C         write(*,'(20I3)')J-ch,(neighb(J-ch,k),k=1,lastneig(J-ch))
C         write(*,*)''
        endif
       endif

C       if(NSTEP.eq.22861) then
C        open(80,FILE='EBErr6.dat')
C        do i=1,N
C         write(80,'(66I4)')i,(neighb(i,k),k=1,20)
C        enddo
C        close(80)
C        write(*,*)J,oe
C       endif


       i=J-2*ch
C       if(NSTEP.eq.75435) write(*,*)'i=',i

       XI=X(1,J)-X(1,I)
       YI=X(2,J)-X(2,I)
       ZI=X(3,J)-X(3,I)
       R=dsqrt(XI*XI+YI*YI+ZI*ZI)/BOND
       if(R.le.RMAX)then
C      podumat, mojet nado -ch
        tmp(1:lastneig(J))=neighb(J,1:lastneig(J))
        k=FINDNBH(lastneig(J),tmp,I)
        if(neighb(J,k).ne.I) then
         tmp(1:lastneig(J)-k+1)=neighb(J,k:lastneig(J))
         neighb(J,k)=I
         neighb(J,k+1:lastneig(J)+1)=tmp(1:lastneig(J)-k+1)
         lastneig(J)=lastneig(J)+1
         if(lastneig(I).gt.0) then
          tmp(1:lastneig(I))=neighb(I,1:lastneig(I))
          k=FINDNBH(lastneig(I),tmp,J)
          tmp(1:lastneig(I)-k+1)=neighb(I,k:lastneig(I))
          neighb(I,k)=J
          neighb(I,k+1:lastneig(I)+1)=tmp(1:lastneig(I)-k+1)
          lastneig(I)=lastneig(I)+1
         else
          lastneig(I)=1
          neighb(I,1)=J
         endif
C         write(*,'(20I3)')J,(neighb(J,k),k=1,lastneig(J))
C         write(*,'(20I3)')I,(neighb(I,k),k=1,lastneig(I))
C         write(*,*)''
        endif
       endif
       return
      endif
C      stRAR(2,4,stENV,stEst)=stRAR(2,4,stENV,stEst)+1

      end
      
      INTEGER FUNCTION FINDNBH(inr,arr,i)
      implicit none
      integer lb,inr,rb,mid,mid2,arr(inr),i,k
C naiti libo index i libo index pervogo 4to bolwe i
      lb=1
      mid=1
      mid2=0
      rb=inr
C      write(*,'(A4,20I4)')'ent:',rb,i,(arr(k),k=1,rb)
      do while(lb.lt.rb)
       mid=(lb+rb)/2
       if(mid.eq.mid2.and.mid.lt.rb)mid=mid+1
       mid2=mid
       if(arr(mid).lt.i)then
        lb=mid
       else
        if(arr(mid).eq.i) goto 50
        if(rb.eq.mid) goto 50
        rb=mid
       endif
      enddo
      if(arr(mid).gt.i) goto 50
      if(arr(mid).lt.i) mid=mid+1
50      FINDNBH=mid
C      write(*,*)FINDNBH
      RETURN
      END


