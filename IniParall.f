      subroutine IniParall(N,NB,X,rand0,CopType,StifKey,CopPhi,E1,E2,
     *                     WLFact,neighb,lastneig,MPKey,parallf)
      use DynamicResize
      implicit none
      include 'var.common'

      integer (kind=8)NSTEP

      integer N,NB,CopType(N),StifKey(N),WLRes,IKey,WinKeytmp,MPKey
      integer E1,E2,Ecur,edum1,edum2,edum3,edum4
      integer i,k,neighb(N,N),lastneig(N),CKey
      integer EofKey

      real*8 X(3,N),V(3,N),rand0,Phi(N),CopPhi(2,N),WLFact
      real*8 rndctrl(100)

      character(Len=256) :: vrmlname,inxf,parallf,xconfname,rep

      if(MPrnk.eq.0) then
       open(45,file='inmc-whlp.dat')
       write(45,*)'Input conformations to next window:'
       write(45,'(8A8)')'Window','E2','E1mx','|','Emax','Emin',
     *   'E2max','E2min'
      endif

      if(KEY.eq.0) then
       open(10,FILE=parallf)
       read(10,'(I3)')nWin
       MPKey=1
       if(nWin.le.MPrnk) return

       MPKey=0

       read(10,*)
       do 15 i=1,nWin
        if(MPrnk.eq.i-1) then
C         read(10,'(5I8)')Emaxbnd,Eminbnd,WinKeytmp,Wkey,CKey
         if(i.ne.1) then
          ngwinbnd(1,1)=edum2
          ngwinbnd(1,2)=1
         endif
         read(10,'(7I8)')edum1,edum2,edum3,edum4,WinKeytmp,Wkey,CKey
         if(i.lt.nWin) then
          read(10,'(7I8)')EDeepbnd,k,k,k,k,k,CDeepKey
          ngwinbnd(2,1)=EDeepbnd
          ngwinbnd(2,2)=1
C Men9t zdes!!! $tobi na4alnie konformazii peredavat po jestkosti
          backspace 10
         endif
C         exit
         Emaxbnd=edum1
         Eminbnd=edum2
         E2maxbnd=edum3
         E2minbnd=edum4
         write(*,*)Mprnk,Emaxbnd,Eminbnd,E2maxbnd,E2minbnd
         if(MPrnk.ne.0)exit
        else
         read(10,'(4I8)')edum1,edum2,edum3,edum4
        endif
        if(MPrnk.eq.0)write(45,'(3I8,A8,4I8)')i-1,0,0,'|',edum1,edum2,
     *  edum3,edum4
15     continue
        if(MPrnk.eq.0)close(45)

       if(nWin.gt.MPsize) then
        nWin=MPSize
        if(MPrnk.eq.0) write(*,*)'Maximum number of windows: ',nWin
        if(MPrnk.eq.MPsize-1)WinKey=-1
       endif
       close(10)

       call BEGIN(N,NB,X,rand0,CopType,Phi,StifKey,CopPhi)

       write(xconfname,'(A6,I0.2,A4)')'x_cnfi',MPrnk,'.dat'

       if(CKey.ne.0) then
        EofKey=0
        open(10,FILE=xconfname,Status='old',IOStat=EofKey)
        if(EofKey.ne.0) goto 500
        do 30 i=1,N
         read(10,*)X(1,I),X(2,I),X(3,I)
30      continue
        close(10)
        WinKey=WinKeytmp
C        close(60)
        goto 600
       endif
       

500    call E1Full(N,X,E1,1,CopType)
       call E2Full(N,X,E2,Phi,CopPhi,CopType,StifKey,1)

       Ecur=E1+E2

       WinKey=2
       IniKey=1
       call SAMC(E1,0,E2,0,-1,0.D0,WLRes,0.D0,0,0,0,N)
C       Write(*,*)MPrnk,'SAMC returned'



       E1cmax=E1
       E1cmin=E1
       E2cmax=E2
       E2cmin=E2
       Efcmax=E1cmax+E2cmax
       Efcmin=E1cmin+E2cmin

       call cpu_time(tpr1)
       IKey=0
       NSTEP=1
       write(inxf,'(A5,I0.2,A4)')'x_cnf',MPrnk,'.dat'
       write(vrmlname,'(A5,I0.2,A4)')'x_cnf',MPrnk,'.wrl'
C       Write(*,*)MPrnk,'Starting step'
       do while (IKey.EQ.0)
        call Step(N,X,V,NSTEP,neighb,lastneig,CopType,E1,WLFact,Phi,
     *            CopPhi,StifKey,E2,rndctrl)
        Ecur=E1+E2
        if(E1.le.Emaxbnd-1.and.E1.ge.Eminbnd+1) then
        if(E2.le.E2maxbnd.and.E2.ge.E2minbnd) then
         open(10,FILE=inxf)
         do 20, k=1,N
          write(10,'(3E18.10)')X(1,k),X(2,k),X(3,k)
20       continue
         close(10)
         call ToVRML(N,NB,X,SIG,CopType,vrmlname)
         IKey=1
        endif
        endif
        NSTEP=NSTEP+1
        
        if(MOD(NSTEP,10000).eq.0.and.CKey.eq.2) then
         EofKey=0
         open(10,FILE=xconfname,Status='old',IOStat=EofKey)
         if(EofKey.ne.0) cycle
         do 31 i=1,N
          read(10,*)X(1,I),X(2,I),X(3,I)
31       continue
         close(10)
         WinKey=WinKeytmp
         exit
        endif
       end do
C       open(60,FILE=rep,position="append")
C       Write(60,*)MPrnk,'Initial steps ended'

       call cpu_time(tpr2)
       Write(*,*)'Pr ',MPrnk,' Prob_Time: ',tpr2-tpr1

       deallocate(WLH)
       deallocate(WLG)
       deallocate(PW)
       deallocate(Mask)
       
C       write(*,*)'Ini_priv'
       if(allocated(gemod)) deallocate(gemod)
       if(allocated(stan)) deallocate(stan)
C       write(*,*)'Ini_pokedava'
       WLFact=2.D0
       WinKey=WinKeytmp
       IniKey=0

600    return
C       close(60)
      endif

      if(KEY.gt.0) then
       MPKey=1
       open(30,FILE=parallf)
       read(30,'(I3)')nWin
C       read(10,'(2I8)')Eglmn,Eglmx
       read(30,*)
       if(nWin.le.MPrnk) return
       MPKey=0
       i=0
       do while (MPrnk.gt.i)
        read(30,*)
        i=i+1
       enddo
       read(30,'(6I8)')Emaxbnd,Eminbnd,E2maxbnd,E2minbnd,WinKey,Wkey
       call BEGIN(N,NB,X,rand0,CopType,Phi,StifKey,CopPhi)
       if(MPrnk.eq.0) then
        write(45,'(3I8,A8,4I8)')0,0,0,'|',Emaxbnd,Eminbnd,E2maxbnd,
     *   E2minbnd
C        write(*,'(3I8,A8,2I8)')0,0,0,'|',Emaxbnd,Eminbnd
C        backspace 30
        do while (nWin-1.gt.i)
         read(30,'(4I8)')edum1,edum2,edum3,edum4
         i=i+1
         write(45,'(3I8,A8,4I8)')i,0,0,'|',edum1,edum2,edum3,edum4
C         write(*,'(3I8,A8,2I8)')i,0,0,'|',edum1,edum2
        enddo
        close(45)
       endif
       close(30)
      endif
      endln
