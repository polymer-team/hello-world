      subroutine OutErrorEnergyNV(mono,J,N,R2,V,neig,ln)
      implicit none
      include 'var.common'
      integer N,mono,J,neig(N,N),ln(N),k
      real*8 R2,V(3,N)
      
      Write(*,*)'Error: EnergyNV: Bad distance'
      Write(*,'(A3,2I5,A5,I5,A12,F12.8)')'Pr',Mprnk,mono,' and ',J,
     *    ' distance',R2
      open(20,FILE='Err_V.dat')
      write(20,'(13E18.10)')V(1,J),V(2,J),V(3,J),V(1,mono),V(2,mono),
     *  V(3,mono),R2
      do 30 k=1,N
       write(20,'(3E14.6)')V(1,k),V(2,k),V(3,k)
30    continue
      close(20)
      open(80,FILE='ENTab.dat')
      do k=1,N
       write(80,'(66I4)')k,(neig(k,j),j=1,ln(k))
      enddo
      close(80)
      stop''
C      EKey=EKey+1
      return
      end
      
      subroutine OutIni(MTrue,N,CopType,StifKey,E1,E2,logfname,prrunf,
     *                  genmtmp,henmtmp)
      implicit none
      include 'var.common'
      
      integer N,CopType(N),StifKey(N),E1,E2,MTrue,i
      character(Len=256) :: logfname,prrunf,genmtmp,henmtmp

      
      if(KEY.gt.0) then
       open(31,FILE=logfname,position="append")
      else
       open(31,FILE=logfname)
      endif

      Write(31,*)'Process ',MPrnk
      Write(31,*)'Bounds: Emin ',Eminbnd,' Emax: ',Emaxbnd,' Key: ',
     *  WinKey
      Write(31,*)'Copolymer sequence: '
      Write(31,'(32I2)') (CopType(i),i=1,N)
      Write(31,*)'StifKeyArray:'
      Write(31,'(32I2)') (StifKey(i),i=1,N)

      Write(*,*)'Rank',MTrue,'Initial E1 ',E1,' E2 ',E2
      Write(31,*)'Initial E1 ',E2,' E2 ',E2

      if (KEY.eq.0) then
       select case(Wkey)
        case(0)
         write(31,*)'Constant Weights'
        case(1)
         write(31,*)'Linear Weights'
        case(2)
         write(31,*)'Quadratic Weights'
        case(3)
         write(31,*)'Exponential Weights'
        case(4)
         write(31,*)'Exponential Weights - two dimensional'
        case(5)
         write(31,*)'Quadratic Weights - two dimensional'
        case(6)
         write(31,*)'Linear Weights - two dimensional'

        case default
         stop 'Error. Incorrect weights key'
       end select
C       call wanglandau(0.D0,0.D0,-1,0.D0,WLRes,WLFact,WLAcRate,WLNBRate)
       open(30,FILE=genmtmp,position="rewind")
       open(32,FILE=henmtmp,position="rewind")
       close(30,STATUS="delete")
       close(32,STATUS="delete")

      else
       if(KEY.eq.1) then
        write(*,*)'Pr_',MTrue,' Continue from step: ',SAMCPST
        open(30,FILE=genmtmp,position="append")
        open(32,FILE=henmtmp,position="append")
        write(30,*)'*********************Continue from step: ',SAMCPST
        write(31,*)'*********************Continue from step: ',SAMCPST
        write(32,*)'*********************Continue from step: ',SAMCPST
       endif
       if(KEY.eq.2) then
        write(*,*)'Pr_',MTrue,'**Continue gam0=',SMgam0,'T0=',SMT0
        open(30,FILE=genmtmp,position="append")
        open(32,FILE=henmtmp,position="append")
        write(30,*)'**Continue gam0=',SMgam0,'T0=',SMT0,'WLfact',SMgam0
        write(31,*)'**Continue gam0=',SMgam0,'T0=',SMT0,'WLfact',SMgam0
        write(32,*)'**Continue gam0=',SMgam0,'T0=',SMT0,'WLfact',SMgam0
       endif
       if(KEY.eq.3) then
        write(*,*)'Pr_',MTrue,'**Umbrella sampling'
        open(30,FILE=genmtmp,position="append")
        open(32,FILE=henmtmp,position="append")
        write(30,*)'**Umbrella sampling**'
        write(31,*)'**Umbrella sampling**'
        write(32,*)'**Umbrella sampling**'
       endif
       if(KEY.eq.4) then
        write(31,*)'*****Productive run*******************************'
       endif
      endif

      if(MPsize.eq.1.or.MPrnk.eq.0) then
      write(*,'(A3,A15,A4,2A6,A4,4A6,A4,10A8)')'Rn','STEP',
     *'Er1','E1mx','E1mn','Er2','E2mx','E2mn','Eflmx','Eflmn',
     *'flt','AR_loc','T_loc,m','AR_ECB','T_ecb,m','AR_CB','T_cb,m',
     * 'AR_EB','T_eb,m','AR_Rpt','T_rpt,m'
      endif
      write(31,'(A15,A10,A4,2A8,A4,4A8,A4,10A8)')'STEP','SSMEAN','Er1',
     * 'E1max','E1min','Er2','E2max','E2min','Eflmx','Eflmn','flt',
     * 'AR_loc','T_loc,m','AR_ECB','T_ecb,m','AR_CB','T_cb,m',
     * 'AR_EB','T_eb,m','AR_Rpt','T_rpt,m'

      close(31)

      if(KEY.ge.4) then
       open(36,FILE=prrunf)
       write(36,'(A15,3A5,20A13,3A5)')'NSTEP','Eful','ENV','Estf',
     * 'SSMEAN','R2END','SSMEANC','K1','K2','Eta1','Eta2','Eta3','SSA',
     * 'K1A','K2A','Eta1A','Eta2A','Eta3A','SSB','K1B','K2B','Eta1B',
     * 'Eta2B','Eta3B','EAA','EAB','EBB'
      endif
      close(36)
      end
      
      subroutine OutPrintCurCnf(N,X,E1,E2,CopType,PrintKey)
      implicit none
      include 'var.common'
      integer N,E1,E2,CopType(N),nxc,PrintKey,i,k
      integer cnfe1(100),cnfe2(100),cnfkey(100),cnfdel(100)
      real*8 X(3,N)
      character(Len=256) :: confname,vrmlname

      data cnfe1,cnfe2,cnfkey,cnfdel,nxc/100*0,100*0,100*0,100*0,0/

      if(PrintKey.eq.1) then
      do 50 i=1,nxc
       if(E1.eq.cnfe1(i).and.E2.eq.cnfe2(i).and.cnfkey(i).ne.0)then
        if(cnfdel(i).eq.0) then
         write(vrmlname,'(A5,I0.4,A1,I0.4,2(A2,I0.2),A4)')'xtri_',E1,
     *        '_',E2,'n_',cnfkey(i),'r_',MPrnk,'.wrl'
         call ToVRML(N,N,X,SIG,CopType,vrmlname)
         write(confname,'(A5,I0.4,A1,I0.4,2(A2,I0.2),A4)')'xtri_',E1,
     *        '_',E2,'n_',cnfkey(i),'r_',MPrnk,'.dat'
C          call ContMatrix(N,X,RCUT,mtrxfname)
         open(51,FILE=confname)
         do 51 k=1,N
          write(51,'(3E18.10)')X(1,k),X(2,k),X(3,k)
51       continue
         close(51)
         cnfkey(i)=cnfkey(i)-1
         cnfdel(i)=5000
        else
         cnfdel(i)=cnfdel(i)-1
        endif
       endif
50    continue
      return
      endif
      
      if(PrintKey.eq.0) then
       open(50,FILE='inmc-econ.dat')
       read(50,'(I4)')nxc
       if(nxc.gt.100.and.MPrnk.eq.0) then
        write(*,*)'PrintCurCnf: Only 100 first confs will be printed'
        nxc=100
       endif
       read(50,*)
       do i=1,nxc
        read(50,'(3I6)')cnfe1(i),cnfe2(i),cnfkey(i)
       enddo
       close(50)
      endif
      end
      
      
      subroutine OutPrintGHTMP(henmtmp,genmtmp,WLRes,NSTEP,WLFact)
      use dynamicResize
      implicit none
      include 'var.common'
      integer i,k,WLRes
      integer (Kind=8) NSTEP
      real*8 WLFact
      character WLResC
      character(Len=256) :: henmtmp,genmtmp
      character(Len=1024) :: fmt,fmt2

      WLResC='n'
      if(WLRes.eq.2) WLResC='y'

      if(KEY.lt.4.or.MOD(NSTEP,50000000).eq.0) then
       open(32,FILE=henmtmp,position="append")
       write(32,'(I15,A18,A4)')NSTEP+SAMCPST,' steps. H(Env,Est)',
     *   WLResC
       write(32,'(6I8,E15.6,I8,E15.8)')Nbin,Emin,Emax,N2bin,E2min,
     *  E2max,WLFact,nWLHf,SumGam
       write(fmt,'(A4,I0,A4)')'(A8,',N2bin,'I18)'
       write(32,fmt)' ENV\Est',(E2min+k-1,k=1,N2bin)
       write(fmt2,'   (A4,I0,A4)')'(I8,',N2bin,'I18)'

       do i=1,Nbin
        write(32,fmt2) Emin+i-1,(WLH(i,k),k=1,N2bin)
       enddo
       close(32)
      endif

      if(KEY.lt.4) then
       open(30,FILE=genmtmp,position="append")
       write(30,'(I15,A18,A4)')NSTEP+SAMCPST,' steps. G(Env,Est)',
     *   WLResC
       write(30,'(6I8,E15.6,I8,E15.8)')Nbin,Emin,Emax,N2bin,E2min,
     *   E2max,WLFact,nWLHf,SumGam
       write(30,fmt)' ENV\Est',(E2min+k-1,k=1,N2bin)
       write(fmt,'(A4,I0,A6)')'(I8,',N2bin,'E18.8)'

       do i=1,Nbin
        write(30,fmt)Emin+i-1,(WLG(i,k),k=1,N2bin)
       enddo
       close(30)
      endif
      end
