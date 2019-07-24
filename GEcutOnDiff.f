      Program GECutOnDiff
C     finds values of dlng(E) higher than specified and cuts
C      use DynamicResize
       implicit none
        integer Nbin,N2bin,N1binmax
      integer Nbinnew,N2binnew
      integer(kind=8),dimension(:,:),allocatable::WLH,WLHtmp,PrWLH
       integer(kind=8),dimension(:,:),allocatable::matrxtmpH
      real*8,dimension(:,:),allocatable::WLG,WLGtmp,matrxtmpG
      real*8,dimension(:,:),allocatable::PW,PWtmp,Mask,MaskTmp
      real*8,dimension(:,:,:),allocatable::gemod,gemodtmp
      integer(kind=8),dimension(:,:,:,:),allocatable::StAn,StAntmp
      integer,dimension(:),allocatable:: E1max,E1min,N1bin


      integer i,k,ENcut,j,ENbin,n
      character(Len=1024) :: fmt,fmt2
      real*8 GEmin,dGEmin,WLFact,SumGam
      integer Emin,Emax,E2min,E2max,nWLHf
      integer (kind=8)SAMCPST
C      include 'var.common'
      
      open(40,FILE='mc-gnm_N64_Nc002_b04_aa1_bb1_ab1_00.dat')
      open(41,FILE='mc-hnm_N64_Nc002_b04_aa1_bb1_ab1_00.dat')
      read(40,*)SAMCPST
C      read(40,'(6I8)')Nbin,Emin,Emax,N2bin,E2min,E2max
      read(40,'(6I8,E15.6,I8,E15.8)')Nbin,Emin,Emax,N2bin,E2min,E2max,
     * WLFact,nWLHf,SumGam

      read(40,*)
      read(41,*)
      read(41,*)
      read(41,*)
      allocate(WLG(Nbin,N2bin))
      allocate(WLH(Nbin,N2bin))
      allocate(Mask(Nbin,N2bin))
      allocate(matrxtmpG(Nbin*N2bin,3))
       allocate(matrxtmpH(Nbin*N2bin,3))
      WLG=0.d0
      Mask=0.d0
      WLH=0
      matrxtmpG=0.d0
      matrxtmpH=0
      write(*,'(6I8)')Nbin,Emin,Emax,N2bin,E2min,E2max

        write(fmt2,'(A4,I0,A6)')'(I8,',N2bin,'E18.8)'
        write(fmt,'(A4,I0,A6)')'(I8,',N2bin,'I18)'
      
      do   i=1,Nbin*N2bin
       read(40,*)(matrxtmpG(i,k),k=1,3)
       read(41,*)(matrxtmpH(i,k),k=1,3)
      end do
      
      n=0
      do 11 i=1,Nbin
      n=n+1
       do 12 k=1,N2bin
      WLG(i,k)=matrxtmpG(k+(n-1)*N2bin,3)
      WLH(i,k)=matrxtmpH(k+(n-1)*N2bin,3)
       if(WLH(i,k).ne.0) Mask(i,k)=1.d0
12    continue
11    continue
      deallocate(matrxtmpG)
      deallocate(matrxtmpH)
C       read(40,fmt2)j,(WLG(i,k),k=1,N2bin)
C       read(41,fmt)j,(WLH(i,k),k=1,N2bin)
C       do 12 k=1,N2bin
C         if(WLH(i,k).ne.0) Mask(i,k)=1.d0
C12     continue
C11     continue
      close(40)
      close(41)

      write(*,*)'Enter maximum value of dlnG(E). Cut lower values of g'
      read(*,*)dGEmin

      GEmin=minval(WLG)
C                inverse order of k and i!!!!
      do k=1,N2bin
       do i=2,Nbin
        if(Mask(i,k).lt.0.5)cycle
        if(Mask(i-1,k).lt.0.5)cycle
        if(WLG(i,k)-WLG(i-1,k).gt.dGEmin) then
         if(GEmin.lt.WLG(i,k)-1.D0) GEmin=WLG(i,k)-1.D0
        endif
       enddo
      enddo

      
      ENcut=0
      do i=1,Nbin
      ENbin=0
      do k=1,N2bin
       if(WLG(i,k).lt.GEmin) then
        WLG(i,k)=0.D0
        Mask(i,k)=0.D0
        WLH(i,k)=0
       else
        if(Mask(i,k).gt.0.5)ENbin=1
       endif
      enddo
      if(ENbin.eq.0)ENcut=ENcut+1
      enddo
      Nbin=Nbin-ENcut
      Emin=Emin+ENcut
      write(*,*)'ENcut',ENcut
      write(*,*)'NbinNew',Nbin

      allocate(matrxtmpG(Nbin*N2bin,3))
      allocate(matrxtmpH(Nbin*N2bin,3))
     
     
        do k=1,Nbin

        do 15 i=1,N2bin
        matrxtmpG(i+(k-1)*N2bin,1)=Emin+k-1
        matrxtmpH(i+(k-1)*N2bin,1)=Emin+k-1
15       continue
        end do 
     
        do k=1,Nbin
        do  i=1,N2bin
           matrxtmpG(i+(k-1)*N2bin,2)=E2min+i-1
           matrxtmpH(i+(k-1)*N2bin,2)=E2min+i-1
        end do
        end do


        do   k=1,Nbin
         do 22 i=1,N2bin
         matrxtmpG(i+(k-1)*N2bin,3)=WLG(k+ENcut,i)
         matrxtmpH(i+(k-1)*N2bin,3)=WLH(k+ENcut,i)
22         end do
         end do
      
      
      write(*,'(6I8)')Nbin,Emin,Emax,N2bin,E2min,E2max
      open(40,FILE='mc-gnm_cut2.dat')
      
      
      write(40,'(I15,A18,A4)')SAMCPST,' steps. G(Env,Est)'
      write(40,'(6I8,E15.6,I8,E15.8)')Nbin,Emin,Emax,N2bin,E2min
     * ,E2max,WLFact,nWLHf,SumGam
     
C      write(fmt,'(A4,I0,A4)')'(A8,',N2bin,'I18)'
C      write(40,fmt)' ENV\Est',(E2min+k-1,k=1,N2bin)
C      write(fmt,'(A4,I0,A6)')'(I8,',N2bin,'E18.8)'

      open(41,FILE='mc-hnm_cut2.dat')
      write(41,'(I15,A18,A4)')SAMCPST,' steps. H(Env,Est)'
      write(41,'(6I8,E15.6,I8,E15.8)')Nbin,Emin,Emax,N2bin,E2min
     * ,E2max,WLFact,nWLHf,SumGam
C      write(fmt2,'(A4,I0,A4)')'(A8,',N2bin,'I18)'
C      write(41,fmt2)' ENV\Est',(E2min+k-1,k=1,N2bin)
C      write(fmt2,'(A4,I0,A4)')'(I8,',N2bin,'I18)'

        do  i=1,Nbin*N2bin
        write(40,*)(matrxtmpG(i,k),k=1,3)
        write(41,*)(matrxtmpH(i,k),k=1,3)
        end do
C     do i=1,Nbin
C       write(40,fmt)Emin+i-1,(WLG(i+ENcut,k),k=1,N2bin)
C       write(41,fmt2)Emin+i-1,(WLH(i+ENcut,k),k=1,N2bin)
C      enddo
      close(40)
      close(41)

      
      deallocate(WLG)
      deallocate(WLH)
      deallocate(Mask)
      deallocate(matrxtmpG)
      deallocate(matrxtmpH)
      
      end program
