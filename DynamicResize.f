      module DynamicResize
      implicit none
      integer Nbin,N2bin,N1binmax
      integer Nbinnew,N2binnew
      integer(kind=8),dimension(:,:),allocatable::WLH,WLHtmp,PrWLH,oWLH
      real*8,dimension(:,:),allocatable::WLG,WLGtmp
      real*8,dimension(:,:),allocatable::PW,PWtmp,Mask,MaskTmp
      real*8,dimension(:,:,:),allocatable::gemod,gemodtmp
      integer(kind=8),dimension(:,:,:,:),allocatable::StAn,StAntmp
      integer,dimension(:),allocatable:: E1max,E1min,N1bin
      End module DynamicResize

