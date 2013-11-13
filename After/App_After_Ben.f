
      program after
c     This program processes PMCAMx results
c     The parameters are HARDCODED for the specific output file
c
c     BE VERY CAREFUL WITH THIS CODE!
c     IT IS VERY HARDCODED!  (BNM 8-23-11)
      integer simdays

      parameter(ifin = 10)
      parameter(ifout = 30)

      parameter(ncol = 150)
      parameter(nrow = 162)
      parameter(nlay = 1)
      parameter(nspc = 161) ! Number of avg species
      parameter(nsrc = 1) !Number of sources
      parameter(nhr  = 24)
      parameter(mspc = 18) ! Number of species to be processed
      parameter(simdays = 7)
      parameter(nday = simdays)
      parameter(naero = 0)
      parameter(UOMOC = 1.4, ROMOC = 2.0)

      character*199 fname1, outdir
      character*99 fname2
      character*198 fname
      character*198 fdir
      character*198 fnamex
      character*4  name(10),note(60),mspec(10,nspc)
      character*10 spcname(nspc)
      character*4  pname(mspc)
      character*3  cdate(simdays)
      character*2  chr(nhr),csrc(nsrc)

      real*4 btime,etime,rdum,xorg,yorg,delx,dely
      real*4 tbtime
      real*4 conh(nspc,ncol,nrow,nhr)
      real*4 sconh(mspc,ncol,nrow)
      real*4 maxh(mspc)
      real Average(ncol*nrow,2*mspc,nsrc)
      real Appave(ncol*nrow*nspc*nsrc)

      integer*4 ione,nspec,ibdate,iedate,idum,izero,nx,ny,nz
      integer*4 itbdate
      integer*4 ihr,msec,icut
      integer*4 iday,ix,iy
      integer i,j,s,n,iloc1,iloc2,indx

      data fdir /'/home/bnmurphy/Research/Output/PSAT/CAMx.NTSOA.PSAT.PITLOC.080911/'/
      data cdate /'193','194','195','196','197','198','199'/  !,'200'/
       !&            '201','202','203','204','205','206','207','208',
       !&            '209'/
c      data cdate /'001','002','003','004','005','006','007','008',
c     &            '009','010','011','012','013','014','015','016',
c     &            '017'/
c      data cdate /'091','092','093','094','095','096','097','098',
c     &            '099','100','101','102','103','104','105','106',
c     &            '107'/
c      data cdate /'274','275','276','277','278','279','280','281',
c     &            '282','283','284','285','286','287','288','289',
c     &            '290'/
c      data cdate /'193','194','195','196','197','198','199','200',
c     &            '201','202','203'/
      data chr /'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &          '21','22','23','24'/
      data pname /'PTOC','SOA_','PEC_','PNH4','PNO3','PSO4',
     &            'PTOT','POC_','SO2_','AVOC','NOX_','NH3_','POA_',
     &            'PTOM','BVOC','ACOG','BCOG','FCOG'/
c      data csrc /'01','02','03','04','05','06','07','08','09','10','11',
c     &           '12','13','14'/
      data csrc /'01','02','03','04','05','06','07','08','09','10','11','12','13'/
c      data csrc /'01','02','03','04'/
c
c Loop for PM 2.5 and PM 10
c      do icut=1,2 
c       msec = 2*icut + 4
c        if (msec.eq.6) then
c          write(6,*)'PM 2.5 and gases'
c        endif
c        if (msec.eq.8) then
c          write(6,*)'PM 10'
c        endif
c
c initialization
c
      do i=1,mspc
        maxh(i) = 0.0
      enddo
c      do i=1,ncol*nrow
c        do j=1,mspc+naero
c         do k=1,nsrc
c            Average(i,j,k) = 0.0
c          enddo
c        enddo
c      enddo
c
c iterate days
c
      do idate = 1,nday
c
c open input file
c
      fname1='/home/bnmurphy/Research/Output/PSAT/'
      outdir='/home/bnmurphy/Research/processed_output/PSAT/'
      isblk=INDEX(fname1,' ')
      fname1=fname1(1:isblk-1)//'CAMx.NTSOA.PSAT.PITLOC.080911/'
      isblk=INDEX(outdir,' ')
      outdir=outdir(1:isblk-1)//'CAMx.NTSOA.PSAT.PITLOC.080911/'
      fname2='4rpos.baseE.'//cdate(idate)//'.CAMx.NTSOA.PSAT.PITLOC.080911.A1'
      isblk=INDEX(fname1,' ')
      fname=fname1(1:isblk-1)//fname2
      isblk=INDEX(fname,' ')
      write(6,*)'INPUT FILE: ',fname(1:isblk-1)
      open(unit=20,file=fname(1:isblk-1),form='UNFORMATTED',
     &     status='old')
c
c read hourly portion
c
      do ihr=1,24
        write(6,*)'ihr = ',ihr
c        write(6,*) ncol,nrow,nsrc,nspc
        read(20) (Appave(i),i=1,ncol*nrow*nspc*nsrc)
c
      do i=1,ncol*nrow
        do j=1,2*mspc
          do k=1,nsrc
            Average(i,j,k) = 0.0
          enddo
        enddo
      enddo
c
c calculate desired quantities
c
        do i=1,ncol
          do j=1,nrow
            iloc1 = i+ncol*(j-1)
            do s=1,nsrc
              n=1 !PTOC
                ist = 65  !Get POC
                do spc=ist,ist+5 !PM2.5
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)/ROMOC
                enddo
                do spc=125,132
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)/UOMOC
                enddo
                do spc=133,158
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)/ROMOC
                enddo

              n=14  !PTOM
                ist = 65  !Get POC
                do spc=ist,ist+5 !PM2.5
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo
                do spc=125,132
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo
                do spc=133,158
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo

              n=2 !SOA
                do spc=133,158
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo

              n=13 !POA
                do spc=125,132
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo

              n=3 !PEC
                ist = 65
                do spc=ist,ist+5 !PM2.5
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo
                do spc=ist,ist+7 !PM10
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n+mspc,s) = Average(iloc1,n+mspc,s)+
     &                                    Appave(iloc2)
                enddo
              n=4 !PNH4
                spc = 160
                iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)

              n=5 !PNO3
                spc = 159
                iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)

              n=6 !PSO4
                ist = 115
                do spc=ist,ist+5 !PM2.5
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo
                do spc=ist,ist+7 !PM10
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n+mspc,s) = Average(iloc1,n+mspc,s)+
     &                                    Appave(iloc2)
                enddo
              n=7 !PTOT
                do ict=65,124,10  !PEC,POC,CRST,PCL,NA,PSO4
                  do ict2=ict,ict+5 !PM2.5
                    spc = ict2
                    iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                    ncol*nrow*nspc*(s-1)
                    Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                  enddo
                enddo
                do ict=125,160    !APO,AOO,ANS,ABS,AAS,PNO3,PNH4
                  spc = ict2
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo

              n=8 !POC
                ist = 75
                do spc=ist,ist+5 !PM2.5
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo
                do spc=ist,ist+7 !PM10
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n+mspc,s) = Average(iloc1,n+mspc,s)+
     &                                    Appave(iloc2)
                enddo

              n=9 !SO2
                spc=15
                iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                ncol*nrow*nspc*(s-1)
                Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)

              n=10 !AVOC
                do spc=5,6  !OLE1,OLE2
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo
                do spc=8,9    !ALK5,ARO1
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo
                spc=21     !ARO2
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                spc=56     !ALK4
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)

              n=15  !BVOC
                spc=7 !TERP
                iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                spc=12 !ISOP
                iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)

              n=16  !ACOG
                do spc=51,55  !CAS
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo

              n=17  !BCOG
                do spc=48,50  !CBS
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo

              n=18  !FCOG
                do spc=22,37  !CPO,COO,CNS
                  iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                  Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                enddo

              n=11 !NOx
                spc=1
                iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
                spc=2
                iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)

              n=12 !NH3
                spc=17
                iloc2 = i+ncol*(j-1)+ncol*nrow*(spc-1)+
     &                  ncol*nrow*nspc*(s-1)
                Average(iloc1,n,s) = Average(iloc1,n,s)+Appave(iloc2)
            enddo
          enddo
        enddo

        do ict=1,mspc
c          do ict2=4,5
c            iday = int((ihr-1)/24)+1
c            ih = ihr - (iday-1)*24
c            if (ict.le.naero) then
c              fname1=outdir(1:INDEX(outdir,' ')-1)
c     &               //'/'//pname(ict)//'25.'//cdate(idate)//
c     &               '.'//chr(ihr)
c            elseif (ict.gt.mspc) then
c              fname1=outdir(1:INDEX(outdir,' ')-1)
c     &               //'/'//pname(ict-mspc)//'10.'//cdate(idate)//
c     &               '.'//chr(ihr)
c              goto 433
c            else
              fname1=outdir(1:INDEX(outdir,' ')-1)
     &               //'/'//pname(ict)//'.'//cdate(idate)//
     &               '.'//chr(ihr)
c            endif
            open(unit=30,file=fname1,form='FORMATTED',status='new')
            do i=1,ncol
            do j = 1,nrow
              indx = (j-1)*ncol + i
              write(30,'I2,1x,I2,1x,30(E10.3,1x)') i,j,(Average(indx,ict,ict2),ict2=1,nsrc)
            enddo
            enddo
            close(30)
        enddo
      enddo
      enddo
      close(20)
c
      end
