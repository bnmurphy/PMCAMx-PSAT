      subroutine Appemmap()
c
c     This subroutine maps the Apportionment species to the species in the
c     Apportionment input emissions files.
c  
c     DEVELOPED BY: 
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     04/12/2007
c
c     CALLED BY:
c     startup
c
c     ROUTINES CALLED:
c     none
c
c     VARIABLES (common):
c     Appnam  - Apportionment species names
c     Appmap  - Mapping values between model and apportionment 
c     MXSPEC  - Maximum number of species
c     MXTRK   - Maximum number of tracked species
c     AppemisP- Mapping between apportionment and emissions species - point
c     AppemisA- Mapping between apportionment and emissions species - area
c     
c     VARIABLES (declared):
c     i,j     - Counters
c     match   - Check whether a match was found
c     num     - Numbers 1-10 in text form
c     
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'filunit.com'
      include 'chmstry.com'
      include 'ahomap.com'
      include 'grid.com'
      include 'pigsty.com'
      include 'flags.com'
      include 'App.com'
c
      character*10 ptspec1(MXTRK)
      character*10 arspec1(MXTRK)
      character*4 ifile1(10),note1(60)
      integer s
      real xloc1(120000),yloc1(120000),hstk1(120000)
      real dstk1(120000),tstk1(120000),vstk1(120000)
      real tim11,tim21,orgx1,orgy1,utmx1,utmy1,dx1,dy1
      integer nseg1,nptspc1,idat11,idat21,izone1,nx1,ny1,nz1,idum1
c
c-----Loop over each source (POINT then AREA)
c
      do s=1,Appnum
c
c-----Reading in Headers for Point File
c
c      Header Line 1 (name,note,ione,nspec,ibdate,btime,iedate,etime)
         read(iAppiP(s)) ifile1,note1,nseg1,nptspc1,
     &                   idat11,tim11,idat21,tim21
         nAppP(s) = nptspc1
c      Header Line 2 (rdum,rdun,iutm,xorg,yorg,delx,dely,nx,
c                      ny,nz,idum,idum,rdum,rdum,rdum)
         read(iAppiP(s)) orgx1,orgy1,izone1,utmx1,
     &                   utmy1,dx1,dy1,nx1,ny1,nz1
c      Header Line 3 (ione,ione,nx,ny)
         read(iAppiP(s)) (idum1,idum1,idum1,idum1,n=1,nseg1)
c
c-----Read Species Names -- Point
c
c      Header Line 4 (mspec(l),l=1,nspec)
c
         read(iAppiP(s)) (ptspec1(n),n=1,nptspc1)
         do i=1,nAppP(s)
           do j=1,MXTRK
              if (ptspec1(i).eq.Appnam(j)) then
c                write(6,*) 'POINT: Found #',i,': ',Appnam(j),j
                AppemisP(s,i)=j    !Emissions --> Modeled
              endif
           enddo
         enddo
c
c-----Read More Header Lines -- Point
c
c      Header Line 5 (ione, nstk)
         read(iAppiP(s)) idum1,nptsrc1
c
c      Header Line 6
         read(iAppiP(s)) (xloc1(n),yloc1(n),hstk1(n),dstk1(n),
     &                    tstk1(n),vstk1(n),n=1,nptsrc1)
c
c-----Read Header Lines
c      Header Line 1 (name,note,ione,nspec,ibdate,btime,iedate,etime)
         read(iAppiA(s)) ifile1,note1,nseg1,narspc1,
     &                   idat11,tim11,idat21,tim21
         nAppA(s) = narspc1
c      Header Line 2 (rdum,rdun,iutm,xorg,yorg,delx,dely,nx,ny,nz,idum,idum,rdum,rdum,rdum)
         read(iAppiA(s)) orgx1,orgy1,izone1,utmx1,
     &                   utmy1,dx1,dy1,nx1,ny1,nz1
c      Header Line 3 (ione,ione,nx,ny)
         read(iAppiA(s)) (idum1,idum1,idum1,idum1,n=1,nseg1)
c      Header Line 4 (mspec(l),l=1,nspec)
         read(iAppiA(s)) (arspec1(n),n=1,narspc1)
c
c-----Match and Map Species
c
         do i=1,nAppA(s)
           do j=1,MXTRK
             if (arspec1(i).eq.Appnam(j)) then
c               write(6,*) 'AREA: Found #',i,': ',Appnam(j),j
               AppemisA(s,i)=j
             endif
           enddo
         enddo
c
       enddo
       end
