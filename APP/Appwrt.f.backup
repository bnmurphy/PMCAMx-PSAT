      subroutine Appwrt (Ave)
c            
c     This writes the averages of the concentrations for each hour
c     to the .App output file.
c
c     DEVELOPED BY:
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     08/20/2006
c     
c     CALLED BY:
c     Appave.f
c
c     ROUTINES CALLED:
c     none
c
c     VARIABLES (common):
c     Appavg  -- Hourly averaged source concentrations
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
      integer i,j,k,s
      integer Ave
c
c     Printing the Average File
      if (Ave.eq.1) then
         write(6,*) 'Writing average apportionment concentrations'
         do i=1,MXCOL1
           do j=1,NROW1
             do n=1,MXTRK
               do s=1,Appnum+3
                 iloc=i+MXCOL1*(j-1)+MXCOL1*MXROW1*MXLAY1*(n-1)+
     &                MXCOL1*MXROW1*MXLAY1*MXTRK*(s-1)
c                 if (n.eq.15.and.s.eq.4.and.Appavg(iloc).ne.0)
c     &               write(6,*) i,j,n,s,Appavg(iloc)
               enddo
             enddo
           enddo
         enddo
         write(iAppoA) (Appavg(i),i=1,MXROW1*MXCOL1*
     &                  MXTRK*(Appnum+3))
         write(6,*) MXROW1,MXCOL1,MXTRK,Appnum+3
      endif
c
c     Printing the Instantaneous File
      if (Ave.eq.0) then
         rewind(iAppo)
         write(iAppo) (Appconc(i),i=1,
     &             MXSOUR*MXVEC3D*MXTRK)
      endif
c
      end
