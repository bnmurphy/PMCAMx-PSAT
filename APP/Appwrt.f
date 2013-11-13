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
         do k=1,1
           istart = MXCOL1*MXROW1*MXTRK*MXSOUR*(k-1)
           istop = MXCOL1*MXROW1*MXTRK*MXSOUR*(k) - 1
           write(iAppoA(k)) (Appavg(i),i=istart,istop)
         enddo
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
