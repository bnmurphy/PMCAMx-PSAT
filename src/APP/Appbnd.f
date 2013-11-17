      subroutine Appbnd (i,j,k,l,n4d)
c
c     This is a subroutine assigned the boundary concentrations to the 
c     appropriate cells and sets the concentration of the remaining sources
c     in those cells to 0.
c
c     DEVELOPED BY:
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     08/18/2006
c
c     CALLED BY:
c     readbnd.f
c
c     ROUTINES CALLED:
c     none
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
      integer nx, ny, nz, loc, m
c
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)
c
c-----Setting concentrations values in the boundary cells
c
      if (Appmap(l).ne.0) then  !BNM Fix
        do m = 1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &        nx*ny*nz*MXTRK*(m-1)
          if (m.eq.2) then
            Appconc(loc) = conc(n4d)
          else
            Appconc(loc) = 0
          endif
        enddo
      endif
c
      end
        
