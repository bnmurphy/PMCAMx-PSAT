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
c     VARIABLES (passed):
c     i     - Column Number
c     j     - Row Number
c     k     - Layer Number
c     l     - Species Number
c     n4d   - Location in conc Array
c
c     VARIABLES (common):
c     ncol  - Number of columns
c     nrow  - Number of rows
c     nlay  - Number of layers
c     nspec - Number of species
c     conc  - Concentrations (umol/m^3)
c     Appnum- Number of sources
c
c     VARIABLES (declared):
c     nx    - Number of columns
c     ny    - Number of rows
c     nz    - Number of layers
c     loc   - Location in Array (Appconc)
c     m     - Source Type Counter
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
      do m = 1,Appnum+3
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &        nx*ny*nz*MXTRK*(m-1)
        if (m.eq.2) then
          Appconc(loc) = conc(n4d)
        else
          Appconc(loc) = 0
        endif
      enddo
c
      end
        
