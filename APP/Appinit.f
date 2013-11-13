      subroutine Appinit
c
c     This subroutine allocates all the concentrations due to the initial
c     counditions to there location in the source array.  The source array
c     is ordered so that it loops over (in order) x, y, z, spec and sources.
c     The order of the source classes is I.C., B.C., unspecified, source #1,
c     source #2, etc.
c
c     DEVELOPED BY:
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     08/18/2006
c
c     CALLED BY:
c     readcnc.f
c
c     ROUTINES CALLED:
c     none
c
c     VARIABLES (common - used):
c     ncol  - Number of columns
c     nrow  - Number of rows
c     nlay  - Number of layers
c     nspec - Number of species
c     conc  - Concentrations (umol/m^3)
c 
c     VARIABLES (declared):
c     nx    - Number of columns
c     ny    - Number of rows
c     nz    - Number of layers
c     i     - Column Counter
c     j     - Row Counter
c     k     - Layer Counter
c     l     - Species Counter
c     loc   - Location in Array (Appconc)
c     n4d   - Location in Array (conc)
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
      integer nx,ny,nz,i,j,k,l,loc
c
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)
      loc = 1
c
c-----Set all App. concentrations to 0 to start
c
      do i=1,MXVEC3D*MXSOUR*MXTRK
        Appconc(i) = 0
      enddo
c
c-----Copying the current concentrations (which only includes the I.C.
c      values) to the I.C. location in the apportionment array.
c
      if(.NOT.lrstrt) then
        do l = 1,nspec
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                if (Appmap(l).ne.0) then
                  n4d = i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(l-1)
                  loc = i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(Appmap(l)-1)               
                  if(Appmap(l).ne.0) Appconc(loc) = conc(n4d)
                  if(i.eq.50.and.j.eq.21.and.k.eq.2.and.
     &               Appmap(l).eq.91) write(6,*) 'IC Reset',conc(n4d)
     &                                           ,loc,Appconc(loc)
                endif
              enddo
            enddo
          enddo
        enddo
      else
         read(iAppinst) (Appconc(i),i=1,MXSOUR*MXVEC3D*MXTRK)
      endif
c
      end
