      subroutine vdepread(ncol,nrow,nspc,vdep,modvdep)
c
c-----Created by Kristina 08/06/08
c
c     VDEPREAD just copies the values from the pointer array for vdep to an
c     actual array - modvdep.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003 
c     ENVIRON International Corporation
c          
c     Modifications:
c        4/17/00   Revised diffusion equations to weight fluxes by density
c       12/07/01   added instructions for OMP
c        1/13/03   added deposited mass array
c
c     Input arguments:
c        ncol              number of columns
c        nrow              number of rows
c        nspc              number of species
c        vdep              deposition velocity (m/s)
c
c
c     Output arguments:
c        modvdep           modeled species deposition velocities (m/s)
c
c     Called by:
c        EMISTRNS
c
      include "camx.prm"
      include "filunit.com"
      include "bndary.com"
      include "chmstry.com"
      include "flags.com"
c
      integer ncol,nrow,nspc
      dimension modvdep(ncol,nrow,nspc)
      real vdep(ncol,nrow,nspc)
c
      do i=1,ncol
        do j=1,nrow
          do k=1,nspc
            modvdep(i,j,k)=vdep(i,j,k)
          enddo
        enddo
      enddo
c
      return
      end
