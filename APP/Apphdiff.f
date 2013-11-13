      subroutine Apphdiff(flx,i,j,k,spc,original,Apptot,Actconc)
c
c     This subroutine handles the changes in apportionment due to 
c     horizontal diffusion. 
c  
c     DEVELOPED BY: 
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     04/13/2007
c
c     CALLED BY:
c     diffus
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
      real flx(8),original(MXCOL1,MXROW1,MXSOUR)
      real Apptot(MXCOL1,MXROW1),total,Actconc
      integer i,j,k,spc,loc,s,nx,ny,nz
                    !spc is already equal to Appmap(ispc) from diffus.f
c
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)  
c
      total = 0.0
      do s=1,Appnum+3
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        Appconc(loc) = Appconc(loc)+flx(1)*
     &                 original(i-1,j,s)/Apptot(i-1,j)
        Appconc(loc) = Appconc(loc)-flx(2)*
     &                 original(i,j,s)/Apptot(i,j)
        Appconc(loc) = Appconc(loc)+flx(3)*
     &                 original(i+1,j,s)/Apptot(i+1,j)
        Appconc(loc) = Appconc(loc)-flx(4)*
     &                 original(i,j,s)/Apptot(i,j)
        Appconc(loc) = Appconc(loc)+flx(5)*
     &                 original(i,j-1,s)/Apptot(i,j-1)
        Appconc(loc) = Appconc(loc)-flx(6)*
     &                 original(i,j,s)/Apptot(i,j)
        Appconc(loc) = Appconc(loc)+flx(7)*
     &                 original(i,j+1,s)/Apptot(i,j+1)
        Appconc(loc) = Appconc(loc)-flx(8)*
     &                 original(i,j,s)/Apptot(i,j)
        if (Appconc(loc).eq.'NaN') 
     &    write(6,*) 'Problem in Apphdiff'
        total = total + Appconc(loc)
      enddo
c
c CHECKS
c
c   90/10 Check
c      loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &      nx*ny*nz*MXTRK*(3-1)
c      loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &       nx*ny*nz*MXTRK*(4-1)
c      if (abs(Appconc(loc)*9-Appconc(loc2)).gt.0.01*Appconc(loc2)
c     &    .and.Appconc(loc).gt.1E-10) then
c        write(6,*) 'ERROR in Apphdiff: 90/10 split incorrect'
c        write(6,*) 'i,j,k,spc: ',i,j,k,spc
c        write(6,*) 'Actual Conc.', Actconc
c        write(6,*) 'total ', total
c        do s=1,Appnum+3
c          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &          nx*ny*nz*MXTRK*(s-1)
c          write(6,*) s,Appconc(loc),original(i,j,s)
c        enddo
c        write(6,*) 'Fluxes:',(flx(i),i=1,8)
c        stop
c      endif
c
c
c-----Check total & Negative Concentrations
      if (Appconc(loc).lt.0) then
        write(6,*) 'Negative Value in Apphdiff'
        write(6,*) 'i,j,k,spc,: ',i,j,k,spc
        write(6,*) 'Appconc(loc): ', Appconc(loc)
        write(6,*) 'Actual Conc.', Actconc
        write(6,*) 'total ', total
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          write(6,*) s,Appconc(loc),original(i,j,s)
        enddo
        stop
      endif
      if (abs(total-Actconc).gt.0.01*Actconc) then
        write(6,*) 'ERROR in Apphdiff: total incorrect'
        write(6,*) 'i,j,k,spc: ',i,j,k,spc
        write(6,*) 'Actual Conc.', Actconc
        write(6,*) 'total ', total
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          write(6,*) s,Appconc(loc),original(i,j,s)
        enddo
        write(6,*) 'Fluxes:',(flx(i),i=1,8)
        stop
      endif
      do s=1,Appnum+3
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        Appconc(loc) = Appconc(loc)*Actconc/total
c        if (i.eq.91.and.j.eq.86.and.k.eq.13.and.spc.eq.5)
c     &     write(6,*) Appconc(loc),total,Actconc
      enddo
c
      end

