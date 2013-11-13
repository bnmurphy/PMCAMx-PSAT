      subroutine Appwetdep(con)
c
c     This subroutine handles the changes in apportionment due to wet deposition. 
c  
c     DEVELOPED BY: 
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     04/13/2007
c
c     CALLED BY:
c     wetdep 
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
      real con(MXCOL1,MXROW1,MXLAY1,MXSPEC),total
      integer i,j,k,s,spc,nx,ny,nz,v
c
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)
c
      do i=1,nx
        do j=1,ny
          do k=1,nz
            do spc=1,MXSPEC
              if (Appmap(spc).ne.0) then
                total = 0
                do s=1,Appnum+3
                  loc = i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(Appmap(spc)-1)+
     &                nx*ny*nz*MXTRK*(s-1)
                  total = total + Appconc(loc)
                enddo
                do s=1,Appnum+3
                  loc = i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(Appmap(spc)-1)+
     &                  nx*ny*nz*MXTRK*(s-1)
                  Appconc(loc) = con(i,j,k,spc)*Appconc(loc)/total
                  if (Appconc(loc).eq.'NaN') 
     &              write(6,*) 'Problem in Appwetdep',i,j,k,spc,s
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c     CHECKS
c
      do i=1,nx
        do j=1,ny
          do k=1,nz
            do spc=1,MXSPEC
              total = 0
              do s=1,Appnum+3
                loc = i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(Appmap(spc)-1)+
     &                nx*ny*nz*MXTRK*(s-1)
                total = total + Appconc(loc)
              enddo
                if (Appmap(spc).ne.0) then
c------90/10 check
c                  loc = i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(Appmap(spc)-1)+
c     &                  nx*ny*nz*MXTRK*(3-1)
c                  loc2 =i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(Appmap(spc)-1)+
c     &                   nx*ny*nz*MXTRK*(4-1)
c                  if (abs(Appconc(loc)*9-Appconc(loc2)).gt.0.01*
c     &                Appconc(loc2).and.Appconc(loc2).gt.1E-10) then
c                    write(6,*) 'ERROR in Appwetdep: 
c     &                          90/10 split incorrect'
c                    write(6,*) 'Other, Source 1: ',Appconc(loc), 
c     &                          Appconc(loc2)
c                    write(6,*) 'i,j,k,spc: ',i,j,k,Appmap(spc)
c                    loc = i + nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(Appmap(spc)
c     &                    -1) + nx*ny*nz*MXTRK*(1-1)
c                    loc2 = i + nx*(j-1) + nx*ny*(k-1)+nx*ny*nz*
c     &                     (Appmap(spc)-1) + nx*ny*nz*MXTRK*(2-1)
c                    write(6,*) 'IC,BC :',Appconc(loc), Appconc(loc2)
c                    stop
c                  endif
c-------Check for total
                  total = 0.0
                  do s=1,Appnum+3
                    loc=i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(Appmap(spc)-1)+
     &                    nx*ny*nz*MXTRK*(s-1)
                    total = total + Appconc(loc)
                    if (Appconc(loc).lt.0) then
                      write(6,*) 'Negative Value in Appxyadv'
                      write(6,*) 'i,j,k,spc,s: ',i,j,k,Appmap(spc),s
                      write(6,*) 'Appconc(loc): ', Appconc(loc)
                      write(6,*) 'Actual Conc.',con(i,j,k,spc)
                      write(6,*) 'total ', total
                      do v=1,Appnum+3
                        loc = i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(spc-1) +
     &                        nx*ny*nz*MXTRK*(v-1)
                        write(6,*) v,Appconc(loc)
                      enddo
                      stop
                     endif
                  enddo
                  if (abs(total-con(i,j,k,spc))
     &               .gt.0.01*con(i,j,k,spc)) then
                    write(6,*) 'ERROR in Appwetdep: total incorrect'
                    write(6,*) 'i,j,k,spc,: ',i,j,k,Appmap(spc)
                    write(6,*) 'Actual Conc.', con(i,j,k,spc)
                    write(6,*) 'total ', total
                    do s=1,Appnum+3
                      loc=i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(Appmap(spc)
     &                    -1)+nx*ny*nz*MXTRK*(s-1)
                      write(6,*) s,Appconc(loc)
                    enddo
                    stop
                  endif
                endif
            enddo
          enddo
        enddo
      enddo
c
      end
