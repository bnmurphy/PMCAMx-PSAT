      subroutine Appxyadv(fp,fm,spc,xy,k,ij,Actconc,dep,scale,old,dx)
c
c     This subroutine handles the changes in apportionment due to emissions. 
c  
c     DEVELOPED BY: 
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     04/13/2007
c
c     CALLED BY:
c     
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
      real fp(MX1D),fm(MX1D),Actconc(MXCOL1,MXROW1,MXLAY1,MXSPEC)
      real old(MXCOL1,MXROW1,MXLAY1,MXSPEC)
      integer spc,xy,k,ij,loc,loc2,loc3,s,src,nx,ny,nz,maxi,maxj
      integer mini,minj
      real Apptot(MXCOL1),Original(MXCOL1,MXSOUR),total
      real dep(MX1D),scale(MX1D),dx(MX1D),convfac
c
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)     
c 
      if (xy.eq.1) then !X-Advection
c
c     Totals
        do i=1,MXCOL1
          Apptot(i) = 0
          do src=1,Appnum+3
            loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &            nx*ny*nz*MXTRK*(src-1)
            Apptot(i) = Appconc(loc)+Apptot(i)
            Original(i,src)=Appconc(loc)
          enddo
          loc2 = i + nx*(ij-1) + nx*ny*(k-1) + 
     &           nx*ny*nz*(Appmaprev(spc)-1)
          loc3 = i + nx*(ij-1) + nx*ny*(k-1)
          if (spc.lt.25) then
            convfac = densfac*(273/tempk(loc3))*(press(loc3)/1013)
          else
            convfac = 1
          endif
          if (abs(Apptot(i)-old(i,ij,k,Appmaprev(spc))).gt.
     &        0.01*MIN(Apptot(i),old(i,ij,k,Appmaprev(spc))).or.
     &        old(i,ij,k,Appmaprev(spc)).eq.bdnl(Appmaprev(spc))) then
            if (abs(old(i,ij,k,Appmaprev(spc))-bdnl(Appmaprev(spc))
     &          *convfac).lt.0.05*bdnl(Appmaprev(spc))
     &          *convfac.or.spc.eq.5) then
              do s = 1,Appnum+3
                loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                if (Apptot(i).eq.0) then
                  if (s.eq.1) Appconc(loc) = old(i,ij,k,Appmaprev(spc))
                  if (s.ne.1) Appconc(loc) = 0.0
                else
                  Appconc(loc) = Appconc(loc)*old(i,ij,k,Appmaprev(spc))
     &                         /Apptot(i)
                endif
                Original(i,s)=Appconc(loc)
              enddo
              Apptot(i) = old(i,ij,k,Appmaprev(spc))
            elseif (i.eq.1.or.i.eq.96.or.j.eq.1.or.j.eq.90) then
              do s = 1,Appnum+3
                loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                if (Apptot(i).eq.0) then
                  if (s.eq.1) Appconc(loc) = old(i,ij,k,Appmaprev(spc))
                  if (s.ne.1) Appconc(loc) = 0.0
                else
                  Appconc(loc) = Appconc(loc)*old(i,ij,k,Appmaprev(spc))
     &                         /Apptot(i)
                endif
                Original(i,s)=Appconc(loc)
              enddo
            else
              write(6,*) 'Totals not same at beginning of Appxyadv',i,ij
     &                    ,k,spc
              write(6,*) 'Total,old,bdnl:',Apptot(i),
     &                   old(i,ij,k,Appmaprev(spc)),bdnl(Appmaprev(spc))
              write(6,*) 'New: ',Actconc(i,ij,k,Appmaprev(spc))
              stop
            endif
          endif
        enddo
c
c
        do s=1,Appnum+3
c
          do i=2,MXCOL1-1
c
c           BOUNDARY CELLS
            if (i.eq.1.and.fp(2).eq.0) then
              loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              loc3 = i + nx*(ij-1) + nx*ny*(k-1) + 
     &               nx*ny*nz*(Appmaprev(spc)-1)
              Appconc(loc)=conc(loc3)*(Original(i,s)/Apptot(i)) 
c
            elseif (i.eq.1.and.fp(2).ne.0) then
              loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              Appconc(loc)=Appconc(loc) - (fp(i)+fm(i))*Original(i,s)/
     &                     Apptot(i) + fm(i+1)*Original(i+1,s)/
     &                     Apptot(i+1)
c
            elseif (i.eq.MX1D.and.fm(MX1D-1).eq.0) then
              loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              loc3 = i + nx*(ij-1) + nx*ny*(k-1) + 
     &               nx*ny*nz*(Appmaprev(spc)-1)
              Appconc(loc)=conc(loc3)*(Original(i,s)/Apptot(i))
c
            elseif (i.eq.MX1D.and.fm(MX1D-1).ne.0) then
              loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              Appconc(loc)=Appconc(loc) + fp(i-1)*Original(i-1,s)/
     &                     Apptot(i-1) - (fp(i)+fm(i))*Original(i,s)/
     &                     Apptot(i)
c
c           NON-BOUNDARY CELLS
            else
              loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              Appconc(loc)=Appconc(loc) + (fp(i-1)*Original(i-1,s)/
     &                     Apptot(i-1) - (fp(i)+fm(i))*Original(i,s)/
     &                     Apptot(i) + fm(i+1)*Original(i+1,s)/
     &                     Apptot(i+1))*scale(i)/dep(i)
              if (Appconc(loc).eq.'NaN') 
     &          write(6,*) 'Problem in Appxyadv - x'
            endif
          enddo
        enddo
c
      elseif (xy.eq.2) then !Y-Advection
c
c       Totals
        do j=1,MXROW1
          Apptot(j) = 0
          do src=1,Appnum+3
            loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &            nx*ny*nz*MXTRK*(src-1)
            Apptot(j) = Appconc(loc)+Apptot(j)
            Original(j,src) = Appconc(loc)
          enddo
          loc2 = ij + nx*(j-1) + nx*ny*(k-1) + 
     &           nx*ny*nz*(Appmaprev(spc)-1)
          loc3 = ij + nx*(j-1) + nx*ny*(k-1)
          if (spc.lt.25) then
            convfac = densfac*(273/tempk(loc3))*(press(loc3)/1013)
          else
            convfac = 1
          endif
          if (abs(Apptot(j)-old(ij,j,k,Appmaprev(spc))).gt.
     &        0.01*MIN(Apptot(j),old(ij,j,k,Appmaprev(spc)))) then
            if (abs(old(ij,j,k,Appmaprev(spc))-bdnl(Appmaprev(spc))
     &          *convfac).lt.0.01*old(ij,j,k,Appmaprev(spc)).or.
     &          spc.eq.5) then
              do s = 1,Appnum+3
                loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                if (Apptot(j).eq.0) then
                  if (s.eq.1) Appconc(loc) = old(ij,j,k,Appmaprev(spc))
                  if (s.ne.1) Appconc(loc) = 0.0
                else
                  Appconc(loc) = Appconc(loc)*old(ij,j,k,Appmaprev(spc))
     &                         /Apptot(j)
                endif
                Original(j,s)=Appconc(loc)
              enddo
              Apptot(j) = old(ij,j,k,Appmaprev(spc))
            else
              write(6,*) 'Totals not same at beginning of Appxyadv',ij,j
     &                    ,k,spc,xy,convfac
              write(6,*) 'Total,old,bdnl:',Apptot(j),
     &                   old(ij,j,k,Appmaprev(spc)),bdnl(Appmaprev(spc))
              write(6,*) 'New: ',Actconc(ij,j,k,Appmaprev(spc))
              do s = 1,Appnum+3
                loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                write(6,*) s,Appconc(loc)
              enddo
              stop
            endif
          else
            do s = 1,Appnum+3
              loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              Appconc(loc) = Appconc(loc)*old(ij,j,k,Appmaprev(spc))
     &                       /Apptot(j)
              Original(j,s)=Appconc(loc)
            enddo
            Apptot(j) = old(ij,j,k,Appmaprev(spc))
          endif
        enddo
c
        do s=1,Appnum+3
c
          do j=2,MXROW1-1
c
c           BOUNDARY CELLS
            if (j.eq.1.and.fp(2).eq.0) then
              loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              loc3 = ij + nx*(j-1) + nx*ny*(k-1) + 
     &               nx*ny*nz*(Appmaprev(spc)-1)
              Appconc(loc)=conc(loc3)*(Original(j,s)/Apptot(j)) 
c
            elseif (j.eq.1.and.fp(2).ne.0) then
              loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              Appconc(loc)=Appconc(loc) - (fp(j)+fm(j))*Original(j,s)/
     &                     Apptot(j) + fm(j+1)*Original(j+1,s)/
     &                     Apptot(j+1)
c
            elseif (j.eq.MX1D.and.fm(MX1D-1).eq.0) then
              loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              loc3 = ij + nx*(j-1) + nx*ny*(k-1) + 
     &               nx*ny*nz*(Appmaprev(spc)-1)
              Appconc(loc)=conc(loc3)*(Original(j,s)/Apptot(j))
c
            elseif (j.eq.MX1D.and.fm(MX1D-1).ne.0) then
              loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              Appconc(loc)=Appconc(loc) + fp(j-1)*Original(j-1,s)/
     &                     Apptot(j-1) - (fp(j)+fm(j))*Original(j,s)/
     &                     Apptot(j)
c
c           NON-BOUNDARY CELLS
            else
              loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              Appconc(loc)=Appconc(loc) + (fp(j-1)*Original(j-1,s)/
     &                     Apptot(j-1) - (fp(j)+fm(j))*Original(j,s)/
     &                     Apptot(j) + fm(j+1)*Original(j+1,s)/
     &                     Apptot(j+1))*scale(j)/dep(j)/dx(j)
            endif
c            if (Appconc(loc).eq.'NaN') then
c              write(6,*) 'NaN Problem in Appxyadv - y:',ij,j,k,spc,s
c     &                    ,Appconc(loc)
c              do s=1,4
c                loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &                nx*ny*nz*MXTRK*(s-1)
c                write(6,*) Appconc(loc),Original(j,s)
c              enddo
c              stop
c            endif
          enddo
        enddo
      endif
c
c----Checks
      if (xy.eq.1) then
        maxi = MXCOL1
        mini = 1
        maxj = ij
        minj = ij
      else
        maxi = ij
        mini = ij
        maxj = MXROW1
        minj = 1
      endif
c
c----Checks for 90/10 split
c      do i=mini,maxi
c        do j=minj,maxj
c            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &            nx*ny*nz*MXTRK*(3-1)
c            loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &             nx*ny*nz*MXTRK*(4-1)
c            if (abs(Appconc(loc)*9-Appconc(loc2)).gt.0.1*Appconc(loc2)
c     &          .and.Appconc(loc).gt.1E-15) 
c     &          then
c              write(6,*) 'ERROR in Appxyadv: 90/10 split incorrect'
c              write(6,*) 'Other, Source 1: ',Appconc(loc), Appconc(loc2)
c              write(6,*) 'i,j,k,spc,xy?,ij: ',i,j,k,spc,xy,ij
c              loc = i + nx*(j-1) + nx*ny*(k-1)+nx*ny*nz*(Appmap(l)-1) +
c     &              nx*ny*nz*MXTRK*(1-1)
c              loc2 = i + nx*(j-1) + nx*ny*(k-1)+nx*ny*nz*(Appmap(l)-1) +
c     &               nx*ny*nz*MXTRK*(2-1)
c              write(6,*) 'IC,BC :',Appconc(loc), Appconc(loc2)
c              stop
c            endif
c        enddo
c      enddo
c
c
c-----Check total
      do i=mini,maxi
        do j=minj,maxj
            total = 0.0
            do s=1,Appnum+3
              loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &            nx*ny*nz*MXTRK*(s-1)
              total = total + Appconc(loc)
              if (Appconc(loc).lt.0) then
                write(6,*) 'Negative Value in Appxyadv'
                write(6,*) 'i,j,k,spc,s,xy?: ',i,j,k,spc,s,xy
                write(6,*) 'Appconc(loc): ', Appconc(loc)
                write(6,*) 'Actual Conc.', Actconc(i,j,k,Appmaprev(spc))
                write(6,*) 'total ', total
                do s=1,Appnum+3
                  loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                  nx*ny*nz*MXTRK*(s-1)
                  write(6,*) s,Appconc(loc),Original(j,s),
     &                       Original(j+1,s)
                enddo
                write(6,*) 'fp(i-1),fp(i),fm(i),fm(i+1):',
     &                      fp(j-1),fp(j),fm(j),fm(j+1)
                write(6,*) 'scale,depth,dx: ',
     &                      scale(j),dep(j),dx(j)
                stop
              endif
            enddo
            if (abs(total-Actconc(i,j,k,Appmaprev(spc)))
     &          .gt.0.05*MIN(Actconc(i,j,k,Appmaprev(spc)),total)) then
              write(6,*) 'ERROR in Appxyadv: total incorrect'
              write(6,*) 'i,j,k,spc,xy?,ij: ',i,j,k,spc,xy,ij
              write(6,*) 'Actual Conc.', Actconc(i,j,k,Appmaprev(spc))
              write(6,*) 'total ', total
              do s=1,Appnum+3
                loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                write(6,*) s,Appconc(loc),Original(j,s),Original(j+1,s)
              enddo
              write(6,*) 'fp(i-1),fp(i),fm(i),fm(i+1):',
     &                    fp(j-1),fp(j),fm(j),fm(j+1)
              write(6,*) 'scale,depth,dx: ',
     &                    scale(j),dep(j),dx(j)
              stop
            else
              do s=1,Appnum+3
                loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                Appconc(loc) = Appconc(loc)*
     &            Actconc(i,j,k,Appmaprev(spc))/total
              enddo
            endif
        enddo
      enddo
c
      end
