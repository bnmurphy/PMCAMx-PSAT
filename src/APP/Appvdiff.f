      subroutine Appvdiff(fcup,fcdn,i,j,spc,Actconc)
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
      real fcup(MXLAY1),fcdn(MXLAY1),Actconc(MXLAY1),remain
      integer i,j,spc,k,s,nx,ny,nz,loc
      real original(MXLAY1+1,MXSOUR),Apptot(MXLAY1+1),rho(MXLAY1)
c
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1) 
c
c     Save original values
      do k=1,nz
        Apptot(k) = 0
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          original(k,s) = Appconc(loc)
          if (Appconc(loc).eq.'NaN') write(6,*) 'conc off beg. Appvdiff'
     &      ,i,j,k,spc
          Apptot(k) = Apptot(k) + Appconc(loc)
          if (k.eq.14.and.s.eq.2) original(15,s) = 1.0
          if (k.eq.14.and.s.ne.2) original(15,s) = 0.0
        enddo
      enddo
      Apptot(15) = 1.0
c
      do s=1,Appnum+3
        do k=1,nz
c
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
c
          if (-1*fcup(k).gt.Apptot(k)) then !If more than original conc is transported out
            remain = Apptot(k) + fcup(k)
            if (fcdn(k).ge.0) then
              Appconc(loc)=(fcdn(k)+remain)*original(k-1,s)/Apptot(k-1)
            elseif (fcdn(k).le.0) then
              write(6,*) 'Negative conc in Appvdiff - fluxes too high'
            endif
c
          else if (-1*fcdn(k).gt.Apptot(k)) then !If more than original conc is transported out
            remain = Apptot(k) + fcdn(k)
            if (fcup(k).ge.0) then
              Appconc(loc)=(fcup(k)+remain)*original(k+1,s)/Apptot(k+1)
            elseif (fcup(k).le.0) then
              write(6,*) 'Negative conc in Appvdiff - fluxes too high'
            endif
c
          else
            if (fcup(k).ge.0.and.fcdn(k).ge.0) then
              Appconc(loc) = Appconc(loc) + fcup(k)*original(k+1,s)/
     &                       Apptot(k+1) + fcdn(k)*original(k-1,s)/
     &                       Apptot(k-1)
            elseif (fcup(k).ge.0.and.fcdn(k).le.0) then
              Appconc(loc) = Appconc(loc) + fcup(k)*original(k+1,s)/
     &                       Apptot(k+1) + fcdn(k)*original(k,s)/
     &                       Apptot(k)
            elseif (fcup(k).le.0.and.fcdn(k).ge.0) then
              Appconc(loc) = Appconc(loc) + fcup(k)*original(k,s)/
     &                       Apptot(k) + fcdn(k)*original(k-1,s)/
     &                       Apptot(k-1)
            elseif (fcup(k).le.0.and.fcdn(k).le.0) then
              Appconc(loc) = Appconc(loc) + fcup(k)*original(k,s)/
     &                       Apptot(k) + fcdn(k)*original(k,s)/
     &                       Apptot(k)
            endif
          endif
          if (Appconc(loc).eq.'NaN') 
     &      write(6,*) 'Problem in Appvdiff',i,j,k,spc,original(k,s)
        enddo
      enddo
c
c
c----Checks
c
c----Checks for 90/10 split
c      do k=1,nz
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &        nx*ny*nz*MXTRK*(3-1)
c        loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &         nx*ny*nz*MXTRK*(4-1)
c        if ((abs(Appconc(loc)*9-Appconc(loc2)).gt.0.01*Appconc(loc2))
c     &       .and.(Appconc(loc).gt.0.00000001)) then
c          write(6,*) 'ERROR in Appzadv: 90/10 split incorrect'
c          write(6,*) 'Other, Source 1: ',Appconc(loc), Appconc(loc2)
c          write(6,*) 'i,j,k,spc: ',i,j,k,spc
c          loc = i + nx*(j-1) + nx*ny*(k-1)+nx*ny*nz*(spc-1) +
c     &          nx*ny*nz*MXTRK*(1-1)
c          loc2 = i + nx*(j-1) + nx*ny*(k-1)+nx*ny*nz*(spc-1) +
c     &           nx*ny*nz*MXTRK*(2-1)
c          write(6,*) 'IC,BC :',Appconc(loc), Appconc(loc2)
c          stop
c        endif
c      enddo
c
c
c-----Check total
      do k=1,nz
 100    total = 0.0
        check = 0.0
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          if (Appconc(loc).lt.0) then
c           CORRECT IF ONE SOURCE IS NEGATIVE
c            if ((abs(fcup(k))+abs(fcdn(k))).gt.Apptot(k)) then
c              fixtot = 0.0
c              do n=1,Appnum+3
c                loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &                   nx*ny*nz*MXTRK*(n-1)
c                if (Appconc(loc2).gt.0)  fixtot = fixtot+Appconc(loc2)
c              enddo
c              do n=1,Appnum+3
c                if (n.ne.s) then
c                  loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &                   nx*ny*nz*MXTRK*(n-1)
c                  Appconc(loc2) = Appconc(loc2) + Appconc(loc)*
c     &                            Appconc(loc2)/fixtot
c                  if (Appconc(loc2).gt.0) check = 1
c                endif
c              enddo
c              Appconc(loc) = 0.0
c              if (check.eq.1) goto 100
c            endif
c           FINISH CORRECTION
c          if (Appconc(loc).lt.0) then
            write(6,*) 'Negative Value in Appvdiff'
            write(6,*) 'i,j,k,spc: ',i,j,k,spc
            write(6,*) 'Actual Conc.', Actconc(k)
            write(6,*) 'total ', total
            write(6,*) 'fcup,fcdn: ',fcup(k),fcdn(k)
            do n=1,Appnum+3
              loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(n-1)
              write(6,*) n,Appconc(loc),original(k,n),original(k-1,n)
            enddo
            stop
           endif
        enddo
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          total = total + Appconc(loc)
        enddo
        if (abs(total-Actconc(k)).gt.0.01*Actconc(k)) then
          write(6,*) 'ERROR in Appvdiff: total incorrect'
          write(6,*) 'i,j,k,spc: ',i,j,k,spc
          write(6,*) 'Actual Conc.', Actconc(k)
          write(6,*) 'total ', total
          write(6,*) 'fcup,fcdn: ',fcup(k),fcdn(k)
          do s=1,Appnum+3
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &            nx*ny*nz*MXTRK*(s-1)
            write(6,*) s,Appconc(loc),original(k,s),original(k-1,s)
          enddo
          do icount=1,MXLAY1
            write(6,*) icount,fcup(icount),fcdn(icount)
          enddo
          stop
        endif
      enddo
c
c
      end
