      subroutine Appzadv(f1,f2,f3,i,j,spc,Actconc)
c
c     This subroutine handles the changes in apportionment due to
c     vertical advection. 
c  
c     DEVELOPED BY: 
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     04/13/2007
c
c     CALLED BY:
c     zadvec
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
      real f1(MXLAY1),f2(MXLAY1),f3(MXLAY1)
      integer i,j,k,spc,s,nx,ny,nz,loc,n,m,check
      real original(MXLAY1+1,MXSOUR),Apptot(MXLAY1),Actconc(MXLAY1)
      real old_conc,fixtot,remain
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
          Apptot(k) = Apptot(k) + Appconc(loc)
        enddo
      enddo
c
      do s=1,Appnum+3
c
c-------Surface (f1 should = 0)
c
          k = 1
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
c
c         ENTRAINMENT
c
          if (-1*f2(k).gt.Apptot(k)) then !If more than original conc is transported out
            remain = Apptot(k) + f2(k)
            if (f1(k).lt.0) then
              write(6,*) 'Negative conc in Appzadv - fluxes too high'
            endif
c
          elseif (-1*f1(k).gt.Apptot(k)) then !If more than original conc is transported out
            remain = Apptot(k) + f1(k)
            if (f2(k).ge.0) then
              Appconc(loc)=(f2(k)+remain)*original(k+1,s)/Apptot(k+1)
            elseif (f2(k).lt.0) then
              write(6,*) 'Negative conc in Appzadv - fluxes too high'
            endif
c
          else
            if (f2(k).ge.0.and.f1(k).ge.0) then
              Appconc(loc) = Appconc(loc) + f2(k)*original(k+1,s)/
     &                       Apptot(k+1)
            elseif (f2(k).ge.0.and.f1(k).lt.0) then
              Appconc(loc) = Appconc(loc) + f2(k)*original(k+1,s)/
     &                       Apptot(k+1) + f1(k)*original(k,s)/
     &                       Apptot(k)
            elseif (f2(k).le.0.and.f1(k).ge.0) then
              Appconc(loc) = Appconc(loc) + f2(k)*original(k,s)/
     &                       Apptot(k)
            elseif (f2(k).le.0.and.f1(k).lt.0) then
              Appconc(loc) = Appconc(loc) + f2(k)*original(k,s)/
     &                       Apptot(k) + f1(k)*original(k,s)/
     &                       Apptot(k)
            endif
          endif
c
c         DILUTION
          Appconc(loc) = Appconc(loc) - f3(k)*original(k,s)/
     &                   Apptot(k)
          if (Appconc(loc).eq.'NaN') 
     &      write(6,*) 'Problem in Appzadv',i,j,k,spc
        
c
c-------Layer 2 to NZ-1
c
        do k=2,nz-1
c
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
c
c         ENTRAINMENT
c         
          if (-1*f2(k).gt.Apptot(k)) then !If more than original conc is transported out
            remain = Apptot(k) + f2(k)
            if (f1(k).ge.0) then
              Appconc(loc)=(f1(k)+remain)*original(k-1,s)/Apptot(k-1)
            elseif (f1(k).lt.0) then
              write(6,*) 'Negative conc in Appzadv - f2 fluxes too high'
            endif
c
          elseif (-1*f1(k).gt.Apptot(k)) then !If more than original conc is transported out
            remain = Apptot(k) + f1(k)
            if (f2(k).ge.0) then
              Appconc(loc)=(f2(k)+remain)*original(k+1,s)/Apptot(k+1)
            elseif (f2(k).lt.0) then
              write(6,*) 'Negative conc in Appzadv - f1 fluxes too high'
            endif
c
          else
            if (f2(k).ge.0.and.f1(k).ge.0) then
              Appconc(loc) = Appconc(loc) + f2(k)*original(k+1,s)/
     &                       Apptot(k+1) + f1(k)*original(k-1,s)/
     &                       Apptot(k-1)
            elseif (f2(k).ge.0.and.f1(k).le.0) then
              Appconc(loc) = Appconc(loc) + f2(k)*original(k+1,s)/
     &                       Apptot(k+1) + f1(k)*original(k,s)/
     &                       Apptot(k)
            elseif (f2(k).le.0.and.f1(k).ge.0) then
              Appconc(loc) = Appconc(loc) + f2(k)*original(k,s)/
     &                       Apptot(k) + f1(k)*original(k-1,s)/
     &                       Apptot(k-1)
            elseif (f2(k).le.0.and.f1(k).le.0) then
              Appconc(loc) = Appconc(loc) + f2(k)*original(k,s)/
     &                       Apptot(k) + f1(k)*original(k,s)/
     &                       Apptot(k)
            endif
          endif
c
c         DILUTION
          Appconc(loc) = Appconc(loc) - f3(k)*original(k,s)/
     &                   Apptot(k)
          if (Appconc(loc).eq.'NaN') 
     &      write(6,*) 'Problem in Appzadv',i,j,k,spc
        enddo
c
c----TOP LAYER 
c    Assign anyting entrained from above to BC source category
c
        k = nz
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
c
c       ENTRAINMENT
c
        if (-1*f2(k).gt.Apptot(k)) then !If more than original conc is transported out
          remain = Apptot(k) + f2(k)
          if (f1(k).ge.0) then
            Appconc(loc)=(f1(k)+remain)*original(k-1,s)/Apptot(k-1)
          elseif (f1(k).lt.0) then
            write(6,*) 'Negative conc in Appzadv - fluxes too high'
            stop
          endif
c
        else if (-1*f1(k).gt.Apptot(k)) then !If more than original conc is transported out
          remain = Apptot(k) + f1(k)
          if (f2(k).ge.0) then
            Appconc(loc) = 0.0
            loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &            nx*ny*nz*MXTRK*(2-1)  !BNM changed s=1 to s=2
            Appconc(loc2)=f2(k)+remain
          elseif (f2(k).lt.0) then
            write(6,*) 'Negative conc in Appzadv - fluxes too high'
            stop
          endif
c
        else
          if (f1(k).ge.0) then
            Appconc(loc) = Appconc(loc) + f1(k)*original(k-1,s)/
     &                     Apptot(k-1)
          elseif (f1(k).le.0) then
            Appconc(loc) = Appconc(loc) + f1(k)*original(k,s)/
     &                     Apptot(k)
          endif
c
          if (f2(k).gt.0.and.s.eq.2) then  !BNM changed s.eq.1 to s.eq.2
            Appconc(loc) = Appconc(loc) + f2(k)
          elseif (f2(k).le.0) then
            Appconc(loc) = Appconc(loc) + f2(k)*original(k,s)/
     &                     Apptot(k)
          endif
        endif
c
c       DILUTION
        Appconc(loc) = Appconc(loc) - f3(k)*original(k,s)/
     &                 Apptot(k)
        if (Appconc(loc).eq.'NaN') 
     &    write(6,*) 'Problem in Appzadv',i,j,k,spc
      enddo
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
          total = total + Appconc(loc)
        enddo
       do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          if (Appconc(loc).lt.0) then
c           CORRECT IF ONE SOURCE IS NEGATIVE
c            if ((abs(f1(k))+abs(f2(k))).gt.Apptot(k)) then
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
c            if (Appconc(loc2).lt.0) then

            !IF THE CONCENTRATION IS LESS THAN THE LOWER BOUND VALUE,
            !SET IT TO ZERO
            !if (abs(Appconc(loc)).lt.bdnl(Appmaprev(spc)) ) then
            if (abs(Appconc(loc)).lt.Actconc(k)/100 .or.
     &          abs(Appconc(loc)).gt.Actconc(k)*100
     &            ) then
              Appconc(loc) = 0.0
            else
              write(6,*) 'Negative Value in Appzadv'
              write(6,*) 'i,j,k,spc: ',i,j,k,spc
              write(6,*) 'Actual Conc.', Actconc(k)
              write(6,*) 'total ', total
              write(6,*) 'fc1,fc2,fc3: ',f1(k),f2(k),f3(k)
              do n=1,Appnum+3
                loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(n-1)
                print *,'loc=',loc
                write(6,*) n,Appconc(loc),original(k,n),original(k-1,n)
              enddo
              stop
             endif
c            endif
          endif
        enddo
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          total = total + Appconc(loc)
        enddo
        if (abs(total-Actconc(k)).gt.0.01*Actconc(k)) then
         if (k.lt.9
     &         ) then
          write(6,*) 'ERROR in Appzadv: total incorrect'
          write(6,*) 'i,j,k,spc: ',i,j,k,spc
          write(6,*) 'Actual Conc.', Actconc(k)
          write(6,*) 'total ', total
          write(6,*) 'fc1,fc2,fc3: ',f1(k),f2(k),f3(k)
          do s=1,Appnum+3
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &            nx*ny*nz*MXTRK*(s-1)
              write(6,*) s,Appconc(loc),original(k,s),original(k-1,s)
            enddo
            stop
          else
            do s=1,Appnum+3
              loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              Appconc(loc) = Appconc(loc)/total * Actconc(k) 
            enddo
          endif
        endif
      enddo
c
c
      end
