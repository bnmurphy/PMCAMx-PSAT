      subroutine Appaero(orig,final,ii,jj,kk,convfac)
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
      real orig(MXSPEC),final(MXSPEC),MW(8),conv,convfac
      integer i,j,k,aspc,gspc,n,s,loc,locA,locG,ii,jj,kk
      integer nx,ny,nz,a,g,coorA(8),coorG(8),count
      real total(MXTRK),frac,allsize(8,MXSOUR),allsizetot(8)
      integer loc2,loc3,loc4,spc,coorP(5)
c
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)
      i=ii
      j=jj
      k=kk
c
c     ASSIGNING MW FOR GAS SPECIES
c
      MW(1) = 98.1 !H2SO4 (SULF)
      MW(2) = 150 !CG1
      MW(3) = 150 !CG2
      MW(4) = 150 !CG3
      MW(5) = 180 !CG4
      MW(6) = 36.5 !HCl
      MW(7) = 17 !NH3
      MW(8) = 63 !HNO3
c
c     MAP GAS TO AEROSOL
      coorA(1) = 145 !PSO4
      coorA(2) = 25 !SOA1
      coorA(3) = 35 !SOA2
      coorA(4) = 45 !SOA3
      coorA(5) = 55 !SOA4
      coorA(6) = 105 !PCL
      coorA(7) = 125 !PNH4
      coorA(8) = 135 !PNO3
      coorG(1) = 15 !SO2
      coorG(2) = 19 !CG1
      coorG(3) = 20 !CG2
      coorG(4) = 21 !CG3
      coorG(5) = 22 !CG4
      coorG(6) = 24 !HCL
      coorG(7) = 17 !NH3
      coorG(8) = 11 !HNO3
      !Primary species
      coorP(1) = 65
      coorP(2) = 75
      coorP(3) = 85
      coorP(4) = 95
      coorP(5) = 115
c
c     GETTING TOTALS
      do n=1,MXTRK
        total(n) = 0
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          total(n) = total(n) + Appconc(loc)
        enddo
        if (n.lt.25) conv = convfac
        if (n.ge.25) conv = 1
c
c     CHECKING TOTALS AT BEGINNING
        if (abs(total(n)-orig(Appmaprev(n))).gt.
     &      0.05*MIN(orig(Appmaprev(n)),total(n)).or.
     &      orig(Appmaprev(n)).eq.bdnl(Appmaprev(n))*conv) then 
          if ((abs(orig(Appmaprev(n))-bdnl(Appmaprev(n))*conv).lt.
     &         0.05*bdnl(Appmaprev(n))*conv).or.
     &        (n.eq.4).or.(n.eq.5).or.
     &        (i.eq.1.or.i.eq.96.or.j.eq.1.or.j.eq.90))
     &         then
            do ict=1,Appnum+3
              loc3 = i + nx*(j-1) + nx*ny*(k-1)+nx*ny*nz*(n-1)+
     &               nx*ny*nz*MXTRK*(ict-1)
              Appconc(loc3) = Appconc(loc3)*(orig(Appmaprev(n))
     &                        /total(n))
            enddo
            total(n) = 0.0
            do s=1,Appnum+3
              loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              total(n) = total(n) + Appconc(loc)
            enddo
          else     
            write(6,*) 'Totals not same at beginning of Appaero',i,j,k,
     &                 n,Appmaprev(n)
            write(6,*) 'Total,old:',total(n),orig(Appmaprev(n))
            write(6,*) 'New: ',final(Appmaprev(n))
            write(6,*) 'Bdnl: ',bdnl(Appmaprev(n))*conv,conv
            stop
          endif
        endif
      enddo
c
      conv = convfac
c     TOTALS FOR AEROSOL SPECIES -- ADD ALL SIZE SECTIONS TOGETHER
      do count=1,8
        allsizetot(count) = 0.0
        do n=1,10
          do s=1,Appnum+3
            if (n.eq.1) allsize(count,s) = 0.0
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*
     &            ((coorA(count)-1+n)-1) + nx*ny*nz*MXTRK*(s-1)
            allsize(count,s) = allsize(count,s) + Appconc(loc)
            allsizetot(count) = allsizetot(count) + Appconc(loc)
          enddo
        enddo
      enddo
c
c     SULFATE
      do n=1,10
        a = coorA(1)-1+n
        g = coorG(1)
        aspc = Appmaprev(a)
        gspc = Appmaprev(g)
        do s=1,Appnum+3
          locA = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(a-1) +
     &            nx*ny*nz*MXTRK*(s-1)
          locG = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(g-1) +
     &            nx*ny*nz*MXTRK*(s-1)
          if ((final(aspc)-orig(aspc)).gt.0) then
            Appconc(locA) = Appconc(locA) + (final(aspc) - orig(aspc))*
     &                      Appconc(locG)/orig(gspc)
          else
            Appconc(locA) = final(aspc)*Appconc(locA)/orig(aspc)
          endif
          if (n.eq.10) then
            Appconc(locG) = final(gspc)*Appconc(locG)/orig(gspc)
            locSulf = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(16-1) +
     &             nx*ny*nz*MXTRK*(s-1)
            Appconc(locSulf) = final(Appmaprev(16))*
     &                         Appconc(locSulf)/orig(Appmaprev(16))
          endif
          if (Appconc(locA).eq.'NaN'.or.Appconc(locG).eq.'NaN') 
     &      write(6,*) 'Problem in Appaero'
        enddo
      enddo
c
c     Remaining Species
      do count=2,8
        do n=1,10
          a = coorA(count)-1+n
          g = coorG(count)
          aspc = Appmaprev(a)
          gspc = Appmaprev(g)
          check = 0.0
          do s=1,Appnum+3
            locA = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(a-1) +
     &             nx*ny*nz*MXTRK*(s-1)
            locG = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(g-1) +
     &             nx*ny*nz*MXTRK*(s-1)
            frac = (Appconc(locG)*MW(count)+allsize(count,s))/
     &             (total(g)*MW(count)+allsizetot(count))
            check = check + frac
c            write(6,*) s,check,Appconc(locG),allsize(count,s),
c     &                 total(g),allsizetot(count)
            if (frac.lt.0.or.frac.gt.1) write(6,*) 'Improper frac: ',
     &        frac,i,j,k,a,g,s
            if (abs(1-check).gt.0.001.and.s.eq.Appnum+3) then
              write(6,*) 'frac Total not 1: ',frac,i,j,k,a,g,check,
     &                   count
              do s=1,Appnum+3
                locA = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(a-1) +
     &                 nx*ny*nz*MXTRK*(s-1)
                locG = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(g-1) +
     &                 nx*ny*nz*MXTRK*(s-1)
                frac = (Appconc(locG)*MW(count)+allsize(count,s))/
     &                 (total(g)*MW(count)+allsizetot(count))
                write(6,*) allsize(count,s),Appconc(locG),frac,total(g),
     &                     total(g),allsizetot(count),MW(count)
              enddo
              stop
            endif
            Appconc(locA) = final(aspc)*frac
            if (n.eq.10) Appconc(locG) = final(gspc)*frac
            if (Appconc(locA).eq.'NaN'.or.Appconc(locG).eq.'NaN') 
     &        write(6,*) 'Problem in Appaero'
          enddo
        enddo
      enddo
c
c     Primary Emissions
c     Two cases: Net condensation, Net evaporation
c     Case 1: Net Condensation
      if ((final(Appmaprev(coorP(1)))-orig(Appmaprev(coorP(1)))).lt.0)
     &    then
        do count = 1,5
          do n = 1,10
            spc = coorP(count) + (11-n) - 1  !(11-n) counts backwards
            do s = 1,Appnum+3
              loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) + 
     &              nx*ny*nz*MXTRK*(s-1)
              loc2 = Appmaprev(spc)
              loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*((spc-1)-1) + 
     &              nx*ny*nz*MXTRK*(s-1)
              if (n.eq.10) then !the smallest size distribution
                Appconc(loc) = final(loc2)*Appconc(loc)/orig(loc2)
              elseif ((final(loc2)-orig(loc2)).gt.0) then
                Appconc(loc) = Appconc(loc)+(final(loc2)-orig(loc2))*
     &                         Appconc(loc3)/orig(loc2-1)
              elseif ((final(loc2)-orig(loc2)).le.0) then
                Appconc(loc) = final(loc2)*Appconc(loc)/orig(loc2)
              endif
            enddo
          enddo
        enddo
c     Case 2: Net evaporation
      else
        do count = 1,5
          do n = 1,10
            spc = coorP(count) + n - 1 
            do s = 1,Appnum+3
              loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) + 
     &              nx*ny*nz*MXTRK*(s-1)
              loc2 = Appmaprev(spc)
              loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*((spc+1)-1) + 
     &              nx*ny*nz*MXTRK*(s-1)
              if (n.eq.10) then !the largest size distribution
                Appconc(loc) = final(loc2)*Appconc(loc)/orig(loc2)
              elseif ((final(loc2)-orig(loc2)).gt.0) then
                Appconc(loc) = Appconc(loc)+(final(loc2)-orig(loc2))*
     &                         Appconc(loc3)/orig(loc2+1)
              elseif ((final(loc2)-orig(loc2)).le.0) then
                Appconc(loc) = final(loc2)*Appconc(loc)/orig(loc2)
              endif
            enddo
          enddo
        enddo
      endif
c
c-------------------------------------------------------------------
c------------------------CHECKS-------------------------------------
c-------------------------------------------------------------------
c   90/10 Check
      do spc=1,MXTRK
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &        nx*ny*nz*MXTRK*(3-1)
c        loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &         nx*ny*nz*MXTRK*(4-1)
c        loc3 = Appmaprev(spc)
c        do n=1,4
c          loc4 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &          nx*ny*nz*MXTRK*(n-1)
c          if (Appconc(loc4).lt.0.and.abs(Appconc(loc4)).lt.
c     &         0.0001*final(Appmaprev(spc))) then
c            Appconc(loc4) = 0
c          endif
c        enddo
c        if (abs(Appconc(loc)*9-Appconc(loc2)).gt.0.01*Appconc(loc2).and.
c     &      Appconc(loc).gt.0.001*final(loc3)) then
c          write(6,*) 'ERROR in Appaero: 90/10 split incorrect'
c          write(6,*) 'i,j,k,spc: ',i,j,k,spc
c          write(6,*) 'Actual Conc.(old, new)',orig(loc3)
c     &               ,final(loc3),loc3
c          do s=1,Appnum+3
c            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &            nx*ny*nz*MXTRK*(s-1)
c            write(6,*) s,Appconc(loc)
c          enddo
c          stop
c        endif
c
c
c-----Check total & Negative Concentrations
        tot = 0.0
        do n=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(n-1)
          loc2 = Appmaprev(spc)
c          if (Appconc(loc).lt.0.or.Appconc(loc).eq.'NaN') then
c            write(6,*) 'Negative (or NaN) Value in Appgas'
c            write(6,*) 'i,j,k,spc,: ',i,j,k,spc,Appmaprev(spc)
c            write(6,*) 'Appconc(loc): ', Appconc(loc)
c            write(6,*) 'Actual Conc.(old,new)', orig(loc2), 
c     &                 final(loc2)
c            do s=1,Appnum+3
c              loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &              nx*ny*nz*MXTRK*(s-1)
c              write(6,*) s,Appconc(loc)
c            enddo
c            stop
c          endif
          tot = tot+Appconc(loc)
        enddo
        if (abs(tot-final(loc2)).gt.0.05*MIN(final(loc2),tot)) 
     &       then
c          if (final(Appmaprev(n)).eq.bdnl(Appmaprev(n))) then
c            do ict=1,Appnum+3
c              loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
c     &               nx*ny*nz*MXTRK*(ict-1)
c              Appconc(loc3) = Appconc(loc3)*(new(Appmaprev(n))
c     &                            /conv/tot)
c            enddo
c          else      
            write(6,*) 'ERROR in Appaero: total incorrect'
            write(6,*) 'i,j,k,spc',i,j,k,spc,Appmaprev(spc)
            write(6,*) 'Actual Conc.(old,new)',orig(loc2),
     &                 final(loc2)
            write(6,*) 'total ', tot
            do s=1,Appnum+3
              loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              write(6,*) s,Appconc(loc)
            enddo
            stop
c          endif
        endif
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          loc2 = Appmaprev(spc)
          Appconc(loc) = Appconc(loc)*final(loc2)/tot
        enddo
c
      enddo
c
      end
