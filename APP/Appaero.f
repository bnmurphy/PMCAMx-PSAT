      subroutine Appaero(orig,final,ii,jj,kk,convfac)
c
c     This subroutine handles the changes in apportionment due to aerosol 
c     processes. 
c  
c     DEVELOPED BY: 
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     04/13/2007
c
c     CALLED BY:
c     chemdriv
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
      real orig(MXSPEC),final(MXSPEC),MW(38),conv,convfac
      integer i,j,k,aspc,gspc,n,s,loc,locA,locG,ii,jj,kk,v
      integer nx,ny,nz,a,g,coorA(38),coorG(38),count
      real total(MXTRK),frac,allsize(38,MXSOUR),allsizetot(38)
      integer loc2,loc3,loc4,spc,coorP(5), nend
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
      MW(2) = 250 !CPO1
      MW(3) = 250 !CPO2
      MW(4) = 250 !CPO3
      MW(5) = 250 !CPO4
      MW(6) = 250 !CPO5
      MW(7) = 250 !CPO6
      MW(8) = 250 !CPO7
      MW(9) = 250 !CPO8
      MW(10) = 250 !COO1
      MW(11) = 250 !COO2
      MW(12) = 250 !COO3
      MW(13) = 250 !COO4
      MW(14) = 250 !COO5
      MW(15) = 250 !COO6
      MW(16) = 250 !COO7
      MW(17) = 250 !COO8
      MW(18) = 250 !CNS1
      MW(19) = 250 !CNS2
      MW(20) = 250 !CNS3
      MW(21) = 250 !CNS4
      MW(22) = 250 !CNS5
      MW(23) = 250 !CNS6
      MW(24) = 250 !CNS7
      MW(25) = 250 !CNS8
      MW(26) = 180 !CBS1
      MW(27) = 180 !CBS2
      MW(28) = 180 !CBS3
      MW(29) = 180 !CBS4
      MW(30) = 180 !CBS5
      MW(31) = 150 !CAS1
      MW(32) = 150 !CAS2
      MW(33) = 150 !CAS3
      MW(34) = 150 !CAS4
      MW(35) = 150 !CAS5
      MW(36) = 36.5 !HCl
      MW(37) = 17 !NH3
      MW(38) = 63 !HNO3
c
c     MAP GAS TO AEROSOL (the gas/aerosol pairs share an index)
      coorA(1) = 115 !PSO4
      coorA(2) = 125 !APO1
      coorA(3) = 126 !APO2
      coorA(4) = 127 !APO3
      coorA(5) = 128 !APO4
      coorA(6) = 129 !APO5
      coorA(7) = 130 !APO6
      coorA(8) = 131 !APO7
      coorA(9) = 132 !APO8
      coorA(10) = 133 !AOO1
      coorA(11) = 134 !AOO2
      coorA(12) = 135 !AOO3
      coorA(13) = 136 !AOO4
      coorA(14) = 137 !AOO5
      coorA(15) = 138 !AOO6
      coorA(16) = 139 !AOO7
      coorA(17) = 140 !AOO8
      coorA(18) = 141 !ANS1
      coorA(19) = 142 !ANS2
      coorA(20) = 143 !ANS3
      coorA(21) = 144 !ANS4
      coorA(22) = 145 !ANS5
      coorA(23) = 146 !ANS6
      coorA(24) = 147 !ANS7
      coorA(25) = 148 !ANS8
      coorA(26) = 149 !ABS1
      coorA(27) = 150 !ABS2
      coorA(28) = 151 !ABS3
      coorA(29) = 152 !ABS4
      coorA(30) = 153 !ABS5
      coorA(31) = 154 !AAS1
      coorA(32) = 155 !AAS2
      coorA(33) = 156 !AAS3
      coorA(34) = 157 !AAS4
      coorA(35) = 158 !AAS5
      coorA(36) = 95  !PCL
      !coorA(36) = 161 !PCL
      coorA(37) = 160 !PNH4
      coorA(38) = 159 !PNO3
      coorG(1) = 15 !SO2
      coorG(2) = 22 !CPO1
      coorG(3) = 23 !CPO2
      coorG(4) = 24 !CPO3
      coorG(5) = 25 !CPO4
      coorG(6) = 26 !CPO5
      coorG(7) = 27 !CPO6
      coorG(8) = 28 !CPO7
      coorG(9) = 29 !CPO8
      coorG(10) = 30 !COO1
      coorG(11) = 31 !COO2
      coorG(12) = 32 !COO3
      coorG(13) = 33 !COO4
      coorG(14) = 34 !COO5
      coorG(15) = 35 !COO6
      coorG(16) = 36 !COO7
      coorG(17) = 37 !COO8
      coorG(18) = 38 !CNS1
      coorG(19) = 39 !CNS2
      coorG(20) = 40 !CNS3
      coorG(21) = 41 !CNS4
      coorG(22) = 42 !CNS5
      coorG(23) = 43 !CNS6
      coorG(24) = 44 !CNS7
      coorG(25) = 45 !CNS8
      coorG(26) = 46 !CBS1
      coorG(27) = 47 !CBS2
      coorG(28) = 48 !CBS3
      coorG(29) = 49 !CBS4
      coorG(30) = 50 !CBS5
      coorG(31) = 51 !CAS1
      coorG(32) = 52 !CAS2
      coorG(33) = 53 !CAS3
      coorG(34) = 54 !CAS4
      coorG(35) = 55 !CAS5
      coorG(36) = 19 !HCL
      coorG(37) = 17 !NH3
      coorG(38) = 11 !HNO3
      !Primary species
      coorP(1) = 65
      coorP(2) = 75
      coorP(3) = 85
      coorP(4) = 105

c
c     GETTING TOTALS
      !print *,'Appaero: Getting totals'
      do n=1,MXTRK
        if (Appmaprev(n).eq.0) goto 200
c        if ( (n.gt.sa_num_gas .and. n.lt.65) .or.n.ge.162.or. 
c     &       (n.ge.95.and.n.lt.104) ) goto 200
        total(n) = 0
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          total(n) = total(n) + Appconc(loc)
        enddo
        if (n.le.sa_num_gas) conv = convfac
        if (n.gt.sa_num_gas) conv = 1
c
c     CHECKING TOTALS AT BEGINNING
        if (abs(total(n)-orig(Appmaprev(n))).gt.
     &      0.10*MIN(orig(Appmaprev(n)),total(n)).or.
     &      orig(Appmaprev(n)).eq.bdnl(Appmaprev(n))*conv) then 
          if ((abs(orig(Appmaprev(n))-bdnl(Appmaprev(n))*conv).lt.
     &         0.05*bdnl(Appmaprev(n))*conv).or.
     &        (n.eq.4).or.(n.eq.5).or.
     &        (i.eq.1.or.i.eq.ncol(1).or.j.eq.1.or.j.eq.nrow(1)).or.
     &        (abs(total(n)-orig(Appmaprev(n))).lt.10*bdnl(Appmaprev(n))*conv)
     &         ) then
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
 200    continue  !Skip the blank arrays between ALK4 and PEC_1
      enddo
c
      conv = convfac
c     TOTALS FOR AEROSOL SPECIES -- ADD ALL SIZE SECTIONS TOGETHER
      do count=1,38
        allsizetot(count) = 0.0
        nend = 1
        if (count.eq.1.or.count.eq.36) nend = 10 !Sulfate and Chloride
        do n=1,nend
          do s=1,Appnum+3
            if (n.eq.1) allsize(count,s) = 0.0
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*
     &            ((coorA(count)-1+n)-1) + nx*ny*nz*MXTRK*(s-1)
            allsize(count,s)  = allsize(count,s)  + Appconc(loc)
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
          if (Appconc(locA).eq.'NaN'.or.Appconc(locG).eq.'NaN') then
            write(6,*) 'Problem in Appaero'
            write(6,*) '   i=',i,' j=',j,' k=',k,' s=',s
            write(6,*) '   a=',a,' Appconc=',Appconc(locA)
            write(6,*) '   g=',g,' Appconc=',Appconc(locG)
          endif
        enddo
      enddo
c
c     Remaining Species (semi-volatiles)
      !print *,'Appaero: Doing semi-volatiles'
      do count=2,38
        nend = 1
        if (count.eq.36) nend=10  !It's PCL, do size resolved partitioning
        do n=1,nend  !BNM remove size from semivolatiles
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
            if (frac.lt.0.or.frac.gt.1) then 
              write(6,*) 'Improper frac: ',frac,i,j,k,a,g,s,count
              write(6,*) '         Appconc(gas)=',Appconc(locG),' total=',total(g)
              write(6,*) '         Allsizetot=',allsizetot(count)
            endif
     
            if (abs(1-check).gt.0.001.and.s.eq.Appnum+3) then
              write(6,*) 'frac Total not 1: ',frac,i,j,k,a,g,check,
     &                   count
              do v=1,Appnum+3
                locA = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(a-1) +
     &                 nx*ny*nz*MXTRK*(v-1)
                locG = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(g-1) +
     &                 nx*ny*nz*MXTRK*(v-1)
                frac = (Appconc(locG)*MW(count)+allsize(count,v))/
     &                 (total(g)*MW(count)+allsizetot(count))
                write(6,*) allsize(count,v),Appconc(locG),frac,total(g),
     &                     total(g),allsizetot(count),MW(count)
              enddo
              stop
            endif
            Appconc(locA) = final(aspc)*frac
            !if (n.eq.10) Appconc(locG) = final(gspc)*frac
            Appconc(locG) = final(gspc)*frac
            if (Appconc(locA).eq.'NaN'.or.Appconc(locG).eq.'NaN') 
     &        write(6,*) 'Problem in Appaero'
          enddo
        enddo
      enddo
c
c     Primary Emissions
c     Two cases: Net condensation, Net evaporation
c     Case 1: Net Condensation
      !print *,'Appaero: Doing primary species'
      if ((final(Appmaprev(coorP(1)))-orig(Appmaprev(coorP(1)))).lt.0)
     &    then
        do count = 1,4
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
        do count = 1,4
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
      !print *,'Appaero: Before Check'
      do spc=1,MXTRK

        if (Appmaprev(spc).eq.0) goto 100
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
        if (abs(tot-final(loc2)).gt.0.10*MIN(final(loc2),tot)) 
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
 100    continue

      enddo
c          
      !print *,'Done with Appaero'
      end
