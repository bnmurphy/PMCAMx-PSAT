      subroutine Appgas(old,new,r,i,j,k,dt,convfac,NO2case,NO2cons,
     &                  NO3conc,NO3new)
c
c     This subroutine handles the changes in apportionment due to gas phase chemistry. 
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
c
      real old(MXSPEC),new(MXSPEC),frac(53,MXSOUR),r(300)
      integer i,j,k,s,loc,n,loc2,loc3,loc4,spc,ict,NO2case,isrc
      real total(53),tot,conv,total2(53),test,NO2cons(3)
      real NXOYtot,oldNO(MXTRK),dt,dummy,NO3conc
      real NOxTotal, NOxRemain, NOconc(MXSOUR),NOxfrac(MXSOUR)
      real gNTR,gPNA,gPAN,gHNO3,gHONO,NOxSplit(MXSOUR),NOxFracs(MXSOUR)
      real lPNA,lPAN,lHNO3,lHONO,NOxTotNew,modelTot,NO3new,NOxEnd
      real cycPAN,cycPNA,dPNA,dPAN,ngPAN,ngPNA,nlPAN,nlPNA
      real cycHONO,dHONO,ngHONO,nlHONO
c
c
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1) 
      loc = i + nx*(j-1) + nx*ny*(k-1)
      conv = 1/convfac
c----------------------------------------------------------------------
c-------------Changes due to chemical production-----------------------
c----------------------------------------------------------------------
c     Calculate Totals
      do n=1,53
        total(n) = 0
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          total(n) = total(n) + Appconc(loc)
c
          if (Appconc(loc).lt.0) then
             write(6,*) 'Negative Concentrations',n
             do is=1,Appnum+3
               loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &               nx*ny*nz*MXTRK*(is-1)
                write(6,*) Appconc(loc)
             enddo
          endif
c
c         Adjust NO and NO2 concentrations
          if (s.eq.Appnum+3.and.n.eq.1) then
            do ict=1,Appnum+3
              loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &               nx*ny*nz*MXTRK*(ict-1)
              oldNO(ict) = Appconc(loc3)
              Appconc(loc3) = Appconc(loc3)*(old(Appmaprev(n))
     &                        /conv/total(n))
              
            enddo
          elseif (s.eq.Appnum+3.and.(n.eq.4.or.n.eq.3)) then
            do ict=1,Appnum+3
              loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &               nx*ny*nz*MXTRK*(ict-1)
              Appconc(loc3) = Appconc(loc3)*(old(Appmaprev(n))
     &                        /conv/total(n))
            enddo
          elseif (s.eq.Appnum+3.and.n.eq.2) then
            NXOYtot = 0
            do ict=1,Appnum+3
              loc4 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(4-1) +
     &               nx*ny*nz*MXTRK*(ict-1)
              NXOYtot = NXOYtot+Appconc(loc4)
            enddo
            if (NO2case.eq.1) then
              do ict=1,Appnum+3
                loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &                 nx*ny*nz*MXTRK*(ict-1)
                loc4 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(4-1) +
     &                 nx*ny*nz*MXTRK*(ict-1)
                Appconc(loc3) = Appconc(loc3)+Appconc(loc4)/NXOYtot*
     &                  NO2cons(2)/conv
              enddo
            elseif (NO2case.eq.2) then
              do ict=1,Appnum+3
                loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &                 nx*ny*nz*MXTRK*(ict-1)
                loc4 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(4-1) +
     &                 nx*ny*nz*MXTRK*(ict-1)
                Appconc(loc3) = Appconc(loc3)+Appconc(loc4)/NXOYtot*
     &                 NO2cons(2)/conv+Appconc(loc3)/total(2)*
     &                 NO2cons(3)/conv
              enddo
            elseif (NO2case.eq.3) then
              do ict=1,Appnum+3
                loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &                 nx*ny*nz*MXTRK*(ict-1)
                loc4 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(4-1) +
     &                 nx*ny*nz*MXTRK*(ict-1)
                loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(1-1) +
     &                 nx*ny*nz*MXTRK*(ict-1)
                Appconc(loc3) = Appconc(loc3)+Appconc(loc4)/NXOYtot*
     &                 NO2cons(2)/conv+oldNO(ict)/total(1)*
     &                 NO2cons(3)/conv
              enddo
            endif
          elseif ((n.ne.4).and.s.eq.Appnum+3.and.
     &            abs(total(n)-old(Appmaprev(n))/conv).gt.
     &            0.01*MIN(old(Appmaprev(n))/conv,total(n))) 
     &            then 
            if (abs(old(Appmaprev(n))-bdnl(Appmaprev(n))).lt.
     &          0.05*bdnl(Appmaprev(n)).or.n.eq.5) then
              do ict=1,Appnum+3
                loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &                 nx*ny*nz*MXTRK*(ict-1)
                Appconc(loc3) = Appconc(loc3)*(old(Appmaprev(n))
     &                          /conv/total(n))
              enddo
            else     
              write(6,*) 'Totals not same at beginning of Appgas',i,j,k,
     &                   n,Appmaprev(n)
              write(6,*) 'Total,old:',total(n),
     &                   old(Appmaprev(n))/conv,old(Appmaprev(n))
              write(6,*) 'New: ',new(Appmaprev(n))/conv
              write(6,*) 'Bdnl: ',bdnl(Appmaprev(n))
              stop
            endif
          endif
        enddo
      enddo
c
c     Calculate Fractions
      do n=1,53
          total(n) = 0
          do s=1,Appnum+3
            loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &             nx*ny*nz*MXTRK*(s-1)
            total(n) = total(n) + Appconc(loc3)
          enddo
          if (n.eq.2) NOxTotal = total(1)+total(2)
        check = 0.0
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          frac(n,s) = Appconc(loc)/total(n)
          check = check + frac(n,s)
c
          if (n.eq.1) NOconc(s) = Appconc(loc)
          if (n.eq.2) NOxfrac(s) = (Appconc(loc)+NOconc(s))/NOxTotal
c
          if (frac(n,s).lt.0.and.(n.ge.3.).and.(n.ne.4)) then
             write(6,*) 'Appgas: fraction < 0',frac(n,s),n,s,i,j,k
             write(6,*) Appconc(loc),total(n)
             stop
          endif
        enddo
        if (abs(check-1).gt.0.00001.and.(n.ge.3.).and.(n.ne.4)) then
           write(6,*) 'Fraction sum NOT EQUAL to 1',check,n
        endif
      enddo
c
c     Calculate Final Gas Concentrations
c     
c-------B is the branching ratio (NOx dependent SOA yields)
c        B = 0.0
        B =(r( 46)+r( 48)+r( 51)+r( 53)+r( 56)+r( 58)+r( 62)+r( 65)) /
     &     ( r( 46)+r( 48)+r( 51)+r( 53)+r( 56)+r( 58)+r( 62)+r( 65) +
     &     r( 47)+r( 49)+r( 50)+r( 52)+r( 54)+r( 55)+r( 57)+r( 59)+
     &     r( 60)+r( 61)+r( 63)+r( 64)+r( 66)+r( 67)+r( 68) )
c
      do s=1,Appnum+3
        !SULF (only produced)
        n=16
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.000)*r( 44))*frac(15,s))/conv
        !SO2 (only consumed)
        n=15
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO1 (only consumed)
        n=20
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO2 (only consumed)
        n=21
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO3 (only consumed)
        n=22
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO4 (only consumed)
        n=23
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO5 (only consumed)
        n=24
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO6 (only consumed)
        n=25
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO7 (only consumed)
        n=26
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO8 (only consumed)
        n=27
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !COO1
        n=28
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(213)*frac(21,s)
     &               ( 1.075)*r(221)*frac(29,s))/conv
        !COO2
        n=29
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(214)*frac(22,s)
     &               ( 1.075)*r(222)*frac(30,s)
     &               (-1.000)*r(221)*frac(29,s))/conv
        !COO3
        n=30
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(215)*frac(23,s)
     &               ( 1.075)*r(223)*frac(31,s)
     &               (-1.000)*r(222)*frac(30,s))/conv
        !COO4
        n=31
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(216)*frac(24,s)
     &               ( 1.075)*r(224)*frac(32,s)
     &               (-1.000)*r(223)*frac(31,s))/conv
        !COO5
        n=32
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(217)*frac(25,s)
     &               ( 1.075)*r(225)*frac(33,s)
     &               (-1.000)*r(224)*frac(32,s))/conv
        !COO6
        n=33
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(218)*frac(26,s)
     &               ( 1.075)*r(226)*frac(34,s)
     &               (-1.000)*r(225)*frac(33,s))/conv
        !COO7
        n=34
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(219)*frac(27,s)
     &               ( 1.075)*r(227)*frac(35,s)
     &               (-1.000)*r(226)*frac(34,s))/conv
        !COO8
        n=35
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(220)*frac(28,s)
     &               ( 1.075)*r(228)*frac(36,s)
     &               (-1.000)*r(227)*frac(35,s))/conv
        !CNS1
        n=36
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(239)*frac(37,s)
     &               (-0.000)*r(238)*frac(36,s))/conv
        !CNS2
        n=37
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(240)*frac(38,s)
     &               (-1.000)*r(239)*frac(37,s))/conv
        !CNS3
        n=38
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(241)*frac(39,s)
     &               (-1.000)*r(240)*frac(38,s))/conv
        !CNS4
        n=39
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(242)*frac(40,s)
     &               (-1.000)*r(241)*frac(39,s))/conv
        !CNS5
        n=40
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(243)*frac(41,s)
     &               (-1.000)*r(242)*frac(40,s))/conv
        !CNS6
        n=41
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(244)*frac(42,s)
     &               (-1.000)*r(243)*frac(41,s))/conv
        !CNS7
        n=42
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(245)*frac(43,s)
     &               (-1.000)*r(244)*frac(42,s))/conv
        !CNS8
        n=43
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(246)*frac(44,s)
     &               (-1.000)*r(245)*frac(43,s))/conv
