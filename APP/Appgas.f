      subroutine Appgas(old,new,r,i,j,k,dt,convfac,NO2case,NO2cons)
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
      real old(MXSPEC),new(MXSPEC),frac(24,MXSOUR),r(100)
      integer i,j,k,s,loc,n,loc2,loc3,loc4,spc,ict,NO2case
      real total(24),tot,conv,total2(24),test,NO2cons(3)
      real NXOYtot,oldNO(MXTRK),dt,dummy
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
      do n=1,24
        total(n) = 0
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          total(n) = total(n) + Appconc(loc)
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
          elseif (s.eq.Appnum+3.and.
     &            abs(total(n)-old(Appmaprev(n))/conv).gt.
     &            0.1*MIN(old(Appmaprev(n))/conv,total(n))) 
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
      do n=1,24
c        if (n.le.4) then
          total(n) = 0
          do s=1,Appnum+3
            loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &             nx*ny*nz*MXTRK*(s-1)
            total(n) = total(n) + Appconc(loc3)
          enddo
c        endif
        check = 0.0
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          frac(n,s) = Appconc(loc)/total(n)
          check = check + frac(n,s)
        if (frac(n,s).lt.0) then
           write(6,*) 'Appgas: fraction < 0',frac(n,s),n,s,i,j,k
           write(6,*) Appconc(loc),total(n)
           stop
        endif
        enddo
        if (abs(check-1).gt.0.00001) then
           write(6,*) 'Fraction sum NOT EQUAL to 1',check,n
        endif
      enddo
c
c     Calculate Final Gas Concentrations
      do s=1,Appnum+3
        !PAN
        n=3 
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)
     &               +dt*(( 1.000)*r( 47)*frac(2,s)))/conv
        !HNO3
        n=11
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)
     &              +dt*(( 2.000)*r( 18)+( 1.000)*r( 26)+( 1.000)*r( 41)
     &               +( 1.000)*r( 44)+( 1.000)*r( 67)+( 0.150)*r( 94))
     &               *frac(2,s))/conv
        !HONO
        n=10
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 2.000)*r( 21)*(frac(1,s)+frac(2,s))/2
     &               +( 1.000)*r( 22)*frac(1,s)))/conv
        !PNA
        n=9
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.000)*r( 29)*frac(2,s)))/conv
        !PAR
        n=6
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &                +( 0.220)*r( 56)*frac(5,s)
     &                +( 1.100)*r( 72)*frac(23,s)
     &                +(( 0.250)*r( 75)+( 0.350)*r( 77)
     &                +( 2.400)*r( 78))*frac(2,s)
     &                +(( 1.565)*r( 92)+( 0.360)*r( 93)
     &                +( 1.282)*r( 94)+( 0.832)*r( 95))*frac(13,s)
     &                +( 2.400)*r( 96)*frac(2,s)
     &                +( 0.220)*r( 97)*frac(18,s)))/conv
        !OLE (only consumed)
        n=5
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !TOL (only consumed)
        n=7
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CRES
        n=8
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 0.360)*r( 63)*frac(7,s)
     &               +( 1.000)*r( 65)*(0.3*r(63)*frac(23,s)+0.56*
     &               r(72)*frac(7,s))/(0.3*r(63)+0.56*r(72))
     &               +( 0.200)*r( 72)*frac(23,s)))/conv
        !XYL (only consumed)
        n=23
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !ISOP (only consumed)
        n=12
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !SO2 (only consumed)
        n=15
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !SULF (only produced)
        n=16
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.000)*r( 82)+(1.000)*r( 83))*frac(15,s))/conv
        !ISPD 
        n=13
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               (( 0.750)*r( 75)+( 0.912)*r( 76)+( 0.650)*r( 77)
     &                +( 0.200)*r( 78)+( 0.200)*r( 96))*frac(12,s)))
     &                /conv
        !OLE2 (only comsumed) 
        n=18
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CG1 (only produced)
        n=19
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 0.070)*r( 63)*frac(7,s)
     &               +( 0.044)*r( 72)*frac(23,s)))/conv
        !CG2 (only produced)
        n=20
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 0.137)*r( 63)*frac(7,s)
     &               +( 0.192)*r( 72)*frac(23,s)))/conv
        !CG3 (only produced) 
        n=21
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( .0024)*r( 52)*frac(6,s)
     &               +(( .0024)*r( 56)+( .0024)*r( 57)
     &               +( .0024)*r( 58)+( .0024)*r( 59))*frac(5,s)
     &               +(( 0.036)*r( 66)+( 0.036)*r( 67))*frac(8,s)))
     &               /conv
        !CG4 (only produced)
        n=22
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               (( 0.136)*r( 97)+( 0.136)*r( 98)+( 0.136)*r( 99)
     &               +( 0.136)*r(100))*frac(18,s)))/conv
        !NO2
        n=2
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
c        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
c     &               -1*(r(1)+r(4)+r(5)+r(7)+r(16)+r(17)+r(21)+r(26)
c     &               +r(29)+r(47)+r(55)+r(68)+r(96))*frac(2,s)
c     &               +(r(3)+r(6)+r(15)+2*r(20)+r(79))*frac(1,s)
c     &               +(0.89*r(14)+r(15)+r(16)+r(19)+.2*r(78))*frac(n,s)
c     &               +(r(24))*frac(10,s)
c     &               +(r(30)+r(31))*frac(9,s)
c     &               +r(48)*frac(3,s)))/conv*new(loc2)/(
c     &               old(loc2)+dt*(
c     &               -1*(r(1)+r(4)+r(5)+r(7)+r(16)+r(17)+r(21)+r(26)
c     &               +r(29)+r(47)+r(55)+r(68)+r(96))
c     &               +(r(3)+r(6)+r(15)+2*r(20)+r(79))
c     &               +(0.89*r(14)+r(15)+r(16)+r(19)+.2*r(78))
c     &               +(r(24))
c     &               +(r(30)+r(31))
c     &               +r(48)))
        Appconc(loc)=new(loc2)/conv*frac(n,s)
        !NO
        n=1
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
c        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
c     &               -1*(r(3)+r(6)+2*r(20)+r(21)+r(22)+r(46)+r(28)+
c     &               r(15)+r(64)+r(79)+r(81))*frac(n,s)
c     &               +(r(1)+r(4)+0.11*r(14)+r(16))*frac(2,s)
c     &               +(r(23)+r(25))*frac(10,s)))/conv*new(loc2)/(
c     &               old(loc2)+dt*(
c     &               -1*(r(3)+r(6)+2*r(20)+r(21)+r(22)+r(46)+r(28)+
c     &               r(15)+r(64)+r(79)+r(81))
c     &               +(r(1)+r(4)+0.11*r(14)+r(16))
c     &               +(r(23)+r(25))))
        Appconc(loc)=new(loc2)/conv*frac(n,s)
        !NXOY
        n=4
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !NTR
        n=14
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
      enddo
c---------------------------------------------------------------
c----------------Changes due to chemical loss-------------------
c---------------------------------------------------------------
c     Calculate Totals
      do n=1,24
        total(n) = 0
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          total(n) = total(n) + Appconc(loc)
        enddo
      enddo
c
c     Calculate Fractions
      do n=1,24
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          frac(n,s) = Appconc(loc)/total(n)
        enddo
      enddo
c
c     Calculate Final Gas Concentrations
      do s=1,Appnum+3
        !PAN
        n=3 
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        test=Appconc(loc)
        Appconc(loc)=Appconc(loc)
     &               +dt*((-1.000)*r( 48)*frac(n,s))/conv
        !HNO3
        n=11
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=Appconc(loc)
     &               +dt*((-1.000)*r( 27)*frac(n,s))/conv
        !HONO
        n=10
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=Appconc(loc)+dt*(
     &               (-1.000)*r( 23)+(-1.000)*r( 24)+
     &               (-2.000)*r( 25))*frac(n,s)/conv
        !PNA
        n=9
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=Appconc(loc)+dt*(
     &               ((-1.000)*r( 30)+(-1.000)*r( 31))*frac(n,s)
     &               )/conv
        !PAR
        n=6
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=Appconc(loc)+dt*(
     &                (-1.000)*r( 52)+(-0.110)*r( 52)+(-2.100)*r( 53)
     &                +(-1.000)*r( 57)+(-1.000)*r( 58)+(-1.000)*r( 59)
     &                +(-1.000)*r( 98)+(-1.000)*r(99)
     &                +(-1.000)*r(100))*frac(n,s)/conv
        !CRES
        n=8
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=Appconc(loc)+dt*(
     &               ((-1.000)*r( 66)+(-1.000)*r( 67))*frac(n,s)
     &               )/conv
        !ISPD 
        n=13
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=Appconc(loc)+dt*(
     &               ((-1.000)*r( 92)+(-1.000)*r( 93)+(-1.000)*r( 94)
     &                +(-1.000)*r( 95))*frac(n,s))/conv
      enddo
c-------------------------------------------------------------------
c------------------------CHECKS-------------------------------------
c-------------------------------------------------------------------
      do spc=1,24
c
c-----Check total & Negative Concentrations
        tot = 0.0
        do n=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(n-1)
          loc2 = Appmaprev(spc)
          tot = tot+Appconc(loc)
        enddo
c
        if (abs(tot-new(loc2)/conv).gt.0.01*MIN(new(loc2)/conv,tot)) 
     &       then
          if (new(Appmaprev(spc)).eq.bdnl(Appmaprev(spc))) then
            do ict=1,Appnum+3
              loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &               nx*ny*nz*MXTRK*(ict-1)
              Appconc(loc3) = Appconc(loc3)*(new(Appmaprev(spc))
     &                            /conv/tot)
            enddo
          else      
            write(6,*) 'ERROR in Appgas: total incorrect'
            write(6,*) 'i,j,k,spc',i,j,k,spc
            write(6,*) 'Actual Conc.(old,new)',old(loc2)/conv,
     &                 new(loc2)/conv
            write(6,*) 'total ', tot
            do s=1,Appnum+3
              loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              write(6,*) s,Appconc(loc),frac(spc,s),loc
            enddo
            stop
          endif
        endif
        do ict=1,Appnum+3
          loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &           nx*ny*nz*MXTRK*(ict-1)
          Appconc(loc3) = Appconc(loc3)*(new(Appmaprev(spc))
     &                    /conv/tot)
        enddo
      enddo
c
      end
