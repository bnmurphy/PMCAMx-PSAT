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
c     old     - Concentrations before Trap call. Gases (ppm). Ptcls (ug/m3)
c     new     - Concentrations after  Trap call. Gases (ppm). Ptcls (ug/m3)
c     Appconc - Array with source apportioned concentrations in it (umol/m3, ug/m3)
c     bdnl    - lower bound concentrations. Gases (ppm), Ptcls (ug/m3)
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
      real old(MXSPEC),new(MXSPEC+1),frac(sa_num_gas,MXSOUR),r(MXRXN)
      integer i,j,k,s,loc,n,loc2,loc3,loc4,spc,ict,NO2case,isrc
      real total(sa_num_gas),tot,conv,total2(sa_num_gas),test,NO2cons(3)
      real NXOYtot,oldNO(MXTRK),dt,dummy,NO3conc
      real NOxTotal, NOxRemain, NOconc(MXSOUR),NOxfrac(MXSOUR)
      real gPBZN,gMPAN,gPAN2,gPAN,gHNO3,gNPHE,gHONO,gHNO4,gXN,gRNO3
      real NOxSplit(MXSOUR),NOxFracs(MXSOUR)
      real lPBZN,lMPAN,lPAN2,lPAN,lHNO3,lHONO,lHNO4,lRNO3
      real NOxTotNew,modelTot,NO3new,NOxEnd
      real cycPAN,cycPAN2,cycMPAN,cycHONO
      real dPAN,dPAN2,dMPAN,dHONO
      real ngPAN,ngPAN2,ngMPAN,ngHONO
      real nlPAN,nlPAN2,nlMPAN,nlHONO
      real total3
      integer locB
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
      do n=1,sa_num_gas
        total(n) = 0
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          total(n) = total(n) + Appconc(loc)   !Appconc (umol/m3)
c
          if (Appconc(loc).lt.0) then
             write(6,*) 'Negative Concentrations',n
             do is=1,Appnum+3
               loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &               nx*ny*nz*MXTRK*(is-1)
                write(6,*) Appconc(loc)
             enddo
          endif
cccc  BNM commented out 12-6-11
cc
cc         Adjust NO and NO2 concentrations
c          if (s.eq.Appnum+3.and.n.eq.1) then
c            do ict=1,Appnum+3
c              loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
c     &               nx*ny*nz*MXTRK*(ict-1)
c              oldNO(ict) = Appconc(loc3)
c              Appconc(loc3) = Appconc(loc3)*(old(Appmaprev(n))
c     &                        /conv/total(n))
c              
c            enddo
c          elseif (s.eq.Appnum+3.and.(n.eq.4.or.n.eq.3)) then
c            do ict=1,Appnum+3
c              loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
c     &               nx*ny*nz*MXTRK*(ict-1)
c              Appconc(loc3) = Appconc(loc3)*(old(Appmaprev(n))
c     &                        /conv/total(n))
c            enddo
c          elseif (s.eq.Appnum+3.and.n.eq.2) then
c            NXOYtot = 0
c            do ict=1,Appnum+3
c              loc4 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(4-1) +
c     &               nx*ny*nz*MXTRK*(ict-1)
c              NXOYtot = NXOYtot+Appconc(loc4)
c            enddo
c            
c            if (NO2case.eq.0) then
c              do ict=1,Appnum+3
c                loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
c     &                 nx*ny*nz*MXTRK*(ict-1)
c                Appconc(loc3) = Appconc(loc3)*(old(Appmaprev(n))
c     &                        /conv/total(n))
c              enddo
c            elseif (NO2case.eq.1) then
c              do ict=1,Appnum+3
c                loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
c     &                 nx*ny*nz*MXTRK*(ict-1)
c                loc4 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(4-1) +
c     &                 nx*ny*nz*MXTRK*(ict-1)
c                Appconc(loc3) = Appconc(loc3)+Appconc(loc4)/NXOYtot*
c     &                  NO2cons(2)/conv
c              enddo
c            elseif (NO2case.eq.2) then
c              do ict=1,Appnum+3
c                loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
c     &                 nx*ny*nz*MXTRK*(ict-1)
c                loc4 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(4-1) +
c     &                 nx*ny*nz*MXTRK*(ict-1)
c                Appconc(loc3) = Appconc(loc3)+Appconc(loc4)/NXOYtot*
c     &                 NO2cons(2)/conv+Appconc(loc3)/total(2)*
c     &                 NO2cons(3)/conv
c              enddo
c            elseif (NO2case.eq.3) then
c              do ict=1,Appnum+3
c                loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
c     &                 nx*ny*nz*MXTRK*(ict-1)
c                loc4 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(4-1) +
c     &                 nx*ny*nz*MXTRK*(ict-1)
c                loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(1-1) +
c     &                 nx*ny*nz*MXTRK*(ict-1)
c                Appconc(loc3) = Appconc(loc3)+Appconc(loc4)/NXOYtot*
c     &                 NO2cons(2)/conv+oldNO(ict)/total(1)*
c     &                 NO2cons(3)/conv
c              enddo
c            endif
c          elseif ((n.ne.4).and.s.eq.Appnum+3.and.
c     &            abs(total(n)-old(Appmaprev(n))/conv).gt.
c     &            0.05*MIN(old(Appmaprev(n))/conv,total(n))) 
c     &            then 
c            if (abs(old(Appmaprev(n))-bdnl(Appmaprev(n))).lt.
c     &          0.05*bdnl(Appmaprev(n)).or.n.eq.5) then
c              do ict=1,Appnum+3
c                loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
c     &                 nx*ny*nz*MXTRK*(ict-1)
c                Appconc(loc3) = Appconc(loc3)*(old(Appmaprev(n))
c     &                          /conv/total(n))
c              enddo
c            else     
c              write(6,*) 'Totals not same at beginning of Appgas',i,j,k,
c     &                   n,Appmaprev(n)
c              write(6,*) 'Total,old:',total(n),
c     &                   old(Appmaprev(n))/conv,old(Appmaprev(n))
c              write(6,*) 'New: ',new(Appmaprev(n))/conv
c              write(6,*) 'Bdnl: ',bdnl(Appmaprev(n))
c              stop
c            endif
c          endif
        enddo
        if (n.eq.8.and.i.eq.2.and.j.eq.2) then
          print *,'Appgas: SO2 Total=',total(n),' old=',old(Appmaprev(n))/conv,' new=',new(Appmaprev(n))/conv
        endif
      enddo
c
c     Calculate Fractions
      do n=1,sa_num_gas
          total(n) = 0
          do s=1,Appnum+3
            loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &             nx*ny*nz*MXTRK*(s-1)
            total(n) = total(n) + Appconc(loc3)
          enddo
c          if (n.eq.2) NOxTotal = total(1)+total(2)
        check = 0.0
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          frac(n,s) = Appconc(loc)/total(n)
          check = check + frac(n,s)
c
c          if (n.eq.1) NOconc(s) = Appconc(loc)
c          if (n.eq.2) NOxfrac(s) = (Appconc(loc)+NOconc(s))/NOxTotal
cc
c          if (frac(n,s).lt.0.and.(n.ge.3.).and.(n.ne.4)) then
c             write(6,*) 'Appgas: fraction < 0',(frac(n,is),is=1,Appnum+3),n,s,i,j,k
c             write(6,*) Appconc(loc),total(n)
c             stop
c          endif
        enddo
        if (abs(check-1).gt.0.00001) then !.and.(n.ge.3.).and.(n.ne.4)) then
           write(6,*) 'Fraction sum NOT EQUAL to 1',check,n
        endif
      enddo
cc
cc     Calculate Final Gas Concentrations
cc
cc
cc     NOx calculations: NO, NO2, N2O5 and NO3
c      NOxTotal = old(Appmaprev(1))+old(Appmaprev(2))+old(Appmaprev(4))+
c     &           NO3conc
cc     Write gain and loss terms for non-NOx species from the
cc     perspective of gain and loss to the NOx family
c      gPBZN = +( 1.000)*r( 90)
c      lPBZN = +( 1.000)*r( 91)
c      gMPAN = +( 1.000)*r(102)
c      lMPAN = +( 1.000)*r(103)
c      gPAN2 = +( 1.000)*r( 79)
c      lPAN2 = +( 1.000)*r( 80)
c      gPAN  = +( 1.000)*r( 69)
c      lPAN  = +( 1.000)*r( 70)
c      gNPHE = +( 1.000)*r(117)+( 1.000)*r(121)+( 1.000)*r(122)
c      gHNO3 = +( 2.000)*r( 13)+( 1.000)*r( 25)+( 0.200)*r( 39)
c     &        +( 1.000)*r(129)+( 1.000)*r(132)+( 1.000)*r(135)
c     &        +( 1.000)*r(148)+( 1.000)*r(151)+( 1.000)*r(154)
c     &        +( 1.000)*r(156)+( 1.000)*r(157)+( 1.000)*r(160)
c     &        +( 0.500)*r(163)+( 0.150)*r(172)
c      lHNO3 = +( 1.000)*r( 27)+( 1.000)*r( 28)
c      gRNO3 = +( 1.000)*r(115)+( 0.572)*r(172)+( 0.310)*r(176)
c     &        +( 0.813)*r(191)+( 0.276)*r(195)+( 0.511)*r(206)
c     &        +( 0.321)*r(210)+( 1.000)*r( 62)
c      lRNO3 = +( 1.000)*r(176)+( 1.000)*r(177)
c      
c      gHONO = +( 1.000)*r( 21)
c      lHONO = +( 1.000)*r( 22)+( 1.000)*r( 23)+( 1.000)*r( 24)
c      gHNO4 = +( 1.000)*r( 32)
c      lHNO4 = +( 1.000)*r( 33)+( 1.000)*r( 34)+( 1.000)*r( 35)
c      gXN   = +( 2.000)*r(120)+( 0.500)*r(163)+( 0.278)*r(172)
c     &        +( 0.352)*r(176)+( 1.000)*r(187)+( 0.250)*r(195)
c     &        +( 0.489)*r(206)+( 0.288)*r(210)
c      cycMPAN = MIN(lMPAN,gMPAN)
c      cycPAN2 = MIN(lPAN2,gPAN2)
c      cycPAN  = MIN(lPAN,gPAN)
c      cycHONO = MIN(lHONO,gHONO)
c      dHONO = gHONO-lHONO
c      dPAN  = gPAN-lPAN
c      dPAN2 = gPAN2-lPAN2
c      dMPAN = gMPAN-lMPAN
c
c      !KRISTINA'S CODE
c      !gPBZN = +( 1.000)*r( 90)
c      !lPBZN = +( 1.000)*r( 91)
c      !gMPAN = +( 1.000)*r(102)
c      !lMPAN = +( 1.000)*r(103)
c      !gPAN2 = +( 1.000)*r( 79)
c      !lPAN2 = +( 1.000)*r( 80)
c      !gPAN  = +( 1.000)*r( 69)
c      !lPAN  = +( 1.000)*r( 70)
c      !gNPHE = +( 1.000)*r(117)+( 1.000)*r(121)+( 1.000)*r(122)
c      !gHNO3 = +( 2.000)*r( 13)+( 1.000)*r( 25)+( 0.200)*r( 39)
c      !&        +( 1.000)*r(129)+( 1.000)*r(132)+( 1.000)*r(135)
c      !&        +( 1.000)*r(148)+( 1.000)*r(151)+( 1.000)*r(154)
c      !&        +( 1.000)*r(156)+( 1.000)*r(157)+( 1.000)*r(160)
c      !&        +( 0.500)*r(163)+( 0.150)*r(172)
c      !lHNO3 = +( 1.000)*r( 27)+( 1.000)*r( 28)
c      !gHONO = +( 1.000)*r( 21)
c      !lHONO = +( 1.000)*r( 22)+( 1.000)*r( 23)+( 1.000)*r( 24)
c      !gHNO4 = +( 1.000)*r( 32)
c      !lHNO4 = +( 1.000)*r( 33)+( 1.000)*r( 34)+( 1.000)*r( 35)
c      !gXN   = +( 2.000)*r(120)+( 0.500)*r(163)+( 0.278)*r(172)
c      !&        +( 0.352)*r(176)+( 1.000)*r(187)+( 0.250)*r(195)
c      !&        +( 0.489)*r(206)+( 0.288)*r(210)
c      
c      if (dPAN.gt.0) then
c         ngPAN = dPAN
c         nlPAN = 0
c      else
c         ngPAN = 0
c         nlPAN = -1*dPAN
c      endif
c      if (dPAN2.gt.0) then
c         ngPAN2 = dPAN2
c         nlPAN2 = 0
c      else
c         ngPAN2 = 0
c         nlPAN2 = -1*dPAN2
c      endif
c      if (dMPAN.gt.0) then
c         ngMPAN = dMPAN
c         nlMPAN = 0
c      else
c         ngMPAN = 0
c         nlMPAN = -1*dMPAN
c      endif
c      if (dHONO.gt.0) then
c         ngHONO = dHONO
c         nlHONO = 0
c      else
c         ngHONO = 0
c         nlHONO = -1*dHONO
c      endif
c
c      !Added Gain and Loss of RNO3 to these calculations: BNM 11-14-12
c      NOxRemain = NOxTotal + (-1*gPBZN-gHNO3-ngHONO-ngPAN-ngPAN2-ngMPAN
c     &            -gHNO4-gNPHE-gXN-gRNO3-0.5*cycPAN-0.5*cycPAN2-0.5*cycHONO
c     &            -0.5*cycMPAN)*dt
c      NOxEnd = NOxTotal + (-1*gPBZN-gHNO3-dHONO-dPAN-dPAN2-dMPAN
c     &            -gHNO4-gNPHE-gXN-gRNO3+lPBZN+lHNO3+lHNO4+lRNO3)*dt
c      modelTot = new(Appmaprev(1))+new(Appmaprev(2))+new(Appmaprev(4))+
c     &           NO3new
c
c      !KRISTINA
c      ! NOxRemain = NOxTotal + (-1*gPBZN-gHNO3-ngHONO-ngPAN-ngPAN2-ngMPAN
c      !&            -gHNO4-gNPHE-gXN-0.5*cycPAN-0.5*cycPAN2-0.5*cycHONO
c      !&            -0.5*cycMPAN)*dt
c      ! NOxEnd = NOxTotal + (-1*gPBZN-gHNO3-dHONO-dPAN-dPAN2-dMPAN
c      !&            -gHNO4-gNPHE-gXN+lPBZN+lHNO3+lHNO4)*dt
c      ! modelTot = new(Appmaprev(1))+new(Appmaprev(2))+new(Appmaprev(4))+
c      !&           NO3new
c      if (abs(modelTot-NOxEnd).gt.0.05*MIN(modelTot,NOxEnd)) then
c         write(6,*) 'Big diff in end concentrations',modelTot,NOxEnd,
c     &              i,j,k,NO3new,NO3conc
c      endif
cc
c      NOxTotNew=0.0
c      do s=1,Appnum+3
c         if (NOxRemain.gt.0) then
c         NOxSplit(s)=NOxTotal*frac(2,s)-dt*(gNPHE+gXN+gPBZN+gHNO4+gHNO3+gRNO3+
c     &               0.5*cycHONO+ngHONO+0.5*cycPAN+0.5*cycPAN2+
c     &               0.5*cycMPAN+ngMPAN+ngPAN+ngPAN2)*frac(2,s)+
c     &               dt*(lHNO3*frac(11,s)+lHNO4*frac(20,s)+
c     &               lPBZN*frac(13,s)+(nlHONO+0.5*cycHONO)*frac(10,s)+
c     &               (nlPAN+0.5*cycPAN)*frac(3,s)+
c     &               (nlMPAN+0.5*cycMPAN)*frac(18,s)+
c     &               (nlPAN2+0.5*cycPAN2)*frac(14,s)+lRNO3*frac(57,s))
c         else
c            write(6,*) 'NOx used up ...',i,j,k,NOxRemain,NOxTotal
c           NOxSplit(s)=(NOxRemain+lPAN2+lMPAN+lPBZN+lHNO4+lPAN+lHNO3+
c     &                  lHONO+lRNO3) * (lPAN2*frac(14,s)+lMPAN*frac(18,s)
c     &                  +lPAN*frac(3,s)+lHNO3*frac(11,s)+lHONO*
c     &                   frac(10,s)+lPBZN*frac(13,s)+lHNO4*frac(20,s)
c     &                  +lRNO3*frac(57,s))/
c     &                  (lPAN2+lMPAN+lPBZN+lHNO4+lPAN+lHNO3+lHONO+lRNO3)
c         endif
cc
c         NOxTotNew = NOxTotNew+NOxSplit(s)
c         if (NOxSplit(s).lt.0) then
c            write(6,*) "NOxSplit is less than 0"
c            stop
c         endif
c      enddo
c      check = 0.0
c      do s=1,Appnum+3
c         NOxFracs(s)=NOxSplit(s)/NOxTotNew
c         check = check+NOxFracs(s)
c      enddo
cc
c      if (abs(check-1).gt.0.001) then
c         write(6,*) 'NOx Fraction not equal to 1'
c         do s=1,Appnum+3
c            write(6,*) NOxFracs(s)
c         enddo
c         stop
c      endif
cccc BNM done commenting out NOx apportionment 12-6-11


c-------B is the branching ratio (NOx dependent SOA yields)
c        B = 0.0
        B =(r( 46)+r( 48)+r( 51)+r( 53)+r( 56)+r( 58)+r( 62)+r( 65)) /
     &     ( r( 46)+r( 48)+r( 51)+r( 53)+r( 56)+r( 58)+r( 62)+r( 65) +
     &     r( 47)+r( 49)+r( 50)+r( 52)+r( 54)+r( 55)+r( 57)+r( 59)+
     &     r( 60)+r( 61)+r( 63)+r( 64)+r( 66)+r( 67)+r( 68) )
      ro = 1.4
c
      do s=1,Appnum+3
c        !NO
c        n=1
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
c     &        nx*ny*nz*MXTRK*(s-1)
c        loc2 = Appmaprev(n)
c        Appconc(loc) = new(loc2)*NOxFracs(s)/conv

c        !NO2
c        n=2
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
c     &        nx*ny*nz*MXTRK*(s-1)
c        loc2 = Appmaprev(n)
c        Appconc(loc) = new(loc2)*NOxFracs(s)/conv
        
c        !PAN
c        n=3
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
c     &        nx*ny*nz*MXTRK*(s-1)
c        loc2 = Appmaprev(n)
c        Appconc(loc) = (old(loc2)*frac(n,s)
c     &                 + dt*(( 1.000)*r( 69)*frac(2,s)
c     &                 +     (-1.000)*r( 70)*frac(n,s) ))/conv
     
c        !NXOY
c        n=4
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
c     &        nx*ny*nz*MXTRK*(s-1)
c        loc2 = Appmaprev(n)
c        Appconc(loc) = new(loc2)*NOXFracs(s)/conv

c        !HONO
c        n=10
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
c     &        nx*ny*nz*MXTRK*(s-1)
c        loc2 = Appmaprev(n)
c        Appconc(loc) = (old(loc2)*frac(n,s)
c     &                 +dt*( ( 1.000)*r( 21)*frac(2,s)
c     &                      +((-1.000)*r( 22)+(-1.00)*r(23)
c     &                      + (-1.000)*r( 24) )*frac(n,s)  ) )/conv
        
c        !HNO3
c        n=11
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
c     &        nx*ny*nz*MXTRK*(s-1)
c        loc2 = Appmaprev(n)
c        Appconc(loc) = (old(loc2)*frac(n,s)
c     &                 +dt*( ( (2.00)*r(13)+(1.00)*r(25)+(0.20)*r(39)
c     &                        +(1.00)*r(129)+(1.00)*r(132)+(1.00)*r(135)
c     &                        +(1.00)*r(148)+(1.00)*r(151)+(1.00)*r(154)
c     &                        +(1.00)*r(156)+(1.00)*r(157)+(1.00)*r(160)
c     &                        +(0.50)*r(163)+(0.15)*r(172) )*frac(2,s)
c     &                       +((-1.00)*r(27)+(-1.00)*r(28))*frac(n,s)   ))/conv
      
c        !PBZN
c        n=13
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
c     &        nx*ny*nz*MXTRK*(s-1)
c        loc2 = Appmaprev(n)
c        Appconc(loc) = (old(loc2)*frac(n,s)
c     &                 +dt*( (1.00)*r(90)*frac(2,s) 
c     &                      +(-1.0)*r(91)*frac(n,s) ))/conv

c        !PAN2
c        n=14
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
c     &        nx*ny*nz*MXTRK*(s-1)
c        loc2 = Appmaprev(n)
c        Appconc(loc) = (old(loc2)*frac(n,s)
c     &                 +dt*( (1.00)*r(79)*frac(2,s)
c     &                      +(-1.0)*r(80)*frac(n,s)))/conv

c        !MPAN
c        n=18
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
c     &        nx*ny*nz*MXTRK*(s-1)
c        loc2 = Appmaprev(n)
c        Appconc(loc) = (old(loc2)*frac(n,s)
c     &                 +dt*( (1.00)*r(102)*frac(2,s)
c     &                      +(-1.0)*r(103)*frac(n,s) ))/conv
      
c        !HNO4
c        n=20
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
c     &        nx*ny*nz*MXTRK*(s-1)
c        loc2 = Appmaprev(n)
c        Appconc(loc) = (old(loc2)*frac(n,s)
c     &                 +dt*( (1.00)*r(32)*frac(2,s)
c     &                      +((-1.0)*r(33)+(-1.0)*r(34)
c     &                      + (-1.0)*r(35))*frac(n,s)  ))/conv

c        !RNO3 - Added by BNM 11-14-12
c        n=57
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
c     &        nx*ny*nz*MXTRK*(s-1)
c        loc2 = Appmaprev(n)
c        Appconc(loc) = (old(loc2)*frac(n,s)
c     &                 +dt*((+( 1.000)*r( 62)
c     &                       +( 1.000)*r(115)+( 0.572)*r(172)+( 0.310)*r(176)
c     &                       +( 0.813)*r(191)+( 0.276)*r(195)+( 0.511)*r(206)
c     &                       +( 0.321)*r(210))*frac(2,s) 
c     &                       +((-1.0)*r(176)+(-1.0)*r(177))*frac(n,s) ))/conv


        !SULF (only produced)
        n=9
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.000)*r( 44))*frac(8,s))/conv

        !SO2 (only consumed)
        n=8
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv

        !CPO1 (only consumed)
        n=13
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO2 (only consumed)
        n=14
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO3 (only consumed)
        n=15
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO4 (only consumed)
        n=16
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO5 (only consumed)
        n=17
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO6 (only consumed)
        n=18
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO7 (only consumed)
        n=19
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO8 (only consumed)
        n=20
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !COO1
        n=21
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(213)*frac(14,s)+
     &               ( 1.075)*r(221)*frac(22,s)))/conv
        !COO2
        n=22
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(214)*frac(15,s)+
     &               ( 1.075)*r(222)*frac(23,s)+
     &               (-1.000)*r(221)*frac(n,s)))/conv
        !COO3
        n=23
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(215)*frac(16,s)+
     &               ( 1.075)*r(223)*frac(24,s)+
     &               (-1.000)*r(222)*frac(n,s)))/conv
        !COO4
        n=24
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(216)*frac(17,s)+
     &               ( 1.075)*r(224)*frac(25,s)+
     &               (-1.000)*r(223)*frac(n,s)))/conv
        !COO5
        n=25
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(217)*frac(18,s)+
     &               ( 1.075)*r(225)*frac(26,s)+
     &               (-1.000)*r(224)*frac(n,s)))/conv
        !COO6
        n=26
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(218)*frac(19,s)+
     &               ( 1.075)*r(226)*frac(27,s)+
     &               (-1.000)*r(225)*frac(n,s)))/conv
        !COO7
        n=27
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(219)*frac(20,s)+
     &               ( 1.075)*r(227)*frac(28,s)+
     &               (-1.000)*r(226)*frac(n,s)))/conv
        !COO8
        n=28
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               (-1.000)*r(227)*frac(n,s)))/conv
        !CNS1
        n=29
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(239)*frac(30,s)+
     &               (-0.000)*r(238)*frac(n,s)))/conv
        !CNS2
        n=30
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(240)*frac(31,s)+
     &               (-1.000)*r(239)*frac(n,s)))/conv
        !CNS3
        n=31
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(241)*frac(32,s)+
     &               (-1.000)*r(240)*frac(n,s)))/conv
        !CNS4
        n=32
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(242)*frac(33,s)+
     &               (-1.000)*r(241)*frac(n,s)))/conv
        !CNS5
        n=33
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(243)*frac(34,s)+
     &               (-1.000)*r(242)*frac(n,s)))/conv
        !CNS6
        n=34
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(244)*frac(35,s)+
     &               (-1.000)*r(243)*frac(n,s)))/conv
        !CNS7
        n=35
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(245)*frac(36,s)+
     &               (-1.000)*r(244)*frac(n,s)))/conv
        !CNS8
        n=36
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               (-1.000)*r(245)*frac(n,s)))/conv
        !CBS1
        n=37
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(229)*frac(38,s)+
     &               (-0.000)*r(228)*frac(n,s)))/conv
        !CBS2
        n=38
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(230)*frac(39,s)+
     &               (-1.000)*r(229)*frac(n,s)+
     &               ro*((0.0023*(1-B) + 0.0001*B )*r(189)+
     &               (0.0023*(1-B) + 0.0001*B )*r(190))*frac(7,s)+
     &               ro*((0.0541*(1-B) + 0.0061*B )*r(193)+
     &               (0.0541*(1-B) + 0.0061*B )*r(194)+
     &               (0.0541*(1-B) + 0.0061*B )*r(195)+
     &               (0.0541*(1-B) + 0.0061*B )*r(196)+
     &               (0.0170 )*r(193)+(0.0170 )*r(194)+
     &               (0.0170 )*r(195)+(0.0170 )*r(196))*frac(3,s)))/conv
        !CBS3
        n=39
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(231)*frac(40,s)+
     &               (-1.000)*r(230)*frac(n,s)+
     &               ro*((0.0076*(1-B) + 0.0057*B )*r(189)+
     &               (0.0076*(1-B) + 0.0057*B )*r(190))*frac(7,s)+
     &               ro*((0.0463*(1-B) + 0.0613*B )*r(193)+
     &               (0.0463*(1-B) + 0.0613*B )*r(194)+
     &               (0.0463*(1-B) + 0.0613*B )*r(195)+
     &               (0.0463*(1-B) + 0.0613*B )*r(196)+
     &               (0.0340 )*r(193)+(0.0340 )*r(194)+
     &               (0.0340 )*r(195)+(0.0340 )*r(196))*frac(3,s)))/conv
        !CBS4
        n=40
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(231)*frac(41,s)+
     &               (-1.000)*r(230)*frac(n,s)+
     &               ro*((0.0038*(1-B) + 0.0038*B )*r(189)+
     &               (0.0038*(1-B) + 0.0038*B )*r(190))*frac(7,s)+
     &               ro*((0.1809*(1-B) + 0.1014*B )*r(193)+
     &               (0.1809*(1-B) + 0.1014*B )*r(194)+
     &               (0.1809*(1-B) + 0.1014*B )*r(195)+
     &               (0.1809*(1-B) + 0.1014*B )*r(196)+
     &               (0.1700 )*r(193)+(0.1700 )*r(194)+
     &               (0.1700 )*r(195)+(0.1700 )*r(196))*frac(3,s)))/conv
        !CBS5
        n=41
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               (-1.000)*r(231)*frac(n,s)+
     &               ro*((0.3065*(1-B) + 0.2558*B )*r(193)+
     &               (0.3065*(1-B) + 0.2558*B )*r(194)+
     &               (0.3065*(1-B) + 0.2558*B )*r(195)+
     &               (0.3065*(1-B) + 0.2558*B )*r(196)+
     &               (0.2040 )*r(193)+(0.2040 )*r(194)+
     &               (0.2040 )*r(195)+(0.2040 )*r(196))*frac(3,s)))/conv
        !CAS1
        n=42
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(234)*frac(43,s)+
     &               (-0.000)*r(233)*frac(n,s)))/conv
        !CAS2
        n=43
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(235)*frac(44,s)+
     &               (-1.000)*r(234)*frac(n,s)+
     &               ro*((0.0048*(1-B) + 0.0048*B )*r(202))*frac(6,s)+   !ARO1
     &               ro*((0.0380*(1-B) + 0.0008*B )*r(203))*frac(12,s)+  !ARO2
     &               ro*((0.0012*(1-B) + 0.0002*B )*r(204)+
     &               (0.0012*(1-B) + 0.0002*B )*r(205)+
     &               (0.0012*(1-B) + 0.0002*B )*r(206)+
     &               (0.0012*(1-B) + 0.0002*B )*r(207))*frac(1,s)+
     &               ro*((0.0079*(1-B) + 0.0011*B )*r(208)+
     &               (0.0079*(1-B) + 0.0011*B )*r(209)+
     &               (0.0079*(1-B) + 0.0011*B )*r(210)+
     &               (0.0079*(1-B) + 0.0011*B )*r(211))*frac(2,s)))/conv
        !CAS3
        n=44
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(236)*frac(45,s)+
     &               (-1.000)*r(235)*frac(n,s)+
     &               ro*((0.1148*(1-B) + 0.1148*B )*r(202))*frac(6,s)+
     &               ro*((0.1519*(1-B) + 0.0987*B )*r(203))*frac(12,s)+
     &               ro*((0.0025*(1-B) + 0.0012*B )*r(204)+
     &               (0.0025*(1-B) + 0.0012*B )*r(205)+
     &               (0.0025*(1-B) + 0.0012*B )*r(206)+
     &               (0.0025*(1-B) + 0.0012*B )*r(207))*frac(1,s)+
     &               ro*((0.0153*(1-B) + 0.0090*B )*r(208)+
     &               (0.0153*(1-B) + 0.0090*B )*r(209)+
     &               (0.0153*(1-B) + 0.0090*B )*r(210)+
     &               (0.0153*(1-B) + 0.0090*B )*r(211))*frac(2,s)+
     &               ro*((0.0244*(1-B)+0.0122*B)*r(200))*frac(4,s)+
     &               ro*((0.1426*(1-B)+0.0713*B)*r(201))*frac(5,s)))/conv
        !CAS4
        n=45
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(237)*frac(46,s)+
     &               (-1.000)*r(236)*frac(n,s)+
     &               ro*((0.3349*(1-B) + 0.2153*B )*r(202))*frac(6,s)+
     &               ro*((0.1899*(1-B) + 0.1519*B )*r(203))*frac(12,s)+
     &               ro*((0.0165*(1-B) + 0.0103*B )*r(204)+
     &               (0.0165*(1-B) + 0.0103*B )*r(205)+
     &               (0.0165*(1-B) + 0.0103*B )*r(206)+
     &               (0.0165*(1-B) + 0.0103*B )*r(207))*frac(1,s)+
     &               ro*((0.0453*(1-B) + 0.0290*B )*r(208)+
     &               (0.0453*(1-B) + 0.0290*B )*r(209)+
     &               (0.0453*(1-B) + 0.0290*B )*r(210)+
     &               (0.0453*(1-B) + 0.0290*B )*r(211))*frac(2,s)))/conv
        !CAS5
        n=46
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               (-1.000)*r(237)*frac(n,s)+
     &               ro*((0.4306*(1-B) + 0.3349*B )*r(202))*frac(6,s)+
     &               ro*((0.2658*(1-B) + 0.2203*B )*r(203))*frac(12,s)+
     &               ro*((0.0617*(1-B) + 0.0411*B )*r(204)+
     &               (0.0617*(1-B) + 0.0411*B )*r(205)+
     &               (0.0617*(1-B) + 0.0411*B )*r(206)+
     &               (0.0617*(1-B) + 0.0411*B )*r(207))*frac(1,s)+
     &               ro*((0.1318*(1-B) + 0.0949*B )*r(208)+
     &               (0.1318*(1-B) + 0.0949*B )*r(209)+
     &               (0.1318*(1-B) + 0.0949*B )*r(210)+
     &               (0.1318*(1-B) + 0.0949*B )*r(211))*frac(2,s)))/conv
        !OLE1 (only consumed)
        n=1
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !OLE2 (only consumed)
        n=2
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !TERP (only consumed)
        n=3
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !ALK4 (only consumed)
        n=4
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !ALK5 (only consumed)
        n=5
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !ARO1 (only consumed)
        n=6
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !ARO2 (only consumed)
        n=12
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !ISOP (only consumed)
        n=7
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv

        !BLANK GAS INDICES (Change if you add species!!!)
        do n=47,49
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
          Appconc(loc)=0.0
        enddo

      enddo
c-------------------------------------------------------------------
c------------------------CHECKS-------------------------------------
c-------------------------------------------------------------------
      do spc=1,sa_num_gas
c
c-----Check total & Negative Concentrations
        tot = 0.0
        do n=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(n-1)
          loc2 = Appmaprev(spc)
          if (Appconc(loc).lt.0) then
             write(6,*) 'Negative Concentrations in Appgas: ',i,j,k,n,
     &                   spc,Appconc(loc)
c            if ( (Appconc(loc).gt.-1e-5 .and. spc.eq.14 ).or.
c     &           (Appconc(loc).gt.-1e-7 .and. spc.eq.18 )    ) then
c              !The value is just barely less than 0. Make it zero. Usually for PAN2
c              !Add the negative value to the source with the most mass in it
c              ismax = 1
c              do is = 2,Appnum+3
c                loc3=i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(spc-1)+nx*ny*nz*MXTRK*(ismax-1)
c                loc4=i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(spc-1)+nx*ny*nz*MXTRK*(is-1)
c                if (Appconc(loc4).gt.Appconc(loc3)) ismax = is
c              enddo
c              loc3=i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(spc-1)+nx*ny*nz*MXTRK*(ismax-1)
c              !print *,'Appgas: Negative Concentrations. ismax=',ismax,' spc=',spc
c              !print *,'Appgas: Negative Concentrations. Appconc(loc3)=',Appconc(loc3)
c              Appconc(loc3) = Appconc(loc3) + Appconc(loc)
c              Appconc(loc) = bdnl(Appmaprev(spc))/conv
c              !print *,'Appgas: Negative Concentrations. Appconc(loc3)=',Appconc(loc3)
c            endif
          endif
        enddo
        do n=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &          nx*ny*nz*MXTRK*(n-1)
          tot = tot+Appconc(loc)
        enddo
c
        if (abs(tot-new(loc2)/conv).gt.0.01*MIN(new(loc2)/conv,tot)) then
          if (new(Appmaprev(spc)).eq.bdnl(Appmaprev(spc)) .or.
     &        tot.lt.100*bdnl(Appmaprev(spc))/conv ) then
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
