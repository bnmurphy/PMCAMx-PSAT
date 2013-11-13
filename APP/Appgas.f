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
            
            if (NO2case.eq.0) then
              do ict=1,Appnum+3
                loc3 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &                 nx*ny*nz*MXTRK*(ict-1)
                Appconc(loc3) = Appconc(loc3)*(old(Appmaprev(n))
     &                        /conv/total(n))
              enddo
            elseif (NO2case.eq.1) then
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
     &            0.05*MIN(old(Appmaprev(n))/conv,total(n))) 
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
      do n=1,sa_num_gas
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
             write(6,*) 'Appgas: fraction < 0',(frac(n,is),is=1,Appnum+3),n,s,i,j,k
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
c
c     NOx calculations: NO, NO2, N2O5 and NO3
      NOxTotal = old(Appmaprev(1))+old(Appmaprev(2))+old(Appmaprev(4))+
     &           NO3conc
c     Write gain and loss terms for non-NOx species from the
c     perspective of gain and loss to the NOx family
      gPBZN = +( 1.000)*r( 90)
      lPBZN = +( 1.000)*r( 91)
      gMPAN = +( 1.000)*r(102)
      lMPAN = +( 1.000)*r(103)
      gPAN2 = +( 1.000)*r( 79)
      lPAN2 = +( 1.000)*r( 80)
      gPAN  = +( 1.000)*r( 69)
      lPAN  = +( 1.000)*r( 70)
      gNPHE = +( 1.000)*r(117)+( 1.000)*r(121)+( 1.000)*r(122)
      gHNO3 = +( 2.000)*r( 13)+( 1.000)*r( 25)+( 0.200)*r( 39)
     &        +( 1.000)*r(129)+( 1.000)*r(132)+( 1.000)*r(135)
     &        +( 1.000)*r(148)+( 1.000)*r(151)+( 1.000)*r(154)
     &        +( 1.000)*r(156)+( 1.000)*r(157)+( 1.000)*r(160)
     &        +( 0.500)*r(163)+( 0.150)*r(172)
      lHNO3 = +( 1.000)*r( 27)+( 1.000)*r( 28)
      gRNO3 = +( 1.000)*r(115)+( 0.572)*r(172)+( 0.310)*r(176)
     &        +( 0.813)*r(191)+( 0.276)*r(195)+( 0.511)*r(206)
     &        +( 0.321)*r(210)+( 1.000)*r( 62)
      lRNO3 = +( 1.000)*r(176)+( 1.000)*r(177)
      
      gHONO = +( 1.000)*r( 21)
      lHONO = +( 1.000)*r( 22)+( 1.000)*r( 23)+( 1.000)*r( 24)
      gHNO4 = +( 1.000)*r( 32)
      lHNO4 = +( 1.000)*r( 33)+( 1.000)*r( 34)+( 1.000)*r( 35)
      gXN   = +( 2.000)*r(120)+( 0.500)*r(163)+( 0.278)*r(172)
     &        +( 0.352)*r(176)+( 1.000)*r(187)+( 0.250)*r(195)
     &        +( 0.489)*r(206)+( 0.288)*r(210)
      cycMPAN = MIN(lMPAN,gMPAN)
      cycPAN2 = MIN(lPAN2,gPAN2)
      cycPAN  = MIN(lPAN,gPAN)
      cycHONO = MIN(lHONO,gHONO)
      dHONO = gHONO-lHONO
      dPAN  = gPAN-lPAN
      dPAN2 = gPAN2-lPAN2
      dMPAN = gMPAN-lMPAN

      !KRISTINA'S CODE
      !gPBZN = +( 1.000)*r( 90)
      !lPBZN = +( 1.000)*r( 91)
      !gMPAN = +( 1.000)*r(102)
      !lMPAN = +( 1.000)*r(103)
      !gPAN2 = +( 1.000)*r( 79)
      !lPAN2 = +( 1.000)*r( 80)
      !gPAN  = +( 1.000)*r( 69)
      !lPAN  = +( 1.000)*r( 70)
      !gNPHE = +( 1.000)*r(117)+( 1.000)*r(121)+( 1.000)*r(122)
      !gHNO3 = +( 2.000)*r( 13)+( 1.000)*r( 25)+( 0.200)*r( 39)
      !&        +( 1.000)*r(129)+( 1.000)*r(132)+( 1.000)*r(135)
      !&        +( 1.000)*r(148)+( 1.000)*r(151)+( 1.000)*r(154)
      !&        +( 1.000)*r(156)+( 1.000)*r(157)+( 1.000)*r(160)
      !&        +( 0.500)*r(163)+( 0.150)*r(172)
      !lHNO3 = +( 1.000)*r( 27)+( 1.000)*r( 28)
      !gHONO = +( 1.000)*r( 21)
      !lHONO = +( 1.000)*r( 22)+( 1.000)*r( 23)+( 1.000)*r( 24)
      !gHNO4 = +( 1.000)*r( 32)
      !lHNO4 = +( 1.000)*r( 33)+( 1.000)*r( 34)+( 1.000)*r( 35)
      !gXN   = +( 2.000)*r(120)+( 0.500)*r(163)+( 0.278)*r(172)
      !&        +( 0.352)*r(176)+( 1.000)*r(187)+( 0.250)*r(195)
      !&        +( 0.489)*r(206)+( 0.288)*r(210)
      
      if (dPAN.gt.0) then
         ngPAN = dPAN
         nlPAN = 0
      else
         ngPAN = 0
         nlPAN = -1*dPAN
      endif
      if (dPAN2.gt.0) then
         ngPAN2 = dPAN2
         nlPAN2 = 0
      else
         ngPAN2 = 0
         nlPAN2 = -1*dPAN2
      endif
      if (dMPAN.gt.0) then
         ngMPAN = dMPAN
         nlMPAN = 0
      else
         ngMPAN = 0
         nlMPAN = -1*dMPAN
      endif
      if (dHONO.gt.0) then
         ngHONO = dHONO
         nlHONO = 0
      else
         ngHONO = 0
         nlHONO = -1*dHONO
      endif

      !Added Gain and Loss of RNO3 to these calculations: BNM 11-14-12
      NOxRemain = NOxTotal + (-1*gPBZN-gHNO3-ngHONO-ngPAN-ngPAN2-ngMPAN
     &            -gHNO4-gNPHE-gXN-gRNO3-0.5*cycPAN-0.5*cycPAN2-0.5*cycHONO
     &            -0.5*cycMPAN)*dt
      NOxEnd = NOxTotal + (-1*gPBZN-gHNO3-dHONO-dPAN-dPAN2-dMPAN
     &            -gHNO4-gNPHE-gXN-gRNO3+lPBZN+lHNO3+lHNO4+lRNO3)*dt
      modelTot = new(Appmaprev(1))+new(Appmaprev(2))+new(Appmaprev(4))+
     &           NO3new

      !KRISTINA
      ! NOxRemain = NOxTotal + (-1*gPBZN-gHNO3-ngHONO-ngPAN-ngPAN2-ngMPAN
      !&            -gHNO4-gNPHE-gXN-0.5*cycPAN-0.5*cycPAN2-0.5*cycHONO
      !&            -0.5*cycMPAN)*dt
      ! NOxEnd = NOxTotal + (-1*gPBZN-gHNO3-dHONO-dPAN-dPAN2-dMPAN
      !&            -gHNO4-gNPHE-gXN+lPBZN+lHNO3+lHNO4)*dt
      ! modelTot = new(Appmaprev(1))+new(Appmaprev(2))+new(Appmaprev(4))+
      !&           NO3new
      if (abs(modelTot-NOxEnd).gt.0.05*MIN(modelTot,NOxEnd)) then
         write(6,*) 'Big diff in end concentrations',modelTot,NOxEnd,
     &              i,j,k,NO3new,NO3conc
      endif
c
      NOxTotNew=0.0
      do s=1,Appnum+3
         if (NOxRemain.gt.0) then
         NOxSplit(s)=NOxTotal*frac(2,s)-dt*(gNPHE+gXN+gPBZN+gHNO4+gHNO3+gRNO3+
     &               0.5*cycHONO+ngHONO+0.5*cycPAN+0.5*cycPAN2+
     &               0.5*cycMPAN+ngMPAN+ngPAN+ngPAN2)*frac(2,s)+
     &               dt*(lHNO3*frac(11,s)+lHNO4*frac(20,s)+
     &               lPBZN*frac(13,s)+(nlHONO+0.5*cycHONO)*frac(10,s)+
     &               (nlPAN+0.5*cycPAN)*frac(3,s)+
     &               (nlMPAN+0.5*cycMPAN)*frac(18,s)+
     &               (nlPAN2+0.5*cycPAN2)*frac(14,s)+lRNO3*frac(57,s))
         else
            write(6,*) 'NOx used up ...',i,j,k,NOxRemain,NOxTotal
           NOxSplit(s)=(NOxRemain+lPAN2+lMPAN+lPBZN+lHNO4+lPAN+lHNO3+
     &                  lHONO+lRNO3) * (lPAN2*frac(14,s)+lMPAN*frac(18,s)
     &                  +lPAN*frac(3,s)+lHNO3*frac(11,s)+lHONO*
     &                   frac(10,s)+lPBZN*frac(13,s)+lHNO4*frac(20,s)
     &                  +lRNO3*frac(57,s))/
     &                  (lPAN2+lMPAN+lPBZN+lHNO4+lPAN+lHNO3+lHONO+lRNO3)
         endif
c
         NOxTotNew = NOxTotNew+NOxSplit(s)
         if (NOxSplit(s).lt.0) then
            write(6,*) "NOxSplit is less than 0"
            stop
         endif
      enddo
      check = 0.0
      do s=1,Appnum+3
         NOxFracs(s)=NOxSplit(s)/NOxTotNew
         check = check+NOxFracs(s)
      enddo
c
      if (abs(check-1).gt.0.001) then
         write(6,*) 'NOx Fraction not equal to 1'
         do s=1,Appnum+3
            write(6,*) NOxFracs(s)
         enddo
         stop
      endif

c-------B is the branching ratio (NOx dependent SOA yields)
c        B = 0.0
        B =(r( 46)+r( 48)+r( 51)+r( 53)+r( 56)+r( 58)+r( 62)+r( 65)) /
     &     ( r( 46)+r( 48)+r( 51)+r( 53)+r( 56)+r( 58)+r( 62)+r( 65) +
     &     r( 47)+r( 49)+r( 50)+r( 52)+r( 54)+r( 55)+r( 57)+r( 59)+
     &     r( 60)+r( 61)+r( 63)+r( 64)+r( 66)+r( 67)+r( 68) )
      ro = 1.4
c
      do s=1,Appnum+3
        !NO
        n=1
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc) = new(loc2)*NOxFracs(s)/conv

        !NO2
        n=2
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc) = new(loc2)*NOxFracs(s)/conv
        
        !PAN
        n=3
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)

        
        Appconc(loc) = (old(loc2)*frac(n,s)
     &                 + dt*(( 1.000)*r( 69)*frac(2,s)
     &                 +     (-1.000)*r( 70)*frac(n,s) ))/conv
        !print *,'Appgas: i=',i,' j=',j,' k=',k,' s=',s
        !if (Appconc(loc).lt.0.and.abs(Appconc(loc)).lt.1e-9) then
        !  print *,'  Appconc= ',Appconc(loc),'  bdnl=',bdnl(loc2),' conv=',conv
        !  Appconc(loc) = bdnl(loc2)/conv
        !  print *,'  Appconc reassigned to: ',Appconc(loc),'  bdnl=',bdnl(loc2),' conv=',conv
        !endif
      !if (i.eq.105.and.j.eq.2) then
      ! print *,'Appgas: PAN. s=',s,' k=',k, ' frac(2,s) = ',frac(2,s)
      ! print *,'   frac(n,s)=',frac(n,s),'  old(PAN)=',old(loc2)
      ! print *,'   Appconc(PAN) =', Appconc(loc),' NOxfracs(s)=',NOxFracs(s)
      !endif
     
        !NXOY
        n=4
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc) = new(loc2)*NOXFracs(s)/conv

        !HONO
        n=10
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc) = (old(loc2)*frac(n,s)
     &                 +dt*( ( 1.000)*r( 21)*frac(2,s)
     &                      +((-1.000)*r( 22)+(-1.00)*r(23)
     &                      + (-1.000)*r( 24) )*frac(n,s)  ) )/conv
        
        !HNO3
        n=11
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc) = (old(loc2)*frac(n,s)
     &                 +dt*( ( (2.00)*r(13)+(1.00)*r(25)+(0.20)*r(39)
     &                        +(1.00)*r(129)+(1.00)*r(132)+(1.00)*r(135)
     &                        +(1.00)*r(148)+(1.00)*r(151)+(1.00)*r(154)
     &                        +(1.00)*r(156)+(1.00)*r(157)+(1.00)*r(160)
     &                        +(0.50)*r(163)+(0.15)*r(172) )*frac(2,s)
     &                       +((-1.00)*r(27)+(-1.00)*r(28))*frac(n,s)   ))/conv
      
        !PBZN
        n=13
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc) = (old(loc2)*frac(n,s)
     &                 +dt*( (1.00)*r(90)*frac(2,s) 
     &                      +(-1.0)*r(91)*frac(n,s) ))/conv

        !PAN2
        n=14
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc) = (old(loc2)*frac(n,s)
     &                 +dt*( (1.00)*r(79)*frac(2,s)
     &                      +(-1.0)*r(80)*frac(n,s)))/conv

        !MPAN
        n=18
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc) = (old(loc2)*frac(n,s)
     &                 +dt*( (1.00)*r(102)*frac(2,s)
     &                      +(-1.0)*r(103)*frac(n,s) ))/conv
        !if (Appconc(loc).lt.0.and.abs(Appconc(loc)).lt.1e-9) then
        !  Appconc(loc) = bdnl(loc2)/conv
        !endif
      
        !HNO4
        n=20
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc) = (old(loc2)*frac(n,s)
     &                 +dt*( (1.00)*r(32)*frac(2,s)
     &                      +((-1.0)*r(33)+(-1.0)*r(34)
     &                      + (-1.0)*r(35))*frac(n,s)  ))/conv

        !RNO3 - Added by BNM 11-14-12
        n=57
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) + 
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc) = (old(loc2)*frac(n,s)
     &                 +dt*((+( 1.000)*r( 62)
     &                       +( 1.000)*r(115)+( 0.572)*r(172)+( 0.310)*r(176)
     &                       +( 0.813)*r(191)+( 0.276)*r(195)+( 0.511)*r(206)
     &                       +( 0.321)*r(210))*frac(2,s) 
     &                       +((-1.0)*r(176)+(-1.0)*r(177))*frac(n,s) ))/conv


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
        n=22
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO2 (only consumed)
        n=23
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO3 (only consumed)
        n=24
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO4 (only consumed)
        n=25
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO5 (only consumed)
        n=26
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO6 (only consumed)
        n=27
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO7 (only consumed)
        n=28
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !CPO8 (only consumed)
        n=29
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !COO1
        n=30
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(213)*frac(23,s)+
     &               ( 1.075)*r(221)*frac(31,s)))/conv
        !COO2
        n=31
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(214)*frac(24,s)+
     &               ( 1.075)*r(222)*frac(32,s)+
     &               (-1.000)*r(221)*frac(31,s)))/conv
        !COO3
        n=32
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(215)*frac(25,s)+
     &               ( 1.075)*r(223)*frac(33,s)+
     &               (-1.000)*r(222)*frac(32,s)))/conv
        !COO4
        n=33
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(216)*frac(26,s)+
     &               ( 1.075)*r(224)*frac(34,s)+
     &               (-1.000)*r(223)*frac(33,s)))/conv
        !COO5
        n=34
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(217)*frac(27,s)+
     &               ( 1.075)*r(225)*frac(35,s)+
     &               (-1.000)*r(224)*frac(34,s)))/conv
        !COO6
        n=35
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(218)*frac(28,s)+
     &               ( 1.075)*r(226)*frac(36,s)+
     &               (-1.000)*r(225)*frac(35,s)))/conv
        !COO7
        n=36
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(219)*frac(29,s)+
     &               ( 1.075)*r(227)*frac(37,s)+
     &               (-1.000)*r(226)*frac(36,s)))/conv
        !COO8
        n=37
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               (-1.000)*r(227)*frac(37,s)))/conv
        !CNS1
        n=38
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(239)*frac(39,s)+
     &               (-0.000)*r(238)*frac(38,s)))/conv
        !CNS2
        n=39
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(240)*frac(40,s)+
     &               (-1.000)*r(239)*frac(39,s)))/conv
        !CNS3
        n=40
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(241)*frac(41,s)+
     &               (-1.000)*r(240)*frac(40,s)))/conv
        !CNS4
        n=41
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(242)*frac(42,s)+
     &               (-1.000)*r(241)*frac(41,s)))/conv
        !CNS5
        n=42
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(243)*frac(43,s)+
     &               (-1.000)*r(242)*frac(42,s)))/conv
        !CNS6
        n=43
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(244)*frac(44,s)+
     &               (-1.000)*r(243)*frac(43,s)))/conv
        !CNS7
        n=44
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(245)*frac(45,s)+
     &               (-1.000)*r(244)*frac(44,s)))/conv
        !CNS8
        n=45
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               (-1.000)*r(245)*frac(45,s)))/conv
        !CBS1
        n=46
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(229)*frac(47,s)+
     &               (-0.000)*r(228)*frac(46,s)))/conv
        !CBS2
        n=47
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(230)*frac(48,s)+
     &               (-1.000)*r(229)*frac(47,s)+
     &               ro*((0.0023*(1-B) + 0.0001*B )*r(189)+
     &               (0.0023*(1-B) + 0.0001*B )*r(190))*frac(12,s)+
     &               ro*((0.0541*(1-B) + 0.0061*B )*r(193)+
     &               (0.0541*(1-B) + 0.0061*B )*r(194)+
     &               (0.0541*(1-B) + 0.0061*B )*r(195)+
     &               (0.0541*(1-B) + 0.0061*B )*r(196)+
     &               (0.0170 )*r(193)+(0.0170 )*r(194)+
     &               (0.0170 )*r(195)+(0.0170 )*r(196))*frac(7,s)))/conv
        !CBS3
        n=48
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(231)*frac(49,s)+
     &               (-1.000)*r(230)*frac(48,s)+
     &               ro*((0.0076*(1-B) + 0.0057*B )*r(189)+
     &               (0.0076*(1-B) + 0.0057*B )*r(190))*frac(12,s)+
     &               ro*((0.0463*(1-B) + 0.0613*B )*r(193)+
     &               (0.0463*(1-B) + 0.0613*B )*r(194)+
     &               (0.0463*(1-B) + 0.0613*B )*r(195)+
     &               (0.0463*(1-B) + 0.0613*B )*r(196)+
     &               (0.0340 )*r(193)+(0.0340 )*r(194)+
     &               (0.0340 )*r(195)+(0.0340 )*r(196))*frac(7,s)))/conv
        !CBS4
        n=49
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(231)*frac(50,s)+
     &               (-1.000)*r(230)*frac(49,s)+
     &               ro*((0.0038*(1-B) + 0.0038*B )*r(189)+
     &               (0.0038*(1-B) + 0.0038*B )*r(190))*frac(12,s)+
     &               ro*((0.1809*(1-B) + 0.1014*B )*r(193)+
     &               (0.1809*(1-B) + 0.1014*B )*r(194)+
     &               (0.1809*(1-B) + 0.1014*B )*r(195)+
     &               (0.1809*(1-B) + 0.1014*B )*r(196)+
     &               (0.1700 )*r(193)+(0.1700 )*r(194)+
     &               (0.1700 )*r(195)+(0.1700 )*r(196))*frac(7,s)))/conv
        !CBS5
        n=50
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               (-1.000)*r(231)*frac(50,s)+
     &               ro*((0.3065*(1-B) + 0.2558*B )*r(193)+
     &               (0.3065*(1-B) + 0.2558*B )*r(194)+
     &               (0.3065*(1-B) + 0.2558*B )*r(195)+
     &               (0.3065*(1-B) + 0.2558*B )*r(196)+
     &               (0.2040 )*r(193)+(0.2040 )*r(194)+
     &               (0.2040 )*r(195)+(0.2040 )*r(196))*frac(7,s)))/conv
        !CAS1
        n=51
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(234)*frac(52,s)+
     &               (-0.000)*r(233)*frac(51,s)))/conv
        !CAS2
        n=52
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(235)*frac(53,s)+
     &               (-1.000)*r(234)*frac(52,s)+
     &               ro*((0.0048*(1-B) + 0.0048*B )*r(202))*frac(9,s)+   !ARO1
     &               ro*((0.0380*(1-B) + 0.0008*B )*r(203))*frac(21,s)+  !ARO2
     &               ro*((0.0012*(1-B) + 0.0002*B )*r(204)+
     &               (0.0012*(1-B) + 0.0002*B )*r(205)+
     &               (0.0012*(1-B) + 0.0002*B )*r(206)+
     &               (0.0012*(1-B) + 0.0002*B )*r(207))*frac(5,s)+
     &               ro*((0.0079*(1-B) + 0.0011*B )*r(208)+
     &               (0.0079*(1-B) + 0.0011*B )*r(209)+
     &               (0.0079*(1-B) + 0.0011*B )*r(210)+
     &               (0.0079*(1-B) + 0.0011*B )*r(211))*frac(6,s)))/conv
        !CAS3
        n=53
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(236)*frac(54,s)+
     &               (-1.000)*r(235)*frac(53,s)+
     &               ro*((0.1148*(1-B) + 0.1148*B )*r(202))*frac(9,s)+
     &               ro*((0.1519*(1-B) + 0.0987*B )*r(203))*frac(21,s)+
     &               ro*((0.0025*(1-B) + 0.0012*B )*r(204)+
     &               (0.0025*(1-B) + 0.0012*B )*r(205)+
     &               (0.0025*(1-B) + 0.0012*B )*r(206)+
     &               (0.0025*(1-B) + 0.0012*B )*r(207))*frac(5,s)+
     &               ro*((0.0153*(1-B) + 0.0090*B )*r(208)+
     &               (0.0153*(1-B) + 0.0090*B )*r(209)+
     &               (0.0153*(1-B) + 0.0090*B )*r(210)+
     &               (0.0153*(1-B) + 0.0090*B )*r(211))*frac(6,s)+
     &               ro*((0.0244*(1-B)+0.0122*B)*r(200))*frac(56,s)+
     &               ro*((0.1426*(1-B)+0.0713*B)*r(201))*frac(8,s)))/conv
        !CAS4
        n=54
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               ( 1.075)*r(237)*frac(55,s)+
     &               (-1.000)*r(236)*frac(54,s)+
     &               ro*((0.3349*(1-B) + 0.2153*B )*r(202))*frac(9,s)+
     &               ro*((0.1899*(1-B) + 0.1519*B )*r(203))*frac(21,s)+
     &               ro*((0.0165*(1-B) + 0.0103*B )*r(204)+
     &               (0.0165*(1-B) + 0.0103*B )*r(205)+
     &               (0.0165*(1-B) + 0.0103*B )*r(206)+
     &               (0.0165*(1-B) + 0.0103*B )*r(207))*frac(5,s)+
     &               ro*((0.0453*(1-B) + 0.0290*B )*r(208)+
     &               (0.0453*(1-B) + 0.0290*B )*r(209)+
     &               (0.0453*(1-B) + 0.0290*B )*r(210)+
     &               (0.0453*(1-B) + 0.0290*B )*r(211))*frac(6,s)))/conv
        !CAS5
        n=55
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=(old(loc2)*frac(n,s)+dt*(
     &               (-1.000)*r(237)*frac(55,s)+
     &               ro*((0.4306*(1-B) + 0.3349*B )*r(202))*frac(9,s)+
     &               ro*((0.2658*(1-B) + 0.2203*B )*r(203))*frac(21,s)+
     &               ro*((0.0617*(1-B) + 0.0411*B )*r(204)+
     &               (0.0617*(1-B) + 0.0411*B )*r(205)+
     &               (0.0617*(1-B) + 0.0411*B )*r(206)+
     &               (0.0617*(1-B) + 0.0411*B )*r(207))*frac(5,s)+
     &               ro*((0.1318*(1-B) + 0.0949*B )*r(208)+
     &               (0.1318*(1-B) + 0.0949*B )*r(209)+
     &               (0.1318*(1-B) + 0.0949*B )*r(210)+
     &               (0.1318*(1-B) + 0.0949*B )*r(211))*frac(6,s)))/conv
        !OLE1 (only consumed)
        n=5
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !OLE2 (only consumed)
        n=6
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !TERP (only consumed)
        n=7
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !ALK4 (only consumed)
        n=56
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !ALK5 (only consumed)
        n=8
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !ARO1 (only consumed)
        n=9
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(n-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        loc2 = Appmaprev(n)
        Appconc(loc)=new(loc2)*frac(n,s)/conv
        !ARO2 (only consumed)
        n=21
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

        !BLANK GAS INDICES (Change if you add species!!!)
        do n=58,64
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
            if ( (Appconc(loc).gt.-1e-5 .and. spc.eq.14 ).or.
     &           (Appconc(loc).gt.-1e-7 .and. spc.eq.18 )    ) then
              !The value is just barely less than 0. Make it zero. Usually for PAN2
              !Add the negative value to the source with the most mass in it
              ismax = 1
              do is = 2,Appnum+3
                loc3=i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(spc-1)+nx*ny*nz*MXTRK*(ismax-1)
                loc4=i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(spc-1)+nx*ny*nz*MXTRK*(is-1)
                if (Appconc(loc4).gt.Appconc(loc3)) ismax = is
              enddo
              loc3=i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz*(spc-1)+nx*ny*nz*MXTRK*(ismax-1)
              !print *,'Appgas: Negative Concentrations. ismax=',ismax,' spc=',spc
              !print *,'Appgas: Negative Concentrations. Appconc(loc3)=',Appconc(loc3)
              Appconc(loc3) = Appconc(loc3) + Appconc(loc)
              Appconc(loc) = bdnl(Appmaprev(spc))/conv
              !print *,'Appgas: Negative Concentrations. Appconc(loc3)=',Appconc(loc3)
            endif
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
