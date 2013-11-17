      subroutine Appave ()
c
c     This rountine averages the concentrations for each hour.
c
c     DEVELOPED BY:
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     08/20/2006
c     
c     CALLED BY:
c     CAMx.f
c
c     ROUTINES CALLED:
c     Appwrt.f
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
      integer i,j,spc,s,loc,loc2,loc3
      real count,convfac
c
      count = amod(time,100.)
c
c-----Each hour
c
      call Appwrt(0)
      if (count.eq.0.0) then
        call Appwrt(1)
        Appct = 1
        do i = 1,MXCOL1
          do j = 1,MXROW1
            do spc = 1,MXTRK
              do s = 1,MXSOUR
                loc = i+MXCOL1*(j-1)+MXCOL1*MXROW1*(spc-1)+
     &                MXCOL1*MXROW1*MXTRK*(s-1)
                loc2 = i+MXCOL1*(j-1)+MXCOL1*MXROW1*MXLAY1*(spc-1)+
     &                MXCOL1*MXROW1*MXLAY1*MXTRK*(s-1)
                loc3 = i+MXCOL1*(j-1)
                if (spc.le.sa_num_gas) then !Change if change PSAT species
                  convfac=densfac*(273./tempk(loc3))*
     &                    (press(loc3)/1013.)    !For Gases, umol/m3->ppm
                  Appavg(loc) = Appconc(loc2)/convfac
c                  if (s.eq.4.and.spc.eq.15.and.Appavg(loc).ne.0)
c     &              write(6,*) i,j,spc,s,Appavg(loc),Appconc(loc)
                else
                  Appavg(loc) = Appconc(loc2)
                endif
              enddo
            enddo
          enddo
        enddo
c
c-----Each time step
c
      else
        do i = 1,MXCOL1
          do j = 1,MXROW1
            do spc = 1,MXTRK
              do s = 1,MXSOUR
                loc = i+MXCOL1*(j-1)+MXCOL1*MXROW1*(spc-1)+
     &                MXCOL1*MXROW1*MXTRK*(s-1)
                loc2 = i+MXCOL1*(j-1)+MXCOL1*MXROW1*MXLAY1*(spc-1)+
     &                MXCOL1*MXROW1*MXLAY1*MXTRK*(s-1)
                loc3 = i+MXCOL1*(j-1)
                if (spc.le.sa_num_gas) then !Change if change PSAT species
                  convfac=densfac*(273./tempk(loc3))*
     &                    (press(loc3)/1013.)    !For Gases, umol/m3->ppm
                  Appavg(loc) = (Appconc(loc2)/convfac+Appavg(loc)*
     &                          Appct)/(Appct+1)
                else
                  Appavg(loc) = Appconc(loc2)
                endif
              enddo
            enddo
          enddo
        enddo
        Appct = Appct + 1
      endif
c
      end
                  
