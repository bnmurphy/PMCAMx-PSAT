      subroutine Appalign(actual,saconc)
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
c     actual  - The model predicted total species concentrations
c     saconc  - The array of source concentrations
c     
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
      real saconc(MXCOL1,MXROW1,MXLAY1,MXTRK*MXSOUR)
      real actual(MXCOL1,MXROW1,MXLAY1,MXSPEC)
      real totals(MXCOL1,MXROW1,MXLAY1,MXTRK)
      integer i,j,k,n,s
c
c     COPY concentrations from saconc to Appconc
c
      do i=1,MXCOL1
        do j=1,MXROW1
          do k=1,MXLAY1
            do n=1,MXTRK
              totals(i,j,k,n) = 0.0
              do s=1,MXSOUR
                loc1=n+MXTRK*(s-1)
                loc2=i + MXCOL1*(j-1) + MXCOL1*MXROW1*(k-1) +
     &               MXCOL1*MXROW1*MXLAY1*(n-1) +
     &               MXCOL1*MXROW1*MXLAY1*MXTRK*(s-1)
                Appconc(loc2) = saconc(i,j,k,loc1)
                totals(i,j,k,n) = totals(i,j,k,n) + Appconc(loc2)
              enddo
              if (totals(i,j,k,n).ne.actual(i,j,k,Appmaprev(n))) then
c                if ((totals(i,j,k,n)-actual(i,j,k,Appmaprev(n))).gt.
c     &               0.99*actual(i,j,k,Appmaprev(n)).and.
c     &               totals(i,j,k,n).gt.1.0E-6) then
c                  write(6,*) 'HUGE difference in totals in Appalign'
c                  write(6,*) 'Location: ',i,j,k,n,s
c                  write(6,*) 'SA concentration total: ',totals(i,j,k,n)
c                  write(6,*) 'Actual concentration: ',
c     &                       actual(i,j,k,Appmaprev(n))
c                  write(6,*) 'Source concentrations: '
c                  do s=1,MXSOUR
c                    loc2=i + MXCOL1*(j-1) + MXCOL1*MXROW1*(k-1) +
c     &                   MXCOL1*MXROW1*MXLAY1*(n-1) +
c     &                   MXCOL1*MXROW1*MXLAY1*MXTRK*(s-1)
c                    write(6,*) s,Appconc(loc2)
c                  enddo
c                  stop
c                endif
                do s=1,MXSOUR
                  loc2=i + MXCOL1*(j-1) + MXCOL1*MXROW1*(k-1) +
     &                 MXCOL1*MXROW1*MXLAY1*(n-1) +
     &                 MXCOL1*MXROW1*MXLAY1*MXTRK*(s-1)
                  Appconc(loc2) = Appconc(loc2) *
     &                            actual(i,j,k,Appmaprev(n))/
     &                            totals(i,j,k,n)
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
      return
      end

