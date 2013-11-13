       subroutine Appsetar(areaemis,nareasp)
c
c     This subroutine reads in the source specific area emissions.
c  
c     DEVELOPED BY: 
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     04/13/2007
c
c     CALLED BY:
c     startup
c
c     ROUTINES CALLED:
c     none
c
c     VARIABLES (common):
c     spname  - Modeled species names
c     Appnam  - Apportionment species names
c     Appmap  - Mapping values between model and apportionment 
c     MXSPEC  - Maximum number of species
c     MXTRK   - Maximum number of tracked species
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
      real areaemis(MXCOL1,MXROW1,MXSPEC)
      real Appemis(MXCOL1,MXROW1,MXTRK)
      character*4 name1(10)
      integer s
      logical match
c
c-----Initialize the fractions array
c
        do i=1,MXCOL1
          do j=1,MXROW1
            do k=1,MXSPEC
              AppemfracA(i,j,k,1)=1.0
              do s=2,Appnum+1
                AppemfracA(i,j,k,s)=0.0
              enddo
            enddo
          enddo
        enddo
c
c-----Loop over each source
c
c      do s=2,2
c
c-----Calculate fractions for the sources
c
        do i=1,MXCOL1
          do j=1,MXROW1
            do ll=1,MXTRK
              idum=AppemisA(s,ll)
              if (idum.ne.0) then
c       
c               PITTSBURGH
                if (i.ge.64.and.i.le.67.and.
     &               j.ge.49.and.j.le.52) then
                       AppemfracA(i,j,idum,2)=1.0
                endif
c
c               CHICAGO
                if (i.ge.44.and.i.le.48.and.
     &               j.ge.49.and.j.le.54) then
                       AppemfracA(i,j,idum,3)=1.0
                endif
c
c               MID-ATLANTIC
                if (i.ge.73.and.i.le.80.and.
     &               j.ge.42.and.j.le.58) then
                       AppemfracA(i,j,idum,4)=1.0
                endif
                if (i.ge.70.and.i.le.72.and.
     &               j.ge.45.and.j.le.51) then
                       AppemfracA(i,j,idum,4)=1.0
                endif
c
c               DETROIT
                if (i.ge.55.and.i.le.58.and.
     &               j.ge.53.and.j.le.57) then
                       AppemfracA(i,j,idum,5)=1.0
                endif
c
c               LITTLE ROCK
                if (i.ge.35.and.i.le.37.and.
     &               j.ge.29.and.j.le.31) then
                       AppemfracA(i,j,idum,6)=1.0
                endif
c
c               NORTH GULF COAST (including New Orleans)
                if (i.ge.35.and.i.le.40.and.
     &               j.ge.12.and.j.le.16) then
                       AppemfracA(i,j,idum,7)=1.0
                endif
                if (i.ge.41.and.i.le.48.and.
     &               j.ge.10.and.j.le.19) then
                       AppemfracA(i,j,idum,7)=1.0
                endif
                if (i.ge.49.and.i.le.59.and.
     &               j.ge.15.and.j.le.19) then
                       AppemfracA(i,j,idum,7)=1.0
                endif
c
c               BOSTON
                if (i.ge.82.and.i.le.84.and.
     &               j.ge.59.and.j.le.62) then
                       AppemfracA(i,j,idum,8)=1.0
                endif
c
c               ST. LOUIS
                if (i.ge.39.and.i.le.42.and.
     &               j.ge.40.and.j.le.44) then
                       AppemfracA(i,j,idum,9)=1.0
                endif
c
c               ATLANTA
                if (i.ge.55.and.i.le.58.and.
     &               j.ge.26.and.j.le.30) then
                       AppemfracA(i,j,idum,10)=1.0
                endif
c
c               TEXAS GULF COAST (including Houston)
                if (i.ge.22.and.i.le.26.and.
     &               j.ge.1.and.j.le.6) then
                       AppemfracA(i,j,idum,11)=1.0
                endif
                if (i.ge.23.and.i.le.28.and.
     &               j.ge.7.and.j.le.9) then
                       AppemfracA(i,j,idum,11)=1.0
                endif
                if (i.ge.24.and.i.le.31.and.
     &               j.ge.10.and.j.le.11) then
                       AppemfracA(i,j,idum,11)=1.0
                endif
                if (i.ge.25.and.i.le.34.and.
     &               j.ge.12.and.j.le.14) then
                       AppemfracA(i,j,idum,11)=1.0
                endif
                if (i.ge.28.and.i.le.34.and.
     &               j.ge.15.and.j.le.17) then
                       AppemfracA(i,j,idum,11)=1.0
                endif
c
c               OHIO RIVER VALLEY
                if (i.eq.45.and.
     &               j.ge.36.and.j.le.39) then
                       AppemfracA(i,j,idum,12)=1.0
                endif
                if (i.ge.46.and.i.le.47.and.
     &               j.ge.37.and.j.le.41) then
                       AppemfracA(i,j,idum,12)=1.0
                endif
                if (i.ge.48.and.i.le.51.and.
     &               j.ge.39.and.j.le.41) then
                       AppemfracA(i,j,idum,12)=1.0
                endif
                if (i.ge.52.and.i.le.53.and.
     &               j.ge.40.and.j.le.44) then
                       AppemfracA(i,j,idum,12)=1.0
                endif
                if (i.ge.54.and.i.le.59.and.
     &               j.ge.42.and.j.le.45) then
                       AppemfracA(i,j,idum,12)=1.0
                endif
                if (i.ge.60.and.i.le.61.and.
     &               j.ge.42.and.j.le.57) then
                       AppemfracA(i,j,idum,12)=1.0
                endif
                if (i.eq.62.and.
     &               j.ge.44.and.j.le.49) then
                       AppemfracA(i,j,idum,12)=1.0
                endif
                if (i.eq.63.and.
     &               j.ge.45.and.j.le.52) then
                       AppemfracA(i,j,idum,12)=1.0
                endif
c
              endif
            enddo
          enddo
        enddo
c
c-----Assigning remaining fraction to "others" (4th dimension = 1)
c
        do i=1,MXCOL1
          do j=1,MXROW1
            do k=1,MXSPEC
              do s=2,Appnum+1
c                AppemfracA(i,j,k,2)=0.9
                AppemfracA(i,j,k,1)=AppemfracA(i,j,k,1)-
     &                    AppemfracA(i,j,k,s)
              enddo
            enddo
          enddo
        enddo
c
c----Check for sum ------------------------
c
        do i=1,MXCOL1
          do j=1,MXROW1
            do k=1,MXSPEC
              total = 0.0
              do s=1,Appnum+1
                total = total + AppemfracA(i,j,k,s)
              enddo
              if (abs(1-total).gt.0.01) then
                write(6,*) 'ERROR in Appsetar.f:'
                write(6,*) 'sum of fractions not equal to 1'
                write(6,*) 'total: ',total
                do s=1,Appnum+1
                  write(6,*) s,AppemfracA(i,j,k,s)
                enddo
              endif
            enddo
          enddo
        enddo
c
      end

              
