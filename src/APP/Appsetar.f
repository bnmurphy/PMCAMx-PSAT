       subroutine Appsetar
c       subroutine Appsetar(areaemis,nareasp)
c
c     This subroutine sets the source specific area emissions.
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
c      real areaemis(MXCOL1,MXROW1,MXSPEC)
c      real Appemis(MXCOL1,MXROW1,MXTRK)
      character*4 name1(10)
      integer s,test,nx,ny,nz
      logical match
c
      nx = MXCOL1
      ny = MXROW1
      nz = MXLAY1
c
c-----Initialize the fractions array
c
        do i=1,MXCOL1
          do j=1,MXROW1
            do k=1,MXTRK
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
                if (date.eq.1193) then
                  AppemfracA(i,j,Appmap(ll),1) = 1.0
                  do n = 2,Appnum+1
                    AppemfracA(i,j,Appmap(ll),n) = 0.0
                  enddo
                endif 
              endif
            enddo
          enddo
        enddo
c
c-----Assigning remaining fraction to "others" (4th dimension = 1)
c
        do i=1,MXCOL1
          do j=1,MXROW1
            do k=1,MXTRK
              do s=2,Appnum+1
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
            do k=1,MXTRK
              total = 0.0
              do s=1,Appnum+1
                total = total + AppemfracA(i,j,k,s)
              enddo
              if (abs(1-total).gt.0.01) then
                write(6,*) 'ERROR in Appsetar.f:'
                write(6,*) 'sum of fractions not equal to 1'
                write(6,*) 'total: ',total
                write(6,*) 'i,j,k: ',i,j,k
                do s=1,Appnum+1
                  write(6,*) s,AppemfracA(i,j,k,s)
                enddo
                stop
              endif
            enddo
          enddo
        enddo
c
      end

              
