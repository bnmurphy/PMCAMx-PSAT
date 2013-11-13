       subroutine Appsetpt
c       subroutine Appsetpt(pointemis,npointsp)
c
c     This subroutine sets the source specific point emissions.
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
c      real pointemis(120000,MXSPEC)
c      real Appemis(120000,MXTRK)
      real flowrat1(120000),effph1(120000)
      character*4 name1(10)
      integer s
c
c-----Initialize the fractions array
c
        do i=1,120000
          do k=1,MXSPEC
            AppemfracP(i,k,1)=1.0
            do s=2,Appnum+1
              AppemfracP(i,k,s)=0.0
            enddo
          enddo
        enddo
c
c-----Loop over each source
c
c      do s=1,Appnum
c
c-----Calculate fractions for the sources
c
c        do i=1,120000
c          do ll=1,MXTRK
c           idum=AppemisP(s,ll)
c            if (idum.ne.0) then
c                AppemfracP(i,idum,s+1)=0.9
c            endif
c          enddo
c        enddo
c      enddo
c
c-----Assigning remaining fraction to "others" (4th dimension = 1)
c
c        do i=1,120000
c          do k=1,MXSPEC
c            do s=1,Appnum
c              AppemfracP(i,k,2)=0.9
c              AppemfracP(i,k,1)=AppemfracP(i,k,1)-
c     &                    AppemfracP(i,k,s+1)
c            enddo
c          enddo
c        enddo
c
c----Check for 90/10 ------------------------
c      do i=1,120000
c        do k=1,MXSPEC
c          if (abs(9*AppemfracP(i,k,1)-AppemfracP(i,k,2)).gt.0.000001)
c     &        then
c            write(6,*) 'Error in Appsetpt: 90/10 split not correct'
c            write(6,*) AppemfracP(i,k,1),AppemfracP(i,k,2),i,k,Appnum
c            stop
c          endif
c        enddo
c      enddo
c
      end


              
