       subroutine Appreadpt(pointemis,npointsp)
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
      real pointemis(120000,MXSPEC)
      real Appemis(120000,MXTRK)
      real flowrat1(120000),effph1(120000)
      character*4 name1(10)
      integer s
c
c-----Initialize the fractions array
c
        write(6,*) 'Reading the Appportionment point emissions'
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
      do s=1,Appnum
c
c-----Time
c
        read(iAppiP(s)) idat11,tim11,idat21,tim21
c
        read(iAppiP(s)) idum1,npts1
c
        read(iAppiP(s)) (idum1,idum1,idum1,flowrat1(n),
     &                   effph1(n),n=1,npts1)
c
c-----Loop over species
c       
        do k=1,nAppP(s)
c
c       Read emissions values
          read(iAppiP(s)) idum1,(name1(i),i=1,10),
     &                    (Appemis(i,k),i=1,npts1)
        enddo
        
c
c-----Calculate fractions for the sources
c
        do i=1,npts1
          do k=1,npointsp
            do ll=1,nAppP(s)
              idum=AppemisP(s,ll)
              if (lptmap(k).eq.Appmaprev(idum).and.
     &            lptmap(k).ne.0) then
c
                if (pointemis(i,k).eq.0) then
                  AppemfracP(i,ll,s+1) = 0
                else
                  AppemfracP(i,ll,s+1)=Appemis(i,ll)/
     &                                   pointemis(i,k)
                endif
c
                if((AppemfracP(i,k,s+1).lt.0.48.or.AppemfracP(i,k,s+1)
     &             .gt.52).and.AppemfracP(i,k,s+1).ne.0.0) then
                  write(6,*) 'Fraction not between 0 and 1 ....',
     &                     AppemfracP(i,k,s+1),i,k,ll,pointemis(i,k)
     &                     ,Appemis(i,ll)
                  write(6,*) 'AppemisP: ',(AppemisP(s,n),n=1,154)
                  write(6,*) 'lptmap: ',(lptmap(n),n=1,npointsp)
                  write(6,*) 'Appmaprev: ',(Appmaprev(n),n=1,154)
                  stop
                endif
c
              endif
            enddo
          enddo
        enddo
      enddo
c
c-----Assigning remaining fraction to "others" (4th dimension = 1)
c
        write(6,*) 'Location 2...'
        do i=1,npts1
          do k=1,nAppP(s)
            do s=1,Appnum
              AppemfracP(i,k,1)=AppemfracP(i,k,1)-
     &                    AppemfracP(i,k,s+1)
            enddo
          enddo
        enddo
c
      end

              
