       subroutine Appreadar(areaemis,nareasp)
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
      integer s, nareasp
      logical match
c
c-----Initialize the fractions array
c
        write(6,*) 'Reading Apportionment Area emissions'
        write(6,*) 'Date: ',date
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
      do s=1,Appnum
c
c-----Time
c
        read(iAppiA(s)) idat11,tim11,idat21,tim21
c
c-----Loop over species
c       
        do k=1,nAppA(s)
c
c       Read emissions values
          read(iAppiA(s)) idum1,(name1(i),i=1,10),
     &                    ((Appemis(i,j,k),i=1,MXCOL1),j=1,MXROW1)
        enddo
        
c
c-----Calculate fractions for the sources
c
        do i=1,MXCOL1
          do j=1,MXROW1
            do k=1,nareasp
              do ll=1,nAppA(s)
                idum=AppemisA(s,ll)
                if (larmap(k,1).eq.Appmaprev(idum).and.
     &              larmap(k,1).ne.0) then
c
                  if (areaemis(i,j,k).eq.0) then
                    AppemfracA(i,j,ll,s+1) = 0
                  else
                    AppemfracA(i,j,ll,s+1)=Appemis(i,j,ll)/
     &                                     areaemis(i,j,k)
                  endif
c
                  if(AppemfracA(i,j,k,s+1).lt.0.or.AppemfracA(i,j,k,s+1)
     &               .gt.1) then
                    write(6,*) 'Fraction not between 0 and 1 ....',
     &                       AppemfracA(i,j,k,1),i,j,k,areaemis(i,j,k)
     &                       ,Appemis(i,j,ll)
                    stop
                  endif

                endif
              enddo
            enddo
          enddo
        enddo
      enddo
c
c-----Assigning remaining fraction to "others" (4th dimension = 1)
c
        do i=1,MXCOL1
          do j=1,MXROW1
            do k=1,nAppA(s)
              do s=1,Appnum
                AppemfracA(i,j,k,1)=AppemfracA(i,j,k,1)-
     &                    AppemfracA(i,j,k,s+1)
              enddo

            enddo
          enddo
        enddo
c
      end

              
