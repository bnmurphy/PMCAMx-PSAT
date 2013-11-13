      subroutine Appalign(actual)
c
c     This program realigns the mass totals between CAMx and PSAT after
c     the source tracers are sent through the vertical diffusion
c     CAMx routines.
c  
c     DEVELOPED BY: 
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     04/13/2007
c
c     CALLED BY:
c     emistrns
c
c     ROUTINES CALLED:
c     none
c
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
      !real saconc(MXCOL1,MXROW1,MXLAY1,MXTRK*MXSOUR)
      real actual(MXCOL1,MXROW1,MXLAY1,MXSPEC)
      real totals, app_align_con
      integer i,j,k,n,s,isize
c
c     COPY concentrations from saconc to Appconc
c
      do i=1,MXCOL1
        do j=1,MXROW1
          do k=1,MXLAY1
            do n=1,MXTRK
              totals = 0.0
              do s=1,MXSOUR
                loc1=n+MXTRK*(s-1)
                loc2=i + MXCOL1*(j-1) + MXCOL1*MXROW1*(k-1) +
     &               MXCOL1*MXROW1*MXLAY1*(n-1) +
     &               MXCOL1*MXROW1*MXLAY1*MXTRK*(s-1)
                Appconc(loc2) = saconc(i,j,k,loc1)
                totals = totals + Appconc(loc2)
              enddo

              if (n.ge.sa_num_sv) then
                app_align_con = 0.0
                do isize = Appmaprev(n),Appmaprev(n)+sv_bin
                  app_align_con = app_align_con + actual(i,j,k,isize)
                enddo
              else
                app_align_con = actual(i,j,k,Appmaprev(n))
              endif             
                
              !if (totals.ne.actual(i,j,k,Appmaprev(n))) then
              if (totals.ne.app_align_con) then
                do s=1,MXSOUR
                  loc2=i + MXCOL1*(j-1) + MXCOL1*MXROW1*(k-1) +
     &                 MXCOL1*MXROW1*MXLAY1*(n-1) +
     &                 MXCOL1*MXROW1*MXLAY1*MXTRK*(s-1)
                  Appconc(loc2) = Appconc(loc2) *
     &                            app_align_con/totals
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
      return
      end

