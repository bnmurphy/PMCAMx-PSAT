      subroutine Appemiss(dconc,i,j,k,l,type,srnum,Actconc,original)
c
c     This subroutine handles the changes in apportionment due to emissions. 
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
      real dconc,Actconc,total,conv
      integer i,j,k,l,srnum,type,spc,s,nx,ny,nz
      integer ring_size(10)
c
c      if(l.eq.146.and.srnum.eq.152) then
c        write(6,*) 'Problem Spot'
c      endif
      receptor_x=65
      receptor_y=51
      ring_size(1)=1
      ring_size(2)=1
      ring_size(3)=2
      ring_size(4)=2
      ring_size(5)=2
      ring_size(6)=3
      ring_size(7)=3
      ring_size(8)=3
      ring_size(9)=4
      ring_size(10)=4
      spc = Appmap(l)
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)
      loc = i+nx*(j-1)+nx*ny*(k-1)
      conv = densfac*(273/tempk(loc))*(press(loc)/1013)
c
c     Checking the totals at the start
      total = 0.0
      do s = 1, Appnum+3
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        total = total + Appconc(loc)
      enddo
      if (abs(total-original).gt.0.01*MIN(original,total).and.
     &    Appmap(l).ne.0) then
        if (abs(original-bdnl(l)*conv).lt.0.05*bdnl(l)*conv.and.
     &      Appmap(l).lt.25.or.Appmap(l).eq.5) then
          do s = 1,Appnum+3
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &            nx*ny*nz*MXTRK*(s-1)
            Appconc(loc) = Appconc(loc)*original/total
          enddo
        elseif (abs(original-bdnl(l)).lt.0.05*bdnl(l)) then
          do s = 1,Appnum+3
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &            nx*ny*nz*MXTRK*(s-1)
            Appconc(loc) = Appconc(loc)*original/total
          enddo
        else
          write(6,*) 'Totals not same at beginning of Appemis',i,j,k,
     &               Appmap(l),l,type
          write(6,*) 'Total,old,bdnl:',total,original,bdnl(l)*conv,
     &               bdnl(l)
          write(6,*) 'New: ',Actconc
          do s=1,Appnum+3
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &            nx*ny*nz*MXTRK*(s-1)
            write(6,*) s,Appconc(loc),loc
          enddo
          stop
        endif
      endif
c
c
      if (Appmap(l).ne.0) then
        do s = 3,Appnum+3
c
c       Point Emissions
          if (type.eq.1) then
            AppemfracP(srnum,Appmap(l),1) = 1.0
            do is = 1,Appnum
              AppemfracP(srnum,Appmap(l),is+1) = 0.0
            enddo
c
          if (i.eq.receptor_x.and.j.eq.receptor_y) then
            AppemfracP(srnum,Appmap(l),2) = 1.0
            AppemfracP(srnum,Appmap(l),1) = 0.0
          endif
c
          sum_n=0
          do n = 3,11
            sum_n = sum_n + ring_size(n-2)
c            
            inner_top = receptor_y + sum_n - ring_size(n-2)
            outer_top = inner_top + ring_size(n-2)
            outer_bottom = receptor_y - sum_n
            inner_bottom = outer_bottom + ring_size(n-2)
            inner_east = receptor_x + sum_n - ring_size(n-2)
            outer_east = inner_east + ring_size(n-2)
            outer_west = receptor_x - sum_n
            inner_west = outer_west + ring_size(n-2)
c       
            if (((i.le.outer_east.and.i.ge.outer_west).AND.
     &           (j.le.outer_top.and.j.ge.outer_bottom)).AND.
     &           (.NOT.
     &          ((i.le.inner_east.and.i.ge.inner_west).AND.
     &           (j.le.inner_top.and.j.ge.inner_bottom)))) then
              AppemfracP(srnum,Appmap(l),1) = 0.0
              AppemfracP(srnum,Appmap(l),n) = 1.0
            endif  
          enddo
c
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &          nx*ny*nz*MXTRK*(s-1)
            Appconc(loc) = Appconc(loc)+AppemfracP(srnum,Appmap(l),s-2)
     &                     *dconc
          endif
c
c       Area Emissions
          if (type.eq.2) then
c
            AppemfracA(i,j,Appmap(l),1) = 1.0
            do is = 1,Appnum
              AppemfracA(i,j,Appmap(l),is+1) = 0.0
            enddo
c       
          if (i.eq.receptor_x.and.j.eq.receptor_y) then
            AppemfracA(i,j,Appmap(l),1)=0.0
            AppemfracA(i,j,Appmap(l),2) = 1.0
          endif
c
          sum_n=0
          do n = 3,11
            sum_n = sum_n + ring_size(n-2)
c       
            inner_top = receptor_y + sum_n - ring_size(n-2)
            outer_top = inner_top + ring_size(n-2)
            outer_bottom = receptor_y - sum_n
            inner_bottom = outer_bottom + ring_size(n-2)
            inner_east = receptor_x + sum_n - ring_size(n-2)
            outer_east = inner_east + ring_size(n-2)
            outer_west = receptor_x - sum_n
            inner_west = outer_west + ring_size(n-2)
c
            if (((i.le.outer_east.and.i.ge.outer_west).AND.
     &           (j.le.outer_top.and.j.ge.outer_bottom)).AND.
     &           (.NOT.
     &          ((i.le.inner_east.and.i.ge.inner_west).AND.
     &           (j.le.inner_top.and.j.ge.inner_bottom)))) then
              AppemfracA(i,j,Appmap(l),1)=0.0
              AppemfracA(i,j,Appmap(l),n) = 1.0
            endif 
          enddo
c
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &          nx*ny*nz*MXTRK*(s-1)
            Appconc(loc) = Appconc(loc) + AppemfracA(i,j,Appmap(l),s-2)
     &          *dconc
          endif
c
        enddo
c
c
c-----Check 90/10 Split
c        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
c     &        nx*ny*nz*MXTRK*(3-1)
c        loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
c     &         nx*ny*nz*MXTRK*(4-1)
c        if (abs(Appconc(loc)*9-Appconc(loc2)).gt.0.01*Appconc(loc2)
c     &      .and.Appconc(loc).gt.1E-15) 
c     &       then
c          write(6,*) 'ERROR in Appemiss: 90/10 split incorrect'
c          write(6,*) Appconc(loc), Appconc(loc2),i,j,k,Appmap(l),type
c          write(6,*) AppemfracP(srnum,spc,1),AppemfracA(i,j,spc,1)
c          write(6,*) AppemfracP(srnum,spc,2),AppemfracA(i,j,spc,2)
c          write(6,*) Actconc,l
c          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
c     &          nx*ny*nz*MXTRK*(1-1)
c          loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
c     &           nx*ny*nz*MXTRK*(2-1)
c          write(6,*) Appconc(loc), Appconc(loc2),dconc
c          stop
c        endif
c
c-----Check total
      total = 0.0
      do s=1,Appnum+3
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        total = total + Appconc(loc)
      enddo
      if (abs(total-Actconc).gt.0.01*MIN(Actconc,total)) then
        if (Actconc.eq.bdnl(l)) then
          do s=1,Appnum+3
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &            nx*ny*nz*MXTRK*(s-1)
            Appconc(loc) = Appconc(loc)*Actconc/total
          enddo
        else
          write(6,*) 'ERROR in Appemiss: total incorrect: type: ', type
          write(6,*) 'i,j,k,Appmap(l),l: ', i,j,k,Appmap(l),l
          write(6,*) 'Actual Conc.', Actconc
          do s=1,Appnum+3
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &            nx*ny*nz*MXTRK*(s-1)
            write(6,*) s,Appconc(loc),loc
          enddo
          stop
        endif
      else
        do s=1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          Appconc(loc) = Appconc(loc)*Actconc/total
        enddo
      endif
c
      endif
      if (Appconc(11010520).eq.0) then
        write(6,*) 'Equal zero',i,j,k,Appmap(l)
      endif
c

      if (i.eq.65.and.j.eq.51.and.Appmap(l).eq.75) then
        print *,'Source Attribution at Pittsburgh; Layer',k,':'
        do s = 1,Appnum+3
          loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &          nx*ny*nz*MXTRK*(s-1)
          print *,'   Source ',s,':',Appconc(loc)
        enddo
      endif

      end
