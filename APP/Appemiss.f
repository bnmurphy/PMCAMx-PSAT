      subroutine Appemiss(dconc,i,j,k,l,type,srnum,conc1,original)
c
c     This subroutine handles the changes in apportionment due to emissions. 
c  
c     DEVELOPED BY: 
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     04/13/2007
c
c     CALLED BY:
c     emiss
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
      real dconc,Actconc,total,conv,original
      integer i,j,k,l,srnum,type,spc,s,nx,ny,nz,n,m,test
      integer isize, size_bin, loc
      real conc1(ncol(1),nrow(1),nlay(1),nspec)

      integer receptor_x, receptor_y, ring_size(20), sum_n
      integer inner_top, inner_bottom, outer_top, outer_bottom
      integer inner_east, outer_east, inner_west, outer_west
c
      spc = Appmap(l)
      if (spc.eq.0) return
      
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)
      loc = i+nx*(j-1)+nx*ny*(k-1)
      conv = densfac*(273/tempk(loc))*(press(loc)/1013)
      if (Appmap(l).gt.sa_num_gas) conv = 1.0
      Actconc = conc1(i,j,k,l)
c
c     Checking the totals at the start
      total = 0.0
      !print *,'Appemiss: i=',i,' j=',j,' k=',k,' conv=',conv
      !print *,'          l=',l
      do s = 1, Appnum+3
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        !print *,'    Appconc=',Appconc(loc),'  total=',total
        if (Appconc(loc).lt.bdnl(l)*conv) Appconc(loc) = bdnl(l)*conv
        !if (Appconc(loc).eq.'NaN') Appconc(loc) = bdnl(l)*conv
        total = total + Appconc(loc)
        !print *,'    Appconc=',Appconc(loc),'  total=',total
      enddo
      if (total.eq.0.0) then  !The species ran out of mass.
                              !Put the lower-bound mass in the s1 bin
        !total = bdnl(l)
        total = bdnl(l)*conv
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &        nx*ny*nz*MXTRK*(1-1)
        !Appconc(loc) = bdnl(l)
        Appconc(loc) = bdnl(l)*conv
      endif
      
      if (Appmap(l).ge.sa_num_sv) then
        !The Species is Semivolatile Aerosol. Sum up the 
        !size bins before making a comparison        

        !Find out what bin this is
        size_bin = l - Appmaprev(Appmap(l)) + 1
        !Sum up before and after
        original = 0.0
        Actconc = 0.0
        !print *,'Appemiss: l=',l,'  Appmap=',Appmap(l)
        do isize = 1,6
          original = original + conc1(i,j,k,l-size_bin+isize)
          Actconc  = Actconc  + conc1(i,j,k,l-size_bin+isize)
          !print *,'    isize=',isize,'  conc1=',conc1(i,j,k,l-size_bin+isize),' bdnl=',bdnl(l-size_bin+isize)
        enddo
        original = original - dconc
        !print *,'   original = ',original,' Actconc=',Actconc,' dconc=',dconc
      else
        original = original
        Actconc = Actconc
        !print *,'Appemiss: l=',l,'  Appmap=',Appmap(l)
        !print *,'  conc1=',conc1(i,j,k,l),' bdnl=',bdnl(l)
      endif
      
      
      if (abs(total-original).gt.0.01*MIN(original,total).and.
     &    Appmap(l).ne.0.and.total.gt.0.0) then
        if (abs(original-bdnl(l)*conv).lt.0.05*bdnl(l)*conv.and.
     &      Appmap(l).lt.25.or.Appmap(l).eq.5) then
          do s = 1,Appnum+3
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &            nx*ny*nz*MXTRK*(s-1)
            Appconc(loc) = Appconc(loc)*original/total
          enddo
        elseif (abs(original-bdnl(l)).lt.0.05*bdnl(l).or.
     &         abs(total-original).lt.0.05*MIN(original,total).or.
     &         abs(total-original).lt.(Appnum+4)*bdnl(l)*conv
     &          ) then
          do s = 1,Appnum+3
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &            nx*ny*nz*MXTRK*(s-1)
            Appconc(loc) = Appconc(loc)*original/total
          enddo
        elseif (original.ne.0) then
          write(6,*) 'Totals not same at beginning of Appemis',i,j,k,
     &                 Appmap(l),l,type
          write(6,*) 'Total,old,bdnl:',total,original,bdnl(l)*conv,
     &               bdnl(l)
          write(6,*) '  Source Apportionment:'
          do s = 1,Appnum+3
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &            nx*ny*nz*MXTRK*(s-1)
            write(6,*)'       s=',s,'  Appconc=',Appconc(loc)
          enddo
          write(6,*) 'New: ',Actconc
          stop
        endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cTHIS CODE IS FOR PSAT APPLIED TO AGE
      !8121 is for May, 2008
      !1192 is for July, 2001
      !do n=1,Appnum+1
      !  AppemfracP(srnum,Appmap(l),n) = 0.0
      !  if (n.eq.date-8121) AppemfracP(srnum,Appmap(l),n) = 1.0
      !  AppemfracA(i,j,Appmap(l),n) = 0.0
      !  if (n.eq.date-8121) AppemfracA(i,j,Appmap(l),n) = 1.0
      !enddo
cAGE CODE!!!!!! BNM

cTHIS CODE IS FOR PSAT APPLIED TO SOURCE LOCATION
c     
c     !Set Receptor Location
c     receptor_x = 65      !Pittsburgh (65,51)
c     receptor_y = 51
c     !receptor_x = 79      !New York (79,55)
c     !receptor_y = 55
c     !receptor_x = 70      !Duke Forest (70,38)
c     !receptor_y = 38
c     !receptor_x = 51      !Bowling Green, KY
c     !receptor_y = 38
c     ring_size(1) = 1
c     ring_size(2) = 1
c     ring_size(3) = 2
c     ring_size(4) = 2
c     ring_size(5) = 2
c     ring_size(6) = 3
c     ring_size(7) = 3
c     ring_size(8) = 3
c     ring_size(9) = 4
c     ring_size(10) = 4

c       AppemfracP(srnum,Appmap(l),Appnum+1) = 1.0
c       AppemfracA(i,j,Appmap(l),Appnum+1) = 1.0
c       do n=1,Appnum
c         AppemfracP(srnum,Appmap(l),n) = 0.0
c         AppemfracA(i,j,Appmap(l),n) = 0.0
c       enddo

c       !Loop Through Source Categories
c       if (i.eq.receptor_x.and.j.eq.receptor_y) then
c         !Local Emissions
c         AppemfracP(srnum,Appmap(l),1) = 1.0
c         AppemfracA(i,j,Appmap(l),1) = 1.0
c         !Erase Attribution from OTHER Category
c         AppemfracP(srnum,Appmap(l),Appnum+1) = 0.0
c         AppemfracA(i,j,Appmap(l),Appnum+1) = 0.0
c       endif

c       sum_n = 0
c       do n = 2,Appnum
c         sum_n = sum_n + ring_size(n-1)
c         outer_top = receptor_y + sum_n
c         inner_top = outer_top - ring_size(n-1)
c         outer_bottom = receptor_y - sum_n
c         inner_bottom = outer_bottom + ring_size(n-1)
c         outer_east = receptor_x + sum_n
c         inner_east = outer_east - ring_size(n-1)
c         outer_west = receptor_x - sum_n
c         inner_west = outer_west + ring_size(n-1)
c           
c         if ((i.le.outer_east.and.i.ge.outer_west.AND.
c    &       j.le.outer_top.and.j.ge.outer_bottom).AND.
c    &      .NOT.
c    &      (i.le.inner_east.and.i.ge.inner_west.AND.
c    &       j.le.inner_top.and.j.ge.inner_bottom)) then
c           !Source is inside this ring
c           AppemfracP(srnum,Appmap(l),n) = 1.0
c           AppemfracA(i,j,Appmap(l),n) = 1.0
c           !Erase Attribution from OTHER Category
c           AppemfracP(srnum,Appmap(l),Appnum+1) = 0.0
c           AppemfracA(i,j,Appmap(l),Appnum+1) = 0.0
c         endif
c       enddo
c   
      !This is for Paris
      if (i.eq.58.and.j.eq.73) then
        n=1   !Local Emissions
        AppemfracP(srnum,Appmap(l),n) = 1.0
        AppemfracA(i,j,Appmap(l),n) = 1.0
      else
        n=2   !Everything Else
        AppemfracP(srnum,Appmap(l),n) = 1.0
        AppemfracA(i,j,Appmap(l),n) = 1.0
      endif
cLOCATION CODE        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c
      if (Appmap(l).ne.0) then
        do s = 3,Appnum+3
c
c       Point Emissions
          if (type.eq.1) then    
c
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &          nx*ny*nz*MXTRK*(s-1)
            Appconc(loc) = Appconc(loc)+AppemfracP(srnum,Appmap(l),s-2)
     &                     *dconc
          endif
c
c       Area Emissions
          if (type.eq.2) then
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &          nx*ny*nz*MXTRK*(s-1)
c            if (k.eq.1.and.i.eq.2.and.j.eq.18) then
c              print *,'Appemiss: Assigning Area. spc=',Appmap(l),' s=',s,
c     &             ' Appconc=',Appconc(loc),'  AppemfracA=',AppemfracA(i,j,Appmap(l),s-2)
c	    endif
            Appconc(loc) = Appconc(loc) + AppemfracA(i,j,Appmap(l),s-2)
     &          *dconc
c            if (k.eq.1.and.i.eq.2.and.j.eq.18) then
c              print *,'Appemiss: Assigning Area. Appconc=',Appconc(loc)
c	    endif
          endif
c
        enddo
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
          write(6,*) 'i,j,k,Appmap(l): ', i,j,k,Appmap(l)
          write(6,*) 'Actual Conc.', Actconc,total,dconc
          do s=1,Appnum+3
            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &            nx*ny*nz*MXTRK*(s-1)
            write(6,*) s,Appconc(loc)
          enddo
          if (type.eq.1) then
            do s=1,Appnum+1
               loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              write(6,*) s,AppemfracP(srnum,Appmap(l),s)
            enddo
          endif        
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
      end
