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
c
      spc = Appmap(l)
      if (spc.eq.0) return
      
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)
      loc = i+nx*(j-1)+nx*ny*(k-1)
      conv = densfac*(273/tempk(loc))*(press(loc)/1013)

      !if (i.eq.2.and.j.eq.18.and.Appmap(l).eq.130.and.k.eq.1) then
      ! loc = 2+nx*(18-1) + nx*ny*(1-1)+nx*ny*nz*(130-1)+nx*ny*nz*MXTRK*(1-1)
      ! print *,'\nBegin Appemiss: Appconc(i=2,j=17,spc=130,s=1)=',Appconc(loc)
      ! loc = 2+nx*(18-1) + nx*ny*(1-1)+nx*ny*nz*(130-1)+nx*ny*nz*MXTRK*(2-1)
      ! print *,'Begin Appemiss: Appconc(i=2,j=17,spc=130,s=2)=',Appconc(loc)
      ! loc = 2+nx*(18-1) + nx*ny*(1-1)+nx*ny*nz*(130-1)+nx*ny*nz*MXTRK*(3-1)
      ! print *,'Begin Appemiss: Appconc(i=2,j=17,spc=130,s=3)=',Appconc(loc)
      ! loc = 2+nx*(18-1) + nx*ny*(1-1)+nx*ny*nz*(130-1)+nx*ny*nz*MXTRK*(4-1)
      ! print *,'Begin Appemiss: Appconc(i=2,j=18,spc=130,s=4)=',Appconc(loc)
      ! print *,'dconc=',dconc,'  actconc=',actconc,'  original=',original 
      !endif
c
c     Checking the totals at the start
      total = 0.0
      do s = 1, Appnum+3
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &        nx*ny*nz*MXTRK*(s-1)
        total = total + Appconc(loc)
      enddo
      if (total.eq.0.0) then  !The species ran out of mass.
                              !Put the lower-bound mass in the s1 bin
        total = bdnl(l)
        loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(Appmap(l)-1) +
     &        nx*ny*nz*MXTRK*(1-1)
        Appconc(loc) = bdnl(l)
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
     &         original.lt.1.0E-7) then
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
          write(6,*) 'New: ',Actconc
          stop
        endif
      endif
c
cTHIS IS FOR THE AGE CODE!!!!!!
      do n=1,Appnum+1
        AppemfracP(srnum,Appmap(l),n) = 0.0
        if (n.eq.date-1192) AppemfracP(srnum,Appmap(l),n) = 1.0
        AppemfracA(i,j,Appmap(l),n) = 0.0
        if (n.eq.date-1192) AppemfracA(i,j,Appmap(l),n) = 1.0
      enddo
cAGE CODE!!!!!! BNM
      
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
      if (k.eq.1.and.i.eq.2.and.j.eq.18) print *,'Appemiss: spc=',Appmap(l),' total=',total
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
      if (i.eq.4.and.j.eq.3.and.k.eq.14) then
       print *,'spc =',l,'  Appspc=',Appmap(l) 
       print *,'i=',i,'  j=',j
       loc = 4+nx*(3-1) + nx*ny*(1-1)+nx*ny*nz*(2-1)+nx*ny*nz*MXTRK*(1-1)
       print *,'Appemiss: Appconc(i=4,j=3,spc=2,s=1)=',Appconc(loc)
       loc = 4+nx*(3-1) + nx*ny*(1-1)+nx*ny*nz*(2-1)+nx*ny*nz*MXTRK*(2-1)
       print *,'Appemiss: Appconc(i=4,j=3,spc=2,s=2)=',Appconc(loc)
       loc = 4+nx*(3-1) + nx*ny*(1-1)+nx*ny*nz*(2-1)+nx*ny*nz*MXTRK*(3-1)
       print *,'Appemiss: Appconc(i=4,j=3,spc=2,s=3)=',Appconc(loc)
       loc = 4+nx*(3-1) + nx*ny*(1-1)+nx*ny*nz*(2-1)+nx*ny*nz*MXTRK*(4-1)
       print *,'Appemiss: Appconc(i=4,j=3,spc=2,s=4)=',Appconc(loc)
      endif
c
      end
