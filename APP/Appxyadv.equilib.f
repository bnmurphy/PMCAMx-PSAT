      subroutine Appxyadv(fp,fm,spc,xy,k,ij,dep,scale,dx,l)
c
c     This subroutine handles the changes in apportionment due to
c     horizontal advection. 
c  
c     DEVELOPED BY: 
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     04/13/2007
c
c     CALLED BY:
c     xyadvec
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
      real fp(MX1D),fm(MX1D)  !,Actconc(MXCOL1,MXROW1,MXLAY1,MXSPEC)
      !real old(MXCOL1,MXROW1,MXLAY1,MXSPEC)
      !BNM Replaced AppAFT and old with AppAFT and AppFORE, respectively
      integer spc,xy,k,ij,loc,loc2,loc3,s,src,nx,ny,nz,maxi,maxj
      integer mini,minj,v,l
      real Apptot(MXCOL1),Original(MXCOL1,MXSOUR),total
      real dep(MX1D),scale(MX1D),dx(MX1D),convfac
      real original2, original3(MX1D), size_bin
c
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)     
c 
      if (xy.eq.1) then !X-Advection
c
c     Totals
        
        do i=1,MXCOL1
          Apptot(i) = 0
          do src=1,Appnum+3
            loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &            nx*ny*nz*MXTRK*(src-1)
            Apptot(i) = Appconc(loc)+Apptot(i)
            Original(i,src)=Appconc(loc)
            print *,'i=',i,' src=',src,'  Apptot(i)=',Apptot(i)
            print *,'   Appconc=',Appconc(loc),'  Original=',Original(i,src)
          enddo
          loc2 = i + nx*(ij-1) + nx*ny*(k-1) + 
     &           nx*ny*nz*(Appmaprev(spc)-1)
          loc3 = i + nx*(ij-1) + nx*ny*(k-1)
          if (spc.lt.sa_num_gas) then
            convfac = densfac*(273/tempk(loc3))*(press(loc3)/1013)
          else
            convfac = 1
          endif

          if (spc.ge.sa_num_sv) then
            !Species is a semivolatile aerosol. Sum up size bins
            !Find out size bin
            size_bin = l - Appmaprev(spc) + 1
            original2    = 0.0
            original3(i) = 0.0
            print *,'i=',i,' j=',ij,' k=',k,'  spc=',l
            do isize = 1,10
              original2    = original2 + AppFORE(i,ij,k,Appmaprev(spc)-1+isize)
              original3(i) = original3(i) + AppAFT(i,ij,k,Appmaprev(spc)-1+isize)
              print *,'  isize=',isize,' original2=',original2,' AppFORE=',AppFORE(i,ij,k,Appmaprev(spc)-1+isize)
              print *,'  isize=',isize,' original3=',original3(i),' AppAFT=',AppAFT(i,ij,k,Appmaprev(spc)-1+isize)
            enddo
            print *,'   Apptot(i)=',Apptot(i)
          else
            original2    = AppFORE(i,ij,k,Appmaprev(spc))
            original3(i) = AppAFT(i,ij,k,Appmaprev(spc))
          endif

          if (abs(Apptot(i)-original2).gt.
     &        0.01*MIN(Apptot(i),original2).or.
     &        original2.eq.bdnl(Appmaprev(spc))*convfac) then
            if (abs(original2-bdnl(Appmaprev(spc))
     &          *convfac).lt.0.05*bdnl(Appmaprev(spc))*convfac.or.spc.eq.5) then
              do s = 1,Appnum+3
                loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                if (Apptot(i).eq.0) then
                  if (s.eq.1) Appconc(loc) = original2
                  if (s.ne.1) Appconc(loc) = 0.0
                else
                  Appconc(loc) = Appconc(loc)*original2
     &                         /Apptot(i)
                endif
                Original(i,s)=Appconc(loc)
              enddo
              Apptot(i) = original2
            elseif (i.eq.1.or.i.eq.97.or.j.eq.1.or.j.eq.90) then  !BNMchanged 96->97
              do s = 1,Appnum+3
                loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                if (Apptot(i).eq.0) then
                  if (s.eq.1) Appconc(loc) = original2
                  if (s.ne.1) Appconc(loc) = 0.0
                else
                  Appconc(loc) = Appconc(loc)*original2
     &                         /Apptot(i)
                endif
                Original(i,s)=Appconc(loc)
                print *,'  Corrected: Appconc=',Appconc(loc),'  Original=',Original(i,s)
                print *,'             Apptot=',Apptot(i)
              enddo
            else
              write(6,*) 'Totals not same at beginning of Appxyadv',i,ij
     &                    ,k,spc
              write(6,*) 'Total,old,bdnl:',Apptot(i),
     &                   original2,bdnl(Appmaprev(spc))
              write(6,*) 'New: ',AppAFT(i,ij,k,Appmaprev(spc))
              stop
            endif
          endif
          print *,'  After total check: Apptot(i)=',Apptot(i)
        enddo
c
c
        do s=1,Appnum+3
c
          !do i=2,MXCOL1-1
          do i=1,MXCOL1
c
c           BOUNDARY CELLS DO NOT CHANGE
            !if (i.eq.1.and.fm(2).eq.0.and.fp(1).eq.0) then
            if (i.eq.1) then
              loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*mxtRK*(s-1)
              loc3 = i + nx*(ij-1) + nx*ny*(k-1) + 
     &               nx*ny*nz*(Appmaprev(spc)-1)
              !Appconc(loc)=conc(loc3)*(Original(i,s)/Apptot(i))
              !BNM editing to make semivolatiles independent of size
              print *,' Boundary Cell: i=',i,' j=',j,' s=',s,' original3=',original3(i)
              print *,'     Original(i,s)=',Original(i,s),' Apptot=',Apptot(i),' Appconc=',Appconc(loc)
              Appconc(loc)=original3(i)*(Original(i,s)/Apptot(i)) 
              print *,'     Appconc=',Appconc(loc)
c
            !elseif (i.eq.1.and.(fm(2).ne.0 .or. fp(1).ne.0) ) then
            !  loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &      !        nx*ny*nz*MXTRK*(s-1)
              !print *,'\n  Appxyadvec: fp(1) = ',fp(1),' fm(2)=',fm(2)
              !print *,'   i=',i,' j=',ij,' k=',k,' spc=',spc,' s=',s
              !print *,'   loc=',loc,' Appconc=',Appconc(loc)
              !print *,'   fp(i)=',fp(i),' fm(i)=',fm(i),' fm(i+1)=',fm(i+1)
              !print *,'   Apptot(i)=',Apptot(i),'  Apptot(i+1)=',Apptot(i+1)
              !print *,'   Original(i,s)=',Original(i,s)
            !  Appconc(loc)=Appconc(loc) - ( fp(i)*Original(i,s)/
     &      !               Apptot(i) + fm(i+1)*Original(i+1,s)/
     &      !               Apptot(i+1) )*scale(i)/dep(i)
c
            elseif (i.eq.MXCOL1) then
              loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              loc3 = i + nx*(ij-1) + nx*ny*(k-1) + 
     &               nx*ny*nz*(Appmaprev(spc)-1)
              !Appconc(loc)=conc(loc3)*(Original(i,s)/Apptot(i))
              !BNM editing to make semivolatiles independent of size
              Appconc(loc)=original3(i)*(Original(i,s)/Apptot(i))
c
            !elseif (i.eq.MX1D.and.fm(MX1D-1).ne.0) then
            !elseif (i.eq.MXCOL1.and.(fp(MXCOL1-1).ne.0 .or. fm(MXCOL1).ne.0) ) then
            !  loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &      !        nx*ny*nz*MXTRK*(s-1)
            !  Appconc(loc)=Appconc(loc) + (fp(i-1)*Original(i-1,s)/
     &      !               Apptot(i-1) - (fp(i)+fm(i))*Original(i,s)/
     &      !               Apptot(i))*scale(i)/dep(i)
c
c           NON-BOUNDARY CELLS
            else
              loc = i + nx*(ij-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              Appconc(loc)=Appconc(loc) + (fp(i-1)*Original(i-1,s)/
     &                     Apptot(i-1) - (fp(i)+fm(i))*Original(i,s)/
     &                     Apptot(i) + fm(i+1)*Original(i+1,s)/
     &                     Apptot(i+1))*scale(i)/dep(i)
              if (Appconc(loc).lt.0.and.
     &            abs(Appconc(loc)).lt.bdnl(Appmaprev(spc))*convfac) then
                 Appconc(loc) = 0.0
              endif
              if (Appconc(loc).eq.'NaN') 
     &          write(6,*) 'Problem in Appxyadv - x'
            endif
          enddo
        enddo
c
      elseif (xy.eq.2) then !Y-Advection
c
      !if (k.eq.1.and.ij.eq.2) then 
      !  loc = 2+nx*(18-1)+nx*ny*(1-1)+nx*ny*nz*(130-1)+nx*ny*nz*MXTRK*(4-1)
      !  print *,'Appyadv: spc=',spc,' Appconc(i=2,j=18,spc=130,s=4)=',Appconc(loc)
      !endif
c       Totals

        do j=1,MXROW1
          Apptot(j) = 0
          do src=1,Appnum+3
            loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &            nx*ny*nz*MXTRK*(src-1)
            Apptot(j) = Appconc(loc)+Apptot(j)
            Original(j,src) = Appconc(loc)
          enddo
          loc2 = ij + nx*(j-1) + nx*ny*(k-1) + 
     &           nx*ny*nz*(Appmaprev(spc)-1)
          loc3 = ij + nx*(j-1) + nx*ny*(k-1)
          if (spc.lt.sa_num_gas) then
            convfac = densfac*(273/tempk(loc3))*(press(loc3)/1013)
          else
            convfac = 1
          endif

          if (spc.ge.sa_num_sv) then
            !Species is a semivolatile aerosol. Sum up size bins
            !Find out size bin
            size_bin = l - Appmaprev(spc) + 1
            original2    = 0.0
            original3(j) = 0.0
            do isize = 1,10
              original2    = original2 + AppFORE(ij,j,k,Appmaprev(spc)-1+isize)
              original3(j) = original3(j) + AppAFT(ij,j,k,Appmaprev(spc)-1+isize)
            enddo
          else
            original2    = AppFORE(ij,j,k,Appmaprev(spc))
            original3(j) = AppAFT(ij,j,k,Appmaprev(spc))
          endif
          
          if (abs(Apptot(j)-original2).gt.
     &        0.02*MIN(Apptot(j),original2)) then
            if (abs(original2-bdnl(Appmaprev(spc))*convfac)
     &          .lt.0.01*original2.or.
     &          spc.eq.5) then
              do s = 1,Appnum+3
                loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                if (Apptot(j).eq.0) then
                  if (s.eq.1) Appconc(loc) = original2
                  if (s.ne.1) Appconc(loc) = 0.0
                else
                  Appconc(loc) = Appconc(loc)*original2/Apptot(j)
                endif
                Original(j,s)=Appconc(loc)
              enddo
              Apptot(j) = original2
            else
              write(6,*) 'Totals not same at beginning of Appxyadv',ij,j
     &                    ,k,spc,xy,convfac
              write(6,*) 'Total,old,bdnl:',Apptot(j),
     &                   original2,bdnl(Appmaprev(spc))
              write(6,*) 'New: ',AppAFT(ij,j,k,Appmaprev(spc))
              do s = 1,Appnum+3
                loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                write(6,*) s,Appconc(loc)
              enddo
              stop
            endif
          else
            do s = 1,Appnum+3
              loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              Appconc(loc) = Appconc(loc)*original2
     &                       /Apptot(j)
              Original(j,s)=Appconc(loc)
            enddo
            Apptot(j) = original2
          endif

        enddo
c
        do s=1,Appnum+3
c
          !do j=2,MXROW1-1
          do j=1,MXROW1
c
c           BOUNDARY CELLS DO NOT CHANGE
            !if (j.eq.1.and.fp(2).eq.0) then
            if (j.eq.1) then
              loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              loc3 = ij + nx*(j-1) + nx*ny*(k-1) + 
     &               nx*ny*nz*(Appmaprev(spc)-1)
              !Appconc(loc)=conc(loc3)*(Original(j,s)/Apptot(j))
              !BNM editing to make semivolatiles independent of size
              Appconc(loc)=original3(j)*(Original(j,s)/Apptot(j))
c
            !elseif (j.eq.1.and.fp(2).ne.0) then
            !  loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &      !        nx*ny*nz*MXTRK*(s-1)
            !  Appconc(loc)=Appconc(loc) - (fp(j)+fm(j))*Original(j,s)/
     &      !               Apptot(j) + fm(j+1)*Original(j+1,s)/
     &      !               Apptot(j+1)
c
            !elseif (j.eq.MX1D.and.fm(MX1D-1).eq.0) then
            !elseif (j.eq.MXROW1.and.fm(MXROW1-1).eq.0) then
            elseif (j.eq.MXROW1) then
              loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              loc3 = ij + nx*(j-1) + nx*ny*(k-1) + 
     &               nx*ny*nz*(Appmaprev(spc)-1)
              !Appconc(loc)=conc(loc3)*(Original(j,s)/Apptot(j))
              !BNM editing to make semivolatiles independent of size
              Appconc(loc)=original3(j)*(Original(j,s)/Apptot(j))
c
            !elseif (j.eq.MX1D.and.fm(MX1D-1).ne.0) then
            !elseif (j.eq.MXROW1.and.fm(MXROW1-1).ne.0) then
            !  loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &      !        nx*ny*nz*MXTRK*(s-1)
            !  Appconc(loc)=Appconc(loc) + fp(j-1)*Original(j-1,s)/
     &      !               Apptot(j-1) - (fp(j)+fm(j))*Original(j,s)/
     &      !               Apptot(j)
c
c           NON-BOUNDARY CELLS
            else
              loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &              nx*ny*nz*MXTRK*(s-1)
              Appconc(loc)=Appconc(loc) + (fp(j-1)*Original(j-1,s)/
     &                     Apptot(j-1) - (fp(j)+fm(j))*Original(j,s)/
     &                     Apptot(j) + fm(j+1)*Original(j+1,s)/
     &                     Apptot(j+1))*scale(j)/dep(j)/dx(j)
              if (Appconc(loc).lt.0.and.
     &           abs(Appconc(loc)).lt.bdnl(Appmaprev(spc))*convfac) then
                 Appconc(loc) = 0.0
              endif
            endif
c            if (Appconc(loc).eq.'NaN') then
c              write(6,*) 'NaN Problem in Appxyadv - y:',ij,j,k,spc,s
c     &                    ,Appconc(loc)
c              do s=1,4
c                loc = ij + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &                nx*ny*nz*MXTRK*(s-1)
c                write(6,*) Appconc(loc),Original(j,s)
c              enddo
c              stop
c            endif
          enddo
        enddo
      endif
c
c----Checks
      if (xy.eq.1) then
        maxi = MXCOL1
        mini = 1
        maxj = ij
        minj = ij
      else
        maxi = ij
        mini = ij
        maxj = MXROW1
        minj = 1
      endif
c
c----Checks for 90/10 split
c      do i=mini,maxi
c        do j=minj,maxj
c            loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &            nx*ny*nz*MXTRK*(3-1)
c            loc2 = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
c     &             nx*ny*nz*MXTRK*(4-1)
c            if (abs(Appconc(loc)*9-Appconc(loc2)).gt.0.1*Appconc(loc2)
c     &          .and.Appconc(loc).gt.1E-15) 
c     &          then
c              write(6,*) 'ERROR in Appxyadv: 90/10 split incorrect'
c              write(6,*) 'Other, Source 1: ',Appconc(loc), Appconc(loc2)
c              write(6,*) 'i,j,k,spc,xy?,ij: ',i,j,k,spc,xy,ij
c              loc = i + nx*(j-1) + nx*ny*(k-1)+nx*ny*nz*(Appmap(l)-1) +
c     &              nx*ny*nz*MXTRK*(1-1)
c              loc2 = i + nx*(j-1) + nx*ny*(k-1)+nx*ny*nz*(Appmap(l)-1) +
c     &               nx*ny*nz*MXTRK*(2-1)
c              write(6,*) 'IC,BC :',Appconc(loc), Appconc(loc2)
c              stop
c            endif
c        enddo
c      enddo
c
c
c-----Check total
      do i=mini,maxi
        do j=minj,maxj
            total = 0.0
            do s=1,Appnum+3
              loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &            nx*ny*nz*MXTRK*(s-1)
              total = total + Appconc(loc)
              if (Appconc(loc).lt.0) then
                write(6,*) 'Negative Value in Appxyadv'
                write(6,*) 'i,j,k,spc,s,xy?: ',i,j,k,spc,s,xy
                write(6,*) 'loc:          ', loc
                write(6,*) 'Appconc(loc): ', Appconc(loc)
                write(6,*) 'Actual Conc.', AppAFT(i,j,k,Appmaprev(spc))
                write(6,*) 'Apptot ', Apptot
                write(6,*) 'total ', total
                do v=1,Appnum+3
                  loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                  nx*ny*nz*MXTRK*(v-1)
                  write(6,*) v,Appconc(loc),Original(j,v),Original(j+1,v)
                enddo
                write(6,*) 'fp(i-1),fp(i),fm(i),fm(i+1):',
     &                      fp(j-1),fp(j),fm(j),fm(j+1)
                write(6,*) 'scale,depth,dx: ',
     &                      scale(j),dep(j),dx(j)
                stop
              endif
            enddo

            if (spc.ge.sa_num_sv) then
              print *,'Appxyadv: i=',i,' j=',j,' k=',k,' spc=',spc
              print *,'          total=',total,'  AppAFT=',AppAFT(i,j,k,Appmaprev(spc))
              !Species is a semivolatile aerosol. Sum up size bins
              !Find out size bin
              size_bin = l - Appmaprev(spc) + 1
              original2 = 0.0
              do isize = 1,10
                original2 = original2 + AppAFT(i,j,k,Appmaprev(spc)-1+isize)
              enddo
            else
              original2 = AppAFT(i,j,k,Appmaprev(spc))
            endif
          
            if (abs(total-original2)
     &          .gt.0.05*MIN(original2,total) ) then
              if (spc.ge.sa_num_sv
     &            .and.original2.le.10.5*bdnl(Appmaprev(spc))*convfac ) then
                  goto 300
              endif
              write(6,*) 'ERROR in Appxyadv: total incorrect'
              write(6,*) 'i,j,k,spc,xy?,ij: ',i,j,k,spc,xy,ij
              write(6,*) 'Actual Conc.', original2
              write(6,*) 'total ', total, '  bdnl ',bdnl(Appmaprev(spc))*convfac
              if (xy.eq.1) then
               do s=1,Appnum+3
                 loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                 write(6,*) s,Appconc(loc),Original(i,s),Original(i+1,s)
               enddo
               write(6,*) 'fp(i-1),fp(i),fm(i),fm(i+1):',
     &                     fp(i-1),fp(i),fm(i),fm(i+1)
               write(6,*) 'scale,depth,dy: ',
     &                     scale(i),dep(i),dx(i)
              else
               do s=1,Appnum+3
                 loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                 write(6,*) s,Appconc(loc),Original(j,s),Original(j+1,s)
               enddo
               write(6,*) 'fp(j-1),fp(j),fm(j),fm(j+1):',
     &                     fp(j-1),fp(j),fm(j),fm(j+1)
               write(6,*) 'scale,depth,dx: ',
     &                     scale(j),dep(j),dx(j)
              endif
              stop
            else
 300          do s=1,Appnum+3
                loc = i + nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(spc-1) +
     &                nx*ny*nz*MXTRK*(s-1)
                Appconc(loc) = Appconc(loc)*
     &            original2/total
              enddo
            endif
        enddo
      enddo
c
      end
