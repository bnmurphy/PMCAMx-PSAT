      subroutine Appprep()
c
c     This subroutine maps the Apportionment species to the modeled species 
c     and prints out the headers for the output file.
c  
c     DEVELOPED BY: 
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     08/19/2006
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
      integer i,j,ISUND, PM_indx
      character*2 num(10)
      logical match
c
      data num /'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ','10'/
c
      count = 0
c
c-----Assigning species to the apportionment name array
c
c      Appnam(1)  = 'NO        '
c      Appnam(2)  = 'NO2       '
c      Appnam(3)  = 'PAN       '
c      Appnam(4)  = 'NXOY      '
      Appnam(1)  = 'OLE1      '
      Appnam(2)  = 'OLE2      '
      Appnam(3)  = 'TERP      '
      Appnam(4) = 'ALK4      '
      Appnam(5)  = 'ALK5      '
      Appnam(6)  = 'ARO1      '
c      Appnam(10) = 'HONO      '
c      Appnam(11) = 'HNO3      '
      Appnam(7) = 'ISOP      '
c      Appnam(13) = 'PBZN      '
c      Appnam(14) = 'PAN2      '
      Appnam(8) = 'SO2       '
      Appnam(9) = 'SULF      '
      Appnam(10) = 'NH3       '
c      Appnam(18) = 'MPAN      '
      Appnam(11) = 'HCL       '
c      Appnam(20) = 'HNO4      '
      Appnam(12) = 'ARO2      '
      Appnam(13) = 'CPO1      '
      Appnam(14) = 'CPO2      '
      Appnam(15) = 'CPO3      '
      Appnam(16) = 'CPO4      '
      Appnam(17) = 'CPO5      '
      Appnam(18) = 'CPO6      '
      Appnam(19) = 'CPO7      '
      Appnam(20) = 'CPO8      '
      Appnam(21) = 'COO1      '
      Appnam(22) = 'COO2      '
      Appnam(23) = 'COO3      '
      Appnam(24) = 'COO4      '
      Appnam(25) = 'COO5      '
      Appnam(26) = 'COO6      '
      Appnam(27) = 'COO7      '
      Appnam(28) = 'COO8      '
      Appnam(29) = 'CNS1      '
      Appnam(30) = 'CNS2      '
      Appnam(31) = 'CNS3      '
      Appnam(32) = 'CNS4      '
      Appnam(33) = 'CNS5      '
      Appnam(34) = 'CNS6      '
      Appnam(35) = 'CNS7      '
      Appnam(36) = 'CNS8      '
      Appnam(37) = 'CBS1      '
      Appnam(38) = 'CBS2      '
      Appnam(39) = 'CBS3      '
      Appnam(40) = 'CBS4      '
      Appnam(41) = 'CBS5      '
      Appnam(42) = 'CAS1      '
      Appnam(43) = 'CAS2      '
      Appnam(44) = 'CAS3      '
      Appnam(45) = 'CAS4      '
      Appnam(46) = 'CAS5      '
c      Appnam(57) = 'RNO3      '
      Appnam(50) = 'PEC_      '
      Appnam(60) = 'POC_      '
      Appnam(70) = 'CRST_     '
      Appnam(80) = 'PCL_      '
      Appnam(90) = 'NA_       '
      Appnam(100) = 'PSO4_     '
      Appnam(110) = 'APO1      '
      Appnam(111) = 'APO2      '
      Appnam(112) = 'APO3      '
      Appnam(113) = 'APO4      '
      Appnam(114) = 'APO5      '
      Appnam(115) = 'APO6      '
      Appnam(116) = 'APO7      '
      Appnam(117) = 'APO8      '
      Appnam(118) = 'AOO1      '
      Appnam(119) = 'AOO2      '
      Appnam(120) = 'AOO3      '
      Appnam(121) = 'AOO4      '
      Appnam(122) = 'AOO5      '
      Appnam(123) = 'AOO6      '
      Appnam(124) = 'AOO7      '
      Appnam(125) = 'AOO8      '
      Appnam(126) = 'ANS1      '
      Appnam(127) = 'ANS2      '
      Appnam(128) = 'ANS3      '
      Appnam(129) = 'ANS4      '
      Appnam(130) = 'ANS5      '
      Appnam(131) = 'ANS6      '
      Appnam(132) = 'ANS7      '
      Appnam(133) = 'ANS8      '
      Appnam(134) = 'ABS1      '
      Appnam(135) = 'ABS2      '
      Appnam(136) = 'ABS3      '
      Appnam(137) = 'ABS4      '
      Appnam(138) = 'ABS5      '
      Appnam(139) = 'AAS1      '
      Appnam(140) = 'AAS2      '
      Appnam(141) = 'AAS3      '
      Appnam(142) = 'AAS4      '
      Appnam(143) = 'AAS5      '
c      Appnam(159) = 'PNO3      '
      Appnam(144) = 'PNH4      '
      !Appnam(161) = 'PCL_      '
c
c-----Assign full names to all the PM species
c
      PM_indx = 50
      do i=1,6
          isblk=INDEX(Appnam((i-1)*10+PM_indx),' ')
        do j=1,10
          Appnam((i-1)*10+(PM_indx-1)+j)=Appnam((i-1)*10+PM_indx)(1:isblk-1)//
     &       num(j)//Appnam((i-1)*10+PM_indx)(isblk+2:10)
        enddo
      enddo
      sa_num_gas = 46
      sa_num_sv  = 110
      sv_bin = 6 - 1
c
c-----Match tracers to modeled species
c
      do i=1,MXTRK
        Appmaprev(i)=0
      enddo
c
      do i=1,nspec
        match = .FALSE.
        Appmap(i)=0
        do j=1,MXTRK
          if (spname(i).eq.Appnam(j)) then 
            Appmap(i)=j
            Appmaprev(j)=i
            match = .TRUE.
          elseif (j.ge.sa_num_sv .and.  !Semivolatile Aerosol
     &          spname(i)(1:4).eq.Appnam(j)(1:4) ) then
            if (spname(i)(5:7).eq.'_1 ') then
              !Match the bulk Source Apportionment of SV 
              !aerosols to the index number of the matching
              !aerosol's lowest size bin
              Appmap(i)=j
              Appmaprev(j)=i
              match = .TRUE.
            elseif (spname(i)(5:7).eq.'_2 '.or.spname(i)(5:7).eq.'_3 '.or.
     &              spname(i)(5:7).eq.'_4 '.or.spname(i)(5:7).eq.'_5 '.or.
     &              spname(i)(5:7).eq.'_6 ') then
              !Match the bulk Source Apportionment of SV 
              !aerosols to the index number of the matching
              !aerosol's lowest size bin
              Appmap(i)=j
              match = .TRUE.
            endif
          endif
        enddo
        !if (.not.match) write(6,*) 'Match not found: ',spname(i)
      enddo
c
c-----Print out the headers for the output file
c     
      end
