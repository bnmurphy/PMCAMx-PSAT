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
      integer i,j,ISUND
      character*2 num(10)
      logical match
c
      data num /'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ','10'/
c
      count = 0
c
c-----Assigning species to the apportionment name array
c
      Appnam(1)  = 'NO        '
      Appnam(2)  = 'NO2       '
      Appnam(3)  = 'PAN       '
      Appnam(4)  = 'NXOY      '
      Appnam(5)  = 'OLE1      '
      Appnam(6)  = 'OLE2      '
      Appnam(7)  = 'TERP      '
      Appnam(8)  = 'ALK5      '
      Appnam(9)  = 'ARO1      '
      Appnam(10) = 'HONO      '
      Appnam(11) = 'HNO3      '
      Appnam(12) = 'ISOP      '
      Appnam(13) = 'PBZN      '
      Appnam(14) = 'PAN2      '
      Appnam(15) = 'SO2       '
      Appnam(16) = 'SULF      '
      Appnam(17) = 'NH3       '
      Appnam(18) = 'MPAN      '
      Appnam(19) = 'HCL       '
      Appnam(20) = 'HNO4      '
      Appnam(21) = 'ARO2      '
      Appnam(22) = 'CPO1      '
      Appnam(23) = 'CPO2      '
      Appnam(24) = 'CPO3      '
      Appnam(25) = 'CPO4      '
      Appnam(26) = 'CPO5      '
      Appnam(27) = 'CPO6      '
      Appnam(28) = 'CPO7      '
      Appnam(29) = 'CPO8      '
      Appnam(30) = 'COO1      '
      Appnam(31) = 'COO2      '
      Appnam(32) = 'COO3      '
      Appnam(33) = 'COO4      '
      Appnam(34) = 'COO5      '
      Appnam(35) = 'COO6      '
      Appnam(36) = 'COO7      '
      Appnam(37) = 'COO8      '
      Appnam(38) = 'CNS1      '
      Appnam(39) = 'CNS2      '
      Appnam(40) = 'CNS3      '
      Appnam(41) = 'CNS4      '
      Appnam(42) = 'CNS5      '
      Appnam(43) = 'CNS6      '
      Appnam(44) = 'CNS7      '
      Appnam(45) = 'CNS8      '
      Appnam(46) = 'CBS1      '
      Appnam(47) = 'CBS2      '
      Appnam(48) = 'CBS3      '
      Appnam(49) = 'CBS4      '
      Appnam(50) = 'CBS5      '
      Appnam(51) = 'CAS1      '
      Appnam(52) = 'CAS2      '
      Appnam(53) = 'CAS3      '
      Appnam(54) = 'CAS4      '
      Appnam(55) = 'CAS5      '
      Appnam(56) = 'ALK4      '
      Appnam(57) = 'RNO3      '
      Appnam(65) = 'PEC_      '
      Appnam(75) = 'POC_      '
      Appnam(85) = 'CRST_     '
      Appnam(95) = 'PCL_      '
      Appnam(105) = 'NA_       '
      Appnam(115) = 'PSO4_     '
      Appnam(125) = 'APO1      '
      Appnam(126) = 'APO2      '
      Appnam(127) = 'APO3      '
      Appnam(128) = 'APO4      '
      Appnam(129) = 'APO5      '
      Appnam(130) = 'APO6      '
      Appnam(131) = 'APO7      '
      Appnam(132) = 'APO8      '
      Appnam(133) = 'AOO1      '
      Appnam(134) = 'AOO2      '
      Appnam(135) = 'AOO3      '
      Appnam(136) = 'AOO4      '
      Appnam(137) = 'AOO5      '
      Appnam(138) = 'AOO6      '
      Appnam(139) = 'AOO7      '
      Appnam(140) = 'AOO8      '
      Appnam(141) = 'ANS1      '
      Appnam(142) = 'ANS2      '
      Appnam(143) = 'ANS3      '
      Appnam(144) = 'ANS4      '
      Appnam(145) = 'ANS5      '
      Appnam(146) = 'ANS6      '
      Appnam(147) = 'ANS7      '
      Appnam(148) = 'ANS8      '
      Appnam(149) = 'ABS1      '
      Appnam(150) = 'ABS2      '
      Appnam(151) = 'ABS3      '
      Appnam(152) = 'ABS4      '
      Appnam(153) = 'ABS5      '
      Appnam(154) = 'AAS1      '
      Appnam(155) = 'AAS2      '
      Appnam(156) = 'AAS3      '
      Appnam(157) = 'AAS4      '
      Appnam(158) = 'AAS5      '
      Appnam(159) = 'PNO3      '
      Appnam(160) = 'PNH4      '
      !Appnam(161) = 'PCL_      '
c
c-----Assign full names to all the PM species
c
      do i=1,6
          isblk=INDEX(Appnam(i*10+55),' ')
        do j=1,10
          Appnam(i*10+54+j)=Appnam(i*10+55)(1:isblk-1)//
     &       num(j)//Appnam(i*10+55)(isblk+2:10)
        enddo
      enddo
      sa_num_gas = 57
      sa_num_sv  = 125
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
            ISUND=INDEX(spname(i),'_')
            if (spname(i)(ISUND:ISUND+2).eq.'_1 ') then
              !Match the bulk Source Apportionment of SV 
              !aerosols to the index number of the matching
              !aerosol's lowest size bin
              Appmap(i)=j
              Appmaprev(j)=i
              match = .TRUE.
            elseif (spname(i)(5:7).eq.'_2 '.or.spname(i)(5:7).eq.'_3 '.or.
     &              spname(i)(5:7).eq.'_4 '.or.spname(i)(5:7).eq.'_5 '.or.
     &              spname(i)(5:7).eq.'_6 ' ) then
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
