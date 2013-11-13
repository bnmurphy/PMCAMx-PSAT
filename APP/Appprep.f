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
      integer i,j
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
      Appnam(65) = 'PEC_      '
      Appnam(75) = 'POC_      '
      Appnam(85) = 'CRST_     '
      Appnam(95) = 'PCL_      '
      Appnam(105) = 'NA_       '
      Appnam(115) = 'PSO4_     '
      Appnam(125) = 'APO1_     '
      Appnam(135) = 'APO2_     '
      Appnam(145) = 'APO3_     '
      Appnam(155) = 'APO4_     '
      Appnam(165) = 'APO5_     '
      Appnam(175) = 'APO6_     '
      Appnam(185) = 'APO7_     '
      Appnam(195) = 'APO8_     '
      Appnam(205) = 'AOO1_     '
      Appnam(215) = 'AOO2_     '
      Appnam(225) = 'AOO3_     '
      Appnam(235) = 'AOO4_     '
      Appnam(245) = 'AOO5_     '
      Appnam(255) = 'AOO6_     '
      Appnam(265) = 'AOO7_     '
      Appnam(275) = 'AOO8_     '
      Appnam(285) = 'ANS1_     '
      Appnam(295) = 'ANS2_     '
      Appnam(305) = 'ANS3_     '
      Appnam(315) = 'ANS4_     '
      Appnam(325) = 'ANS5_     '
      Appnam(335) = 'ANS6_     '
      Appnam(345) = 'ANS7_     '
      Appnam(355) = 'ANS8_     '
      Appnam(365) = 'ABS1_     '
      Appnam(375) = 'ABS2_     '
      Appnam(385) = 'ABS3_     '
      Appnam(395) = 'ABS4_     '
      Appnam(405) = 'ABS5_     '
      Appnam(415) = 'AAS1_     '
      Appnam(425) = 'AAS2_     '
      Appnam(435) = 'AAS3_     '
      Appnam(445) = 'AAS4_     '
      Appnam(455) = 'AAS5_     '
      Appnam(465) = 'PNO3_     '
      Appnam(475) = 'PNH4_     '
c
c-----Assign full names to all the PM species
c
      do i=1,42
          isblk=INDEX(Appnam(i*10+55),' ')
        do j=1,10
          Appnam(i*10+54+j)=Appnam(i*10+55)(1:isblk-1)//
     &       num(j)//Appnam(i*10+55)(isblk+2:10)
        enddo
      enddo
      sa_num_gas = 56
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
          endif
        enddo
c        if (.not.match) write(6,*) 'Match not found: ',spname(i)
      enddo
c
c-----Print out the headers for the output file
c     
      end
