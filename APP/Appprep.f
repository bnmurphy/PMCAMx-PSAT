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
      Appnam(5)  = 'OLE       '
      Appnam(6)  = 'PAR       '
      Appnam(7)  = 'TOL       '
      Appnam(8)  = 'CRES      '
      Appnam(9)  = 'PNA       '
      Appnam(10) = 'HONO      '
      Appnam(11) = 'HNO3      '
      Appnam(12) = 'ISOP      '
      Appnam(13) = 'ISPD      '
      Appnam(14) = 'NTR       '
      Appnam(15) = 'SO2       '
      Appnam(16) = 'SULF      '
      Appnam(17) = 'NH3       '
      Appnam(18) = 'OLE2      '
      Appnam(19) = 'CG1       '
      Appnam(20) = 'CG2       '
      Appnam(21) = 'CG3       '
      Appnam(22) = 'CG4       '
      Appnam(23) = 'XYL       '
      Appnam(24) = 'HCL       '
      Appnam(25) = 'SOA1_     '
      Appnam(35) = 'SOA2_     '
      Appnam(45) = 'SOA3_     '
      Appnam(55) = 'SOA4_     '
      Appnam(65) = 'POC_      '
      Appnam(75) = 'PEC_      '
      Appnam(85) = 'CRST_     '
      Appnam(95) = 'PH2O_     '
      Appnam(105)= 'PCL_      '
      Appnam(115)= 'NA_       '
      Appnam(125)= 'PNH4_     '
      Appnam(135)= 'PNO3_     '
      Appnam(145)= 'PSO4_     '
c
c-----Assign full names to all the PM species
c
      do i=1,13
          isblk=INDEX(Appnam(i*10+15),' ')
        do j=1,10
          Appnam(i*10+14+j)=Appnam(i*10+15)(1:isblk-1)//
     &       num(j)//Appnam(i*10+15)(isblk+2:10)
        enddo
      enddo
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
