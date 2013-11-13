c
c     App.com contains the basic variables that will be used in the 
c     apportionment algorithm.
c
c     DEVELOPED BY:
c     Kristina Wagstrom
c     Carnegie Mellon University (Chemical Engineering)
c     08/16/2006
c
c-------------------------------------------------------------------------
c
c     Appconc      -- Contains the individual source concentrations
c     Appnum       -- Number of sources to be tracks (not including B.C. & I.C.)
c     Appmap       -- Location of the apportionment vs. modeled species
c     Appnam       -- Names of the apportionment species
c     Appct        -- Counter for averaging
c     Appavg       -- Hourly average concentration
c     AppemisA     -- Mapping between app. species and emissions species - area
c     AppemisP     -- Mapping between app. species and emissions species - point
c     nAppA        -- Number of area emissions species
c     nAppP        -- Number of point emissions species
c     Appmaprev    -- Location of the modeled species vs. apportionment
c     AppemfracA   -- Fraction of area of emissions from each source - area
c     AppemfracP   -- Fraction of area of emissions from each source - point
c
c--------------------------------------------------------------------------
c
      real Appconc(MXVEC3D*MXSOUR*MXTRK)
      real AppFORE(MXCOL1,MXROW1,MXLAY1,MXSPEC)
      real AppAFT(MXCOL1,MXROW1,MXLAY1,MXSPEC)
      real Appavg(MXCOL1*MXROW1*MXLAY1*MXSOUR*MXTRK)
      integer Appnum, Appmap(MXSPEC), Appct
      integer Appmaprev(MXTRK), sa_num_gas, sa_num_sv
      character*10 Appnam(MXTRK)
      integer AppemisP(MXSOUR-3,MXTRK)
      integer AppemisA(MXSOUR-3,MXTRK)
      integer nAppA(MXSOUR-3),nAppP(MXSOUR-3)
      real AppemfracA(MXCOL1,MXROW1,MXSPEC,MXSOUR-2)
      real AppemfracP(120000,MXSPEC,MXSOUR-2)
      real modvdep(MXCOL1,MXROW1,MXSPEC)
      real saconc(MXCOL1,MXROW1,MXLAY1,MXTRK*MXSOUR)
      integer AppEmissNum, AppEmissMap(40)
      character*200 AppEmissFile(40)
c
      common /App/ Appconc,Appnum,Appnam,Appmap,Appct,
     &             Appavg, AppFORE, AppAFT
      common /App2/ nAppA,nAppP,Appmaprev
      common /App3/ AppemfracA,AppemfracP
      common /App4/ AppemisP, AppemisA
      common /App5/ modvdep, sa_num_gas, sa_num_sv
      common /App6/ saconc, sv_bin
      common /App6/ AppEmissNum, AppEmissFile, AppEmissMap
c     
c     END
