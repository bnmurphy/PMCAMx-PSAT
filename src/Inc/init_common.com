c----BNM
c  
c     init_common.COM initializes common block variables that are causing problems
c       during linking because they are initialized more than once
c                            
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c            
c     Modifications:  
c        none  
c 
c-----------------------------------------------------------------------

      data crads /'OH','NO3','HO2'/ 