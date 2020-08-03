C ======================================================================
      Subroutine populateControl(nE, control)
C ======================================================================
      include 'status.fd'
      Integer nE, control(6)

C size of the basket (old name is MaxSkipE)
      control(1) = max(100, nE / 15)     

C number of local mesh modifications (MaxQItr)
      control(2) = min(1000, nE * 10)       

c advanced control of mesh generation (status)
      control(3) = ANIForbidBoundaryElements  

c positive = recovery of missing objects such as boundary edges (flagAuto)
      control(4) = 1    

c verbosity level
      control(5) = 3 

c reserved
      control(6) = 0    

      Return
      End


