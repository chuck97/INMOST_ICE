C ======================================================================
      Subroutine errMesg(iERR, routine, message)
C ======================================================================
C Generic routine for terminating the code with an error message.
C ======================================================================
      Integer iERR
      Character*(*) routine, message
C ======================================================================

      Write(*,5000) iERR, routine, message
      Stop 911
 5000 Format(/,'[MBA] Error:', I7,/, 
     &         '[MBA] Routine: ', A,/, '[MBA] Message: ', A,/)

      Return
      End




