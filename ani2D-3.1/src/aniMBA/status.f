C ======================================================================
      Subroutine setStatus(flagAuto, status, iPrint)
C ======================================================================
      include 'status.fd'
C ======================================================================
C Routine adds additional properties to the variable status
C ======================================================================
      Logical  flagAuto
      Integer  status, iPrint

      Logical  ifXnode

C ======================================================================
      status = max(0, status)

c ... remove obsolete and not-implemented input features 
      If(status.GT.0) Then
         If(ifXnode(status, ANISmoothMesh)) Then
            Call delXnode(status, ANISmoothMesh)
            If(iPrint.GE.1) Write(*, 5001) 
         End if

         If(ifXnode(status, ANIMultiConnectedGeometry)) Then
           Call delXnode(status, ANIMultiConnectedGeometry)
            If(iPrint.GE.1) Write(*, 5002) 
         End if
      End if

      
c ... inform the user about requested features  
      If(iPrint.GE.1) Then
         If(ifXnode(status, ANIForbidBoundaryElements)) Write(*, 5003) 
         If(ifXnode(status, ANIFixBoundaryEdges))       Write(*, 5004)
         If(ifXnode(status, ANIFixBoundaryPoints))      Write(*, 5007)

         If(ifXnode(status, ANIDeleteTemporaryEdges))
     &      Write(*, 5006) '[user]'
         If(ifXnode(status, ANIUntangleMesh))
     &      Write(*, 5008) '[user]'
      End if


c ... set up default features
      If(flagAuto) Then
         If(.NOT.ifXnode(status, ANIDeleteTemporaryEdges)) Then
            If(iPrint.GE.2) Write(*, 5006) '[system]'

            Call addXnode(status, ANIDeleteTemporaryEdges)
         End if
      End if

      If(iPrint.GE.1) Write(*,*) 

      Return

 5001 Format('status.fd: -256  [ANISmoothMesh]             [obsolete]')
 5002 Format('status.fd: -64   [ANIMultiConnectedGeometry] [not used]')
 5003 Format('status.fd: +1    [ANIForbidBoundaryElements] [user]') 
 5004 Format('status.fd: +4    [ANIFixBoundaryEdges]       [user]')
 5006 Format('status.fd: +8    [ANIDeleteTemporaryEdges]   ', A)
 5007 Format('status.fd: +16   [ANIFixBoundaryPoints]      [user]')
 5008 Format('status.fd: +32   [ANIUntangleMesh]           ', A)

      Return
      End


