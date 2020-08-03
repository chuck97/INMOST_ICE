#
#  Ani2D 3-X
#  Finite Element Method
#

add_library(fem2D-${VERSION} STATIC
            "template.f" "assemble.f" 
            "fem2Dtri.f" "fem2Derr.f" "fem2Dedge.f" "fem2Dsub.f"
            "opIDEN.f"   "opGRAD.f"   "opDIV.f"     "opCURL.f"   "opDUDX.f" 
            "opTRACE.f"  "opDUDN.f"
            "bc.f"       "algebra.f"  "forlibfem.f" "isys.f"     "dof.f"
            "utils.f"    "error.f"    "DG.f"
            "../aniIO/io.f")

target_link_libraries(fem2D-${VERSION} ${LAPACK_LIBS})

install(TARGETS fem2D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES "assemble.fd" "fem2Dtri.fd" DESTINATION ${INCLUDE_OUTPUT_DIRECTORY})





