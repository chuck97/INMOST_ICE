#
#  Ani2D 3-X
#  VIzualization 
#

add_library(view2D-${VERSION} STATIC
            "demo.f"       "isoline.f"  "mesh.f"  
            "mesh_full.f"  "streamlines.f"
            "matrix.f")

install(TARGETS view2D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




