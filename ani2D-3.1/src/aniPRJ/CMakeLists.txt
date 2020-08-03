#
#  Ani2D 3-X
#  Compatible ReMap
#

add_library(prj2D-${VERSION} STATIC
            "metamesh.f"  "rhs.f"  "tri2tri.f"
            "../aniIO/io.f")

install(TARGETS prj2D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




