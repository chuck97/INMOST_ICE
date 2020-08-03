#
#  Ani2D 3-X
#  Interface to C-code
#

add_library(c2f2D-${VERSION} STATIC
            "ani2D.c"  "ani2D.h")

install(TARGETS c2f2D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES "ani2D.h" DESTINATION ${INCLUDE_OUTPUT_DIRECTORY})





