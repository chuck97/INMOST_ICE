#
#  Ani2D 3-X
#  Inexact Newton-Ktylov Backtracking
#

add_library(inb-${VERSION} STATIC
            "bcg_inb.f"  "bt_inb.f"  "jv_inb.f"  "solver_inb.f")

install(TARGETS inb-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




