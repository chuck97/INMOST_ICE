#
#  Ani2D 3-X
#  Advance Front Technique
#

add_library(aft2D-${VERSION} STATIC
            "memory2.c"  "aft2d.c"    "refine2.c" 
            "region2.c"  "struct2.c"  "tree2.c" 
            "tria2.c"    "user2.c")

install(TARGETS aft2D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




