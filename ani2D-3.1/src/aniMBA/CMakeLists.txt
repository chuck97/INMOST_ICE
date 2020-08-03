#
#  Ani2D 3-X
#  Metric Based Adaptation
#

set(CORE "mba_analytic.f"  "mba_nodal.f"  "mba_fixshape.f" 
         "ani2.f"  "forlibmba.f"  "calcrv.f" 
         "makM.f"  "mesh_audit.f"  "makQ_linear.f"  "module_hessian.f"  
         "list_new.f" "makS.f"  "prjcrv.f"   "lintrp2D.f"   "dsort.f" 
         "nlnfnc.f"  "minim.f" 
         "colapse_edgeF1.f"  "colapse_edgeF2.f"  "insrt_point.f"  "swap_edge.f"  
         "move_point.f"  "split_triangle.f" 
         "auxSP.f"  "auxSF.f"  "tangled.f"  "check.f"  "status.f" "control.f"
         "triangle.f" "edge.f" 
         "update.f"  "maps.f"  "utils.f"
         "refine.f"  "smoothing.f"  "untangle.f" 
         "ZZ.f"  "delaunay.f"  "datacopy.f"
         "statistics.f"  "error.f"  "time.f" "isnan.c"  "random.c" 
         "../aniIO/io.f")

set(IO  "load_mesh.f"  "loadM_other.f"  "saveM.f"  "saveM_other.f")

set(DEBUG  "debug.f")


add_library(mba2D-${VERSION} STATIC ${CORE} ${IO} ${DEBUG})

install(TARGETS mba2D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES "status.fd" DESTINATION ${INCLUDE_OUTPUT_DIRECTORY})





