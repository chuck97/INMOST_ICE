#
#  Ani2D 3-X
#  Local Metric Recovery
#

add_library(lmr2D-${VERSION} STATIC
           "Func2MetricMAX.f"     "Func2GradMetricMAX.f"
           "Nodal2MetricVAR.f"    "Nodal2MetricZZ.f"       
	   "EdgeEst2MetricMAX.f"  "EdgeEst2GradMetricMAX.f" 
           "EdgeEst2MetricLS.f"   "CellEst2MetricZZ.f"     
	   "metric_tri.f"  
           "Lp_norm.f"  "metric_norm.f"  "utils.f"  "error.f"
           "../aniIO/io.f")

install(TARGETS lmr2D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




