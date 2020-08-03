#
#  Ani2D 3.X
#  Basic linear algebra package LAPACK-3.2
#

add_library(lapack-3.2 STATIC
            "dsyev.f"   "ilaenv.f"  "dlansy.f" 
            "dlascl.f"  "dsytrd.f"  "dsterf.f"  "dorgtr.f" 
            "dsteqr.f"  "ieeeck.f"  "dlassq.f"  "dlatrd.f" 
            "dsytd2.f"  "dlanst.f"  "dlae2.f"   "dlapy2.f" 
            "dlasrt.f"  "dorgql.f"  "dorgqr.f"  "dlaset.f" 
            "dlaev2.f"  "dlasr.f"   "dlartg.f"  "dlarfg.f" 
            "dorg2l.f"  "dlarft.f"  "dorg2r.f"  "dlarfb.f" 
            "dlarf.f"   "dsysv.f"   "dsytrf.f"  "dsytrs.f" 
	    "dsytf2.f"  "dlasyf.f"  "disnan.f"  
	    "dgesv.f"   "dgetrf.f"  "dgetrs.f"  "iparmq.f"
	    "dgetf2.f"  "dlaswp.f"  "iladlr.f"  "iladlc.f" 
	    "dpotrf.f"  "dpotf2.f"  "dposv.f"   "dpotrs.f"
            "dlaisnan.f")

install(TARGETS lapack-3.2 DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




