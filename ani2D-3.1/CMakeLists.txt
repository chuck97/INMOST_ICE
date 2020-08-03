cmake_minimum_required(VERSION 2.6)
project(ani2D)
enable_language(Fortran)


# default parameters
set(CMAKE_BUILD_TYPE Release)
set(VERSION 3.1)


set(LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../lib)
set(INCLUDE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../include)
set(RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../bin)
set(PROJECT_DIR ${PROJECT_BINARY_DIR}/..)

set(ADDITIONAL_MAKE_CLEAN_FILES build)


# process input parameters
if (ENABLE_SYSTEM_LAPACK)
  set(LAPACK_LIBS "lapack;blas")
else ()
  set(LAPACK_LIBS "lapack-3.2;lapack_ext-3.2;blas-3.2")
endif()

if (DISABLE_C2F_TESTS)
  set(FLAG_C2F_TESTS 0)
else ()
  set(FLAG_C2F_TESTS 1)
endif()

if (FIX_CMAKE_BUGS)
  set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
endif()


# detect Platform
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows") 
  set(SUFFICS_EXE "")
else ()
  set(SUFFICS_EXE ".exe")
endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "SunOS")
# set(LIBRARY_G2C "g2c") Solaris 11 and earlier
  set(LIBRARY_G2C "")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Windows") 
  set(LIBRARY_G2C "gfortran")
else ()
  set(LIBRARY_G2C "")
endif()


# version of TPLs
set(UMFPACK_VER "5.1")
set(LAPACK_VER "3.2")


###########################################################
# include directories
###########################################################
include_directories(${PROJECT_BINARY_DIR}/../src/aniFEM)
include_directories(${PROJECT_BINARY_DIR}/../src/aniMBA)
include_directories(${PROJECT_BINARY_DIR}/../src/aniC2F)
include_directories(${PROJECT_BINARY_DIR}/../src/aniINB)
include_directories(${PROJECT_BINARY_DIR}/../src/aniLU/UMFPACK/Include)
include_directories(${PROJECT_BINARY_DIR}/../src/aniLU/AMD/Include)
include_directories(${PROJECT_BINARY_DIR}/../src/aniLU/UFconfig)


# list libraries
add_subdirectory(src/aniAFT)
add_subdirectory(src/aniFEM)
add_subdirectory(src/aniILU)
add_subdirectory(src/aniLU)
add_subdirectory(src/aniINB)
add_subdirectory(src/aniLMR)
add_subdirectory(src/aniMBA)
add_subdirectory(src/aniPRJ)
add_subdirectory(src/aniRCB)
add_subdirectory(src/aniC2F)
add_subdirectory(src/aniVIEW)
add_subdirectory(src/blas)
add_subdirectory(src/lapack)
add_subdirectory(src/lapack/extension)


###########################################################
# test packages
###########################################################
enable_testing()

add_subdirectory(src/Tutorials/PackageAFT)
if (FLAG_C2F_TESTS)
  add_subdirectory(src/Tutorials/PackageC2F)
endif()
add_subdirectory(src/Tutorials/PackageFEM)
add_subdirectory(src/Tutorials/PackageILU)
add_subdirectory(src/Tutorials/PackageLU)
add_subdirectory(src/Tutorials/PackageINB)
add_subdirectory(src/Tutorials/PackageLMR)
add_subdirectory(src/Tutorials/PackageMBA)
add_subdirectory(src/Tutorials/PackagePRJ)
add_subdirectory(src/Tutorials/PackageRCB)
add_subdirectory(src/Tutorials/PackageVIEW)

add_subdirectory(python)

add_subdirectory(src/Tutorials/MultiPackage/Interpolation)
add_subdirectory(src/Tutorials/MultiPackage/ConvectionDiffusion)
add_subdirectory(src/Tutorials/MultiPackage/HybridMixedFem)
add_subdirectory(src/Tutorials/MultiPackage/IsotropicBVP)
add_subdirectory(src/Tutorials/MultiPackage/NonisotropicBVP1)
add_subdirectory(src/Tutorials/MultiPackage/NonisotropicBVP2)
add_subdirectory(src/Tutorials/MultiPackage/PosterioriEstimates)
add_subdirectory(src/Tutorials/MultiPackage/Stokes)
add_subdirectory(src/Tutorials/MultiPackage/StokesNavier)
add_subdirectory(src/Tutorials/MultiPackage/UnsteadyConDif)
add_subdirectory(src/Tutorials/MultiPackage/DiscontinuousGalerkin)

add_test(AFT_front_method1 ${RUNTIME_OUTPUT_DIRECTORY}/aniAFT_front1.exe)
add_test(AFT_front_method2 ${RUNTIME_OUTPUT_DIRECTORY}/aniAFT_front2.exe)
add_test(AFT_square        ${RUNTIME_OUTPUT_DIRECTORY}/aniAFT_boundary_square.exe)
add_test(AFT_wing          ${RUNTIME_OUTPUT_DIRECTORY}/aniAFT_boundary_wing.exe)

if (FLAG_C2F_TESTS)
  add_test(C2F_nodal    ${RUNTIME_OUTPUT_DIRECTORY}/aniC2F_nodal.exe)
  add_test(C2F_analytic ${RUNTIME_OUTPUT_DIRECTORY}/aniC2F_analytic.exe)
  add_test(C2F_audit    ${RUNTIME_OUTPUT_DIRECTORY}/aniC2F_audit.exe)
endif()

add_test(FEM_BCs      ${RUNTIME_OUTPUT_DIRECTORY}/aniFEM_BC.exe)
add_test(FEM_Template ${RUNTIME_OUTPUT_DIRECTORY}/aniFEM_template.exe)
add_test(FEM_Triangle ${RUNTIME_OUTPUT_DIRECTORY}/aniFEM_triangle.exe)

add_test(ILU_bcg0    ${RUNTIME_OUTPUT_DIRECTORY}/aniILU_bcg0.exe)
add_test(ILU_bcg2    ${RUNTIME_OUTPUT_DIRECTORY}/aniILU_bcg2.exe)
add_test(ILU_gmares0 ${RUNTIME_OUTPUT_DIRECTORY}/aniILU_gmr0.exe)
add_test(LU          ${RUNTIME_OUTPUT_DIRECTORY}/aniLU.exe)

add_test(INB_simple ${RUNTIME_OUTPUT_DIRECTORY}/aniINB_simple.exe)
add_test(INB_bratu  ${RUNTIME_OUTPUT_DIRECTORY}/aniINB_bratu.exe)

add_test(LMR_estimate ${RUNTIME_OUTPUT_DIRECTORY}/aniLMR_edgeest.exe)
add_test(LMR_gradient ${RUNTIME_OUTPUT_DIRECTORY}/aniLMR_gradient.exe)
add_test(LMR_nodal    ${RUNTIME_OUTPUT_DIRECTORY}/aniLMR_nodal.exe)

add_test(MBA_nodal    ${RUNTIME_OUTPUT_DIRECTORY}/aniMBA_nodal.exe)
add_test(MBA_analytic ${RUNTIME_OUTPUT_DIRECTORY}/aniMBA_analytic.exe)
add_test(MBA_fixShape ${RUNTIME_OUTPUT_DIRECTORY}/aniMBA_fixshape.exe)
add_test(MBA_untangle ${RUNTIME_OUTPUT_DIRECTORY}/aniMBA_untangle.exe)
add_test(MBA_triangle ${RUNTIME_OUTPUT_DIRECTORY}/aniMBA_triangle.exe)

add_test(PRJ_remap ${RUNTIME_OUTPUT_DIRECTORY}/aniPRJ.exe)
add_test(RCB ${RUNTIME_OUTPUT_DIRECTORY}/aniRCB.exe)
add_test(Vizualization_mesh ${RUNTIME_OUTPUT_DIRECTORY}/aniVIEW_mesh.exe)
add_test(Vizualization_matrix ${RUNTIME_OUTPUT_DIRECTORY}/aniVIEW_matrix.exe)

add_test(multi_Interpolation        ${RUNTIME_OUTPUT_DIRECTORY}/multi_Interpolation.exe)
add_test(multi_InterpolationGrad    ${RUNTIME_OUTPUT_DIRECTORY}/multi_InterpolGrad.exe)
add_test(multi_ConvectionDiffusion_long  ${RUNTIME_OUTPUT_DIRECTORY}/multi_ConDifAdap.exe)
add_test(multi_HybridMixedFEM       ${RUNTIME_OUTPUT_DIRECTORY}/multi_HybridMFEM.exe)
add_test(multi_IsotropicHessian     ${RUNTIME_OUTPUT_DIRECTORY}/multi_IsoAdapHessian.exe)
add_test(multi_IsotropicAPosteriori ${RUNTIME_OUTPUT_DIRECTORY}/multi_IsoAdapPosteriori.exe)
add_test(multi_Anisotropic_1        ${RUNTIME_OUTPUT_DIRECTORY}/multi_AnisoAdap1.exe)
add_test(multi_Anisotropic_2_long   ${RUNTIME_OUTPUT_DIRECTORY}/multi_AnisoAdap2.exe)
add_test(multi_APosterioriEstimates ${RUNTIME_OUTPUT_DIRECTORY}/multi_APostEst.exe)
add_test(multi_Stokes               ${RUNTIME_OUTPUT_DIRECTORY}/multi_Stokes.exe)
add_test(multi_StokesNavier_long    ${RUNTIME_OUTPUT_DIRECTORY}/multi_StokesNavier.exe)
add_test(multi_UnsteadyConvDiff     ${RUNTIME_OUTPUT_DIRECTORY}/multi_UnsteadyConDif.exe)
add_test(multi_DG                   ${RUNTIME_OUTPUT_DIRECTORY}/multi_DG.exe)

include (CTest)


###########################################################
# production packages
###########################################################
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
set(CPACK_GENERATOR "TGZ;ZIP")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Advanced Numerical Methods 2D")
set(CPACK_PACKAGE_NAME "ani2D")
set(CPACK_PACKAGE_VERSION "3.1")
set(CPACK_PACKAGE_VERSION_MAJOR "3")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_BINARY_DIR}/../dist")

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_BINARY_DIR}/../COPYRIGHT")
set(CPACK_RESOURCE_FILE_README  "${PROJECT_BINARY_DIR}/../README")

set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "ani2D-3.1")
set(CPACK_SOURCE_IGNORE_FILES "testing;book;lib/lib;doc/;build/;bin/;dist;apps/;TODO;\\\\.svn;\\\\.list;\\\\.ps")

set(ARCHIVE_TGZ "${PROJECT_BINARY_DIR}/_CPack_Packages/Linux-Source/TGZ/${CPACK_SOURCE_PACKAGE_FILE_NAME}") 
set(ARCHIVE_ZIP "${PROJECT_BINARY_DIR}/_CPack_Packages/Linux-Source/ZIP/${CPACK_SOURCE_PACKAGE_FILE_NAME}") 
set(CPACK_INSTALL_COMMANDS 
  "${CMAKE_COMMAND} -E make_directory ${ARCHIVE_TGZ}/bin/;" 
  "${CMAKE_COMMAND} -E touch ${ARCHIVE_TGZ}/bin/README;"
  "${CMAKE_COMMAND} -E make_directory ${ARCHIVE_TGZ}/doc/;" 
  "${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/../doc/user_guide.pdf ${ARCHIVE_TGZ}/doc/;"
  "${CMAKE_COMMAND} -E make_directory ${ARCHIVE_ZIP}/bin/;" 
  "${CMAKE_COMMAND} -E touch ${ARCHIVE_ZIP}/bin/README;"
  "${CMAKE_COMMAND} -E make_directory ${ARCHIVE_ZIP}/doc/;" 
  "${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/../doc/user_guide.pdf ${ARCHIVE_ZIP}/doc/")

include(CPack)

