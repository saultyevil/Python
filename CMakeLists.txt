#
# Python CMakeLists
# -----------------
#
# One should also be able to build Python using CMakeLists, if that is the
# poison which they prefer. This CMakeLists will also allow once to import
# Python into an IDE which uses CMake as its build system.
#
# Note that this CMakeLists program is nowhere near as complex as the standard
# Makefile is, as an example it CANNOT indent files automatically nor does it
# create libatomic or change the compiler flags depending on the GCC version
#
# TODO: make building from scratch via CMake possible
#       - get version from Makefile
#       - build GSL libraries if required
#       - create source files as libraries to re-use objects
#       - check minimum cmake version required
#       - check if MPI library is installed, otherwise use gcc or CC
#       - check correct GSL is installed in system and use that instead
#

project(Python C)
cmake_minimum_required(VERSION 3.0)

# Set the Python version

set(PYTHON_VERSION 84g)
message(STATUS "Python version : ${PYTHON_VERSION}")

# The C standard is being set to 99, although if we follow this standard is
# somewhat a mystery. The compiler is, by default, set to mpicc as we want to
# exploit MPI to reduce runtimes. The compilation flags used here work both for
# clang and gcc

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_COMPILER mpicc)
add_definitions(-Wall -Wextra -O3 -D MPI_ON)

# Setup any include directories, as well as any library directories. This is
# used to link the GSL directories

include_directories(include)
link_directories(lib)

# Create the version.h file. During the configure process, the file version.h
# should be created and whenever the program is build again, the information
# in version.h should be updated.

set(GIT_HASH "")
set(GIT_NUM_CHANGES 0)
set(VERSION_FILE source/version.h)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
    find_package(Git)
    find_program(BASH bash HINTS /bin)

    if(GIT_FOUND)
        execute_process(
                COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                OUTPUT_VARIABLE GIT_HASH
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        execute_process(
                COMMAND ${BASH} "-c" "${GIT_EXECUTABLE} status --porcelain | grep '^ M' | wc -l"
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                OUTPUT_VARIABLE GIT_NUM_CHANGES
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    else(GIT_FOUND)
        set(GIT_HASH 0)
        set(GIT_NUM_CHANGES 0)
    endif(GIT_FOUND)

    message(STATUS "Git hash       : ${GIT_HASH}")
    message(STATUS "Git changes    : ${GIT_NUM_CHANGES}")

    file(WRITE ${VERSION_FILE} "#define VERSION \"${PYTHON_VERSION}\"\n\
#define GIT_COMMIT_HASH \"${GIT_HASH}\"\n\
#define GIT_DIFF_STATUS \"${GIT_NUM_CHANGES}\"\n"
     )
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)

# Define the source for the various different programs contained with the
# Python repository -- that is Python, Py_wind and Windsave2table

set(PYTHON_SOURCE
        source/agn.c source/atomicdata.c source/atomicdata_init.c
        source/atomicdata_sub.c source/anisowind.c source/bands.c source/bb.c
        source/bilinear.c source/brem.c source/cdf.c source/charge_exchange.c
        source/compton.c source/continuum.c source/cooling.c source/corona.c
        source/cv.c source/cylind_var.c source/cylindrical.c source/density.c
        source/diag.c source/dielectronic.c source/direct_ion.c source/disk.c
        source/disk_init.c source/emission.c source/estimators.c source/estimators_simple.c
        source/extract.c source/frame.c source/get_models.c source/gradv.c
        source/gridwind.c source/homologous.c source/hydro_import.c
        source/import.c source/import_calloc.c source/import_cylindrical.c
        source/import_rtheta.c source/import_spherical.c source/ionization.c
        source/knigge.c source/levels.c source/lines.c source/macro_gov.c
        source/matom.c source/matom_diag.c source/matrix_ion.c
        source/para_update.c source/parse.c source/partition.c source/paths.c
        source/phot_util.c source/photo_gen_matom.c source/photon2d.c
        source/photon_gen.c source/pi_rates.c source/radiation.c source/random.c
        source/rdpar.c source/rdpar_init.c source/recipes.c source/recomb.c
        source/reposition.c source/resonate.c source/reverb.c source/roche.c
        source/rtheta.c source/run.c source/saha.c source/setup.c
        source/setup_disk.c source/setup_domains.c source/setup_files.c
        source/setup_line_transfer.c source/setup_reverb.c
        source/setup_star_bh.c source/shell_wind.c source/signal.c
        source/spectra.c source/spectral_estimators.c source/spherical.c
        source/stellar_wind.c source/sv.c source/synonyms.c source/time.c
        source/trans_phot.c source/vvector.c source/walls.c source/wind.c
        source/wind2d.c source/wind_sum.c source/wind_updates2d.c
        source/wind_util.c source/windsave.c source/windsave2table_sub.c
        source/xlog.c source/zeta.c
)

set(PY_WIND_SOURCE
        source/py_wind_ion.c source/py_wind_macro.c
        source/py_wind_sub.c source/py_wind_write.c
)

set(WINDSAVE2TABLE_SOURCE
        source/windsave2table_sub.c
)

set(RAD_HYDRO_SOURCE
        source/rad_hydro_files.c
)

# Create the executables for each program, and link the required libraries

# Py_wind
add_executable(py_wind source/py_wind.c ${PYTHON_SOURCE} ${PY_WIND_SOURCE})
target_link_libraries(py_wind gsl gslcblas m)

# Windsave2table
add_executable(windsave2table source/windsave2table.c ${PYTHON_SOURCE} ${WINDSAVE2TABLE_SOURCE})
target_link_libraries(windsave2table gsl gslcblas m)

# Rad_hydro_files
add_executable(rad_hydro_files ${PYTHON_SOURCE} ${RAD_HYDRO_SOURCE})
target_link_libraries(rad_hydro_files gsl gslcblas m)

# Python
add_executable(python source/python.c ${PYTHON_SOURCE})
target_link_libraries(python gsl gslcblas m)
