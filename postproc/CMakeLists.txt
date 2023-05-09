####################################################################################
# Binary target: SW
####################################################################################

# Grab all the *.c files to build the library:
FILE(GLOB postprocSrc *.c)

add_executable (postproc ${COMMONsrc} ${postprocSrc})
target_link_libraries(postproc ${LIBM} ${EXTRA_LIBS})

# The library will be installed to lib/ directory:
install (TARGETS postproc DESTINATION .)
