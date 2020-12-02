####################################################################################
# Binary target: SW
####################################################################################

# Grab all the *.c files to build the library:
FILE(GLOB casegenSrc *.c)

add_executable (casegen ${COMMONsrc} ${casegenSrc})
target_link_libraries(casegen ${LIBM} ${EXTRA_LIBS})

# The library will be installed to lib/ directory:
install (TARGETS casegen DESTINATION .)
