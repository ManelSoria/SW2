####################################################################################
# Binary target: SW
####################################################################################

# Grab all the *.c files to build the library:
FILE(GLOB mainSrc *.c)

add_executable (SW ${COMMONsrc} ${mainSrc})
target_link_libraries(SW ${LIBM} ${EXTRA_LIBS})

# The library will be installed to lib/ directory:
install (TARGETS SW DESTINATION .)
