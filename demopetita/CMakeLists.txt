####################################################################################
# Binary target: SW
####################################################################################

# Grab all the *.c files to build the library:
FILE(GLOB demopetitaSrc *.c)

add_executable (demopetita ${COMMONsrc} ${demopetitaSrc})
target_link_libraries(demopetita ${LIBM} ${EXTRA_LIBS})

# The library will be installed to lib/ directory:
install (TARGETS demopetita DESTINATION .)

if(BUILD_TESTING)

	add_test (demopetitaRuns demopetita)
	set_tests_properties(demopetitaRuns PROPERTIES PASS_REGULAR_EXPRESSION
                                                "Please run with 4 procs")
	add_test (demopetitaRuns4procs mpirun -np 4 demopetita)
	set_tests_properties(demopetitaRuns4procs PROPERTIES PASS_REGULAR_EXPRESSION
                                                "ite=90 max")
endif(BUILD_TESTING)

