####################################################################################
# libsppde.a CMakeFileList
####################################################################################

# Grab all the *.c files to build the library:
FILE(GLOB sppdesrc *.c)

# We will build an STATIC version:
add_library(sppde STATIC ${sppdesrc})

# The library will be installed to lib/ directory:
install (TARGETS sppde DESTINATION lib)
