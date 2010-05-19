##  Skeleton Makefile for Vamp plugin builds using command-line tools.
##
##  Rename this to Makefile, and edit as appropriate.
##  This Makefile WILL NOT WORK until you have edited it as described
##  below -- the Makefile as supplied does nothing useful at all!
##
##  Various sets of options are provided, commented out -- just uncomment
##  (remove the '#' characters for) the set that most closely resembles
##  your own situation, and adjust to taste.  Then run "make".
##
##  (For Windows builds using MS Visual Studio, start instead with the
##  VampExamplePlugins project found in the build directory of the SDK.)


# Edit this to the base name of your plugin library
#
PLUGIN_LIBRARY_NAME = matthiasm

# Edit this to list one .o file for each .cpp file in your plugin project
#
PLUGIN_CODE_OBJECTS = NNLSChroma.o plugins.o nnls.o

# Edit this to the location of the Vamp plugin SDK, relative to your
# project directory
#
VAMP_SDK_DIR = ../vamp-plugin-sdk
# LAPACK_DIR = ../lapack
QMDSP_DIR = ../qm-dsp/build/osx/20091028
FFT_DIR = ../qm-dsp/dsp/transforms
NNLS_DIR = ../nnls_suvrit



##  Uncomment these for an OS/X native build using command-line tools:
CXXFLAGS = -arch i386 -I$(VAMP_SDK_DIR) -I$(FFT_DIR) -I$(NNLS_DIR) -Wall -fPIC -g
PLUGIN_EXT = .dylib
PLUGIN = $(PLUGIN_LIBRARY_NAME)$(PLUGIN_EXT)
LDFLAGS = -g -m32 -dynamiclib -install_name $(PLUGIN) $(VAMP_SDK_DIR)/libvamp-sdk.a  -exported_symbols_list vamp-plugin.list -framework Accelerate


##  Uncomment these for an OS/X universal binary using command-line tools:

# CXXFLAGS = -isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch i386 -arch ppc -I$(VAMP_SDK_DIR) -Wall -fPIC
# PLUGIN_EXT = .dylib
# PLUGIN = $(PLUGIN_LIBRARY_NAME)$(PLUGIN_EXT)
# LDFLAGS = -dynamiclib -install_name $(PLUGIN) $(VAMP_SDK_DIR)/libvamp-sdk.a -exported_symbols_list vamp-plugin.list


##  Uncomment these for Linux using the standard tools:

# CXXFLAGS = -I$(VAMP_SDK_DIR) -Wall -fPIC
# PLUGIN_EXT = .so
# PLUGIN = $(PLUGIN_LIBRARY_NAME)$(PLUGIN_EXT)
# LDFLAGS = -shared -Wl,-soname=$(PLUGIN) $(VAMP_SDK_DIR)/libvamp-sdk.a -Wl,--version-script=vamp-plugin.map


##  Uncomment these for a cross-compile from Linux to Windows using MinGW:

# CXX = i586-mingw32msvc-g++
# CXXFLAGS = -I$(VAMP_SDK_DIR) -Wall 
# PLUGIN_EXT = .dll
# PLUGIN = $(PLUGIN_LIBRARY_NAME)$(PLUGIN_EXT)
# LDFLAGS = --static-libgcc -Wl,-soname=$(PLUGIN) -shared $(VAMP_SDK_DIR)/libvamp-sdk.a


##  Uncomment these for OpenSolaris using SunStudio compiler and GNU make:

# CXX = CC
# CXXFLAGS = -G -I$(VAMP_SDK_DIR) +w -KPIC
# PLUGIN_EXT = .so
# PLUGIN = $(PLUGIN_LIBRARY_NAME)$(PLUGIN_EXT)
# LDFLAGS = -G -h$(PLUGIN) $(VAMP_SDK_DIR)/libvamp-sdk.a -Qoption ld -Mvamp-plugin.map



##  All of the above

$(PLUGIN): $(PLUGIN_CODE_OBJECTS)
	   $(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o

