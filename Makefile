PLUGIN_LIBRARY_NAME = matthiasm

# Edit this to list one .o file for each .cpp file in your plugin project
#
PLUGIN_CODE_OBJECTS = NNLSChroma.o plugins.o nnls.o

# Edit this to the location of the Vamp plugin SDK, relative to your
# project directory
#
VAMP_SDK_DIR = ../vamp-plugin-sdk
QMDSP_DIR = ../qm-dsp/build/osx/20091028
BOOST_ROOT = ../boost_1_44_0


##  Uncomment these for an OS/X native build using command-line tools:
ARCHFLAGS = -isysroot /Developer/SDKs/MacOSX10.4u.sdk -mmacosx-version-min=10.4 -arch ppc -arch i386
CFLAGS = $(ARCHFLAGS) -Wall -fPIC -g
CXXFLAGS = $(CFLAGS) -I$(VAMP_SDK_DIR) -I$(BOOST_ROOT)
PLUGIN_EXT = .dylib
PLUGIN = $(PLUGIN_LIBRARY_NAME)$(PLUGIN_EXT)
LDFLAGS = $(ARCHFLAGS) -dynamiclib -install_name $(PLUGIN) $(VAMP_SDK_DIR)/libvamp-sdk.a  -exported_symbols_list vamp-plugin.list -framework Accelerate


$(PLUGIN): $(PLUGIN_CODE_OBJECTS)
	   $(CXX) -o $@ $^ $(LDFLAGS)

NNLSChroma.o: NNLSChroma.h
plugins.o: NNLSChroma.h

nnls.o:	nnls.c		# not nnls.f

clean:
	rm -f *.o

