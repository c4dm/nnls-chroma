PLUGIN_LIBRARY_NAME = matthiasm

# Edit this to list one .o file for each .cpp file in your plugin project
#
PLUGIN_CODE_OBJECTS = NNLSBase.o NNLSChroma.o Chordino.o Tuning.o plugins.o nnls.o chromamethods.o viterbi.o

# Edit this to the location of the Vamp plugin SDK, relative to your
# project directory
#
VAMP_SDK_DIR = ../vamp-plugin-sdk
BOOST_ROOT = ../boost_1_44_0


##  Uncomment these for an OS/X native build using command-line tools:
ARCHFLAGS = -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5 -arch i386
CFLAGS = $(ARCHFLAGS) -Wall -fPIC -g
CXXFLAGS = $(CFLAGS) -I$(VAMP_SDK_DIR) -I$(BOOST_ROOT)
PLUGIN_EXT = .dylib
PLUGIN = $(PLUGIN_LIBRARY_NAME)$(PLUGIN_EXT)
LDFLAGS = $(ARCHFLAGS) -dynamiclib -install_name $(PLUGIN) $(VAMP_SDK_DIR)/libvamp-sdk.a  -exported_symbols_list vamp-plugin.list -framework Accelerate


$(PLUGIN): $(PLUGIN_CODE_OBJECTS)
	   $(CXX) -o $@ $^ $(LDFLAGS)

nnls.o:	nnls.c		# not nnls.f

clean:
	rm -f *.o

# DO NOT DELETE

nnls.o: nnls.h
Chordino.o: Chordino.h NNLSBase.h chromamethods.h nnls.h
chromamethods.o: chromamethods.h nnls.h chorddict.cpp
NNLSBase.o: NNLSBase.h chromamethods.h nnls.h
NNLSChroma.o: NNLSChroma.h NNLSBase.h chromamethods.h nnls.h
plugins.o: NNLSChroma.h NNLSBase.h Chordino.h Tuning.h
Tuning.o: Tuning.h NNLSBase.h chromamethods.h nnls.h
Chordino.o: NNLSBase.h
chromamethods.o: nnls.h
NNLSChroma.o: NNLSBase.h
Tuning.o: NNLSBase.h
viterbi.o: viterbi.h
