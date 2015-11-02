
PLUGIN_LIBRARY_NAME = nnls-chroma

PLUGIN_CODE_OBJECTS = chromamethods.o NNLSBase.o NNLSChroma.o Chordino.o Tuning.o plugins.o nnls.o viterbi.o

CFLAGS += $(ARCHFLAGS) $(OPTFLAGS) -I$(VAMP_SDK_DIR) -Wall -fPIC
CXXFLAGS += $(ARCHFLAGS) $(OPTFLAGS) -I$(VAMP_SDK_DIR) -I$(BOOST_ROOT) -Wall -fPIC

PLUGIN = $(PLUGIN_LIBRARY_NAME)$(PLUGIN_EXT)

$(PLUGIN): $(PLUGIN_CODE_OBJECTS)
	   $(CXX) -o $@ $^ $(LDFLAGS)

nnls.o:	nnls.c		# not nnls.f

test:	$(PLUGIN)
	bash regression/regression.sh

clean:
	rm -f *.o

# DO NOT DELETE

Chordino.o: Chordino.h NNLSBase.h chromamethods.h nnls.h viterbi.h
chromamethods.o: chromamethods.h nnls.h
NNLSBase.o: NNLSBase.h chromamethods.h nnls.h
NNLSChroma.o: NNLSChroma.h NNLSBase.h chromamethods.h nnls.h
plugins.o: NNLSChroma.h NNLSBase.h Chordino.h Tuning.h
Tuning.o: Tuning.h NNLSBase.h chromamethods.h nnls.h
viterbi.o: viterbi.h
