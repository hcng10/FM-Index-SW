XCC=g++

HEADERS = ../fm-sa-def.h file_op.h format.h buildFM.h
HEADER_CPPS = file_op.cpp format.cpp buildFM.cpp
LINKDIVSUF=-Wl,-R divsufsort64/lib/  -L divsufsort64/lib/ -ldivsufsort64
INCDIVSUF=-I divsufsort64/include -ldivsufsort64 -L divsufsort64/lib/ -ldivsufsort64

FLAGS = -fopenmp

FM_SA_LIST = fm-sa-build
	
all: $(FM_SA_LIST)

fm-sa-build: fm-sa-build.cpp $(HEADER_CPPS) $(HEADERS)
	$(XCC) -o $@ $< \
	$(HEADER_CPPS) \
	$(HEADERS) \
	$(FLAGS) \
	$(INCDIVSUF) 
	
.PHONY: all clean

clean:
	rm -f $(FM_SA_LIST)
