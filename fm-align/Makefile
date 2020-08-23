XCC=g++

HEADERS = ../fm-sa-def.h file_op.h readFM.h reads.h writeSAM.h
HEADER_CPPS = file_op.cpp readFM.cpp reads.cpp align.hpp writeSAM.cpp


FM_SA_LIST = fm-sa-align

FLAGS = -fopenmp  -std=c++11
	
all: $(FM_SA_LIST)

fm-sa-align: fm-sa-align.cpp $(HEADER_CPPS) $(HEADERS)
	$(XCC) -o $@ $< \
	$(HEADER_CPPS) \
	$(HEADERS) \
	$(FLAGS)
	
.PHONY: all clean

clean:
	rm -f $(FM_SA_LIST)
