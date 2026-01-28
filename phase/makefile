#COMPILER MODE C++11
CXX=g++ -std=c++17


#create folders
dummy_build_folder_bin := $(shell mkdir -p bin)
dummy_build_folder_obj := $(shell mkdir -p obj)

#COMPILER & LINKER FLAGS
CXXFLAG=-O3 -mavx2 -mfma -msse3 -Wno-ignored-attributes
LDFLAG=-O3
LDFLAG += -Wl,-rpath,'$$ORIGIN/../../boost/lib'
#CXXFLAG=-O0 -g -mavx2 -mfma -Wno-ignored-attributes
#LDFLAG=-O0

#DYNAMIC LIBRARIES
DYN_LIBS=-lz -lpthread -lbz2 -llzma -lcurl -lcrypto -fopenmp

HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

NAME=$(shell basename $(CURDIR))
BFILE=bin/GLIMPSE2_$(NAME)
EXEFILE=bin/GLIMPSE2_$(NAME)_static

COMMIT_VERS=$(shell git rev-parse --short HEAD)
COMMIT_DATE=$(shell git log -1 --format=%cd --date=short)
CXXFLAG+= -D__COMMIT_ID__=\"$(COMMIT_VERS)\"
CXXFLAG+= -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"

# xSqueezeIt Support [YES/NO]
ifeq ($(XSI_SUPPORT),)
	XSI_SUPPORT=NO
endif
ifeq ($(XSI_SUPPORT),YES)
	XSI_INC=xsqueezeit/include
	XSI_LIB=xsqueezeit/libxsqueezeit.a
	DYN_LIBS+= -lzstd -lhts
	CXXFLAG+= -D__XSI__ -I$(XSI_INC)
endif

# BGEN Support [YES/NO]
ifeq ($(BGEN_SUPPORT),)
	BGEN_SUPPORT=NO
endif

ifeq ($(BGEN_SUPPORT),)
	BGEN_SUPPORT=NO
	BGEN_LIB=""
endif
ifeq ($(BGEN_SUPPORT),YES)
	BGEN_INC=../3rd_party/bgen/genfile/include/
	BGEN_LIB=../3rd_party/bgen/build/libbgen.a ../3rd_party/zstd/lib/libzstd.a
	#DYN_LIBS+=../3rd_party/zstd/lib/libzstd.a
	CXXFLAG+= -D__BGEN__ -I$(BGEN_INC)
endif

BOOST_INC=../boost/include
BOOST_LIB_IO=../boost/lib/libboost_iostreams.so
BOOST_LIB_PO=../boost/lib/libboost_program_options.so
BOOST_LIB_SE=../boost/lib/libboost_serialization.so
HTSLIB_INC=../htslib
SDSL_INC=../sdsl/include
SDSL_LIB=../sdsl/lib


#COMPILATION RULES
all: $(BFILE)
	
$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ -I${HTSLIB_INC} -L${HTSLIB_INC} -lhts -I${SDSL_INC} -L${SDSL_LIB} -lsdsl -fopenmp $(BOOST_LIB_IO) $(BOOST_LIB_PO) $(BOOST_LIB_SE) $(BGEN_LIB) -o $@ $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I${HTSLIB_INC} -L${HTSLIB_INC} -lhts -I${SDSL_INC} -L${SDSL_LIB} -lsdsl -fopenmp -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE)
