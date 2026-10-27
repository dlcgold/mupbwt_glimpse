PROJECTS = phase split_reference build_ser_from_bin build_ser_from_bcf
ORIGINAL = chunk ligate concordance phase split_reference
HTSLIB = htslib/libhts.a
BOOST_LIBS = \
	boost/lib/libboost_iostreams.so \
	boost/lib/libboost_program_options.so \
	boost/lib/libboost_serialization.so

ifeq ($(DNANEXUS),1)
DEPS =
ORIGINAL_TARGET = system
ORIGINAL_ENV =
else
DEPS = $(HTSLIB) $(BOOST_LIBS)
ORIGINAL_TARGET = docker
ORIGINAL_ENV = \
	HTSLIB_INC=$(CURDIR)/htslib \
	HTSLIB_LIB=$(CURDIR)/htslib/libhts.a \
	BOOST_INC=$(CURDIR)/boost/include \
	BOOST_LIB_IO=$(CURDIR)/boost/lib/libboost_iostreams.a \
	BOOST_LIB_PO=$(CURDIR)/boost/lib/libboost_program_options.a \
	BOOST_LIB_SE=$(CURDIR)/boost/lib/libboost_serialization.a
endif

.PHONY: all clean deps original $(PROJECTS) $(addprefix original-,$(ORIGINAL))

all: deps $(PROJECTS) original

deps: $(DEPS)

$(PROJECTS):
	$(MAKE) -C $@ $(COMPILATION_ENV)

original: $(addprefix original-,$(ORIGINAL))

$(addprefix original-,$(ORIGINAL)):
	$(MAKE) -C original/$(patsubst original-%,%,$@) $(ORIGINAL_TARGET) $(ORIGINAL_ENV)

$(HTSLIB):
	wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
	tar -xf htslib-1.16.tar.bz2
	mv htslib-1.16 htslib
	cd htslib && make
	rm -f htslib-1.16.tar.bz2

$(BOOST_LIBS):
	wget https://archives.boost.io/release/1.73.0/source/boost_1_73_0.tar.bz2
	tar -xjf boost_1_73_0.tar.bz2
	cd boost_1_73_0 && \
		./bootstrap.sh \
			--with-libraries=iostreams,program_options,serialization \
			--prefix=../boost && \
		./b2 install
	rm -f boost_1_73_0.tar.bz2
	rm -rf boost_1_73_0

clean:
	for dir in $(PROJECTS); do \
		$(MAKE) clean -C $$dir; \
	done
	for dir in $(ORIGINAL); do \
		$(MAKE) clean -C original/$$dir; \
	done
