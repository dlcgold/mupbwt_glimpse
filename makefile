PROJECTS = phase split_reference build_ser_from_bin build_ser_from_bcf
UPSTREAM = chunk ligate concordance
HTSLIB = htslib/libhts.a
BOOST_LIBS = \
	boost/lib/libboost_iostreams.so \
	boost/lib/libboost_program_options.so \
	boost/lib/libboost_serialization.so

ifeq ($(DNANEXUS),1)
DEPS =
UPSTREAM_TARGET = system
UPSTREAM_ENV =
else
DEPS = $(HTSLIB) $(BOOST_LIBS)
UPSTREAM_TARGET = docker
UPSTREAM_ENV = \
	HTSLIB_INC=$(CURDIR)/htslib \
	HTSLIB_LIB=$(CURDIR)/htslib/libhts.a \
	BOOST_INC=$(CURDIR)/boost/include \
	BOOST_LIB_IO=$(CURDIR)/boost/lib/libboost_iostreams.a \
	BOOST_LIB_PO=$(CURDIR)/boost/lib/libboost_program_options.a \
	BOOST_LIB_SE=$(CURDIR)/boost/lib/libboost_serialization.a
endif

.PHONY: all clean deps upstream $(PROJECTS) $(addprefix upstream-,$(UPSTREAM))

all: deps $(PROJECTS) upstream

deps: $(DEPS)

$(PROJECTS) $(addprefix upstream-,$(UPSTREAM)): | deps

$(PROJECTS):
	$(MAKE) -C $@ $(COMPILATION_ENV)

upstream: $(addprefix upstream-,$(UPSTREAM))

$(addprefix upstream-,$(UPSTREAM)):
	$(MAKE) -C $(patsubst upstream-%,%,$@) $(UPSTREAM_TARGET) $(UPSTREAM_ENV)

$(HTSLIB):
	wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
	tar -xf htslib-1.16.tar.bz2
	mv htslib-1.16 htslib
	cd htslib && make
	rm -f htslib-1.16.tar.bz2

$(BOOST_LIBS) &:
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
	for dir in $(UPSTREAM); do \
		$(MAKE) clean -C $$dir; \
	done
