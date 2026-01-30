PROJECTS = chunk concordance ligate phase split_reference

HTSLIB = htslib/libhts.a
BOOST_LIBS = \
	boost/lib/libboost_iostreams.so \
	boost/lib/libboost_program_options.so \
	boost/lib/libboost_serialization.so
SDSL = sdsl/lib/libsdsl.a


ifeq ($(DNANEXUS),1)
DEPS = $(SDSL)
else
DEPS = $(HTSLIB) $(BOOST_LIBS) $(SDSL)
endif



.PHONY: all clean deps $(PROJECTS)

all: deps $(PROJECTS)

deps: $(DEPS)

$(PROJECTS):
	$(MAKE) -C $@ $(COMPILATION_ENV)


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


$(SDSL):
	git clone https://github.com/simongog/sdsl-lite.git
	cd sdsl-lite && ./install.sh ../sdsl
	rm -rf sdsl-lite


clean:
	for dir in $(PROJECTS); do \
		$(MAKE) clean -C $$dir; \
	done

