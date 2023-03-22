all:
	 @mkdir -p $(ARCH)
	 cd $(ARCH) && make SRC=.. -f ../Makefile.build spcm.x

clean:
	 cd $(ARCH) && make SRC=.. -f ../Makefile.build clean
