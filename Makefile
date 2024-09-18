#
SHELL=/bin/sh
MAKE=make
MAKEFILE=Makefile

LIBDIRS = common-v4.5 \
          numeric-util-v1.0 \
          tables-v8.0 \
          regex-v2.2 \
          cif-file-v1.0 \
          cifparse-obj-v7.0 \
          cif-file-util-v1.0 \
          utillib-v1.1 \
          vflib-2.0.6 \
          connect-v3.3 \
          filterlib-v10.1 \
          rolodist-util \
          validation-v10.1 \
          annotation-v1.0 \
          maxit-v10.1

all:	compile

binary:	compile
	@sh -c './binary.sh'
	@./bin/DictToSdb -ddlFile ./data/ascii/mmcif_ddl.dic -dictFile ./data/ascii/mmcif_pdbx.dic -dictSdbFile mmcif_pdbx.sdb
	@mv mmcif_pdbx.sdb ./data/binary
	@rm -f ./bin/DictToSdb ./bin/cif2bin ./bin/connect_main
	@sh -c 'if [ -e ./mmcif_pdbx.dic-parser.log ]; then rm -rf ./mmcif_pdbx.dic-parser.log; fi'

compile:
	@sh -c 'cd ./etc; ./platform.sh'
	@for libdir in $(LIBDIRS); do \
		echo " "; \
		echo "------------------------------------------------------------"; \
		echo "**** Making $$libdir ****"; \
		echo "------------------------------------------------------------"; \
		(cd $$libdir && $(MAKE) -f $(MAKEFILE) "OPT=-O" install) || exit 1; \
	done
#
debug:
	@sh -c 'cd ./etc; ./platform.sh'
	@for libdir in $(LIBDIRS); do \
		echo " "; \
		echo "------------------------------------------------------------"; \
		echo "**** Making (debug) $$libdir ****"; \
		echo "------------------------------------------------------------"; \
		(cd $$libdir && $(MAKE) -f $(MAKEFILE)  "OPT=-g" install) || exit 1; \
	done
#
clean:  
	@for libdir in $(LIBDIRS); do \
		echo cleaning $$libdir; \
		(cd $$libdir && $(MAKE) -f $(MAKEFILE) clean) || exit 1; \
	done
#
