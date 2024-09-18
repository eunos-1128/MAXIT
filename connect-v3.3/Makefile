#
#        connectlib module makefile
#
#----------------------------------------------------------------------------
# Project specific path defintions.
#----------------------------------------------------------------------------
M_INCL_DIR  = ../include
M_LIB_DIR   = ../lib
M_BIN_DIR   = ../bin

PROJ_DIR    = .

L_INCL_DIR  = $(PROJ_DIR)/include
SRC_DIR     = $(PROJ_DIR)/src
OBJ_DIR     = $(PROJ_DIR)/obj
L_LIB_DIR   = $(PROJ_DIR)/lib
L_BIN_DIR   = $(PROJ_DIR)/bin

VPATH = $(OBJ_DIR) $(L_BIN_DIR)

#----------------------------------------------------------------------------
# LINCLUDES and LDEFINES are appended to CFLAGS and C++FLAGS
#----------------------------------------------------------------------------
LDEFINES  =
LINCLUDES = -I$(L_INCL_DIR) -I$(M_INCL_DIR) 

#----------------------------------------------------------------------------
# Include the appropriate compiler/platform definitions ...
#----------------------------------------------------------------------------
include ../etc/Makefile.platform

# Dependent libraries
UTILLIB_LIB       = $(M_LIB_DIR)/utillib.a
CIF_FILE_UTIL_LIB = $(M_LIB_DIR)/cif-file-util.a
CIF_FILE_LIB      = $(M_LIB_DIR)/cif-file.a
CIFPARSE_LIB      = $(M_LIB_DIR)/cifparse-obj.a
TABLES_LIB        = $(M_LIB_DIR)/tables.a
COMMON_LIB        = $(M_LIB_DIR)/common.a
REGEX_LIB         = $(M_LIB_DIR)/regex.a

ALL_DEP_LIBS = $(UTILLIB_LIB) $(CIF_FILE_UTIL_LIB) $(CIF_FILE_LIB) $(CIFPARSE_LIB) \
               $(TABLES_LIB) $(COMMON_LIB) $(REGEX_LIB)

# Module libraries
MOD_LIB = connectlib.a

# Agregate library
AGR_LIB = all.a

# Temporary library. Used to obtain the agregate library.
TMP_LIB = tmp.a

L_MOD_LIB = $(L_LIB_DIR)/$(MOD_LIB)
M_MOD_LIB = $(M_LIB_DIR)/$(MOD_LIB)
L_AGR_LIB = $(L_LIB_DIR)/$(AGR_LIB)
M_AGR_LIB = $(M_LIB_DIR)/$(AGR_LIB)


# Base file names. Must have ".ext" at the end of the file.

BASE_FILES = AtomFormat.ext \
             ConnectDic.ext \
             ConnectFormat.ext

MAIN_FILES = connect_main.ext

# Source files.
SRC_FILES = ${MAIN_FILES:.ext=.C} ${BASE_FILES:.ext=.C}

# Object files. Replace ".ext" with ".o"
OBJ_BASE_FILES = ${BASE_FILES:.ext=.o}

OBJ_FILES = $(OBJ_BASE_FILES) ${MAIN_FILES:.ext=.o}

HEADER_FILES = ${BASE_FILES:.ext=.h}

# Executables. Remove ".ext"
TARGETS = ${MAIN_FILES:.ext=}

# Scripts
TARGET_SCRIPTS = binary.csh


ALL_OBJ_FILES = *.o

.PHONY: ../etc/Makefile.platform all install export clean clean_build

all: install


install: $(M_MOD_LIB) $(TARGETS)

clean: clean_build

# Rule for making executables
%: $(OBJ_DIR)/%.o $(M_MOD_LIB) $(ALL_DEP_LIBS)
	$(CCC) $(LDFLAGS) $< $(M_MOD_LIB) $(ALL_DEP_LIBS) $(MALLOCLIB) -lm -o $(L_BIN_DIR)/$@
	@cp -f $(L_BIN_DIR)/$@ $(M_BIN_DIR)/$@


$(M_MOD_LIB): $(L_MOD_LIB)
#       Install header files
	@cd $(L_INCL_DIR); \
          ../$(INSTALL) $(INSTALLOPTS) $(HEADER_FILES) ../$(M_INCL_DIR)

#       Install module library
	$(INSTALL) $(INSTALLOPTS) $(L_MOD_LIB) $(M_LIB_DIR)

#       Create agregate library

	@cd $(M_LIB_DIR); ../etc/initlib.sh $(MOD_LIB)

	@cd $(L_LIB_DIR); cp ../$(M_AGR_LIB) $(TMP_LIB)
	@cd $(L_LIB_DIR); $(AR) $(AR_GETFLAGS) $(TMP_LIB)
	@cd $(L_LIB_DIR); rm -f $(TMP_LIB)

	@cd $(L_LIB_DIR); cp $(MOD_LIB) $(TMP_LIB)
	@cd $(L_LIB_DIR); $(AR) $(AR_GETFLAGS) $(TMP_LIB)
	@cd $(L_LIB_DIR); rm -f $(TMP_LIB)

	@cd $(L_LIB_DIR); $(AR) $(AR_PUTFLAGS) $(AGR_LIB) $(ALL_OBJ_FILES)
	@cd $(L_LIB_DIR); rm -f $(ALL_OBJ_FILES)

	$(INSTALL) $(INSTALLOPTS) $(L_AGR_LIB) $(M_LIB_DIR)
	@rm -f $(L_AGR_LIB)


clean_build:
	@cd $(M_INCL_DIR); rm -f $(HEADER_FILES)
	@rm -f $(OBJ_DIR)/*.o
	@rm -rf $(OBJ_DIR)/ii_files
	@rm -f $(L_MOD_LIB)
	@rm -f $(M_MOD_LIB)
	@rm -f $(M_AGR_LIB)
	@cd $(L_BIN_DIR); rm -f $(TARGETS)
	@cd $(M_BIN_DIR); rm -f $(TARGETS)


$(L_MOD_LIB): $(OBJ_BASE_FILES)
#       Create module library
	@cd $(OBJ_DIR); $(AR) $(AR_PUTFLAGS) ../$@ $(OBJ_BASE_FILES)
	$(RANLIB) $@
	@echo $@ " is up to date."

# Rule for making object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.C
	$(CCC) $(C++FLAGS) -c $< -o $@

%.o: $(SRC_DIR)/%.C
	$(CCC) $(C++FLAGS) -c $< -o $(OBJ_DIR)/$@

