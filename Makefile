##################################
#ViSUS Visualization Project                    
# Copyright (c) 2010 University of Utah          
# Scientific Computing and Imaging Institute     
# 72 S Central Campus Drive, Room 3750           
# Salt Lake City, UT 84112                       
#                                                
#For information about this project see:        
#http://www.pascucci.org/visus/                 
#                                                
#      or contact: pascucci@sci.utah.edu         
###################################

PIDX_NEW_OBJS = PIDX.o
PIDX_RST_OBJS = PIDX_rst.o
PIDX_BLOCK_RST_OBJS = PIDX_block_restructure.o
PIDX_COMPRESSION_OBJS = PIDX_compression.o
PIDX_HZ_ENCODE_OBJS = PIDX_hz_encode.o
PIDX_AGG_OBJS = PIDX_agg.o
PIDX_IO_OBJS = PIDX_io.o
PIDX_DATA_STRUCTS = ../pidx/PIDX_memory_layout_data_structs.h ../pidx/PIDX_idx_data_structs.h

MPICC = mpicc -Wall 
MPICXX = mpic++ -Wall 
MPI_LDFLAGS_PIDX = -L. -lPIDX -lm -lzfp -lfpzip -I ../pidx -I../zfp-0.2.1/inc -I../fpzip-1.1.0/inc
MPI_CFLAGS = -Wno-write-strings -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -lm -lzfp -lfpzip -I ../pidx -I../zfp-0.2.1/inc -I../fpzip-1.1.0/inc

FPZIP_FP = FPZIP_FP_FAST
FPZIP_BLOCK_SIZE = 0x1000
FPZIP_CONV = -DWITH_UNION
CXX = mpic++
CXXFLAGS += -ansi -Wall -g
CXXFLAGS += -O3
DEFS += -DFPZIP_BLOCK_SIZE=$(FPZIP_BLOCK_SIZE) -DFPZIP_FP=$(FPZIP_FP) $(FPZIP_CONV)


ZFPCXX = mpic++
ZFPCXXFLAGS = -g -O3 -ansi -Wall -I../zfp-0.2.1/inc $(ZFPDEFS)

all: error.o rcdecoder.o rcencoder.o rcqsmodel.o read.o write.o libfpzip.a libzfp.o libzfp.a PIDX_error_codes.o PIDX_file_access_modes.o PIDX_data_types.o PIDX_data_layout.o PIDX_comm.o PIDX_blocks.o PIDX_utils.o PIDX_point.o PIDX_file_name.o PIDX_header_io.o PIDX_rst.o PIDX_hz_encode.o PIDX_block_restructure.o PIDX_compression.o PIDX_io.o PIDX_agg.o PIDX.o libPIDX.a pidxtest idx-verify

error.o: ../fpzip-1.1.0/src/error.cpp
	$(CXX) $(CXXFLAGS) $(DEFS) -I../fpzip-1.1.0/inc -c $<

rcdecoder.o: ../fpzip-1.1.0/src/rcdecoder.cpp
	$(CXX) $(CXXFLAGS) $(DEFS) -I../fpzip-1.1.0/inc -c $<

rcencoder.o: ../fpzip-1.1.0/src/rcencoder.cpp
	$(CXX) $(CXXFLAGS) $(DEFS) -I../fpzip-1.1.0/inc -c $<
	
rcqsmodel.o: ../fpzip-1.1.0/src/rcqsmodel.cpp
	$(CXX) $(CXXFLAGS) $(DEFS) -I../fpzip-1.1.0/inc -c $<
	
read.o: ../fpzip-1.1.0/src/read.cpp
	$(CXX) $(CXXFLAGS) $(DEFS) -I../fpzip-1.1.0/inc -c $<
	
write.o: ../fpzip-1.1.0/src/write.cpp
	$(CXX) $(CXXFLAGS) $(DEFS) -I../fpzip-1.1.0/inc -c $<

libfpzip.a: error.o rcdecoder.o rcencoder.o rcqsmodel.o read.o write.o
	@rm -f $@
	ar rc $@ $^

libzfp.o: ../zfp-0.2.1/src/libzfp.cpp
	$(ZFPCXX) $(ZFPCXXFLAGS) -c $<

libzfp.a: libzfp.o
	rm -f $@
	ar rc $@ $^
	
PIDX_error_codes.o: ../pidx/PIDX_error_codes.c ../pidx/PIDX_error_codes.h ../pidx/PIDX_inc.h  ../pidx/PIDX_memory_layout_data_structs.h ../pidx/PIDX_idx_data_structs.h
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@

PIDX_file_access_modes.o: ../pidx/PIDX_file_access_modes.c  ../pidx/PIDX_file_access_modes.h ../pidx/PIDX_inc.h  ../pidx/PIDX_memory_layout_data_structs.h ../pidx/PIDX_idx_data_structs.h
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@
	
PIDX_data_types.o: ../pidx/PIDX_data_types.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@

PIDX_data_layout.o: ../pidx/PIDX_data_layout.c ../pidx/PIDX_data_layout.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@
	
PIDX_comm.o: ../pidx/PIDX_comm.c ../pidx/PIDX_comm.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@
	
PIDX_blocks.o: ../pidx/PIDX_blocks.c ../pidx/PIDX_blocks.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@	
	
PIDX_utils.o: ../pidx/PIDX_utils.c ../pidx/PIDX_utils.h $(PIDX_DATA_STRUCTS) ../pidx/PIDX_inc.h
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@

PIDX_point.o: ../pidx/PIDX_point.c ../pidx/PIDX_point.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@
	
PIDX_file_name.o: ../pidx/PIDX_file_name.c ../pidx/PIDX_file_name.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@

PIDX_header_io.o: ../pidx/PIDX_header_io.c ../pidx/PIDX_header_io.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@
	
$(PIDX_RST_OBJS): %.o: ../pidx/PIDX_rst.c ../pidx/PIDX_rst.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c  -o $@

$(PIDX_BLOCK_RST_OBJS): %.o: ../pidx/PIDX_block_restructure.c ../pidx/PIDX_block_restructure.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c  -o $@

$(PIDX_HZ_ENCODE_OBJS): %.o: ../pidx/PIDX_hz_encode.c ../pidx/PIDX_hz_encode.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c  -o $@

$(PIDX_COMPRESSION_OBJS): %.o: ../pidx/PIDX_compression.c ../pidx/PIDX_compression.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c  -o $@

$(PIDX_AGG_OBJS): %.o: ../pidx/PIDX_agg.c ../pidx/PIDX_agg.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c  -o $@

$(PIDX_IO_OBJS): %.o: ../pidx/PIDX_io.c ../pidx/PIDX_io.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c  -o $@

$(PIDX_NEW_OBJS): %.o: ../pidx/PIDX.c ../pidx/PIDX.h ../pidx/PIDX_io.c ../pidx/PIDX_io.h ../pidx/PIDX_agg.c ../pidx/PIDX_agg.h ../pidx/PIDX_compression.c ../pidx/PIDX_compression.h ../pidx/PIDX_hz_encode.c ../pidx/PIDX_hz_encode.h ../pidx/PIDX_block_restructure.c ../pidx/PIDX_block_restructure.h ../pidx/PIDX_rst.c ../pidx/PIDX_rst.h ../pidx/PIDX_header_io.c ../pidx/PIDX_header_io.h ../pidx/PIDX_inc.h  $(PIDX_DATA_STRUCTS)
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@

libPIDX.a: $(PIDX_NEW_OBJS) ../pidx/PIDX_io.c ../pidx/PIDX_io.h ../pidx/PIDX_agg.c ../pidx/PIDX_agg.h ../pidx/PIDX_compression.c ../pidx/PIDX_compression.h ../pidx/PIDX_hz_encode.c ../pidx/PIDX_hz_encode.h ../pidx/PIDX_block_restructure.c ../pidx/PIDX_block_restructure.h ../pidx/PIDX_rst.c ../pidx/PIDX_rst.h ../pidx/PIDX_header_io.c ../pidx/PIDX_header_io.h
	ar rcs $@ $(PIDX_NEW_OBJS) PIDX_error_codes.o PIDX_file_access_modes.o PIDX_data_types.o PIDX_data_layout.o PIDX_comm.o PIDX_blocks.o PIDX_utils.o PIDX_point.o PIDX_file_name.o PIDX_header_io.o PIDX_rst.o PIDX_hz_encode.o PIDX_io.o PIDX_compression.o PIDX_block_restructure.o PIDX_agg.o PIDX.o

pidxtest:  ../test/pidxtest.c ../test/testdefs.h ../test/pidxtest.h ../test/test-PIDX-multi-patch.c ../test/test-PIDX-multi-var-writer.c ../test/test-PIDX-one-var-writer.c ../test/test-PIDX-reader.c ../test/test-PIDX-writer.c ../test/test-serial-writer.c PIDX_data_layout.o  PIDX_error_codes.o PIDX_file_access_modes.o PIDX_data_types.o PIDX_data_layout.o PIDX_comm.o PIDX_blocks.o PIDX_utils.o PIDX_point.o PIDX_file_name.o PIDX_header_io.o PIDX_rst.o PIDX_hz_encode.o PIDX_io.o PIDX_agg.o PIDX.o libPIDX.a
	$(MPICXX) $(MPI_CFLAGS) $< ../test/test-PIDX-multi-idx-writer.c -o $@ $(MPI_LDFLAGS_PIDX)
	rm *.o

idx-verify: ../profile/idx-verify.c
	cc   ../profile/idx-verify.c   -o idx-verify -lm

clean::
	rm -f *.o libPIDX* pidxtest idx-verify