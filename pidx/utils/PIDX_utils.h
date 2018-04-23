/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2018 ViSUS L.L.C., 
 * Scientific Computing and Imaging Institute of the University of Utah
 * 
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *  
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 * 
 */
 
#ifndef __PIDX_UTILS_H
#define __PIDX_UTILS_H


#define Min2ab(a,b)      (((a)<=(b))?(a):(b))
#define Max2ab(a,b)      (((a)> (b))?(a):(b))
#define OffsetFor(_D_,_From_,_Off_) for((_D_)=(_From_);(_D_)<(PIDX_MAX_DIMENSIONS+(_Off_));(_D_)++)
#define For(_D_) for((_D_)=0;(_D_)<PIDX_MAX_DIMENSIONS;(_D_)++)
#define PGET(_Point_,_Coordinate_) ((&((_Point_).x))[(_Coordinate_)])

typedef struct {int x,y,z;} Point3D;

unsigned int getNumBits ( unsigned int v );

uint64_t getPowerOf2(int x);

unsigned int getLevelFromBlock (uint64_t block, int bits_per_block);

unsigned int getLeveL (uint64_t index);

int isValidBox(int** box);

void Deinterleave(const char* bitmask, int maxh, uint64_t zaddress, int* point);

uint64_t ZBitmask(const char* bitmask,int maxh);

uint64_t ZStart(const char* bitmask,int maxh,int BlockingH);

uint64_t ZEnd(const char* bitmask,int maxh,int BlockingH);

void ZDelta(const char* bitmask, int maxh, int BlockingH, int* point);

void GetBoxIntersection(int** a, int** b, int** c);

int** AlignEx(int** box, int* p0, int* delta);

void revstr(char* str);

void GuessBitmaskPattern(char* _bits, Point3D dims);

void Align(int maxh, int H, const char* bitmask, int** userBox, int** a_offset, int** a_count, int** nsamples);

int RegExBitmaskBit(const char* bitmask_pattern,int N);

uint64_t xyz_to_HZ(const char* bitmask, int maxh, Point3D xyz);

void Hz_to_xyz(const char* bitmask,  int maxh, uint64_t hzaddress, uint64_t* xyz);

int VisusSplitFilename(const char* filename,char* dirname,char* basename);

void guess_bit_string(char* bit_string, const Point3D dims);

void guess_bit_string_X(char* bit_string, const Point3D dims);

void guess_bit_string_Y(char* bit_string, const Point3D dims);

void guess_bit_string_Z(char* bit_string, const Point3D dims);

void guess_bit_string_ZYX(char* bit_string, const Point3D dims);

void guess_bit_string_YXZ(char* bit_string, const Point3D dims);

void guess_bit_string_XZY(char* bit_string, const Point3D dims);

int pow_greater_equal(int base, int num);

double PIDX_get_time();

Point3D get_num_samples_per_block(const char* bit_string, int bs_len, int hz_level, int bits_per_block);

Point3D get_inter_block_strides(const char* bit_string, int bs_len, int hz_level, int bits_per_block);

Point3D get_intra_block_strides(const char* bit_string, int bs_len, int hz_level);

Point3D get_strides(const char* bit_string, int bs_len, int len);

void get_grid( Point3D sub_vol_from, Point3D sub_vol_to, int hz_level, const char* bit_string, int bs_len, Point3D* from, Point3D* to, Point3D* stride);

Point3D get_first_coord(const char* bit_string, int bs_len, int hz_level);

Point3D get_last_coord(const char* bit_string, int bs_len, int hz_level);

void intersect_grid(Point3D vol_from, Point3D vol_to, Point3D from, Point3D to, Point3D* stride, Point3D* output_from, Point3D* output_to);

PIDX_return_code PIDX_get_datatype_details(PIDX_data_type type, int* values, int* bits);
#endif
