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

#include "../PIDX_inc.h"

unsigned int getNumBits ( unsigned int v )
{
  return (unsigned int)floor((log2(v))) + 1;
}

uint64_t getPowerOf2(int x)
{
  /*  find the power of 2 of an integer value (example 5->8) */
  int n = 1;
  while (n < x) n <<= 1;
  return n;
}

unsigned int getLevelFromBlock (uint64_t block, int bits_per_block)
{
  if (block == 0)
    return 0;
  else
    return (unsigned int)floor((log2(getPowerOf2(block)))) + bits_per_block;

  return 0;
}

unsigned int getLeveL (uint64_t index)
{
  if (index)
    return (unsigned int)floor((log2(index))) + 1;
  /* To deal with index 0 - level 0 */
  return 0;
}

int isValidBox(int** box)
{
  int D;
  for (D = 0 ; D < PIDX_MAX_DIMENSIONS ; D++)
  {
    //fprintf(stderr, "VFY: %d %d\n", box[0][D], box[1][D]);
    if (! (box[0][D]>=0 && box[0][D]<=box[1][D]))
      return 0;
  }
  return 1;
}

void Deinterleave(const char* bitmask, int maxh, uint64_t zaddress, int* point)
{
  //Z deinterleave (see papers!)
  int n = 0, bit;
  int* cnt_point = (int*)malloc(PIDX_MAX_DIMENSIONS*sizeof(int));
  memset(cnt_point, 0, PIDX_MAX_DIMENSIONS*sizeof(int));

  for (;zaddress;zaddress >>= 1,++n,maxh--)
  {
    bit=bitmask[maxh];
    point[bit] |= (zaddress & 1) << cnt_point[bit];
    ++cnt_point[bit];
  }
  free(cnt_point);
  cnt_point = 0;
}

uint64_t ZBitmask(const char* bitmask,int maxh)
{
  return ((uint64_t)1)<<maxh;
}

uint64_t ZStart(const char* bitmask,int maxh,int BlockingH)
{
  if (!BlockingH)
    return 0;
  assert(BlockingH>=1 && BlockingH<=maxh);
  return ((uint64_t)1)<<(maxh-BlockingH);
}

uint64_t ZEnd(const char* bitmask,int maxh,int BlockingH)
{
  if (!BlockingH)
    return 0;
  assert(BlockingH>=1 && BlockingH<=maxh);
  return (ZBitmask(bitmask,maxh)-1)-(ZStart(bitmask,maxh,BlockingH)-1);
}

void ZDelta(const char* bitmask, int maxh, int BlockingH, int* point)
{
  int K, bit;
  for (K = 0; K < PIDX_MAX_DIMENSIONS; K++)
    point[K] = 1;

  if (!BlockingH) return;
  assert(BlockingH>=1 && BlockingH<=maxh);
  for (K=maxh;K>=BlockingH;K--)
  {
    bit=bitmask[K];
    point[bit] <<= 1;
  }
}

void GetBoxIntersection(int** inputBox1, int** inputBox2, int** outputBox)
{
  //returns the intersection of two boxes
  int i;
  for (i = 0 ; i < PIDX_MAX_DIMENSIONS ; i++)
  {
    outputBox[0][i]=Max2ab(inputBox1[0][i], inputBox2[0][i]);
    outputBox[1][i]=Min2ab(inputBox1[1][i], inputBox2[1][i]);
  }
}

int** AlignEx(int** box, int* p0, int* delta)
{
  int mod,i;
  if (!isValidBox(box))
    return box;
  for ( i = 0 ; i < PIDX_MAX_DIMENSIONS ; i++)
  {
    mod=(box[0][i]-p0[i]) % delta[i];
    if (mod) box[0][i]+=(delta[i]-mod);
    mod=(box[1][i]-p0[i]) % delta[i];
    if (mod) box[1][i]-=mod;
  }
  return box;
}

void revstr(char* str)
{
  uint64_t i;
  char* cpstr = (char*)malloc(strlen(str)+1);
  //char cpstr[strlen(str)+1];
  for (i=0; i < (int)strlen(str); i++)
    cpstr[i] = str[strlen(str)-i-1];

  cpstr[i] = '\0';
  strcpy(str, cpstr);
  free(cpstr);
}

void GuessBitmaskPattern(char* _bits, Point3D dims)
{
  int D,N,ordered;
  int dim = 1;
  char* p=_bits;

  Point3D id,sorted_id;

  *p++='V';

  For (D)
  {
    PGET(dims,D)=( int)getPowerOf2(PGET(dims,D));
    PGET(id,D)=D;
  }

  //order is ASC order (from smaller dimension to bigger)
  for (ordered=0;!ordered;)
  {
    ordered=1;
    OffsetFor (D,0,-1)
    {
      int ref0=PGET(id,D  ),dim0=PGET(dims,ref0);
      int ref1=PGET(id,D+1),dim1=PGET(dims,ref1);
      if (!(dim0<dim1 || (dim0==dim1 && ref0<ref1)))
      {
    int _temp=PGET(id,D);
    PGET(id,D)=PGET(id,D+1);
    PGET(id,D+1)=_temp;
    ordered=0;
      }
    }
  }

  For (D)
  {
    //order in DESC order
    for (ordered=0,sorted_id=id;!ordered;)
    {
      ordered=1;
      OffsetFor (N,D,-1)
      {
    if (PGET(sorted_id,N+0)<PGET(sorted_id,N+1))
    {
      int _temp=PGET(sorted_id,N);
      PGET(sorted_id,N)=PGET(sorted_id,N+1);
      PGET(sorted_id,N+1)=_temp;
      ordered=0;
    }
      }
    }
    //while dim is not consumed
    for (;dim<PGET(dims,PGET(id,D));dim<<=1)
    {
      OffsetFor (N,D,0)
      {
    *p++='0'+PGET(sorted_id,N);
      }
    }
  }
  *p++=0;
  revstr(_bits+1);
  //strrev(_bits+1)
}

/** Find the first power of, say, 2 that is greater than or equal to the given number. */
int pow_greater_equal(int base, int num)
{
  assert(base > 1);
  assert(num > 0);

  int result = 1;
  while (result < num)
    result *= base;

  return result;
}

/** Guess a bit string given the dimensions. Z has the highest priority, then Y, then X. */
void guess_bit_string(char* bit_string, const Point3D dims)
{
  Point3D power_2_dims;
  power_2_dims.x = pow_greater_equal(2, dims.x);
  power_2_dims.y = pow_greater_equal(2, dims.y);
  power_2_dims.z = pow_greater_equal(2, dims.z);

  uint64_t size = 0;
  char buffer[65];
  while (power_2_dims.x > 1 || power_2_dims.y > 1 || power_2_dims.z > 1)
  {
    int max = Max2ab(power_2_dims.z, Max2ab(power_2_dims.x, power_2_dims.y));
    if (max == power_2_dims.z)
    {
      power_2_dims.z /= 2;
      buffer[size++] = '2';
    }
    else if (max == power_2_dims.y)
    {
      power_2_dims.y /= 2;
      buffer[size++] = '1';
    }
    else
    {
      power_2_dims.x /= 2;
      buffer[size++] = '0';
    }
  }

  //bit_string.size = size;
  assert(size >= 0);
  uint64_t i = 0;
  bit_string[0] = 'V';
  for ( i = 1; i < size + 1; ++i)
  {
    bit_string[i] = buffer[i - 1];
  }
  bit_string[size + 1] = '\0';
}



void guess_bit_string_ZYX(char* bit_string, const Point3D dims)
{
  Point3D power_2_dims;
  power_2_dims.x = pow_greater_equal(2, dims.x);
  power_2_dims.y = pow_greater_equal(2, dims.y);
  power_2_dims.z = pow_greater_equal(2, dims.z);

  uint64_t size = 0;
  char buffer[65];
  while (power_2_dims.x > 1 || power_2_dims.y > 1 || power_2_dims.z > 1)
  {
    int max = Max2ab(power_2_dims.z, Max2ab(power_2_dims.x, power_2_dims.y));
    if (max == power_2_dims.z)
    {
      power_2_dims.z /= 2;
      buffer[size++] = '2';
    }
    else if (max == power_2_dims.y)
    {
      power_2_dims.y /= 2;
      buffer[size++] = '1';
    }
    else
    {
      power_2_dims.x /= 2;
      buffer[size++] = '0';
    }
  }

  //bit_string.size = size;
  assert(size >= 0);
  uint64_t i = 0;
  bit_string[0] = 'V';
  for ( i = 1; i < size + 1; ++i)
  {
    bit_string[i] = buffer[i - 1];
  }
  bit_string[size + 1] = '\0';
}


void guess_bit_string_YXZ(char* bit_string, const Point3D dims)
{
  Point3D power_2_dims;
  power_2_dims.x = pow_greater_equal(2, dims.x);
  power_2_dims.y = pow_greater_equal(2, dims.y);
  power_2_dims.z = pow_greater_equal(2, dims.z);

  uint64_t size = 0;
  char buffer[65];
  while (power_2_dims.x > 1 || power_2_dims.y > 1 || power_2_dims.z > 1)
  {
    int max = Max2ab(power_2_dims.y, Max2ab(power_2_dims.z, power_2_dims.x));
    if (max == power_2_dims.y)
    {
      power_2_dims.y /= 2;
      buffer[size++] = '1';
    }
    else if (max == power_2_dims.x)
    {
      power_2_dims.x /= 2;
      buffer[size++] = '0';
    }
    else
    {
      power_2_dims.z /= 2;
      buffer[size++] = '2';
    }
  }

  //bit_string.size = size;
  assert(size >= 0);
  uint64_t i = 0;
  bit_string[0] = 'V';
  for ( i = 1; i < size + 1; ++i)
  {
    bit_string[i] = buffer[i - 1];
  }
  bit_string[size + 1] = '\0';
}


void guess_bit_string_XZY(char* bit_string, const Point3D dims)
{
  Point3D power_2_dims;
  power_2_dims.x = pow_greater_equal(2, dims.x);
  power_2_dims.y = pow_greater_equal(2, dims.y);
  power_2_dims.z = pow_greater_equal(2, dims.z);

  uint64_t size = 0;
  char buffer[65];
  while (power_2_dims.x > 1 || power_2_dims.y > 1 || power_2_dims.z > 1)
  {
    int max = Max2ab(power_2_dims.x, Max2ab(power_2_dims.y, power_2_dims.z));
    if (max == power_2_dims.x)
    {
      power_2_dims.x /= 2;
      buffer[size++] = '0';
    }
    else if (max == power_2_dims.z)
    {
      power_2_dims.z /= 2;
      buffer[size++] = '2';
    }
    else
    {
      power_2_dims.y /= 2;
      buffer[size++] = '1';
    }
  }

  //bit_string.size = size;
  assert(size >= 0);
  uint64_t i = 0;
  bit_string[0] = 'V';
  for ( i = 1; i < size + 1; ++i)
  {
    bit_string[i] = buffer[i - 1];
  }
  bit_string[size + 1] = '\0';
}



void guess_bit_string_Z(char* bit_string, const Point3D dims)
{
  Point3D power_2_dims;
  power_2_dims.x = pow_greater_equal(2, dims.x);
  power_2_dims.y = pow_greater_equal(2, dims.y);
  power_2_dims.z = pow_greater_equal(2, dims.z);

  uint64_t size = 0;
  char buffer[65];
  while (power_2_dims.z > 1)
  {
    power_2_dims.z /= 2;
    buffer[size++] = '2';
  }

  while (power_2_dims.y > 1)
  {
    power_2_dims.y /= 2;
    buffer[size++] = '1';
  }

  while (power_2_dims.x > 1)
  {
    power_2_dims.x /= 2;
    buffer[size++] = '0';
  }

  //bit_string.size = size;
  assert(size >= 0);
  uint64_t i = 0;
  bit_string[0] = 'V';
  for ( i = 1; i < size + 1; ++i)
  {
    bit_string[i] = buffer[i - 1];
  }
  bit_string[size + 1] = '\0';
}

void guess_bit_string_Y(char* bit_string, const Point3D dims)
{
  Point3D power_2_dims;
  power_2_dims.x = pow_greater_equal(2, dims.x);
  power_2_dims.y = pow_greater_equal(2, dims.y);
  power_2_dims.z = pow_greater_equal(2, dims.z);

  uint64_t size = 0;
  char buffer[65];

  while (power_2_dims.y > 1)
  {
    power_2_dims.y /= 2;
    buffer[size++] = '1';
  }

  while (power_2_dims.z > 1)
  {
    power_2_dims.z /= 2;
    buffer[size++] = '2';
  }

  while (power_2_dims.x > 1)
  {
    power_2_dims.x /= 2;
    buffer[size++] = '0';
  }

  //bit_string.size = size;
  assert(size >= 0);
  uint64_t i = 0;
  bit_string[0] = 'V';
  for ( i = 1; i < size + 1; ++i)
  {
    bit_string[i] = buffer[i - 1];
  }
  bit_string[size + 1] = '\0';
}

void guess_bit_string_X(char* bit_string, const Point3D dims)
{
  Point3D power_2_dims;
  power_2_dims.x = pow_greater_equal(2, dims.x);
  power_2_dims.y = pow_greater_equal(2, dims.y);
  power_2_dims.z = pow_greater_equal(2, dims.z);

  uint64_t size = 0;
  char buffer[65];

  while (power_2_dims.x > 1)
  {
    power_2_dims.x /= 2;
    buffer[size++] = '0';
  }

  while (power_2_dims.y > 1)
  {
    power_2_dims.y /= 2;
    buffer[size++] = '1';
  }

  while (power_2_dims.z > 1)
  {
    power_2_dims.z /= 2;
    buffer[size++] = '2';
  }

  //bit_string.size = size;
  assert(size >= 0);
  uint64_t i = 0;
  bit_string[0] = 'V';
  for ( i = 1; i < size + 1; ++i)
  {
    bit_string[i] = buffer[i - 1];
  }
  bit_string[size + 1] = '\0';
}


static void freeBox(int** box)
{
  free(box[0]);
  box[0] = 0;
  free(box[1]);
  box[1] = 0;

  free(box);
  box = 0;
}

void Align(int maxh, int H, const char* bitmask, int** userBox, int** a_offset, int** a_count, int** nsamples)
{
  int** alignedBox;
  int* h_delta;
  int** h_box;

  if (!isValidBox(userBox))
    return;

  if (!(H>=0 && H<=maxh))
    return;

  h_delta = (int*)malloc(PIDX_MAX_DIMENSIONS *sizeof(int));
  memset(h_delta, 0, PIDX_MAX_DIMENSIONS *sizeof(int));

  ZDelta(bitmask,maxh,H, h_delta);

  //ZBox
  h_box = (int**)malloc(2* sizeof(int*));
  memset(h_box, 0, 2* sizeof(int*));

  h_box[0] = (int*)malloc(PIDX_MAX_DIMENSIONS * sizeof(int));
  h_box[1] = (int*)malloc(PIDX_MAX_DIMENSIONS * sizeof(int));
  memset(h_box[0], 0, PIDX_MAX_DIMENSIONS * sizeof(int));
  memset(h_box[1], 0, PIDX_MAX_DIMENSIONS * sizeof(int));

  if (!H)
  {
    h_box[0][0] = h_box[0][1] = h_box[0][2] = 0;
    h_box[1][0] = h_box[1][1] = h_box[1][2] = 0;
  }
  else
  {
    assert(H>=1 && H<=maxh);
    Deinterleave(bitmask,maxh,ZStart(bitmask,maxh,H), h_box[0]);
    Deinterleave(bitmask,maxh,ZEnd  (bitmask,maxh,H), h_box[1]);
  }

  alignedBox = (int**)malloc(2* sizeof(int*));
  memset(alignedBox, 0, 2* sizeof(int*));
  alignedBox[0] = (int*)malloc(PIDX_MAX_DIMENSIONS * sizeof(int));
  alignedBox[1] = (int*)malloc(PIDX_MAX_DIMENSIONS * sizeof(int));
  memset(alignedBox[0], 0, PIDX_MAX_DIMENSIONS * sizeof(int));
  memset(alignedBox[1], 0, PIDX_MAX_DIMENSIONS * sizeof(int));

  //calculate intersection of the query with current H box
  GetBoxIntersection(userBox, h_box, alignedBox);

  //the box is not valid
  if (!isValidBox(alignedBox))
  {
    freeBox(h_box);
    freeBox(alignedBox);
    free(h_delta);
    return;
  }

  alignedBox = AlignEx(alignedBox, h_box[0], h_delta);

  //invalid box
  if (!isValidBox(alignedBox))
  {
    freeBox(h_box);
    freeBox(alignedBox);
    free(h_delta);
    return;
  }

  int i;
  for ( i = 0 ; i < PIDX_MAX_DIMENSIONS ; i++)
    nsamples[H][i]=1 + (alignedBox[1][i]-alignedBox[0][i])/h_delta[i];

  memcpy(a_offset[H], alignedBox[0], PIDX_MAX_DIMENSIONS * sizeof(int));
  memcpy(a_count[H], alignedBox[1], PIDX_MAX_DIMENSIONS * sizeof(int));

  freeBox(h_box);
  freeBox(alignedBox);
  free(h_delta);
  return;
}

int RegExBitmaskBit(const char* bitmask_pattern,int N)
{
  const char *OpenRegEx;
  int S, L;
  assert(bitmask_pattern[0]=='V');

  if (!N)
    return bitmask_pattern[0];

  if ((OpenRegEx=strchr(bitmask_pattern,'{')))
  {
    S = 1+OpenRegEx-bitmask_pattern;
    L = strchr(bitmask_pattern,'}')-bitmask_pattern-S;

    if ((N+1)<S)
      return bitmask_pattern[N]-'0';
    else
      return bitmask_pattern[S+((N+1-S)%L)]-'0';
  }
  return bitmask_pattern[N]-'0';
}

void Hz_to_xyz(const char* bitmask,  int maxh, uint64_t hzaddress, uint64_t* xyz)
{
  uint64_t lastbitmask=((uint64_t)1)<<maxh;

  hzaddress <<= 1;
  hzaddress  |= 1;
  while ((lastbitmask & hzaddress) == 0) hzaddress <<= 1;
    hzaddress &= lastbitmask - 1;

  Point3D cnt;
  Point3D p  ;
  int n = 0;

  memset(&cnt,0,sizeof(Point3D));
  memset(&p  ,0,sizeof(Point3D));

  for (;hzaddress; hzaddress >>= 1,++n, maxh--)
  {
    int bit= bitmask[maxh];
    PGET(p,bit) |= (hzaddress & 1) << PGET(cnt,bit);
    ++PGET(cnt,bit);
  }
  xyz[0] = p.x;
  xyz[1] = p.y;
  xyz[2] = p.z;
}

uint64_t xyz_to_HZ(const char* bitmask, int maxh, Point3D xyz)
{
  uint64_t zaddress=0;
  int cnt   = 0;
  Point3D zero;
  int temp_maxh = maxh;
  memset(&zero,0,sizeof(Point3D));

  for (cnt=0 ; memcmp(&xyz, &zero, sizeof(Point3D)) ; cnt++, maxh--)
  {
    int bit= bitmask[maxh];
    zaddress |= ((uint64_t)PGET(xyz,bit) & 1) << cnt;
    PGET(xyz,bit) >>= 1;
  }

  uint64_t lastbitmask=((uint64_t)1)<<temp_maxh;
  zaddress |= lastbitmask;
  while (!(1 & zaddress)) zaddress >>= 1;
    zaddress >>= 1;

  return zaddress;
}

int VisusSplitFilename(const char* filename,char* dirname,char* basename)
{
  int i;
  int N=strlen(filename);

  if (!N)
    return 0;

  //find the last separator
  for (i=N-1;i>=0;i--)
  {
    if (filename[i]=='/' || filename[i]=='\\')
    {
      strncpy(dirname,filename,i);
      dirname[i]=0;
      strcpy(basename,filename+i+1);
      return 1;
    }
  }
  //assume is all filename (without directory name)
  dirname [0]=0;
  strcpy(basename,filename);
  return 1;
}


/// Returns elapsed time
double PIDX_get_time()
{
  return MPI_Wtime();

//  struct timeval temp;
//  gettimeofday(&temp, NULL);
//  return (double)(temp.tv_sec) + (double)(temp.tv_usec)/1000000.0;

}


#undef max
#define max(a,b) ((a) > (b) ? (a) : (b))
Point3D get_strides(const char* bit_string, int bs_len, int len)
{
  uint64_t i = 0;
  assert(len >= 0);
  Point3D stride = { 0, 0, 0 };
  uint64_t start = max(bs_len - len, 0);
  for (i = start; i < bs_len; ++i)
  {
    if      (bit_string[i] == '0') { ++stride.x; }
    else if (bit_string[i] == '1') { ++stride.y; }
    else if (bit_string[i] == '2') { ++stride.z; }
  }
  if (len > bs_len) { ++stride.x; ++stride.y; ++stride.z; }
  stride.x = 1 << stride.x;
  stride.y = 1 << stride.y;
  stride.z = 1 << stride.z;
  return stride;
}
#undef max


Point3D get_intra_block_strides(const char* bit_string, int bs_len, int hz_level)
{
  // count the number of x, y, z in the least significant (z_level + 1) bits
  // in the bit_string
  int z_level = bs_len - hz_level;
  int len = z_level + 1;
  return get_strides(bit_string, bs_len, len);
}


/* Return the strides (in terms of the first sample) of idx blocks, in x, y, and z. */
Point3D get_inter_block_strides(const char* bit_string, int bs_len, int hz_level, int bits_per_block)
{
  assert(bs_len >= hz_level);
  // count the number of x, y, z in the least significant
  // (z_level + bits_per_block + 1) bits in the bit_string
  int len = bs_len - hz_level + bits_per_block + 1;
  // len can get bigger than bit_string.size if the input hz_level is smaller
  // than the mininum hz level
  return get_strides(bit_string, bs_len, len);
}


/* Get the number of samples in each dimension for a block at the given hz level */
Point3D get_num_samples_per_block(const char* bit_string, int bs_len, int hz_level, int bits_per_block)
{
  Point3D intra_stride = get_intra_block_strides(bit_string, bs_len, hz_level);
  Point3D inter_stride = get_inter_block_strides(bit_string, bs_len, hz_level, bits_per_block);
  Point3D block_nsamples;
  block_nsamples.x = inter_stride.x / intra_stride.x;
  block_nsamples.y = inter_stride.y / intra_stride.y;
  block_nsamples.z = inter_stride.z / intra_stride.z;
  return block_nsamples;
}


void get_grid( Point3D sub_vol_from, Point3D sub_vol_to, int hz_level, const char* bit_string, int bs_len, Point3D* from, Point3D* to, Point3D* stride)
{
  *stride = get_intra_block_strides(bit_string, bs_len, hz_level);
  Point3D start = get_first_coord(bit_string, bs_len, hz_level);
  Point3D end = get_last_coord(bit_string, bs_len, hz_level);

  intersect_grid(sub_vol_from, sub_vol_to, start, end, stride, from, to);

  return;
}


Point3D get_first_coord(const char* bit_string, int bs_len, int hz_level)
{
  if (hz_level == 0)
  {
    Point3D zero = {0, 0, 0};
    return zero;
  }

  uint64_t i = 0;
  int pos = hz_level - 1; // the position of the right-most 1 bit in the bit string
  // count the number of "bits" that is the same with the one at position pos
  int count = 0;
  char c = bit_string[pos];
  for (i = pos + 1; i < bs_len; ++i)
  {
    if (bit_string[i] == c)
      ++count;
  }

  // raise the corresponding coordinate to the appropriate power of 2 (the other
  // 2 coordinates are 0)
  Point3D coord = {0, 0, 0};
  if (c == '0')
    coord.x = (int)pow(2, count);

  else if (c == '1')
    coord.y = (int)pow(2, count);

  else if (c == '2')
    coord.z = (int)pow(2, count);

  return coord;
}


Point3D get_last_coord(const char* bit_string, int bs_len, int hz_level)
{

  int i = 0;
  if (hz_level == 0)
  {
    Point3D zero = {0, 0, 0};
    return zero;
  }


  int pos = hz_level - 1; // the position of the right-most 1 bit in the bit string
  int size = bs_len;
  Point3D count = {0, 0, 0};
  for (i = size - 1; i > pos; --i)
  {
    if (bit_string[i] == '0')
      ++count.x;

    else if (bit_string[i] == '1')
      ++count.y;

    else if (bit_string[i] == '2')
      ++count.z;

  }

  Point3D coord = {0, 0, 0};
  for (i = pos; i >= 0; --i)
  {
    if (bit_string[i] == '0')
      coord.x += (int)pow(2, count.x++);

    else if (bit_string[i] == '1')
      coord.y += (int)pow(2,count.y++);

    else if (bit_string[i] == '2')
      coord.z += (int)pow(2, count.z++);

  }
  return coord;
}


#undef min
#define min(a,b) ((a) < (b) ? (a) : (b))
void intersect_grid(Point3D vol_from, Point3D vol_to, Point3D from, Point3D to, Point3D* stride, Point3D* output_from, Point3D* output_to)
{
    Point3D min_to = vol_to;
    min_to.x = min(min_to.x, to.x);
    min_to.y = min(min_to.y, to.y);
    min_to.z = min(min_to.z, to.z);

    (*output_from).x = from.x + ((vol_from.x + (*stride).x - 1 - from.x) / (*stride).x) * (*stride).x;
    (*output_to).x = from.x + ((min_to.x - from.x) / (*stride).x) * (*stride).x;

    (*output_from).y = from.y + ((vol_from.y + (*stride).y - 1 - from.y) / (*stride).y) * (*stride).y;
    (*output_to).y = from.y + ((min_to.y - from.y) / (*stride).y) * (*stride).y;

    (*output_from).z = from.z + ((vol_from.z + (*stride).z - 1 - from.z) / (*stride).z) * (*stride).z;
    (*output_to).z = from.z + ((min_to.z - from.z) / (*stride).z) * (*stride).z;

    // we need to do the following corrections because the behavior of integer
    // division with negative integers are not well defined...
    if (vol_from.x < from.x) { (*output_from).x = from.x; }
    if (vol_from.y < from.y) { (*output_from).y = from.y; }
    if (vol_from.z < from.z) { (*output_from).z = from.z; }
    if (min_to.x < from.x) { (*output_to).x = from.x - (*stride).x; }
    if (min_to.y < from.y) { (*output_to).y = from.y - (*stride).y; }
    if (min_to.z < from.z) { (*output_to).z = from.z - (*stride).z; }

    return;// output_from <= output_to;
}
#undef min


PIDX_return_code PIDX_get_datatype_details(PIDX_data_type type, int* values, int* bits)
{
    if (strcmp(type, INT8) == 0)
    {
      *bits = 8;
      *values = 1;
    }
    else if (strcmp(type, INT8_GA) == 0)
    {
      *bits = 8;
      *values = 2;
    }
    else if (strcmp(type, INT8_RGB) == 0)
    {
      *bits = 8;
      *values = 3;
    }
    else if (strcmp(type, INT8_RGBA) == 0)
    {
      *bits = 8;
      *values = 4;
    }

    else if (strcmp(type, UINT8) == 0)
    {
      *bits = 8;
      *values = 1;
    }
    else if (strcmp(type, UINT8_GA) == 0)
    {
      *bits = 8;
      *values = 2;
    }
    else if (strcmp(type, UINT8_RGB) == 0)
    {
      *bits = 8;
      *values = 3;
    }
    else if (strcmp(type, UINT8_RGBA) == 0)
    {
      *bits = 8;
      *values = 4;
    }

    else if (strcmp(type, INT16) == 0)
    {
      *bits = 16;
      *values = 1;
    }
    else if (strcmp(type, INT16_GA) == 0)
    {
      *bits = 16;
      *values = 2;
    }
    else if (strcmp(type, INT16_RGB) == 0)
    {
      *bits = 16;
      *values = 3;
    }
    else if (strcmp(type, INT16_RGBA) == 0)
    {
      *bits = 16;
      *values = 4;
    }

    else if (strcmp(type, UINT16) == 0)
    {
      *bits = 16;
      *values = 1;
    }
    else if (strcmp(type, UINT16_GA) == 0)
    {
      *bits = 16;
      *values = 2;
    }
    else if (strcmp(type, UINT16_RGB) == 0)
    {
      *bits = 16;
      *values = 3;
    }
    else if (strcmp(type, UINT16_RGBA) == 0)
    {
      *bits = 16;
      *values = 4;
    }

    else if (strcmp(type, INT32) == 0)
    {
      *bits = 32;
      *values = 1;
    }
    else if (strcmp(type, INT32_GA) == 0)
    {
      *bits = 32;
      *values = 2;
    }
    else if (strcmp(type, INT32_RGB) == 0)
    {
      *bits = 32;
      *values = 3;
    }
    else if (strcmp(type, INT32_RGBA) == 0)
    {
      *bits = 32;
      *values = 4;
    }

    else if (strcmp(type, UINT32) == 0)
    {
      *bits = 32;
      *values = 1;
    }
    else if (strcmp(type, UINT32_GA) == 0)
    {
      *bits = 32;
      *values = 2;
    }
    else if (strcmp(type, UINT32_RGB) == 0)
    {
      *bits = 32;
      *values = 3;
    }
    else if (strcmp(type, UINT32_RGBA) == 0)
    {
      *bits = 32;
      *values = 4;
    }

    else if (strcmp(type, INT64) == 0)
    {
      *bits = 64;
      *values = 1;
    }
    else if (strcmp(type, INT64_GA) == 0)
    {
      *bits = 64;
      *values = 2;
    }
    else if (strcmp(type, INT64_RGB) == 0)
    {
      *bits = 64;
      *values = 3;
    }
    else if (strcmp(type, INT64_RGBA) == 0)
    {
      *bits = 64;
      *values = 4;
    }

    else if (strcmp(type, UINT64) == 0)
    {
      *bits = 64;
      *values = 1;
    }
    else if (strcmp(type, UINT64_GA) == 0)
    {
      *bits = 64;
      *values = 2;
    }
    else if (strcmp(type, UINT64_RGB) == 0)
    {
      *values = 3;
      *bits = 64;
    }
    else if (strcmp(type, UINT64_RGBA) == 0)
    {
      *values = 4;
      *bits = 64;
    }

    else if (strcmp(type, FLOAT32) == 0)
    {
      *values = 1;
      *bits = 32;
    }
    else if (strcmp(type, FLOAT32_GA) == 0)
    {
      *values = 2;
      *bits = 32;
    }
    else if (strcmp(type, FLOAT32_RGB) == 0)
    {
      *values = 3;
      *bits = 32;
    }
    else if (strcmp(type, FLOAT32_RGBA) == 0)
    {
      *values = 4;
      *bits = 32;
    }
    else if (strcmp(type, FLOAT32_7STENCIL) == 0)
    {
      *values = 7;
      *bits = 32;
    }
    else if (strcmp(type, FLOAT32_9TENSOR) == 0)
    {
      *values = 9;
      *bits = 32;
    }

    else if (strcmp(type, FLOAT64) == 0)
    {
      *values = 1;
      *bits = 64;
    }
    else if (strcmp(type, FLOAT64_GA) == 0)
    {
      *values = 2;
      *bits = 64;
    }
    else if (strcmp(type, FLOAT64_RGB) == 0)
    {
      *values = 3;
      *bits = 64;
    }
    else if (strcmp(type, FLOAT64_RGBA) == 0)
    {
      *values = 4;
      *bits = 64;
    }
    else if (strcmp(type, FLOAT64_7STENCIL) == 0)
    {
      *values = 7;
      *bits = 64;
    }
    else if (strcmp(type, FLOAT64_9TENSOR) == 0)
    {
      *values = 9;
      *bits = 64;
    }
    else
      *values = 0;

    return PIDX_success;
}
