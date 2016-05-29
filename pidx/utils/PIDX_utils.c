/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       **
 **                                                 **
 **  PIDX is licensed under the Creative Commons    **
 **  Attribution-NonCommercial-NoDerivatives 4.0    **
 **  International License. See LICENSE.md.         **
 **                                                 **
 **  For information about this project see:        **
 **  http://www.cedmav.com/pidx                     **
 **  or contact: pascucci@sci.utah.edu              **
 **  For support: PIDX-support@visus.net            **
 **                                                 **
 *****************************************************/

#include "../PIDX_inc.h"

unsigned int getNumBits ( unsigned int v )
{
  return (unsigned int)floor((log2(v))) + 1;
}

unsigned long long getPowerOf2(int x)
{
  /*  find the power of 2 of an integer value (example 5->8) */
  int n = 1;
  while (n < x) n <<= 1;
  return n;
}

unsigned int getLevelFromBlock (unsigned long long block, int bits_per_block)
{
  if (block == 0)
    return 0;
  else
    return (unsigned int)floor((log2(getPowerOf2(block)))) + bits_per_block;
  
  return 0;
}

unsigned int getLeveL (unsigned long long index)
{
  if (index)
    return (unsigned int)floor((log2(index))) + 1;
  /* To deal with index 0 - level 0 */
  return 0;
}

int isValidBox(int** box)
{
  int D;
  for(D = 0 ; D < PIDX_MAX_DIMENSIONS ; D++)
  {
    //printf("VFY: %d %d\n", box[0][D], box[1][D]);
    if (! (box[0][D]>=0 && box[0][D]<=box[1][D]))
      return 0;
  }
  return 1;
}

void Deinterleave(const char* bitmask, int maxh, unsigned long long zaddress, int* point)
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

unsigned long long ZBitmask(const char* bitmask,int maxh)
{
  return ((unsigned long long)1)<<maxh;
}

unsigned long long ZStart(const char* bitmask,int maxh,int BlockingH)
{
  if (!BlockingH) 
    return 0;
  assert(BlockingH>=1 && BlockingH<=maxh);
  return ((unsigned long long)1)<<(maxh-BlockingH);
}

unsigned long long ZEnd(const char* bitmask,int maxh,int BlockingH)
{
  if (!BlockingH) 
    return 0;
  assert(BlockingH>=1 && BlockingH<=maxh);
  return (ZBitmask(bitmask,maxh)-1)-(ZStart(bitmask,maxh,BlockingH)-1);
}

void ZDelta(const char* bitmask, int maxh, int BlockingH, int* point)
{
  int K, bit;
  for(K = 0; K < PIDX_MAX_DIMENSIONS; K++)
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
  for(i = 0 ; i < PIDX_MAX_DIMENSIONS ; i++)
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
  for( i = 0 ; i < PIDX_MAX_DIMENSIONS ; i++)
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
  unsigned long long i;
  char cpstr[strlen(str)+1];
  for(i=0; i < (int)strlen(str); i++)
    cpstr[i] = str[strlen(str)-i-1];
  
  cpstr[i] = '\0';
  strcpy(str, cpstr);
}

void GuessBitmaskPattern(char* _bits, PointND dims)
{
  int D,N,ordered;
  int dim = 1;
  char* p=_bits;
	      
  PointND id,sorted_id;

  *p++='V';

  For(D)
  {
    PGET(dims,D)=( int)getPowerOf2(PGET(dims,D));
    PGET(id,D)=D;
  }

  //order is ASC order (from smaller dimension to bigger)
  for (ordered=0;!ordered;)
  {
    ordered=1;
    OffsetFor(D,0,-1)
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
  
  For(D)
  {
    //order in DESC order
    for (ordered=0,sorted_id=id;!ordered;)
    {
      ordered=1;
      OffsetFor(N,D,-1)
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
      OffsetFor(N,D,0)
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

  size_t size = 0;
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
  size_t i = 0;
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

  size_t size = 0;
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
  size_t i = 0;
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

  size_t size = 0;
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
  size_t i = 0;
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

  size_t size = 0;
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
  size_t i = 0;
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
  
  if(!H)
  {
    h_box[0][0] = h_box[0][1] = h_box[0][2] = h_box[0][3] = h_box[0][4] = 0;
    h_box[1][0] = h_box[1][1] = h_box[1][2] = h_box[1][3] = h_box[1][4] = 0;
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
  for( i = 0 ; i < PIDX_MAX_DIMENSIONS ; i++)
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

void Hz_to_xyz(const char* bitmask,  int maxh, unsigned long long hzaddress, unsigned long long* xyz)
{
  unsigned long long lastbitmask=((unsigned long long)1)<<maxh;
  
  hzaddress <<= 1;
  hzaddress  |= 1;
  while ((lastbitmask & hzaddress) == 0) hzaddress <<= 1;
    hzaddress &= lastbitmask - 1;
  
  PointND cnt;
  PointND p  ;
  int n = 0;

  memset(&cnt,0,sizeof(PointND));
  memset(&p  ,0,sizeof(PointND));

  for (;hzaddress; hzaddress >>= 1,++n, maxh--) 
  {
    int bit= bitmask[maxh];
    PGET(p,bit) |= (hzaddress & 1) << PGET(cnt,bit);
    ++PGET(cnt,bit);
  }
  xyz[0] = p.x;
  xyz[1] = p.y;
  xyz[2] = p.z;
  xyz[3] = p.u;
  xyz[4] = p.v;
}

unsigned long long xyz_to_HZ(const char* bitmask, int maxh, PointND xyz)
{
  unsigned long long zaddress=0;
  int cnt   = 0;
  PointND zero;
  int temp_maxh = maxh;
  memset(&zero,0,sizeof(PointND));

  for (cnt=0 ; memcmp(&xyz, &zero, sizeof(PointND)) ; cnt++, maxh--)
  {
    int bit= bitmask[maxh];
    zaddress |= ((unsigned long long)PGET(xyz,bit) & 1) << cnt;
    PGET(xyz,bit) >>= 1;
  }
  
  unsigned long long lastbitmask=((unsigned long long)1)<<temp_maxh;
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
