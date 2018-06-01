#ifndef TYPES_H
#define TYPES_H

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;

#if __STDC_VERSION__ >= 199901L
  #include <stdint.h>
  typedef int8_t zfp_int8;
  typedef uint8_t zfp_uint8;
  typedef int16_t zfp_int16;
  typedef uint16_t zfp_uint16;
  typedef int32_t zfp_int32;
  typedef uint32_t zfp_uint32;
  typedef uint64_t zfp_int64;
  typedef uint64_t zfp_uint64;
#else
  /* assume common integer types in C89 */
  typedef signed char zfp_int8;
  typedef unsigned char zfp_uint8;
  typedef signed short zfp_int16;
  typedef unsigned short zfp_uint16;
  typedef signed int zfp_int32;
  typedef unsigned int zfp_uint32;
  typedef signed long long zfp_int64; /* not ANSI C89 compliant */
  typedef unsigned long long zfp_uint64; /* not ANSI C89 compliant */
#endif

#endif
