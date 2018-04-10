#ifndef __PIDX_WINDOWS_DEFINE_H
#define __PIDX_WINDOWS_DEFINE_H

/** Macros to swap bytes in a multi-byte value, to convert from big-endian data
to little-endian data and vice versa. These are taken from the Boost library.*/
#ifndef __has_builtin
#define __has_builtin(x) 0 // Compatibility with non-clang compilers
#endif
#define htonl(x) _byteswap_ulong(x)
#define ntohl(x) htonl(x)

#include <io.h>
#include <stdio.h>
#include <string.h>

typedef long int __off_t;
typedef long int __off64_t;
typedef __off_t off_t;

#define ssize_t size_t
#define PATH_MAX 256
#define mkdir(x,y) mkdir(x) // This returns error if directory exists already
#define PIDX_HAVE_MPI 1  // TODO: remove this later cause it is going to be deprecated

inline int pwrite(int fd, const void *buf, size_t nbytes, off_t offset)
{
  long ret = _lseek(fd, offset, SEEK_SET);

  if (ret == -1) {
    return(-1);
  }
  return(_write(fd, buf, nbytes));
}

inline int pread(int fd, void *buf, size_t nbytes, off_t offset)
{
  if (_lseek(fd, offset, SEEK_SET) != offset) {
    return -1;
  }
  return read(fd, buf, nbytes);
}

#endif
