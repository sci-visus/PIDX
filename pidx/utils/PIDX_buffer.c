#include "../PIDX_inc.h"

PIDX_buffer PIDX_buffer_create_empty()
{
  PIDX_buffer b = { NULL, 0, 0 };
  return b;
}
PIDX_buffer PIDX_buffer_create_with_capacity(size_t capacity)
{
  PIDX_buffer b = { malloc(capacity), 0, capacity };
  return b;
}
void PIDX_buffer_free(PIDX_buffer *b)
{
  free(b->buffer);
  b->size = 0;
  b->capacity = 0;
}
void PIDX_buffer_append(PIDX_buffer *b, const unsigned char *data, const size_t size)
{
  if (b->capacity - b->size < size) {
    b->capacity = Max2ab(b->capacity + size, b->capacity * 1.5);
    b->buffer = realloc(b->buffer, b->capacity);
  }
  memcpy(b->buffer, data, size);
  b->size += size;
}
void PIDX_buffer_resize(PIDX_buffer *b, const size_t size)
{
  if (b->capacity < size) {
    b->capacity = size;
    b->buffer = realloc(b->buffer, b->capacity);
  }
  b->size = size;
}

