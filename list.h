#ifndef LIST_H
#define LIST_H

#include <stddef.h>

#define LIST_GROWTH_RATE 1.5

typedef struct
{
    void *array;
    size_t n_bytes;
    size_t len;
    size_t used;
}list; // Dynamic array

void init_list(list *l, size_t n_bytes, size_t len);
void append_to_list(list *l, void *val);
void *index_list(list *l, size_t index);

#endif
