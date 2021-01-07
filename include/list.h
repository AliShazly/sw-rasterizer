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

void list_init(list *l, size_t n_bytes, size_t len);
void *list_index(list *l, size_t index);
void list_append(list *l, void *val);
void list_set(list *l, size_t index, void *val);
void list_clear(list *l);
void *list_array(list *l);
size_t list_used(list *l);
void list_free(list *l);

#endif
