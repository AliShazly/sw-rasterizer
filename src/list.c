#include "list.h"

#include <string.h>
#include <assert.h>
#include <stdlib.h>

void init_list(list *l, size_t n_bytes, size_t len)
{
    l->n_bytes = n_bytes;
    l->array=malloc(n_bytes * len);
    assert(l->array != NULL);
    l->len = len;
    l->used = 0;
}

void append_to_list(list *l, void *val)
{
    if (l->len == l->used)
    {
        l->len *= LIST_GROWTH_RATE;
        l->array = realloc(l->array, l->len * l->n_bytes);
        assert(l->array != NULL);
    }
    void *dest = index_list(l, l->used);
    memcpy(dest, val, l->n_bytes);
    l->used++;
}

void *index_list(list *l, size_t index)
{
    // casting to char* for pointer arithmetic on a void*
    return (char*)l->array + (l->n_bytes * index);
}

