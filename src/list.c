#include "list.h"

#include <string.h>
#include <assert.h>
#include <stdlib.h>

void list_init(list *l, size_t n_bytes, size_t len)
{
    l->n_bytes = n_bytes;
    l->array = malloc(n_bytes * len);
    assert(l->array != NULL);
    l->len = len;
    l->used = 0;
}

void *list_index(list *l, size_t index)
{
    assert(index >= 0);
    // casting to char* for pointer arithmetic on a void*
    return (char*)l->array + (l->n_bytes * index);
}

void list_set(list *l, size_t index, void *val)
{
    assert(index < l->used);
    void *dest = list_index(l, index);
    memcpy(dest, val, l->n_bytes);
}

void list_append(list *l, void *val)
{
    if (l->len == l->used)
    {
        l->len *= LIST_GROWTH_RATE;
        l->array = realloc(l->array, l->len * l->n_bytes);
        assert(l->array != NULL);
    }

    void *dest = list_index(l, l->used);
    memcpy(dest, val, l->n_bytes);
    l->used++;
}

void list_clear(list *l)
{
    l->used = 0;
}

void list_free(list *l)
{
    free(l->array);
}
