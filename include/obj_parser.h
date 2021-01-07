#ifndef OBJ_PARSER_H
#define OBJ_PARSER_H

#define OBJ_LIST_SIZE 2048
#define OBJ_FILE_BUF_SIZE (128 * 128)

#include <stdio.h>
#include "list.h"

void parse_obj(char *filename, size_t *out_size,
        double (**out_verts)[3],
        double (**out_texcoords)[2],
        double (**out_normals)[3]);

#endif
