#ifndef OBJ_PARSER_H
#define OBJ_PARSER_H

#define OBJ_LIST_SIZE 2048
#define OBJ_FILE_BUF_SIZE (128 * 128)

#include <stdio.h>
#include "list.h"

char *skip_whitespace(char *str);

int backtrack_fileptr(FILE *fp, char *buffer, size_t buf_size, char key);

void parse_coordinate(char *line, double ret[3], int n_values);
void parse_idx_cluster(char *start, int ret[3]);
void parse_idx_line(char *line_start_ptr, list *out_values);

void triangulate(list *face_verts, list *out_list);

void parse_obj(char *filename, size_t *out_size,
        double (**out_verts)[3],
        double (**out_texcoords)[2],
        double (**out_normals)[3]);

#endif
