#include <stdint.h>
#include <GL/glut.h>
#include <time.h>

#include "renderer.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define CTX global_render_context_ptr

RenderCtx *CTX;

void exit_func();
void keyboard_exit(unsigned char c, int x, int y);
void window_reshape(int width, int height);
void render_2d_texture();

int main(int argc, char** argv)
{
    srand(time(NULL));

    RenderCtx ctx = init_renderer();
    CTX = &ctx;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(1024, 1024);
    glutCreateWindow("renderer");

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glutReshapeFunc(window_reshape);
    glutDisplayFunc(render_2d_texture);
    glutKeyboardFunc(keyboard_exit);
    atexit(exit_func);

    glutMainLoop();
}

// glutMainLoop doesn't return control, need to cleanup this way
void exit_func()
{
    destroy_renderer(CTX);
}

void keyboard_exit(unsigned char c, int x, int y)
{
    if (c == 27)
    {
        exit(0);
    }
}

void window_reshape(int width, int height)
{
    int min = MIN(width, height);

    glViewport(0, 0, min, min);

    glEnable(GL_TEXTURE_2D);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
}


// draws a quad over the screen, writes the buffer to it as a texture
void render_2d_texture()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_TEXTURE_2D);

    draw_object(CTX->mesh->verts, CTX->mesh->size, CTX);

    glTexImage2D(GL_TEXTURE_2D,
               0,                    // level 0
               3,                    // use only R, G, and B components
               CTX->cols, CTX->rows, // texture width x height
               0,                    // no border
               GL_RGB,               // texels are in RGB format
               GL_UNSIGNED_BYTE,     // color components are unsigned bytes
               CTX->buffer);

    // texcoords are arranged to fit 2D arrays
    glBegin(GL_QUADS);
        /* glColor3f(1, 0, 0); */
        glTexCoord2f(0.0, 1.0);
        glVertex2f(-1, -1);

        /* glColor3f(rand() / (double)RAND_MAX,rand() / (double)RAND_MAX, rand() / (double)RAND_MAX); */
        glTexCoord2f(1.0, 1.0);
        glVertex2f(1, -1);

        /* glColor3f(rand() / (double)RAND_MAX,rand() / (double)RAND_MAX, rand() / (double)RAND_MAX); */
        glTexCoord2f(1.0, 0.0);
        glVertex2f(1, 1);

        /* glColor3f(rand() / (double)RAND_MAX,rand() / (double)RAND_MAX, rand() / (double)RAND_MAX); */
        glTexCoord2f(0.0, 0.0);
        glVertex2f(-1, 1);
    glEnd();

    glutSwapBuffers();
    glutPostRedisplay();
}

