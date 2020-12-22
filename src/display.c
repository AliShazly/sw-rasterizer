#include "context.h"
#include "rasterize.h"

#include <GL/glut.h>
#include <sys/time.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define CTX global_render_context_ptr

RenderCtx *CTX;

static unsigned long long int global_frame_count = 1;
static unsigned long long int global_usec_sum = 0;

void exit_func();
void keyboard_exit(unsigned char c, int x, int y);
void window_reshape(int width, int height);
void render_2d_texture();

int main(int argc, char** argv)
{
    RenderCtx ctx = init_renderer();
    CTX = &ctx;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(CTX->cols, CTX->rows);
    glutCreateWindow("rasterizer");

    glutReshapeFunc(window_reshape);
    glutDisplayFunc(render_2d_texture);
    glutKeyboardFunc(keyboard_exit);
    atexit(exit_func);
    glutMainLoop();
}

// glutMainLoop doesn't return control, need to cleanup this way
void exit_func()
{
    long double usec_avg = global_usec_sum / (long double)global_frame_count;
    double sec = (usec_avg / 1000000.);
    printf("avg FPS: %f\n\n", 1/sec);
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
}

// draws a quad over the screen, writes the buffer to it as a texture
void render_2d_texture()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_TEXTURE_2D);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);


    struct timeval tval_before, tval_after, tval_result;
    gettimeofday(&tval_before, NULL);

    // drawing the frame
    /* clear_buffers(CTX); */
    /* draw_grid(CTX); */
    /* draw_object_wireframe(CTX); */
    draw_object_threads(CTX);

    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);

    global_frame_count+=1;
    global_usec_sum+=tval_result.tv_usec; // usec/1000000 for seconds


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
        glTexCoord2f(0.0, 1.0);
        glVertex2f(-1, -1);

        glTexCoord2f(1.0, 1.0);
        glVertex2f(1, -1);

        glTexCoord2f(1.0, 0.0);
        glVertex2f(1, 1);

        glTexCoord2f(0.0, 0.0);
        glVertex2f(-1, 1);
    glEnd();

    glutSwapBuffers();
    glutPostRedisplay();
}

