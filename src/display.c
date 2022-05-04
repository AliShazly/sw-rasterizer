#include "context.h"
#include "rasterize.h"
#include "utils.h"
#include <MiniFB.h>
#include <stdlib.h>

void resize_cb(struct mfb_window *window, int width, int height);
void keyboard_cb(struct mfb_window *window, mfb_key key, mfb_key_mod mod,
                 bool pressed);

// clang-format off
static vec3 lowercase_to_vec3[26][3] = {
    {-0.1, 0, 0},  // 'a'
    {-1}, {-1},
    {0.1, 0, 0}, {0, 0.1, 0,},  // 'd' and 'e'
    {-1}, {-1}, {-1}, {-1}, {-1}, {-1},
    {-1}, {-1}, {-1}, {-1}, {-1},
    {0, -0.1, 0},  // 'q'
    {-1},
    {0, 0, -0.1},  // 's'
    {-1}, {-1}, {-1},
    {0, 0, 0.1},  // 'w'
    {-1}, {-1}, {-1}
};
// clang-format on

int main(int argc, char **argv) {
    if (argc != 2) {
        printf("Usage: %s <input_file>\n", argv[0]);
        return 1;
    }
    RenderCtx ctx = init_renderer(argv[1]);

    struct mfb_window *window =
        mfb_open_ex("Rasterizer", ctx.cols, ctx.rows, WF_RESIZABLE);
    if (!window) {
        fprintf(stderr, "Failed to open window\n");
        return 1;
    }

    mfb_set_user_data(window, &ctx);
    mfb_set_resize_callback(window, resize_cb);
    mfb_set_keyboard_callback(window, keyboard_cb);

    uint32_t *mfb_buffer = calloc(ctx.cols * ctx.rows, sizeof(uint32_t));

    do {
        clear_buffers(&ctx);
        draw_object(&ctx);

        for (int i = 0; i < ctx.cols * ctx.rows; i++) {
            unsigned char *pixel = ctx.buffer[i];
            mfb_buffer[i] = MFB_RGB(pixel[0], pixel[1], pixel[2]);
        }

        int state = mfb_update_ex(window, mfb_buffer, ctx.cols, ctx.rows);
        if (state < 0) {
            window = NULL;
            break;
        }
    } while (mfb_wait_sync(window));

    free(mfb_buffer);
    destroy_renderer(&ctx);
}

void resize_cb(struct mfb_window *window, int width, int height) {
    int dim = imin(width, height);
    mfb_set_viewport(window, 0, 0, dim, dim);
}

void keyboard_cb(struct mfb_window *window, mfb_key key, mfb_key_mod mod,
                 bool pressed) {
    if (pressed && key >= 'A' && key <= 'Z') {
        unsigned char idx_by_key = key - 'A';
        double *offset = *lowercase_to_vec3[idx_by_key];

        if (*offset != -1) {
            RenderCtx *ctx = mfb_get_user_data(window);
            move_camera(ctx, offset);
        }
    }
}
