/* Runs Lomax's FilterPicker5 C function on a trace, packet by packet as in
 * picker_func_test_memory.c, for tests/test_filterpicker.cpp. The other files
 * in this directory are unchanged from picker_20130212.tgz (README.md); they
 * are compiled here as one unit because ew_bridge.h defines message_str. */

#include <limits.h>
#include <math.h>
#include <string.h>

#include "FilterPicker5.c"
#include "FilterPicker5_Memory.c"
#include "PickData.c"

#include "fp_reference.h"

/* Picks, with indices from the first sample, into out[i]: time, begin, end,
 * polarity, polarity weight, amplitude, period. The first packet holds one
 * sample, the rest `packet` samples. Returns the number of picks, or -1 if
 * there are more than max_picks. With cf, also writes the clipped
 * characteristic function of every sample there (a second pass). */
int fp_reference(double dt, const float *x, int n, int packet,
                 double filter_window, double long_term_window,
                 double threshold1, double threshold2, double t_up,
                 double (*out)[7], int max_picks, double *cf) {
    int pass, found = 0;
    float *buf = malloc(sizeof(float) * (size_t)n);
    for (pass = 0; pass < (cf ? 2 : 1); pass++) {
        FilterPicker5_Memory *mem = NULL;
        int start = 0, len, i;
        resultType = pass == 0 ? PICKS : CHAR_FUNC;
        memcpy(buf, x, sizeof(float) * (size_t)n);
        while (start < n) {
            PickData **list = NULL;
            int num = 0;
            len = start == 0 ? 1 : packet;
            if (start + len > n)
                len = n - start;
            Pick(dt, buf + start, len, filter_window, long_term_window,
                 threshold1, threshold2, t_up, &mem, TRUE_INT, &list, &num,
                 "test");
            for (i = 0; i < num && pass == 0; i++) {
                PickData *p = list[i];
                if (found == max_picks) {
                    free_PickList(list, num);
                    free_FilterPicker5_Memory(&mem);
                    free(buf);
                    return -1;
                }
                out[found][0] = start + (p->indices[0] + p->indices[1]) / 2.0;
                out[found][1] = start + p->indices[0];
                out[found][2] = start + p->indices[1];
                out[found][3] = p->polarity;
                out[found][4] = p->polarityWeight;
                out[found][5] = p->amplitude;
                out[found][6] = p->period;
                found++;
            }
            if (list)
                free_PickList(list, num);
            start += len;
        }
        if (pass == 1)
            for (i = 0; i < n; i++)
                cf[i] = buf[i];
        free_FilterPicker5_Memory(&mem);
    }
    free(buf);
    return found;
}
