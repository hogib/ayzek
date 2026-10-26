#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int fp_reference(double dt, const float *x, int n, int packet,
                 double filter_window, double long_term_window,
                 double threshold1, double threshold2, double t_up,
                 double (*out)[7], int max_picks, double *cf);

#ifdef __cplusplus
}
#endif
