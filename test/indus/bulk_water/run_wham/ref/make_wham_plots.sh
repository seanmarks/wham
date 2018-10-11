#!/bin/bash

# Run from the directory containing the WHAM output files

T=273.0
# x = lambda-tilde_CH
avg_x_0=29.8188
std_dev_x_0=2.36484
avg_n_v_0=29.8188
std_dev_n_v_0=2.36484

scripts_dir="${HOME}/scripts/bias_lambda/wham"

python "${scripts_dir}/make_wham_plots.py" -T $T \
       -avg_x_0   $avg_x_0   -std_dev_x_0   $std_dev_x_0 \
       -avg_n_v_0 $avg_n_v_0 -std_dev_n_v_0 $std_dev_n_v_0

# FIXME
exit 0

python "${scripts_dir}/make_wham_reweighted_plots.py" -T $T \
       -avg_x_0   $avg_x_0   -std_dev_x_0   $std_dev_x_0 \
       -avg_n_v_0 $avg_n_v_0 -std_dev_n_v_0 $std_dev_n_v_0
