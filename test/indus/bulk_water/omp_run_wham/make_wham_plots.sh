#!/bin/bash

# Run from the directory containing the WHAM output files

T=273.0
# x = Ntilde_v
avg_x_0=29.9657
std_dev_x_0=2.12469
avg_n_v_0=29.954
std_dev_n_v_0=2.258

scripts_dir="${HOME}/scripts/bias_lambda/wham"

python "${scripts_dir}/make_wham_plots.py" -T $T \
       -avg_x_0   $avg_x_0   -std_dev_x_0   $std_dev_x_0 \
       -avg_n_v_0 $avg_n_v_0 -std_dev_n_v_0 $std_dev_n_v_0

python "${scripts_dir}/make_wham_reweighted_plots.py" -T $T \
       -avg_x_0   $avg_x_0   -std_dev_x_0   $std_dev_x_0 \
       -avg_n_v_0 $avg_n_v_0 -std_dev_n_v_0 $std_dev_n_v_0
