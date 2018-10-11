#!/bin/bash -l

# Program
wham="${HOME}/source/wham/bin/wham"

# Input
wham_options="wham_options.input"
data_summary="../data/data_summary.log"
biasing_parameters="biasing_parameters.log"
aux_var_files="../data/data_summary.log"

#valgrind \
$wham  $wham_options  $data_summary \
       $biasing_parameters  $aux_var_files
