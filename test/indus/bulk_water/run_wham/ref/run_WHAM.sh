#!/bin/bash -l

# Program
wham="${HOME}/source/wham/build/bin/wham"

# Input
wham_options="wham_options.input"

#valgrind \
$wham  $wham_options  $data_summary \
       $biasing_parameters  $aux_var_files
