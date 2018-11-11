#!/bin/bash

# Program
wham="${HOME}/source/wham/bin/wham"

# Input
wham_options="wham_options.input"

#valgrind \
$wham  $wham_options  $data_summary \
       $biasing_parameters  $aux_var_files
