#!/bin/bash

declare -a files=( \
	# Optimal biasing free energies from WHAM
	"f_bias_WHAM.out" \
	# F_0(Ntilde)
	"F_Ntilde_WHAM.out" \
	"F_Ntilde_biased.out" \
	"F_Ntilde_unbiased.out" \
	"F_Ntilde_rebiased.out" \
	# F_0(Ntilde,N)
	"F_Ntilde_N_WHAM.out" \
	"samples_Ntilde_N.out" \
	# F_0(N)
	"F_N_WHAM.out" \
	"F_N_biased.out" \
	"F_N_unbiased.out" \
	"F_N_rebiased.out" \
)

echo_failed_diffs=1

for output_file in ${files[@]}; do
	diff_out=$( diff $output_file ref/$output_file )
	if [[ ! -z $diff_out ]]; then
		echo "FAILED"
		echo "  diff $output_file ref/$output_file"

		if [[ $echo_failed_diffs -eq 1 ]]; then
			diff $output_file ref/$output_file
			echo ""
		fi
	fi
done
