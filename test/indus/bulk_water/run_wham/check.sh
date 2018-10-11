#!/bin/bash

declare -a files=( \
	"F_x_WHAM.out" \
	"F_x_biased.out" \
	"GlobalSampleCounts_x.out" \
	"F_x_unbiased.out" \
	"F_x_rebiased.out" \
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
