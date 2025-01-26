#!/bin/bash

# sudo apt install dos2unix
# dos2unix unipept_pipeline.bash

help(){
cat << EOF

Usage: bash $0 -i [Input directory] -o [Output directory]
       mandatory arguments:
       -i [Input directory]
       -o [Output directory]

       optional arguments:
       -r [remove duplicated peptides: yes (default)/no]
       -t [lca analysis: yes (default)/no]
       -f [functional analysis: yes (default)/no]
       -fo [functional option: go (default)/ec/interpro/all]

EOF
exit
}

# Parse arguments
while getopts "i:o:r:t:f:fo:h" option; do
	case $option in
		i) indir=$OPTARG;;
		o) outdir=$OPTARG;;
		r) rmdup=$OPTARG;;
		t) taxonomy=$OPTARG;;
		f) functional=$OPTARG;;
		fo) functional_option=$OPTARG;; 
		h) help;;
	esac
done

# Default values for arguments not defined
if [[ -z $rmdup ]]; then
	rmdup="yes"
fi
if [[ -z $taxonomy ]]; then
	taxonomy="yes"
fi
if [[ -z $functional ]]; then
	functional="yes"
fi
if [[ $functional == "yes" ]] && [[ -z $functional_option ]]; then
	functional_option="go"
fi

# Check mandatory arguments
if [[ -z $indir ]] || [[ -z $outdir ]]; then
	help
fi
if [[ $functional == "yes" ]] && [[ -z $functional_option ]]; then
	echo ""
	echo "argument -fo [functional_option] is mandatory for functional analysis"
	help
fi

# Create output directory
mkdir -p "$outdir"
echo "$outdir created"

# Process each file in the input directory
for file in "$indir"/*; do
	echo "$file"
	# Get the base name of the file
	filename=$(basename $file | cut -d'.' -f1 )

	echo "Processing $filename ..."

	# Taxonomic analysis
	if [[ $taxonomy == "yes" ]]; then
		echo "LCA analysis ..."
		if [[ $rmdup == "yes" ]]; then
			cat "$file" | prot2pept | peptfilter | tr I L | sort -u | unipept pept2lca -a > "$outdir/${filename}_lca_output.csv"
		else
			cat "$file" | prot2pept | peptfilter | tr I L | unipept pept2lca -a > "$outdir/${filename}_lca_output.csv"
		fi

		# Check if the command was successful
		if [[ $? -eq 0 ]]; then
			echo "LCA analysis: SUCCESS"
			echo "Output saved to: $outdir/${filename}_lca_output.csv"
		else
			echo "LCA analysis: FAIL"
			exit 1
		fi
	fi

	# Functional analysis (if enabled)
done
