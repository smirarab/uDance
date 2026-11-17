#!/usr/bin/env bash
export SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# create a conda environment including the tools benchmarked.
# channel order matters here: conda-forge must come before bioconda, since this repo
# has been seen paired with a user ~/.condarc that sets channel_priority: strict - with
# that setting, listing bioconda first can make the solver ignore conda-forge entirely
# for ordinary (non-bioinformatics) packages like python/urllib3/requests, producing
# bogus "unsolvable" errors that have nothing to do with uDance's actual pins.
# mafft is there for WITCH (pip-installed below), which uses the mafft found on PATH.
if conda env list | awk '{print $1}' | grep -xq "udance"; then
        echo "conda environment udance exists"
else
	basecomm=$(echo conda create -y -c conda-forge -c bioconda --channel smirarab --name udance python=3.9  pip newick_utils=1.6 setuptools seqkit=2.1.0 scipy dendropy=4.5.2 pandas=1.3.0 snakemake raxml=8.2.12 iqtree=2.1.2 treeshrink=1.3.9 fasttree=2.1.10 julia=1.7.1 gappa=0.7.1 trimal=1.4.1 mafft)
	condaplat=$(conda info | grep "platform" | awk '{print $3}')
	if [ "$condaplat" == "osx-arm64" ] ; then
		echo "Installing conda environment on Mac OS X Apple Chip platform"
		echo "tqdist will not be installed. Installing tqdist package is optional for uDance."
		CONDA_SUBDIR=osx-64 $basecomm
	elif [ "$condaplat" == "osx-64" ] ; then
		echo "Installing conda environment on Mac OS X platform"
                echo "tqdist will not be installed. Installing tqdist package is optional for uDance."
		$basecomm
	elif [ "$condaplat" == "linux-64" ] ; then
		echo "Installing conda environment on Linux 64-bit platform"
                $basecomm tqdist
	else
		echo "Your platform is not supported. The following conda platforms are supported: osx-arm64, osx-64, and linux-64."
		exit 1
	fi

        source activate udance
        conda activate udance
        pip install apples==2.0.11 kmeans1d==0.3.1 witch-msa
        # julia needs to download packages during the first run. Let's do it now while we are online :)
        julia $SCRIPTS_DIR/uDance/correction_multi.jl datasmall/alignments/p0309.fasta > /dev/null
        # install astral
        pushd $SCRIPTS_DIR/uDance/tools/ASTRAL
       	./make.sh > .astralInstallation.log 2>&1
        java -D"java.library.path=lib/" -jar native_library_tester.jar
        popd
fi

# raxml-ng 2.0 is installed in its own, separate conda environment rather than into "udance"
# itself: through htslib it needs zlib>=1.3, which conflicts with the zlib<1.3 pins of
# gappa=0.7.1, iqtree=2.1.2 and julia=1.7.1 in the main environment - a genuine ABI-generation
# conflict that can't be resolved without changing the versions of those tools.
if conda env list | awk '{print $1}' | grep -xq "udance_raxmlng"; then
	echo "conda environment udance_raxmlng exists"
else
	conda create -y -c conda-forge -c bioconda --name udance_raxmlng "raxml-ng=2.0"
fi
