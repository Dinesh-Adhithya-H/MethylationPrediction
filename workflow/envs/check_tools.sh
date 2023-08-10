#!/bin/bash

# Specify required versions
required_python_version="3.10"
required_snakemake_version="7.25.2"
required_samtools_version="1.8"
required_bamtools_version="2.5.1"

# Check Python version
if command -v python3 &>/dev/null; then
    python_version=$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")

    if [[ $(echo "$python_version >= $required_python_version" | bc -l) -eq 1 ]]; then
        echo "Python version is sufficient: $python_version"
    else
        echo "Python version $required_python_version or higher is required, but found $python_version"
        exit 1
    fi
else
    echo "Python 3 is required but not found"
    exit 1
fi

# Check Snakemake version
if snakemake --version &>/dev/null; then
    snakemake_version=$(snakemake --version | awk '{print $2}')

    if [[ "$snakemake_version" == "$required_snakemake_version" ]]; then
        echo "Snakemake version is sufficient: $snakemake_version"
    else
        echo "Snakemake version $required_snakemake_version is required, but found $snakemake_version"
        exit 1
    fi
else
    echo "Snakemake is required but not found"
    exit 1
fi

# Check Samtools version
if samtools --version &>/dev/null; then
    samtools_version=$(samtools --version | awk '{print $2}')

    if [[ "$samtools_version" == "$required_samtools_version" ]]; then
        echo "Samtools version is sufficient: $samtools_version"
    else
        echo "Samtools version $required_samtools_version is required, but found $samtools_version"
        exit 1
    fi
else
    echo "Samtools is required but not found"
    exit 1
fi

# Check Bamtools version
if bamtools --version &>/dev/null; then
    bamtools_version=$(bamtools --version | awk '{print $2}')

    if [[ "$bamtools_version" == "$required_bamtools_version" ]]; then
        echo "Bamtools version is sufficient: $bamtools_version"
    else
        echo "Bamtools version $required_bamtools_version is required, but found $bamtools_version"
        exit 1
    fi
else
    echo "Bamtools is required but not found"
    exit 1
fi

echo "All required tools with the correct versions are available!"
