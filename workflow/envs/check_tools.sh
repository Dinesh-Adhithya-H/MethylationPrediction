#!/bin/bash

# Specify required versions
required_python_version="3.10.10"
required_snakemake_version="7.25.2"
required_samtools_version="1.18"
required_bedtools_version="2.30.0"


# Function to compare versions
version_compare() {
    local required_version=$1
    local installed_version=$2

    if [[ "$(printf '%s\n' "$installed_version" "$required_version" | sort -V | head -n1)" == "$required_version" ]]; then
        return 0  # Installed version is equal or greater
    else
        return 1  # Installed version is older
    fi
}

# Check Python version
if command -v python3 &>/dev/null; then
    python_version=$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")

    if version_compare "$required_python_version" "$python_version"; then
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
    snakemake_version=$(snakemake --version | awk '{print $1}')

    if version_compare "$required_snakemake_version" "$snakemake_version"; then
        echo "Snakemake version is sufficient: $snakemake_version"
    else
        echo "Snakemake version $required_snakemake_version or higher is required, but found $snakemake_version"
        exit 1
    fi
else
    echo "Snakemake is required but not found"
    exit 1
fi



# Check Bedtools version
if bedtools --version &>/dev/null; then
    bedtools_version=$(bedtools --version | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')

    if version_compare "$required_bedtools_version" "$bedtools_version"; then
        echo "Bedtools version is sufficient: $bedtools_version"
    else
        echo "Bedtools version $required_bedtools_version or higher is required, but found $bedtools_version"
        exit 1
    fi
else
    echo "Bedtools is required but not found"
    exit 1
fi



# Check Samtools version
if samtools --version &>/dev/null; then
    samtools_version=$(samtools --version | grep -oE 'samtools [0-9]+\.[0-9]+' | awk '{print $2}')


    if version_compare "$required_samtools_version" "$samtools_version"; then
        echo "Samtools version is sufficient: $samtools_version"
    else
        echo "Samtools version $required_samtools_version or higher is required, but found $samtools_version"
        exit 1
    fi
else
    echo "Samtools is required but not found"
    exit 1
fi


echo "All required tools with the correct versions are available!"
