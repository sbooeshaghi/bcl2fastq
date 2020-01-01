#!/bin/bash
# Get major.minor.patch version from src/configure.
B2F_MAJOR=`sed -r '/bcl2fastq_version_major=/I!d;s/.*"([0-9]*)"/\1/' src/configure `
B2F_MINOR=`sed -r '/bcl2fastq_version_minor=/I!d;s/.*"([0-9]*)"/\1/' src/configure `
B2F_PATCH=`sed -r '/bcl2fastq_version_patch=/I!d;s/.*"([0-9]*)"/\1/' src/configure `
# Set build version in src/configure.
sed -i.bak 's:bcl2fastq_version_build="[0-9]*":bcl2fastq_version_build="422":' src/configure
"##teamcity[buildNumber '$B2F_MAJOR.$B2F_MINOR.$B2F_PATCH.422']"
chmod +x src/configure
