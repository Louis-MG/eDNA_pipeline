#!/bin/bash

echo -e "
This script will install the mamba environment necessary to run the pipeline.
"

echo -e "Checking mamba install ..."

mamba_version=$(mamba --version)

if [ "$mamba_version" == "" ];
then
	echo -e "\033[0;31mMamba not found, exiting.\033[0m"
	exit 1
else
	echo -e "\033[0;32mMamba found, continue ...\033[0m"
fi

for yaml in ./*.yml;
do
	echo -e "\033[0;32mInstalling environment ${yaml} ...\033[0m"
	mamba create -f "$yaml" --yes &>> mamba_installation.logs
done


