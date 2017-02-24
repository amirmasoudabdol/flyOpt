#!/bin/sh
# Convert directory structure with lots of info to simann-style output

# Use default simulation folder or one given by the user
if [ -z "$1" ]
then
    full_path="$HOME/projects/2016_scat_vs_simann/dm_p_sss_1/dm_hkgn58_sss"
else
	full_path=$1
fi

# Decompose simulation folder
base_fname=`basename ${full_path}`
scenario_path=`dirname ${full_path}`
project_path=`dirname ${scenario_path}`
new_scenario_path=${scenario_path}_simple

echo "Simplifying folder structure for:"
echo "- project path: ${project_path}"
echo "- scenario path: ${scenario_path}"
echo "- base simulation file: ${base_fname}"


# Make the new directory with a base config file
mkdir -p ${new_scenario_path}
cp ${full_path} ${new_scenario_path}

# For each directory, assuming the names are 000, 001, etc.
shopt -s dotglob
find ${scenario_path}/??? -prune -type d | while read run; do 

	run_idx=`basename ${run}`
	base_runfname=${base_fname}_${run_idx}

    # Copy best individual, log and times files to new scenario folder
    cp --remove-destination ${run}/${base_runfname}_ref_00 ${new_scenario_path}/${base_runfname}
    cp --remove-destination ${run}/${base_runfname}.log -t ${new_scenario_path}/
    cp --remove-destination ${run}/${base_runfname}.times -t ${new_scenario_path}/

    best=$(grep -A1 "Best Solution:" ${run}/terminal_output.txt | tail -n1) 
    echo "Score${best:7: -1}" > ${new_scenario_path}/${base_runfname}.best
done
# eof
