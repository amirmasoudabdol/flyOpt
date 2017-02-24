#!/bin/sh
# Run sequential Scatter Search (sSS) on a single workstation

fly="$HOME/code/flyOpt/fly/fly_ss"

opt="-s kr -i 0.2 -a 0.001"
fname="dm_hkgn53_sss"
outpath="$HOME/projects/2016_scat_vs_simann/dm_r_sss_0"

i=201
mkdir "$outpath/$i"
cp "$outpath/$fname" "$outpath/$i/$fname"_"$i"

#echo  ${opt} ${outpath}/${i}/${fname}_${i}
#gdb ${fly}
${fly} ${opt} ${outpath}/${i}/${fname}_${i}
exit

# outpath="$HOME/projects/2016_scat_vs_simann/dm_p_sss_limls"

# for i in $(seq -f "%03g" 0 9)
# do
# 	# Create folders with config file
# 	mkdir "$outpath/$i"
# 	cp "$outpath/$fname" "$outpath/$i/$fname"_"$i"

# 	# Start simulation 
# 	nohup ${fly} ${opt} ${outpath}/${i}/${fname}_${i} > ${outpath}/${i}/terminal_output.txt &
# done
# eof
