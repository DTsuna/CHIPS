mass=18 # FIXME CHANGE THIS
name=tsuna # FIXME CHANGE THIS

today=$(date "+%m%d")
for f in 0.3 0.4 0.5 0.6 0.8; do
	for t in 03 10 30; do
		cp -r ~/mesa-r12778/star/CHIPS ./
		mv CHIPS ${mass}Msun_finj${f}_t${t}yr_${today}_${name}
		cd ${mass}Msun_finj${f}_t${t}yr_${today}_${name}
		cp /SET/PATH/TO/intermediate${t}yr.txt EruptionFiles/ #FIXME CHANGE THIS
		python3 after_eruption.py --zams-m ${mass} --zams-z 1 --profile-at-cc EruptionFiles/intermediate${t}yr.txt --analytical-CSM
		cd ../
	done
done
