
instances=(data/Ireland/GCC_07_Langside_ballots.txt)

for instance in "${instances[@]}"; do 
    bn="${instance%%.*}"
    voters=`cat "${bn}.voters"`
    quota=`cat "${bn}.quota"`
    rl=0.10
    er=0.002
    
    python3 audit_general.py -d ${instance} -outcome ${bn}.otc_D2 -quota ${quota} -voters ${voters} -e $er -r $rl -log "${bn}.${rl}.${er}.429.logExtraArgs" -maxtv 0.429 
done
