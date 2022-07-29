
instances=(data/US/2014_OAK_MAYOR.us data/US/BK_CC_D8_2014.us data/US/test_Aspen09CC_ballots.us data/US/test_Pierce08CAS_ballots.us)

for instance in "${instances[@]}"; do 
    bn="${instance%%.*}"
    voters=`cat "${bn}.voters"`
    quota=`cat "${bn}.quota"`
    rl=0.10
    er=0.002
    
    python3 audit_general.py -d ${instance} -outcome ${bn}.otc_D2 -quota ${quota} -voters ${voters} -e $er -r $rl -log "${bn}.${rl}.${er}.429.TEST.log" -maxtv 0.429 
done
