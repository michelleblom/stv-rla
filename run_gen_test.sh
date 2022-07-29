
instances=(data/US/2014_OAK_MAYOR.us data/US/BK_CC_D8_2014.us data/US/test_Aspen09CC_ballots.us data/US/test_Pierce08CAS_ballots.us data/Ireland/DublinNorth2002_ballots.txt data/Ireland/DublinWest2002_ballots.txt data/Ireland/GCC_07_Anderson_ballots.txt data/Ireland/GCC_07_Baillieston_ballots.txt data/Ireland/GCC_07_Canal_ballots.txt data/Ireland/GCC_07_Craigton_ballots.txt data/Ireland/GCC_07_Drumchapel_ballots.txt data/Ireland/GCC_07_EastCentre_ballots.txt data/Ireland/GCC_07_Garscadden_ballots.txt data/Ireland/GCC_07_Govan_ballots.txt data/Ireland/GCC_07_GreaterPollock_ballots.txt data/Ireland/GCC_07_Hillhead_ballots.txt data/Ireland/GCC_07_Langside_ballots.txt data/Ireland/GCC_07_Linn_ballots.txt data/Ireland/GCC_07_Maryhill_ballots.txt data/Ireland/GCC_07_Newlands_ballots.txt data/Ireland/GCC_07_NorthEast_ballots.txt data/Ireland/GCC_07_Partick_ballots.txt data/Ireland/GCC_07_Pollockshields_ballots.txt data/Ireland/GCC_07_Shettleston_ballots.txt data/Ireland/GCC_07_SouthsideCentral_ballots.txt data/Ireland/GCC_07_Springburn_ballots.txt)

for instance in "${instances[@]}"; do 
    bn="${instance%%.*}"
    voters=`cat "${bn}.voters"`
    quota=`cat "${bn}.quota"`
    rl=0.10
    er=0.002
    
    python3 audit_general.py -d ${instance} -outcome ${bn}.otc_D2 -quota ${quota} -voters ${voters} -e $er -r $rl -log "${bn}.${rl}.${er}.429.log" -maxtv 0.429 
done
