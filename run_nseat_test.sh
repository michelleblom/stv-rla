arr=()

while read -r line; do
    arr+=($line)
done < 3seat_fwc_Scotland17_files.txt

while read -r line; do
    arr+=($line)
done < 4seat_fwc_Scotland17_files.txt

while read -r line; do
    arr+=($line)
done < 3seat_fwc_Scotland22_files.txt

while read -r line; do
    arr+=($line)
done < 4seat_fwc_Scotland22_files.txt



for f in "${arr[@]}"; do
    bn=`basename $f`
    bn="${bn%.*}"
    fe="${f%.*}"

    quota=`cat "${fe}.quota"`
    voters=`cat "${fe}.voters"`

    python3.9 audit_n_seats_fwc_simplified.py -d $f -outcome "${fe}.otc" -quota $quota -voters $voters -log "${fe}.ass_tv_0.log" 
done

