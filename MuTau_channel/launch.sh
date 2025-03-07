voms-proxy-init --rfc --voms cms --valid 172:00

for i in {87..354}
do
  	echo "$i" | ./nanotry_ttjets >> cortes_ttjets_2018_UL_ETau_1.txt
        clear
done
