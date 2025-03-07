voms-proxy-init --rfc --voms cms --valid 172:00

g++ -o nanotry_ttJets nanotry_ttJets.cpp `root-config --cflags --glibs`

for i in {354..354}
do
  	echo "$i" | ./nanotry_ttJets
        clear
done
