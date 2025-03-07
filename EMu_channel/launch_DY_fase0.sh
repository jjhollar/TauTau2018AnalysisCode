voms-proxy-init --rfc --voms cms --valid 172:00

g++ -o nanotry_DY nanotry_DY.cpp `root-config --cflags --glibs`

for i in {158..204}
do
  	echo "$i" | ./nanotry_DY
        clear
done
