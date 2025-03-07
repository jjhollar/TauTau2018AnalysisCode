voms-proxy-init --rfc --voms cms --valid 172:00

g++ -o nanotry_data nanotry_data.cpp `root-config --cflags --glibs`

for i in {120..138}
do
  	echo "$i" | ./nanotry_data
        clear
done
