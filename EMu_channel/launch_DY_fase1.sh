voms-proxy-init --rfc --voms cms --valid 172:00

g++ -o nanotry_fase1_DY  nanotry_fase1_DY.cpp  `root-config --cflags --glibs`

for i in {1..204}
do
  	echo "$i" | ./nanotry_fase1_DY >> cortes_cinematicos_fase1_DY.txt
        clear
done
