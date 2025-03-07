voms-proxy-init --rfc --voms cms --valid 172:00

g++ -o nanotry_fase1_dados nanotry_fase1_dados.cpp  `root-config --cflags --glibs`

for i in {136..138}
do
  	echo "$i" | ./nanotry_fase1_dados >> cortes_cinematicos_fase1_dados.txt
        clear
done
