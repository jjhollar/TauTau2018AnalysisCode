voms-proxy-init --rfc --voms cms --valid 172:00

g++ -o sinal  sinal.cpp  `root-config --cflags --glibs`

for i in {1..1}
do
  	echo "$i" | ./sinal >> cortes_sinal_JHtest_1.txt
        clear
done

