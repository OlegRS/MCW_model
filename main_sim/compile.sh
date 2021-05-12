#!/bin/bash

rm ../bin/*
[ -d "../data" ] || mkdir "../data"
[ -d "../bin" ] || mkdir "../bin"
[ -d "../data/CW_model" ] || mkdir "../data/CW_model"
[ -d "../data/CW_model/h_scanning" ] || mkdir "../data/CW_model/h_scanning"
[ -d "../data/CW_model/T_scanning" ] || mkdir "../data/CW_model/T_scanning"
[ -d "../data/2CW_model" ] || mkdir "../data/2CW_model"
[ -d "../data/2CW_model/h1_scanning" ] || mkdir "../data/2CW_model/h1_scanning"
[ -d "../data/2CW_model/h1=h2_scanning" ] || mkdir "../data/2CW_model/h1=h2_scanning"
[ -d "../data/2CW_model/T_scanning" ] || mkdir "../data/2CW_model/T_scanning"

g++ -std=c++11 CW_model_h_scan.cpp ../MCW_lib/*.cpp -Ofast -o ../bin/CW_model_h_scan
g++ -std=c++11 2CW_crit.cpp ../MCW_lib/*.cpp -Ofast -o ../bin/2CW_crit
g++ -std=c++11 2CW_bistability.cpp ../MCW_lib/*.cpp -Ofast -o ../bin/2CW_bistability
g++ -std=c++11 2CW_temperature_dependence.cpp ../MCW_lib/*.cpp -Ofast -o ../bin/2CW_temperature_dependence
g++ -std=c++11 CW_temperature_dependence.cpp ../MCW_lib/*.cpp -Ofast -o ../bin/CW_temperature_dependence
