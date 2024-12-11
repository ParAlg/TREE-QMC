echo "COMPILING..."
# g++ -std=c++17 -g -I external/MQLib/include -I external/toms743 -I external/parlaylib/include -o tree-qmc -pthread -mcx16 -march=native src/*.cpp external/toms743/toms743.cpp external/MQLib/bin/MQLib.a -lm -DVERSION=\"$(cat version.txt)\"
g++ -std=c++17 -O3 -I external/MQLib/include -I external/toms743 -I external/parlaylib/include -o tree-qmc -pthread -mcx16 -march=native src/*.cpp external/toms743/toms743.cpp external/MQLib/bin/MQLib.a -lm -DVERSION=\"$(cat version.txt)\"

echo "RUNNING..."
# gdb --args ./tree-qmc -i tutorial/gene-trees/avian_uce_trees_3679.tre
./tree-qmc -i tutorial/gene-trees/avian_uce_trees_3679.tre