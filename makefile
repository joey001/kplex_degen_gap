all : 
	g++ -flto -Ofast -DNDEBUG -march=native -o KPLEX main.cpp Graph.cpp -w
clean:
	rm -rf KPLEX