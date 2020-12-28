all:
	g++ -I/usr/include/python3.8 -o fem main.cpp -lquadmath -lpython3.8

clear:
	rm fem
