all:
	g++ -I/usr/include/python3.8 -o mes main.cpp -lquadmath -lpython3.8

clear:
	rm mes
