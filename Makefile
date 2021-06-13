gcc_options = -std=c++17 -Wall -O2 --pedantic-errors

iss_sgp4_teme: iss_sgp4_teme.o sgp4.o tle.o time.o
	g++ $(gcc_options) -o $@ $^

iss_sgp4_teme.o : iss_sgp4_teme.cpp
	g++ $(gcc_options) -c $<

sgp4.o : sgp4.cpp
	g++ $(gcc_options) -c $<

tle.o : tle.cpp
	g++ $(gcc_options) -c $<

time.o : time.cpp
	g++ $(gcc_options) -c $<

run : iss_sgp4_teme
	./iss_sgp4_teme

clean :
	rm -f ./iss_sgp4_teme
	rm -f ./*.o

.PHONY : run clean

