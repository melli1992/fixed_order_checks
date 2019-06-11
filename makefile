#
#test-DY-num-v01: test-DY-num-v01.o DY-num-v01.o

#	g++ -o test-DY-num-v01 test-DY-num-v01.o DY-num-v01.o -lcuba -L/home/mbeekveld/Cuba/lib  -Wreturn-local-addr -L/home/mbeekveld/LHAPDF/lib -lLHAPDF
	
main: main.o deriv_pdf.o monte_carlo.o parameters.o k_factors_dy.o k_factors_higgs.o polygamma.o k_factors_nnlo_dy.o
	g++ -o DY_num main.o deriv_pdf.o monte_carlo.o parameters.o k_factors_dy.o k_factors_higgs.o polygamma.o k_factors_nnlo_dy.o -Wall -Wextra -lgsl -lgslcblas -lm -lcuba -L/home/mbeekveld/Cuba/lib  -Wreturn-local-addr -L/home/mbeekveld/LHAPDF/lib -lLHAPDF

deriv_pdf.o: deriv_pdf.cpp
	g++ -c deriv_pdf.cpp -I/home/mbeekveld/Cuba/include  -Wreturn-local-addr -O3 -I/home/mbeekveld/LHAPDF/include -std=c++11

monte_carlo.o: monte_carlo.cpp
	g++ -c monte_carlo.cpp -Wreturn-local-addr -O3 -I/home/mbeekveld/LHAPDF/include -std=c++11


polygamma.o: polygamma.cpp
	g++ -c polygamma.cpp -Wreturn-local-addr -O3 -std=c++11


parameters.o: parameters.cpp
	g++ -c parameters.cpp -I/home/mbeekveld/Cuba/include  -Wreturn-local-addr -O3 -I/home/mbeekveld/LHAPDF/include -std=c++11

	
k_factors_nnlo_dy.o: k_factors_nnlo_dy.cpp
	g++ -c k_factors_nnlo_dy.cpp -I/home/mbeekveld/Cuba/include  -Wreturn-local-addr -O3 -I/home/mbeekveld/LHAPDF/include -std=c++11

	
k_factors_dy.o: k_factors_dy.cpp
	g++ -c k_factors_dy.cpp -I/home/mbeekveld/Cuba/include  -Wreturn-local-addr -O3 -I/home/mbeekveld/LHAPDF/include -std=c++11
	
k_factors_higgs.o: k_factors_higgs.cpp
	g++ -c k_factors_higgs.cpp -I/home/mbeekveld/Cuba/include  -Wreturn-local-addr -O3 -I/home/mbeekveld/LHAPDF/include -std=c++11

main.o: main.cpp
	g++ -c main.cpp -I/home/mbeekveld/Cuba/include  -Wreturn-local-addr -O3 -I/home/mbeekveld/LHAPDF/include 

#-Wall -Wextra  voor als je echt niet weet wat er mis gaat, misschien lost het op
