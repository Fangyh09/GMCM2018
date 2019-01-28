LIBPATH = -L/home/yinghong/tmp -L/usr/local/matlab/2016b/runtime/glnxa64 -lmwmclmcrrt -lMyLinprog
INCLUDEPATH = -I/usr/local/matlab/2016b/extern/include
 

LD_LIBRARY_PATH = /home/yinghong/tmp:/usr/local/matlab/2016b/runtime/glnxa64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
 
MainSin.o:main.cpp
	g++ -std=c++11 -c main.cpp $(INCLUDEPATH)
Genetic.o: Genetic.cpp
	g++ -std=c++11 -c Genetic.cpp $(INCLUDEPATH)
MainSinApp:main.o
	g++ -std=c++11   -o MainSinApp main.o Genetic.o $(LIBPATH)
	./MainSinApp
 
clean:
	rm -f *.o
