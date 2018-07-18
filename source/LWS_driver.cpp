#include "simulation.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;


int main(int argc, char** argv)
{
	if (argc < 2 || argc % 2 != 0)
	{
		cout << "Usage: ./simulation inputFile [outputVar outputFile] ..." << endl;
		cout <<	"outputVar options: pos" << endl;
		throw invalid_argument("Check number of arguments.");
	}


	char* inputFile = argv[1];
	Simulation simul(inputFile);


	simul.run();


	int argc_read = 2;
	while (argc_read < argc)
	{
		string outputVar = argv[argc_read];
		char* outputFilename = argv[argc_read+1];

		if (outputVar  == "pos")
		{
			ofstream outFile;
			outFile.open(outputFilename);
			simul.print_pos(outFile, 100);
			outFile.close();

		} else {
			cout << "Usage: ./simulation inputFile [outputVar outputFile] ..." << endl;
			cout <<	"outputVar options: pos" << endl;
			throw invalid_argument("Check outputVar names");
		}

		argc_read += 2;
	}

	return 0;
}