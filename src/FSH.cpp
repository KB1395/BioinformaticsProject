//============================================================================
// Name        : HashTest.cpp
// Author      : Samuele Girotto
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "FSH/FSH.h"

int main(int argc, char* argv[]) {
	string dir_output = "output/";
	FileParameter param;

	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i], "-si") == 0)
		{
			i++;
			if(!param.init(argv[i], ""))
			{
				cerr<<endl<<"Please enter an input filename single-end: -si <AbsPathFile>"<<flush;
				return 0;
			}
		}
		else if(strcmp(argv[i], "-pi") == 0)
		{
			i++;
			if(!param.init(argv[i], argv[i+1]))
			{
				cerr<<endl<<"Please enter an input filenames paired-end: -pi <AbsPathFile1> <AbsPathFile2>"<<flush;
				return 0;
			}
		}
		else if(strcmp(argv[i], "-q") == 0)
		{
			i++;
			regex rgx("^1(0|1)*1$");
			string spaced_qmer(argv[i]);
			if(regex_match(spaced_qmer, rgx))
			{
				cerr<<endl<<"Please enter -q parameter with 1 at begin and end of the string. "
						"Ex. 1**1*11*1. 1 is the simbol considered, any others are like *"<<flush;
				return 0;
			}
			else
				param.addSpacedQmer(spaced_qmer, spaced_qmer);
		}
		else if(strcmp(argv[i], "-dirO") == 0)
		{
			i++;
			if(argv[i] == "")
			{
				cerr<<endl<<"Enter valid path for output directory"<<flush;
				return 0;
			}
			else
				dir_output = argv[i];
		}
	}
	//Creo cartella output se non presente
	createDirAndSubDir(dir_output);

	bool single_test_equals = false;
	string dir_output_1 = dir_output + "single/";
	FSH test;
	if(test.load_sequences(param))
	{
		for(size_t j = 0; j < param.spaced_seed.size(); ++j)
		{
			string dir_output_tmp = dir_output_1 + param.spaced_seed[j].first +  param.spaced_seed[j].second.toString() + "/";

			test.single_run(param.spaced_seed[j].second, single_test_equals);
			test.single_save(param, dir_output_tmp);
		}
	}

	////////////////////////////////////////////////////////////////7

	bool multi_test_equals = false;
	string dir_output_2 = dir_output + "multi/";
	vector<SpacedQmer> multi_spaced;
	for(size_t i = 0; i < param.spaced_seed.size(); ++i)
		multi_spaced.push_back(param.spaced_seed[i].second);

	FSH test;
	if(test.load_sequences(param))
	{
		test.multi_run(multi_spaced, multi_test_equals);
		test.multi_save(param, multi_spaced, dir_output_2);
	}
}
