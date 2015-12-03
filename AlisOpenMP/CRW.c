//============================================================================
// Name        : DataRetriever.cpp
// Author      : Utkarsh
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <ifstream.h>
#include <string>

using namespace std;

#ifdef WINDOWS
	#include<direct.h>
	#define getCWD _getcwd
#else
	#include<unistd.h>
	#define getCWD getcwd
#endif

int main()
{
	char wd[1024]="";
	DIR *dirp;
	struct dirent *dp;
	std::ifstream fp;

								//cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	if (NULL==getcwd(wd,1024))
		cout<<"Error opening curent path"<<getcwd(wd,1024);
	else{cout<<"Current path is"<< wd;}

	cout<<" Enter data source path:"<<endl;
	cin>>wd;
								//cout<<wd;
	if (NULL==wd) return -1;

int a=0;
dirp=opendir(wd);
std::vector<long> row(17700,0);
std::vector<long> col(177000,0);
std::vector<bool> rating(177000,0);
std::string line;
	while (dirp!=NULL && (dp = readdir(dirp)) != NULL)
	{
		a++;
		if (!strcmp (dp->d_name, "."))
		            continue;
		if (!strcmp (dp->d_name, ".."))
		            continue;

//		cout<<dp->d_name;
	//	if (a==50) return 0;

		fp.open(dp->d_name);

		char output[100];
		char unfilteredrow[20];
		while (std::getline(fp, line))
		{
		    std::istringstream iss(line);

		    double row;
		    char rating[2];
		    char date[20];

		    if (!(iss >> row >> rating >> date)) { break; } // error
		    storeCRW(row,col,rating);
		    // process pair (a,b)
		}

	}
	myfile.close();

	return 1;
}

