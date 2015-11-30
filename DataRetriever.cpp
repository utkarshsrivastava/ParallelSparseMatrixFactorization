//============================================================================
// Name        : DataRetriever.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <string>
//#include <sys/types>
#include <sstream>
#include <dirent.h>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
using namespace std;

#ifdef WINDOWS
	#include<direct.h>
	#define getCWD _getcwd
#else
	#include<unistd.h>
	#define getCWD getcwd
#endif

std::vector<long> rows(17700,0);
std::vector<long> cols(100000000,0);
std::vector<bool> ratings(100000000,0);
std::string line="";
ifstream fp;
ifstream rowptrfile;

int main();
void CRW_colm_rate(int row_num, double col,bool rating);

 int main()
{
	char wd[1024]="";
	DIR *dirp;
	struct dirent *dp;

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

char output[100];
char unfilteredrow[20];

//update row_ptr////////////
bool rowptr_done;
//if (this is UNIX)
rowptr_done=system("rm abcd.txt;wc -l *.txt> abcd.txt;"
		" sed -i 's/mv_//g' abcd.txt; sed -i 's/.txt//g' abcd.txt;"
		" sed -i '$ d'  abcd.txt");

rowptrfile. open("abcd.txt");
if (!rowptr_done){return -1;}

int total_lines,rownum=0;
while (std::getline(rowptrfile, line))
{
    std::istringstream iss(line);
    if (!(iss>>total_lines>>rownum)){break;}
    rows[rownum]=total_lines-1;
}
long tmp=0;
long rowitr =1;
long offset=1;
while(rowitr<rows.size())
{
tmp=rows[rowitr];
rows[rowitr]=offset+rows[rowitr-1];
offset=tmp;
rowitr++;
}
rows[0]=1;

//column update now///////////////////////////

int col_num;
int rating[2];
char date[20];
bool isRow;

	while (dirp!=NULL && (dp = readdir(dirp)) != NULL)
	{
		a++;
		if (! (dp->d_name.compare(".")))
		            continue;
		if (! (dp->d_name.compare("..")))
		            continue;

//		cout<<dp->d_name;
	//	if (a==50) return 0;

		fp.open(dp->d_name);
		std::string rowunfiltered;
	    isRow=1;
	    col_num=0;
		while (std::getline(fp, line))
		{
		    std::istringstream iss(line);
		    if(isRow)
		    {
		    	isRow=false;
		    	iss >> rowunfiltered;
		    	rowunfiltered[rowunfiltered.size()-1]='\0';
		    	row_num=rowunfiltered; //convert char to int
		    	continue;
		    }
		    if (!(iss >> col >> rating >> date)) { break; } // error
		    col_num++;
		    CRW_colm_rate(rownum,col,col_num,rating/4);
		    // process pair (a,b)
		}


	}
	myfile.close();

	return 1;
}

void CRW_colm_rate(int row_num, double col,double col_num,bool rating)
{
cols[rows[row_num]+col_num]=col;
ratings[rows[row_num]+col_num]=rating;
}

