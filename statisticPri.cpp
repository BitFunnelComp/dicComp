#include<iostream>
#include<fstream>
#include<string>
#include<cstdio>
#include<cstdint>
#include <algorithm>
#include<math.h>
#include<string.h>
#include<vector>
#include<map>
#include<nmmintrin.h>

#include <sstream>  
#include <string> 

using namespace std;
const unsigned blocksize = 32;
const unsigned blockbyte = blocksize / 8;
const unsigned mapsize = 16;
const unsigned block_map_size = 65536;
const unsigned mapbyte = mapsize / 8;
const unsigned shardcount = 9;

unsigned maxrow[shardcount];
unsigned maxcolold[shardcount];
unsigned maxcol[shardcount];

int result_count = 0;
map<string, vector<unsigned>>term_to_rows;
unsigned termcount = 0;
unsigned long long sum_short = 0;
unsigned query_term_count = 0;
vector<string>query_term(1000050);
vector<string>term;
unsigned myshardid = 0;


typedef pair<unsigned, unsigned long>FRE;
map<unsigned, unsigned long>row_fre;
vector<FRE>vfre;
double thresh;
unsigned char*rowflag;
void read_term_to_rows()
{
	string filename = "";
	ifstream fin(filename);
	string line = "", term;
	unsigned count = 0;
	while (getline(fin, line))  
	{
		count = 0;
		istringstream sin(line);
		string field;
		while (getline(sin, field, ',')) 
		{  
			if (count == 0)
				term = field;
			else if (count % 3 == 0)
				term_to_rows[term].push_back(atoi(field.c_str()));
			count++;
		}
	}
	fin.close();
	cout << "shard " << myshardid << " has " << term_to_rows.size() << " terms" << endl;
}

void read_query()
{
	ifstream fin("");
	
	char temp_line[2000];
	string str, term;
	int pos = 0;
	unsigned linecount = 0;
	while (fin.getline(temp_line, 2000))
	{
		str = temp_line;
		query_term[linecount++] = str;
	}
	query_term_count = linecount;
	fin.close();
}

int find_rowid(vector<unsigned>v, unsigned value)
{
	for (unsigned i = 0; i < v.size(); i++)
	{
		if (v[i] == value)
			return i;
	}
	return -1;
}

void parse_query(unsigned qid)
{
	term.clear();
	string str = query_term[qid];
	int pos = str.find(" ");
	string sterm;
	while (pos != -1)
	{
		sterm = str.substr(0, pos);
		term.push_back(sterm);
		str = str.substr(pos + 1, str.length() - pos - 1);
		pos = str.find(" ");
	}
 term.push_back(str);
}
bool cmp(FRE &b1, FRE &b2)
{
	if (b1.second> b2.second)
		return true;
	return false;
}

void statistic()
{
	fill(rowflag, rowflag + maxrow[myshardid] + 10, 1);
	for (unsigned q = 0; q < query_term_count; q++)
	{
		parse_query(q);
		vector<unsigned>rowset;
		for (unsigned i = 0; i < term.size(); i++)
		{
			if (term_to_rows.find(term[i]) == term_to_rows.end()){ continue; }
			for (unsigned j = 0; j < term_to_rows[term[i]].size(); j++)
			{
				if (find_rowid(rowset, term_to_rows[term[i]][j]) == -1)
				{
					rowset.push_back(term_to_rows[term[i]][j]); 
				}
			}
		}
		for (unsigned i = 0; i < rowset.size(); i++)
		{
			row_fre[rowset[i]]++;
		}
	}
	vfre = vector<FRE>(row_fre.begin(), row_fre.end());
	sort(vfre.begin(), vfre.end(), cmp);
	unsigned long sumfre = 0, tmpfre = 0;
	for (unsigned i = 0; i < vfre.size(); i++)
	{
		sumfre += vfre[i].second;
	}
	cout << "sumfre=" << sumfre << endl;
	sumfre *= thresh;
	for (unsigned i = 0; i < vfre.size(); i++)
	{
		tmpfre += vfre[i].second;
		if (tmpfre>sumfre)
		{
			cout << "T rows count=" << i << endl;
			break;
		}
		rowflag[vfre[i].first] = 0;
	}
 cout<<"all rows "<<vfre.size()<<endl;
}

void init_data()
{
	
	maxrow[0] = 1032 - 3;
	maxrow[1] = 1269 - 3;
	maxrow[2] = 2357 - 3;
	maxrow[3] = 4480 - 3;
	maxrow[4] = 6510 - 3;
	maxrow[5] =10719 - 3;
	maxrow[6] = 23330 - 3;
	maxrow[7] = 43666 - 3;
	maxrow[8] = 117399 - 3;

	maxcolold[0] = 2647295;
	maxcolold[1] = 2614477;
	maxcolold[2] = 5519536;
	maxcolold[3] = 7373033;
	maxcolold[4] = 4503823;
	maxcolold[5] = 1872882;
	maxcolold[6] = 497135;
	maxcolold[7] = 157717;
	maxcolold[8] = 19285;
 
	unsigned col_max = 0, colblock = 0;

	maxcol[myshardid] = ((maxcolold[myshardid] + blocksize - 1) / blocksize)*blocksize;
	colblock = maxcol[myshardid] / blocksize;
	col_max = max(col_max, colblock);
	
	rowflag = new unsigned char[maxrow[myshardid] + 10];
}
void write_file()
{
	string filename = "";
	FILE *filem = fopen(filename.c_str(), "wb");
	unsigned char temp;
	for (unsigned j = 0; j < maxrow[myshardid]; j += 8)
	{
		temp = 0;
		for (unsigned t = j, k = 128; t < j + 8 && t < maxrow[myshardid]; t++, k /= 2)
			temp += k*rowflag[t];
	
		fwrite(&temp, sizeof(unsigned char), 1, filem);
	}
	fclose(filem);
}
int main(int argc, char *argv[])
{

	myshardid = atoi(argv[1]);
	thresh = atof(argv[2]);
	init_data();
	cout << "init over" << endl;
	read_term_to_rows();
	cout << "read_term_to_rows over" << endl;
	read_query();
	cout << "read query over querycount=" << query_term_count << endl;
	statistic();
	write_file();
	cout << "over" << endl;

}
