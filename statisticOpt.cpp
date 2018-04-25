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

const unsigned shardcount = 9;
const unsigned rankcount = 6;
unsigned maxrow[shardcount][rankcount];
map<string, vector<unsigned>>term_to_rows[rankcount];
unsigned termcount = 0;
unsigned query_term_count = 0;
vector<string>query_term(1000000);
vector<string>term;
unsigned myshardid = 0;


typedef pair<unsigned, unsigned long>FRE;
map<unsigned, unsigned long>row_fre[rankcount];
vector<FRE>vfre;
double thresh;
unsigned char*rowflag[rankcount];
void read_term_to_rows()
{
	string filename = "";
	ifstream fin(filename);
	string line = "", term;
	unsigned count = 0, rankid = 0;
	while (getline(fin, line))   
	{
		count = 0;
		rankid = 0;
		istringstream sin(line); 
		string field;
		while (getline(sin, field, ','))
		{ 
			if (count == 0)
				term = field;
			else if ((count + 1) % 3 == 0)
				rankid = field[field.length() - 1] - '0';
			else if (count % 3 == 0)
				term_to_rows[rankid][term].push_back(atoi(field.c_str()));
			count++;
		}
	}
	fin.close();
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

void statistic(unsigned rankid)
{
  vfre.clear();
	fill(rowflag[rankid], rowflag[rankid] + maxrow[myshardid][rankid] + 10, 1);
	for (unsigned q = 0; q < query_term_count; q++)
	{
		parse_query(q);
		vector<unsigned>rowset;
		for (unsigned i = 0; i < term.size(); i++)
		{
			if (term_to_rows[rankid].find(term[i]) == term_to_rows[rankid].end()){ continue; }
			for (unsigned j = 0; j < term_to_rows[rankid][term[i]].size(); j++)
			{
				if (find_rowid(rowset, term_to_rows[rankid][term[i]][j]) == -1)
				{
					rowset.push_back(term_to_rows[rankid][term[i]][j]); 
				}
			}
		}
		for (unsigned i = 0; i < rowset.size(); i++)
		{
			row_fre[rankid][rowset[i]]++;
		}
	}
	vfre = vector<FRE>(row_fre[rankid].begin(), row_fre[rankid].end());
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
		rowflag[rankid][vfre[i].first] = 0;
	}
   cout<<"all visited rows "<<vfre.size()<<endl;
	cout << "all rows " << maxrow[myshardid][rankid] << endl;
}

void init_data()
{

	maxrow[0][0] = 1310-3; maxrow[0][1] = 1000; maxrow[0][2] = 1027; maxrow[0][3] = 1000; maxrow[0][4] = 1077; maxrow[0][5] = 1010;
	maxrow[1][0] = 1616-3; maxrow[1][1] = 1000; maxrow[1][2] = 1095; maxrow[1][3] = 1000; maxrow[1][4] = 1232; maxrow[1][5] = 1961;
	maxrow[2][0] = 1995-3; maxrow[2][1] = 1000; maxrow[2][2] = 1154; maxrow[2][3] = 1000; maxrow[2][4] = 1460; maxrow[2][5] = 3560;
	maxrow[3][0] = 2967-3; maxrow[3][1] = 1000; maxrow[3][2] = 1667; maxrow[3][3] = 2107; maxrow[3][4] = 1801; maxrow[3][5] = 5701;
	maxrow[4][0] = 4752-3; maxrow[4][1] = 1179; maxrow[4][2] = 1497; maxrow[4][3] = 1502; maxrow[4][4] = 2318; maxrow[4][5] = 8965;
	maxrow[5][0] = 7659-3; maxrow[5][1] = 2014; maxrow[5][2] = 2403; maxrow[5][3] = 2632; maxrow[5][4] = 4390; maxrow[5][5] = 13666;
	maxrow[6][0] = 15820-3; maxrow[6][1] = 4696; maxrow[6][2] = 6297; maxrow[6][3] = 7252; maxrow[6][4] = 9454; maxrow[6][5] = 32485;
	maxrow[7][0] = 29314-3; maxrow[7][1] = 9022; maxrow[7][2] = 14588; maxrow[7][3] = 16300; maxrow[7][4] = 19810; maxrow[7][5] = 36242;
	maxrow[8][0] = 69714-3; maxrow[8][1] = 34880; maxrow[8][2] = 52364; maxrow[8][3] = 74024; maxrow[8][4] = 12052; maxrow[8][5] = 98188;
	for (unsigned i = 0; i < rankcount;i++)
		rowflag[i] = new unsigned char[maxrow[myshardid][i] + 10];
}
void write_file(unsigned rankid)
{
	string filename = "";
	FILE *filem = fopen(filename.c_str(), "wb");
	unsigned char temp;
	for (unsigned j = 0; j < maxrow[myshardid][rankid]; j += 8)
	{
		temp = 0;
		for (unsigned t = j, k = 128; t < j + 8 && t < maxrow[myshardid][rankid]; t++, k /= 2)
			temp += k*rowflag[rankid][t];
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
	for (unsigned i = 0; i < rankcount; i++)
	{
		cout << "rank " << i << endl;
		statistic(i);
		write_file(i);
	}
	
	cout << "over" << endl;

}
