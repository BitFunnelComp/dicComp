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
#include <sstream>  
#include <string> 
#include<sys/time.h> 
#include <climits>
#include <xmmintrin.h>
#include <nmmintrin.h>
#include <chrono>
#include <ctime>
using namespace std;
unsigned myshardid = 0;
const unsigned shardcount = 9;
unsigned long **bit;
unsigned maxrow[shardcount];
unsigned maxcolold[shardcount];
unsigned col64;
map<string, vector<unsigned>>term_to_rows;
unsigned query_term_count = 0;
vector<string>query_term(1000000);
vector<string>term;
unsigned result_count = 0;
unsigned long *result;
std::chrono::time_point<std::chrono::system_clock> Start;
std::chrono::time_point<std::chrono::system_clock> End;
double time_baseline = 0;
unsigned long fixbit[64];
unsigned *density;
void readfile()
{
	string filename = "";
	FILE *file1 = fopen(filename.c_str(), "rb");
	unsigned colbyte = (maxcolold[myshardid] + 7) / 8;
	for (unsigned i = 0; i < maxrow[myshardid]; i++)
	{
		for (unsigned j = 0; j < colbyte; j += 8)
		{
			unsigned char temp[8];
			fill(temp, temp + 8, 0);
			for (unsigned t = 0, k = j; t < 8 && k < j + 8 && k < colbyte; t++, k++)
			{
				fread(&temp[t], sizeof(unsigned char), 1, file1);
			}
			for (unsigned t = 0; t < 8; t++)
			{
				bit[i][j / 8] <<= 8;
				bit[i][j / 8] |= temp[t];
			}
      density[i] += _mm_popcnt_u64(bit[i][j / 8]);
		}
	}
	fclose(file1);
}
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
bool cmpdensity(unsigned &a, unsigned &b)
{
	if (density[a] < density[b])
		return true;
	return false;
}
void do_query( unsigned qid)
{
   result_count = 0;
	vector<unsigned>rowset;
	for (unsigned i = 0; i < term.size(); i++)
	{
		if (term_to_rows.find(term[i]) == term_to_rows.end()){ return; }
		for (unsigned j = 0; j < term_to_rows[term[i]].size(); j++)
		{
			if (find_rowid(rowset, term_to_rows[term[i]][j]) == -1)
			{
				rowset.push_back(term_to_rows[term[i]][j]);
			}
		}
	}
 
	unsigned long fixnum = ULLONG_MAX,temp=0;
	fill(result, result + col64, fixnum);
	unsigned rowid = 0;

Start = std::chrono::system_clock::now();
 for (unsigned j = 0; j < col64; ++j)
 {
   temp=fixnum;
   for (unsigned i = 0; i < rowset.size(); ++i)
   {
     rowid = rowset[i];
     temp &= bit[rowid][j];
   }
   result[j]=temp;
 }
	for (unsigned i = 0; i < col64; ++i)
		result_count += _mm_popcnt_u64(result[i]);
  End = std::chrono::system_clock::now();
 std::chrono::duration<double> elapsed_seconds = End-Start;
  time_baseline += elapsed_seconds.count()*1000;
}
void init_data()
{
	for (unsigned i = 0; i<64; i++)
	{
		fixbit[i] = (unsigned long)pow(2, 63 - i);
	}
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
 

	unsigned col_max = 0;
	density= new unsigned [maxrow[myshardid]];
	fill(density, density + maxrow[myshardid], 0);
	
	col64 = ((maxcolold[myshardid] + 7)/8+7) / 8;
	col_max = max(col_max, col64);
	bit = new unsigned long *[maxrow[myshardid]];
	for (unsigned row = 0; row < maxrow[myshardid]; row++)
	bit[row] = new unsigned long[col64];
	
	result = new unsigned long[col_max];
}
void get_error_rate()
{
	ifstream fino("");
	ifstream finn("");
	string stmpo, stmpn;
	unsigned tmpo, tmpn, qcount = 0;
	double error_rate = 0, error = 0;
	while (getline(fino, stmpo) && getline(finn, stmpn))
	{
		tmpo = atoi(stmpo.c_str());
		tmpn = atoi(stmpn.c_str());
		if (tmpo > tmpn)
		{
			cout << "wrong " << qcount << endl;
      return;
		}
		error = 0;
		if (tmpo == 0 && tmpn != 0)error = 1;
		else if (tmpo == 0 && tmpn == 0)error = 0;
		else error = (double)(tmpn - tmpo) / tmpn;
		if (error > 1)
			error = 1;
		error_rate += error;
		qcount++;
	}
	fino.close();
	finn.close();
	cout << "query count=" << qcount << endl;
	cout << "shard " << myshardid << "\terror rate=" << error_rate / qcount << endl;
}
int main(int argc, char *argv[])
{
  cpu_set_t mask;
	CPU_ZERO(&mask);
	CPU_SET(0, &mask);
	myshardid = atoi(argv[1]);
	init_data();
	readfile();
	cout << "read file over" << endl;
	read_term_to_rows();
	cout << "read_term_to_rows over" << endl;
	read_query(); cout << "read query over querycount=" << query_term_count << endl;

	ofstream fout("");
	for (unsigned i = 0; i <query_term_count; i++)
	{
		
		parse_query(i);
		do_query(i);
		fout << result_count << endl;
	}
	fout.close();
	cout << "shard id=" << myshardid << endl;
	cout << "all match time test =" << time_baseline << "ms\t" << "time per query=" << time_baseline / query_term_count << "ms" << endl;
 get_error_rate();
 time_baseline=0;
   for(unsigned turn=0;turn <10;turn++)
 { for (unsigned i = 0; i <query_term_count; i++)
  	{
  		parse_query(i);
  		do_query(i);
  	}
  }
  cout << "match time=" <<time_baseline  << "ms\ttime per query=" <<time_baseline/ query_term_count/10 << "ms" << endl;
 cout << "**********************************************************************" << endl;
}

