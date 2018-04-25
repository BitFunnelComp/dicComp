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
#include<set>
#include<nmmintrin.h>
#include <sstream>  
#include <string> 
#include<sys/time.h>  
#include<climits>
#include <chrono>
#include <ctime>
using namespace std;

typedef pair<unsigned, unsigned >Rank_Row;
map<string, unsigned>termindoc;
const unsigned shard_count = 9;
const unsigned rank_count = 6;
map<string, vector<Rank_Row>>term_to_rows;

unsigned long **matrix[rank_count];
unsigned Row[shard_count][rank_count];
unsigned Col[shard_count][rank_count];
unsigned myshardid;
unsigned query_term_count = 0;
unsigned long *midresult[rank_count];
unsigned long *result;
unsigned resultcount = 0;
vector<string>query_term(1000000);
vector<string>term;
unsigned rankflag[rank_count];
std::chrono::time_point<std::chrono::system_clock> Start;
std::chrono::time_point<std::chrono::system_clock> End;
double mytime = 0;
unsigned *density[rank_count];
unsigned denrankid = 0;
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
			{
				Rank_Row rr;
				rr.first = rankid; rr.second = atoi(field.c_str());
				term_to_rows[term].push_back(rr);
			}
			count++;
		}
	}
	fin.close();
}
void read_file(unsigned rankid)
{
	string filename = "";
	FILE *file = fopen(filename.c_str(), "rb");
	unsigned rows = Row[myshardid][rankid], cols = Col[myshardid][rankid] / 64;
	
	for (unsigned r = 0; r < rows; r++)
	{
		for (unsigned c = 0; c < cols; c++)
		{
			fread(&matrix[rankid][r][c], sizeof(unsigned long), 1, file);
			density[rankid][r] += _mm_popcnt_u64(matrix[rankid][r][c]);
		}
	}
	fclose(file);
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
int find_rowid(vector<unsigned>v, unsigned value)
{
	for (unsigned i = 0; i < v.size(); i++)
	{
		if (v[i] == value)
			return i;
	}
	return -1;
}
bool cmpdensity(unsigned &a, unsigned &b)
{
	if (density[denrankid][a] < density[denrankid][b])
		return true;
	return false;
}
void do_query(unsigned qid)
{
	resultcount = 0;
	fill(rankflag, rankflag + rank_count, 0);
	vector<unsigned>rowset[rank_count];
	for (unsigned i = 0; i < term.size(); ++i)
	{
		if (term_to_rows.find(term[i]) == term_to_rows.end())return;
		Rank_Row rr;
		for (unsigned j = 0; j < term_to_rows[term[i]].size(); ++j)
		{
			rr.first = term_to_rows[term[i]][j].first; rr.second = term_to_rows[term[i]][j].second;
			if (find_rowid(rowset[rr.first], rr.second) == -1)
			{
				rowset[rr.first].push_back(rr.second);
			}
		}

	}
	for (unsigned denrankid = 0; denrankid < rank_count; ++denrankid)
	{
		if (rowset[denrankid].size() != 0)
		{
			rankflag[denrankid] = 1;
		}
	}
	unsigned long rowid = 0, fixone = ULLONG_MAX, docinbit = 0, shortr = 0, temp = 0;
	unsigned col64 = 0,col64_2=0;
  Start = std::chrono::system_clock::now();
	for (int r = rank_count - 1; r >= 0; --r)
	{
		if (rankflag[r] == 0)
			;
		else
		{
			col64 = (Col[myshardid][r] >> 6);
			for (unsigned i = 0; i < col64; ++i)
			{
				temp = matrix[r][rowset[r][0]][i];
				for (unsigned j = 1; j < rowset[r].size(); ++j)
				{
					rowid = rowset[r][j];
					temp &= matrix[r][rowid][i];
				}
				midresult[r][i] = temp;
			}
		}
	}
	col64_2 = 0;
	col64 = 0;
	for (int r = rank_count - 1; r >= 0; --r)
	{
		if (rankflag[r] == 0)
			;
		else
		{
			if (col64 == 0)
			{
				col64 = (Col[myshardid][r] >> 6);
				memcpy(result, midresult[r], col64*sizeof(unsigned long));
			}
			else
			{
				col64_2 = (Col[myshardid][r] >> 6);
				for (unsigned i = 0; i < col64_2; i += col64)
				{
					for (unsigned j = i, t = 0; t< col64; ++j, ++t)
					{
						midresult[r][j] &= result[t];
					}
				}
				swap(midresult[r], result);
				col64 = col64_2;
			}
		}
	}
	col64_2 = (Col[myshardid][0] >> 6);
	
	for (unsigned i = 0; i < col64; ++i)
	{
		resultcount += _mm_popcnt_u64(result[i]);
	}
	resultcount *= (col64_2 / col64);
  End = std::chrono::system_clock::now();
 std::chrono::duration<double> elapsed_seconds = End-Start;
	mytime += elapsed_seconds.count()*1000;

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
			cout << tmpo << " " << tmpn << endl;
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
	cout << tmpo << " " << tmpn << " " << error << endl;
	fino.close();
	finn.close();
	cout << "query count=" << qcount << endl;
	cout << "shard " << myshardid << "\terror rate=" << error_rate / qcount << endl;
}
void init_data()
{
	Row[0][0] = 1310 - 3; Row[0][1] = 1000; Row[0][2] = 1027; Row[0][3] = 1000; Row[0][4] = 1077; Row[0][5] = 1010;
	Row[1][0] = 1616 - 3; Row[1][1] = 1000; Row[1][2] = 1095; Row[1][3] = 1000; Row[1][4] = 1232; Row[1][5] = 1961;
	Row[2][0] = 1995 - 3; Row[2][1] = 1000; Row[2][2] = 1154; Row[2][3] = 1000; Row[2][4] = 1460; Row[2][5] = 3560;
	Row[3][0] = 2967 - 3; Row[3][1] = 1000; Row[3][2] = 1667; Row[3][3] = 2107; Row[3][4] = 1801; Row[3][5] = 5701;
	Row[4][0] = 4752 - 3; Row[4][1] = 1179; Row[4][2] = 1497; Row[4][3] = 1502; Row[4][4] = 2318; Row[4][5] = 8965;
	Row[5][0] = 7659 - 3; Row[5][1] = 2014; Row[5][2] = 2403; Row[5][3] = 2632; Row[5][4] = 4390; Row[5][5] = 13666;
	Row[6][0] = 15820 - 3; Row[6][1] = 4696; Row[6][2] = 6297; Row[6][3] = 7252; Row[6][4] = 9454; Row[6][5] = 32485;
	Row[7][0] = 29314 - 3; Row[7][1] = 9022; Row[7][2] = 14588; Row[7][3] = 16300; Row[7][4] = 19810; Row[7][5] = 36242;
	Row[8][0] = 69714 - 3; Row[8][1] = 34880; Row[8][2] = 52364; Row[8][3] = 74024; Row[8][4] = 12052; Row[8][5] = 98188;

	Col[0][0] = 2647295;
	Col[1][0] = 2614477;
	Col[2][0] = 5519536;
	Col[3][0] = 7373033;
	Col[4][0] = 4503823;
	Col[5][0] = 1872882;
	Col[6][0] = 497135;
	Col[7][0] = 157717;
	Col[8][0] = 19285;
	cout << "shard id=" << myshardid << endl;
	Col[myshardid][5] = (Col[myshardid][0] + pow(2, 5) - 1) / pow(2, 5);
	Col[myshardid][5] = (Col[myshardid][5] + 63) / 64 * 64;
	for (unsigned r = 0; r < rank_count; r++)
	{
		cout << "rank " << r << endl;
		Col[myshardid][r] = Col[myshardid][5] * pow(2, 5 - r);
		cout << "maxrow=" << Row[myshardid][r] << "\tmaxcol(ali)=" << Col[myshardid][r] << endl;
	}

	for (unsigned r = 0; r < rank_count; r++)
	{
		matrix[r] = new unsigned long*[Row[myshardid][r]];
		density[r] = new unsigned [Row[myshardid][r]];
		fill(density[r], density[r] + Row[myshardid][r], 0);
		for (unsigned row = 0; row < Row[myshardid][r]; row++)
		{
			matrix[r][row] = new unsigned long[Col[myshardid][r] / 64];
			fill(matrix[r][row], matrix[r][row] + Col[myshardid][r] / 64, 0);
		}
	}
	result = new unsigned long[25205183 / 64];
	for (unsigned i = 0; i<rank_count; i++)
		midresult[i] = new unsigned long[25205183 / 64];
}
int main(int argc, char *argv[])
{
	cpu_set_t mask;
	CPU_ZERO(&mask);
	CPU_SET(0, &mask);
	myshardid = atoi(argv[1]);
	init_data();
	read_term_to_rows();
	cout << "read term to rows over" << endl;
	read_query();
	cout << "read query over" << endl;
	cout << "query count=" << query_term_count << endl;
	for (unsigned r = 0; r < rank_count; r++)
	read_file(r);
	ofstream fout("");
	for (unsigned i = 0; i <query_term_count; i++)
	{
		parse_query(i);
		do_query(i);
		
		fout << resultcount << endl;
	}
	fout.close();
	cout << "test match time=" << mytime << "ms\ttime per query=" << mytime / query_term_count << "ms" << endl;
	get_error_rate();
   mytime=0;
   for(unsigned turn=0;turn <10;turn++)
 { for (unsigned i = 0; i <query_term_count; i++)
  	{
  		parse_query(i);
  		do_query(i);
  	}
  }
  	cout << "match time=" << mytime << "ms\ttime per query=" << mytime / query_term_count/10 << "ms" << endl;
   cout << "**********************************************************************" << endl;
}
