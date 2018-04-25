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
#include<sys/time.h>  
#include <sstream>  
#include <string> 
#include<climits>
#include <chrono>
#include <ctime>
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
unsigned pattern[block_map_size] __attribute__((aligned(64)));
unsigned short *ind2;
unsigned long *Row_offset;
unsigned *Row_size;

int result_count = 0;
unsigned  long *result;
unsigned  *midresult;
unsigned *result_map;

map<string, vector<unsigned>>term_to_rows;
unsigned termcount = 0;
unsigned long long sum_short = 0;
unsigned query_term_count = 0;
vector<string>query_term(1000000);
vector<string>term;
unsigned myshardid = 0;
unsigned long fixbit[64];
unsigned long andfix[4];
vector<unsigned>resultdoc(25205183);

std::chrono::time_point<std::chrono::system_clock> Start;
std::chrono::time_point<std::chrono::system_clock> End;
double time_notspread;
unsigned char *rowflag;
unsigned thresh = 0;
double rthresh = 0;
void read_rowflag()
{
	string filename = ""; 
	FILE *filem = fopen(filename.c_str(), "rb");
	unsigned char temp;
	unsigned u_temp = 0;
	unsigned T = 0, NT = 0;
	for (unsigned i = 0; i <maxrow[myshardid]; i += 8)
	{
		fread(&temp, sizeof(unsigned char), 1, filem);
		u_temp = (unsigned)(unsigned char)temp;
		for (unsigned t = i, k = 128; t < i + 8 && t <maxrow[myshardid]; t++, k /= 2)
		{
			rowflag[t] = u_temp / k;
			u_temp %= k;
			if (rowflag[t] == 0)T++;
			else  NT++;
		}

	}
	fclose(filem);
	cout << "T rows count=" << T << endl;
	cout << "all rows count=" << T + NT << endl;
}
void read_pattern()
{
	string filenamep = "";
	FILE *filep = fopen(filenamep.c_str(), "rb");
	for (unsigned j = 0; j < block_map_size; j++)
	{
		fread(&pattern[j], sizeof(unsigned), 1, filep);
	}
	fclose(filep);
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
void read_ind2()
{
	sum_short = 0;
	string filename2 = "";
	FILE *file = fopen(filename2.c_str(), "rb");
	while (fread(&ind2[sum_short++], sizeof(unsigned short), 1, file));
	cout << "shard " << myshardid << " ind2=" << sum_short << " shorts data" << endl;

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
void get_row_offset_size()
{
	long long offset = 0;
	unsigned colblock = maxcol[myshardid] / blocksize;
	for (unsigned i = 0; i < maxrow[myshardid]; i++)
	{
		if (rowflag[i] == 0)
			Row_size[i] = colblock * 2;
		else
			Row_size[i] = colblock;
		Row_offset[i] = offset;
		offset += Row_size[i];
	}
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
void do_notspread(unsigned qid)
{
	result_count = 0;
	vector<unsigned>rowset0;
	vector<unsigned>rowset1;
	for (unsigned i = 0; i < term.size(); i++)
	{
		if (term_to_rows.find(term[i]) == term_to_rows.end()){ return; }
		for (unsigned j = 0; j < term_to_rows[term[i]].size(); j++)
		{
			unsigned rowid = term_to_rows[term[i]][j];
			if (rowflag[rowid] == 0)
			{
				if (find_rowid(rowset0, rowid) == -1)
				{
					rowset0.push_back(rowid); 
				}
			}
			else
			{
				if (find_rowid(rowset1, rowid) == -1)
				{
					rowset1.push_back(rowid);
				}
			}
		}
	}
	unsigned colblock = maxcol[myshardid] / blocksize;
	unsigned col64 = (colblock + 1) / 2;
	unsigned long rowid, offset;
	unsigned long fixnumb = (unsigned long)pow(2, 32) - 1, fixnum = ULLONG_MAX;
   unsigned long temp0,temp1,fixnum1=(unsigned long)pow(2,32)-1;
	fill(midresult, midresult + colblock, fixnumb);
	fill(result, result + col64, fixnum);

  Start = std::chrono::system_clock::now();
	for (unsigned i = 0; i < rowset1.size(); ++i)
	{
		rowid = rowset1[i]; offset = Row_offset[rowid];
		for (unsigned j = 0; j < colblock; ++j)
		{
			midresult[j] &= pattern[ind2[offset++]];
		}
	}
  /* for (unsigned j = 0; j < colblock; ++j)
	{
			temp1=fixnum1;
      for (unsigned i = 0; i < rowset1.size()&&temp1!=0; ++i)
      {
        rowid = rowset1[i]; offset = Row_offset[rowid]+j;
        temp1&=pattern[ind2[offset]];
      }
      midresult[j]=temp1;
	}*/
	unsigned long *tmp = NULL;
	for (unsigned i = 0; i < rowset0.size(); ++i)
	{
		rowid = rowset0[i]; offset = Row_offset[rowid];
		for (unsigned j = 0; j < col64; ++j, offset += 4)
		{
			tmp = reinterpret_cast<unsigned long *>(ind2 + offset);
			result[j] &= *tmp;
		}
	}
   /*for (unsigned j = 0; j < col64; ++j)
   {
       temp0=fixnum;
     	for (unsigned i = 0; i < rowset0.size()&&temp0!=0; ++i)
      {
          rowid = rowset0[i]; offset = Row_offset[rowid]+j*4;
          temp0&=*reinterpret_cast<unsigned long *>(ind2 + offset);
      }
      result[j]=temp0;
   }*/
	if (rowset1.size() != 0)
	{
		for (unsigned i = 0, j = 0; i < col64; ++i, j += 2)
		{
			tmp = reinterpret_cast<unsigned long *>(midresult + j);
			*tmp = (((*tmp&andfix[0]) >> 16) | ((*tmp&andfix[1]) << 16) | ((*tmp&andfix[2]) >> 16) | ((*tmp&andfix[3]) << 16));
			result[i] &= *tmp;
		}
	}
	for (unsigned i = 0; i < col64; ++i)
		result_count += _mm_popcnt_u64(result[i]);
 End = std::chrono::system_clock::now();
 std::chrono::duration<double> elapsed_seconds = End-Start;
 time_notspread += elapsed_seconds.count()*1000;
}
void init_data()
{
	for (unsigned i = 0; i<64; i++)
	{
		fixbit[i] = (unsigned long)pow(2,i);
	}
	andfix[0] = 0xffff000000000000;
	andfix[1] = 0x0000ffff00000000;
	andfix[2] = 0x00000000ffff0000;
	andfix[3] = 0x000000000000ffff;

	maxrow[0] = 1032 - 3;
	maxrow[1] = 1269 - 3;
	maxrow[2] = 2357 - 3;
	maxrow[3] = 4480 - 3;
	maxrow[4] = 6510 - 3;
	maxrow[5] = 10719 - 3;
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
	ind2 = new unsigned short[(unsigned long)maxrow[myshardid] * maxcol[myshardid] / 16];
	Row_offset = new unsigned long[maxrow[myshardid] + 10];
	Row_size = new unsigned[maxrow[myshardid] + 10];
	result = new unsigned long[col_max + 10];
	midresult = new unsigned[col_max + 10];
	result_map = new unsigned[col_max + 10];
	rowflag = new unsigned char[maxrow[myshardid] + 10];
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
	cout << tmpo << " " << tmpn << " " << error << endl;
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
	_mm_prefetch(pattern, _MM_HINT_T0);
	myshardid = atoi(argv[1]);
	
	rthresh = atof(argv[2]);
	init_data();
	cout << "init over" << endl;
	read_rowflag();
	cout << "read row flag over" << endl;
	read_pattern();
	cout << "read_pattern over" << endl;
	read_term_to_rows();
	cout << "read_term_to_rows over" << endl;
	read_ind2();
	cout << "read_ind2 over" << endl;
	read_query();
	cout << "read query over querycount=" << query_term_count << endl;
	get_row_offset_size();
	cout << "get_row_offset_size over" << endl;
 for (unsigned i = 0; i < 65536;i+=16)
			_mm_prefetch(pattern+i, _MM_HINT_T0);
	ofstream fout("");
	for (unsigned i = 0; i <query_term_count; i++)
	{
		parse_query(i);
		
		do_notspread(i);
	
		fout << result_count << endl;
	}
	fout.close();
	cout << "shard " << myshardid << " \tall match time=" << time_notspread << "ms\t" << " time per query=" << time_notspread / query_term_count << "ms" << endl;
	get_error_rate();
  time_notspread=0;
   for(unsigned turn=0;turn <10;turn++)
 { for (unsigned i = 0; i <query_term_count; i++)
  	{
  		parse_query(i);
  		do_notspread(i);
  	}
  }
  cout << "match time=" <<time_notspread  << "ms\ttime per query=" <<time_notspread/ query_term_count/10 << "ms" << endl;
 cout << "**********************************************************************" << endl;
}
