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
#include <climits>
#include <xmmintrin.h>
#include <nmmintrin.h>
using namespace std;
const unsigned blocksize = 32;
const unsigned blockbyte = blocksize / 8;
const unsigned mapsize = 16;
const unsigned block_map_size = 65536;
const unsigned mapbyte = mapsize / 8;
const unsigned shardcount = 9;
const unsigned alig = 32;
unsigned maxcolold[shardcount];
unsigned maxrow[shardcount];
unsigned maxcol[shardcount];
typedef pair<unsigned long, unsigned long> BLOCK;
unsigned char **Row;
unsigned *Rowsize;
unsigned **reduce_block_flag;

map<unsigned long, unsigned long>block_fre;
vector<BLOCK>vblock;
unsigned long T_count = 0;
unsigned long NT_count = 0;
unsigned long fix_count = 0;
unsigned long S_count = 0;
map<unsigned long, unsigned>pattern_to_id;
unsigned myshardid = 0;
double rthresh = 0;
unsigned char *rowflag;
void read_file()
{
	string filename = "";
	FILE *file1 = fopen(filename.c_str(), "rb");
	unsigned colbyte = maxcol[myshardid] / 8;
	for (unsigned i = 0; i < maxrow[myshardid]; i++)
	{
		for (unsigned j = 0; j < colbyte; j++)
		{
			fread(&Row[i][j], sizeof(unsigned char), 1, file1);
		}
		Rowsize[i] = colbyte;
	}
	fclose(file1);

}
void read_rowflag()
{
	string filename = "";
	FILE *filem = fopen(filename.c_str(), "rb");
	unsigned char temp;
	unsigned u_temp = 0;
	unsigned T = 0, NT = 0;
	for (unsigned i = 0; i <maxrow[myshardid]; i+=8)
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
	cout << "all rows count=" << T+NT << endl;
}
void get_block_frequency(unsigned rows, unsigned startcol)
{
	unsigned char temp[blockbyte];
	for (unsigned i = 0; i < blockbyte; i++)
		temp[i] = 0;
	for (unsigned i = startcol; i < startcol + blockbyte&&i<Rowsize[rows]; i++)
		temp[i - startcol] = Row[rows][i];
	unsigned long  ltemp = 0;
	for (unsigned i = 0; i < blockbyte; i++)
	{
		ltemp *= pow(2, 8);
		ltemp |= (unsigned)temp[i];
	}
	block_fre[ltemp]++;
}
bool cmp(BLOCK &b1, BLOCK &b2)
{
	if (b1.second> b2.second)
		return true;
	return false;
}
void get_block_map()
{
	unsigned ltemp = (unsigned long)pow(2, 32) - 1;
	block_fre[ltemp] = (unsigned long)pow(2, 60);
	for (unsigned i = 0; i < maxrow[myshardid]; i++)
	{
		if (rowflag[i] == 0)continue;
		for (unsigned j = 0; j < Rowsize[i]; j += blockbyte)
		{
			get_block_frequency(i, j);
		}
	}
	vblock = vector<BLOCK>(block_fre.begin(), block_fre.end());
	sort(vblock.begin(), vblock.end(), cmp);
}
long find_map(unsigned rows, unsigned startcol)
{
	unsigned char temp[blockbyte];
	for (unsigned i = 0; i < blockbyte; i++)
		temp[i] = 0;
	for (unsigned i = startcol; i < startcol + blockbyte&&i<Rowsize[rows]; i++)
		temp[i - startcol] = Row[rows][i];
	unsigned long ltemp = 0;
	for (unsigned i = 0; i < blockbyte; i++)
	{
		ltemp *= pow(2, 8);
		ltemp |= (unsigned)temp[i];
	}
	map<unsigned long, unsigned>::iterator it;
	it = pattern_to_id.find(ltemp);
	if (it != pattern_to_id.end())
	{
		NT_count++;
		return it->second;
	}
	else
	{
		long  pos = -1, count1 = INT_MAX, popl = 0, popv = 0;
		for (unsigned i = 0; i < block_map_size&&i < vblock.size(); i++)
		{
			if ((vblock[i].first | ltemp) == vblock[i].first)
			{
				popv = _mm_popcnt_u32((unsigned)vblock[i].first);
				popl = _mm_popcnt_u32((unsigned)ltemp);
				if (count1>popv)
				{
					count1 = popv;
					pos = i;
				}
			}
		}
		if (pos == 0)fix_count++;
		if (pos > 0)S_count++;
		if (pos < 0)T_count++;
		return pos;
	}
}
double map_all_rows()
{
	for (unsigned i = 0; i < block_map_size; i++)
	{
		pattern_to_id[vblock[i].first] = i;
	}
	long long sumrowsize = 0;
	unsigned colbyte = maxcol[myshardid] / 8;
	for (unsigned i = 0; i < maxrow[myshardid]; i++)
	{
		if (rowflag[i] == 0)
		{ 
			sumrowsize += Rowsize[i];
			continue; 
		}
		unsigned char *temprow;
		temprow = new unsigned char[colbyte];
		fill(temprow, temprow + colbyte, 0);
		unsigned rowsize = 0;
		unsigned blockcount = 0;
		for (unsigned j = 0; j < Rowsize[i]; j += blockbyte, blockcount++)
		{
			long pos = find_map(i, j);
			if (pos == -1)
			{
				for (unsigned t = j; t < j + blockbyte&&t<Rowsize[i]; t++)
					temprow[rowsize++] = Row[i][t];
			}
			else
			{
				temprow[rowsize++] = pos / 256;
				temprow[rowsize++] = pos % 256;
				reduce_block_flag[i][blockcount] = 1;
			}
		}
		for (int j = 0; j < rowsize; j++)
			Row[i][j] = temprow[j];
		Rowsize[i] = rowsize;
		sumrowsize += rowsize;
		delete[] temprow;
	}
	return (double)sumrowsize / maxrow[myshardid];
}
void write_file()
{
	string filename2 = "";
	string filenamep = "";

	FILE *file2 = fopen(filename2.c_str(), "wb");
	FILE *file_P = fopen(filenamep.c_str(), "wb");
	
	for (unsigned i = 0; i < maxrow[myshardid]; i++)
	{
	
		for (unsigned j = 0; j < Rowsize[i] / 2; j++)
		{
			unsigned short tmp = (unsigned short)Row[i][j * 2] * 256 + Row[i][j * 2 + 1];
			fwrite(&tmp, sizeof(unsigned short), 1, file2);
		}
	}
	for (unsigned j = 0; j < block_map_size&&j < vblock.size(); j++)
	{
		unsigned tmp = (unsigned)vblock[j].first;
		fwrite(&tmp, sizeof(unsigned), 1, file_P);
	}
	fclose(file2);
	fclose(file_P);

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
 

	maxcol[myshardid] = ((maxcolold[myshardid] + blocksize - 1) / blocksize)*blocksize;
	unsigned colblock = maxcol[myshardid] / blocksize, colbyte = maxcol[myshardid] / 8;
	Row = new unsigned char*[maxrow[myshardid]];
	Rowsize = new unsigned[maxrow[myshardid]];
	reduce_block_flag = new unsigned*[maxrow[myshardid]];
	rowflag = new unsigned char[maxrow[myshardid] + 10];
	for (unsigned i = 0; i < maxrow[myshardid]; i++)
	{
		Row[i] = new unsigned char[colbyte];
		reduce_block_flag[i] = new unsigned[colblock];
     fill(Row[i],Row[i]+colbyte,0);
	}
	cout << "maxrow=" << maxrow[myshardid] << "\tmaxcol(ali)=" << maxcol[myshardid] << endl;
	cout << "colbyte=" << colbyte << endl;
	cout << "colblock=" << colblock << endl;
}
int main(int argc, char *argv[])
{
	myshardid = atoi(argv[1]);
	rthresh = atof(argv[2]);
	init_data();
	read_file();
	cout << "read over" << endl;
	read_rowflag();
	cout << "read rowflagover" << endl;
	get_block_map();
	cout << "get_block_map and sort" << endl;
	long long sumfre = 0;
	sumfre = 0;
	for (int i = 0; i < block_map_size&&i < vblock.size(); i++)
	{
		if (i != 0)
			sumfre += vblock[i].second;
	}
	cout << "all we have " << min(block_map_size, (unsigned)vblock.size()) << " patterns" << endl;
	cout << "average NT fre=" << (double)sumfre / min(block_map_size, (unsigned)vblock.size()) << endl;
	double avrowsize = map_all_rows();
	cout << "map_all_rows" << endl;
	cout << "average row size=" << avrowsize << endl;
	write_file();
	cout << "write over" << endl;
	cout << "T_count=" << T_count << "\tNT_count=" << NT_count << endl;
	cout << "S_count=" << S_count << "fixed count=" << fix_count << endl;
}
