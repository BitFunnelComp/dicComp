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
const unsigned rankcount = 6;
const unsigned alig = 64;
unsigned maxrow[shardcount][rankcount];
unsigned maxcol[shardcount][rankcount];
typedef pair<unsigned long, unsigned long> BLOCK;
unsigned char **Row[rankcount];
unsigned *Rowsize[rankcount];
map<unsigned long, unsigned long>block_fre[rankcount];
vector<BLOCK>vblock[rankcount];
unsigned long T_count = 0;
unsigned long NT_count = 0;
unsigned long fix_count = 0;
unsigned long S_count = 0;
map<unsigned long, unsigned>pattern_to_id[rankcount];
unsigned myshardid = 0;
unsigned thresh = 32;
double rthresh = 0;
unsigned char *rowflag[rankcount];
void read_file(unsigned rankid)
{
	string filename = "";
	FILE *file1 = fopen(filename.c_str(), "rb");
	unsigned colbyte = maxcol[myshardid][rankid] / 8, collong = maxcol[myshardid][rankid]/64;
	for (unsigned i = 0; i < maxrow[myshardid][rankid]; i++)
	{
		for (unsigned j = 0; j < collong; j++)
		{
			unsigned long temp=0;
			fread(&temp, sizeof(unsigned long), 1, file1);
			for (unsigned long t = j * 8,k=(unsigned long)pow(2,56); t < (j + 1) * 8&&k>0; t++,k/=256)
			{
				Row[rankid][i][t] = temp/k;
				temp %= k;
			}
		}
		Rowsize[rankid][i] = colbyte;
	}
	fclose(file1);

}
void read_rowflag(unsigned rankid)
{
	string filename = "";
	FILE *filem = fopen(filename.c_str(), "rb");
	unsigned char temp;
	unsigned u_temp = 0;
	unsigned T = 0, NT = 0;
	for (unsigned i = 0; i <maxrow[myshardid][rankid]; i += 8)
	{
		fread(&temp, sizeof(unsigned char), 1, filem);
		u_temp = (unsigned)(unsigned char)temp;
		for (unsigned t = i, k = 128; t < i + 8 && t <maxrow[myshardid][rankid]; t++, k /= 2)
		{
			rowflag[rankid][t] = u_temp / k;
			u_temp %= k;
			if (rowflag[rankid][t] == 0)T++;
			else  NT++;
		}

	}
	fclose(filem);
	cout << "T rows count=" << T << endl;
	cout << "all rows count=" << T + NT << endl;
}
void get_block_frequency(unsigned rankid,unsigned rows, unsigned startcol)
{
	unsigned char temp[blockbyte];
	for (unsigned i = 0; i < blockbyte; i++)
		temp[i] = 0;
	for (unsigned i = startcol; i < startcol + blockbyte&&i<Rowsize[rankid][rows]; i++)
		temp[i - startcol] = Row[rankid][rows][i];
	unsigned long  ltemp = 0;
	for (unsigned i = 0; i < blockbyte; i++)
	{
		ltemp <<=8;
		ltemp |= (unsigned)temp[i];
	}
	block_fre[rankid][ltemp]++;
}
bool cmp(BLOCK &b1, BLOCK &b2)
{
	if (b1.second> b2.second)
		return true;
	return false;
}
void get_block_map(unsigned rankid)
{
	unsigned ltemp = (unsigned long)pow(2, 32) - 1;
	block_fre[rankid][ltemp] = (unsigned long)pow(2, 60);
	for (unsigned i = 0; i < maxrow[myshardid][rankid]; i++)
	{
		if (rowflag[rankid][i] == 0)continue;
		for (unsigned j = 0; j < Rowsize[rankid][i]; j += blockbyte)
		{
			get_block_frequency(rankid,i, j);
		}
	}
	vblock[rankid] = vector<BLOCK>(block_fre[rankid].begin(), block_fre[rankid].end());
	sort(vblock[rankid].begin(), vblock[rankid].end(), cmp);
}
long find_map(unsigned rankid,unsigned rows, unsigned startcol)
{
	unsigned char temp[blockbyte];
	for (unsigned i = 0; i < blockbyte; i++)
		temp[i] = 0;
	for (unsigned i = startcol; i < startcol + blockbyte&&i<Rowsize[rankid][rows]; i++)
		temp[i - startcol] = Row[rankid][rows][i];
	unsigned long ltemp = 0;
	for (unsigned i = 0; i < blockbyte; i++)
	{
		ltemp <<=8;
		ltemp |= (unsigned)temp[i];
	}
	map<unsigned long, unsigned>::iterator it;
	it = pattern_to_id[rankid].find(ltemp);
	if (it != pattern_to_id[rankid].end())
	{
		NT_count++;
		return it->second;
	}
	else
	{
		long  pos = -1, count1 = INT_MAX, popl = 0, popv = 0;
		for (unsigned i = 0; i < block_map_size&&i < vblock[rankid].size(); i++)
		{
			if ((vblock[rankid][i].first | ltemp) == vblock[rankid][i].first)
			{
				popv = _mm_popcnt_u32((unsigned)vblock[rankid][i].first);
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
double map_all_rows(unsigned rankid)
{
	for (unsigned i = 0; i < block_map_size&&i<vblock[rankid].size(); i++)
	{
		pattern_to_id[rankid][vblock[rankid][i].first] = i;
	}
	long long sumrowsize = 0;
	unsigned colbyte = maxcol[myshardid][rankid] / 8;
	for (unsigned i = 0; i < maxrow[myshardid][rankid]; i++)
	{
		if (rowflag[rankid][i] == 0)
		{
			sumrowsize += Rowsize[rankid][i];
			continue;
		}
		unsigned char *temprow;
		temprow = new unsigned char[colbyte];
		fill(temprow, temprow + colbyte, 0);
		unsigned rowsize = 0;
		unsigned blockcount = 0;
		for (unsigned j = 0; j < Rowsize[rankid][i]; j += blockbyte, blockcount++)
		{
			long pos = find_map(rankid,i, j);
			if (pos == -1)
			{
				for (unsigned t = j; t < j + blockbyte&&t<Rowsize[rankid][i]; t++)
					temprow[rowsize++] = Row[rankid][i][t];
				cout << "Wrong !!!!!!" << endl; return 0;
			}
			else
			{
				temprow[rowsize++] = pos / 256;
				temprow[rowsize++] = pos % 256;
			}
		}
		for (int j = 0; j < rowsize; j++)
			Row[rankid][i][j] = temprow[j];
		Rowsize[rankid][i] = rowsize;
		sumrowsize += rowsize;
		delete[] temprow;
	}
	return (double)sumrowsize / maxrow[myshardid][rankid];
}
void write_file(unsigned rankid)
{
	string filename2 = "";
	string filenamep = "";
	FILE *file2 = fopen(filename2.c_str(), "wb");
	FILE *file_P = fopen(filenamep.c_str(), "wb");

	for (unsigned i = 0; i < maxrow[myshardid][rankid]; i++)
	{
		for (unsigned j = 0; j < Rowsize[rankid][i] / 2; j++)
		{
			unsigned short tmp = (unsigned short)Row[rankid][i][j * 2] * 256 + Row[rankid][i][j * 2 + 1];
			fwrite(&tmp, sizeof(unsigned short), 1, file2);
		}
	}
	for (unsigned j = 0; j < block_map_size&&j < vblock[rankid].size(); j++)
	{
		unsigned tmp = (unsigned)vblock[rankid][j].first;
		fwrite(&tmp, sizeof(unsigned), 1, file_P);
	}
	fclose(file2);
	fclose(file_P);

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

	maxcol[0][0] = 2647295;
	maxcol[1][0] = 2614477;
	maxcol[2][0] = 5519536;
	maxcol[3][0] = 7373033;
	maxcol[4][0] = 4503823;
	maxcol[5][0] = 1872882;
	maxcol[6][0] = 497135;
	maxcol[7][0] = 157717;
	maxcol[8][0] = 19285;

	cout << "shard id=" << myshardid << endl;
	maxcol[myshardid][5] = (maxcol[myshardid][0] + pow(2, 5) - 1) / pow(2, 5);
	maxcol[myshardid][5] = (maxcol[myshardid][5] + 63) / 64 * 64;
	for (unsigned r = 0; r < rankcount; r++)
	{
		cout << "rank " << r << endl;
		maxcol[myshardid][r] = maxcol[myshardid][5] * pow(2, 5 - r);
		cout << "maxrow=" << maxrow[myshardid][r] << "\tmaxcol(ali)=" << maxcol[myshardid][r] << endl;
	}
	
	
	for (unsigned i = 0; i < rankcount; i++)
	{
		Row[i] = new unsigned char*[maxrow[myshardid][i]];
		Rowsize[i] = new unsigned [maxrow[myshardid][i]];
		rowflag[i] = new unsigned char[maxrow[myshardid][i] + 10];
		for (unsigned j = 0; j < maxrow[myshardid][i]; j++)
		{
			Row[i][j] = new unsigned char[maxcol[myshardid][i] / 8];
			fill(Row[i][j], Row[i][j] + maxcol[myshardid][i] / 8, 0);
		}
	}

	
}
int main(int argc, char *argv[])
{
	myshardid = atoi(argv[1]);
	rthresh = atof(argv[2]);
	init_data();
	for (unsigned r = 0; r < rankcount; r++)
	{
		cout << "****************************************************************" << endl;
		cout << "rankid=" << r << endl;
		read_file(r);
		cout << "read over" << endl;
		read_rowflag(r);
		cout << "read rowflagover" << endl;
		get_block_map(r);
		cout << "get_block_map and sort" << endl;
		long long sumfre = 0;
		sumfre = 0;
		for (int i = 0; i < block_map_size&&i < vblock[r].size(); i++)
		{
			if (i != 0)
				sumfre += vblock[r][i].second;
		}
		cout << "all we have " << min(block_map_size, (unsigned)vblock[r].size()) << " patterns" << endl;
		cout << "average NT fre=" << (double)sumfre / min(block_map_size, (unsigned)vblock[r].size()) << endl;
		double avrowsize = map_all_rows(r);
		cout << "map_all_rows" << endl;
		cout << "average row size=" << avrowsize << endl;
		write_file(r);
		cout << "write over" << endl;
		cout << "T_count=" << T_count << "\tNT_count=" << NT_count << endl;
		cout << "S_count=" << S_count << "fixed count=" << fix_count << endl;
		T_count = 0; NT_count = 0; S_count = 0; fix_count = 0;
	}
}
