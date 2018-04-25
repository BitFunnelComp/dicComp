#include<iostream>
#include<fstream>
#include<string>
#include<cstdio>
#include<cstdint>
#include <algorithm>
#include<math.h>
#include<string.h>
#include<vector>
#include<set>
using namespace std;
const unsigned maxrow = 7413-3;
const unsigned maxcol =1069958;
const unsigned alig=32;
const unsigned alicol=((maxcol+alig-1)/alig)*alig;
const unsigned colbyte = (maxcol + 7) / 8;
unsigned char bit[maxrow][maxcol];
unsigned char temp[maxrow][alicol];

const unsigned blocksize = 8;
const unsigned blockcount = alicol / blocksize;
struct col_1{ unsigned long count; unsigned long old_col; };
col_1 col_1_count[maxcol];

unsigned long block1_num[9];


void read_file()
{
	FILE *file = fopen("");
	for (unsigned i = 0; i<maxrow; i++)
	{
		for (unsigned j = 0; j<colbyte; j++)
		{
			unsigned char t = 0;
			fread(&t, sizeof(unsigned char), 1, file);
			for (unsigned k = 128,col=j*8; k>0&&col<(j+1)*8&&col<maxcol; k /= 2,col++)
			{
				bit[i][col] = t / k;
				t %= k;
			}
		}
	}
	fclose(file);
}
bool cmp(col_1 &a, col_1& b)
{
	return a.count>b.count;
}
void get_density_sort_interval()
{
	for (unsigned i = 0; i < maxcol; i++)
	{
		col_1_count[i].old_col = i; col_1_count[i].count = 0;
		for (unsigned j = 0; j < maxrow; j++)
		{
			if ((unsigned)bit[j][i] == 1)
				col_1_count[i].count++;
		}
	}
	
	sort(col_1_count, col_1_count + maxcol, cmp);
	unsigned oldcolid = 0, newcolid = 0,colid=0;
	for (unsigned turn = 0; turn < blocksize&&colid<maxcol; turn++)
	{
		for (unsigned b = 0; b < blockcount&&colid<maxcol; b++)
		{
			newcolid = b*blocksize + turn;
			if (newcolid >= alicol){cout<<"wrong"<<endl;continue;}
			oldcolid = col_1_count[colid++].old_col; 
			for (unsigned row = 0; row < maxrow; row++)
			{
				temp[row][newcolid] = bit[row][oldcolid];
			}
		}
	}

}
void write_file()
{
	FILE *file = fopen("", "wb");
	for (unsigned i = 0; i < maxrow; i++)
	{
		for (unsigned j = 0; j < alicol; j += 8)
		{
			unsigned char tmp = 0;
			for (unsigned t = j,k=128;k>0&& t < j+8 && t < alicol; t++,k/=2)
			{
				tmp += k*temp[i][t];
			}
			fwrite(&tmp, sizeof(unsigned char), 1, file);
		}
	}
	fclose(file);
}
unsigned check_block_1(unsigned start)
{
	unsigned  count = 0;
	for (unsigned i = 0; i < maxrow; i++)
	{
		int ttemp = 0;
		for (unsigned j = start*blocksize; j < (start + 1)*blocksize&&j<alicol; j++)
		{
			if ((unsigned)temp[i][j] == 1)
				ttemp++;
		}
		block1_num[ttemp]++;
	}
	return count; 
}

int main()
{
	cout<<"maxcol="<<maxcol<<"\talicol="<<alicol<<endl;
	read_file();
 cout<<"read over"<<endl;
	get_density_sort_interval();
 cout<<"sort over"<<endl;
	
	write_file();
 cout<<"write over"<<endl;
 for (unsigned i = 0; i < blockcount; i++)
	{
		check_block_1(i);
	}
 unsigned long sum=0;
 for(unsigned i=0;i<9;i++)
 {
     cout<<"1 count="<<i<<"\tblock count="<<block1_num[i]<<endl;
     sum+=block1_num[i];
 }
 cout<<"all block count="<<sum<<endl;
 cout<<(unsigned long)maxrow*blockcount<<endl;
 cout<<"over"<<endl;
}
