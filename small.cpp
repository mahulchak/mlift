#include "sv.h"
#include<iostream>
#include "seqIO.h"

void readUniq(ifstream & fin,vector<mI> & cm, map<int,vector<qord> > & umRef) //records one to one mapping
{
	string refName,qName,indexAln,line;
	size_t pos1,pos2,namePos;
	int refLen=0,qLen=0,count =0,indelPos =0, refStart =0, qStart =0, refEnd =0, qEnd = 0;
	vector<int> vi;
	
	mI tempmi;
//	fin.open(argv[1]);
	while(getline(fin,line))
	{
		if(line.find('>') != string::npos)//start of an aligning chromosome description
		{
			refName = line.substr(1,line.find(' ')-1);
			pos1 = line.find(' '); //position of the first space
			pos2 = line.find(' ',pos1+1);//position of the second space
			qName = line.substr(pos1+1, pos2-pos1); //up to the second space
			pos1 = line.find(' ',pos2+1); //recycling pos1 to find pos3
			refLen = stoi(line.substr(pos2+1,pos1-pos2));//reference length
			qLen = stoi(line.substr(pos1));//from last space till end 
			indexAln = refName + qName;
			count = -1;
		}
		if((line.size() <10) && (refName != "") && (count > -1))
		{
			indelPos = abs(stoi(line));
			refStart = refStart + indelPos;
			if(indelPos <0)
			{
				refStart = refStart * (-1);
			}
			vi.push_back(refStart);
			if(indelPos ==0) //reached the end of the indel description
			{
				tempmi.mv = vi;
				if(find(cm.begin(),cm.end(),tempmi) != cm.end()) //if tempmi is present within the unique alignments
				{
					storeCords(umRef,tempmi); //add the new mRef or umRef here
				}
				vi.clear();//reset it once its values are used
			}
			count++;
		}
		if((line.find('>') == string::npos) && (line.size() >10) && (refName != "")) //when describing alignment segments
		{
			tempmi.rn = refName;
			tempmi.qn = qName;
			refStart = stoi(line,&pos1);
			refEnd = stoi(line.substr(pos1),&pos2);
			qStart = stoi(line.substr(pos1+pos2), &namePos);
			qEnd = stoi(line.substr(pos1+pos2+namePos));
			tempmi.x1 = refStart;
			tempmi.x2 = refEnd;
			tempmi.y1 = qStart;
			tempmi.y2 = qEnd;
			count = 0;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void readfasta(ifstream & fin,map<string, string> & fastaseq) //reading fasta files
{
	string str,index;
	while(getline(fin,str))
	{
		if(str[0] == '>')
		{
			index = str.substr(1);
		}
		if(str[0] != '>')
		{
			fastaseq[index].append(str);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void callSmall(map<int,vector<qord> > & umRef, map<string,string> & refseq, map<string, string> & qseq) //calls SNPs and indels contained within aligned regions
//{
	
