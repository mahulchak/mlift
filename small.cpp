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

			indelPos = stoi(line);
			refStart = refStart + abs(indelPos);
			if(indelPos <0)
			{
				//refStart = refStart * (-1);
				vi.push_back(refStart*-1);
			}
			if(indelPos > 0)
			{
				vi.push_back(refStart);
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
void callSmall(string & refName,map<int,vector<qord> > & umRef, string & refseq, string & qseq) //just get the individual sequences passed
{
	int refPos =0,pos =0,refGap =0,qGap =0;
	qord lq,delStart;//lq is last qord :P
	vector<qord> tq; //creating one element
	vector<int> ti;
	
	for(map<int,vector<qord> >::iterator it= umRef.begin();it != umRef.end();it++)
	{
		pos = it->first;
		if(umRef[pos].size() >1)
		{
			sort(umRef[pos].begin(),umRef[pos].end());//sort the coordinates to remain consistent
		}
		if(umRef[pos].size() == 1) //a mapping is present
		{
			refGap = (pos - refPos);
			qGap = abs(lq.cord - umRef[pos][0].cord);
			if((refGap == 1) && (qGap == 1))
			{
				if(tq.size()>0)
				{
					cout<<"DEL "<<refName<<" "<<ti[0]<<" "<<ti[ti.size()-1]<<" "<<tq[0].name<<" "<<tq[0].cord<<" "<<tq[tq.size()-1].cord<<endl;
					tq.clear();
					ti.clear();
				}
				cout<<"SNP "<<refName<<" "<<pos<<" "<<umRef[pos][0].name<<" "<<umRef[pos][0].cord<<endl;
			}
			if((refGap == 1) && (qGap ==0))
			{
				tq.push_back(lq);
				ti.push_back(pos);
			}
			//if((refGap == 0) && (qGap == 1))
			if((refGap == 1) && (qGap>1))
			{
				cout<<"INS "<<refName<<" "<<pos<<" "<<pos<<" "<<umRef[pos][0].name<<" "<<lq.cord<<" "<<umRef[pos][0].cord<<endl;
				//gq.push_back(umRef[pos][0]);
				//gi.push_back(refPos);
				//if(tq.size() > 0)
				//{
				//	cout<<"DEL "<<refName<<" "<<ti[0]<<" "<<ti[ti.size()-1]<<" "<<tq[0].name<<" "<<tq[0].cord<<" "<<tq[tq.size()-1].cord<<endl;
				//	tq.clear();
				//	ti.clear();
				//}
					
			}	
			lq = umRef[pos][0];
			refPos = pos;
		}
	}
}				
