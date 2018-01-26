#include<iostream>
#include "ml.h"
using namespace std;
using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;

/////////////////////////////////////////////////////////
bool qusort(mI mi1, mI mi2)
{
	return (min(mi1.y1,mi1.y2) < min(mi2.y1,mi2.y2)) ||((min(mi1.y1,mi1.y2) == min(mi2.y1,mi2.y2)) && (max(mi1.y1,mi1.y2)<max(mi2.y1,mi2.y2)));
}
ccov makeChromBucket(int refLen)
{
	ccov v;
	for(int i=0;i<refLen;i++)
	{
		v.push_back(0);
	}
return v;
}	

///////////////////////////////////////////////////////////
void storeCords(ccov & masterRef,ccov & masterQ, mI & mi)
{

	int ty1 = 0, ty2 =0;
	
	if(mi.y1 > mi.y2)//if reverse oriented
	{
		ty1 = mi.y2;
		ty2 = mi.y1;
	}
	if(mi.y1 < mi.y2)//forward oriented
	{
		ty1 = mi.y1;
		ty2 = mi.y2;
	}
	for(int i = mi.x1-1; i<mi.x2;i++)
	{
		masterRef[i]++;
	}
	
	for(int j = ty1-1; j<ty2;j++)
	{
		masterQ[j]++;
	}
}
///////////////////////////////////////////////////////
void storeCords(map<int,vq> & mRef, mI & mi)
{	
	int refC = mi.x1;
	int ci = refC * (-1); //ci is minus i
	qord temp;
	if(mi.y1 < mi.y2 ) //both are on the same strand
	{
		int qC = mi.y1;
		while( refC<mi.x2+1)
		{
		
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //if this position does not have a indel
			{
				refC++;
				qC++;
				temp.name = mi.qn;
				temp.cord = qC-1;
				mRef[refC-1].push_back(temp);
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) != mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //position has insertion
			{
				refC++;
				temp.name = mi.qn;
				temp.cord = qC-1;
				mRef[refC-1].push_back(temp);
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())) //position has deletion
			{
				qC++;
			}
		}
	}
	if(mi.y1 > mi.y2 )//if two are on different strands
	{
		int qC = mi.y1; //y1 is bigger than y2
	
		while(refC<mi.x2+1)
		{
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //if this position does not have a indel
			{
				refC++;
				qC--;
				temp.name = mi.qn;
				temp.cord = qC+1;	
				mRef[refC-1].push_back(temp);
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) != mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //position has insertion
			{
				refC++;
				temp.name = mi.qn;
				temp.cord = qC+1;
				mRef[refC-1].push_back(temp);				
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())) //position has deletion
			{
				qC--;
			}
		}
	}
}
		
	
////////////////////////////////////////////////////
void findQuery(map<int,vq> & mRef, mI & mi,string & refName)
{
	
	if((mRef[mi.x1]. size() >0) && (mRef[mi.x2].size()>0)) //if the position has something
	{
		for(unsigned int j=0; j<min(mRef[mi.x1].size(),mRef[mi.x2].size());j++)//traverse all items at position i
		{
			if(mRef[mi.x1][j].name == mRef[mi.x2][j].name)
			{
				cout<<refName<<"\t"<<mi.x1<<"\t"<<mi.x2<<"\t"<<mRef[mi.x1][j].name<<"\t"<<mRef[mi.x1][j].cord<<"\t"<<mRef[mi.x2][j].cord<<endl;
			}
		}
	}
	else //if the position is empty
	{
	//dunno what to do yet		
	}	
	
}
/////////////////////////////////////////////////
vector<string> splitField(string & str, char c)
{
	size_t pos =1, pos1 =0;
	vector<string> vs;
	string tempstr;
	while(pos1 <str.size())
	{
		pos1 = str.find(c,pos);
		if(pos1 < str.size())
		{
			tempstr = str.substr(pos-1,pos1-pos+1);
			if(tempstr[0] == c)
			{
				tempstr = tempstr.substr(1); //remove the preceding delimiter
			}
			pos = pos1+1;
			vs.push_back(tempstr);
		}
	}
	tempstr = str.substr(pos);
	vs.push_back(tempstr);
	return vs;
}
