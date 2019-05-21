#include<iostream>
#include "ml.h"
using namespace std;
using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;


///////////////////////////////////////////////////////
mI liftCords(mI & cm,mI & mi)
{	
	int refC = mi.x1;
	int ci = refC * (-1); //ci is minus i
	qord temp;
	mI tempmi,liftmi;//stores the positions as mI
	tempmi.rn = cm.rn;
	liftmi.rn = cm.rn;
	liftmi.x1 = cm.x1;
	liftmi.x2 = cm.x2;
	liftmi.qn = mi.qn;
	liftmi.y1 = 0;
	liftmi.y2 = 0;
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
				tempmi.x1 = refC-1;
				tempmi.x2 = refC-1;
				tempmi.y1 = temp.cord;
				tempmi.y2 = temp.cord;
				tempmi.qn = mi.qn;
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) != mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //position has insertion
			{
				refC++;
				temp.name = mi.qn;
				temp.cord = qC-1;
				tempmi.x1 = refC-1;
				tempmi.x2 = refC-1;
				tempmi.y1 = temp.cord;
				tempmi.y2 = temp.cord;
				tempmi.qn = mi.qn;
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())) //position has deletion
			{
				qC++;
			}
			//cout<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<endl;
			if((cm.x1 == tempmi.x1) && (liftmi.y1 ==0))//if liftmi hasn't been changed yet
			{
				liftmi.y1 = tempmi.y1;
			}
			if((cm.x2 == tempmi.x2) && (liftmi.y2 ==0))		
			{
				liftmi.y2 = tempmi.y2;
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
				tempmi.x1 = refC-1;
				tempmi.x2 = refC-1;
				tempmi.y1 = temp.cord;
				tempmi.y2 = temp.cord;
				tempmi.qn = mi.qn;
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) != mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //position has insertion
			{
				refC++;
				temp.name = mi.qn;
				temp.cord = qC+1;
				tempmi.x1 = refC-1;
				tempmi.x2 = refC-1;
				tempmi.y1 = temp.cord;
				tempmi.y2 = temp.cord;
				tempmi.qn = mi.qn;
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())) //position has deletion
			{
				qC--;
			}
			//cout<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<endl;
			if((cm.x1 == tempmi.x1) && (liftmi.y1 ==0))//if liftmi hasn't been changed yet
			{
				liftmi.y1 = tempmi.y1;
			}
			if((cm.x2 == tempmi.x2) && (liftmi.y2 ==0))
			{
				liftmi.y2 = tempmi.y2;
			}
		}
	}
	return liftmi;
}
		
	
////////////////////////////////////////////////////
vector<mI> findMum(vector<mI> & mums,mI & cm)
{
	unsigned int i =0;
	mI tempmi;
	vector<mI> vmi;
	while(!(cm.x2 < mums[i].x1))
	{
//cout<<"all\t"<<cm.rn<<'\t'<<cm.x1<<'\t'<<cm.x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<'\t'<<mums[i].mv.size()<<endl;
		//if(!((cm.x2 < mums[i].x1) && (cm.x1 > mums[i].x2)))
		if(!(cm.x2 < mums[i].x1) && !(cm.x1 > mums[i].x2))
		{
			tempmi = mums[i];
			vmi.push_back(tempmi);
cout<<"query\t"<<cm.rn<<'\t'<<cm.x1<<'\t'<<cm.x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<'\t'<<mums[i].mv.size()<<endl;
		}
		++i;
	}
	return vmi;
}
//////////////////////////////////////////////////////
void writeLift(vector<mI> & vmi,mI & cm,ofstream & fout)
{
	mI tempmi;
	for(unsigned int i=0;i<vmi.size();i++)
	{
		tempmi = liftCords(cm,vmi[i]);
		if((tempmi.y1 != 0) && (tempmi.y2 !=0))
		{
			fout<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<endl;
		}
		if(tempmi.y1 ==0)
		{
			fout<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<"NA"<<'\t'<<tempmi.y2<<endl;
		}
		if(tempmi.y2 == 0)
		{
			fout<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<"NA"<<endl;
		}
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
