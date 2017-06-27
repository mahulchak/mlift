#include "sv.h"
#include<iostream>

using namespace std;
using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;

void annotGaps(vector<mI> & cm, map<int,vq> & mRef, ccov & masterRef, ccov & masterQ, vector<mI> & cnv, map<int,vector<qord> > & umRef, string & refseq, string & qseq, vector<int> & seqLen)
{
	int refOvl =0; //overlap between reference intervals
	int qOvl = 0; //overlap between query intervals
	vector<double> cov(2);
	mI gapmi,cnvmi,tempmi;
	vector<mI> storedCNV;
	ofstream fout;
	fout.open("sv.txt");
	for(unsigned int i=1; i< cm.size();i++)
	{
		tempmi = cm[i];
		callSmall(tempmi,umRef,refseq,qseq,seqLen);
		refOvl = cm[i].x1 - cm[i-1].x2;
		qOvl = min(cm[i].y1,cm[i].y2) - max(cm[i-1].y1,cm[i-1].y2);
		cnvmi = findDup(cm[i-1],cm[i]);
		if(refOvl < 0)
		{
			if(qOvl<0)
			{
				//cout<<"BD "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
			}
			else
			{
				//cout<<"INS "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
				fout<<"INS "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
				gapmi.y1 = min(cm[i-1].y2,cm[i].y1);
				gapmi.y2 = max(cm[i].y1,cm[i-1].y2);
fout<<"DEBUG "<<"previous: "<<gapmi.y1<<" "<<gapmi.y2;
				findCnvOverlap(cnv,gapmi,storedCNV);
fout<<" now:"<<gapmi.y1<<" "<<gapmi.y2<<endl;
			}
		}
		if(!(refOvl<0))
		{
			if(refOvl > qOvl)
			{
				//cout<<"DEL "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
				fout<<"DEL "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
			}
			if(refOvl < qOvl)
			{
				//cout<<"INS "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
				fout<<"INS "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
				gapmi.y1 = min(cm[i-1].y2,cm[i].y1);
                                gapmi.y2 = max(cm[i].y1,cm[i-1].y2);
                                findCnvOverlap(cnv,gapmi,storedCNV);
			}
		}
		if(cnvmi.x1 != 0)
		{
fout<<"cnv "<<cnvmi.rn<<" "<<cnvmi.x1<<" "<<cnvmi.x2<<" "<<cnvmi.qn<<" "<<cnvmi.y1<<" "<<cnvmi.y2<<endl;
		}
		if(cm[i].y1>cm[i].y2)//inversion
		{
fout<<"inv "<<cm[i].rn<<" "<<cm[i].x1<<" "<<cm[i].x2<<" "<<cm[i].y2<<" "<<cm[i].y1<<endl;
		}

	}
	fout.close();
} 
//////////////////////////////////////////////////////////////////////////////////
	
	
void findCnvOverlap(vector<mI> & cnv,mI & mi,vector<mI> & storedCNV) //returns a cnv if it overlaps a gap
{
	mI tempmi;
	int ovl =0;
	unsigned int count = 0;
	for(unsigned int i=0;i<cnv.size();i++)
	{
		tempmi = cnv[i];
		if(tempmi.y1 <tempmi.y2) //forward strand. assumes that the gap is also forward oriented
		{
			if(!(tempmi.y2<mi.y1) && !(tempmi.y1>mi.y2))
			{
				ovl = min(tempmi.y2,mi.y2) - max(tempmi.y1,mi.y1);
			
				if((double(ovl)/double(tempmi.y2 -tempmi.y1)) >0.9)
				{
					//return tempmi;
					//found = true;
					if((find(storedCNV.begin(),storedCNV.end(),tempmi) == storedCNV.end()) && (mi.y1 != mi.y2))
					{
//cout<<mi.y1<<" "<<mi.y2<<" "<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.qn<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;
						storedCNV.push_back(tempmi);
						mi.y1 = tempmi.y2 + 1;//reduce the gap
						//mi.y2 = tempmi.y2 + 1;
						findCnvOverlap(cnv,mi,storedCNV);			
						break;
					}
				}
				
			}
		}
		if(tempmi.y1 > tempmi.y2) // reverse strand. assumes that the gap is also reverse oriented
		{
			if(!(tempmi.y2 > mi.y1) && !(tempmi.y1 < mi.y2))
			{
				ovl = min(tempmi.y1,mi.y1) - max(tempmi.y2,mi.y2);
			
				if((double(ovl)/double(tempmi.y1 -tempmi.y2))>0.9) 
				{
					//return tempmi;
					//found = true;
					if((find(storedCNV.begin(),storedCNV.end(),tempmi) == storedCNV.end()) && (mi.y1 != mi.y2))
					{
//cout<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.qn<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;
						storedCNV.push_back(tempmi);
						mi.y1 = tempmi.y2 - 1;
						//mi.y2 = tempmi.y2 -1;
						findCnvOverlap(cnv,mi,storedCNV);
						break;
					}
				}
			}
		}
		count++;
	}
if(count != cnv.size())
{
//	cout<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.qn<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;
}

}			
/////////////////////////////////////////////////////////////////
mI findDup(mI & mi1, mI & mi2)
{
	mI tempmi;
	tempmi.x1= 0;//to use it as a filter for later
//cout<<"candidate "<<mi1.x1<<" "<<mi1.x2<<" "<<mi2.x1<<" "<<mi2.x2<<endl;
	if(mi2.x1 < mi1.x2) //i.e. mi2.x1 falls within the previous interval. this works because mums are sorted by reference cord
	{
//cout<<"candidate "<<mi1.x1<<" "<<mi1.x2<<" "<<mi2.x1<<" "<<mi2.x2<<endl;
		if(mi2.x2 > mi1.x2) //mi1 amd mi2 just overlaps but one is not contained within another
		{
			tempmi.rn = mi1.rn;
			tempmi.x1 = mi2.x1;
			tempmi.x2 = mi1.x2;
			tempmi.qn = mi1.qn;
			if(mi2.y1 <mi2.y2) //forward oriented
			{
				tempmi.y1 = mi2.y1;
				tempmi.y2 = mi2.y1 + (mi1.x2 - mi2.x1);
			}
			if(mi2.y1 > mi2.y2)//reverse oriented
			{
				tempmi.y1 = mi2.y1;
				tempmi.y2 = mi2.y1 - (mi1.x2 - mi2.x1);
			}
//cout<<"cnv "<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.qn<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;				
		}
		if(!(mi2.x2>mi1.x2))//mi2 is contained within mi1
		{
			tempmi.rn = mi1.rn;
			tempmi.x1 = mi2.x1;
			tempmi.x2 = mi1.x2;
			tempmi.qn = mi1.qn;
			tempmi.y1 = mi2.y1;
			tempmi.y2 = mi2.y2;
//cout<<"cnv "<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.qn<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;
		}
	}
	return tempmi;
}
/////////////////////////////////////////////////////////////////////////////	
	
