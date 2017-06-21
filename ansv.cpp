#include "sv.h"
#include<iostream>

using namespace std;
using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;

void annotGaps(vector<mI> & cm,map<int,vq> & mRef,ccov & masterRef, ccov & masterQ, vector<mI> & cnv)
{
	int refOvl =0; //overlap between reference intervals
	int qOvl = 0; //overlap between query intervals
	vector<double> cov(2);
	mI gapmi;
	for(unsigned int i=1; i< cm.size();i++)
	{
		refOvl = cm[i].x1 - cm[i-1].x2;
		qOvl = min(cm[i].y1,cm[i].y2) - max(cm[i-1].y1,cm[i-1].y2);
		if(refOvl < 0)
		{
			if(qOvl<0)
			{
				//cout<<"BD "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
			}
			else
			{
				//cout<<"INS "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
				gapmi.y1 = min(cm[i-1].y2,cm[i].y1);
				gapmi.y2 = max(cm[i].y1,cm[i-1].y2);
				findCnvOverlap(cnv,gapmi);
				
			}
		}
		if(!(refOvl<0))
		{
			if(refOvl > qOvl)
			{
				//cout<<"DEL "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
			}
			if(refOvl < qOvl)
			{
				//cout<<"INS "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
				gapmi.y1 = min(cm[i-1].y2,cm[i].y1);
                                gapmi.y2 = max(cm[i].y1,cm[i-1].y2);
                                findCnvOverlap(cnv,gapmi);
			}
		}

	}
} 
//////////////////////////////////////////////////////////////////////////////////
	
	
mI findCnvOverlap(vector<mI> & cnv,mI & mi) //returns a cnv if it overlaps a gap
{
	mI tempmi;
	int ovl =0;
	unsigned int count = 0;
	//bool found = 0;	
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
cout<<mi.y1<<" "<<mi.y2<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.y1<<" "<<tempmi.y2<<ovl<<" "<<(double(tempmi.y2 -tempmi.y1))/double(ovl)<<endl;					
					//return tempmi;
					//found = true;
					break;
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
					break;
cout<<mi.y1<<" "<<mi.y2<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.y1<<" "<<tempmi.y2<<ovl<<" "<<(double(tempmi.y2 -tempmi.y1))/double(ovl)<<endl;
				}
			}
		}
		count++;
//cout<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.qn<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;
	}
if(count != cnv.size())
{
	cout<<count<<" "<<cnv.size()<<" "<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.qn<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;
}
return tempmi;
}			
		
