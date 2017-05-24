#include "sv.h"
#include<iostream>

using namespace std;
using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;

void annotGaps(vector<mI> & cm,map<int,vq> & mRef)
{
	int refOvl =0; //overlap between reference intervals
	int qOvl = 0; //overlap between query intervals
	for(unsigned int i=1; i< cm.size();i++)
	{
		refOvl = cm[i].x1 - cm[i-1].x2;
		qOvl = min(cm[i].y1,cm[i].y2) - max(cm[i-1].y1,cm[i-1].y2);
		if(refOvl < 0)
		{
			if(qOvl<0)
			{
				cout<<"BD "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
			}
			else
			{
				cout<<"INS "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
			}
		}
		if(!(refOvl<0))
		{
			//if(qOvl<0)
			//{
			//	cout<<"DEL "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;	
			//}
			if(refOvl > qOvl)
			{
				cout<<"DEL "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
			}
			if(refOvl < qOvl)
			{
				cout<<"INS "<<cm[i-1].x2<<" "<<cm[i].x1<<" "<<max(cm[i-1].y1,cm[i-1].y2)<<" "<<min(cm[i].y1,cm[i].y2)<<endl;
			}
		}

	}
} 
		
