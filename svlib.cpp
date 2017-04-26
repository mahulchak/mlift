#include<iostream>
#include "sv.h"
using namespace std;
using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;


ccov makeChromBucket(int refLen)
{
	ccov v;
	for(unsigned int i=0;i<refLen;i++)
	{
		v.push_back(0);
	}
return v;
}	
///////////////////////////////////////////////////////////
int findDist(int & x1,int & y1, int & c) //c is the intercept of the absolute diagonal
{
	int c1 = 0, dist =0;
	
	c1 = abs(x1 - y1);
	
	dist = abs(c1 -c);
	dist = int(dist/sqrt(2));

	return dist;
}
///////////////////////////////////////////////////////////
//bool chkOverlap(int & x1, int & x2,int & y1, int &y2) 
bool detectShadow(mI & mum, vector<mI> & mums, unsigned int i) //return whether mum is shadow mum or not
{
	unsigned int count = 0;
	bool found = false;
	while((found == false) && (count<mums.size())) // end of mum reference cannot be beyond mums[i] start if mum is a shadow
	{
		if((!(mum.x1<mums[count].x1)) && (!(mum.x2>mums[count].x2)) && (count != i))
		{
			found =true;
		}
		if((!(mum.y1 < mums[count].y1)) && (!(mum.y2 > mums[count].y2)) && (count != i))
		{
			found = true;
		}
	count++;
	}

return found;
}
/////////////////////////////////////////////////////////
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
		
//////////////////////////////////////////////////////////
vector<double> getCoverage(mI & mi, ccov & masterRef, ccov & masterQ)
{
	int d = 0, cov = 0;
	double c;
	vector<double> cc;
//cout<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t";
	d = mi.x2 - mi.x1;
	for(int i = mi.x1-1;i<mi.x2;i++)
	{
		cov = cov + masterRef[i];	
		
	}
c = cov/double(d);
cc.push_back(c);

	cov = 0;
	d = abs(mi.y1-mi.y2);
	for(int i = min(mi.y1,mi.y2)-1;i<max(mi.y1,mi.y2);i++)
	{
		cov = cov + masterQ[i];
	}
c = cov/double(d);
cc.push_back(c);	
return cc;
//cout<<c<<endl;
}
		
		
		
//////////////////////////////////////////////////////
mI findClosest(mI & mi, vector<mI> & mums, unsigned int i,ccov & masterRef, ccov & masterQ)
{
	map<double,mI> storeDist;
	map<double,mI>::iterator it;
	double d1 =0, d2 =0, d = 0;
	int ty1 =0, ty2 =0; //to switch the coordinates of the inverted MUMs
	vector<double> vd;
	ty1 = mi.y2;
	
	if(mi.y1 > mi.y2) // if reverse oriented
	{
		ty1 = mi.y1;
	}

	for(unsigned int j=i+1; j <mums.size();j++)
	{
		vd = getCoverage(mums[j],masterRef,masterQ);
		if((mums[j].y1 > mums[j].y2) && (vd[0]<2)) //on the other strand
		{
			ty2 = mums[j].y2; //swap the values
			d1 = pow(abs(mi.x2 - mums[j].x1),2);
			d2 = pow(abs(ty1 - ty2),2);
              		d = sqrt(d1+d2);
		}
		if((mums[j].y1 < mums[j].y2) && (vd[0]<2))
		{
			d1 = pow(abs(mi.x2 -mums[j].x1),2);
			d2 = pow(abs(ty1 - mums[j].y1),2);
			d = sqrt(d1+d2);
		}
//cout << mi.x2<<"\t"<<mi.y2<<"\t"<<mums[j].x1<<"\t"<< mums[j].y1<<endl;
		storeDist[d] = mums[j];
		
	}
	
	it = storeDist.begin();
	 		
return it->second;
}
////////////////////////////////////////////////////
unsigned int reSet(mI & mi, vector<mI> & mums, unsigned int i)
{
	bool found = false;
	while(found == false)//could use overloaded = operator to compare the mI objects	
	{
		i++;
		//if((mi.x1 == mums[i].x1) && (mi.x2 == mums[i].x2) && (mi.y1 == mums[i].y1) && (mi.y2 == mums[i].y2))
		if( mi == mums[i])
		{
			found = true;
		}
	}
return i;
}
/////////////////////////////////////////////////
void recordShadow(unsigned int i, unsigned int j, vector<mI> & mums, vector<mI> & sm)
{
	for(unsigned int k = i; k<j;k++)
	{
		
		if(find(sm.begin(),sm.end(),mums[k]) == sm.end()) // if it does not exist within the sm
		{
			sm.push_back(mums[k]);
//cout<<mums[k].x1<<"\t"<<mums[k].x2<<"\t"<<mums[k].y1<<"\t"<<mums[k].y2<<endl;
		}
	}
}

//////////////////////////////////////////////
void splitByCoverage(chromPair & cp, ccov & chrom,vector<mI> & mums, ccov & masterRef, ccov & masterQ) // returns the percentage of gap filled by the mi in mums
{
	int cov=0, lastcov=0, nextcov =0;
	//vector<mI> mum;
	mI mi;
	vector<double> vd;
	mi.x1 =1;
	for(unsigned int i =1; i<chrom.size()-1;i++)
	{
		cov = chrom[i];
		lastcov = chrom[i-1];
		nextcov = chrom[i+1];
		if((cov != lastcov) && (cov == nextcov))
		{
			mi.x1 = i+1;
			
		}
		if((cov == lastcov) && (cov != nextcov))
		{
			mi.x2 = i+1;
			//mum.push_back(mi);
			//findPartnerCord(mi,mums,'R'); //mi is the current mem, mums is the original vector of mems
			findPartnerCord(mi,cp.mums,'R');
			vd = getCoverage(mi,masterRef,masterQ);//reset mi.y1, and mi.y2 as 1 because we won't use them here
			if(vd[0] >1)
			{
				cp.cc.push_back(mi);
			}
			if(vd[0] == 1)
			{
				cp.cm.push_back(mi);
			}
			if((vd[0] <0.01) && (vd[1] >0.1))
			{
				cp.in.push_back(mi);
			}
			if((vd[0] > 0.1) && (vd[1] <0.01))
			{
				cp.del.push_back(mi);
			}	
				
//cout<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t"<<vd[0]<<"\t"<<vd[1]<<endl;
		}
//cout<<i<<"\t"<<chrom[i]<<endl;
	}
	
//return mum;	
}						
/////////////////////////////////////////////////
void findPartnerCord(mI & mi, vector<mI> & mums,char c)//finds corresponding query or reference coordinates for supplied reference or query					
{
	unsigned int i =0;
	int span =0, endDiff =0, ovl =0;
	
	if(c == 'R')//provided reference coordinates
	{
		while(i<mums.size()-1) //when no overlap is found
		{
			endDiff = abs(min(mums[i].x1,mi.x1) - max(mums[i].x2,mi.x2));//distance between the two ends of the mems
			span = abs(mums[i].x2 - mums[i].x1) + abs(mi.x2- mi.x1); // total length of the two mems
						
			if((span -endDiff) > ovl) //if the overlap is greater than the overlap stored in the system
			{
				if(mums[i].y1 < mums[i].y2)//the query is one the same strand as the 
				{
					mi.y1 = mums[i].y1 + abs(mums[i].x1 - max(mums[i].x1,mi.x1));
					mi.y2 = mums[i].y2 - (mums[i].x2 -  min(mums[i].x2,mi.x2));
				}
				if(mums[i].y1 > mums[i].y2)//the query is on the reverse strand
				{
					mi.y1 = mums[i].y1 - abs(mums[i].x1 - max(mums[i].x1,mi.x1));
					mi.y2 = mums[i].y2 + (mums[i].x2 - max(mums[i].x1,mi.x1));
				} 
			}
			i++;
		}
	}
}		
