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
	int ci = 0; //ci is minus i
	int found = 0;
	qord temp;
	bool oneFound = 0;//whether at least one if statement has been hit
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
		while(refC<mi.x2+1)
		{
//cout<<cm.rn<<'\t'<<cm.x1<<'\t'<<cm.x2<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mi.y1<<'\t'<<mi.y2<<endl;
			ci = refC * (-1);
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //if this position does not have a indel
			{
				refC++;
				qC++;
				temp.cord = qC-1;
				tempmi.x1 = refC-1;
				tempmi.x2 = refC-1;
				tempmi.y1 = qC-1;
				tempmi.y2 = qC-1;
				tempmi.qn = mi.qn;
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) != mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //position has insertion
			{
				refC++;
				tempmi.x1 = refC-1;
				tempmi.x2 = refC-1;
				tempmi.y1 = qC-1;
				tempmi.y2 = qC-1;
				tempmi.qn = mi.qn;
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())) //position has deletion
			{
				refC++;
				qC++;
				tempmi.x1 = refC-1;
				tempmi.x2 = refC-1;
				tempmi.y1 = qC-1;
				tempmi.y2 = qC-1;
				tempmi.qn = mi.qn;
				qC++;
				found = returnIndex(mi.mv,ci);//get the index
				mi.mv[found] = 0;//reset the deletion start
				found++;
				while((mi.mv[found] == -1) && (found<mi.mv.size()))
				{
					qC++;
					found++;
				}
				//refC++;//increment refC to move to the next	
			}
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
			ci = refC * (-1);
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //if this position does not have a indel
			{
				refC++;
				qC--;
				tempmi.x1 = refC-1;
				tempmi.x2 = refC-1;
				tempmi.y1 = qC+1;
				tempmi.y2 = qC+1;
				tempmi.qn = mi.qn;
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) != mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) == mi.mv.end())) //position has insertion
			{
				refC++;
				tempmi.x1 = refC-1;
				tempmi.x2 = refC-1;
				tempmi.y1 = qC+1;
				tempmi.y2 = qC+1;
				tempmi.qn = mi.qn;
			}
			if((find(mi.mv.begin(),mi.mv.end(),refC) == mi.mv.end()) && (find(mi.mv.begin(),mi.mv.end(),ci) != mi.mv.end())) //position has deletion
			{
				refC++;
				qC--;
				tempmi.x1 = refC-1;
				tempmi.x2 = refC-1;
				tempmi.y1 = qC+1;
				tempmi.y2 = qC+1;
				tempmi.qn = mi.qn;
				qC--;
				found = returnIndex(mi.mv,ci);//get the index
				mi.mv[0] = 0;//reset the deletion start
				found++;
				while((mi.mv[found] == -1) && (found<mi.mv.size()))
				{
					qC--;
					found++;
				}
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
int returnIndex(vector<int> & vi, int delPos)
{
	int i =0;
	while(( vi[i] != delPos) && (i<vi.size()))
	{
		i++;
	}
	return	i;
}
///////////////////////////////////////////////////
void findMum(vector<mI> & mums,vector<mI> & cm,ofstream & fout, char & c, map<string,string> & rseq, map<string,string> & qseq)
{
	unsigned int i =0;
	mI tempmi;
	//vector<mI> vmi;
	for(unsigned int j=0;j<cm.size();j++)
	{
		while((!(cm[j].x2 < mums[i].x1)) && (i<mums.size()))
		{
//cout<<"scan\t"<<cm[j].rn<<'\t'<<cm[j].x1<<'\t'<<cm[j].x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<'\t'<<i<<endl;			
			if(!(cm[j].x2 < mums[i].x1) && !(cm[j].x1 > mums[i].x2))
			{
				tempmi = mums[i];
				//vmi.push_back(tempmi);
				writeLift(tempmi,cm[j],fout,c,rseq,qseq);
//cout<<"query\t"<<cm[j].rn<<'\t'<<cm[j].x1<<'\t'<<cm[j].x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<'\t'<<i<<endl;
			}
			++i;
		}
		i = 0;//reset i
	}
//	return vmi;
}
//////////////////////////////////////////////////////
void writeLift(mI & vmi,mI & cm,ofstream & fout, char & c, map<string,string> & rseq, map<string,string> & qseq)
{
	mI tempmi,mi;
	if(c=='s')//if sensitive mode is used
	{
		cm.x1 = max(cm.x1,vmi.x1);
		cm.x2 = min(cm.x2,vmi.x2);
	}
	tempmi = liftCords(cm,vmi);
	if((tempmi.y1 != 0) && (tempmi.y2 !=0))
	{
		if(vmi.y1 < vmi.y2)//forward strand
		{
			fout<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<'\t'<<rseq[tempmi.rn].substr(tempmi.x1-1,max(tempmi.x2-tempmi.x1,1))<<'\t'<<qseq[tempmi.qn].substr(tempmi.y1-1,max(tempmi.y2-tempmi.y1,1))<<endl;
		}
		if(vmi.y1 > vmi.y2)//reverse strand
		{
			fout<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<'\t'<<rseq[tempmi.rn].substr(tempmi.x1-1,max(tempmi.x2-tempmi.x1,1))<<'\t'<<revcom(qseq[tempmi.qn].substr(tempmi.y2-1,max(tempmi.y1-tempmi.y2,1)))<<endl;
		}
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
/////////////////////////////////////////////////
void readFasta(ifstream & fin, map<string,string> & seq)
{
	string str,name;
	while(getline(fin,str))
	{
		if(str[0] == '>')//fasta header
		{
			name = str.substr(1);
		}
		else
		{
			seq[name].append(str);
		}
	}
}	
//////////////////////////////////////////////////
string revcom(string seq) //return the complementary sequence
{
	string revseq;
	char N,c;
	for(int i=seq.size()-1;i>-1;i--)
	{
		N = seq[i];
		if(N == 'N')
		{
			c = 'N';
		}
		if(N == 'n')
		{
			c = 'n';
		}
		if(N == 'A')
		{
			c = 'T';
		}
		if(N == 'a')
		{
			c = 't';
		}
		if(N == 'T')
		{
			c = 'A';
		}
		if(N == 't')
		{
			c = 'a';
		}
		if(N == 'G')
		{
			c = 'C';
		}
		if(N == 'g')
		{
			c = 'c';
		}	
		if(N == 'C')
		{
			c = 'G';
		}
		if(N == 'c')
		{
			c = 'g';
		}
		revseq.push_back(c);
	}
	return revseq;
}	
