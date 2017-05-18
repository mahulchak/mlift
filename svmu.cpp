#include<iostream>
#include "sv.h"


using namespace std;

using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;
int main(int argc, char *argv[])
{
	if(argc <2)
	{
		cerr<<"Usage: "<<argv[0]<<" foo.delta"<<endl;
		exit(EXIT_FAILURE);
	}

	chroms allChrom;
	
	map<string,ccov> masterRef; //stores sequence coverage but it can also be used to find reference chromosome lengths
	map<string,ccov>masterQ; //stores sequence coverage but it can also be used to find query chromosome lengths
	map<string,vector<string> > cp; //cp is an alias for Chromosome partner. Each reference name index has a vector of unqiue alignments which are part of these
	map<string,vector<string> > hcp;//hcp stands for homologous cp
	map<string,map<int,vq> > mRef; //stores the coordinates of query on reference chromosomes
	
	mI tempmi,gapmi;

	string foo = string(argv[1]);
	string line, chromName,refName,qName,indexAln;
	int refStart = 0, refEnd = 0, qStart = 0, qEnd = 0, refLen =0, qLen =0, count = -1,indelPos =0;
	
	vector<double> vd;
	vector<int> vi;
	vector<mI> vmi;
	size_t pos1,pos2,namePos;
	ifstream fin;
	
	fin.open(argv[1]);
	
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
			cp[refName].push_back(indexAln); //adding the alignment to the list of refName alignments
			if(masterRef[refName].size() == 0)//if they have not been created
			{
				masterRef[refName] = makeChromBucket(refLen);
			}
			if(masterQ[qName].size() == 0)//if they have not been created
			{
				masterQ[qName] = makeChromBucket(qLen);
			}
//cout<<indexAln<<"\t"<<refName<<"\t"<<qName<<"\t"<<refLen<<"\t"<<qLen<<"\t"<<refName.size()<<"\t"<<qName.size()<<endl;			
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
				allChrom[indexAln].mums.push_back(tempmi);
				storeCords(masterRef[refName],masterQ[qName],tempmi);
				storeCords(mRef[refName],tempmi);
//cout<<refName<<"\t"<<refStart<<"\t"<<refEnd<<"\t"<<qName<<"\t"<<qStart<<"\t"<<qEnd<<"\t"<<allChrom[indexAln].mums.size()<<endl;
				vi.clear();//reset it once its values are used
			}
				
			count++;
			
		}
		if((line.find('>') == string::npos) && (line.size() >10) && (refName != "")) //when describing alignment segments
		{
//cout<<line<<endl;		
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

				//allChrom[indexAln].mums.push_back(tempmi);
				//storeCords(masterRef[refName],masterQ[qName],tempmi);
				
//cout<<refName<<"\t"<<refStart<<"\t"<<refEnd<<"\t"<<qName<<"\t"<<qStart<<"\t"<<qEnd<<"\t"<<allChrom[indexAln].mums.size()<<endl;
		}
	}
	fin.close();
//	splitByCoverage(masterRef[refName],allChrom[indexAln].mums,masterRef[refName],masterQ[qName]);
	for(chroms::iterator it = allChrom.begin();it!= allChrom.end();it++)
	{
		indexAln = it->first;
		sort(allChrom[indexAln].mums.begin(),allChrom[indexAln].mums.end());
		//for(int j =15000;j<20000;j++)
		//{
		//	cout<<"2L"<<"\t"<<j;
		//	for(unsigned int ct=0;ct<mRef["2L"][j].size();ct++)
		//	{
		//		cout<<"\t"<<mRef["2L"][j][ct].name<<"\t"<<mRef["2L"][j][ct].cord;
		//	}
		//	cout<<endl;
		//}
		for(unsigned int i= 0; i<allChrom[indexAln].mums.size();i++)
		{
			tempmi = allChrom[indexAln].mums[i];
			vd = getCoverage(tempmi,masterRef[tempmi.rn],masterQ[tempmi.qn]);
			if((vd[0] <1.3) && (vd[1]<1.3))
			{
				if(allChrom[indexAln].cm.size() == 0)
				{
					gapmi.rn = tempmi.rn;
					gapmi.qn = tempmi.qn;
					gapmi.x1 = 1;
					gapmi.x2 = tempmi.x2;
					gapmi.y1 = 1;
					gapmi.y2 = max(tempmi.y1,tempmi.y2);
					allChrom[indexAln].gap.push_back(gapmi);
				}
					
				if(allChrom[indexAln].cm.size() >0)//once more than one element has been entered
				{	
					refStart = allChrom[indexAln].cm[allChrom[indexAln].cm.size() -1].x2;
					qStart = max(allChrom[indexAln].cm[allChrom[indexAln].cm.size() -1].y2,allChrom[indexAln].cm[allChrom[indexAln].cm.size() -1].y1);
					if(refStart < tempmi.x1)
					{
						gapmi.rn = tempmi.rn;
						gapmi.qn = tempmi.qn;
						gapmi.x1 = refStart;
						gapmi.x2 = tempmi.x1;
						gapmi.y1 = qStart;
						gapmi.y2 = min(tempmi.y1,tempmi.y2);
						allChrom[indexAln].gap.push_back(gapmi);
				//cout<<"gap\t"<<indexAln<<"\t"<<gapmi.x1<<"\t"<<gapmi.x2<<"\t"<<gapmi.y1<<"\t"<<gapmi.y2<<"\t"<<vd[0]<<"\t"<<vd[1]<<endl;
					}
				}
				allChrom[indexAln].cm.push_back(tempmi);
				count = count + tempmi.x2 - tempmi.x1; //keeping a count of total unique alignment
				//cout<<"cm\t"<<indexAln<<"\t"<<tempmi.x1<<"\t"<<tempmi.x2<<"\t"<<tempmi.y1<<"\t"<<tempmi.y2<<"\t"<<vd[0]<<"\t"<<vd[1]<<endl;
			}
			//if(!(vd[0] <1.3) && !(vd[1]<1.3))	
			else
			{
				 allChrom[indexAln].ncm.push_back(tempmi);
				//cout<<"ncm\t"<<indexAln<<"\t"<<tempmi.x1<<"\t"<<tempmi.x2<<"\t"<<tempmi.y1<<"\t"<<tempmi.y2<<"\t"<<vd[0]<<"\t"<<vd[1]<<endl;
			}
		}
//cout<<"size of cm is "<<allChrom[indexAln].cm.size()<<endl;
//cout<<"Beginning gap filling for "<<indexAln<<endl;
		if(allChrom[indexAln].cm.size()>10) //if less than 10 unique regions map to an alignment
		{
			hcp[allChrom[indexAln].cm[0].rn].push_back(indexAln);//homologous alignment
			
			for(unsigned int i=0;i<allChrom[indexAln].gap.size();i++)
			{
				tempmi = allChrom[indexAln].gap[i];
//				gapCloser(allChrom[indexAln].gap[i], allChrom[indexAln].ncm, allChrom[indexAln].cm);
				gapCloser(tempmi, allChrom[indexAln].ncm, allChrom[indexAln].cm);
//cout<<"up "<<allChrom[indexAln].gap.size()<<" "<<i<<" "<<allChrom[indexAln].gap[i].x1<<" "<<allChrom[indexAln].gap[i].x2<<" "<<allChrom[indexAln].gap[i].y1<<" "<<allChrom[indexAln].gap[i].y2<<endl;

			}
		}
		allChrom[indexAln].mums.clear();
//cout<<"size of cm is "<<allChrom[indexAln].cm.size()<<endl;	 
		count = 0; //reset count for the next alignment
	}
//cout<<"Done with gap filling "<<endl;	
	//for(map<string,vector<string> >::iterator it = cp.begin(); it != cp.end();it++)
	for(map<string,vector<string> >::iterator it = hcp.begin(); it != hcp.end();it++)
	{
		refName = it->first;
		//for(unsigned int i =0; i<cp[refName].size(); i++)
		for(unsigned int i = 0; i<hcp[refName].size();i++)
		{		
			//indexAln = cp[refName][i];
			indexAln = hcp[refName][i];
			qName = allChrom[indexAln].mums[i].qn;
			sort(allChrom[indexAln].cm.begin(),allChrom[indexAln].cm.end());
			splitByCoverage(allChrom[indexAln],masterRef[refName],masterQ[qName]);
//cout<<"Completed splitting coverage for "<<refName<<endl;
//cout<<"came here"<<endl;
//cout<<"size of cm is "<<allChrom[indexAln].cm.size()<<endl;
//cout<<"size of cm is "<<allChrom[indexAln].cm.size()<<endl;
			for(unsigned int j=0;j<allChrom[indexAln].cm.size();j++)
			{
				tempmi = allChrom[indexAln].cm[j];
				//cout<<tempmi.rn<<" "<<tempmi.x1<<" "<<tempmi.x2<<" "<<tempmi.qn<<" "<<tempmi.y1<<" "<<tempmi.y2<<endl;
			}
			for(unsigned int j=0; j<allChrom[indexAln].cc.size();j++)
			{
				vmi = findQuery(mRef[refName],allChrom[indexAln].cc[j],masterRef[refName],masterQ[qName]);
				if(vmi.size() >0) //if it is not empty
				{
					for(unsigned int k=0;k<vmi.size();k++)
					{
						cout<<vmi[k].rn<<"\t"<<vmi[k].x1<<"\t"<<vmi[k].x2<<"\t"<<vmi[k].qn<<"\t"<<vmi[k].y1<<"\t"<<vmi[k].y2<<endl;
					}
				}
			}
		}
	}

		


return 0;
}
			


