#include<iostream>
#include "ml.h"

using namespace std;

using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;
int main(int argc, char *argv[])
{
	if(argc <2)
	{
		cerr<<"Usage: "<<argv[0]<<" foo.delta foo.bed"<<endl;
		exit(EXIT_FAILURE);
	}

	chroms allChrom;
	
	map<string,ccov> masterRef; //stores sequence coverage but it can also be used to find reference chromosome lengths
	map<string,ccov>masterQ; //stores sequence coverage but it can also be used to find query chromosome lengths
	map<string,ccov>masterHQ; //same as masterQ but records coverage only for homologous pairs
	map<string,vector<string> > cp; //cp is an alias for Chromosome partner. Each reference name index has a vector of unqiue alignments which are part of these
		
	map<string,map<int,vq> > mRef; //stores the coordinates of query on reference chromosomes
	map<string,vector<int> > seqLen;//length of sequences.first element is ref and second is query
	map<string,bool> qStrand; //stores whether query strand is forward strand or reverse strand
	mI tempmi;
	string line, chromName,refName,qName,indexAln;
	int refStart = 0, refEnd = 0, qStart = 0, qEnd = 0, refLen =0, qLen =0, count = -1,indelPos =0;
	vector<int> vi;
	vector<string> vstr;
	vector<mI> cm,vmi;
	size_t pos1,pos2,namePos;
	
	ifstream fin;
	ofstream fout;
	fin.open(argv[2]);//open the bed file with intervals
	while(getline(fin,line))
	{
		vstr = splitField(line,'\t');
//		cout<<vstr[0]<<"\t"<<vstr[1]<<"\t"<<vstr[2]<<endl;
		tempmi.rn = vstr[0];
		tempmi.x1 = stoi(vstr[1]);
		tempmi.x2 = stoi(vstr[2]);
		cm.push_back(tempmi);
//cout<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<endl;		
	}
	fin.close();
	sort(cm.begin(),cm.end());
	fin.open(argv[1]);	
	while(getline(fin,line))
	{
		if(line.find('>') != string::npos)//start of an aligning chromosome description
		{			
			refName = line.substr(1,line.find(' ')-1);
			pos1 = line.find(' '); //position of the first space
			pos2 = line.find(' ',pos1+1);//position of the second space
			qName = line.substr(pos1+1, pos2-pos1-1); //up to the second space
//cout<<qName<<endl;
			pos1 = line.find(' ',pos2+1); //recycling pos1 to find pos3
			refLen = stoi(line.substr(pos2+1,pos1-pos2));//reference length
			qLen = stoi(line.substr(pos1));//from last space till end 
			indexAln = refName + qName;
			count = -1;
			//seqLen[indexAln].push_back(refLen);
			//seqLen[indexAln].push_back(qLen);
			cp[refName].push_back(indexAln); //adding the alignment to the list of refName alignments
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
//cout<<refName<<"\t"<<indelPos<<" " <<refStart<<"\t"<<refEnd<<"\t"<<qName<<"\t"<<qStart<<"\t"<<qEnd<<endl;
			if(indelPos ==0) //reached the end of the indel description
			{
				tempmi.mv = vi;
				allChrom[indexAln].mums.push_back(tempmi);
				tempmi.mv.clear();
//cout<<refName<<"\t"<<refStart<<"\t"<<refEnd<<"\t"<<qName<<"\t"<<qStart<<"\t"<<qEnd<<"\t"<<allChrom[indexAln].mums.size()<<endl;
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
	//			--refStart;//to count the mutation distance
		}
	}
	fin.close();

	fout.open("lifted.txt");
	fout<<"REF_CHROM\tREF_START\tREF_END\tQ_CHROM\tQ_START\tQ_END"<<endl;
	for(unsigned int i= 0;i<cm.size();i++) // cm is the user provided list of intervals
	{
		refName = cm[i].rn;
		for(unsigned int j=0; j<cp[refName].size();j++)
		{
			indexAln = cp[refName][j];
//cout<<refName<<'\t'<<cm[i].x1<<'\t'<<cm[i].x2<<'\t'<<indexAln<<'\t'<<allChrom[indexAln].mums[0].x1<<'\t'<<allChrom[indexAln].mums[0].x2<<'\t'<<allChrom[indexAln].mums[allChrom[indexAln].mums.size()-1].x1<<'\t'<<allChrom[indexAln].mums[allChrom[indexAln].mums.size()-1].x2<<endl;
			sort(allChrom[indexAln].mums.begin(),allChrom[indexAln].mums.end());
			vmi = findMum(allChrom[indexAln].mums,cm[i]);//find the mum corresponding to cm
			writeLift(vmi,cm[i],fout);
		}
	}
	fout.close();
	return 0;
}
			


