#include<iostream>
#include "ml.h"

using namespace std;

using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;
int main(int argc, char *argv[])
{
	if((argc <3)||(string(argv[1]) == "-h")||(string(argv[1]) == "--help"))
	{
		cerr<<"Usage: "<<argv[0]<<" foo.delta foo.bed mode(s/n) ref.fasta query.fasta"<<endl;
		exit(EXIT_FAILURE);
	}

	chroms allChrom;
	char mode = 'n';
	mode = argv[3][0];
	map<string,string> refFasta;
	map<string,string> qFasta;
	map<string,vector<string> > cp; //cp is an alias for Chromosome partner. Each reference name index has a vector of unqiue alignments which are part of these
	map<string,vector<int> > seqLen;//length of sequences.first element is ref and second is query
	map<string,bool> qStrand; //stores whether query strand is forward strand or reverse strand
	mI tempmi;
	string line, chromName,refName,qName,indexAln;
	int refStart = 0, refEnd = 0, qStart = 0, qEnd = 0,refLen =0,qLen=0,  count = -1,indelPos =0, insPos=0;
	vector<int> vi;
	vector<string> vstr;
	vector<mI> cm,vmi,sameRef;
	vector<string> uniqChrom,uniqRef;
	size_t pos1,pos2,namePos;
	
	ifstream fin;
	ofstream fout;
	fin.open(argv[2]);//open the bed file with intervals
	while(getline(fin,line))
	{
		vstr = splitField(line,'\t');
		tempmi.rn = vstr[0];
		tempmi.x1 = stoi(vstr[1]);
		tempmi.x2 = stoi(vstr[2]);
		cm.push_back(tempmi);
		if(find(uniqChrom.begin(),uniqChrom.end(),tempmi.rn) == uniqChrom.end())
		{
			uniqChrom.push_back(tempmi.rn);
		}
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
			pos1 = line.find(' ',pos2+1); //recycling pos1 to find pos3
			refLen = stoi(line.substr(pos2+1,pos1-pos2));//reference length
			qLen = stoi(line.substr(pos1));//from last space till end 
			indexAln = refName + qName;
			count = -1;	
			cp[refName].push_back(indexAln); //adding the alignment to the list of refName alignments
		}
		if((line.size() <10) && (refName != "") && (count > -1))
		{
			indelPos = abs(stoi(line));
			insPos = stoi(line);
//cout<<refName<<'\t'<<refStart<<'\t'<<indelPos<<'\t'<<insPos<<'\t'<<qName<<'\t'<<qStart<<endl;
		//	refStart = refStart + indelPos;
			if(insPos < -1)//deletion starts
			{	
				refStart = refStart + abs(insPos)-1;
				refStart = (refStart * (-1));
				vi.push_back(refStart);
//cout<<"ins start\t"<<refName<<'\t'<<refStart<<'\t'<<indelPos<<'\t'<<insPos<<'\t'<<qName<<'\t'<<qStart<<endl;
			}
			if(insPos == -1)//increments of deletion
			{
				vi.push_back(-1);
//cout<<"ins\t"<<refName<<'\t'<<refStart<<'\t'<<"-1"<<'\t'<<insPos<<'\t'<<qName<<'\t'<<qStart<<endl;
			}
			if(insPos > 0)//0 or positive
			{
				refStart = refStart + insPos;
				vi.push_back(refStart);
//cout<<"del\t"<<refName<<'\t'<<refStart<<'\t'<<indelPos<<'\t'<<insPos<<'\t'<<qName<<'\t'<<qStart<<endl;
			}
//cout<<refName<<'\t'<<refStart<<'\t'<<indelPos<<'\t'<<insPos<<'\t'<<qName<<'\t'<<qStart<<endl;
			insPos = 0;
			refStart = abs(refStart);
			if(indelPos ==0) //reached the end of the indel description
			{
				tempmi.mv = vi;
				allChrom[indexAln].mums.push_back(tempmi);
				tempmi.mv.clear();
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
				refStart = refStart -1;//because alignment is 1 based
	//			--refStart;//to count the mutation distance
		}
	}
	fin.close();
	fin.open(argv[4]);
	readFasta(fin,refFasta);
	fin.close();
	fin.open(argv[5]);
	readFasta(fin,qFasta);
	fin.close();
	fout.open("lifted.txt");
	fout<<"REF_CHROM\tREF_START\tREF_END\tQ_CHROM\tQ_START\tQ_END"<<endl;
	//for(unsigned int i= 0;i<cm.size();i++) // cm is the user provided list of intervals
	sort(cm.begin(),cm.end());
	for(unsigned int i=0;i<uniqChrom.size();i++)
	{
		refName = uniqChrom[i];
	//	refName = cm[i].rn;
		for(unsigned int k =0;k<cm.size();k++)
		{
			if(cm[k].rn == refName)//if same reference names
			{
				sameRef.push_back(cm[k]);
			}
		}
//		sort(sameRef.begin(),sameRef.end());
		for(unsigned int j=0; j<cp[refName].size();j++)
		{
			indexAln = cp[refName][j];
			sort(allChrom[indexAln].mums.begin(),allChrom[indexAln].mums.end());
//			sort(sameRef.begin(),sameRef.end());
			findMum(allChrom[indexAln].mums,sameRef,fout,mode,refFasta,qFasta);//find the mum corresponding to cm
			//writeLift(vmi,cm[i],fout,mode,refFasta,qFasta);
//			vmi.clear();
		}
		sameRef.clear();
	}
	fout.close();
	return 0;
}
			


