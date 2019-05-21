#ifndef ML_H_
#define ML_H_

#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<vector>
#include<map>
#include<algorithm>
#include<list>
#include<climits>
#include<cmath>
#include<utility>
#include<iomanip>

using namespace std;

//to store coordinates at alignment level
struct mI {
	string rn;
	string qn;
        int x1;//reference start
        int x2;//reference end
        int y1;//query start
        int y2;//query end
        char c;//qualifier. special comments or information can be added to this:i=inversion;
	vector<int> mv;        
        bool operator < (const mI& mum1) const
        {
                //return(x1 < mum1.x1) || ((x1 == mum1.x1) && (x2 < mum1.x2));
                return(x2 < mum1.x2) || ((x2 == mum1.x2) && (x1 < mum1.x1));
        }
        bool operator == (const mI& mum1) const
        {
                return x1 == mum1.x1 && x2 == mum1.x2 && y1 == mum1.y1 && y2 == mum1.y2;
        }
        };
//to store coordinates at base pair level
struct qord {
	string name;
	int cord;
	bool operator < (const qord& q1) const
	{
		return(name<q1.name) || ((name == q1.name) && (cord<q1.cord));
	}
	};

class chromPair {
	public:
	vector<mI> mums;	
	vector<mI> cm; //conserved mems
	vector<mI> ncm; //conserved mems from reverse side
	vector<mI> gap; //gaps are represented as mums
	vector<mI> cc; //cnv candidates
//	vector<mI> in;//stores insertion mums in reference
//	vector<mI> del; //stores deletion mums in query
};

bool qusort(mI mi1, mI mi2); //to sort the mI based on query coordinates
mI liftCords(mI & cm,mI & mi); //extract the lifted coordinates
vector<mI> findMum(vector<mI> & mums,mI & cm);
void writeLift(vector<mI> & vmi,mI cm,ofstream & fout, char & c);
vector<string> splitField(string & str, char c);
#endif
