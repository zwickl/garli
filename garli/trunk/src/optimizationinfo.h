#ifndef OPTIMIZATION_INFO
#define OPTIMIZATION_INFO

typedef pair<FLOAT_TYPE, FLOAT_TYPE> pd;

bool IsEvalLess(pd lhs, pd rhs){
	return lhs.second < rhs.second;
	}


class OptimizationInfo{
	
	int node;
	FLOAT_TYPE initLen;
	FLOAT_TYPE precision;
	bool goodGuess;
	
	FLOAT_TYPE initBracket[3];
	
	vector<pd> brakEvals;
	vector<pd> brentEvals;
	
	public:
	OptimizationInfo(){};
	void Setup(int n, FLOAT_TYPE len, FLOAT_TYPE prec, bool gg, FLOAT_TYPE low, FLOAT_TYPE mid, FLOAT_TYPE high){
		
		brakEvals.clear();
		brentEvals.clear();
		
		node=n;
		initLen=len;
		precision=prec;
		goodGuess=gg;
		
		initBracket[0]=low;
		initBracket[1]=mid;
		initBracket[2]=high;
		}
	void Report(ofstream &out){
		out.precision(12);
		out << "node\t" << node << "\tlen\t" << initLen << "\tprecision\t" << precision;
		if(goodGuess==true) out << "\t(good guess)\n";
		else out << "\t(weak guess)\n";
		out << "init Bracket\t" << initBracket[0] << "\t" << initBracket[1] << "\t" << initBracket[2] << "\n";
		out << "brak";
		for(vector<pd>::iterator it=brakEvals.begin();it!=brakEvals.end();it++){
			out << "\t"<< (*it).first << "\t"  << (*it).second << "\n";
			}
		if(brentEvals.empty() == false)		{
			out << "brent";
			for(vector<pd>::iterator it=brentEvals.begin();it!=brentEvals.end();it++){
				out << "\t"<< (*it).first << "\t" << (*it).second << "\n";
				}		
			}
		}
	void BrakAdd(FLOAT_TYPE val, FLOAT_TYPE score){
		brakEvals.push_back(make_pair(val, score));
		}

	void BrentAdd(FLOAT_TYPE val, FLOAT_TYPE score){
		brentEvals.push_back(make_pair(val, score));
		}
	
	bool IsMinAtMinAllowableLength(){
		vector<pd>::iterator minEval = min_element(brakEvals.begin(),brakEvals.end(), IsEvalLess);
		return FloatingPointEquals((*minEval).first, 0.01, 1e-10);
		}
	};
	

#endif

