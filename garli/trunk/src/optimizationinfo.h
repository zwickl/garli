#ifndef OPTIMIZATION_INFO
#define OPTIMIZATION_INFO

typedef pair<double, double> pd;

bool IsEvalLess(pd lhs, pd rhs){
	return lhs.second < rhs.second;
	}


class OptimizationInfo{
	
	int node;
	double initLen;
	double precision;
	bool goodGuess;
	
	double initBracket[3];
	
	vector<pd> brakEvals;
	vector<pd> brentEvals;
	
	public:
	OptimizationInfo(){};
	void Setup(int n, double len, double prec, bool gg, double low, double mid, double high){
		
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
	void BrakAdd(double val, double score){
		brakEvals.push_back(make_pair(val, score));
		}

	void BrentAdd(double val, double score){
		brentEvals.push_back(make_pair(val, score));
		}
	
	bool IsMinAtMinAllowableLength(){
		vector<pd>::iterator minEval = min_element(brakEvals.begin(),brakEvals.end(), IsEvalLess);
		return (*minEval).first == 0.01;
		}
	};
	

#endif

