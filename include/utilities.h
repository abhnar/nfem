#pragma once
#include<iostream>
#include<vector>
#include<list>
#include<string>
#include<math.h>
#include "matrix.h"
using namespace std;


class Util {
public: 
	static double integratePlynomialxy(dmatrix<int>& poly, double x1, double x2, double x3, double y1, double y2, double y3) {
		//poly.print();
		double xw[6][3] = { {0.44594849091597 ,   0.44594849091597  ,  0.22338158967801},
							{0.44594849091597 ,   0.10810301816807  ,  0.22338158967801},
							{0.10810301816807  ,  0.44594849091597  ,  0.22338158967801},
							{0.09157621350977  ,  0.09157621350977   , 0.10995174365532},
							{0.09157621350977  ,  0.81684757298046  ,  0.10995174365532},
							{0.81684757298046  ,  0.09157621350977  ,  0.10995174365532} };
		double A = A = abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0;
		double z = 0, x, y;
		for (int j = 0; j < 6; j++) {
			x = x1 * (1 - xw[j][0] - xw[j][1]) + x2 * xw[j][0] + x3 * xw[j][1];
			y = y1 * (1 - xw[j][0] - xw[j][1]) + y2 * xw[j][0] + y3 * xw[j][1];
			// cout<<"y: "<<y<<endl;
			z = z + (poly(0, 0) + poly(0, 1) * x + poly(0, 2) * x * x + poly(0, 3) * x * x * x +
				poly(1, 0) + poly(1, 1) * y + poly(1, 2) * y * y + poly(1, 3) * y * y * y + poly(2, 0) * x * y) * xw[j][2];
			// z=z+poly(0,0)*1.0;
		   //  cout<<x <<" "<<y<<" "<<z<<endl;
		}
		//cout<<z<<endl;
		return z * A;

	}
	static  int binarySearch(int arr[], int l, int r, int x)
	{
		if (r >= l)
		{
			int mid = l + (r - l) / 2;
			if (arr[mid] == x)
				return mid;
			if (arr[mid] > x)
				return binarySearch(arr, l, mid - 1, x);
			return binarySearch(arr, mid + 1, r, x);
		}
		return -1;
	}


	// trim from start (copying)
	static inline std::string ltrim(std::string s) {
		s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
			return !std::isspace(ch);
			}));
		return s;
	}
	static void apendToFile(string line) {
		ofstream outfile;

		outfile.open("d:\\test.txt", std::ios_base::app);
		outfile << line << endl;
		outfile.close();
	}
	// trim from end (copying)
	static inline std::string rtrim(std::string s) {
		s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
			return !std::isspace(ch);
			}).base(), s.end());
		return s;
	}

	// trim from both ends (copying)
	static inline std::string trim(std::string s) {
		s = ltrim(s);
		s = rtrim(s);
		return s;
	}
	static vector<string> split(const string& str, const string& delim)
	{
		vector<string> tokens;
		size_t prev = 0, pos = 0;
		do
		{
			pos = str.find(delim, prev);
			if (pos == string::npos) pos = str.length();
			string token = str.substr(prev, pos - prev);
			if (!token.empty()) tokens.push_back(token);
			prev = pos + delim.length();
		} while (pos < str.length() && prev < str.length());
		return tokens;
	}
};

class  pipeline{
    vector <string> pipe;
    int sz;
    int pos;
    
    public:
    pipeline(int s){
        sz=s;
        for(int i=0;i<sz;i++){
            pipe.push_back(" ");
        }
        pos=-1;
    }
    void push(string item){
        pos=pos+1;
        if(pos==sz) pos=0;
        pipe[pos]=item;
       
        
    }
    string fetch(){
        int idx;
        if(pos==sz-1){
            idx=0;
        }else{
            idx=pos+1;
        }
        return pipe[idx];
       
        
    }
    void print(){

        for(int i=0;i<sz;i++){
            cout<<pipe[i]<<endl;
            
        }
        cout<<endl;
    } 
    ~pipeline(){
        
    }
};
