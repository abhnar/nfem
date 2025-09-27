#define _HAS_STD_BYTE 0

#include<iostream>
#include "matrix.h"
#include "EFEM.h"


using namespace std;


int main(){
	

	EFEM fem;
	fem.load("Multi_IRIS");
	fem.start(4);

	


	//SetConsoleTextAttribute(hConsole, 23);
	cout<<"Done"<<endl;
	
	
	return 0;
}

