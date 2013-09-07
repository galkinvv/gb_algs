/**
\file
Реализация не-inline методов CModular
*/

#include "cmodular.h"
#include <cstdio>
#include <cstring>
using namespace std;

namespace F4MPI{
int CModular::MOD = 2;
bool CModular::NEEDMUL64 = true;

int CModular::extendedEuclid(int a, int b, int& x, int& y){
	if (a < b)
		return extendedEuclid(b, a, y, x);
	else if (b == 0)
	{
		x = 1;
		y = 0;
		return a;
	}
	else
	{            
		int dd = extendedEuclid(b, a % b, y, x);
		y = y - (a / b) * x;
		return dd;
	}
}

string CModular::toString() const {
	 char buf[20];
	 sprintf(buf, "%d", val);	 
	 string tmp = buf;
	 return tmp;	 
}
 
void CModular::setValue(string strValue){
	int tmp;
	if(sscanf(strValue.c_str(), "%d", &tmp)==1){
		val = tmp % MOD;
		while(val<0)val+=MOD;
	}	 
}
} //namespace F4MPI
