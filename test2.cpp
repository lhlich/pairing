#include "header/pairing.h"
#include <iostream>
#define HEXMODE 0
using namespace std;

void makepoint(point &P){


	P.y = P.x.pow(apint(3)) + P.x;
	P.y = P.y.pow(apint((MODULO + 1) / 4));

	complex one(apint(1), apint(0), 1);

	if (!P.on_curve()){
		P.x = P.x + one;
		makepoint(P);
	}

}

ostream& operator<<(ostream& os, const apint & obj)
{
	if (HEXMODE)os << "0x";
	for (int i = 0; i<obj.length; i++){
		if (HEXMODE)os << hex;
		os << obj.data[i];
		if (i != obj.length - 1)
			os << ',';
	}
	return os;
}

ostream& operator<<(ostream& os, const complex& obj)
{
	os << obj.a;

	if (obj.ext && !(obj.b.zero()))
		os << '+' << obj.b << 'i';

	return os;
}

ostream& operator<<(ostream& os, const point & obj)
{
	if (obj.id){
		os << "O";
		return os;
	}
	os << '[' << obj.x << ',' << obj.y << ']';

	return os;
}

void test_pairing();
int main(){
	test_pairing();
	return 0;
}



void test_pairing(){
	pairing E;
	point P;

	P.x = complex(apint(30), apint(0), 1);
	makepoint(P);

	P = E.h * P;

	complex one(apint(1), apint(0), 1);


	while (!P.on_curve()){

		P.x = P.x + one;
		makepoint(P);
		P = E.h * P;
	}

	point Q(2, 2);
	makepoint(Q);
	Q=E.h*Q;
	while (!Q.on_curve()){

		Q.x = Q.x + one;
		makepoint(Q);
		Q = E.h * Q;
	}

	complex res;

cout<<"------Evaluating pairing 1000 times-----"<<endl;
	for(int i=0;i<1000;i++)
		res=E.bi(P,Q);
cout<<"Done!"<<endl;


/*
	cout << "M:" << M << endl;
	cout << "N:" << N << endl;
	cout << "tr(M,N):" << E.tr(M,N) << endl;
	cout<<"Reduced:"<<E.tr(M,N).pow(E.qk1r)<<endl;
*/

}
