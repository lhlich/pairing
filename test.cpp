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
void test_apint();
void test_ecc();
void test_pairing();
int main(){



	test_apint();
	test_ecc();
	test_pairing();


	return 0;
}


void test_apint(){
	apint x(234567890);
	apint y(456789012);

	modulo_arithmic P(123121);

	P.reduce(x);
	//P.reduce(y);

	cout << "---------apint test-----------" << endl;
	cout << "x=" << x << endl;
	cout << "y=" << y << endl;

	apint r = x + y;
	cout << "x+y=" << r << endl;

	r = x*y;
	cout << "x*y=" << r << endl;

	apint s = y*y;

	cout << "y^2-x*y=" << s - r << endl;
	cout << "(y-x)*y=" << (y - x)*y << endl;
	cout << "y<<56=" << (y << 56) << endl;

	y.bits();

	apint rem;
	apint q = y.div(apint(123121), rem);

	modulo_arithmic Q(26);
	x = Q.inv(apint(15));
}

void test_ecc(){
	complex x(apint(25), apint(0), 1);
	complex y(apint(30), apint(0), 1);

	point P(x, y);
	makepoint(P);


	point Q(y, y);
	makepoint(Q);

	bool t = P.on_curve();
	t = Q.on_curve();

	cout << "----------ECC test------------" << endl;
	cout << "P:" << P << endl;
	cout << "Q:" << Q << endl;
	cout << "-Q:" << Q.inv() << endl;
	cout << "P+Q:" << P + Q << endl;
	cout << "P-Q:" << P - Q << endl;
	cout << "P+(-Q):" << P + Q.inv() << endl;
	cout << "[2]P:" << (2 * P) << endl;
	cout << "[4]P:" << (4 * P) << endl;
	cout << "(P+Q)+(P-Q):" << (P + Q) + (P - Q) << endl;
	cout << "[3]P+[4]P:" << (3*P)+(4*P) << endl;
	cout << "[7]P:" << (7 * P) << endl;

	cout << "------------------------------" << endl;

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


	P.on_curve();
	Q.on_curve();


	complex res = E.bi(P, Q);

	cout << "----------pairing test------------" << endl;
	cout << "P:" << P << endl;
	cout << "Q:" << Q << endl;
	cout << "e(P,Q):" << res << endl;
	cout << "e(P,Q)^3:" << res*res*res << endl;
	cout << "[3]P:" << (3 * P) << endl;
	cout << "e([3]P,Q):" << E.bi(3 * P, Q) << endl;
	cout << "------------------------------" << endl;

	point M(25,30);
	point N(34,0);
	N.y.b=apint(30);

/*
	cout << "M:" << M << endl;
	cout << "N:" << N << endl;
	cout << "tr(M,N):" << E.tr(M,N) << endl;
	cout<<"Reduced:"<<E.tr(M,N).pow(E.qk1r)<<endl;
*/

}
