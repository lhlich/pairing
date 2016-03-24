#include "header/pairing.h"
//#include <iostream>
//using namespace std;

complex pairing::bi(const point &a,const  point &b){
	point temp;
	temp=dist_map(b);
	complex res=tr(a,temp);

	return res.pow(qk1r);
}
complex pairing::line(const point &M, const point &N,const point &Q){//line passing A,B, value fitting Q.
	
	complex three(apint(3), apint(0), 1);
	complex two(apint(2), apint(0), 1);
	complex paraA(apint(A), apint(0), 1);
	complex lambda, niu;	

	if(M==N){//tangent
		if((two*M.y).zero())//tangent is vertical
			return Q.x-M.x;

		lambda = (three*M.x*M.x + paraA) / (two*M.y);
		niu = M.y - lambda*M.x;
	}
	else if(M.x==M.x){//vertical
		return Q.x-M.x;
	}
	else{
		lambda = (M.y - N.y) / (M.x -N.x);
		niu = M.y - lambda*M.x;		
	}
	return lambda*Q.x + niu -Q.y;

}
complex pairing::tr(const point &P,const point &Q){
	complex x;
	x.a = apint(1);
	x.ext = 1;

	bool vert=0;//vertical between Z and P

	point Z(P);

	for(int i = 2; i <= r.bits(); i++){
		
		x=x*x*line(Z,Z,Q)/line(2*Z,(2*Z).inv(),Q);
		Z = 2*Z;

		if (r.checkbit(i)){
			x=x*line(Z,P,Q)/line(Z+P,(Z+P).inv(),Q);
			Z=Z+P;
		}
	}
	return x;
}