#ifndef PAIRING_H
#define PAIRING_H
#include "ecc.h"

class pairing{
public:
	pairing(){

		r=apint(5);
		
//		r = apint(11);
		h = apint(12);
		k = apint(2);
		q = modulo_arithmic(MODULO);
//		qk1r = apint(1560);
		qk1r=apint(696);
	}
	point dist_map(const point &obj){
		point res;
		res.x = obj.x.neg();

		complex i(apint(0), apint(1), 1);

		res.y = i*obj.y;
		return res;
	}

	complex line(const point &M, const point &N,const point &Q);
	complex bi(const point &a, const  point &b);
	complex tr(const point &P, const  point &Q);

	apint r;
	apint k;
	apint h;	
	modulo_arithmic q;
	apint qk1r;
//preparation stuff for Tate pairing, miller loop
};

//other class definition for extension field consisting Tate pairing

#endif
