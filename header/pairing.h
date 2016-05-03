#ifndef PAIRING_H
#define PAIRING_H
#include "ecc.h"


struct linetoken{
	linetoken(const complex &l,const complex &n,const bool &v){
		lambda = l;
		niu = n;
		vert = v;
	}
	linetoken(){
		lambda = complex();
		niu = complex();
		vert = 0;
	}

	complex lambda;
	complex niu;
	bool vert;
};

struct Node{
	Node(const linetoken rn, const linetoken rd, const linetoken sn, const linetoken sd){
		reg_n = rn;
		reg_d = rd;
		step_n = sn;
		step_d = sd;
		next = nullptr;
	}

	linetoken reg_n;
	linetoken reg_d;
	linetoken step_n;
	linetoken step_d;

	Node *next;
};
class pairing_pp;

complex eval_line(const linetoken &L, const point& Q);

class pairing{
public:
	pairing(){
		
		r=apint(5);
		
//		r = apint(11);
		h = apint(12);
		k = apint(2);
//		q = modulo_arithmic(MODULO);
//		qk1r = apint(1560);
		qk1r=apint(696);

		if(!LIGHT) load();
	}
	point dist_map(const point &obj){
		point res;
		res.x = obj.x.neg();

		complex i(apint(0), apint(1), 1);

		res.y = i*obj.y;
		return res;
	}

	complex line(const point &M, const point &N,const point &Q);

	void line(const point &M, const point &N, complex &LAMBDA, complex &NIU, bool &vert);
	complex bi(const point &a, const  point &b);
	complex bi(const point &Q, pairing_pp &pp);
	complex tr(const point &P, const  point &Q);
	complex tr(const point &Q, pairing_pp &pp);
	void load();//load huge parameter

	apint r;
	apint k;
	apint h;	
	modulo_arithmic q;
	apint qk1r;
//preparation stuff for Tate pairing, miller loop
};
struct pairing_pp{
	pairing_pp(pairing &E,const point &P);

	Node *head;
};

//other class definition for extension field consisting Tate pairing

#endif
