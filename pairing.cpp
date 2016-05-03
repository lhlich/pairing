#include "header/pairing.h"
//#include <iostream>
//using namespace std;

complex pairing::bi(const point &a,const  point &b){
	point temp;
	temp=dist_map(b);
	complex res=tr(a,temp);

	return res.pow(qk1r);
}

complex pairing::bi(const  point &b,pairing_pp &pp){
	point temp;
	temp = dist_map(b);
	complex res = tr(temp,pp);

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

void pairing::line(const point &M, const point &N, complex &lambda, complex &niu, bool &vert){//line passing A,B, value fitting Q.

	complex three(apint(3), apint(0), 1);
	complex two(apint(2), apint(0), 1);
	complex paraA(apint(A), apint(0), 1);
	vert = 0;

	if (M == N){//tangent
		if ((two*M.y).zero()){//tangent is vertical
			vert = 1;
			lambda = M.x;
			return;
		}

		lambda = (three*M.x*M.x + paraA) / (two*M.y);
		niu = M.y - lambda*M.x;
	}
	else if (M.x == M.x){//vertical
		vert = 1;
		lambda = M.x;
		return;
	}
	else{
		lambda = (M.y - N.y) / (M.x - N.x);
		niu = M.y - lambda*M.x;
	}
}

complex eval_line(const linetoken &L, const point& Q){
	if (!L.vert)
		return (L.lambda)*Q.x + L.niu - Q.y;
	else
		return Q.x - L.lambda;
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


complex pairing::tr(const point &Q, pairing_pp &pp){
	complex x;
	x.a = apint(1);
	x.ext = 1;

	Node *temp=pp.head->next;


	for (int i = 2; i <= r.bits(); i++){

		x = x*x*eval_line(temp->reg_n, Q) / eval_line(temp->reg_d, Q);


		if (r.checkbit(i)){
			x = x*eval_line(temp->step_n, Q) / eval_line(temp->step_d, Q);

		}
		temp = temp->next;
	}
	return x;
}

pairing_pp::pairing_pp(pairing &E, const point &P){

	point Z(P);
	linetoken rn,rd,sn,sd;
	complex l, n;
	bool v;

	head = new Node(rn,rd,sn,sd);
	Node *temp = head;
	Node * token;

	for (int i = 2; i <= E.r.bits(); i++){

		E.line(Z, Z, l, n, v);

		rn = linetoken(l, n, v);


		E.line(2 * Z, (2 * Z).inv(), l, n, v);
		rd = linetoken(l, n, v);

		Z = 2 * Z;

		if (E.r.checkbit(i)){
			E.line(Z, P, l, n, v);

			sn = linetoken(l, n, v);


			E.line( Z+P, (Z+P).inv(), l, n, v);
			sd = linetoken(l, n, v);
			Z = Z + P;
		}

		token=new Node(rn, rd, sn, sd);
		temp->next = token;
		temp = temp->next;
	}
}


void pairing::load(){

	//cout<<"loading hval and rval...";
	const char*	hval="12016012264891146079388821366740534204802954401251311822919615131047207289359704531102844802183906537786776";//len 107
	const char* rval="730750818665451621361119245571504901405976559617"; //len 48
	
	//cout<<"Done!"<<endl;
	
	h=apint(0);
	r=apint(0);
	apint ten(10);
	int i;

	//cout<<"importing hval ...";
	for(i=0;i<107;i++){
		h=h*ten;
		h=h+apint(hval[i]-'0');
	}
	//cout<<"Done!"<<endl;

	//cout<<"importing rval ...";
	for(i=0;i<48;i++){
	//	cout<<"the "<<i<<"th bit"<<endl;
		r=r*ten;
		r=r+apint(rval[i]-'0');
	}
	//cout<<"Done!"<<endl;

	//cout<<"calculating qk1r...";
	apint qq=getQ().p;
	qk1r=(qq-apint(1))*h;
//	cout<<"Done!"<<endl;

	q=getQ();
	
}
