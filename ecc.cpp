#include "header/ecc.h"


complex complex::operator+(const complex  &right){
	modulo_arithmic Q = getQ();

	complex res;

	res.a=Q.add(a, right.a);
	res.b = Q.add(b, right.b);
	res.ext = ext;

	return res;
}
complex complex::operator-(const complex &right) const{
	modulo_arithmic Q = getQ();

	complex res;

	res.a = Q.sub(a, right.a);
	res.b = Q.sub(b, right.b);
	res.ext = ext;

	return res;
}

complex complex::operator*(const complex &right) const{
	modulo_arithmic Q = getQ();

	complex res;

	apint temp;

	res.a = Q.mul(a, right.a);

	if (ext){
		temp = Q.mul(b, right.b);
		res.a = Q.sub(res.a, temp);

		res.b = Q.mul(a, right.b);
		temp = Q.mul(b, right.a);
		res.b = Q.add(res.b, temp);
	}

	res.ext = ext;

	return res;
}
complex complex::pow(const apint &n){


	complex res(apint(1),apint(0),ext);

	for (int i = 1; i <= n.bits(); i++){
		res = res*res;

		if (n.checkbit(i)){
			res = res* (*this);
		}
	}

	return res;

}

complex complex::neg() const{
	modulo_arithmic Q = getQ();

	complex res;

	res.a = Q.neg(a);

	res.b = Q.neg(b);
	res.ext = ext;

	return res;

}
complex complex::inv() const{
	modulo_arithmic Q = getQ();

	complex res;

	apint temp, denom;

	if (!ext){
		res.ext = ext;
		res.a = Q.inv(a);

		return res;
	}

	denom = Q.mul(a, a);
	temp = Q.mul(b, b);
	denom = Q.add(denom, temp);

	res.a = Q.div(a, denom);

	temp=Q.neg(b);
	res.b = Q.div(temp, denom);
	res.ext = ext;

	return res;

}

complex complex::operator/(const complex &b){
	complex temp=b.inv();
	return (*this)* temp;
}

point point::operator+(const point &right) const{
	if (id) return right;
	if (right.id)return *this;

	point res;
	apint two(2);
	complex temp1,temp2;

	if ((*this) == right){
		res=two*right;
		return res;
	}

	if (x == right.x){
		res.id = 1;
		return res;
	}


	complex lambda =  (right.y-y)/(right.x-x);
	
	complex niu = y-lambda*x;

	res.x = lambda*lambda - x - right.x;

	res.y = (lambda*res.x + niu).neg();

	return res;

}
point point::inv() const{
	point res(*this);

	res.y = res.y.neg();
	return res;
}
point point::operator-(const point &right) const{
	point res(*this);

	return res+right.inv();
}

point operator*(const apint &s,const point &P){
	point res;
	complex three(apint(3), apint(0), P.x.ext);
	complex two(apint(2), apint(0), P.x.ext);
	complex paraA(apint(A), apint(0), P.x.ext);
	complex temp1,temp2;
	apint TWO(2);

	if (s.length == 1 && s.data[0] == 2){//double a point

		if (P.id)
			res.id = 1;
		else if ((two * P.y).zero())
			res.id = 1;
		else{

			temp1=three*P.x;
			temp1=temp1*P.x;
			temp1=temp1+paraA;
			temp2=two*P.y;

			complex lambda = temp1/ temp2;

			temp1=lambda*P.x;
			complex niu = P.y - temp1;


			temp1=lambda*lambda;
			temp1=temp1-P.x;
			res.x=temp1-P.x;

			res.y = (lambda*res.x + niu).neg();			
		}
	}
	else{
		res = P;

		for (int i = 2; i <= s.bits(); i++){
			res = TWO*res;
			//res.on_curve();

			if (s.checkbit(i)){
				res = res+P;
			//	res.on_curve();
			}
		}
	}
	return res;
}

point operator*(int s,const point &P){
	apint S(s);
	point res=S*P;
	return res;
}

modulo_arithmic getQ(){
	apint m(0);
	apint ten(10);
	int i;

	const char* qval="87807107996633125224377819847540498158068831994142082110286";//len 59
	const char* qval2="5339926647563088022295707862517942266222142315585876958231745";//len 61
	const char* qval3="9277713367317481324925129998224791";//len 34
	
	for(i=0;i<59;i++){
		m=m*ten;
		m=m+apint(qval[i]-'0');
	}
	for(i=0;i<61;i++){
		m=m*ten;
		m=m+apint(qval2[i]-'0');
	}
	for(i=0;i<34;i++){
		m=m*ten;
		m=m+apint(qval3[i]-'0');
	}
	
	if(!LIGHT)return modulo_arithmic(m);
	return modulo_arithmic(MODULO);

}
