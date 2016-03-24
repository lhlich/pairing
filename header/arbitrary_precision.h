#ifndef ARBITRARY_PRECISION_H
#define ARBITRARY_PRECISION_H 
#define high(x) x>>16
#define low(x)  (x<<16)>>16
#define to_high(x)   (((1<<16)-1)&x)<<16
#define to_low(x) ((1<<16)-1)&(x>>16)



struct apint{// int class of arbitrary precision.

public:
	apint();
	apint(unsigned int x);
	apint(const apint &obj);
	apint(unsigned int* buf, int len);
	apint& operator=(const apint& obj) {
		delete[] data;

		data = new unsigned int[obj.length];
		length = obj.length;
		for (int i = 0; i<length; i++){
			data[i] = obj.data[i];
		}
		return *this;
	}

	~apint(){
		delete[] data;
	}

	bool checkbit(int &pos) const;
//	apint square();
	apint operator+(const apint& right) const;
	apint operator*(const apint& right) const;
	apint div(const apint& right, apint &rem);
	apint operator-(const apint& right);//design for the case a>=b
	bool operator>=(const apint& right) const;
	bool operator==(const apint& right) const;
	bool operator!=(const apint& right) const;

	apint operator<<(const int& n) const;
	int bits() const;

	void remove_prefix();

	bool zero() const{
		return (length==0 || length==1 && data[0]==0);
	}

	unsigned int* data;   // "big endian style": LSB at data[length-1]
	int length;   // this integer contains length bytes
};

class modulo_arithmic{// class holding modulo operation

public:
	modulo_arithmic(unsigned int x);
	modulo_arithmic(apint &x);
	modulo_arithmic(){
		p = apint(2333);
	}

	apint add(const apint &a, const apint &b); //(a+b)mod p
	apint sub(const apint &a, const apint &b); //(a-b) mod p
	apint mul(const apint &a, const apint &b);  //(a*b) mod p
	apint pow(const apint &x, const apint &n);   //(x^n) mod p
	apint inv(const apint &a);           //a^(-1) in Zp
	apint div(const apint &a, const apint &b);  //(a/b) mod p
	apint neg(const apint &x);   //-a in Zp

	void reduce(apint &x);

private:
	apint p;
	//... and some stuff for precomputation.
};

#endif
