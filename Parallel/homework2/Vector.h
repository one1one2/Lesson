#ifndef __Vector_h__
#define __Vector_h__
#include<iostream>
#include<cstdlib>
template<class T>
class Vector{
private:
	T *data;
	int length;
public:
	Vector<T>():data(NULL),length(0){}
	Vector<T>(int length1):data(new T[length1]),length(length1){}
	Vector<T>(const Vector<T> &b){ *this = b; }
	~Vector<T>(){
		delete data;
	}
	int size() const{
		return length;
	}
	T& operator[](int i){
		if (i >= 0 && i < length) 
			return data[i];
		else{
			std::cerr<<"数组越界";
			std::abort();
		}
	}
	const T& operator[](int i) const{
		if (i >= 0 && i < length) 
			return data[i];
		else{
			std::cerr<<"数组越界";
			std::abort();
		}
	}
	void resize(int n){ //单纯改变长度，并且直接将所有元素赋值为0.
		data = new T[n];
		length = n;
	}
	Vector<T>& operator+=(const Vector<T> &b);
	Vector<T>& operator*=(T b);
	Vector<T>& operator=(const Vector<T> &b);
};

template<class T>
Vector<T>& Vector<T>::operator+=(const Vector<T> &b){
	if (length != b.size()){
		std::cerr<<"Length is not the same!";
		std::abort();
	} 
	for (int i = 0; i < length; i++){
		data[i] += b[i]; 
	}
	return *this;
}
template<class T>
Vector<T>& Vector<T>::operator*=(T b){
	for (int i = 0; i < length; i++){
		data[i] *= b; 
	}
	return *this;
}
template<class T>
Vector<T>& Vector<T>::operator=(const Vector<T> &b){
	length = b.size();
	data = new T[length];
	for (int i = 0; i < length; i++){
		data[i] = b[i]; 
	}
	return *this;	
}
#endif
