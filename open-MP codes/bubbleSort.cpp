#include<iostream>
// #include<bits/stdc++.h>
using namespace std;

void bsSerial(int a[], int n){
	int tmp;

	for(int l=n; l>=2; l--){
		for(int i=0; i<l-1; i++){
			if (a[i]>a[i+1]){
				tmp=a[i];
				a[i]=a[i+1];
				a[i+1]=tmp;
			}
		}
	}
	for(int i=0; i<n; i++){
		cout<<a[i]<<" ";
	}
	cout<<endl;	

}

// so theres loop carried dependency here, not easy to parallelise

int main(){
	int a[10]={1,4,2,4,9,3,8,7,-1,5};
	bsSerial(a,10);
	return 1;
}