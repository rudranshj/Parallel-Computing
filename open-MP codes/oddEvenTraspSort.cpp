#include<iostream>
#include<omp.h>
using namespace std;


// similar to bubbleSort
// different pass, operations depend whether pass is even or odd
// after n passes the array will be sorted
// loop carried dependency for outer loop
// however inner loop doesnt have dependency because i+=2 :)

void oddEvenTSortSerial(int a[], int n){
	int temp=0;
	for(int pass=0; pass<n; pass++){
		if(pass%2 > 0){ 
		// odd pass
			for(int i=1; i<n-1; i+=2){
				if(a[i]>a[i+1]){
					temp=a[i];
					a[i]=a[i+1];
					a[i+1]=temp;
				}
			}
		}

		else{
			for(int i=1; i<n ; i+=2){
				if(a[i-1]>a[i]){
					temp=a[i-1];
					a[i-1]=a[i];
					a[i]=temp;
				}
			}
		}
	}


	for(int i=0; i<n; i++){
		cout<<a[i]<<" ";
	}
	cout<<endl;
}

void oddEvenTSortParallel(int a[], int n, int thread_count){
	// cout<<"no of threads:"<<thread_count<<endl;
	int i=0,pass=0,temp=0;
	# pragma omp parallel num_threads(thread_count) default(none) private(i,temp,pass) shared(a,n)
	for(int pass=0; pass<n; pass++){
		if(pass%2 > 0){ 
		// odd pass
			# pragma omp for
			for(int i=1; i<n-1; i+=2){
				if(a[i]>a[i+1]){
					temp=a[i];
					a[i]=a[i+1];
					a[i+1]=temp;
				}
			}
		}

		else{
			# pragma omp for
			for(int i=1; i<n ; i+=2){
				if(a[i-1]>a[i]){
					temp=a[i-1];
					a[i-1]=a[i];
					a[i]=temp;
				}
			}
		}
	}


	for(int i=0; i<n; i++){
		cout<<a[i]<<" ";
	}
	cout<<endl;
}


// below version forks and combine threads multiple times and therefore increases the overhead
void oddEvenTSortParallel2(int a[], int n, int thread_count){
	// cout<<"no of threads:"<<thread_count<<endl;
	int i=0,pass=0,temp=0;
	for(int pass=0; pass<n; pass++){
		if(pass%2 > 0){ 
		// odd pass
			# pragma omp parallel for num_threads(thread_count) default(none) private(i,temp) shared(a,n)
			for(int i=1; i<n-1; i+=2){
				if(a[i]>a[i+1]){
					temp=a[i];
					a[i]=a[i+1];
					a[i+1]=temp;
				}
			}
		}

		else{
			# pragma omp parallel for num_threads(thread_count) default(none) private(i,temp) shared(a,n)
			for(int i=1; i<n ; i+=2){
				if(a[i-1]>a[i]){
					temp=a[i-1];
					a[i-1]=a[i];
					a[i]=temp;
				}
			}
		}
	}


	for(int i=0; i<n; i++){
		cout<<a[i]<<" ";
	}
	cout<<endl;
}



int main(int argc, char* argv[]){
	int a[10]={1,4,2,4,9,3,8,7,-1,5};
	int thread_count=1;

	if(argc==2){
		thread_count=strtol(argv[1], NULL, 10);
	}
	else{
		cout<<"No cmd line arg given. huh. quittind"<<endl;
		return 1;
	}
	cout<<"serial version \n";
	oddEvenTSortSerial(a,10);
	cout<<"parallel version (initial) \n";
	oddEvenTSortParallel2(a,10,thread_count);
	cout<<"parallel version (final) \n";
	oddEvenTSortParallel(a,10,thread_count);
	return 1;
}


