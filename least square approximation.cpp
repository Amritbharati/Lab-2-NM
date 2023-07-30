#include<iostream>
#include<cmath>
using namespace std;
	int main(){
		double x[5],y[5],sum_x,sum_y,sum_x2,sum_xy,n=5;
		cout<<"Enter elements for x : ";
			for(int i=0;i<5;i++){
				cin>>x[i];
			}	
		cout<<"Enter elements for y : ";
			for(int i=0;i<5;i++){
				cin>>y[i];
			}
		for(int i=0;i<5;i++){
			sum_x+=x[i];
			sum_y+=y[i];
			sum_x2+=pow(x[i],2);
			sum_xy+=(x[i]+y[i]);
		}
		double A,B;
		B=((n*sum_xy)-(sum_x*sum_y))/((n*sum_x2)-pow(sum_x,2));
		A=((sum_y)-(B*sum_x))/n;
		cout<<"B = "<<B<<endl;
		cout<<"A = "<<A<<endl;
		return 0;
	}
