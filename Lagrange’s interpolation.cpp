#include<iostream>
using namespace std;
int main(){
    int n,i,j;
    float ax[100],gx[100],nu,de,x,y=0.0;
    cout<<"Enter the number of data: ";
    cin>>n;
    cout<<"Enter the elements in x: ";
    for (int i=1;i<=n;i++){
        cin>>ax[i];
    }
    cout<<"ENter the elemets in fx: ";
    for(int i=1;i<=n;i++){
        cin>>gx[i];
    }
    cout<<"Enter the value of x: ";
    cin>>x;
    for (int i=1;i<=n;i++){
        nu=1.0;
        de=1.0;
        for(int j=1;j<=n;j++){
            if(i!=j){
                nu*=x-ax[j];
                de*=(ax[i]-ax[j]);

        }
    }
    y=y+(nu/de)*gx[i];
}
cout<<y;
}
