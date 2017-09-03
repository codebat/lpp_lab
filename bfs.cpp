#include<iostream>
#include<vector>
#include<math.h>
using namespace std;
typedef vector<double> vd;
typedef vector<vd> vvd;
typedef vector<int> vi;


//global variable declaration
int n,m,diff,nbfs=0,nbs=0;
vvd a(100,vd(100,0)),b(100,vd(100,0));
vi c(100,0),p(100,0);

void exch(int p,int q)
{
    for(int i=0;i<=m;i++)
    {
        int t=b[q][i];
        b[q][i]=b[p][i];
        b[p][i]=t;
    }    
} 


void solve_equation(){
	
	//create a m*m  matrix from m*n matrix using the position of non basis variable
	c[diff] = n;
	for(int i=0;i<m;i++){
		int n_ind = 0,o_ind = 0;
		for(int j = 0;j<=diff;j++){
			for(int k=o_ind;k<c[j];k++){
				b[i][n_ind] = a[i][k];
				//cout<<i<<n_ind<<k<<b[i][n_ind]<<a[i][k]<<"\n";
				if(i==0)p[n_ind]=k;n_ind++;
			}
			o_ind = c[j]+1;
		}
		b[i][m] = a[i][n];
	}
	
/* 	for(int i=0;i<m;i++){
		for (int j=0;j<=m;j++){
			cout<<b[i][j]<<" ";
		}
		cout<<"\n";
	} */
	//reduce to row-echelon form
	for(int i=1;i<m;i++){
		for(int j=0;j<i;j++){
			if(b[i][j]==0)continue;
			else if(b[j][j]==0)exch(i,j);
			else{
				double f = b[i][j]/b[j][j];
				for(int k=j;k<m;k++){
					b[i][k]-=(b[j][k]*f);
				}
			}
		}
	}
	vd sol(10,0);
	for(int i=0;i<diff;i++)sol[c[i]]=0;
	sol[p[m-1]] = b[m-1][m]/b[m-1][m-1];
	for(int i=2;i<m+1;i++){
		double ans = b[m-i][m];
		for(int j=1;j<i;j++){
			ans-=sol[p[m-i]]*b[m-i][m-j];
		}
		sol[p[m-i]] = ans/b[m-i][m-i];
	}
	bool bs = true,bfs = true;
	for(int i=0;i<n;i++){
		if(sol[i]!=sol[i]|| isinf(sol[i]))bs=false;
		if(sol[i]!=sol[i]|| isinf(sol[i])||sol[i]<0)bfs = false;
			}
	if(bs){
		nbs++;
		for(int i=0;i<n;i++){
			cout<<"x"<<i+1<<": "<<sol[i]<<" ";
		}
		cout<<"\n";
	}
	if(bfs){
		nbfs++;
		cout<<"This is a basic feasible solution\n";
	}
}


void solve_all_basic(int pos,int var){
	if(pos==diff){
		solve_equation();
		return;
	}
	for(int i=var+1;i<=m+pos;i++){
		c[pos] = i;
		solve_all_basic(pos+1,i);
	}
	return;
}

int main(){

	cout<<"Enter no of equation and no of variables respectively\n";
	cin>>m>>n;
	
	diff = n-m;
	
	cout<<"enter the coefficient of the matrix\n";
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			cout<<"x"<<j<<": ";
			cin>>a[i][j];
		}
		cout<<"rhs: ";
		cin>>a[i][n];
	}
/* 	for(int i=0;i<m;i++){
		for (int j=0;j<=n;j++){
			cout<<a[i][j]<<" ";
		}
		cout<<"\n";
	} */
	for(int i=0;i<=m;i++){
		c[0] = i;
		solve_all_basic(1,i);
	}
	cout<<"Total no of basic solution is "<<nbs<<"\n";
	cout<<"Total no of basic feasible solution is "<<nbfs<<"\n";
	return 0;
}