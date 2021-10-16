#include<iostream>
#include<cmath>

using namespace std;
double eigenValues[100];
double eigenVectors[100][100];
double im[100][100];
int num=0;
int nhi=0;
void matrixInput(double a[][100], int rows)
{
	for(int i=1;i<=rows;i++)
	{
		for(int j=1;j<=rows;j++)
		{
			cout << "a[" <<i<<"]["<<j<<"] = "; cin >> a[i][j];
			 cout << endl;
		}
	}
}
void matrixOutput(double a[][100], int rows)
{
	for(int i=1;i<=rows;i++)
	{
		for(int j=1;j<=rows;j++)
		{
			if(abs(a[i][j])<0.000001) cout << 0 <<"  " ;
			else cout << a[i][j] << "  " ;
		}
		cout << endl;
	}
}
void matrixOutput2(double a[][100], int rows,int cols)
{
	for(int i=1;i<=rows;i++)
	{
		for(int j=1;j<=cols;j++)
		{
			if(abs(a[i][j])<0.000001) cout << 0 <<" " ;
			else cout << a[i][j] << " " ;
		}
		cout << endl;
	}
}
void matrixMultiplication(double X[][100], int a, int b, double Y[][100], int c, int d, double RES[][100])
{
	if(b!=c) return;
	for(int i=1;i<=a;i++)
	{
		for(int j=1;j<=d;j++)
		{
			RES[i][j]=0;
			for(int k=1;k<=b;k++)
			{
				RES[i][j]+=X[i][k]*Y[k][j];
			}
		}
	}
}
bool matrixEmptyCheck(double x[][100], int n1)
{
	for(int i=1;i<=n1;i++)
	{
		for(int j=1;j<=n1;j++)
		{
			if (x[i][j]!=0) return false;
		}
	}
	return true;
}
void nhiIsCloningAMatrix(double im[][100], double y[][100], int n1)
{
	for(int i=1;i<=n1;i++)
	{
		for(int j=1;j<=n1;j++)
		{
			im[i][j]=y[i][j];
		}
	}
}
void squareMatrixMultiplication(double A[100][100], double B[100][100], double res[100][100], int n)
{
	for (int i=1;i<=n;i++)
	{
		for(int j=1;j<=n;j++)
		{
			res[i][j] = 0;
			for(int t=1;t<=n;t++)
			{
				res[i][j]+=A[i][t]*B[t][j];
			}
		}
	}
}
void swpCols(double M[][100], int step, int swpidx, int n1)
{
	for (int i=1;i<=n1;i++)
	{
		int tem = M[i][step];
		M[i][step] = M[i][swpidx];
		M[i][swpidx]=tem;
	}
}
void deleteRowCol(double A[][100], int row, int col, int n)
{
	double tem[100][100],tem2[100][100];
	//delete row
	for (int i=1;i<=(n-1);i++)
	{
		for(int j=1;j<=n;j++)
		{
			if(i>=row)A[i][j]= A[i+1][j];
		}
	}
	//delete col
	for(int i=1;i<=n-1;i++)
	{
		for(int j=1;j<=n-1;j++)
		{
			if(j>=col) A[i][j]=A[i][j+1];
		}
	}
	
}
//////////////////////////////////////////////////////////////
int Danielewski(double A[][100], int n)
{
	double M[100][100], iM[100][100], tem[100][100],sM[100][100];
	int n1=n;
	for(int step=n-1;step>=1;step--)
	{
		int check =0; int swpidx=0;
		//================================
		if(A[step+1][step]==0)
				{
					for (swpidx=1;swpidx<step;swpidx++)
					{
						if(A[step+1][swpidx]!=0)
						{
							check = 1;
							break;
						}
					}
					if (check==0)
					{
						++num;
						++nhi;
						eigenValues[num]=A[step+1][step+1];
						deleteRowCol(A,(step+1),(step+1),n1);
						//************************
						//************************
						deleteRowCol(sM,(step+1),(step+1),n1);
						n1--;
						continue;
					}
					
				}
		//================================
		for(int i=1;i<=n1;i++)
		{
			for(int j=1;j<=n1;j++)
			{
				if (check == 1)
				{
					if(i==j) 
					{
						M[i][j]=1;
						iM[i][j]=1;
					}
					else
					{
						M[i][j]=0;
						iM[i][j]=0;
					}
				}
				else if(A[step+1][step]!=0)
				{
					if(i==step)
				{
					M[i][j]=A[i+1][j];
					if(j==step) iM[i][j]=1.0/A[step+1][step];
					else iM[i][j]=(double)-A[i+1][j]/A[step+1][step];
				}
				else if (i!=step)
				{
					if(i==j)
					{
						M[i][j]=1;
						iM[i][j]=1;
					}
					else
					{
						M[i][j]=0;
						iM[i][j]=0;
					}
				}
				}
			}
		}
		if(check==1)
		{
			swpCols(M,step,swpidx,n1);
			swpCols(iM,step,swpidx,n1);
			step++;
		}
		if (matrixEmptyCheck(sM,n1)==1) 
		{
			nhiIsCloningAMatrix(sM,iM,n1);
			
		}
		else
		{
			squareMatrixMultiplication(sM,iM,tem,n1);
			nhiIsCloningAMatrix(sM,tem,n1);
		}
		squareMatrixMultiplication(A,iM,tem,n1);
		squareMatrixMultiplication(M,tem,A,n1);
	}
	//matrixOutput(sM,n1);
	nhiIsCloningAMatrix(im,sM,n1);
	return n1;
}
void polynomialInitialization(double polynomial[][100], double A[][100], int n1)
{
	int power = n1;
	for(int i=0;i<=n1;i++)
	{
		polynomial[1][i]=power--; //power
		if(i==0)
		polynomial[0][i]=pow(-1,n1); //scalar
		else polynomial[0][i]=pow(-1,n1)*(-A[1][i]);
		//cout << polynomial[0][i] <<"  "<< polynomial[1][i] << endl;
		
	}
}
long double equationValue(double polynomial[][100], int n1,double x)
{
	long double s=0;
	for(int i=0;i<=n1;i++)
	{
		s+=polynomial[0][i]*pow(x,polynomial[1][i]);
		//cout << polynomial[0][i] <<" " <<polynomial[1][i] <<"---" <<s<<endl;
	}
	return s;
}
double nhiIsAccurisingTheRoots(double polynomial[][100], int n1, double c1, double c2)
{
	double mid;
	while(abs(equationValue(polynomial,n1,c1))>0.0000001)
	{
		mid=(c1+c2)/2.0;
		if(equationValue(polynomial,n1,c1)*equationValue(polynomial,n1,mid)<0) 
		{
			c2=mid;
		}
		else if(equationValue(polynomial,n1,c2)*equationValue(polynomial,n1,mid)<0) 
		{
			c1=mid;
		}
	}
	++num;
	eigenValues[num]=c1;
	return c1;
}
void nhiIsFindingTheRoots(double polynomial[][100], int n, int n1)
{
	if(num>=n) return;
	else
	{
		double m1=abs(polynomial[0][1]),m2=abs(polynomial[0][0]);
		for(int i=1;i<=n1;i++)
		{
			if(m1<abs(polynomial[0][i]))
			m1=abs(polynomial[0][i]);
		}
		for(int i=0;i<=n1-1;i++)
		{
			if(m1<abs(polynomial[0][i]))
			m1=abs(polynomial[0][i]);
		}
		double hi,lo;
		lo=(double)abs(polynomial[0][n1])/(abs(polynomial[0][n1])+m2);
		hi=(double)1+(m1/abs(polynomial[0][0]));
		//---------------------------------------
		double step=(hi-lo)/100.0;
		while (num<n1)
		{
			for(double x=-hi;x<=-lo;x+=step)
		{
			if (equationValue(polynomial,n1,x)==0)
			{
				++num;
				eigenValues[num]=x;
			}
			double c1, c2; 
			c1=x;
			c2=x+step;
			if(c2<=-lo&&equationValue(polynomial,n1,c1)*equationValue(polynomial,n1,c2)<0)
			{
				//cout << "day la" <<c1<<" " <<c2<<endl;
				double hihi = nhiIsAccurisingTheRoots(polynomial,n1,c1,c2);
				//cout << hihi <<endl;
			}
		}
		for(double x=lo;x<=hi;x+=step)
		{
			if (equationValue(polynomial,n1,x)==0)
			{
				++num;
				eigenValues[num]=x;
			}
			double c1, c2; 
			c1=x;
			c2=x+step;
			if(c2<=hi&&equationValue(polynomial,n1,c1)*equationValue(polynomial,n1,c2)<0)
			{
				//cout << "day la" <<c1<<" " <<c2<<endl;
				double hihi = nhiIsAccurisingTheRoots(polynomial,n1,c1,c2);
				//cout << endl<< "day la fx: " <<equationValue(polynomial,n1,hihi) << " cua " << hihi <<endl;
			
			}
		}
		}
		
	}
	
}
void eigenValuesOutput()
{
	cout << endl<< "Eigenvalues of the given matrix:" << endl;
	for(int i=1;i<=num;i++)
	{
		cout << eigenValues[i];
		if (i<num) cout << ", ";
	}
	cout << endl;
}
/*void matrixMultiplication(double X[][100], int a, int b, double Y[][100], int c, int d, double RES[][100])
{
	if(b!=c) return;
	for(int i=1;i<=a;i++)
	{
		for(int j=1;j<=d;j++)
		{
			for(int k=1;k<=b;k++)
			{
				RES[i][j]+=X[i][k]*Y[k][j];
			}
		}
	}
}*/
void nhiIsFindingEigenVects(int n1)
{
	cout << "EIGENDECOMPOSITION OF THE GIVEN MATRIX:" << endl;
	for(int ev=nhi+1;ev<=num;ev++)
	{
		double y[100][100]; double res[100][100];
		int t = n1-1;
		int i=0;
		while(t>=0) 
		{
			y[++i][1]=pow(eigenValues[ev],t);
			t--;
		}
		
		matrixMultiplication(im,n1,n1,y,n1,1,res);
		//matrixOutput2(res,n1,1); cout << endl;
		//cout << "[";
		for (int i=1;i<=n1;i++)
		{
			//if(abs(res[i][1])<0.00001) cout << 0<< "  ";
			//else cout << res[i][1] << "  ";
			eigenVectors[i][ev]=res[i][1];
		}
		//cout << "]";
		
		cout << endl;
		
	}
	//matrixOutput(eigenVectors,n1);
}
void eigenDecomposition(int n1)
{
	double eigenMatrix[100][100];
	for(int i=1;i<=num-nhi;i++)
	{
		for(int j=1;j<=num-nhi;j++)
		{
			if(i!=j) eigenMatrix[i][j]=0;
			else if(i==j) 
			{
				eigenMatrix[i][j]=eigenValues[nhi+i];
			}
		}
	}
	for(int i=1;i<=nhi;i++)
	cout << eigenValues[i] << " * "<<endl;
	matrixOutput(eigenMatrix,n1);
	cout << endl << " * " << endl;
	matrixOutput(eigenVectors,n1);
	cout << endl <<" * " <<endl;
	cout << "inverse of:" << endl;
	matrixOutput(eigenMatrix,n1);
	cout << "(is the reciprocal of its diagonal elements if all of them are non-zero, the zero ones stay the same)" << endl;
}
main()
{
	int count=2;
	while (count>=1)
	{
	double A[100][100]; int n;
	if (count==2)
	{
		n=3;
		A[1][1]=2;A[1][2]=1;A[1][3]=0;
		A[2][1]=1;A[2][2]=3;A[2][3]=1;
		A[3][1]=0;A[3][2]=1;A[3][3]=2;
		matrixOutput(A,n);
		cout <<endl; 
		string str;
		getline(cin, str);
		
	}
	
	else
	{
		cout << "n = "; cin >> n; cout << endl;
		matrixInput(A,n);
	}
	
	int n1=Danielewski(A,n);
	cout << "Danielewski matrix:"<<endl;
	matrixOutput(A,n1);
	double polynomial[2][100];
	polynomialInitialization(polynomial,A,n1);
	nhiIsFindingTheRoots(polynomial,n,n1);
	eigenValuesOutput();
	cout <<endl<< "inverse of m:"<<endl;
	matrixOutput(im,n1);
	cout << endl;
	nhiIsFindingEigenVects(n1);
	eigenDecomposition(n1);
	count--;
	}
	system("pause");
	
	
}
