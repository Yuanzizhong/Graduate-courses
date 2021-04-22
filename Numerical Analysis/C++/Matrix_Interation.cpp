#include<iostream>
#include<cmath>
using namespace std;
#define m 3
#define n 3
#define w 1.2

//打印x
void PrintSolve(double x[])
{
	//打印解x
	cout << "解x为：" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << x[i] << endl;
	}
}
int Jcobi(double matrix[][n], double b[], double* x)
{
	int k = 0;
	double y[n];
	while (1)
	{
		double error = 0;   //每次轮到这里都可以变为0
		for (int i = 0; i < m; i++)
		{
			double sum = 0; //每次轮到这里都可以变为0
			for (int j = 0; j < n; j++) {
				if (i != j)
				{
					sum = sum + matrix[i][j] * x[j];
				}
			}
			y[i] = (b[i] - sum) / matrix[i][i];
			error = error + abs(y[i] - x[i]);  //一范数
		}
		k = k + 1; //迭代步数
		if (error < 1e-4)
		{
			break;
		}
		else
		{
			for (int i = 0; i < n; i++)
			{
				x[i] = y[i];//更新向量x
			}
		}
	}
	return k;
}
int GaussSeidel(double matrix[][n],double b[],double *x)
{
	int k = 0;
	while (1) 
	{	
		double error = 0;   //每次轮到这里都可以变为0
		for (int i = 0; i < m; i++)
		{	
			double sum =0 ; //每次轮到这里都可以变为0
			for (int j = 0; j < n; j++) {
				if (i != j)
				{
					sum = sum + matrix[i][j] * x[j]; 
				}
			}
			double temp = x[i];
			x[i] = (b[i] - sum) / matrix[i][i];
			error = error + abs(temp - x[i]);  //一范数

			//不能这样算，因为前面的x[i]已经更新了--sum不准确
			//x[i] = (b[i] - sum - matrix[i][i] * x[i]) / matrix[i][i];
			
		}
		k = k + 1; //迭代步数
		if (error < 1e-4)
		{
			break;
		}
	}
	return k;
}

int SOR(double matrix[][n], double b[], double* x)
{
	int k = 0;
	while (1)
	{
		double error = 0;   //每次轮到这里都可以变为0
		for (int i = 0; i < m; i++)
		{
			double sum = 0; //每次轮到这里都可以变为0
			for (int j = 0; j < n; j++) {

					sum = sum + matrix[i][j] * x[j];
			}
			double temp = x[i];
			x[i] = x[i]+(b[i] - sum)*w / matrix[i][i];
			error = error + abs(temp - x[i]);  //一范数

		}
		k = k + 1; //迭代步数
		if (error < 1e-4)
		{
			break;
		}
	}
	return k;
}
//最速下降法
int Speed(double matrix[][n], double b[], double* x)
{
	int k = 0;
	double error;
	do
	{
		//计算残差r[i]
		double r[n];
		double sum1 =0;
		for (int i = 0; i < m; i++)
		{
			double sum = 0;
			for (int j = 0; j < n; j++)
			{
				sum =sum+ matrix[i][j] * x[j];
			}
			//计算a的分子
			r[i] = b[i] - sum;
			sum1 = sum1 + (r[i])* (r[i]);
		}

		//计算a的分母
		double sum2 = 0;
		for (int i = 0; i < m; i++)
		{
			double q = 0;
			for (int j = 0; j < n; j++)
			{
				q = q+ matrix[i][j] * r[j]; //好像可以合并
			}
			sum2 = sum2 + q * r[i];
		}

		//更新x和误差
		error = 0;
		for (int i = 0;i<n;i++)
		{
			error = error +abs(r[i]); //误差
			x[i] = x[i] + r[i] * (sum1 / sum2);
		}
		k = k + 1;
	} while (error>1e-4);//注意迭代条件--是大于

	return k;
}


//共轭梯度法
int ConGradient(double matrix[][n], double b[], double* x)
{

	//计算残差r[0]、p[0]
	double r[n];
	double p[n];
	
	for (int i = 0; i < m; i++)
	{
		double sum = 0;
		for (int j = 0; j < n; j++)
		{
			sum = sum + matrix[i][j] * x[j];
		}
		r[i] = b[i] - sum;
		p[i] = r[i];
		
	}
	
	int k = 0;
	double error;
	do
	{
		//计算alpha
		double sum1 = 0;
		double sum2 = 0;
		for (int i = 0; i < m; i++)
		{
			double q = 0;
			for (int j = 0; j < n; j++)
			{
				q = q + matrix[i][j] * p[j]; //好像可以合并
			}
			sum2 = sum2 + q * p[i];

			sum1 = sum1 + r[i] * r[i];
		}

		//更新x和误差
		error = 0;
		for (int i = 0; i < n; i++)
		{
	
			error = error + abs(p[i]); //误差
			//cout << error << endl;
			x[i] = x[i] + p[i] * sum1 / sum2;
			
		}
		
		//更新r[k+1]
		double sum3 = 0;
		for (int i = 0; i < m; i++)
		{
			double sum4=0;
			for (int j = 0; j < n; j++)
			{
				sum4 =  sum4 +matrix[i][j] * p[j];//是乘以p[i],不是乘以x[i]
			}
		
			//计算a的分子
			sum3 = sum3 + r[i] * r[i];
			//cout << sum4 << endl;
			r[i] = r[i] - sum4 * sum1/ sum2;
			
		}
		
		double beta =0;
		for (int i = 0; i < n; i++)
		{
			beta = beta + r[i] * r[i];
		}
		
		beta = beta / sum3;
	
		for (int i = 0; i < n; i++)
		{
			p[i] = r[i] + beta * p[i];
		}

		k = k + 1;

	} while (error > 1e-4);//注意迭代条件--是大于

	return k;
}
int main()
{
	double matrix[m][n] = {
		{4,3,0},//打错了一个数，写成4.3
		{3,4,-1},
		{0,-1,4}
	};
	double b[] = {24,30,-24};
	double x[] = {1,1,1};


	//调用GaussSeidel
	//cout << "gaussseidel "<<endl<<GaussSeidel(matrix, b, x)<< endl;
	//PrintSolve(x);

	//调用Jacobo
	//cout <<"jacobi  "<< Jcobi(matrix, b, x) << endl;
	//PrintSolve(x);

	//调用SOR
	cout <<"sor  " <<SOR(matrix, b, x) << endl;
	PrintSolve(x);

	//调用最速下降法
	//cout<<Speed(matrix, b, x)<<endl;
	//PrintSolve(x);
	//system("pause");

	//共轭梯度法
	//cout<<ConGradient(matrix, b, x)<<endl;
	//PrintSolve(x);

	return 0;
}
