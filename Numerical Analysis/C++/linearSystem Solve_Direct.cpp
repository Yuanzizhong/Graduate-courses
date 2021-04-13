#include<iostream>
using namespace std;
#include<cmath>
#define m 3   //矩阵的行数
#define n 3   //矩阵的列数


//高斯消去法--返回一个数组--故要用一个指针------double matrix[][n]==double (*matrix)[n]
double *gauss(double matrix[][n],double *b,double x[])//double x[] ==double *x  不改变值如下 void cit(const int b[],int n)
{	
	
	//遍历矩阵的每一列
	for (int j = 0; j < n-1;j++) {

		//消元
		for (int i = j + 1; i < m;i++) {

			double p = matrix[i][j] / matrix[j][j];//这一步要认真看看
			for (int k = 0; k < n; k++) {

				matrix[i][k] = matrix[i][k] - matrix[i-1][k] *p;
			}
			b[i] = b[i] - b[i-1]*p;//这个p也要替换
		
		}

	}
	//回代
	x[n - 1] = b[n - 1] / matrix[n - 1][n - 1];

	for (int j = n - 2; j >= 0; j--) {
		double sum = 0; 

		for (int i = j+1; i < n; i++) {
			sum = x[i] * matrix[j][i] +sum;
			
		}
		x[j] = (b[j] - sum) / matrix[j][j];
		
	}
	
	return x;
}

void PrintSolve(double x[])
{
	//打印解x
	cout << "解x为：" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << x[i] << endl;
	}
}

void PrintMatrix(double matrix[][n])
{
	//打印解matrix
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

//列主元消去法
double *ColGauss(double (*matrix)[n],double b[],double *x)
{

	//遍历矩阵的每一列
	for (int j = 0; j < n - 1; j++) {
		double max = matrix[j][j];
		int idnex = j;

		////记录每一列的最大值
		for (int i = j + 1; i < m; i++) 
		{	
			if (abs(matrix[i][j])>=max) {
				max = matrix[i][j];
				idnex = i;
			}
		}
		//交换第i-1行与，最大值的一行
		for (int k=0;k<n;k++)
		{
			double temp = matrix[j][k];
			matrix[j][k] = matrix[idnex][k];
			matrix[idnex][k] = temp;
		}

		//消元
		for (int i = j + 1; i < m; i++) 
		{

			double p = matrix[i][j] / matrix[j][j];//这一步要认真看看
			for (int k = 0; k < n; k++) {

				matrix[i][k] = matrix[i][k] - matrix[i - 1][k] * p;
			}
			b[i] = b[i] - b[i - 1] * p;//这个p也要替换

		}
	}
	//回代
	x[n - 1] = b[n - 1] / matrix[n - 1][n - 1];

	for (int j = n - 2; j >= 0; j--) {
		double sum = 0;

		for (int i = j + 1; i < n; i++) {
			sum = x[i] * matrix[j][i] + sum;

		}
		x[j] = (b[j] - sum) / matrix[j][j];

	}

	return x;
}

//LU分解
void LU(double matrix[][n])
{
	for (int k =0;k<m;k++)
	{
	
		for (int j = k; j < n; j++)
		{
			for (int i = 0; i <= k-1 ; i++) //这一步关键呀
			{
				matrix[k][j] = matrix[k][j] - matrix[k][i] * matrix[i][j];
			}
		}
		for (int i = k + 1; i < m; i++) {
			for (int j = 0; j <= k - 1;j++)
			{
				matrix[i][k] = matrix[i][k] - matrix[i][j] * matrix[j][k];
			}
			matrix[i][k] = matrix[i][k] /matrix[k][k];
		}
	
	}

}

//LU1分解--Dolitle
void LU1(double matrix[][n])
{
	for (int i =0;i<m;i++)
	{
		for (int j = i; j < n; j++)
		{
			double sum1 = 0;
			for (int k=0;k<i;k++)
			{
				sum1 = sum1 + matrix[i][k] * matrix[k][j];
			}
			matrix[i][j] = matrix[i][j] - sum1;
		}
		
		for (int j = i + 1; j < m; j++)
		{
			double sum2 = 0;
			for (int k=0;k<i;k++)
			{
				sum2=sum2+matrix[j][k]*matrix[k][i]; //中间变量sum2可以略去
			}
			matrix[j][i] = (matrix[j][i] - sum2)/matrix[i][i];
		}
	

	}
}

//LU
double* LU_Solve(double matrix[][n],double *b,double *x)
{
	x[0] = b[0];
	//求解Ly=b
	for (int j = 1; j < n; j++)
	{
		double sum = 0;
		for (int i = 0; i < j; i++)
		{
			sum = sum + matrix[j][i] * x[i];
		}
		x[j] = b[j] - sum;
	}
	//回代
	x[n - 1] = x[n - 1] / matrix[n - 1][n - 1];

	for (int j = n - 2; j >= 0; j--) {
		double sum = 0;

		for (int i = j + 1; i < n; i++) {
			sum = x[i] * matrix[j][i] + sum;

		}
		x[j] = (x[j] - sum) / matrix[j][j];

	}

	
	return x;
}
//Crout分解
void crout(double matrix[][n])
{
	for (int i = 0; i < m; i++)
	{

		for (int j = i ; j < m; j++)
		{
			double sum2 = 0;
			for (int k = 0; k < i; k++)
			{
				sum2 = sum2 + matrix[j][k] * matrix[k][i]; //中间变量sum2可以略去
			}
			matrix[j][i] = (matrix[j][i] - sum2);
		}

		for (int j = i+1; j < n; j++)
		{
			double sum1 = 0;
			for (int k = 0; k < i; k++)
			{
				sum1 = sum1 + matrix[i][k] * matrix[k][j];
			}
			matrix[i][j] = (matrix[i][j] - sum1)/matrix[i][i];
		}

	}
}

double* crout_Solve(double matrix[][n], double* b, double* x)
{
	x[0] = b[0]/matrix[0][0];
	//求解Ly=b
	for (int j = 1; j < n; j++)
	{
		double sum = 0;
		for (int i = 0; i < j; i++)
		{
			sum = sum + matrix[j][i] * x[i];
		}
		x[j] = (b[j] - sum)/matrix[j][j];
	}
	//回代
	x[n - 1] = x[n - 1];

	for (int j = n - 2; j >= 0; j--) {
		double sum = 0;

		for (int i = j + 1; i < n; i++) {
			sum = x[i] * matrix[j][i] + sum;

		}
		x[j] = (x[j] - sum) ;

	}


	return x;
}

void Cholesky(double matrix[][n])
{
	for (int j = 0; j < n; j++)
	{
		double sum1 = 0;
		for (int i=j;i<m;i++)
		{
			if (i == j)
			{
				//先算对角线元素的值
				for (int k = 0; k < j; k++) {
					sum1 = sum1 + matrix[i][k] * matrix[i][k];

				}

				matrix[j][j] = sqrt(matrix[j][j] - sum1);
			}
			else
			{
				//再算其他的元素
				double sum2 = 0;
				for (int k = 0; k <j ; k++)
				{
					sum2 = sum2 + matrix[i][k] * matrix[j][k];
				}
				matrix[i][j] = (matrix[i][j] - sum2) / matrix[j][j];
			}
		}
	}
}
double* Cholesky_Solve(double matrix[][n], double* b, double* x)
{
	x[0] = b[0] / matrix[0][0];
	//求解Ly=b
	for (int j = 1; j < n; j++)
	{
		double sum = 0;
		for (int i = 0; i < j; i++)
		{
			sum = sum + matrix[j][i] * x[i];
		}
		x[j] = (b[j] - sum) / matrix[j][j];
	}
	//回代
	x[n - 1] = x[n - 1]/matrix[m-1][n-1];

	for (int j = n - 2; j >= 0; j--) {
		double sum = 0;

		for (int i = j + 1; i < n; i++) {
			sum = x[i] * matrix[i][j] + sum; //这里改一下i,j的位置就好（与上面不同）

		}
		x[j] = (x[j] - sum)/matrix[j][j];

	}


	return x;
}
void P_Cholesky(double matrix[][n],double b[],double x[])
{
	for (int j=0;j<n;j++)
	{
		for (int k = 0; k < j;k++) {
			matrix[j][j] = matrix[j][j] -matrix[j][k] * matrix[j][k] * matrix[k][k];
		}
		for (int i = j+1; i < n; i++)
		{
			for (int k = 0; k < j;k++) {//这里不要写成j++
				matrix[i][j] = matrix[i][j]-matrix[i][k] * matrix[k][k] * matrix[j][k];
			}
			matrix[i][j] = matrix[i][j] / matrix[j][j];
		}
	}
	x[0] = b[0];
	//求解Ly=b
	for (int j = 1; j < n; j++)
	{
		double sum = 0;
		for (int i = 0; i < j; i++)
		{
			sum = sum + matrix[j][i] * x[i];
		}
		x[j] = (b[j] - sum);
	}
	x[n - 1] = x[n - 1] / matrix[n - 1][n - 1];//这里不是b/matrix


	for (int i=n-2;i>=0;i--)
	{
		x[i] = x[i] / matrix[i][i];
		for (int j=i+1;j<n;j++)
		{
			x[i] = x[i] - matrix[j][i] * x[j];
		}
	}

}
void RunSolveMatrix(double matrix[][n],double b[],double x[])
{
	matrix[0][1] = matrix[0][1] / matrix[0][0];
	x[0] = b[0] / matrix[0][0];
	for (int i = 1; i < n - 1; i++)
	{
			matrix[i][i] = matrix[i][i] - matrix[i][i-1] * matrix[i-1][i];
			matrix[i][i +1] = matrix[i][i + 1] / matrix[i][i];
			x[i] = (b[i] - matrix[i][i-1] * x[i - 1]) / matrix[i][i];
	}
	matrix[m - 1][n - 1] = matrix[m-1][n-1] - matrix[m - 1][n - 2] * matrix[n - 2][m - 1];
	x[n - 1] = (b[n - 1] - matrix[m - 1][n - 2] * x[n - 2]) / matrix[n - 1][n - 1];
	for (int i=n-2;i>=0;i--)
	{
		x[i] = x[i] - matrix[i][i+1] * x[i + 1];
	}

}


int main() {
	//定义矩阵，向量x，向量b
	double matrix[m][n] = {
		{2,1,0},
		{1,4,2},
		{0,2,6}
	};
	double x[n] = {0,0,0};
	double b[n] = {3,7,8};


	//调用高斯消去法函数,接受向量x
	//gauss(matrix, b, x);
	//PrintSolve(x);

	//调用列主元消去法函数
	//ColGauss(matrix, b, x);
	//PrintSolve(x);
	
	
	//调用LU函数
	//LU(matrix);
	//PrintMatrix(matrix);
	//PrintSolve(x);
	//LU1(matrix);
	//PrintMatrix(matrix);
	//LU_Solve(matrix, b, x);
	//PrintSolve(x);

	//调用court 函数
	//crout(matrix);
	//PrintMatrix(matrix);
	//crout_Solve(matrix, b, x);
	//PrintSolve(x);
	// 
	// 
	//doolittle函数
	//Cholesky(matrix);
	//PrintMatrix(matrix);
	//Cholesky_Solve(matrix, b, x);
	//PrintSolve(x);

	//改进的Cholesky
	//P_Cholesky(matrix, b, x);
	//PrintMatrix(matrix);
	//PrintSolve(x);

	//追赶法
	RunSolveMatrix(matrix, b, x);
	PrintMatrix(matrix);
	PrintSolve(x);

	system("pause");
	return 0;
}
