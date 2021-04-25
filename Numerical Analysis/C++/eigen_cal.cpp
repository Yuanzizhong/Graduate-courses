#include<iostream>
using namespace std;
#define m 3 //定义矩阵的行数和列数
#define n 3

//打印最大特征值的方向
double Print_x(double y[])
{
	double max = abs(y[0]);
	for (int i = 1; i < n; i++) {
		if (abs(y[i]) >= max) {
			max = y[i];
		}
	}
	cout << "最大特征值为：" << max << endl;
	return max;
}
void Print_Matrix(double A[][3])
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << A[i][j] << "  ";
		}
		cout << endl;
	}
}
void Print_vector(double y[])
{
	double max =Print_x(y);
	cout << "其向量为:" << endl;
	for (int i = 0; i < n; i++) {
		cout << y[i] / max << endl;
	}
}

//不能用LU分解--不满足条件
//列主元消去法----有缺陷
//列主元消去法
void ColGauss(double(*matrix)[n], double b[], double* x)
{

	//遍历矩阵的每一列
	for (int j = 0; j < n - 1; j++) {
		double max = matrix[j][j];
		int idnex = j;

		////记录每一列的最大值
		for (int i = j + 1; i < m; i++)
		{
			if (abs(matrix[i][j]) >= max) {
				max = matrix[i][j];
				idnex = i;
			}
		}

		//b也要替换
		double temp1 = b[j];
		b[j] = b[idnex];
		b[idnex] = temp1;

		//交换第i-1行与，最大值的一行
		for (int k = 0; k < n; k++)
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

				matrix[i][k] = matrix[i][k] - matrix[j][k] * p;//这里改了！！！！！！！！！！！[j]
				//cout << matrix[i][k] << endl;

			}
			b[i] = b[i] - b[j] * p;//这个p也要替换,这里改了！！！！！！！！！！！[j]

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
}
//有缺陷--------不收敛，因为模最大的特征值大于1
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
			y[i] = (b[i] - sum) / matrix[i][i]; //-----------y太大了
			
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

				cout << x[i] << endl;
				
			}
		}
	}
	return k;
}
//注意原点平移的假设

//幂法
int eigen(double matrix[][3],double *y,double *x)
{
	double error = 0;
	int k=0;

	//更新向量x
	while(1)
	{
		//矩阵乘向量,更新向量x
		for (int i = 0; i < m; i++)
		{
			x[i]= 0;
			for (int j = 0; j < n; j++){
				x[i] = x[i] + matrix[i][j] * y[j];
			}

		}

		//记录最大的特征值max
		double max =x[0];

		for (int i = 1; i < n; i++){
			if (abs(x[i]) >= abs(max)){
				max = x[i];  //取模最大的分量
			}
		}
		cout << max << endl;
		k = k + 1;

		//更新向量y[i]
		for (int i = 0; i < n; i++){
			y[i] = x[i]/max;
		}


		if (abs(error -max) < 1e-4){//注意要加绝对值
			cout << max << endl;
			break;
		}
		else
		{
			error = max;
		}

	}
	return k;
}

//幂法加速
int acce_eigen(double matrix[][3], double* y, double* x)
{
	double p = 1.5;
	double error;

	//更新矩阵
	for (int i = 0; i < 3; i++)
	{
		matrix[i][i] = matrix[i][i] - p;
	}

	int k = 0;
	//更新向量x
	while (1)
	{


		//矩阵乘向量
		for (int i = 0; i < m; i++)
		{
			x[i] = 0;
			for (int j = 0; j < n; j++) {
				x[i] = x[i] + matrix[i][j] * y[j];
			}
			

		}
		//归一化
		double max = x[0];
		for (int i = 1; i < n; i++) {
			if (x[i] >= max) {
				max = x[i];
			}
		}
		k = k + 1;
		//迭代误差
		error = 0;
		for (int i = 0; i < n; i++) {
			
			
			error = error + abs(x[i] / max - y[i]);//注意迭代条件
			y[i] = x[i] / max; //这是特征向量
			
		}
	
		if (error < 1e-4) {
			cout << max << endl;
			break;
		}
	}
	return k;
}

//反幂法
int inv_eigen(double matrix[][3], double* y, double* x)
{
	double error = 0;
	int k = 0;
	double A[3][3] = {
   {1,1,0.5},
   {1,1,0.25},
   {0.5,0.25,2}
	};
	//更新向量x
	while (1)
	{
		ColGauss(matrix, y, x);//得到解x
				//更新矩阵
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++) {
				matrix[i][j] = A[i][j];
			}
		}
		//记录最大的特征值max
		double max = x[0];

		for (int i = 1; i < n; i++) {
			if (abs(x[i]) >= abs(max)) {
				max = x[i];  //取模最大的分量
				
			}
		}
		k = k + 1;

		//更新向量y[i]
		for (int i = 0; i < n; i++) {
			y[i] = x[i] / max;
		}
		if (abs(error - max) < 1e-4) {//注意要加绝对值
			break;
		}
		else
		{
			error = max;
		}

	}
	return k;
}

//反幂法加速
int acc_inv_eigen(double matrix[][3], double* y, double* x)
{
	double error = 0;
	int k = 0, p = 1.5;
	double A[3][3] = {
   {1,1,0.5},
   {1,1,0.25},
   {0.5,0.25,2}
	};

	//更新矩阵
	for (int i = 0; i < 3; i++)
	{
		matrix[i][i] = matrix[i][i] - p;
	}
	//更新向量x
	while (1)
	{
		ColGauss(matrix, y, x);//得到解x
				//更新矩阵
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++) {
				matrix[i][j] = A[i][j];
			}
		}
		//记录最大的特征值max
		double max = x[0];

		for (int i = 1; i < n; i++) {
			if (abs(x[i]) >= abs(max)) {
				max = x[i];  //取模最大的分量

			}
		}
		k = k + 1;

		//更新向量y[i]
		for (int i = 0; i < n; i++) {
			y[i] = x[i] / max;
		}
		if (abs(error - max) < 1e-4) {//注意要加绝对值
			break;
		}
		else
		{
			error = max;
		}

	}
	return k;
}
void Jacobi_in(double matrix[][3])
{
	double C[3][3];
	double I[3][3] = {//这里的I不能放进去里面，因为每次都会被赋值为单位矩阵
		{1,0,0},
		{0,1,0},
		{0,0,1}
		};
	while (1) {
		double max = 0;
		int x_i;
		int y_j;
		

		//寻找矩阵的最大元素
		for (int i = 0; i < 3; i++)
		{

			for (int j = 0; j < 3; j++)
			{
				if (abs(max) < abs(matrix[i][j]) && i != j )//可能有0出现，其他为负数，就会返回0了,加绝对值--是大于0,而且是两个绝对值
				{
					max = matrix[i][j];
					y_j = j;
					x_i = i;
					
				}
			}
		}
		//cout << max << endl;
		if (abs(max) <1e-5) //迭代停止条件不是等于0，而是绝对值趋于0
		{
			break;
		}

		//矩阵P
		double P[][3] = {
			{1,0,0},
			{0,1,0},
			{0,0,1}
		};

		//对（i，j）进行计算
		double r = (matrix[x_i][x_i] - matrix[y_j][y_j]) / (2 * matrix[x_i][y_j]); //注意，除以要加括号
		double t = (r > 0 ? 1 : -1) / (abs(r) + sqrt(1 + r * r));
		double c = 1 / sqrt(1 + t * t); //cos
		double s = c * t;   //sin


		//赋值
		P[x_i][x_i] = c;
		P[x_i][y_j] = -s;
		P[y_j][x_i] = s;
		P[y_j][y_j] = c;
		


		//记录特征向量--C=IP
		
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				double sum = 0;
				for (int k = 0; k < 3; k++)
				{
					sum = sum + I[i][k] * P[k][j]; //注意矩阵乘法的规则
				}
				C[i][j] = sum;
			}
		}


		
		double A[3][3];
		//矩阵相乘--AP
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				double sum = 0;
				for (int k = 0; k < 3; k++)
				{
					sum = sum + matrix[i][k] * P[k][j];
				}
				A[i][j] = sum;
				I[i][j] = C[i][j];
			}
		}
		//创建矩阵P^t

		P[x_i][y_j] = s;
		P[y_j][x_i] = -s;
		
		//矩阵相乘--PA
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				double sum = 0;
				for (int k = 0; k < 3; k++)
				{
					sum = sum + P[i][k] * A[k][j];
				}
				matrix[i][j] = sum;
			}
		}

	}
	Print_Matrix(C);

}



int main()
{
	double A[3][3] = {
		{1,1,0.5},
		{1,1,0.25},
		{0.5,0.25,2}
	};
	double B[3][3] = {
	{-0.5,1,0.5},
	{1,-0.5,0.25},
	{0.5,0.25,0.5}
	};
	double x0[3] = { 1,1,1 };//----b
	double x[3] = { 1,1,1 };
	//cout << "k = " << eigen(B, x0, x) << endl;
	//cout << "k = " << eigen(B, x0, x) << endl;
	//cout << "k = " << inv_eigen(A, x0, x) << endl;
	//Print_x(x);
	//Print_vector(x);
	Jacobi_in(A);
	Print_Matrix(A);


	system("pause");
	return 0;
}
