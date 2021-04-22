#include<iostream>
using namespace std;

//加载opencv头文件
#include <highgui.h>
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
using namespace cv;

//Cubic权重函数
float cal_Cubic_coeff(float x)
{

	float a = -0.5; //权重参数
	if (abs(x) <= 1.0) {
		return (a + 2) * pow(abs(x), 3) - (a + 3) * pow(abs(x), 2) + 1;
	}
	else if (abs(x) < 2.0)
	{
		return a * pow(abs(x), 3) - 5 * a * pow(abs(x), 2) + 8 * a * abs(x) - 4 * a;
	}
	else
		return 0.0;
}

Mat BiCubicInter(Mat& src, double sx, double sy)
{
	//获取输出图像的分辨率
	int nRows = cvRound(src.rows * sx);
	int nCols = cvRound(src.cols * sy);

	int n = src.channels();

	Mat resultImage(nRows, nCols, src.type());
	if ( n== 1)
	{
		//遍历图像
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				//获取目标图像(i,j)在原图中的坐标
				int xm = i / sx;
				int ym = j / sy;

				//取出映射到原图的整数部分
				int xi = (int)xm;
				int yi = (int)ym;

				//取出当前点映射到原图中的点的四周的16个点的坐标
				int x0 = xi - 1;
				int y0 = yi - 1;
				int x1 = xi;
				int y1 = yi;
				int x2 = xi + 1;
				int y2 = yi + 1;
				int x3 = xi + 2;
				int y3 = yi + 2;

				//防止当前位置在原图像中越界
				if ((x0 >= 0) && (x3 < src.rows) && (y0 >= 0) && (y3 < src.cols))
				{
					//f(i+u,j+v)=A*B*C
					//求出向量A
					float dist_x0 = cal_Cubic_coeff(xm - x0);
					float dist_x1 = cal_Cubic_coeff(xm - x1);
					float dist_x2 = cal_Cubic_coeff(xm - x2);
					float dist_x3 = cal_Cubic_coeff(xm - x3);

					//求出向量B
					float dist_y0 = cal_Cubic_coeff(ym - y0);
					float dist_y1 = cal_Cubic_coeff(ym - y1);
					float dist_y2 = cal_Cubic_coeff(ym - y2);
					float dist_y3 = cal_Cubic_coeff(ym - y3);

					//计算离当前点的像素最近的16个点的权重
					float dist_x0y0 = dist_x0 * dist_y0;
					float dist_x0y1 = dist_x0 * dist_y1;
					float dist_x0y2 = dist_x0 * dist_y2;
					float dist_x0y3 = dist_x0 * dist_y3;
					float dist_x1y0 = dist_x1 * dist_y0;
					float dist_x1y1 = dist_x1 * dist_y1;
					float dist_x1y2 = dist_x1 * dist_y2;
					float dist_x1y3 = dist_x1 * dist_y3;
					float dist_x2y0 = dist_x2 * dist_y0;
					float dist_x2y1 = dist_x2 * dist_y1;
					float dist_x2y2 = dist_x2 * dist_y2;
					float dist_x2y3 = dist_x2 * dist_y3;
					float dist_x3y0 = dist_x3 * dist_y0;
					float dist_x3y1 = dist_x3 * dist_y1;
					float dist_x3y2 = dist_x3 * dist_y2;
					float dist_x3y3 = dist_x3 * dist_y3;

					//由公式知f(i+u,j+v)=A*B*C知，即计算当前点的像素值
					resultImage.at<Vec3b>(i, j) = (src.at<Vec3b>(x0, y0) * dist_x0y0 +
						src.at<Vec3b>(x0, y1) * dist_x0y1 +
						src.at<Vec3b>(x0, y2) * dist_x0y2 +
						src.at<Vec3b>(x0, y3) * dist_x0y3 +
						src.at<Vec3b>(x1, y0) * dist_x1y0 +
						src.at<Vec3b>(x1, y1) * dist_x1y1 +
						src.at<Vec3b>(x1, y2) * dist_x1y2 +
						src.at<Vec3b>(x1, y3) * dist_x1y3 +
						src.at<Vec3b>(x2, y0) * dist_x2y0 +
						src.at<Vec3b>(x2, y1) * dist_x2y1 +
						src.at<Vec3b>(x2, y2) * dist_x2y2 +
						src.at<Vec3b>(x2, y3) * dist_x2y3 +
						src.at<Vec3b>(x3, y0) * dist_x3y0 +
						src.at<Vec3b>(x3, y1) * dist_x3y1 +
						src.at<Vec3b>(x3, y2) * dist_x3y2 +
						src.at<Vec3b>(x3, y3) * dist_x3y3);
				}
			}
		}
	}
	else
	{
		//遍历图像
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				//获取目标图像(i,j)在原图中的坐标
				int xm = i / sx;
				int ym = j / sy;

				//取出映射到原图的整数部分
				int xi = (int)xm;
				int yi = (int)ym;

				//取出当前点映射到原图中的点的四周的16个点的坐标
				int x0 = xi - 1;
				int y0 = yi - 1;
				int x1 = xi;
				int y1 = yi;
				int x2 = xi + 1;
				int y2 = yi + 1;
				int x3 = xi + 2;
				int y3 = yi + 2;

				//防止当前位置在原图像中越界
				if ((x0 >= 0) && (x3 < src.rows) && (y0 >= 0) && (y3 < src.cols))
				{
					//f(i+u,j+v)=A*B*C
					//求出向量A
					float dist_x0 = cal_Cubic_coeff(xm - x0);
					float dist_x1 = cal_Cubic_coeff(xm - x1);
					float dist_x2 = cal_Cubic_coeff(xm - x2);
					float dist_x3 = cal_Cubic_coeff(xm - x3);

					//求出向量B
					float dist_y0 = cal_Cubic_coeff(ym - y0);
					float dist_y1 = cal_Cubic_coeff(ym - y1);
					float dist_y2 = cal_Cubic_coeff(ym - y2);
					float dist_y3 = cal_Cubic_coeff(ym - y3);

					//计算离当前点的像素最近的16个点的权重
					float dist_x0y0 = dist_x0 * dist_y0;
					float dist_x0y1 = dist_x0 * dist_y1;
					float dist_x0y2 = dist_x0 * dist_y2;
					float dist_x0y3 = dist_x0 * dist_y3;
					float dist_x1y0 = dist_x1 * dist_y0;
					float dist_x1y1 = dist_x1 * dist_y1;
					float dist_x1y2 = dist_x1 * dist_y2;
					float dist_x1y3 = dist_x1 * dist_y3;
					float dist_x2y0 = dist_x2 * dist_y0;
					float dist_x2y1 = dist_x2 * dist_y1;
					float dist_x2y2 = dist_x2 * dist_y2;
					float dist_x2y3 = dist_x2 * dist_y3;
					float dist_x3y0 = dist_x3 * dist_y0;
					float dist_x3y1 = dist_x3 * dist_y1;
					float dist_x3y2 = dist_x3 * dist_y2;
					float dist_x3y3 = dist_x3 * dist_y3;

					for (int k = 0; k < 3; k++) {
						//由公式知f(i+u,j+v)=A*B*C知，即计算当前点的像素值
						resultImage.at<Vec3b>(i, j)[k] = (src.at<Vec3b>(x0, y0)[k] * dist_x0y0 +
							src.at<Vec3b>(x0, y1)[k] * dist_x0y1 +
							src.at<Vec3b>(x0, y2)[k] * dist_x0y2 +
							src.at<Vec3b>(x0, y3)[k] * dist_x0y3 +
							src.at<Vec3b>(x1, y0)[k] * dist_x1y0 +
							src.at<Vec3b>(x1, y1)[k] * dist_x1y1 +
							src.at<Vec3b>(x1, y2)[k] * dist_x1y2 +
							src.at<Vec3b>(x1, y3)[k] * dist_x1y3 +
							src.at<Vec3b>(x2, y0)[k] * dist_x2y0 +
							src.at<Vec3b>(x2, y1)[k] * dist_x2y1 +
							src.at<Vec3b>(x2, y2)[k] * dist_x2y2 +
							src.at<Vec3b>(x2, y3)[k] * dist_x2y3 +
							src.at<Vec3b>(x3, y0)[k] * dist_x3y0 +
							src.at<Vec3b>(x3, y1)[k] * dist_x3y1 +
							src.at<Vec3b>(x3, y2)[k] * dist_x3y2 +
							src.at<Vec3b>(x3, y3)[k] * dist_x3y3);
					}
				}
			}
		}
	}

	return resultImage;
}
int main()
{
	//图像沿x轴和y轴放大的倍数
	double sx = 0.5;
	double sy = 0.5;

	Mat src = imread("E://flower.jpg");
	if (!src.data)
	{
		printf("image could not load...\n");
		return -1;
	}

	//创建自适应的窗口--读取原图像---展示原图像
	namedWindow("src", CV_WINDOW_NORMAL);
	imshow("src", src);

	Mat dst = BiCubicInter(src, sx, sy);//调用函数、输出目标图像

	//创建自适应的窗口--读取原图像---展示原图像
	namedWindow("dst", CV_WINDOW_NORMAL);
	imshow("src", dst);

	imwrite("E://flower_hermite_low1.jpg", dst); // 保存图像dst到..

	waitKey(0);
	return 0;

}
