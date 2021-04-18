#include<iostream>
using namespace std;

#include <highgui.h>

#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
using namespace cv; //调用imshow


//最近邻插值算法
//sx,xy -- 缩放因子
void nearest(cv::Mat& src, cv::Mat& dst, float sx, float sy) {
	//由缩放因子计算输出图像的尺寸（四舍五入）
	int dst_cols = round(src.cols * sx);
	int dst_rows = round(src.rows * sy);
	//创建输出图像
	dst = cv::Mat(dst_rows, dst_cols, src.type());
	//灰度图像处理
	if (src.channels() == 1) {
		for (int i = 0; i < dst.rows; i++) {
			for (int j = 0; j < dst.cols; j++) {
				//插值计算，输出图像的像素点由原图像对应的最近的像素点得到（四舍五入）
				int i_index = round(i / sy);
				int j_index = round(j / sx);
				if (i_index > src.rows - 1) i_index = src.rows - 1;//防止越界
				if (j_index > src.cols - 1) j_index = src.cols - 1;//防止越界
				dst.at<uchar>(i, j) = src.at<uchar>(i_index, j_index);
			}
		}
	}
	//彩色图像处理
	else {
		for (int i = 0; i < dst.rows; i++) {
			for (int j = 0; j < dst.cols; j++) {
				//插值计算，输出图像的像素点由原图像对应的最近的像素点得到（四舍五入）
				int i_index = round(i / sy);
				int j_index = round(j / sx);
				if (i_index > src.rows - 1) i_index = src.rows - 1;//防止越界
				if (j_index > src.cols - 1) j_index = src.cols - 1;//防止越界
				//B
				dst.at<cv::Vec3b>(i, j)[0] = src.at<cv::Vec3b>(i_index, j_index)[0];
				//G
				dst.at<cv::Vec3b>(i, j)[1] = src.at<cv::Vec3b>(i_index, j_index)[1];
				//R
				dst.at<cv::Vec3b>(i, j)[2] = src.at<cv::Vec3b>(i_index, j_index)[2];
			}
		}
	}
}

//双线性插值
////sx,xy -- 缩放因子
void Inter_Linear(cv::Mat& src, cv::Mat& dst, double sx, double sy) {
	int dst_rows = round(sx * src.rows);
	int dst_cols = round(sy * src.cols);
	dst = cv::Mat(dst_rows, dst_cols, src.type());
	for (int i = 0; i < dst.rows; i++) {
		//几何中心对齐
		double index_i = (i + 0.5) / sx - 0.5;
		//防止越界
		if (index_i < 0) index_i = 0;
		if (index_i > src.rows - 1) index_i = src.rows - 1;
		//相邻4*4像素的行（坐标）
		int i1 = floor(index_i);
		int i2 = ceil(index_i);
		//u为得到浮点型坐标行的小数部分
		double u = index_i - i1;
		for (int j = 0; j < dst.cols; j++) {
			//几何中心对齐
			double index_j = (j + 0.5) / sy - 0.5;
			//防止越界
			if (index_j < 0) index_j = 0;
			if (index_j > src.cols - 1) index_j = src.cols - 1;
			//相邻4*4像素的列（坐标）
			int j1 = floor(index_j);
			int j2 = ceil(index_j);
			//v为得到浮点型坐标列的小数部分
			double v = index_j - j1;
			if (src.channels() == 1) {
				//灰度图像
				dst.at<uchar>(i, j) = (1 - u) * (1 - v) * src.at<uchar>(i1, j1) + (1 - u) * v * src.at<uchar>(i1, j2) + u * (1 - v) * src.at<uchar>(i2, j1) + u * v * src.at<uchar>(i2, j2);
			}
			else {
				//彩色图像
				dst.at<cv::Vec3b>(i, j)[0] = (1 - u) * (1 - v) * src.at<cv::Vec3b>(i1, j1)[0] + (1 - u) * v * src.at<cv::Vec3b>(i1, j2)[0] + u * (1 - v) * src.at<cv::Vec3b>(i2, j1)[0] + u * v * src.at<cv::Vec3b>(i2, j2)[0];
				dst.at<cv::Vec3b>(i, j)[1] = (1 - u) * (1 - v) * src.at<cv::Vec3b>(i1, j1)[1] + (1 - u) * v * src.at<cv::Vec3b>(i1, j2)[1] + u * (1 - v) * src.at<cv::Vec3b>(i2, j1)[1] + u * v * src.at<cv::Vec3b>(i2, j2)[1];
				dst.at<cv::Vec3b>(i, j)[2] = (1 - u) * (1 - v) * src.at<cv::Vec3b>(i1, j1)[2] + (1 - u) * v * src.at<cv::Vec3b>(i1, j2)[2] + u * (1 - v) * src.at<cv::Vec3b>(i2, j1)[2] + u * v * src.at<cv::Vec3b>(i2, j2)[2];
			}
		}
	}
}


int main()
{
	double sx =0.5;
	double sy = 0.5;
	Mat src, dst;

	src = imread("E://4.png");

	if (src.empty())
	{
		std::cout << "图片读取失败！" << "\n";
		return -1;
	}

	namedWindow("src", CV_WINDOW_NORMAL);
	imshow("src", src);

	//nearest(src, dst, sx, sy);
	Inter_Linear(src,dst, sx, sy);
	if (dst.empty())
	{
		std::cout << "图片读取失败！" << "\n";
		return -1;
	}
	imwrite("E://big_xiao.jpg", dst);
	
	namedWindow("dst", CV_WINDOW_NORMAL);
	imshow("dst", dst);
	waitKey(0);
	//system("pause");//这个函数有问题
	return 0;
}
