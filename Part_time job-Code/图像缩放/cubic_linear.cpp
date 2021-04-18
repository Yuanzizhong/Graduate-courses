
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
using namespace cv;
using namespace std;
#define PI 3.14159265


float BiCubicPoly(float x);
void MyScaleBiCubicInter(Mat& src, Mat& dst, float TransMat[3][3]);
/**
 * @function main
 */
int main(int argc, char** argv)
{
	// load image
	Mat image;
	image = imread("E://4.PNG");

	if (!image.data)
	{
		cout << "No image data" << endl;
		return -1;
	}
	// show image
	namedWindow("image", CV_WINDOW_AUTOSIZE);
	imshow("image", image);
	Mat dst;
	float transMat[3][3] = { {3.0, 0, 0}, {0, 3.0, 0}, {0, 0, 1} };

	MyScaleBiCubicInter(image, dst, transMat);
	namedWindow("out_image", CV_WINDOW_AUTOSIZE);
	imshow("out_image", dst);
	imwrite("E://big_hermite.jpg", dst);
	waitKey(0);
	return 0;
}
float BiCubicPoly(float x)
{
	float abs_x = abs(x);
	float a = -1.0;
	if (abs_x <= 1.0)
	{
		return (a + 2) * pow(abs_x, 3) - (a + 3) * pow(abs_x, 2) + 1;
	}
	else if (abs_x < 2.0)
	{
		return a * pow(abs_x, 3) - 5 * a * pow(abs_x, 2) + 8 * a * abs_x - 4 * a;
	}
	else
		return 0.0;
}

void MyScaleBiCubicInter(Mat& src, Mat& dst, float TransMat[3][3])
{
	CV_Assert(src.data);
	CV_Assert(src.depth() != sizeof(uchar));

	// calculate margin point of dst image
	float left = 0;
	float right = 0;
	float top = 0;
	float down = 0;

	float x = src.cols * 1.0f;
	float y = 0.0f;
	float u1 = x * TransMat[0][0] + y * TransMat[0][1];
	float v1 = x * TransMat[1][0] + y * TransMat[1][1];
	x = src.cols * 1.0f;
	y = src.rows * 1.0f;
	float u2 = x * TransMat[0][0] + y * TransMat[0][1];
	float v2 = x * TransMat[1][0] + y * TransMat[1][1];
	x = 0.0f;
	y = src.rows * 1.0f;
	float u3 = x * TransMat[0][0] + y * TransMat[0][1];
	float v3 = x * TransMat[1][0] + y * TransMat[1][1];

	left = min(min(min(0.0f, u1), u2), u3);
	right = max(max(max(0.0f, u1), u2), u3);
	top = min(min(min(0.0f, v1), v2), v3);
	down = max(max(max(0.0f, v1), v2), v3);

	// create dst image
	dst.create(int(abs(right - left)), int(abs(down - top)), src.type());

	CV_Assert(dst.channels() == src.channels());
	int channels = dst.channels();

	int i, j;
	uchar* p;
	uchar* q0;
	uchar* q1;
	uchar* q2;
	uchar* q3;
	for (i = 0; i < dst.rows; ++i)
	{
		p = dst.ptr<uchar>(i);
		for (j = 0; j < dst.cols; ++j)
		{
			// 
			x = (j + left) / TransMat[0][0];
			y = (i + top) / TransMat[1][1];

			int x0 = int(x) - 1;
			int y0 = int(y) - 1;
			int x1 = int(x);
			int y1 = int(y);
			int x2 = int(x) + 1;
			int y2 = int(y) + 1;
			int x3 = int(x) + 2;
			int y3 = int(y) + 2;

			if ((x0 >= 0) && (x3 < src.cols) && (y0 >= 0) && (y3 < src.rows))
			{
				q0 = src.ptr<uchar>(y0);
				q1 = src.ptr<uchar>(y1);
				q2 = src.ptr<uchar>(y2);
				q3 = src.ptr<uchar>(y3);

				float dist_x0 = BiCubicPoly(x - x0);
				float dist_x1 = BiCubicPoly(x - x1);
				float dist_x2 = BiCubicPoly(x - x2);
				float dist_x3 = BiCubicPoly(x - x3);
				float dist_y0 = BiCubicPoly(y - y0);
				float dist_y1 = BiCubicPoly(y - y1);
				float dist_y2 = BiCubicPoly(y - y2);
				float dist_y3 = BiCubicPoly(y - y3);

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

				switch (channels)
				{
				case 1:
				{
					break;
				}
				case 3:
				{
					p[3 * j] = (uchar)(q0[3 * x0] * dist_x0y0 +
						q1[3 * x0] * dist_x0y1 +
						q2[3 * x0] * dist_x0y2 +
						q3[3 * x0] * dist_x0y3 +
						q0[3 * x1] * dist_x1y0 +
						q1[3 * x1] * dist_x1y1 +
						q2[3 * x1] * dist_x1y2 +
						q3[3 * x1] * dist_x1y3 +
						q0[3 * x2] * dist_x2y0 +
						q1[3 * x2] * dist_x2y1 +
						q2[3 * x2] * dist_x2y2 +
						q3[3 * x2] * dist_x2y3 +
						q0[3 * x3] * dist_x3y0 +
						q1[3 * x3] * dist_x3y1 +
						q2[3 * x3] * dist_x3y2 +
						q3[3 * x3] * dist_x3y3);

					p[3 * j + 1] = (uchar)(q0[3 * x0 + 1] * dist_x0y0 +
						q1[3 * x0 + 1] * dist_x0y1 +
						q2[3 * x0 + 1] * dist_x0y2 +
						q3[3 * x0 + 1] * dist_x0y3 +
						q0[3 * x1 + 1] * dist_x1y0 +
						q1[3 * x1 + 1] * dist_x1y1 +
						q2[3 * x1 + 1] * dist_x1y2 +
						q3[3 * x1 + 1] * dist_x1y3 +
						q0[3 * x2 + 1] * dist_x2y0 +
						q1[3 * x2 + 1] * dist_x2y1 +
						q2[3 * x2 + 1] * dist_x2y2 +
						q3[3 * x2 + 1] * dist_x2y3 +
						q0[3 * x3 + 1] * dist_x3y0 +
						q1[3 * x3 + 1] * dist_x3y1 +
						q2[3 * x3 + 1] * dist_x3y2 +
						q3[3 * x3 + 1] * dist_x3y3);

					p[3 * j + 2] = (uchar)(q0[3 * x0 + 2] * dist_x0y0 +
						q1[3 * x0 + 2] * dist_x0y1 +
						q2[3 * x0 + 2] * dist_x0y2 +
						q3[3 * x0 + 2] * dist_x0y3 +
						q0[3 * x1 + 2] * dist_x1y0 +
						q1[3 * x1 + 2] * dist_x1y1 +
						q2[3 * x1 + 2] * dist_x1y2 +
						q3[3 * x1 + 2] * dist_x1y3 +
						q0[3 * x2 + 2] * dist_x2y0 +
						q1[3 * x2 + 2] * dist_x2y1 +
						q2[3 * x2 + 2] * dist_x2y2 +
						q3[3 * x2 + 2] * dist_x2y3 +
						q0[3 * x3 + 2] * dist_x3y0 +
						q1[3 * x3 + 2] * dist_x3y1 +
						q2[3 * x3 + 2] * dist_x3y2 +
						q3[3 * x3 + 2] * dist_x3y3);

					float thre = 198.0f;
					if ((abs(p[3 * j] - q1[3 * x1]) > thre) || (abs(p[3 * j + 1] - q1[3 * x1 + 1]) > thre) ||
						(abs(p[3 * j + 2] - q1[3 * x1 + 2]) > thre))
					{
						p[3 * j] = q1[3 * x1];
						p[3 * j + 1] = q1[3 * x1 + 1];
						p[3 * j + 2] = q1[3 * x1 + 2];
					}
					break;
				}
				}
			}
		}
	}
}
