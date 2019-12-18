///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


//////////////////////////  CONSTRUCTORS/DESTRUCTOR ////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}

///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage



///////////////////////////// IN PROGRESS //////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
	/*
	int even_even[3][3] = {
		{1, 2, 1},
		{2, 4, 2},
		{1, 2, 1}
	}; int even_even_div = 16;

	int even_odd[4][3] = {
		{1, 2, 1},
		{3, 6, 3},
		{3, 6, 3},
		{1, 2, 1}
	}; int even_odd_div = 32;

	int odd_even[3][4] = {
		{1, 3, 3, 1},
		{2, 6, 6, 2},
		{1, 3, 3, 1},
	}; int odd_even_div = 32;

	int odd_odd[4][4] = {
		{1, 3, 3, 1},
		{3, 9, 9, 3},
		{3, 9, 9, 3},
		{1, 3, 3, 1}
	};  int odd_odd_div = 64;

	int nheight = height * 2;
	int nwidth = width * 2;
	unsigned char * ndata = new unsigned char[nheight * nwidth * 4];

	//stepping through new image.
	for (int i = 0; i < nheight; i++) {
		int os = i * nwidth * 4;
		for (int j = 0; j < nwidth; j++) {
			int x = j * 4;
			for (int k = 0; k < 4; k++) {
				int row, col;
				int pixel = 0;
				int which_case;

				if (i % 2 == 0) {
					if (j % 2 == 0) {
						which_case = 1;
					}
					else {
						which_case = 2;
					}
				}
				else {
					if (j % 2 == 0) {
						which_case = 3;
					}
					else {
						which_case = 4;
					}
				}
				
				switch (which_case) {
				case 1: { //EVEN_EVEN

				}
				case 2: { //EVEN_ODD

				}
				case 3: { //ODD_EVEN

				}
				case 4: { //ODD_ODD

				}
				}


			}
		}
	}
	*/

	ClearToBlack();
	return false;
}// Double_Size



/////////////////////////////////// TO-DO //////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    ClearToBlack();
    return false;
}// Dither_Random

///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    ClearToBlack();
    return false;
}// Dither_Cluster

///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    ClearToBlack();
    return false;
}// Dither_Color

///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
	int nwidth = (int)floor(width * scale);
	int nheight = (int)floor(width * scale);
	unsigned char * ndata = new unsigned char[nwidth * nheight * 4];

	double new_x, new_y;


	int i, j, k, l, m, os, x;

	for (i = 0; i < nheight; i++) {
		os = i * nwidth * 4;
		for (j = 0; j < nwidth; j++) {
			x = j * 4;
			for (k = 0; k < 3; k++) {
				int pixel = 0;
				//name the point we're trying to calculate;
				new_x = j;
				new_y = i;
				//we will intersect exactly 16 points in the old image with the bartlett filter;
				//so we need to start stepping through the points and combining their values;
				for (l = 0; l < 4; l++) {
					for (m = 0; m < 4; m++) {
						double coeff = F((new_x / scale) - 2 + m, (new_y / scale) - 2 + l);
						int row = (int)floor(new_y / scale) - 2 + l;
						int col = (int)floor(new_x / scale) - 2 + m;
						if (row < 0) {
							row = -row;
						}
						else if (row >= height) {
							row = height - (row - height) - 1;
						}
						if (col < 0) {
							col = -col;
						}
						else if (col >= width) {
							col = width - (col - width) - 1;
						}
						row = row * width * 4;
						col = col * 4;
						int intensity = data[row + col + k];
						intensity *= coeff;
						pixel += intensity;
					}
				}
				ndata[os + x + k] = pixel;
			}
			ndata[os + x + 3] = 255;
		}
	}
	delete[] data;
	data = ndata;
	width = nwidth;
	height = nheight;

	return true;
}// Resize

//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


///////////////////////////// COMPLETE //////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
//	NOTE:	Currently does not look at the border of the image.
//			Will need to fix this later. 
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
	unsigned char * filtered = new unsigned char[width * height * 4];
	int i, j, k, l, m, os, x, pixel, row, col;
	//for each row
	for (i = 0; i < height; i++) {
		os = i * width * 4;
		//for each column
		for (j = 0; j < width; j++) {
			x = j * 4;
			//for each color
			for (k = 0; k < 3; k++) {
				pixel = 0;
				//for each subrow
				for (l = -2; l < 3; l++) {
					//for each subcolumn
					for (m = -2; m < 3; m++) {
						row = (i + l);
						col = (j + m);
						if (row < 0) {
							row = -row;
						}
						else if (row >= height) {
							row = height - (row - height) - 1;
						}
						if (col < 0) {
							col = -col;
						}
						else if (col >= width) {
							col = width - (col - width) - 1;
						}
						row = row * width * 4;
						col = col * 4;
						pixel += data[row + col + k];
					}
				}
				pixel /= 25;
				filtered[os + x + k] = pixel;
			}
			filtered[os + x + 3] = data[os + x + 3];
		}
	}
	delete[] data;
	data = filtered;
	return true;
}// Filter_Box

///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
	unsigned char * filtered = new unsigned char[width*height * 4];
	int filter[5][5] = {
		{1, 2, 3, 2, 1 },
		{2, 4, 6, 4, 2 },
		{3, 6, 9, 6, 3 },
		{2, 4, 6, 4, 2 },
		{1, 2, 3, 2, 1 }
	};

	int i, j, k, l, m, os, x, pixel, row, col;

	for (i = 0; i < height; i++) {
		os = i * width * 4;
		for (j = 0; j < width; j++) {
			x = j * 4;
			for (k = 0; k < 3; k++) {
				pixel = 0;
				for (l = -2; l < 3; l++) {
					for (m = -2; m < 3; m++) {
						row = (i + l);
						col = (j + m);
						if (row < 0) {
							row = -row;
						}
						else if (row >= height) {
							row = height - (row - height) - 1;
						}
						if (col < 0) {
							col = -col;
						}
						else if (col >= width) {
							col = width - (col - width) - 1;
						}
						row = row * width * 4;
						col = col * 4;
						pixel += data[row + col + k] * filter[l + 2][m + 2];
					}
				}
				pixel /= 81;
				filtered[os + x + k] = pixel;
			}
			filtered[os + x + 3] = data[os + x + 3];
		}
	}
	delete[] data;
	data = filtered;
	return true;
}// Filter_Bartlett

///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
	this->Filter_Gaussian_N(5);
	return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian_N(unsigned int N)
{
	if (N % 2 == 0) {
		return false;
	}

	unsigned char * filtered = new unsigned char[width * height * 4];

	int* one_d = NULL;
	find_binomial(N, one_d);
	int i, j, k, l, m, os, x, pixel, row, col, lb, ub;
	int mfactor = 1;
	for (i = 0; i < (N - 1); i++) {
		mfactor *= 2;
	}
	mfactor *= mfactor;

	int* mtx = new int[N*N];
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			mtx[i*N + j] = one_d[i] * one_d[j];
		}
	}
	delete[] one_d;
	lb = -1 * ((int)N) / 2;
	ub = (N / 2) + 1;

	for (i = 0; i < height; i++) {
		os = i * width * 4;
		for (j = 0; j < width; j++) {
			x = j * 4;
			for (k = 0; k < 3; k++) {
				pixel = 0;
				for (l = lb; l < ub; l++) {
					for (m = lb; m < ub; m++) {
						row = (i + l);
						col = (j + m);
						if (row < 0) {
							row = -row;
						}
						else if (row >= height) {
							row = height - (row - height) - 1;
						}
						if (col < 0) {
							col = -col;
						}
						else if (col >= width) {
							col = width - (col - width) - 1;
						}
						row = row * width * 4;
						col = col * 4;
						pixel += data[row + col + k] * mtx[((l + (N / 2))*N) + (m + (N / 2))];
					}
				}
				pixel /= mfactor;
				filtered[os + x + k] = pixel;
			}
			filtered[os + x + 3] = data[os + x + 3];
		}
	}
	delete[] data;
	delete[] mtx;
	data = filtered;
	return true;

}// Filter_Gaussian_N

///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
	unsigned char * filtered = new unsigned char[width * height * 4];
	int filter[5][5] = {
		{1, 4, 6, 4, 1},
		{4, 16, 24, 16, 4},
		{6, 24, -220, 24, 6},
		{4, 16, 24, 16, 4},
		{1, 4, 6, 4, 1},
	};
	int filterdiv = -256;

	int i, j, k, l, m, os, x, pixel, row, col;

	for (i = 0; i < height; i++) {
		os = i * width * 4;
		for (j = 0; j < width; j++) {
			x = j * 4;
			for (k = 0; k < 3; k++) {
				pixel = 0;
				for (l = -2; l < 3; l++) {
					for (m = -2; m < 3; m++) {
						row = (i + l);
						col = (j + m);
						if (row < 0) {
							row = -row;
						}
						else if (row >= height) {
							row = height - (row - height) - 1;
						}
						if (col < 0) {
							col = -col;
						}
						else if (col >= width) {
							col = width - (col - width) - 1;
						}
						row = row * width * 4;
						col = col * 4;
						pixel += data[row + col + k] * filter[l + 2][m + 2];
					}
				}
				pixel /= filterdiv;
				if (pixel < 0) {
					pixel = 0;
				}
				else if (pixel > 255) {
					pixel = 255;
				}
				filtered[os + x + k] = pixel;
			}
			filtered[os + x + 3] = data[os + x + 3];
		}
	}
	delete[] data;
	data = filtered;
	return true;

}// Filter_Edge

///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
	TargaImage smoothed(*this);
	smoothed.Filter_Edge();
	unsigned char * filtered = smoothed.data;

	for (int i = 0; i < height; i++) {
		int os = i * width * 4;
		for (int j = 0; j < width; j++) {
			int x = j * 4;
			for (int k = 0; k < 3; k++) {
				int pixel = data[os + x + k] + filtered[os + x + k];
				if (pixel < 0) {
					pixel = 0;
				}
				else if (pixel > 255) {
					pixel = 255;
				}
				data[os + x + k] = pixel;
			}
		}
	}
	return true;
}// Filter_Enhance

///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
	if (width != pImage->width || height != pImage->height)
	{
		cout << "Comp_Over: Images not the same size\n";
		return false;
	}

	unsigned char * fg = this->To_RGB();
	unsigned char * bg = pImage->To_RGB();

	int i, j, k, os3, os4, x3, x4;
	float arr[8];

	for (i = 0; i < height; i++) {
		os3 = i * width * 3;
		os4 = i * width * 4;
		for (j = 0; j < width; j++) {
			x3 = j * 3;
			x4 = j * 4;
			for (k = 0; k < 3; k++) {
				arr[k + 0] = fg[os3 + x3 + k] / 255.0;
				arr[k + 4] = bg[os3 + x3 + k] / 255.0;
			}
			arr[3] = this->data[os4 + x4 + 3] / 255.0;
			arr[7] = pImage->data[os4 + x4 + 3] / 255.0;
			for (k = 0; k < 4; k++) {
				data[os4 + x4 + k] = (int)((arr[k] + (1 - arr[3]) * arr[k + 4]) * 255);
			}
		}
	}
	delete[] fg;
	delete[] bg;
	return true;
}// Comp_Over

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
	if (width != pImage->width || height != pImage->height)
	{
		cout << "Comp_In: Images not the same size\n";
		return false;
	}

	unsigned char * fg = this->To_RGB();
	unsigned char * bg = pImage->To_RGB();
	int i, j, k, os3, os4, x3, x4;
	float arr[8];

	for (i = 0; i < height; i++) {
		os3 = i * width * 3;
		os4 = i * width * 4;
		for (j = 0; j < width; j++) {
			x3 = j * 3;
			x4 = j * 4;
			for (k = 0; k < 3; k++) {
				arr[k + 0] = fg[os3 + x3 + k] / 255.0;
				arr[k + 4] = bg[os3 + x3 + k] / 255.0;
			}
			arr[3] = this->data[os4 + x4 + 3] / 255.0;
			arr[7] = pImage->data[os4 + x4 + 3] / 255.0;
			for (k = 0; k < 4; k++) {
				data[os4 + x4 + k] = (int)((arr[7] * arr[k]) * 255);
			}
		}
	}

	delete[] fg;
	delete[] bg;
	return true;
}// Comp_In

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
	if (width != pImage->width || height != pImage->height)
	{
		cout << "Comp_Out: Images not the same size\n";
		return false;
	}

	unsigned char * fg = this->To_RGB();
	unsigned char * bg = pImage->To_RGB();
	int i, j, k, os3, os4, x3, x4;
	float arr[8];

	for (i = 0; i < height; i++) {
		os3 = i * width * 3;
		os4 = i * width * 4;
		for (j = 0; j < width; j++) {
			x3 = j * 3;
			x4 = j * 4;
			for (k = 0; k < 3; k++) {
				arr[k + 0] = fg[os3 + x3 + k] / 255.0;
				arr[k + 4] = bg[os3 + x3 + k] / 255.0;
			}
			arr[3] = this->data[os4 + x4 + 3] / 255.0;
			arr[7] = pImage->data[os4 + x4 + 3] / 255.0;
			for (k = 0; k < 4; k++) {
				data[os4 + x4 + k] = (int)((1 - arr[7]) * arr[k] * 255);
			}
		}
	}

	delete[] fg;
	delete[] bg;
	return true;
}// Comp_Out

///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
	if (width != pImage->width || height != pImage->height)
	{
		cout << "Comp_Atop: Images not the same size\n";
		return false;
	}

	unsigned char * fg = this->To_RGB();
	unsigned char * bg = pImage->To_RGB();
	int i, j, k, os3, os4, x3, x4;
	float arr[8];

	for (i = 0; i < height; i++) {
		os3 = i * width * 3;
		os4 = i * width * 4;
		for (j = 0; j < width; j++) {
			x3 = j * 3;
			x4 = j * 4;
			for (k = 0; k < 3; k++) {
				arr[k + 0] = fg[os3 + x3 + k] / 255.0;
				arr[k + 4] = bg[os3 + x3 + k] / 255.0;
			}
			arr[3] = this->data[os4 + x4 + 3] / 255.0;
			arr[7] = pImage->data[os4 + x4 + 3] / 255.0;
			for (k = 0; k < 4; k++) {
				data[os4 + x4 + k] = (int)(((arr[7] * arr[k]) + ((1 - arr[3]) * arr[k + 4])) * 255);
			}
		}
	}

	delete[] fg;
	delete[] bg;
	return true;
}// Comp_Atop

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
	if (width != pImage->width || height != pImage->height)
	{
		cout << "Comp_Xor: Images not the same size\n";
		return false;
	}

	unsigned char * fg = this->To_RGB();
	unsigned char * bg = pImage->To_RGB();
	int i, j, k, os3, os4, x3, x4;
	float arr[8];

	for (i = 0; i < height; i++) {
		os3 = i * width * 3;
		os4 = i * width * 4;
		for (j = 0; j < width; j++) {
			x3 = j * 3;
			x4 = j * 4;
			for (k = 0; k < 3; k++) {
				arr[k + 0] = fg[os3 + x3 + k] / 255.0;
				arr[k + 4] = bg[os3 + x3 + k] / 255.0;
			}
			arr[3] = this->data[os4 + x4 + 3] / 255.0;
			arr[7] = pImage->data[os4 + x4 + 3] / 255.0;
			for (k = 0; k < 4; k++) {
				data[os4 + x4 + k] = (int)((((1 - arr[7])*arr[k]) + ((1 - arr[3]) * arr[k + 4])) * 255);
			}
		}
	}

	delete[] fg;
	delete[] bg;
	return true;
}// Comp_Xor

///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
	int i, j;
	for (i = 0; i < height; i++) {
		int os = i * width * 4;
		for (j = 0; j < width; j++) {
			int x = j * 4;
			data[os + x] = data[os + x + 1] = data[os + x + 2]
				= (int)(0.299*data[os + x])
				+ (int)(0.587*data[os + x + 1])
				+ (int)(0.114*data[os + x + 2]);
		}
	}
	return true;
}// To_Grayscale

///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
//		Use the uniform quantization algorithm to convert the current image 
//  from a 24 bit color image to an 8 bit color image. Use 4 levels of blue, 
//	8 levels of red, and 8 levels of green in the quantized image.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
	int i, j;
	for (i = 0; i < height; i++) {
		int os = i * width * 4;
		for (j = 0; j < width; j++) {
			int x = j * 4;
			data[os + x    ] = ((((int)data[os + x    ]) >> 5) << 5);
			data[os + x + 1] = ((((int)data[os + x + 1]) >> 5) << 5);
			data[os + x + 2] = ((((int)data[os + x + 2]) >> 6) << 6);
		}
	}

	return true;
}// Quant_Uniform

///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
//		Use the populosity algorithm to convert the current 24 bit color image 
//	to an 8 bit color image. Before building the color usage histogram, do a 
//	uniform quantization step down to 32 levels of each primary. This gives 
//	32 x 32 x 32 = 32768 possible colors. Then find the 256 most popular colors, 
//	then map the original colors onto their closest chosen color. To find the 
//	closest color, use the euclidean (L2) distance in RGB space. If (r1,g1,b1) 
//	and (r2,g2,b2) are the colors, use sqrt((r1-r2)^2 + (g1-g2)^2 + (b1-b2)^2) 
//	suitably converted into C++ code.
//
//		After quantizing down to 32 levels of each primary, we use bit-shifting
//	to combine them all into a single integer in the range [0, 32767]. We use 
//	the structure "quant_pop_color" to store these colors (along with a counter
//	for each) in an array. We then sort the array, take the 256 most common 
//	colors, and step through the image - comparing each pixel in the image to 
//	each of the 256 most common colors and ultimately assigning each pixel's 
//  color to the closest value from the most common colors. 
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
	quant_pop_color colormap[32768];
	for (int z = 0; z < 32768; z++) {
		colormap[z].value = z;
		colormap[z].count = 0;
	}
	unsigned int colors[256] = { 0 };

	int i, j, r, g, b;
	int value;
	for (i = 0; i < height; i++) {
		int os = i * width * 4;
		for (j = 0; j < width; j++) {
			int x = j * 4;
			r = data[os + x] = ((((int)data[os + x]) >> 3) << 3);
			g = data[os + x + 1] = ((((int)data[os + x + 1]) >> 3) << 3);
			b = data[os + x + 2] = ((((int)data[os + x + 2]) >> 3) << 3);
			value = (r << 7) + (g << 2) + (b >> 3);
			colormap[value].count++;
		}
	}

	sort(colormap, colormap + 32768, [](quant_pop_color a, quant_pop_color b) {return (a.count > b.count); });
	for (i = 0; i < 256; i++) {
		colors[i] = colormap[i].value;
	}

	int cr, cg, cb;
	double frd, fgd, fbd;

	for (i = 0; i < height; i++) {
		int os = i * width * 4;
		for (j = 0; j < width; j++) {
			int x = j * 4;

			r = data[os + x];
			g = data[os + x + 1];
			b = data[os + x + 2];

			unsigned int curr_color;
			double best_dist = numeric_limits<double>::max();
			int best_color;
			double curr_dist;

			for (int k = 0; k < 256; k++) {
				curr_color = colors[k];
				cr = ((curr_color & 0x00007C00) >> 7);
				cg = ((curr_color & 0x000003E0) >> 2);
				cb = ((curr_color & 0x0000001F) << 3);

				frd = (r - cr);
				fgd = (g - cg);
				fbd = (b - cb);

				curr_dist = sqrt(pow(frd, 2) + pow(fgd, 2) + pow(fbd, 2));
				if (curr_dist < best_dist) {
					best_color = curr_color;
					best_dist = curr_dist;
				}
			}
			r = (best_color & 0x00007C00) >> 7;
			g = (best_color & 0x000003E0) >> 2;
			b = (best_color & 0x0000001F) << 3;

			data[os + x] = r;
			data[os + x + 1] = g;
			data[os + x + 2] = b;
		}
	}
	return true;
}// Quant_Populosity

///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
	float *arr = new float[width*height];
	int i, j, k, os, x, t;
	this->To_Grayscale();
	for (i = 0; i < height; i++) {
		os = i * width * 4;
		for (j = 0; j < width; j++) {
			x = j * 4;
			arr[i*width + j] = data[os + x] / (float)256;
			if (arr[i*width + j] < 0.5) {
				t = 0;
			}
			else {
				t = 255;
			}
			for (k = 0; k < 3; k++) {
				data[os + x + k] = t;
			}
		}
	}
	delete[] arr;
	return true;
}// Dither_Threshold
 
///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
	long int intensity_fst = 0;
	long int intensity_sec = 0;

	this->To_Grayscale();
	int i, j, k, os, x, t;
	double intensity = 0;
	double threshold;
	int dimension = width * height;
	double* pixel_intensities = new double[dimension];
	double* sorted_intensities = new double[dimension];

	for (i = 0; i < height; i++) {
		os = i * width * 4;
		for (j = 0; j < width; j++) {
			x = j * 4;
			pixel_intensities[i * width + j] = data[os + x] / 256.0;
			sorted_intensities[i * width + j] = data[os + x] / 256.0;
			intensity_fst += data[os + x];
		}
	}

	for (i = 0; i < dimension; i++) {
		intensity += pixel_intensities[i];
	}

	intensity /= (dimension);
	sort(sorted_intensities, sorted_intensities + dimension);

	int index = (int)((1 - intensity)*dimension);
	threshold = sorted_intensities[index];

	for (i = 0; i < height; i++) {
		os = i * width * 4;
		for (j = 0; j < width; j++) {
			x = j * 4;
			if (pixel_intensities[i * width + j] <= threshold) {
				t = 0;
			}
			else {
				t = 255;
			}
			for (k = 0; k < 3; k++) {
				data[os + x + k] = t;
			}
			intensity_sec += data[os + x];

		}
	}
	delete[] pixel_intensities;
	delete[] sorted_intensities;
	return true;
}// Dither_Bright

///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
	this->To_Grayscale();
	float * fdata = new float[width * height];
	for (int i = 0; i < width*height; i++) {
		fdata[i] = (float)(data[i * 4] / 255.0);
	}
	float thresh = 0.5;
	float error;

	for (int row = 0; row < height; row++) {
		int os = row * width;
		//first we move to the right
		for (int right = 0; right < width; right++) {
			int x = right;
			//calculate error
			if (fdata[os + x] < thresh) {
				error = fdata[os + x];
				fdata[os + x] = 0.0;
			}
			else {
				error = fdata[os + x] - 1.0;
				fdata[os + x] = 1.0;
			}
			//if we're at the bottom row OR rightmost pixel OR leftmost pixel
			if (row == height - 1 || right == width - 1 || right == 0) {
				//this will fire if we're at the bottom row
				if (row == height - 1) {
					//bottom row AND leftmost pixel
					if (right == 0) {
						fdata[os + x + 1] += (7.0 / 16.0) * error;	//same row, next col
					//	fdata[os + width + x + 1]	+= (1.0 / 16.0) * error;	//next row, next col
					//	fdata[os + width + x]		+= (5.0 / 16.0) * error;	//next row, same col
					//	fdata[os + width + x - 1]	+= (3.0 / 16.0) * error;	//next row, prev col
					}

					//bottom row AND rightmost pixel (last pixel)
					else if (right == width - 1) {
						//	fdata[os + x + 1]			+= (7.0 / 16.0) * error;	//same row, next col
						//	fdata[os + width + x + 1]	+= (1.0 / 16.0) * error;	//next row, next col
						//	fdata[os + width + x]		+= (5.0 / 16.0) * error;	//next row, same col
						//	fdata[os + width + x - 1]	+= (3.0 / 16.0) * error;	//next row, prev col
					}

					//just bottom row
					else {
						fdata[os + x + 1] += (7.0 / 16.0) * error;	//same row, next col
					//	fdata[os + width + x + 1]	+= (1.0 / 16.0) * error;	//next row, next col
					//	fdata[os + width + x]		+= (5.0 / 16.0) * error;	//next row, same col
					//	fdata[os + width + x - 1]	+= (3.0 / 16.0) * error;	//next row, prev col
					}
				}
				//NOT at the bottom row
				//leftmost pixel
				else if (right == 0) {
					fdata[os + x + 1] += (7.0 / 16.0) * error;	//same row, next col
					fdata[os + width + x + 1] += (1.0 / 16.0) * error;	//next row, next col
					fdata[os + width + x] += (5.0 / 16.0) * error;	//next row, same col
				//	fdata[os + width + x - 1]	+= (3.0 / 16.0) * error;	//next row, prev col
				}
				//rightmost pixel
				else {
					//	fdata[os + x + 1]			+= (7.0 / 16.0) * error;	//same row, next col
					//	fdata[os + width + x + 1]	+= (1.0 / 16.0) * error;	//next row, next col
					fdata[os + width + x] += (5.0 / 16.0) * error;	//next row, same col
					fdata[os + width + x - 1] += (3.0 / 16.0) * error;	//next row, prev col
				}
			}
			//not bottom row, NOR leftmost pixel, NOR rightmost pixel
			else {
				fdata[os + x + 1] += (7.0 / 16.0) * error;	//same row, next col
				fdata[os + width + x + 1] += (1.0 / 16.0) * error;	//next row, next col
				fdata[os + width + x] += (5.0 / 16.0) * error;	//next row, same col
				fdata[os + width + x - 1] += (3.0 / 16.0) * error;	//next row, prev col
			}
		}

		row++;
		os = row * width;
		if (row == height - 1) {
			break;
		}

		//next we move to the left
		for (int left = width - 1; left >= 0; left--) {
			int x = left;
			if (fdata[os + x] < thresh) {
				error = fdata[os + x];
				fdata[os + x] = 0.0;
			}
			else {
				error = fdata[os + x] - 1.0;
				fdata[os + x] = 1.0;
			}
			//if we're at the bottom row OR the leftmost pixel
			if (row == height - 1 || left == 0 || left == width - 1) {
				//if we're at the bottom row
				if (row == height - 1) {
					//AND the leftmost pixel (last pixel)
					if (left == 0) {
						//	fdata[os + x - 1]			+= (7.0 / 16.0) * error;	//same row, next col
						//	fdata[os + width + x - 1]	+= (1.0 / 16.0) * error;	//next row, next col
						//	fdata[os + width + x]		+= (5.0 / 16.0) * error;	//next row, same col
						//	fdata[os + width + x + 1]	+= (3.0 / 16.0) * error;	//next row, prev col
					}
					//bottom row AND rightmost pixel
					else if (left == width - 1) {
						fdata[os + x - 1] += (7.0 / 16.0) * error;	//same row, next col
					//	fdata[os + width + x - 1]	+= (1.0 / 16.0) * error;	//next row, next col
					//	fdata[os + width + x]		+= (5.0 / 16.0) * error;	//next row, same col
					//	fdata[os + width + x + 1]	+= (3.0 / 16.0) * error;	//next row, prev col
					}
					//otherwise we're just at the last row, not the last pixel;
					else {
						fdata[os + x - 1] += (7.0 / 16.0) * error;	//same row, next col
					//	fdata[os + width + x - 1]	+= (1.0 / 16.0) * error;	//next row, next col
					//	fdata[os + width + x]		+= (5.0 / 16.0) * error;	//next row, same col
					//	fdata[os + width + x + 1]	+= (3.0 / 16.0) * error;	//next row, prev col
					}
				}
				//NOT bottom row, but rightmost pixel
				else if (left == width - 1) {
					fdata[os + x - 1] += (7.0 / 16.0) * error;	//same row, next col
					fdata[os + width + x - 1] += (1.0 / 16.0) * error;	//next row, next col
					fdata[os + width + x] += (5.0 / 16.0) * error;	//next row, same col
				//	fdata[os + width + x + 1]	+= (3.0 / 16.0) * error;	//next row, prev col
				}
				//we're at the leftmost pixel on a row that isn't the bottom row;
				else {
					//	fdata[os + x - 1]			+= (7.0 / 16.0) * error;	//same row, next col
					//	fdata[os + width + x - 1]	+= (1.0 / 16.0) * error;	//next row, next col
					//	fdata[os + width + x]		+= (5.0 / 16.0) * error;	//next row, same col
					fdata[os + width + x + 1] += (3.0 / 16.0) * error;	//next row, prev col
				}
			}
			//neither leftmost pixel NOR bottom row;
			else {
				fdata[os + x - 1] += (7.0 / 16.0) * error;	//same row, next col
				fdata[os + width + x - 1] += (1.0 / 16.0) * error;	//next row, next col
				fdata[os + width + x] += (5.0 / 16.0) * error;	//next row, same col
				fdata[os + width + x + 1] += (3.0 / 16.0) * error;	//next row, prev col
			}
		}
	}

	for (int i = 0; i < height; i++) {
		int os = i * width * 4;
		for (int j = 0; j < width; j++) {
			int x = j * 4;
			for (int k = 0; k < 3; k++) {
				data[os + x + k] = (int)(fdata[i * width + j] * 255.0);
			}
			data[os + x + 3] = 255;
		}
	}

	delete[] fdata;
	return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
	int nwidth = width / 2;
	int nheight = height / 2;
	unsigned char * ndata = new unsigned char[nwidth * nheight * 4];

	int filter[3][3] = {
		{ 1, 2, 1},
		{ 2, 4, 2},
		{ 1, 2, 1}
	};

	for (int i = 0; i < nheight; i++) {
		int os = i * nwidth * 4;
		for (int j = 0; j < nwidth; j++) {
			int x = j * 4;
			for (int k = 0; k < 4; k++) {
				int row, col;
				int pixel = 0;
				for (int l = -1; l < 2; l++) {
					for (int m = -1; m < 2; m++) {
						row = (i * 2) + l;
						col = (j * 2) + m;
						if (row < 0) {
							row = -row;
						}
						else if (row >= height) {
							row = height - (row - height) - 1;
						}
						if (col < 0) {
							col = -col;
						}
						else if (col >= width) {
							col = width - (col - width) - 1;
						}
						row = row * width * 4;
						col = col * 4;
						pixel += data[row + col + k] * filter[l + 1][m + 1];
					}
				}
				pixel /= 16;
				ndata[os + x + k] = pixel;
			}
		}
	}

	delete[] data;
	data = ndata;
	width = nwidth;
	height = nheight;
	return true;
}// Half_Size

 
 
 
 
 
 
 
 
/////////////////////// NOT IMPLEMENTING //////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
	ClearToBlack();
	return false;
}



///////////////////////// HELPER FUNCTIONS //////////////////////////////////////

// Computes n choose s, efficiently
double Binomial(int n, int s)
{
	double        res;

	res = 1;
	for (int i = 1; i <= s; i++)
		res = (n - i + 1) * res / i;

	return res;
}// Binomial

///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
	unsigned char   *rgb = new unsigned char[width * height * 3];
	int		    i, j;

	if (!data)
		return NULL;

	// Divide out the alpha
	for (i = 0; i < height; i++)
	{
		int in_offset = i * width * 4;
		int out_offset = i * width * 3;

		for (j = 0; j < width; j++)
		{
			RGBA_To_RGB(data + (in_offset + j * 4), rgb + (out_offset + j * 3));
		}
	}

	return rgb;
}// To_RGB

///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
	TargaImage	*out_image = Reverse_Rows();

	if (!out_image)
		return false;

	if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
	{
		cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
		return false;
	}

	delete out_image;

	return true;
}// Save_Image

///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
	unsigned char   *temp_data;
	TargaImage	    *temp_image;
	TargaImage	    *result;
	int		        width, height;

	if (!filename)
	{
		cout << "No filename given." << endl;
		return NULL;
	}// if

	temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
	if (!temp_data)
	{
		cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
		width = height = 0;
		return NULL;
	}
	temp_image = new TargaImage(width, height, temp_data);
	free(temp_data);

	result = temp_image->Reverse_Rows();

	delete temp_image;

	return result;
}// Load_Image

///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
	if (!pImage)
		return false;

	if (width != pImage->width || height != pImage->height)
	{
		cout << "Difference: Images not the same size\n";
		return false;
	}// if

	for (int i = 0; i < width * height * 4; i += 4)
	{
		unsigned char        rgb1[3];
		unsigned char        rgb2[3];

		RGBA_To_RGB(data + i, rgb1);
		RGBA_To_RGB(pImage->data + i, rgb2);

		data[i] = abs(rgb1[0] - rgb2[0]);
		data[i + 1] = abs(rgb1[1] - rgb2[1]);
		data[i + 2] = abs(rgb1[2] - rgb2[2]);
		data[i + 3] = 255;
	}

	return true;
}// Difference

//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB

///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows

///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack

///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

bool TargaImage::is_greater(const quant_pop_color& lhs, const quant_pop_color& rhs) {
	if (lhs.count > rhs.count)
		return true;
	return false;
}

//this will give me the bartlett coefficient for the value at (x, y).
double TargaImage::F(double x, double y) {
	double w = 4.0;
	return ((4.0 / (w * w)) * H(floor(x) - x) * H(floor(y) - y));
}

double TargaImage::H(double s) {
	double w = 4.0;
	return (2.0 / w) * (1 - ((2 * abs(s)) / w));
}

int TargaImage::find_binomial(int N, int*& arr) {
	if (N < 3) {
		return 0;
	}

	if (arr) {
		delete[] arr;
	}
	int curr = 3;
	arr = new int[curr];
	arr[0] = 1;
	arr[1] = 2;
	arr[2] = 1;
	return find_binomial_recursive(N, curr, arr);
}

int TargaImage::find_binomial_recursive(int N, int curr, int*& arr) {
	if (N == curr) {
		return N;
	}
	int next_curr = curr + 1;
	int * temp = new int[next_curr];
	temp[0] = temp[curr] = 1;
	for (int i = 1; i < curr; i++) {
		temp[i] = arr[i - 1] + arr[i];
	}
	delete[] arr;
	arr = temp;
	return find_binomial_recursive(N, next_curr, arr);
}
