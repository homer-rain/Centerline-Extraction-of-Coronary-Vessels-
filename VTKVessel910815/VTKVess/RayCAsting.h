#pragma once
#include "itkRGBToLuminanceImageFilter.h"
// Software Guide : EndCodeSnippet


#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNormalizeImageFilter.h"
typedef   float    InputPixelType;
typedef   float    OutputPixelType;
typedef itk::Image< InputPixelType,  3 >   InputImageType;
//typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

#include "cmath"
struct point3D
{
  double _x, _y, _z;
  point3D() : _x(0.0), _y(0.0), _z(0.0) {}
  point3D(double x, double y, double z) : _x(x), _y(y), _z(z) {}
  point3D(const point3D &p) : _x(p._x), _y(p._y), _z(p._z) {}
  point3D &operator=(const point3D &p) { 
    if (this != &p) _x=p._x, _y=p._y, _z=p._z;
	return *this;
  }
  double operator[](int i) const { 
    if (i==0) return _x;
    else if (i==1) return _y; 
    else if (i==2) return _z; 
    else return 0;
  }
  point3D operator+(const point3D &p) const { return point3D(_x+p._x, _y+p._y, _z+p._z); }
  point3D operator-(const point3D &p) const { return point3D(_x-p._x, _y-p._y, _z-p._z); }
  point3D operator*(double f) const { return point3D(_x*f, _y*f, _z*f); }
  point3D operator/(double f) const { return point3D(_x/f, _y/f, _z/f); }
  point3D &operator+=(const point3D &p) { _x+=p._x, _y+=p._y, _z+=p._z; return *this; }
  double dot(const point3D &p) const { return _x*p._x+_y*p._y+_z*p._z; }
  double length() const { return std::sqrt(_x*_x+_y*_y+_z*_z); }
  void normalize() {
	  double f = 1.0/length();
	  _x *= f;
	  _y *= f;
	  _z *= f;
  }
//  point3D operator*(double f, const point3D &p)
//{
//	return p*f;
//}
};



class RayCAsting
{

	
struct light {
	point3D pos;
	float color;
};

struct ray {
	point3D start;
	point3D dir;
};
struct material {
	float reflection;
	float color;
};
public:
	RayCAsting(InputImageType::Pointer image);
	~RayCAsting(void);
	InputImageType::Pointer ChangeGrayToRGB(InputImageType::Pointer image);
	void SphereRayCasting(InputImageType::IndexType center,InputImageType::IndexType VoxelIndex,float radius);
	float dotProduce(float* A,float* B);
	float* normalize(float* A);
	void draw(InputImageType::IndexType center,InputImageType::IndexType VoxelIndex,float radius);

private:
	float Ka;
	float Kd;
	float Ks;
	InputImageType::Pointer image;
	InputImageType::Pointer rescaledImage;
	float color;
	float lambert;
};
