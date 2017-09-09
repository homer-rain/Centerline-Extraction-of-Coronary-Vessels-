#pragma once
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkCastImageFilter.h"
// Software Guide : EndCodeSnippet
class Evaluation
{
	

  typedef  float   InputPixelType;

  typedef   float           InternalPixelType;
  typedef   float   OutputPixelType;

  typedef itk::Image< InputPixelType,    3 >   InputImageType;
  typedef itk::Image< InternalPixelType, 3 >   InternalImageType;
  typedef itk::Image< OutputPixelType,   3 >   OutputImageType;

public:
	Evaluation(void);
	~Evaluation(void);
	int ImageUpsample(InputImageType::Pointer image, double spacing[3]);
};
