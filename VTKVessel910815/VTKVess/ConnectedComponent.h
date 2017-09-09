#pragma once

#include "itkBinomialBlurImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkScalarConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkFilterWatcher.h"
#include "itkChangeInformationImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
typedef   unsigned short           IntPixelType;

  typedef itk::Image< IntPixelType, 3 >  IntImageType;
  // Software Guide : EndCodeSnippet

	typedef   float  InternalPixelType;
	typedef   bool  MaskPixelType;


	typedef itk::Image< IntPixelType, 3 >  InternalImageType1;
	typedef itk::RGBPixel<unsigned char>   RGBPixelType;
	typedef itk::Image<RGBPixelType, 3>    RGBImageType;

class ConnectedComponent
{

public:
	ConnectedComponent(void);
	~ConnectedComponent(void);
	RGBImageType::Pointer filterImage(InternalImageType1::Pointer image);
	RGBImageType::Pointer  itkRelabelComponentImageFilterTest(InternalImageType1::Pointer image );
};
