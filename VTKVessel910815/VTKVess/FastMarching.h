#pragma once

// Software Guide : BeginCodeSnippet
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
// Software Guide : EndCodeSnippet

//  Software Guide : BeginLatex
//  
//  The headers of the GradientMagnitudeRecursiveGaussianImageFilter and
//  SigmoidImageFilter are included below. Together, these two filters will
//  produce the image potential for regulating the speed term in the
//  differential equation describing the evolution of the level set.
//
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
// Software Guide : EndCodeSnippet


//  Software Guide : BeginLatex
//  
//  Of course, we will need the \doxygen{Image} class and the
//  FastMarchingImageFilter class. Hence we include their headers.
//
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include "itkImage.h"
#include "itkFastMarchingImageFilter.h"
// Software Guide : EndCodeSnippet


//  Software Guide : BeginLatex
//  
//  The time-crossing map resulting from the FastMarchingImageFilter
//  will be thresholded using the BinaryThresholdImageFilter. We
//  include its header here.
//
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include "itkBinaryThresholdImageFilter.h"
// Software Guide : EndCodeSnippet


//  Software Guide : BeginLatex
//  
//  Reading and writing images will be done with the \doxygen{ImageFileReader}
//  and \doxygen{ImageFileWriter}.
//
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
// Software Guide : EndCodeSnippet


//  The \doxygen{RescaleIntensityImageFilter} is used to renormailize the
//  output of filters before sending them to files.
// 
#include "itkRescaleIntensityImageFilter.h"


class FastMarching
{
	typedef   float           InternalPixelType;
  
  typedef itk::Image< InternalPixelType, 3 >  InternalImageType;
  
  typedef float                           OutputPixelType;
  typedef itk::Image< OutputPixelType, 3 > OutputImageType;
  
  typedef itk::BinaryThresholdImageFilter< InternalImageType, 
                        OutputImageType    >    ThresholdingFilterType;
public:
	FastMarching(void);
	~FastMarching(void);

	int  Segment(InternalImageType::Pointer image, InternalImageType::IndexType  seedPosition  );
};
