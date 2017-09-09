//  Software Guide : BeginLatex
//
//  This example illustrates the use of the \doxygen{Similarity2DTransform}. A
//  similarity transform involves rotation, translation and scaling. Since the
//  parameterization of rotations is difficult to get in a generic $ND$ case, a
//  particular implementation is available for $2D$.
//
//
//  Software Guide : EndLatex 



#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"



// Software Guide : BeginLatex
//
// The most important headers to include here are the ones corresponding to the
// resampling image filter, the transform, the interpolator and the smoothing
// filter.
//
// Software Guide : EndLatex 


// Software Guide : BeginCodeSnippet
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"
// Software Guide : EndCodeSnippet



#include "itkCastImageFilter.h"

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

//  Software Guide : BeginLatex
//
//  The header file of the transform is included below.
//
//  \index{itk::Similarity2DTransform!header}
//
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include "itkSimilarity2DTransform.h"
// Software Guide : EndCodeSnippet
#pragma once

class Resampler
{
	
public:
	Resampler(void);
	~Resampler(void);
	//const unsigned char Dim ;
	typedef unsigned short PixelType;
	// Declare the types of the images
	typedef itk::Image<PixelType,3> ImageType;
	 
	 int itkResampleImageTest( );
};

