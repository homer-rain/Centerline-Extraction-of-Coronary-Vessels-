#pragma once
#include "itkImage.h"
// Software Guide : BeginCodeSnippet
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
// Software Guide : EndCodeSnippet
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

typedef float        InputPixelType;
typedef float       OutputPixelType;


typedef itk::Image< InputPixelType,  3 >    InputImageType;


class CoronaryTracking
{
public:
	CoronaryTracking(void);
	~CoronaryTracking(void);
	void track(InputImageType::IndexType *m_CircleIndex,InputImageType::Pointer image,int indexCounter);
};
