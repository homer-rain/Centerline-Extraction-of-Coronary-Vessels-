#pragma once

//#include "itkImage.h"
//#include "itkImageFileReader.h"
 #include "itkHoughTransform2DCirclesImageFilter.h"
//// Software Guide : EndCodeSnippet
//
//#include "itkExtractImageFilter.h"
//#include "itkImageFileWriter.h"
//#include "itkImageRegionIterator.h"
#include "itkThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include <list>
#include "itkCastImageFilter.h"
#include "vnl/vnl_math.h"
//#include "QuickView.h"
#include "RegionGrowing.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "itkAnalyzeImageIO.h"
#include "itkOrientImageFilter.h"
#include<sys/stat.h>
//#include <vtkSmartPointer.h>
//#include <vtkImageBlend.h>
#include "itkOrImageFilter.h"
#include "RayCAsting.h" 

#include "itkThresholdSegmentationLevelSetImageFilter.h"
// Software Guide : EndCodeSnippet
#include "itkFastMarchingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "FastMarchingLevelSetBase.h"

#include "itkZeroCrossingImageFilter.h"
#include "vtkRaycaster.h"
#include "fastmarching.h"
typedef float        InputPixelType;
typedef float       OutputPixelType;


//typedef itk::Image< InputPixelType,  3 >    InputImageType;
typedef itk::Image< OutputPixelType, 2 >    OutputImageType;

typedef itk::ImageFileReader< InputImageType  >  ReaderType;
typedef itk::ImageFileWriter< OutputImageType >  WriterType;
typedef itk::ImageFileWriter< InputImageType >  WriterType3D;
const     unsigned int    Dimension = 2;

typedef   float           AccumulatorPixelType;
//	typedef itk::Image<unsigned int, 2>  UnsignedintImageType;
typedef itk::Image<float, 2>  FloatImageType;

typedef itk::Image< AccumulatorPixelType, Dimension > AccumulatorImageType;



class AortaDetection
{
	

public:
//struct CircleIndex
//{
//int x;
//int y;
//}m_CircleIndex[500];
	InputImageType::IndexType *m_CircleIndex ;
	AortaDetection(InputImageType::Pointer image);
	~AortaDetection(void);
	int indexCounter1;


	int ExtractSlice(  InputImageType::IndexType *m_CircleIndex );
	OutputImageType::Pointer  houghCircleTransform(OutputImageType::Pointer image,OutputImageType::Pointer Vesselimage,int i);
	InputImageType::IndexType  FindLastPoint();
	int AortaLevelsetSegmentation(InputImageType::Pointer image,int x,int y);
	OutputImageType::IndexType localIndex;
	OutputImageType::IndexType prevLocalIndex;
	OutputImageType::IndexType radiusIndex;
	OutputImageType::IndexType lastRadiusIndex;

	double sumOfIntensity;
	double lastSumOfIntensity;
	double lastRadius;
	BOOL   FINISHED1;
	BOOL   FINISHED2;
	RayCAsting* m_RayCasting;
	VTKRayCaster* m_VTKRayCaster;
	 RegionGrowing* m_RegionGrawing;
	 InputImageType::Pointer image2;
	 InputImageType::IndexType center;
	 InputImageType::IndexType coronaryOstium;
	 InputImageType::IndexType* Centerline_LCX;
	 InputImageType::IndexType* Centerline_RCX;
	  InputImageType::IndexType* Centerline_LCA;
	 int WriteEachVessel(InputImageType::IndexType* LCX,char* fileName);
	 InputImageType::Pointer VesselOutput;
	 FastMarching* m_FastMarching;
	 InputImageType::IndexType* queue_CenterlinePoint;


	

};
