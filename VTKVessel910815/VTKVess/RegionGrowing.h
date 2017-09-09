#pragma once
#include "itkNeighborhoodConnectedImageFilter.h"
// Software Guide : EndCodeSnippet


#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
// Software Guide : EndCodeSnippet


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"

#include "itkBinaryBallStructuringElement.h" 
#include "itkPasteImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
typedef float                            OutputPixelType;
typedef itk::Image< OutputPixelType, 2 > OutputImageType;
typedef   float           InternalPixelType;
typedef itk::Image<OutputPixelType,3> ImageType;
typedef itk::Image< InternalPixelType, 2 >  InternalImageType;
typedef itk::CastImageFilter< InternalImageType, OutputImageType >    CastingFilterType;


typedef itk::BinaryBallStructuringElement<  OutputPixelType,  3  >             StructuringElementType;
typedef itk::GrayscaleMorphologicalOpeningImageFilter< ImageType,  ImageType, StructuringElementType >  ErodeFilterType;
typedef itk::PasteImageFilter <ImageType, ImageType >
    PasteImageFilterType;
class RegionGrowing
{

public:
	RegionGrowing(ImageType::Pointer image);
	~RegionGrowing(void);
	OutputImageType::Pointer AortaSegmentation(OutputImageType::Pointer image,OutputImageType::IndexType radiusIndex,int axialnum);
	ImageType::Pointer image;
	PasteImageFilterType::Pointer pasteFilter;
	ImageType::Pointer Opening();
	void CreateImage();
	ImageType::Pointer VesselSegmentation(ImageType::Pointer image3d,ImageType::IndexType radiusIndex,int axialnum);
};
