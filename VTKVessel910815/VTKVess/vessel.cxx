/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiScaleHessianBasedMeasureImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007/07/27 21:19:46 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "AortaDetection.h"
#include "ConvertingGrayScaleToBinary.h"
#include "Resampler.h"
#include "CoronaryTracking.h"
#include "VTKRayCaster.h"
#include "FastMarching.h"
#include "Evaluation.h"
struct CircleIndex
{
int x;
int y;
};

	const unsigned char Dim = 3;

	typedef float PixelType;

	// Declare the types of the images
	typedef itk::Image<PixelType,Dim> ImageType;

	typedef itk::ImageFileReader<ImageType> FileReaderType;
	typedef itk::ImageFileWriter<ImageType> FileWriterType;

	typedef itk::RescaleIntensityImageFilter<ImageType> RescaleFilterType;

	// Declare the type of enhancement filter
	typedef itk::HessianToObjectnessMeasureImageFilter<double,Dim> ObjectnessFilterType;

	// Declare the type of multiscale enhancement filter
	typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType,ObjectnessFilterType> 
	MultiScaleEnhancementFilterType;


void CreateImage(ImageType::Pointer image)
{
  // Create an image with 2 connected components
  ImageType::RegionType region;
  ImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
 
  ImageType::SizeType size;
  unsigned int NumRows = 256;
  unsigned int NumCols = 256;
  size[0] = NumRows;
  size[1] = NumCols;
  size[2] = 136;
 
  region.SetSize(size);
  region.SetIndex(start);
 
  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer( itk::NumericTraits<ImageType::PixelType>::Zero);
 
  // Make a rectangle, centered at (100,150) with sides 160 & 240
  // This provides a 20 x 30 border around the square for the crop filter to remove
  for(unsigned int r = 20; r < 180; r++)
    {
    for(unsigned int c = 30; c < 270; c++)
      {
      ImageType::IndexType pixelIndex;
      pixelIndex[0] = r;
      pixelIndex[1] = c;
 
      image->SetPixel(pixelIndex, 200);
      }
    }
}
int main( int argc, char *argv[] )
{
  /*if ( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Input_Image"
              << " Enhanced_Output_Image" 
              << " Scales_Output_Image"
              << std::endl;
    return EXIT_FAILURE;
    }*/

  // Define the dimension of the images
 
	//RayCAsting* m_RayCaster = new RayCAsting();
	//m_RayCaster->ChangeGrayToRGB();

char* filePath = "F:\\CTImage\\training\\dataset02\\vesselness.mhd";
	
  typedef itk::ImageRegionIterator< ImageType> ImageRegionIterator;
 // AortaDetection::CircleIndex* m_CircleIndex = new AortaDetection::CircleIndex[500];
  ImageType::IndexType *m_CircleIndex = new ImageType::IndexType[500];

  // Read the image
  FileReaderType::Pointer imageReader = FileReaderType::New();
  //char* fileName = filePath+"";
  /*Resampler* m_Resampler = new Resampler();
  m_Resampler->itkResampleImageTest();*/
  //imageReader->SetFileName("g:\\rawreader\\downsampledImage1.mhd");
 //  imageReader->SetFileName("g:\\CTImage\\downsampledimage1.mhd");
 //   imageReader->SetFileName("g:\\CTImage\\vesselness.mhd");
  //imageReader->SetFileName("g:\\Rawreader\\output_en2.mhd");
  //imageReader->SetFileName(filePath+" vesselness.mhd");
  imageReader->SetFileName(filePath);
 // imageReader->SetFileName("g:\\CTImage\\FastMarchingFilterOutput2.mhd");
 // imageReader->SetFileName("F:\\CTImage\\training\\dataset00\\output_en1.mhd");
  try
  { 
    imageReader->Update();
  }
  catch (itk::ExceptionObject &ex)
  { 
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

 // m_VTKRayCaster->fixedpointRaycasting();
 // m_VTKRayCaster->vtkrenderer(imageReader->GetOutput());
  //m_VTKRayCaster->SphereCreating(imageReader->GetOutput());
  double spacing[3];
 /* spacing[0] = 0.015;
  spacing[1] = 0.0259;
  spacing[2] = 0.0019;*/
  spacing[0] = 0.03;
  spacing[1] = 0.03;
  spacing[2] = 0.03;
 /* Evaluation* m_Evaluation = new Evaluation();
  m_Evaluation->ImageUpsample(imageReader->GetOutput(),spacing);*/
  AortaDetection* m_AortaDetector = new AortaDetection(imageReader->GetOutput());

  int indexcounter = m_AortaDetector->ExtractSlice(m_CircleIndex);
  

  imageReader->SetFileName("G:\\RawReader\\image00_2.mhd");

	try
	{ 
	imageReader->Update();
	}
	catch (itk::ExceptionObject &ex)
	{ 
	std::cout << ex << std::endl;
	return EXIT_FAILURE;
	}

	CoronaryTracking* m_CoronaryTracking = new CoronaryTracking();
	m_CoronaryTracking->track(m_CircleIndex,imageReader->GetOutput(),m_AortaDetector->indexCounter1);

 //double maxValue1 = imageReader->GetOutput()->GetPixel(m_CircleIndex[0]);
 //double maxValue2 = imageReader->GetOutput()->GetPixel(m_CircleIndex[0]);
 //ImageType::IndexType indexOfMaxValue1;
 //ImageType::IndexType indexOfMaxValue2;
 //int index = 0;
 //for(int i=1;i<m_AortaDetector->indexCounter1;i++)
 // {
	//  if(maxValue1<imageReader->GetOutput()->GetPixel(m_CircleIndex[i]))
	//  {
	//	  maxValue1 = imageReader->GetOutput()->GetPixel(m_CircleIndex[i]);
	//	  indexOfMaxValue1 = m_CircleIndex[i];
	//	  index =i;

	//  }

	// 
 // }
 //std::cerr <<imageReader->GetOutput()->GetPixel(indexOfMaxValue1) << std::endl;;

 //for(int i=1;i<m_AortaDetector->indexCounter1;i++)
 // {
	//  if(maxValue2<imageReader->GetOutput()->GetPixel(m_CircleIndex[i])&&i!=index)
	//	 // if((m_CircleIndex[i].m_Index[0]!=indexOfMaxValue1.m_Index[0])&&(m_CircleIndex[i].m_Index[1]!=indexOfMaxValue1.m_Index[1])&&(m_CircleIndex[i].m_Index[2]!=indexOfMaxValue1.m_Index[2]))
	//  {
	//	  maxValue2 = imageReader->GetOutput()->GetPixel(m_CircleIndex[i]);
	//	  indexOfMaxValue2 = m_CircleIndex[i];

	//  }

	//  std::cout <<imageReader->GetOutput()->GetPixel(m_CircleIndex[i]) << std::endl;;
 // }

 //std::cerr <<imageReader->GetOutput()->GetPixel(indexOfMaxValue2) << std::endl;;

 
	//ImageType::RegionType insideRegion;
	//ImageType::SizeType insideSize = {{ 256, 256, 1 }};
	//ImageType::IndexType insideIndex = m_CircleIndex[0];/*{{ 10, 10, 10 }};*/
	//insideRegion.SetSize( insideSize );
	//insideRegion.SetIndex( insideIndex );

	//
	//ImageRegionIterator it( imageReader->GetOutput(), insideRegion );
	//it.GoToBegin();
	//while( !it.IsAtEnd() )
	//{
	//it.Set( itk::NumericTraits< PixelType >::max() );
	//++it;
	//std::cerr <<it.GetIndex() << std::endl;;
	//}
//    std::cerr <<it.GetIndex() << std::endl;;
	
	

	ConvertingGrayScaleToBinary* m_GrayscaleToBinaryConvertor = new ConvertingGrayScaleToBinary();
	m_GrayscaleToBinaryConvertor->ConvertGraytoBinary(imageReader->GetOutput());;

	imageReader->GetOutput()->GetSpacing();
	imageReader->GetOutput()->GetOrigin();


  // Instantiate the multiscale filter and set the input image
  MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter = 
  MultiScaleEnhancementFilterType::New();
  multiScaleEnhancementFilter->SetInput(imageReader->GetOutput());
  multiScaleEnhancementFilter->SetSigmaMin(0.5);
  multiScaleEnhancementFilter->SetSigmaMax(4.0);
  multiScaleEnhancementFilter->SetNumberOfSigmaSteps(10);

  // Get the objectness filter and set the parameters
  ObjectnessFilterType* objectnessFilter = 
  	multiScaleEnhancementFilter->GetHessianToMeasureFilter();
  objectnessFilter->SetScaleObjectnessMeasure(false);
  objectnessFilter->SetBrightObject(true);
  objectnessFilter->SetAlpha(0.5);
  objectnessFilter->SetBeta(0.5);
  objectnessFilter->SetGamma(5.0);
  objectnessFilter->SetObjectDimension(1);

  

  // The above is equivalent to vesselness

  // Now run the multiscale filter
  try
  {
    multiScaleEnhancementFilter->Update();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
  }
 
  // Write the enhanced image 
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName("G:\\CTImage\\output_en.mhd");
  writer->SetInput(multiScaleEnhancementFilter->GetOutput());



  /* ImageRegionIterator it(multiScaleEnhancementFilter->GetOutput(),
    multiScaleEnhancementFilter->GetOutput()->GetLargestPossibleRegion());*/

	
	
	
	try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
  }
  
  // Write the image containing the best response scales
  FileWriterType::Pointer writer2 = FileWriterType::New();
  writer2->SetFileName("G:\\output_vessel\\image02.mhd");
  writer2->SetInput(multiScaleEnhancementFilter->GetScalesOutput());
  try
  {
    writer2->Update();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
  }
}

