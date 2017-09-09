#pragma once
#include "itkBinaryThresholdImageFilter.h"
// Software Guide : EndCodeSnippet

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

	typedef  float  InputPixelType;
	typedef  float  OutputPixelType;



	typedef itk::Image< InputPixelType,  3 >   InputImageType;
	typedef itk::Image< OutputPixelType, 3 >   OutputImageType3d;



	typedef itk::BinaryThresholdImageFilter<
	InputImageType, OutputImageType3d >  FilterType;


	typedef itk::ImageFileReader< InputImageType >  ReaderType;



	typedef itk::ImageFileWriter< OutputImageType3d >  WriterType_Binary;

class ConvertingGrayScaleToBinary
{
public:
	ConvertingGrayScaleToBinary(void);
	~ConvertingGrayScaleToBinary(void);
	int ConvertGraytoBinary(OutputImageType3d::Pointer image);
};
