#include "ConvertingGrayScaleToBinary.h"

ConvertingGrayScaleToBinary::ConvertingGrayScaleToBinary(void)
{
}

ConvertingGrayScaleToBinary::~ConvertingGrayScaleToBinary(void)
{
}



int ConvertingGrayScaleToBinary::ConvertGraytoBinary(OutputImageType3d::Pointer image)
{
  /*if( argc < 7 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImageFile outputImageFile ";
    std::cerr << " lowerThreshold upperThreshold ";
    std::cerr << " outsideValue insideValue   "  << std::endl;
    return EXIT_FAILURE;
    }*/

  
  
  


  ReaderType::Pointer reader = ReaderType::New();
  FilterType::Pointer filter = FilterType::New();
	//imageReader->SetFileName("g:\\rawreader\\output_en.mhd");

	

  WriterType_Binary::Pointer writer = WriterType_Binary::New();
  writer->SetInput( filter->GetOutput() );
 // reader->SetFileName( "g:\\rawreader\\output_en.mhd" );

	//try
	//{ 
	//imageReader->Update();
	//}
	//catch (itk::ExceptionObject &ex)
	//{ 
	//std::cout << ex << std::endl;
	//return EXIT_FAILURE;
	//}
 
  filter->SetInput( image );
 

  const OutputPixelType outsideValue = atoi( "0" );
  const OutputPixelType insideValue  = atoi( "255" );


  filter->SetOutsideValue( outsideValue );
  filter->SetInsideValue(  insideValue  );
  

  const InputPixelType lowerThreshold = atoi( "80" );
  const InputPixelType upperThreshold = atoi( "900" );


  filter->SetLowerThreshold( lowerThreshold );
  filter->SetUpperThreshold( upperThreshold );
 
  filter->Update();
  
 

  writer->SetFileName( "g:\\CTImage\\BinaryVesselOutput.mhd" );
   try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
  }
 

  return EXIT_SUCCESS;
}
