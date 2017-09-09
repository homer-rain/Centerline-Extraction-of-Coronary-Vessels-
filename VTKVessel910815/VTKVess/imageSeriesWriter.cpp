#include "imageSeriesWriter.h"

imageSeriesWriter::imageSeriesWriter(void)
{
}

imageSeriesWriter::~imageSeriesWriter(void)
{
}
BOOL imageSeriesWriter::write()
{
	 typedef itk::Image< unsigned short, 3 >      ImageType;
  typedef itk::ImageFileReader< ImageType >   ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( "" );
  typedef itk::Image< unsigned short, 2 >     Image2DType;

  typedef itk::ImageSeriesWriter< ImageType, Image2DType > WriterType;

  WriterType::Pointer writer = WriterType::New();

  writer->SetInput( reader->GetOutput() );
  typedef itk::NumericSeriesFileNames    NameGeneratorType;

  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();

   std::string format = "outputaorta";
  format += "%03d.";
  format += "mhd";   // filename extension

  nameGenerator->SetSeriesFormat( format.c_str() );

   try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading the image" << std::endl;
    std::cerr << excp << std::endl;
    }
   ImageType::ConstPointer inputImage = reader->GetOutput();
  ImageType::RegionType   region     = inputImage->GetLargestPossibleRegion();
  ImageType::IndexType    start      = region.GetIndex(); 
  ImageType::SizeType     size       = region.GetSize(); 
   const unsigned int firstSlice = start[2];
  const unsigned int lastSlice  = start[2] + size[2] - 1;

  nameGenerator->SetStartIndex( firstSlice );
  nameGenerator->SetEndIndex( lastSlice );
  nameGenerator->SetIncrementIndex( 1 );

  writer->SetFileNames( nameGenerator->GetFileNames() );

   try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading the image" << std::endl;
    std::cerr << excp << std::endl;
    }

   return EXIT_SUCCESS;
}