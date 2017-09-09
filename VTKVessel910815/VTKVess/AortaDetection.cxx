#include "AortaDetection.h"


	int index1 = 0;
	int index2 = 0;
	float max_Vesselness_LCX[100];// = 1320;
	double angle_max_Vesselness_LCX[100];
int maximumValue(float *array)
{
     int length = index1;  // establish size of array
     double max = array[0];       // start with max = first element
	int ind = 0;

     for(int i = 1; i<length; i++)
     {
          if(array[i] > max)
		  {
                max = array[i];
				ind = i;}
     }
     return ind;                // return highest value in array
}

AortaDetection::AortaDetection(InputImageType::Pointer imageInit)
{
	FINISHED1 = false;
	FINISHED2 = false;
	lastRadius = 0;
	m_VTKRayCaster = new VTKRayCaster(imageInit,imageInit);
	Centerline_LCX = new InputImageType::IndexType[5000];
	Centerline_RCX = new InputImageType::IndexType[5000];
	Centerline_LCA = new InputImageType::IndexType[5000];
	//VesselOutput = imageInit;
	m_FastMarching = new FastMarching();
	queue_CenterlinePoint = new InputImageType::IndexType[500000];
}

AortaDetection::~AortaDetection(void)
{
}



 int indexCounter = 0;
 int circleIndex = 0;
OutputImageType::Pointer AortaDetection::  houghCircleTransform(OutputImageType::Pointer image,OutputImageType::Pointer Vesselimage, int i)
{



 int num = 0;
 std::cout << "Computing Hough Map" << std::endl;

  typedef itk::HoughTransform2DCirclesImageFilter<OutputPixelType,AccumulatorPixelType> HoughTransformFilterType;
  HoughTransformFilterType::Pointer houghFilter = HoughTransformFilterType::New();

  houghFilter->SetInput( image );

  houghFilter->SetNumberOfCircles( atoi("2") );
  houghFilter->SetMinimumRadius(   atof("20") );
  houghFilter->SetMaximumRadius(   atof("30") );
  double radius = 0;

  /*if( argc > 6 )
    {
    houghFilter->SetSweepAngle( atof(argv[6]) );
    }
  if( argc > 7 )
    {
    houghFilter->SetSigmaGradient( atoi(argv[7]) );
    }
  if( argc > 8 )
    {
    houghFilter->SetVariance( atof(argv[8]) );
    }
  if( argc > 9 )
    {
    houghFilter->SetDiscRadiusRatio( atof(argv[9]) );
    }*/

  houghFilter->Update();
  AccumulatorImageType::Pointer localAccumulator = houghFilter->GetOutput();




float max_Vesselness_RCX [100];
double angle_max_Vesselness_RCX[100];
float LCX_dist = 0;
float RCX_dist = 0;
for(int i=0;i<100;i++)
	max_Vesselness_LCX[i] = 0;
for(int i=0;i<100;i++)
	angle_max_Vesselness_LCX[i] = 0;
for(int i=0;i<100;i++)
	max_Vesselness_RCX[i] = 0;
for(int i=0;i<100;i++)
	angle_max_Vesselness_RCX[i] = 0;

 HoughTransformFilterType::CirclesListType circles;
  circles = houghFilter->GetCircles( atoi("2") );
  std::cout << "Found " << circles.size() << " circle(s)." << std::endl;


  typedef  float                           OutputPixelType;
  typedef  itk::Image< OutputPixelType, Dimension > OutputImageType;


  OutputImageType::Pointer  localOutputImage = OutputImageType::New();

  OutputImageType::RegionType region;
  region.SetSize(image->GetLargestPossibleRegion().GetSize());
  region.SetIndex(image->GetLargestPossibleRegion().GetIndex());
 
  localOutputImage->SetRegions( region );
  localOutputImage->SetOrigin(image->GetOrigin());
  localOutputImage->SetSpacing(image->GetSpacing());
  localOutputImage->Allocate();
  localOutputImage->FillBuffer(0);
  InputImageType::IndexType EndPoint_LCX;
  OutputImageType::IndexType EndPoint_RCX;
  float distThershold = 10;


   typedef HoughTransformFilterType::CirclesListType CirclesListType;
  CirclesListType::const_iterator itCircles = circles.begin();
  circleIndex =0;
  if(i>130)
  {
  lastRadiusIndex[0] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
  lastRadiusIndex[1] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
  }
  while( itCircles != circles.end() )
    {
    std::cout << "Center: ";
    std::cout << (*itCircles)->GetObjectToParentTransform()->GetOffset()
              << std::endl;
    std::cout << "Radius: " << (*itCircles)->GetRadius()[0] << std::endl;
	radiusIndex[0] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
	radiusIndex[1] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
	double distance = vcl_sqrt(pow((float)( lastRadiusIndex[0] - radiusIndex[0]),2)+pow((float)( lastRadiusIndex[1] - radiusIndex[1]),2));
	/*if(distance> 10)
			FINISHED2 = true;*/
	
	//////////////////////////////////////////////////////////////////////////////////
	/*if((*itCircles)->GetRadius()[0]>radius&&radius>0)
	{
		circleIndex ++;
	}*/
	if(distance<distThershold/*circleIndex == 0*/)
	{
		circleIndex++;
		//indexCounter = 0;
		radiusIndex[0] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
		radiusIndex[1] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
	
		radius = (*itCircles)->GetRadius()[0];
		
			std::cout << "Found circle \n" ;
	// AortaLevelsetSegmentation(image,radiusIndex[0],radiusIndex[1]);
		m_RegionGrawing->AortaSegmentation(image,radiusIndex,i);
			num = 0;
			sumOfIntensity = 0;
			//std::cout<<radiusIndex[0]<<" "<<radiusIndex[1]<< " "<<"\n";
			InputImageType::IndexType NewCenter  ;
	double maxangle = 0;
	InputImageType::IndexType lastCenter;
	double lastAngle = 0;	
	index1 = 0;
	for(double angle = 0;angle <= 2*vnl_math::pi; angle += vnl_math::pi/60.0 )
      {
      localIndex[0] =
         (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0]
                                  + (*itCircles)->GetRadius()[0]*vcl_cos(angle));
      localIndex[1] =
         (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[1]
                                  + (*itCircles)->GetRadius()[0]*vcl_sin(angle));
								
								  if(((angle>0&&angle<vnl_math::pi/4)/*||(angle>3*vnl_math::pi/2&&angle<7*vnl_math::pi/4)*/)&&(i==115))
									{
										angle_max_Vesselness_LCX[index1] = angle;
										max_Vesselness_LCX[index1++] =(double) image->GetPixel(localIndex);
											Vesselimage->SetPixel(localIndex,255);
											Vesselimage->SetPixel(localIndex,255);
											LCX_dist  =LCX_dist+ vcl_sqrt(pow((float)( localIndex[0] - prevLocalIndex[0]),2)+pow((float)( localIndex[1] - prevLocalIndex[1]),2));
											std::cout<<"LCX_dist"<<" "<<LCX_dist<< " max_Vesselness"<<max_Vesselness_LCX[index1-1]<< "\n";
											EndPoint_LCX.m_Index[0] = localIndex.m_Index[0]+5;
											EndPoint_LCX.m_Index[1] = localIndex.m_Index[1]+5;
											EndPoint_LCX.m_Index[2] = i;
											//Vesselimage->SetPixel(EndPoint_LCX,255);
											prevLocalIndex = localIndex ;
											center.m_Index[0] = localIndex.m_Index[0];
											center.m_Index[1] = localIndex.m_Index[1];
											center.m_Index[2] = i;
											lastCenter = center;
											
											//m_VTKRayCaster->SetTwoDimImage(Vesselimage);
													
								  }
								
								 
								
								  if((angle>5*vnl_math::pi/4&&angle<7*vnl_math::pi/4)&&i==96)
									{
											angle_max_Vesselness_RCX[index1] = angle;
											max_Vesselness_RCX[index1++] =(double) image->GetPixel(localIndex);
											Vesselimage->SetPixel(localIndex,255);
											Vesselimage->SetPixel(localIndex,255);
											LCX_dist  =LCX_dist+ vcl_sqrt(pow((float)( localIndex[0] - prevLocalIndex[0]),2)+pow((float)( localIndex[1] - prevLocalIndex[1]),2));
											std::cout<<"LCX_dist"<<" "<<LCX_dist<< " max_Vesselness"<<max_Vesselness_LCX[index1-1]<< "\n";
											EndPoint_LCX.m_Index[0] = localIndex.m_Index[0]+5;
											EndPoint_LCX.m_Index[1] = localIndex.m_Index[1]+5;
											EndPoint_LCX.m_Index[2] = i;
											//Vesselimage->SetPixel(EndPoint_LCX,255);
											prevLocalIndex = localIndex ;
											center.m_Index[0] = localIndex.m_Index[0];
											center.m_Index[1] = localIndex.m_Index[1];
											center.m_Index[2] = i;
											lastCenter = center;
											
											/*max_Vesselness_RCX = image->GetPixel(localIndex);
											Vesselimage->SetPixel(localIndex,255);
											EndPoint_RCX.m_Index[0] = localIndex.m_Index[0]-5;
											EndPoint_RCX.m_Index[1] = localIndex.m_Index[1]-5;
											RCX_dist  = RCX_dist + vcl_sqrt(pow((float)( localIndex[0] - prevLocalIndex[0]),2)+pow((float)( localIndex[1] - prevLocalIndex[1]),2));
											std::cout<<"RCX_dist"<<" "<<RCX_dist<< " max_Vesselness"<<max_Vesselness_RCX<< "\n";
											prevLocalIndex = localIndex ;
											Vesselimage->SetPixel(EndPoint_RCX,255);
											center.m_Index[0] = localIndex.m_Index[0];
											center.m_Index[1] = localIndex.m_Index[1];
											center.m_Index[2] = i;*/
											//if(i<100&&i>90)
											//{
											//	m_VTKRayCaster->lastMaxOpacity = 0;
											//	m_VTKRayCaster->Circle_RayCasting(center,5*vnl_math::pi/4,7*vnl_math::pi/4,2,image2);
											///*	center.m_Index[0] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])
											//	+ (10*vcl_cos(maxangle));
											//	center.m_Index[1] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])
											//	+ (10*vcl_sin(maxangle));*/
											//	center.m_Index[0] = EndPoint_RCX.m_Index[0];
											//		center.m_Index[1] = EndPoint_RCX.m_Index[1];
											//	center.m_Index[2] = i;
											//	m_VTKRayCaster->Intersect(center,2,&NewCenter);
											//}
										
								  }
	
								  m_CircleIndex[indexCounter].m_Index[0] = localIndex[0];
								  m_CircleIndex[indexCounter].m_Index[1] = localIndex[1];
								  m_CircleIndex[indexCounter++].m_Index[2] = i;
								 

	
      OutputImageType::RegionType outputRegion =
                                  localOutputImage->GetLargestPossibleRegion();

      if( outputRegion.IsInside( localIndex ) )
        {
        localOutputImage->SetPixel( localIndex, 255 );
		image->SetPixel(localIndex,(0,255,0));
        }
	 
	    ////////////////////////////////////////////finding min and max of intensity////////////////////////////////////////////////
								
								  for(int r = 1;r<(*itCircles)->GetRadius()[0];r++)
								  {
									   localIndex[0] = 	 (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0]+r*vcl_cos(angle));
									   localIndex[1] =   (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[1]+ r*vcl_sin(angle));
								//  image->SetPixel(localIndex,(255,0,0));
								sumOfIntensity += image->GetPixel(localIndex);
								num++;
							//	std::cout<<localIndex[0]<<" "<<localIndex[1]<<" "<< image->GetPixel(localIndex)<<"\n";
								  }
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      }
double rad = 2;
							 if(i==115)
											{
												maxangle =angle_max_Vesselness_LCX[ maximumValue(max_Vesselness_LCX)];
												center[0] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0]+ (*itCircles)->GetRadius()[0]*vcl_cos(maxangle));
												center[1] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[1]+ (*itCircles)->GetRadius()[0]*vcl_sin(maxangle));
												center[2] = i;
											//	coronaryOstium = center;
 											for(int j=0;j<30;j++)
											{
												maxangle = m_VTKRayCaster->Circle_RayCasting(center,-1*vnl_math::pi/6,2*vnl_math::pi/3,rad,image2);
												//if(abs(((maxangle+lastAngle)-vnl_math::pi))<0.1)
												//{
												//	center.m_Index[0] = center.m_Index[0]/*(long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])*/
											 //+ (2*vcl_cos(maxangle))+1;
											 //center.m_Index[1] = center.m_Index[1]/*(long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])*/
											 //+ (2*vcl_sin(maxangle));
												//center.m_Index[2] = i;
												//}
												//else{
												center.m_Index[0] = center.m_Index[0]/*(long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])*/
											 + (2*vcl_cos(maxangle));
											 center.m_Index[1] = center.m_Index[1]/*(long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])*/
											 + (2*vcl_sin(maxangle));
												center.m_Index[2] = i;
												//}
												if(center.m_Index[0]==lastCenter.m_Index[0]&&center.m_Index[1]==lastCenter.m_Index[1])
												{
													
													/*center.m_Index[0]++;
													center.m_Index[1];*/
													 rad--;
												}
												else 
													rad=2;
												//if(j%2==0)
													lastCenter = center;
													lastAngle = maxangle;


											}
									/*		center.m_Index[0] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])
										 + (10*vcl_cos(maxangle));
											center.m_Index[1] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])
										 + (10*vcl_sin(maxangle));
											center.m_Index[2] = i;*/
											//m_RegionGrawing->VesselSegmentation(image2,center,i);
										//	m_VTKRayCaster->Intersect(center,2,&center,image2,0);
											
										/*	NewCenter.m_Index[0] = 0;
											NewCenter.m_Index[1] = 0;
											NewCenter.m_Index[2] = 0;*/
											//if(i==109)
											//for(int j = i;j>2;j--)
											int k=0;
											rad = 2;
											while(k<2000)
											
											//while(center.m_Index[2]>50)
											{
												if(m_VTKRayCaster->Intersect(center,rad,&NewCenter,image2,1,queue_CenterlinePoint)==1)
												{
												center.m_Index[0] = NewCenter.m_Index[0];
												center.m_Index[1] = NewCenter.m_Index[1];
												center.m_Index[2] = NewCenter.m_Index[2];
												rad =2;
												k++;
												}
												else
													break;
												

											}
											for(int g = 0;g < m_VTKRayCaster->index_queue; g++)
											{
												k=0;
												center = queue_CenterlinePoint[g];
												while(k<20)
											
											//while(center.m_Index[2]>50)
											{
												if(m_VTKRayCaster->Intersect(center,rad,&NewCenter,image2,1,queue_CenterlinePoint)==1)
												{
												center.m_Index[0] = NewCenter.m_Index[0];
												center.m_Index[1] = NewCenter.m_Index[1];
												center.m_Index[2] = NewCenter.m_Index[2];
												rad =2;
												k++;
												}
												else
													break;
												

											}

											}

							
											
											
											m_VTKRayCaster->Intersect(center,rad,&NewCenter,image2,1,queue_CenterlinePoint);
												center.m_Index[0] = NewCenter.m_Index[0];
												center.m_Index[1] = NewCenter.m_Index[1];
												center.m_Index[2] = NewCenter.m_Index[2];

												/*FastMarching* m_FastMarching = new FastMarching();
												m_FastMarching->Segment(coronaryOstium);
												AortaLevelsetSegmentation(image2,20,255);*/
													/*while(center.m_Index[2]>50)
											{
												
												if(m_VTKRayCaster->Intersect(center,rad,&NewCenter,image2,1)==1)
												{
												center.m_Index[0] = NewCenter.m_Index[0];
												center.m_Index[1] = NewCenter.m_Index[1];
												center.m_Index[2] = NewCenter.m_Index[2];
												rad =2;
												k++;
												}
												else
													rad++;			
											}*/
												
							 }
							 else if(i == 96)
							 {
												m_VTKRayCaster->numOfCenterlinePoint = 0;
												m_VTKRayCaster->Centerline = NULL;
											    m_VTKRayCaster->Centerline = new InputImageType::IndexType[5000];
												maxangle =angle_max_Vesselness_RCX[ maximumValue(max_Vesselness_RCX)];
												center[0] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0]+ (*itCircles)->GetRadius()[0]*vcl_cos(maxangle));
												center[1] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[1]+ (*itCircles)->GetRadius()[0]*vcl_sin(maxangle));
 												center[2] = i;
												for(int j=0;j<10;j++)
											{
												maxangle = m_VTKRayCaster->Circle_RayCasting(center,5*vnl_math::pi/4,2*vnl_math::pi,rad,image2);
												//if(abs(((maxangle+lastAngle)-vnl_math::pi))<0.1)
												//{
												//	center.m_Index[0] = center.m_Index[0]/*(long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])*/
											 //+ (2*vcl_cos(maxangle))+1;
											 //center.m_Index[1] = center.m_Index[1]/*(long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])*/
											 //+ (2*vcl_sin(maxangle));
												//center.m_Index[2] = i;
												//}
												//else{
												center.m_Index[0] = center.m_Index[0]/*(long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])*/
											 + (rad*vcl_cos(maxangle));
											 center.m_Index[1] = center.m_Index[1]/*(long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])*/
											 + (rad*vcl_sin(maxangle));
												center.m_Index[2] = i;
												//}
												
												//if(j%2==0)
													


											}
									/*		center.m_Index[0] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])
										 + (10*vcl_cos(maxangle));
											center.m_Index[1] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0])
										 + (10*vcl_sin(maxangle));
											center.m_Index[2] = i;*/
											//m_RegionGrawing->VesselSegmentation(image2,center,i);
											m_VTKRayCaster->Intersect(center,3,&center,image2,0,queue_CenterlinePoint);
											NewCenter.m_Index[0] = 0;
											NewCenter.m_Index[1] = 0;
											NewCenter.m_Index[2] = 0;
											//if(i==109)
											//for(int j = i;j>2;j--)
											rad = 2;
											int j=0;
											m_VTKRayCaster->index_queue = 0;
											

											while(j++<2500)
											//while(center.m_Index[2]>20)
											{
												if(m_VTKRayCaster->Intersect(center,rad,&NewCenter,image2,0,queue_CenterlinePoint)==1/*&&j%3!=0*/)
												{
												
												center.m_Index[0] = NewCenter.m_Index[0];
												center.m_Index[1] = NewCenter.m_Index[1];
												center.m_Index[2] = NewCenter.m_Index[2];
												InputImageType::IndexType NewTemp;
											//	m_VTKRayCaster-> FindBifurcation(center,4,&NewTemp,image2,1);
												if(center.m_Index[2]<50)
												rad =2;

												}
												/*else if(j>=2000&&m_VTKRayCaster->Intersect(center,rad,&NewCenter,image2,0)==1&&j%3!=0)
												{
												center.m_Index[0] = NewCenter.m_Index[0];
												center.m_Index[1] = NewCenter.m_Index[1];
												center.m_Index[2] = NewCenter.m_Index[2];
												}*/
											
												else
													break;


											}
											for(int g = 0;g < m_VTKRayCaster->index_queue; g++)
											{
												j=0;
												center = queue_CenterlinePoint[g];
												while(j<2000)
											
											//while(center.m_Index[2]>50)
											{
												if(m_VTKRayCaster->Intersect(center,rad,&NewCenter,image2,1,queue_CenterlinePoint)==1)
												{
												center.m_Index[0] = NewCenter.m_Index[0];
												center.m_Index[1] = NewCenter.m_Index[1];
												center.m_Index[2] = NewCenter.m_Index[2];
												rad =2;
												j++;
												}
												else
													break;
												

											}

											}
											WriteEachVessel(m_VTKRayCaster->Centerline,"g:\\CTImage\\RCX.mhd");
											Centerline_RCX = m_VTKRayCaster->Centerline;
							 }
											//m_RayCasting->draw(center,EndPoint_LCX,8);
	


	lastRadiusIndex = radiusIndex;
	if(abs(lastSumOfIntensity-(sumOfIntensity/num))>0.05*lastSumOfIntensity&&lastSumOfIntensity>100)
		FINISHED1 = true;
	std::cout<<lastSumOfIntensity<<" "<<sumOfIntensity/num<<"\n";
	lastSumOfIntensity = sumOfIntensity/num;
	if((lastRadius-(*itCircles)->GetRadius()[0]>5)&&lastRadius>0)
		FINISHED2 = true;
	lastRadius = (*itCircles)->GetRadius()[0];
	//circleIndex++;
	distThershold = 10;
}
	else if(distance>distThershold&&circleIndex==2)
	{
		FINISHED2 = true;
		  return image;
	}
	else if(distance>10&&circleIndex<2&&i<96)
	{
		
	//	distThershold = 80;
	//	max_Vesselness_RCX = 900;
		//FINISHED2 = true;
	//	circles.begin();
	//	  lastRadiusIndex[0] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
	//	  lastRadiusIndex[1] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
	//	
	//	  
	}
	

	//////////////////////////////////////////////////////////////////////////////////
	
    // Software Guide : EndCodeSnippet

    //  Software Guide : BeginLatex
    //
    //  We draw white pixels in the output image to represent each circle.
    //
    //  Software Guide : EndLatex

    // Software Guide : BeginCodeSnippet
  
    itCircles++;
    }
	if(circleIndex == 0)
	{
		FINISHED2 = true;FINISHED1 = true;
	}
  //return localOutputImage;
  
  //return image;
	return image;
  // Software Guide : EndCodeSnippet

  //  Software Guide : BeginLatex
  //
  //  We setup a writer to write out the binary image created.
  //
  //  Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  //typedef  itk::ImageFileWriter< ImageType  > WriterType;
  //WriterType::Pointer writer = WriterType::New();

 /* writer->SetFileName( argv[2] );
  writer->SetInput(localOutputImage );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }*/
//  ////////////////////////////////////////////////////////////
//
//
// 
//  QuickView viewer;
// viewer.AddImage<ImageType2d>(img);
//  viewer.Visualize();
// 
//  return EXIT_SUCCESS;
}


 
//const char * inputFilename  = "g:\\CTImage\\downsampledImage1.mhd";
int AortaDetection:: ExtractSlice(  InputImageType::IndexType *m_CircleIndex)
{
	
		
  //const char * inputFilename  = "g:\\CTImage\\downsampledImage1.mhd";
	const char * inputFilename  = "F:\\CTImage\\training\\dataset02\\closingimage.mhd";
  const char * outputFilename = "g:\\CTImage\\image2d.mhd";
  OutputImageType::Pointer image2d_1;
  OutputImageType::Pointer image2d_2;

    ReaderType::Pointer reader = ReaderType::New();
	  ReaderType::Pointer reader2 = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  this->m_CircleIndex = new InputImageType::IndexType[500];

  reader->SetFileName( inputFilename  );
  writer->SetFileName( outputFilename );
  typedef itk::ExtractImageFilter< InputImageType, OutputImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
   FilterType::Pointer filter2 = FilterType::New();
 // filter->InPlaceOn();
  //filter->SetDirectionCollapseToSubmatrix();
   reader->UpdateOutputInformation();
  InputImageType::RegionType inputRegion = reader->GetOutput()->GetLargestPossibleRegion();
  reader->Update();
  InputImageType::Pointer image1 =  reader->GetOutput();
  //reader2->SetFileName( "g:\\rawreader\\output_en2.mhd"  );
  //reader2->SetFileName( "g:\\CTImage\\vesselness.mhd"  );
  reader2->SetFileName( "F:\\CTImage\\training\\dataset02\\output_en1.mhd"  );
  //reader2->SetFileName( "g:\\CTImage\\output_en_c.mhd"  );
    //reader2->SetFileName("g:\\CTImage\\FastMarchingFilterOutput2.mhd");
  reader2->Update();
  image2 =  reader2->GetOutput();
  coronaryOstium[0] = 133;
  coronaryOstium[1] = 131;
  coronaryOstium[2] = 111;
  
 // m_FastMarching->Segment(image2,coronaryOstium);
  
  m_RegionGrawing = new RegionGrowing(image2);
   InputImageType::SizeType size = inputRegion.GetSize();
   //m_RayCasting = new RayCAsting(image2);
   int s = size[2]-1;
  size[2] = 0;
   InputImageType::IndexType start = inputRegion.GetIndex();

   for(int i=s-1;i>0/*&&(!FINISHED1&&!FINISHED2)*/;i--)
   {
	   indexCounter = 0;
		const unsigned int sliceNumber = i;/*atoi( "283" );*/
		
			std::cout<<i<<"\n";
		if(i==97)
		//	break;
			std::cout<<"Found";
		start[2] = sliceNumber;
		InputImageType::RegionType desiredRegion;
		desiredRegion.SetSize(  size  );
		desiredRegion.SetIndex( start );
		filter->SetExtractionRegion( desiredRegion );
		filter->SetInput( image1);
		filter->SetDirectionCollapseToIdentity();
		filter->Update();
		image2d_1 = filter->GetOutput();
		filter2->SetExtractionRegion( desiredRegion );
		filter2->SetInput( image2);
		filter2->SetDirectionCollapseToIdentity();
		filter2->Update();
		image2d_2 = filter2->GetOutput();
		
		
		
		OutputImageType::Pointer outImage = houghCircleTransform( image2d_1,image2d_2,i);

		//this->m_CircleIndex[indexCounter++].m_Index[2] = i;
		this->indexCounter1 = indexCounter;
		
		
		OutputImageType::Pointer InputImage = image2d_1;
		writer->SetInput( outImage );
		 try
		{
		writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
		}
		//indexCounter++;
	////////////////////////////////////////////////////
	//vtkSmartPointer<vtkImageBlend> blend =
	//vtkSmartPointer<vtkImageBlend>::New();
	//blend->AddInputConnection((vtkAlgorithmOutput*)filter->GetOutput());
	//blend->AddInputConnection(outImage);
	//blend->SetOpacity(0,.5);
	//blend->SetOpacity(1,.5);
	/////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////


 

 
//  typedef itk::CastImageFilter< FloatImageType, UnsignedintImageType > CastFilterType;
//  CastFilterType::Pointer castFilter = CastFilterType::New();
// castFilter->SetInput(InputImage);
// castFilter->Update();
// UnsignedintImageType ::Pointer imageA = castFilter->GetOutput();
//
//QuickView viewer;
//	viewer.AddImage<UnsignedintImageType>(imageA);
//  viewer.Visualize();
//
//castFilter->SetInput(outImage);
//	 castFilter->Update();
//UnsignedintImageType ::Pointer imageB = castFilter->GetOutput();
//
//  // Connect the input images
//typedef itk::OrImageFilter<	 UnsignedintImageType,	 UnsignedintImageType,	 UnsignedintImageType  >       myFilterType;
// myFilterType::Pointer Orfilter = myFilterType::New();
//
// Orfilter->SetInput1( imageA );
// 
// 
// Orfilter->SetInput2( imageB );
// UnsignedintImageType::Pointer outputImage = Orfilter->GetOutput();
//
//
//  // Execute the filter
//  Orfilter->Update();
// // // Get the Smart Pointer to the Filter Output
// //OutputImageType::Pointer outputImage = Orfilter->GetOutput();
//	//////////////////////////////////////////////////////////
//
//
//	viewer.AddImage<UnsignedintImageType>(outputImage);
//  viewer.Visualize();
	/*QuickView viewer;
	viewer.AddImage<OutputImageType>(outImage);
  viewer.Visualize();*/
   }
   ImageType::Pointer aorta = m_RegionGrawing->Opening();
	cout<<"FINISHED \n"<<
   memcpy(m_CircleIndex ,(this->m_CircleIndex),indexCounter*(sizeof(InputImageType::IndexType)));
   try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

    WriterType3D::Pointer writer3D = WriterType3D::New();
	writer3D->SetInput(aorta);
	writer3D->SetFileName("g:\\CTImage\\Aorta.mhd");
	try
    {
    writer3D->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

}

InputImageType::IndexType AortaDetection:: FindLastPoint()
{
	VesselOutput =  InputImageType::New() ;
	double maxangle =angle_max_Vesselness_LCX[ maximumValue(max_Vesselness_LCX)];
	
	for(int i=0;i<m_VTKRayCaster->numOfCenterlinePoint;i++)
	{
		VesselOutput->SetPixel(m_VTKRayCaster->Centerline[i],255);
	}

	return coronaryOstium;
}
int AortaDetection::WriteEachVessel(InputImageType::IndexType* LCX,char* fileName)
{
	InputImageType::Pointer VesselOutput = InputImageType::New();
	VesselOutput->SetRegions( image2->GetLargestPossibleRegion() );

const InputImageType::SpacingType& spacing = image2->GetSpacing();
const InputImageType::PointType& inputOrigin = image2->GetOrigin();
double outputOrigin[ 3 ];
for(unsigned int i=0; i< 3; i++)
{
outputOrigin[i] = inputOrigin[i];
}
VesselOutput->SetSpacing( spacing );
VesselOutput->SetOrigin( outputOrigin );
VesselOutput->Allocate();


	VesselOutput->FillBuffer(0);
	for (int i = 0;i<m_VTKRayCaster->numOfCenterlinePoint;i++)
	{
		VesselOutput->SetPixel(m_VTKRayCaster->Centerline[i],255);
	}
	WriterType3D::Pointer writer3D = WriterType3D::New();
	writer3D->SetInput(VesselOutput);
	writer3D->SetFileName(fileName);
	try
    {
    writer3D->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }


}

int AortaDetection:: AortaLevelsetSegmentation(InputImageType::Pointer image,int x,int y)
{
	/*double maxangle =angle_max_Vesselness_LCX[ maximumValue(max_Vesselness_LCX)];
	center[0] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0]+ (*itCircles)->GetRadius()[0]*vcl_cos(maxangle));
	center[1] = (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[1]+ (*itCircles)->GetRadius()[0]*vcl_sin(maxangle));
	center[2] = 109;*/
	FastMarchingLevelSetBase* m_FastMarchingLevelSetBase = new FastMarchingLevelSetBase();
	m_FastMarchingLevelSetBase->LoadInputImage("g:\\CTImage\\aorta.mhd");
	m_FastMarchingLevelSetBase->AddSeed(coronaryOstium);
	m_FastMarchingLevelSetBase->m_SeedValue = 20;
	m_FastMarchingLevelSetBase->RunFastMarching();
	m_FastMarchingLevelSetBase->SetStoppingValue(50);
	m_FastMarchingLevelSetBase->SaveOutputImage("g:\\CTImage\\fast.mhd");




  return 0;


}