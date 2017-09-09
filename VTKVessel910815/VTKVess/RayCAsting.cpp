#include "RayCAsting.h"





RayCAsting::RayCAsting(InputImageType::Pointer image)
{
	Ka = 0.2;
	Kd = 0.6;
	Ks = 0.2;
	this->image = image;
	rescaledImage = ChangeGrayToRGB(this->image);

}

RayCAsting::~RayCAsting(void)
{
}

InputImageType::Pointer RayCAsting::ChangeGrayToRGB(InputImageType::Pointer image)
{
  

  // Software Guide : EndCodeSnippet


  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  typedef itk::ImageFileWriter< InputImageType >  WriterType;
   WriterType::Pointer writer = WriterType::New();
  
  typedef itk::CastImageFilter<
               InputImageType, InputImageType >  CastFilterType;

  typedef itk::RescaleIntensityImageFilter<
               InputImageType, InputImageType >  RescaleFilterType;

  typedef itk::ShiftScaleImageFilter<
               InputImageType, InputImageType >  ShiftScaleFilterType;

  typedef itk::NormalizeImageFilter<
               InputImageType, InputImageType >  NormalizeFilterType;
  // Software Guide : EndCodeSnippet

  ReaderType::Pointer reader = ReaderType::New();





  // Software Guide : BeginCodeSnippet
  CastFilterType::Pointer       castFilter       = CastFilterType::New();
  RescaleFilterType::Pointer    rescaleFilter    = RescaleFilterType::New();
  ShiftScaleFilterType::Pointer shiftFilter      = ShiftScaleFilterType::New();
  NormalizeFilterType::Pointer  normalizeFilter = NormalizeFilterType::New();
  // Software Guide : EndCodeSnippet


  reader->SetFileName( "g:\\rawreader\\output_en2.mhd" );


  
  //castFilter->SetInput(       reader->GetOutput() );
  //shiftFilter->SetInput(      reader->GetOutput() );
  rescaleFilter->SetInput(    image );
  normalizeFilter->SetInput( rescaleFilter->GetOutput() );
  
  // Software Guide : BeginCodeSnippet
  rescaleFilter->SetOutputMinimum(  0 );
  rescaleFilter->SetOutputMaximum( 255 );
  // Software Guide : EndCodeSnippet


  //
  //  Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  //shiftFilter->SetScale( 1.2 );
  //shiftFilter->SetShift( 25 );
 

  // Software Guide : BeginCodeSnippet
  //castFilter->Update();
 // shiftFilter->Update();
  rescaleFilter->Update();
  normalizeFilter->Update();
  // Software Guide : EndCodeSnippet
    writer->SetFileName("G:\\CTImage\\RescaledIntensityImage.mhd");
  writer->SetInput( normalizeFilter->GetOutput() );
		 try
		{
		writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
//		return EXIT_FAILURE;
		}
		return normalizeFilter->GetOutput();
 // return EXIT_SUCCESS;
}


void RayCAsting::  SphereRayCasting(InputImageType::IndexType center,InputImageType::IndexType VoxelIndex,float radius)
{
	InputImageType::PixelType LightIntensity ;
	InputImageType::PixelType Is;
	InputImageType::IndexType LightVector;
	InputImageType::IndexType NormalVector;
	InputImageType::IndexType V;
	InputImageType::IndexType R;
	LightIntensity = 0.75;


	LightVector.m_Index[0] =  VoxelIndex.m_Index[0] - center.m_Index[0];
	LightVector.m_Index[1] =  VoxelIndex.m_Index[1] - center.m_Index[1];
	LightVector.m_Index[2] =  VoxelIndex.m_Index[2] - center.m_Index[2];
	float* LightV;
	LightV[0] = LightVector.m_Index[0];
	LightV[1] = LightVector.m_Index[1];
	LightV[2] = LightVector.m_Index[2];
	LightV = normalize(LightV);
	float L[3];
	L[0] = LightV[0];
	L[1] = LightV[1];
	L[2] = LightV[2];

	float* Gradient;
	InputImageType::IndexType tempIndex1,tempIndex2;
	tempIndex1 = VoxelIndex;
	tempIndex2 = VoxelIndex;
	tempIndex1.m_Index[0] = tempIndex1.m_Index[0]-1;
	tempIndex2.m_Index[0] = tempIndex2.m_Index[0]+1;
	Gradient[0] = rescaledImage->GetPixel(tempIndex2)-rescaledImage->GetPixel(tempIndex1);
	std::cout<<rescaledImage->GetPixel(tempIndex2)<<" "<<rescaledImage->GetPixel(tempIndex1)<<" "<<rescaledImage->GetPixel(center)<<"\n";
	tempIndex1 = VoxelIndex;
	tempIndex2 = VoxelIndex;
	tempIndex1.m_Index[1] = tempIndex1.m_Index[1]-1;
	tempIndex2.m_Index[1] = tempIndex2.m_Index[1]+1;
	Gradient[1] = rescaledImage->GetPixel(tempIndex2)-rescaledImage->GetPixel(tempIndex1);
	tempIndex1 = VoxelIndex;
	tempIndex2 = VoxelIndex;
	tempIndex1.m_Index[2] = tempIndex1.m_Index[2]-1;
	tempIndex2.m_Index[2] = tempIndex2.m_Index[2]+1;
	Gradient[2] = rescaledImage->GetPixel(tempIndex2)-rescaledImage->GetPixel(tempIndex1);
	Gradient = normalize(Gradient);
	float H[3] ;
	//H[0] = 1*LightV[0];
	//H[1] = 1*LightV[1];
	//H[2] = 1*LightV[2];
	double fSpecular = vcl_pow( dotProduce(L,Gradient),5);
	double color = (Ka+(Kd*dotProduce(L,Gradient))+(Ks*fSpecular)/*dot_product(R,V)*Is)*/);

	

	
	/*float fDiffuse=saturate(Gradient, LightVector));
	float3 H = normalize(LightVector+DirectionVector) ;
	float fSpecular = exp(saturate(Gradient, H), 5)
	Current.rgb *= (0.2+0.6*fDiffuse+0.2*fSpecular);*/


}

float RayCAsting::dotProduce(float* A,/*InputImageType::PixelType*/float* B)
{
	float result = 0;
	result = (A[0]*B[0])+(A[1]*B[1])+(A[2]*B[2]);
	if(result<0)
		result = 0;
	else if(result >1)
		result = 1;


	return result;
}

float* RayCAsting::normalize(float* A)
{
	float d = vcl_sqrt((float)(A[0]*A[0])+(A[1]*A[1])+(A[2]*A[2]));
	float result[3] ;
	result[0] = A[0]/d;
	result[1] = A[1]/d;
	result[2] = A[2]/d;
	return result;
}

void RayCAsting:: draw(InputImageType::IndexType center,InputImageType::IndexType pointEnd,float radius ) 
 {
	 InputImageType::IndexType tempIndex;
    //ofstream imageFile(outputName,ios_base::binary);
    //if (!imageFile)
    //    return false; 
    //// Addition of the TGA header
    //imageFile.put(0).put(0);
    //imageFile.put(2);        /* RGB not compressed */

    //imageFile.put(0).put(0);
    //imageFile.put(0).put(0);
    //imageFile.put(0);

    //imageFile.put(0).put(0); /* origin X */ 
    //imageFile.put(0).put(0); /* origin Y */

    //imageFile.put((unsigned char)(myScene.sizex & 0x00FF)).put((unsigned char)((myScene.sizex & 0xFF00) / 256));
    //imageFile.put((unsigned char)(myScene.sizey & 0x00FF)).put((unsigned char)((myScene.sizey & 0xFF00) / 256));
    //imageFile.put(24);       /* 24 bit bitmap */ 
    //imageFile.put(0); 
    // end of the TGA header 

    // Scanning 

	point3D pos;
	pos._x = pointEnd.m_Index[0];
	pos._y = pointEnd.m_Index[1];
	pos._z = pointEnd.m_Index[2];

	 point3D Gradient;
	InputImageType::IndexType tempIndex1,tempIndex2;
	tempIndex1 = pointEnd;
	tempIndex2 = pointEnd;
	tempIndex1.m_Index[0] = tempIndex1.m_Index[0]-1;
	tempIndex2.m_Index[0] = tempIndex2.m_Index[0]+1;
	Gradient._x = abs(rescaledImage->GetPixel(tempIndex2)-rescaledImage->GetPixel(tempIndex1))/255;
	std::cout<<rescaledImage->GetPixel(tempIndex2)<<" "<<rescaledImage->GetPixel(tempIndex1)<<" "<<rescaledImage->GetPixel(center)<<"\n";
	tempIndex1 = pointEnd;
	tempIndex2 = pointEnd;
	tempIndex1.m_Index[1] = tempIndex1.m_Index[1]-1;
	tempIndex2.m_Index[1] = tempIndex2.m_Index[1]+1;
	Gradient._y = abs(rescaledImage->GetPixel(tempIndex2)-rescaledImage->GetPixel(tempIndex1))/255;
	tempIndex1 = pointEnd;
	tempIndex2 = pointEnd;
	tempIndex1.m_Index[2] = tempIndex1.m_Index[2]-1;
	tempIndex2.m_Index[2] = tempIndex2.m_Index[2]+1;
	Gradient._z = abs(rescaledImage->GetPixel(tempIndex2)-rescaledImage->GetPixel(tempIndex1))/255;
	Gradient.normalize();


   /* for (int y = 0; y < myScene.sizey; ++y) { 
    for (int x = 0; x < myScene.sizex; ++x) {*/
        float red = 0, green = 0, blue = 0;
        float coef = 1.0f;
        int level = 0; 
		lambert = 0;
		color =0;
        // Cast the ray 
        // Because we are not in conic perspective (point of origin)
        // but in orthographic perspective, there is no natural starting point
        // We have to put it far enough to enclose the whole scene
        // but not too far to avoid floating point precision problems (acne and other)
        // 1000.0f seems like a good compromise for now..
        ray viewRay;
		viewRay.start._x = center.m_Index[0];
		viewRay.start._y = center.m_Index[1];
		viewRay.start._z = center.m_Index[2];
		//viewRay.start.normalize();
		viewRay.dir._x = pointEnd.m_Index[0] - center.m_Index[0];
		viewRay.dir._y = pointEnd.m_Index[1] - center.m_Index[1];
		viewRay.dir._z = pointEnd.m_Index[2] - center.m_Index[2];
		//viewRay.dir.normalize();


		
      //  do 
        { 
            // Looking for the closest intersection
            float t = 1.0f;
           /* int currentSphere= -1;

            for (unsigned int i = 0; i < myScene.sphereContainer.size(); ++i) 
            { 
                if (hitSphere(viewRay, myScene.sphereContainer[i], t)) 
                {
                    currentSphere = i;
                }
            }

            if (currentSphere == -1)
                break;*/

			point3D newStart;
			newStart = viewRay.start +  viewRay.dir*t; 

            // What is the normal vector at the point of intersection ?
            // It's pretty simple because we're dealing with spheres
            point3D n;
			
			n =  pos - newStart/*myScene.sphereContainer[currentSphere].pos*/;
			n.normalize();
            /*float temp = (n.dot( n));
            if (temp == 0.0f) 
                break; */
			
           /* temp = 1.0f / sqrtf(temp); 
            n =  n*temp; */

            //material currentMat = myScene.materialContainer[myScene.sphereContainer[currentSphere].materialId]; 
			material currentMat;
			tempIndex.m_Index[0] = newStart._x;
			tempIndex.m_Index[1] = newStart._y;
			tempIndex.m_Index[2] = newStart._z;
			currentMat.color= image->GetPixel(pointEnd);
			currentMat.reflection = 0.6;
            // calcul de la valeur d'éclairement au point 
         //   for (unsigned int j = 0; j < myScene.lightContainer.size(); ++j) 
			//{
                light current ;
				current.pos._x = center.m_Index[0];
				current.pos._y = center.m_Index[1];
				current.pos._z = center.m_Index[2];
				current.color = 255;/*myScene.lightContainer[j];*/
                point3D dist;
				dist = current.pos - newStart;
				/*dist[0] = current.pos[0] - newStart[0];
				dist[1] = current.pos[1] - newStart[1];
				dist[2] = current.pos[2] - newStart[2];*/
				//dist.normalize();
               /* if ((n.dot( dist)) <= 0.0f)
                    continue;*/
                 t = sqrtf((dist.dot( dist)));
               /* if ( t <= 0.0f )
                    continue;*/
                ray lightRay;
				lightRay.start = viewRay.start;
				lightRay.dir =   viewRay.dir/*dist *(1/t)*/;
                // computation of the shadows
                bool inShadow = false; 
            /*    for (unsigned int i = 0; i < myScene.sphereContainer.size(); ++i) {
                    if (hitSphere(lightRay, myScene.sphereContainer[i], t)) {
                        inShadow = true;
                        break;
                    }
                }*/
                if (!inShadow) {
                    // lambert
                    lambert = (lightRay.dir.dot( n)) * coef;
                    color += lambert * current.color * currentMat.color;
                   /* green += lambert * current.green * currentMat.green;
                    blue += lambert * current.blue * currentMat.blue;*/
                }
            
			
			//float fDiffuse=saturate(Gradient, LightVector));
			//float3 H = normalize(LightVector+DirectionVector) ;
			//float fSpecular = exp(saturate(Gradient, H), 5)
			//Current.rgb *= (0.2+0.6*fDiffuse+0.2*fSpecular);
			light Current;

			float fDiffuse=Gradient. dot(lightRay.dir);
			point3D H =  lightRay.dir+viewRay.dir ;
			H.normalize();
			float fSpecular = vcl_pow((Gradient.dot(H)), 5);
			//current.color =
			Current.color = (0.2+0.6*fDiffuse+0.2*fSpecular);

			
			/*float reflet = 2.0f * (lightRay.dir * n);
			point3D phongDir = lightRay.dir - reflet * n;
			 float phongTerm = _MAX(phongDir * viewRay.dir, 0.0f) ;
			 phongTerm = currentMat.specvalue * powf(phongTerm, currentMat.specpower) * coef;
			 color = phongTerm * current.color;*/


            // We iterate on the next reflection
            coef = coef * currentMat.reflection;
			double reflet = 2.0f * (viewRay.dir.dot( n));
            viewRay.start = newStart;
            viewRay.dir = viewRay.dir -   n*reflet;

            level++;
        } 
	//		}
     //   while ((coef > 0.0f) && (level < 10));   

        //imageFile.put((unsigned char)min(blue*255.0f,255.0f)).put((unsigned char)min(green*255.0f, 255.0f)).put((unsigned char)min(red*255.0f, 255.0f));
    }
//    }
   // return true;
	
 //}
