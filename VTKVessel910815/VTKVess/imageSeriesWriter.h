#pragma once
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"
class imageSeriesWriter
{
public:
	imageSeriesWriter(void);
	~imageSeriesWriter(void);
	BOOL write();
};
