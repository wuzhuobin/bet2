#ifndef __VTK_BET_FILTER_H__
#define __VTK_BET_FILTER_H__
#pragma once
// vtk 
#include <vtkImageAlgorithm.h>
class vtkNIFTIImageReader;
class vtkNIFTIImageWriter;
class vtkBet2Filter : public vtkImageAlgorithm
{
public:
	static vtkBet2Filter *New();
	vtkTypeMacro(vtkImageAlgorithm, vtkBet2Filter);
	virtual void PrintSelf(ostream &os, vtkIndent indent) override;
protected:
	vtkBet2Filter();
	virtual ~vtkBet2Filter() override;
	virtual int RequestData(vtkInformation *request,
		vtkInformationVector **inputVector,
		vtkInformationVector *outputVector) override;
	vtkNIFTIImageReader *Reader;
	vtkNIFTIImageWriter *Writer;
private:
	vtkBet2Filter(const vtkBet2Filter&) VTK_DELETE_FUNCTION;
	void operator=(const vtkBet2Filter&) VTK_DELETE_FUNCTION;
};
#endif // !__VTK_BET_FILTER_H__
