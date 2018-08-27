#ifndef __VTK_BET_FILTER_H__
#define __VTK_BET_FILTER_H__
#pragma once
// vtk 
#include <vtkImageAlgorithm.h>
class vtkBet2Filter : public vtkImageAlgorithm
{
public:
	static vtkBet2Filter *New();
	vtkTypeMacro(vtkBet2Filter, vtkImageAlgorithm);
	virtual void PrintSelf(ostream &os, vtkIndent indent) VTK_OVERRIDE;
protected:
	vtkBet2Filter();
	virtual ~vtkBet2Filter() VTK_OVERRIDE;
	virtual int RequestData(vtkInformation *request,
		vtkInformationVector **inputVector,
		vtkInformationVector *outputVector) VTK_OVERRIDE;
private:
	vtkBet2Filter(const vtkBet2Filter&) VTK_DELETE_FUNCTION;
	void operator=(const vtkBet2Filter&) VTK_DELETE_FUNCTION;
};
#endif // !__VTK_BET_FILTER_H__
