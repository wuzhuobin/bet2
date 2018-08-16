// me
#include "vtkBet2Filter.h"
#include "bet2.h"
// vtk
#include <vtkNIFTIImageReader.h>
#include <vtkNIFTIImageWriter.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkImageData.h>
#include <vtkObjectFactory.h>
vtkStandardNewMacro(vtkBet2Filter);
void vtkBet2Filter::PrintSelf(ostream & os, vtkIndent indent)
{
	Superclass::PrintSelf(os, indent);
}

vtkBet2Filter::vtkBet2Filter()
{
	this->Reader = vtkNIFTIImageReader::New();
	this->Writer = vtkNIFTIImageWriter::New();
}

vtkBet2Filter::~vtkBet2Filter()
{
	this->Reader->Delete();
	this->Writer->Delete();
}

int vtkBet2Filter::RequestData(vtkInformation * request, vtkInformationVector ** inputVector, vtkInformationVector * outputVector)
{
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and ouptut
	vtkImageData *input = vtkImageData::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	if (input == nullptr) {
		vtkErrorMacro(<< "The input is not vtkImageData type. ");
		return 0;
	}
	vtkImageData *output = vtkImageData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));
	std::string fileName = std::to_string((unsigned long long)this) + ".nii.gz";
	this->Writer->SetInputData(input);
	this->Writer->SetFileName(fileName.c_str());
	this->Writer->Write();

	volume<float> testvol;
	if (read_volume(testvol, fileName) < 0) {
		vtkErrorMacro(<< "Bet read_volume falied. ");
		return 0;
	}
	if (std::remove(fileName.c_str()) != 0) {
		vtkErrorMacro(<< "Removing file " << fileName << " failed. ");
		return 0;
	}
	double xarg = 0, yarg = 0, zarg = 0;
	const double bet_main_parameter = pow(0.5, .275);
	// 2D kludge (worked for bet, but not here in bet2, hohum)
	if (testvol.xsize()*testvol.xdim()<20) testvol.setxdim(200);
	if (testvol.ysize()*testvol.ydim()<20) testvol.setydim(200);
	if (testvol.zsize()*testvol.zdim()<20) testvol.setzdim(200);
	using namespace bet2;
	Mesh m;
	make_mesh_from_icosa(5, m);
	bet_parameters bp = adjust_initial_mesh(testvol, m, 0.0, xarg, yarg, zarg);
	Mesh moriginal = m;

	const double rmin = 3.33 * 1.0;
	const double rmax = 10 * 1.0;
	const double E = (1 / rmin + 1 / rmax) / 2.;
	const double F = 6. / (1 / rmin - 1 / rmax);
	const int nb_iter = 1000;
	const double self_intersection_threshold = 4000;

	double l = 0;
	for (int i = 0; i<nb_iter; i++)
	{
		step_of_computation(testvol, m, bet_main_parameter, 0, 0, i, l, bp.t2, bp.tm, bp.t, E, F, bp.cog.Z, bp.radius, 0.0);
	}

	double tmp = m.self_intersection(moriginal);
	bool self_intersection;
	self_intersection = m.real_self_intersection();
	int pass = 0;
	//self-intersection treatment
	while (m.real_self_intersection())
	{
		m = moriginal;
		l = 0;
		pass++;
		for (int i = 0; i<nb_iter; i++)
		{
			double incfactor = pow(10.0, (double)pass + 1);
			if (i > .75 * (double)nb_iter)
				incfactor = 4.*(1. - i / (double)nb_iter) * (incfactor - 1.) + 1.;
			step_of_computation(testvol, m, bet_main_parameter, pass, incfactor, i, l, bp.t2, bp.tm, bp.t, E, F, bp.cog.Z, bp.radius, 0.0);
		}
		if (pass == 10) // give up
			self_intersection = 0;
	}
	//display
	volume<short> brainmask = make_mask_from_mesh(testvol, m);
	if (save_volume(static_cast<short>(1) - brainmask, fileName) < 0) {
		vtkErrorMacro(<< "Bet2 save_volume failed. ");
		return 0;
	}
	this->Reader->SetFileName(fileName.c_str());
	this->Reader->Update();
	output->ShallowCopy(this->Reader->GetOutput());
	output->SetSpacing(input->GetSpacing());
	output->SetOrigin(input->GetOrigin());
	if (std::remove(fileName.c_str()) != 0) {
		vtkErrorMacro(<< "Removing file " << fileName << " failed. ");
		return 0;
	}
	return 1;
}