// me
#include "vtkBet2Filter.h"
#include "bet2.h"
// vtk
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
}

vtkBet2Filter::~vtkBet2Filter()
{
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
	// convert vtkImageData to NEWIMAGE::volume<float>
	volume<float> testvol;
//////////////////////////////////////// Copy Pointer content////////////////////////////////////////
	// allocate memory
	testvol.reinitialize(
		input->GetDimensions()[0],
		input->GetDimensions()[1],
		input->GetDimensions()[2],
		0, true);
	// doing copy voxel by voxel. 
	// should be more effective way.
	for (int x = input->GetExtent()[0]; x <= input->GetExtent()[1]; ++x) {
		for (int y = input->GetExtent()[2]; y <= input->GetExtent()[3]; ++y) {
			for (int z = input->GetExtent()[4]; z <= input->GetExtent()[5]; ++z) {
				testvol(x, y, z) = input->GetScalarComponentAsFloat(x, y, z, 0);
			}
		}
	}
//////////////////////////////////////// Copy Pointer content////////////////////////////////////////
//////////////////////////////////////// Set Pointer Address //////////////////////////////////////////////////
	//testvol.reinitialize(
	//	input->GetDimensions()[0],
	//	input->GetDimensions()[1],
	//	input->GetDimensions()[2],
	//	reinterpret_cast<float*>(input->GetScalarPointer()), false);
	//float *data = new float[
	//		input->GetDimensions()[0] *
	//		input->GetDimensions()[1] *
	//		input->GetDimensions()[2]];
	//vtkImageExport *imageExport = vtkImageExport::New();
	//imageExport->SetInputData(input);
	//imageExport->Update();
	//imageExport->Export(data);
	//testvol.reinitialize(
	//	imageExport->GetDataDimensions()[0],
	//	imageExport->GetDataDimensions()[1],
	//	imageExport->GetDataDimensions()[2],
	//	data, true);
	//imageExport->Delete();
	//cout << "valueImage: " << input->GetScalarComponentAsFloat(996, 512, 13, 0) << '\n';
	//cout << "array: " << data[13 * input->GetDimensions()[0] * input->GetDimensions()[1] + 512 * input->GetDimensions()[0] + 996];
	//cout << "ValueVolume: " << testvol.value(996, 512, 13) << '\n';
//////////////////////////////////////// Set Pointer Address //////////////////////////////////////////////////
	// setting spacing 
	testvol.setdims(
		input->GetSpacing()[0],
		input->GetSpacing()[1],
		input->GetSpacing()[2]);
	// setting q-form, which is similiar to roientation 
	// since vtkImageData have not orientation info, only setting its origin to 
	// its translation component and setting orientation matrix to indentity.
	NEWMAT::Matrix matrix(NEWMAT::IdentityMatrix(4));
	matrix(1, 4) = input->GetOrigin()[0];
	matrix(2, 4) = input->GetOrigin()[1];
	matrix(3, 4) = input->GetOrigin()[2];
	testvol.set_qform(NIFTI_XFORM_ALIGNED_ANAT, matrix);
//////////////////////////////////////// computation ////////////////////////////////////////
	double xarg = 0, yarg = 0, zarg = 0;
	const double bet_main_parameter = pow(0.5, .275);
	// 2D kludge (worked for bet, but not here in bet2, hohum)
	if (testvol.xsize()*testvol.xdim() < 20) testvol.setxdim(200);
	if (testvol.ysize()*testvol.ydim() < 20) testvol.setydim(200);
	if (testvol.zsize()*testvol.zdim() < 20) testvol.setzdim(200);
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
	for (int i = 0; i < nb_iter; i++)
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
		for (int i = 0; i < nb_iter; i++)
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
	// coverting NEWIMAGE::volume<short> to vtkImageData
	volume<short> outputFSL = static_cast<short>(1) - brainmask;
	// since vtkImageData have no orientation but origin 
//////////////////////////////////////// computation ////////////////////////////////////////
	// only copy translation components. 
	output->SetOrigin(
		brainmask.qform_mat()(1, 4),
		brainmask.qform_mat()(2, 4),
		brainmask.qform_mat()(3, 4));
	// allocate memory.
	output->SetDimensions(
		brainmask.xsize(), 
		brainmask.ysize(), 
		brainmask.zsize());
	output->AllocateScalars(VTK_SHORT, 1);
	output->SetSpacing(outputFSL.xdim(), outputFSL.ydim(), outputFSL.zdim());
	for (int x = output->GetExtent()[0]; x <= output->GetExtent()[1]; ++x) {
		for (int y = output->GetExtent()[2]; y <= output->GetExtent()[3]; ++y) {
			for (int z = output->GetExtent()[4]; z <= output->GetExtent()[5]; ++z) {
				*static_cast<short*>(output->GetScalarPointer(x, y, z)) = outputFSL(x, y, z);
			}
		}
	}
	return 1;
}