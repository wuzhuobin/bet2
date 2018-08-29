// qt
#include <QTest>
// vtk
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkSphere.h>
#include <vtkSampleFunction.h>
#include <vtkImageShiftScale.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkNIFTIImageReader.h>
#include <vtkNIFTIImageWriter.h>
#include <vtkImageExport.h>
// bet2 
#include <newimage.h>
#include <newimageio.h>
// me 
#include "vtkBet2Filter.h"
void createSphereImage(vtkImageData *input)
{
	// Create a spherical implicit function.

	vtkSmartPointer<vtkSphere> sphere =
		vtkSmartPointer<vtkSphere>::New();
	sphere->SetRadius(0.1);
	sphere->SetCenter(0.0, 0.0, 0.0);
	vtkSmartPointer<vtkSampleFunction> sampleFunction =
		vtkSmartPointer<vtkSampleFunction>::New();
	sampleFunction->SetImplicitFunction(sphere);
	sampleFunction->SetOutputScalarTypeToDouble();
	sampleFunction->SetSampleDimensions(100, 100, 100); // intentional NPOT dimensions.
	sampleFunction->SetModelBounds(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
	sampleFunction->SetCapping(false);
	sampleFunction->SetComputeNormals(false);
	sampleFunction->SetScalarArrayName("values");
	sampleFunction->Update();
	vtkDataArray* a = sampleFunction->GetOutput()->GetPointData()->GetScalars("values");
	double range[2];
	a->GetRange(range);
	vtkSmartPointer<vtkImageShiftScale> imageShiftScale =
		vtkSmartPointer<vtkImageShiftScale>::New();
	imageShiftScale->SetInputConnection(sampleFunction->GetOutputPort());
	imageShiftScale->SetShift(-range[0]);
	double magnitude = range[1] - range[0];
	if (magnitude == 0.0)
	{
		magnitude = 1.0;
	}
	imageShiftScale->SetScale(255.0 / magnitude);
	imageShiftScale->SetOutputScalarTypeToDouble();
	imageShiftScale->Update();
	input->ShallowCopy(imageShiftScale->GetOutput());
}
class TestvtkBet2Filter : public QObject
{
	Q_OBJECT;
private Q_SLOTS:
	void initTestCase()
	{
	}
	void cleanupTestCase()
	{
	}
	void init()
	{
	}
	void cleanup()
	{
	}
	void vtkImageDataToNEWIMAGE_volume()
	{
		vtkSmartPointer<vtkImageData> imageData =
			vtkSmartPointer<vtkImageData>::New();
		createSphereImage(imageData);
		imageData->SetSpacing(3, 2, 1);
		imageData->SetOrigin(4, 5, 6);
		int sx = imageData->GetDimensions()[0];
		int sy = imageData->GetDimensions()[1];
		int sz = imageData->GetDimensions()[2];
		int *extent = imageData->GetExtent();
		NEWMAT::Matrix qform(NEWMAT::IdentityMatrix(4));
		qform(1, 4) = imageData->GetOrigin()[0];
		qform(2, 4) = imageData->GetOrigin()[1];
		qform(3, 4) = imageData->GetOrigin()[2];
		NEWIMAGE::volume<double> volume;
//////////////////////////////////////// Set Pointer Address //////////////////////////////////////////////////
		volume.reinitialize(
			sx,	sy,	sz,
			reinterpret_cast<double*>(imageData->GetScalarPointer()), false);
//////////////////////////////////////// Set Pointer Address //////////////////////////////////////////////////
//////////////////////////////////////// Copy Pointer content////////////////////////////////////////
		//volume.reinitialize(sx, sy, sz, 0, true);
		//for (int x = imageData->GetExtent()[0]; x <= imageData->GetExtent()[1]; ++x) {
		//	for (int y = imageData->GetExtent()[2]; y <= imageData->GetExtent()[3]; ++y) {
		//		for (int z = imageData->GetExtent()[4]; z <= imageData->GetExtent()[5]; ++z) {
		//			volume(x, y, z) = imageData->GetScalarComponentAsDouble(x, y, z, 0);
		//		}
		//	}
		//}
//////////////////////////////////////// Copy Pointer content////////////////////////////////////////
		volume.setdims(
			imageData->GetSpacing()[0],
			imageData->GetSpacing()[1],
			imageData->GetSpacing()[2]);
		volume.set_qform(NIFTI_XFORM_ALIGNED_ANAT, qform);
		for (int x = imageData->GetExtent()[0]; x <= imageData->GetExtent()[1]; ++x) {
			for (int y = imageData->GetExtent()[2]; y <= imageData->GetExtent()[3]; ++y) {
				for (int z = imageData->GetExtent()[4]; z <= imageData->GetExtent()[5]; ++z) {
					//QCOMPARE(imageData->GetScalarComponentAsDouble(x, y, z, 0), data[x*sx*sy + y*sy + z]);
					QCOMPARE(imageData->GetScalarComponentAsDouble(x, y, z, 0), volume.value(x, y, z));
				}
			}
		}
		QCOMPARE(imageData->GetSpacing()[0], (double)(volume.xdim()));
		QCOMPARE(imageData->GetSpacing()[1], (double)(volume.ydim()));
		QCOMPARE(imageData->GetSpacing()[2], (double)(volume.zdim()));
		QCOMPARE(imageData->GetOrigin()[0], (double)(volume.qform_mat()(1, 4)));
		QCOMPARE(imageData->GetOrigin()[1], (double)(volume.qform_mat()(2, 4)));
		QCOMPARE(imageData->GetOrigin()[2], (double)(volume.qform_mat()(3, 4)));
	}
	void Update() 
	{
		QSKIP("vtkBet2Filter is not yet completed for foo input");
		vtkSmartPointer<vtkImageData> imageData = 
			vtkSmartPointer<vtkImageData>::New();
		createSphereImage(imageData);
		imageData->SetSpacing(3, 2, 1);
		vtkSmartPointer<vtkBet2Filter> f =
			vtkSmartPointer<vtkBet2Filter>::New();
		f->SetInputData(imageData);
		f->Update();
		vtkImageData *output = f->GetOutput();
		QCOMPARE(output->GetExtent()[0], 0);
		QCOMPARE(output->GetExtent()[1], 100);
		QCOMPARE(output->GetExtent()[2], 0);
		QCOMPARE(output->GetExtent()[3], 100);
		QCOMPARE(output->GetExtent()[4], 0);
		QCOMPARE(output->GetExtent()[5], 100);
		//QCOMPARE(output->GetOrigin()[0], 1.0);
		//QCOMPARE(output->GetOrigin()[1], 2.0);
		//QCOMPARE(output->GetOrigin()[2], 3.0);
		QCOMPARE(output->GetSpacing()[0], 3.0);
		QCOMPARE(output->GetSpacing()[1], 2.0);
		QCOMPARE(output->GetSpacing()[2], 1.0);
	}
	void Update2()
	{
		//QSKIP("Since the input data is hard coded. Please change it to your local path");
		vtkSmartPointer<vtkNIFTIImageReader> reader =
			vtkSmartPointer<vtkNIFTIImageReader>::New();
		reader->SetFileName("C:/Users/jieji/Desktop/T2_RTHANDMOTOR_BOLD/20130610_144057T2AXTE80SENSEs301a1003.nii");
		reader->Update();
		vtkSmartPointer<vtkNIFTIImageWriter> writer =
			vtkSmartPointer<vtkNIFTIImageWriter>::New();
		writer->SetInputConnection(reader->GetOutputPort());
		writer->SetFileName("vtkNiftiImage.nii.gz");
		writer->Write();
		vtkSmartPointer<vtkBet2Filter> bet2 =
			vtkSmartPointer<vtkBet2Filter>::New();
		bet2->SetInputConnection(reader->GetOutputPort());
		bet2->Update();
		writer->SetInputConnection(bet2->GetOutputPort());
		writer->SetFileName("mask.nii.gz");
		writer->Write();
		std::cin.get();
	}
};
// testing
#include "TestvtkBet2Filter.moc"
QTEST_APPLESS_MAIN(TestvtkBet2Filter);