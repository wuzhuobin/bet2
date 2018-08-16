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
// me 
#include "vtkBet2Filter.h"
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
	{}
	void cleanup()
	{}
	void Update123() 
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
		sampleFunction->SetSampleDimensions(127, 127, 127); // intentional NPOT dimensions.
		sampleFunction->SetModelBounds(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
		sampleFunction->SetCapping(false);
		sampleFunction->SetComputeNormals(false);
		sampleFunction->SetScalarArrayName("values");
		sampleFunction->Update();

		vtkDataArray* a = sampleFunction->GetOutput()->GetPointData()->GetScalars("values");
		double range[2];
		a->GetRange(range);

		vtkSmartPointer<vtkImageShiftScale> t =
			vtkSmartPointer<vtkImageShiftScale>::New();
		t->SetInputConnection(sampleFunction->GetOutputPort());

		t->SetShift(-range[0]);
		double magnitude = range[1] - range[0];
		if (magnitude == 0.0)
		{
			magnitude = 1.0;
		}
		t->SetScale(255.0 / magnitude);
		t->SetOutputScalarTypeToUnsignedChar();

		t->Update();
		vtkImageData *imageData = t->GetOutput();
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
};
#include "TestvtkBet2Filter.moc"
QTEST_APPLESS_MAIN(TestvtkBet2Filter);