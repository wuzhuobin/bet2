// qt
#include <QTest>
#include <QDebug>
// vtk
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkStripper.h>
#include <vtkTriangleFilter.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
// me 
#include "vtkPolyDataNormalsCentroids.h"
// std
#include <sstream>
#include <iostream>
class TestvtkPolyDataNormalsCentroids : public QObject
{
	Q_OBJECT;
private Q_SLOTS:
	void initTestCase()
	{
	}
	void cleanupTestCase()
	{
		std::cin.get();
	}
	void init()
	{
	}
	void cleanup()
	{
	}
	void test() {
		vtkSmartPointer<vtkSphereSource> sphereSource =
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetThetaResolution(3);
		sphereSource->SetPhiResolution(3);
		sphereSource->SetCenter(0.0, 0.0, 0.0);
		sphereSource->Update();

		vtkSmartPointer<vtkStripper> stripper =
			vtkSmartPointer<vtkStripper>::New();
		stripper->SetInputConnection(sphereSource->GetOutputPort());
		stripper->Update();
		
		vtkSmartPointer<vtkTriangleFilter> triangleFilter =
			vtkSmartPointer<vtkTriangleFilter>::New();
		triangleFilter->SetInputConnection(stripper->GetOutputPort());
		triangleFilter->Update();
		vtkPolyData *sphere = triangleFilter->GetOutput();

		vtkSmartPointer<vtkPolyDataNormalsCentroids> normalsCentroids =
			vtkSmartPointer<vtkPolyDataNormalsCentroids>::New();
		normalsCentroids->SetInputData(sphere);
		normalsCentroids->SetComputePointCentroids(true);
		normalsCentroids->SetComputePointNormals(true);
		normalsCentroids->SetSplitting(false);
		normalsCentroids->Update();
		sphere = normalsCentroids->GetOutput();
		vtkSmartPointer<vtkPolyDataWriter> writer =
			vtkSmartPointer<vtkPolyDataWriter>::New();
		writer->SetInputData(sphere);
		writer->SetFileName("Sphere.vtk");
		writer->Write();
		std::stringstream ss;
		ss << *normalsCentroids;
		qDebug() << ss.str().c_str();
		// Please use float precision
		float centroid[3]{ 0, 0, 0 };
		double dPoint[3];
		float fPoint[3];
		sphere->GetPoint(2, dPoint);
		std::copy(dPoint, dPoint + 3, fPoint);
		vtkMath::Add(centroid, fPoint, centroid);
		sphere->GetPoint(3, dPoint);
		std::copy(dPoint, dPoint + 3, fPoint);
		vtkMath::Add(centroid, fPoint, centroid);
		sphere->GetPoint(4, dPoint);
		std::copy(dPoint, dPoint + 3, fPoint);
		vtkMath::Add(centroid, fPoint, centroid);
		vtkMath::MultiplyScalar(centroid, 1.0f / 3.0f);
		vtkDataArray *sphereArray = sphere->GetPointData()->GetArray("Centroids");
		QCOMPARE((float)sphereArray->GetTuple3(0)[0] + 1, centroid[0] + 1);
		QCOMPARE((float)sphereArray->GetTuple3(0)[1] + 1, centroid[1] + 1);
		QCOMPARE((float)sphereArray->GetTuple3(0)[2] + 1, centroid[2] + 1);
	}
};

QTEST_APPLESS_MAIN(TestvtkPolyDataNormalsCentroids);
#include "TestvtkPolyDataNormalsCentroids.moc"