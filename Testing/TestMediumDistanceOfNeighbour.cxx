#include <QTest>
// vtk
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkStripper.h>
#include <vtkTriangleFilter.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
// std
#include <iostream>
class TestMediumDistanceOfNeighbour : public QObject
{
	Q_OBJECT;
private Q_SLOTS:
	void initTestCase(){}
	void cleanupTestCase()
	{
		//std::cin.get();
	}
	void init() {}
	void cleanup() {}
	void test(){
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
		vtkFloatArray *points = vtkFloatArray::FastDownCast(sphere->GetPoints()->GetData());
		float *fPoints = points->WritePointer(0, 3 * sphere->GetNumberOfPoints());
		vtkSmartPointer<vtkFloatArray> mediumDistanceOfNeighbour =
			vtkSmartPointer<vtkFloatArray>::New();
		mediumDistanceOfNeighbour->Allocate(sphere->GetNumberOfPoints());
		mediumDistanceOfNeighbour->SetNumberOfComponents(1);
		mediumDistanceOfNeighbour->SetNumberOfTuples(sphere->GetNumberOfPoints());
		mediumDistanceOfNeighbour->SetName("mediumDistanceOfNeighbour");
		float *fMediumDistanceOfNeighbour = 
			mediumDistanceOfNeighbour->WritePointer(0, sphere->GetNumberOfPoints());
		std::fill_n(
			fMediumDistanceOfNeighbour,
			sphere->GetNumberOfPoints(),
			0.0f);
		vtkIdType *counter = new vtkIdType[sphere->GetNumberOfPoints()]();
		sphere->BuildLinks();
		sphere->BuildCells();
		//vtkPoints *points = sphere->GetPoints();
		vtkCellArray *polys = sphere->GetPolys();
		polys->InitTraversal();
		for (vtkIdType cid = 0; cid < polys->GetNumberOfCells(); ++cid) {
			vtkIdType npts;
			vtkIdType *pointIds;
			polys->GetNextCell(npts, pointIds);
			for (vtkIdType pid0 = 0; pid0 < npts; ++pid0) {
				counter[pointIds[pid0]] += npts - 1;
				for (vtkIdType pid1 = 0; pid1 < npts; ++pid1) {
					if (pointIds[pid0] != pointIds[pid1]) {
						float diff[3];
						vtkMath::Subtract(fPoints + 3 * pointIds[pid0], fPoints+ 3 * pointIds[pid1], diff);
						float norm = vtkMath::Norm(diff);
						fMediumDistanceOfNeighbour[pointIds[pid0]] += norm;
					}
				}
			}
		}
		for (vtkIdType id = 0; id < sphere->GetNumberOfPoints(); ++id) {
			fMediumDistanceOfNeighbour[id] /= counter[id];
		}
		delete[]counter;
		sphere->GetPointData()->AddArray(mediumDistanceOfNeighbour);
		float medium = 0;
		double point[3];
		sphere->GetPoint(0, point);
		vtkMath::Subtract(point, sphere->GetPoint(2), point);
		medium += vtkMath::Norm(point);
		sphere->GetPoint(0, point);
		vtkMath::Subtract(point, sphere->GetPoint(3), point);
		medium += vtkMath::Norm(point);
		sphere->GetPoint(0, point);
		vtkMath::Subtract(point, sphere->GetPoint(4), point);
		medium += vtkMath::Norm(point);
		medium /= 3;
		QCOMPARE(fMediumDistanceOfNeighbour[0], medium);
		//vtkSmartPointer<vtkPolyDataWriter> writer =
		//	vtkSmartPointer<vtkPolyDataWriter>::New();
		//writer->SetInputData(sphere);
		//writer->SetFileName("Sphere.vtk");
		//writer->Write();
		//for (vtkIdType pid = 0; pid < points->GetNumberOfPoints(); ++pid) {
		//	unsigned short ncells = 0;
		//	vtkIdType *cellIDs;
		//	sphere->GetPointCells(pid, ncells, cellIDs);
		//	for (vtkIdType cid = 0; cid < ncells; ++cid) {
		//		vtkCell *cell = sphere->GetCell(cellIDs[cid]);
		//		vtkIdList *pointIDs = cell->GetPointIds();
		//		for (vtkIdType id = 0; id< pointIDs->GetNumberOfIds(); ++id) {
		//			
		//		}
		//	}
		//}
	}
};

QTEST_APPLESS_MAIN(TestMediumDistanceOfNeighbour);
#include "TestMediumDistanceOfNeighbour.moc"