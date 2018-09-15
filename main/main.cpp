//#include "bet2.h"
#include <vtkNIFTIImageReader.h>
#include <vtkNIFTIImageWriter.h>
#include "vtkBet2Filter.h"
#include <vtkSmartPointer.h>
#include <vtkPlatonicSolidSource.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointData.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkTransform.h>
#include <vtkImageHistogramStatistics.h>
#include <vtkImageData.h>
#include <vtkImageThreshold.h>
#include <vtkImageIterator.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPolyDataNormals.h>
#include <vtkFloatArray.h>
#include <vtkPolygon.h>
#include <vtkCellData.h>
#include <vtkCleanPolyData.h>
// std 
#include <array>
//using namespace bet2;
struct BET_Parameters{
	double min;
	double max;
	double t98;
	double t2;
	double t;
	double tm;
	double radius;
	typedef std::array<double, 3> CenterOfMass;
	CenterOfMass centerOfMass;
};

std::ostream& operator<<(std::ostream &os, BET_Parameters &bp)
{
	os << '\n' << "BET Parameters" << '\n';
	os << "min " << bp.min << '\n';
	os << "max " << bp.max << '\n';
	os << "t98 " << bp.t98 << '\n';
	os << "t2 " << bp.t2 << '\n';
	os << "t " << bp.t << '\n';
	os << "tm " << bp.tm << '\n';
	os << "radius " << bp.radius << '\n';
	os << "centerOfMass " << bp.centerOfMass[0] << ' ' << bp.centerOfMass[1] << ' '<< bp.centerOfMass[2] << ' ' << '\n';
	return os;
}


void generateSphere(size_t n, vtkPolyData *data)
{
	// create Icosahedron, whose number of faces is twenty
	vtkSmartPointer<vtkPlatonicSolidSource> source =
		vtkSmartPointer<vtkPlatonicSolidSource>::New();
	source->SetSolidTypeToIcosahedron();
	source->Update();
	vtkSmartPointer<vtkPolyData> sphere = source->GetOutput();
	// iterate to re-tesselates.
	for (size_t i = 0; i < n; ++i) {
		vtkSmartPointer<vtkLinearSubdivisionFilter> subdivision = 
			vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
		subdivision->SetInputConnection(source->GetOutputPort());
		subdivision->SetNumberOfSubdivisions(2);
		subdivision->Update();
		sphere = subdivision->GetOutput();
		
		for (vtkIdType id = 0; id < sphere->GetNumberOfPoints(); ++id) {
			double points[3];
			sphere->GetPoint(id, points);
			double norm = vtkMath::Norm(points);
			vtkMath::MultiplyScalar(points, 1/norm);
			sphere->GetPoints()->SetPoint(id, points);
		}
	}
	data->ShallowCopy(sphere);
}

BET_Parameters adjustInitialPolyData(vtkImageData *image, vtkPolyData *polyData)
{
	BET_Parameters bp;

	vtkSmartPointer<vtkImageHistogramStatistics> imageHistogram =
		vtkSmartPointer<vtkImageHistogramStatistics>::New();
	imageHistogram->SetInputData(image);
	imageHistogram->SetAutoRangePercentiles(2, 98);
	imageHistogram->SetGenerateHistogramImage(false);
	imageHistogram->Update();
	imageHistogram->GetAutoRange(bp.t2, bp.t98);
	bp.max = imageHistogram->GetMaximum();
	bp.min = imageHistogram->GetMinimum();
	bp.t = bp.t2 + .1*(bp.t98 - bp.t2);

	const double *spacing = image->GetSpacing();
	double counter = 0;
	size_t numIsBrain = 0;
	double centerOfMass[3]{ 0, 0, 0 };
	for (vtkIdType i = image->GetExtent()[0]; i <= image->GetExtent()[1]; ++i) {
		for (vtkIdType j = image->GetExtent()[2]; j <= image->GetExtent()[3]; ++j) {
			for (vtkIdType k = image->GetExtent()[4]; k <= image->GetExtent()[5]; ++k) {
				double voxel = image->GetScalarComponentAsDouble(i, j, k, 0);
				if (voxel > bp.t) {
					voxel = vtkMath::Min(voxel, bp.t98) - bp.t2;
					counter += voxel;
					centerOfMass[0] += voxel * i * spacing[0];
					centerOfMass[1] += voxel * j * spacing[1];
					centerOfMass[2] += voxel * k * spacing[2];
					numIsBrain++;
				}
			}
		}
	}
	vtkMath::MultiplyScalar(centerOfMass, 1 / counter);
	std::copy(centerOfMass, centerOfMass + 3, bp.centerOfMass.begin());
	// sphere volume equation: 
	// volume = 4/3 * PI * radius^3
	double brainVolume = numIsBrain * spacing[0] * spacing[1] * spacing[2];
	bp.radius = pow(brainVolume * 0.75 / vtkMath::Pi(), 1.0 / 3.0);

	vtkSmartPointer<vtkTransform> transform =
		vtkSmartPointer<vtkTransform>::New();
	transform->Identity();
	transform->Translate(centerOfMass);
	transform->Scale(bp.radius, bp.radius, bp.radius);
	vtkSmartPointer<vtkTransformPolyDataFilter> transformPolyData =
		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformPolyData->SetInputData(polyData);
	transformPolyData->SetTransform(transform);
	transformPolyData->Update();
	polyData->ShallowCopy(transformPolyData->GetOutput());

	vtkSmartPointer<vtkImageThreshold> imageThreshold =
		vtkSmartPointer<vtkImageThreshold>::New();
	imageThreshold->SetInputData(image);
	imageThreshold->ThresholdBetween(bp.t2, bp.t98);
	imageThreshold->SetOutValue(0);
	
	vtkSmartPointer<vtkPolyDataToImageStencil> polyDataToImageStencil =
		vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	polyDataToImageStencil->SetInputData(polyData);
	polyDataToImageStencil->SetOutputOrigin(image->GetOrigin());
	polyDataToImageStencil->SetOutputSpacing(image->GetSpacing());
	polyDataToImageStencil->SetOutputWholeExtent(image->GetExtent());
	polyDataToImageStencil->Update();

	vtkSmartPointer<vtkImageStencil> imageStencil =
		vtkSmartPointer<vtkImageStencil>::New();
	imageStencil->SetInputConnection(imageThreshold->GetOutputPort());
	imageStencil->SetStencilConnection(polyDataToImageStencil->GetOutputPort());
	imageStencil->ReverseStencilOff();
	imageStencil->SetBackgroundValue(0.0);
	imageStencil->Update();

	vtkSmartPointer<vtkImageHistogramStatistics> imageHistogram2 =
		vtkSmartPointer<vtkImageHistogramStatistics>::New();
	imageHistogram2->SetInputConnection(imageStencil->GetOutputPort());
	imageHistogram2->SetStencilConnection(polyDataToImageStencil->GetOutputPort());
	imageHistogram2->SetGenerateHistogramImage(false);
	imageHistogram2->Update();
	bp.tm = imageHistogram2->GetMedian();
	return bp;
}

void stepOfComputation(vtkImageData *image, vtkPolyData *polyData, const BET_Parameters &bp) {
	vtkSmartPointer<vtkPolyDataNormals> polyDataNormals =
		vtkSmartPointer<vtkPolyDataNormals>::New();
	polyDataNormals->SetInputData(polyData);
	polyDataNormals->SetComputePointNormals(true);
	polyDataNormals->SetComputeCellNormals(true);
	polyDataNormals->SetNonManifoldTraversal(false);
	polyDataNormals->SetConsistency(true);
	polyDataNormals->Update();
	polyData->ShallowCopy(polyDataNormals->GetOutput());

	vtkSmartPointer<vtkCleanPolyData> cleanPolyData =
		vtkSmartPointer<vtkCleanPolyData>::New();
	cleanPolyData->SetInputData(polyData);
	cleanPolyData->SetPointMerging(true);
	cleanPolyData->Update();
	polyData->ShallowCopy(cleanPolyData->GetOutput());

	vtkPoints *points = polyData->GetPoints();
	vtkIdType numPoints = points->GetNumberOfPoints();
	vtkCellArray *polys = polyData->GetPolys();
	vtkIdType numPolys = polys->GetNumberOfCells();
	vtkSmartPointer<vtkFloatArray> polyCentroid =
		vtkSmartPointer<vtkFloatArray>::New();
	polyCentroid->Allocate(3 * numPolys);
	polyCentroid->SetNumberOfComponents(3);
	polyCentroid->SetNumberOfTuples(numPolys);
	polyCentroid->SetName("Centroids");

	polys->InitTraversal();
	for (vtkIdType polyId = 0; polyId < numPolys; ++polyId) {
		vtkIdType numOfPointsOfCell;
		vtkIdType *pointsOfCell = nullptr;
		if (!polys->GetNextCell(numOfPointsOfCell, pointsOfCell)) {
			break;
		}
		double centroid[3];
		vtkPolygon::ComputeCentroid(points, numOfPointsOfCell, pointsOfCell, centroid);
		polyCentroid->SetTuple(polyId, centroid);
	}
	
	vtkSmartPointer<vtkFloatArray> pointsCentroid =
		vtkSmartPointer<vtkFloatArray>::New();
	pointsCentroid->Allocate(3 * numPoints);
	pointsCentroid->SetNumberOfComponents(3);
	pointsCentroid->SetNumberOfTuples(numPoints);
	pointsCentroid->SetName("Centroids");
	
	vtkSmartPointer<vtkFloatArray> counterOfPointsCentroid =
		vtkSmartPointer<vtkFloatArray>::New();
	counterOfPointsCentroid->Allocate(numPoints);
	counterOfPointsCentroid->SetNumberOfComponents(1);
	counterOfPointsCentroid->SetNumberOfTuples(numPoints);
	polys->InitTraversal();
	float *fPolyCentroid = polyCentroid->WritePointer(0, 3 * numPolys);
	float *fPointsCentroid = pointsCentroid->WritePointer(0, 3 * numPoints);
	float *fCounterOfPointsCentroid = counterOfPointsCentroid->WritePointer(0, numPoints);
	std::fill(fPointsCentroid, fPointsCentroid + 3 * numPoints, 0.0f);
	std::fill(fCounterOfPointsCentroid, fCounterOfPointsCentroid + numPoints, 0.0f);
	for (vtkIdType polyId = 0; polyId < numPolys; ++polyId) {
		vtkIdType numOfPointsOfCell;
		vtkIdType *pointsOfCell = nullptr;
		if (!polys->GetNextCell(numOfPointsOfCell, pointsOfCell)) {
			break;
		}
		for (vtkIdType i = 0; i < numOfPointsOfCell; ++i)
		{
			fPointsCentroid[3 * pointsOfCell[i]] += fPolyCentroid[3 * polyId];
			fPointsCentroid[3 * pointsOfCell[i] + 1] += fPolyCentroid[3 * polyId + 1];
			fPointsCentroid[3 * pointsOfCell[i] + 2] += fPolyCentroid[3 * polyId + 2];
			fCounterOfPointsCentroid[pointsOfCell[i]]++;
		}
	}

	for (vtkIdType i = 0; i < numPoints; ++i) {
		fPointsCentroid[3 * i] /= fCounterOfPointsCentroid[i];
		fPointsCentroid[3 * i] -= (points->GetPoint(i)[0] / 3);
		fPointsCentroid[3 * i] /= 2;
		fPointsCentroid[3 * i] *= 3;
		fPointsCentroid[3 * i + 1] /= fCounterOfPointsCentroid[i];
		fPointsCentroid[3 * i + 1] -= (points->GetPoint(i)[1] / 3);
		fPointsCentroid[3 * i + 1] /= 2;
		fPointsCentroid[3 * i + 1] *= 3;
		fPointsCentroid[3 * i + 2] /= fCounterOfPointsCentroid[i];
		fPointsCentroid[3 * i + 2] -= (points->GetPoint(i)[2] / 3);
		fPointsCentroid[3 * i + 2] /= 2;
		fPointsCentroid[3 * i + 2] *= 3;
	}

	polyData->GetCellData()->AddArray(polyCentroid);
	polyData->GetPointData()->AddArray(pointsCentroid);

}
void calculate()
{
	auto reader =
		vtkSmartPointer<vtkNIFTIImageReader>::New();
	reader->SetFileName("T2.nii");
	reader->Update();

	auto sphere =
		vtkSmartPointer<vtkPolyData>::New();
	auto image =
		vtkSmartPointer<vtkImageData>::New();
	image->ShallowCopy(reader->GetOutput());
	generateSphere(0, sphere);
	auto bp = adjustInitialPolyData(image, sphere);
	std::cout << bp << '\n';
	stepOfComputation(image, sphere, bp);
	auto polyDataWriter =
		vtkSmartPointer<vtkPolyDataWriter>::New();
	polyDataWriter->SetInputData(sphere);
	polyDataWriter->SetFileName("sphere.vtk");
	polyDataWriter->Write();
}
int main(int argc, char **argv) {

	calculate();
	cin.get();
	return 0;
}
