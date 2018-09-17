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
#include <vtkPolyDataNormalsCentroids.h>
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

void stepOfComputation_legacy(
	vtkImageData *image, 
	vtkPolyData *polyData, 
	const BET_Parameters &bp, 
	double smoothArg) {
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

	//vtkPoints *points = polyData->GetPoints();
	vtkFloatArray *points = vtkFloatArray::FastDownCast(polyData->GetPoints()->GetData());
	vtkIdType numPoints = polyData->GetNumberOfPoints();
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
		vtkPolygon::ComputeCentroid(polyData->GetPoints(), numOfPointsOfCell, pointsOfCell, centroid);
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
	float *fPolysCentroid = polyCentroid->WritePointer(0, 3 * numPolys);
	float *fPointsCentroid = pointsCentroid->WritePointer(0, 3 * numPoints);
	float *fCounterOfPointsCentroid = counterOfPointsCentroid->WritePointer(0, numPoints);
	float *fPoints = points->WritePointer(0, 3 * numPoints);
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
			fPointsCentroid[3 * pointsOfCell[i]] += fPolysCentroid[3 * polyId];
			fPointsCentroid[3 * pointsOfCell[i] + 1] += fPolysCentroid[3 * polyId + 1];
			fPointsCentroid[3 * pointsOfCell[i] + 2] += fPolysCentroid[3 * polyId + 2];
			fCounterOfPointsCentroid[pointsOfCell[i]]++;
		}
	}

	for (vtkIdType i = 0; i < numPoints; ++i) {
		fPointsCentroid[3 * i] /= fCounterOfPointsCentroid[i];
		fPointsCentroid[3 * i] -= (fPoints[3 * i] / 3);
		fPointsCentroid[3 * i] /= 2;
		fPointsCentroid[3 * i] *= 3;
		fPointsCentroid[3 * i + 1] /= fCounterOfPointsCentroid[i];
		fPointsCentroid[3 * i + 1] -= (fPoints[3 * i + 1]/ 3);
		fPointsCentroid[3 * i + 1] /= 2;
		fPointsCentroid[3 * i + 1] *= 3;
		fPointsCentroid[3 * i + 2] /= fCounterOfPointsCentroid[i];
		fPointsCentroid[3 * i + 2] -= (fPoints[3 * i + 2] / 3);
		fPointsCentroid[3 * i + 2] /= 2;
		fPointsCentroid[3 * i + 2] *= 3;
	}

	polyData->GetCellData()->AddArray(polyCentroid);
	polyData->GetPointData()->AddArray(pointsCentroid);

	vtkFloatArray *pointsNormal = vtkFloatArray::FastDownCast(polyData->GetPointData()->GetNormals());
	float *fPointsNormal = pointsNormal->WritePointer(0, 3 * numPoints);
	vtkSmartPointer<vtkFloatArray> s =
		vtkSmartPointer<vtkFloatArray>::New();
	s->Allocate(3 * numPoints);
	s->SetNumberOfComponents(3);
	s->SetNumberOfTuples(numPoints);
	s->SetName("s");
	float *fS = s->WritePointer(0, 3 * numPoints);
	vtkSmartPointer<vtkFloatArray> sn =
		vtkSmartPointer<vtkFloatArray>::New();
	sn->Allocate(3 * numPoints);
	sn->SetNumberOfComponents(3);
	sn->SetNumberOfTuples(numPoints);
	sn->SetName("sn");
	float *fSn = sn->WritePointer(0, 3 * numPoints);
	memcpy(fSn, fPointsNormal, sizeof(float) * 3 * numPoints);
	vtkSmartPointer<vtkFloatArray> st =
		vtkSmartPointer<vtkFloatArray>::New();
	st->Allocate(3 * numPoints);
	st->SetNumberOfComponents(3);
	st->SetNumberOfTuples(numPoints);
	st->SetName("st");
	float *fSt = st->WritePointer(0, 3 * numPoints);
	for (vtkIdType id = 0; id < numPoints; ++id) {
		vtkMath::Subtract(fPointsCentroid + 3 * id, fPoints + 3 * id, fS + 3 * id);
		float tmp = vtkMath::Dot(fS + 3 * id, fPointsNormal + 3 * id);
		vtkMath::MultiplyScalar(fSn + 3 * id, tmp);
		vtkMath::Subtract(fS + 3 * id, fSn + 3 * id, fSt + 3 * id);
	}
	constexpr float f1 = 0.5f;
	vtkSmartPointer<vtkFloatArray> u1 =
		vtkSmartPointer<vtkFloatArray>::New();
	u1->Allocate(3 * numPoints);
	u1->SetNumberOfComponents(3);
	u1->SetNumberOfTuples(numPoints);
	u1->SetName("u1");
	float *fU1 = u1->WritePointer(0, 3 * numPoints);
	memcpy(fU1, fSt, 3 * numPoints * sizeof(float));
	for (vtkIdType id = 0; id < numPolys; ++id) {
		fU1[id * 3] *= f1;
		fU1[id * 3 + 1] *= f1;
		fU1[id * 3 + 2] *= f1;
	}
	const double rmin = 3.33 * smoothArg;
	const double rmax = 10 * smoothArg;
	const double E = (1 / rmin + 1 / rmax) / 2.;
	const double F = 6. / (1 / rmin - 1 / rmax);

}

void stepOfComputation(
	vtkImageData *image,
	vtkPolyData *polyData,
	const BET_Parameters &bp,
	double smoothArg) {
	vtkSmartPointer<vtkPolyDataNormalsCentroids> normalsCentroid =
		vtkSmartPointer<vtkPolyDataNormalsCentroids>::New();
	normalsCentroid->SetInputData(polyData);
	normalsCentroid->SetComputePointNormals(true);
	normalsCentroid->SetComputePointCentroids(true);
	normalsCentroid->Update();

	polyData->ShallowCopy(normalsCentroid->GetOutput());
	const vtkIdType numPoints = polyData->GetNumberOfPoints();
	vtkFloatArray *points = vtkFloatArray::FastDownCast(polyData->GetPoints()->GetData());
	float *points_f = points->WritePointer(0, 3 * numPoints);
	vtkFloatArray *normals = vtkFloatArray::FastDownCast(polyData->GetPointData()->GetNormals());
	float *normals_f = normals->WritePointer(0, 3 * numPoints);
	vtkFloatArray *centroids = vtkFloatArray::FastDownCast(polyData->GetPointData()->GetArray("Centroids"));
	float *centroids_f = centroids->WritePointer(0, 3 * numPoints);
//////////////////////////////////////// s, st, sn ////////////////////////////////////////
	vtkSmartPointer<vtkFloatArray> s =
		vtkSmartPointer<vtkFloatArray>::New();
	s->Allocate(3 * numPoints);
	s->SetNumberOfComponents(3);
	s->SetNumberOfTuples(numPoints);
	s->SetName("s");
	float *s_f = s->WritePointer(0, 3 * numPoints);
	vtkSmartPointer<vtkFloatArray> sn =
		vtkSmartPointer<vtkFloatArray>::New();
	sn->Allocate(3 * numPoints);
	sn->SetNumberOfComponents(3);
	sn->SetNumberOfTuples(numPoints);
	sn->SetName("sn");
	float *sn_f = sn->WritePointer(0, 3 * numPoints);
	memcpy(sn_f, normals_f, sizeof(float) * 3 * numPoints);
	vtkSmartPointer<vtkFloatArray> st =
		vtkSmartPointer<vtkFloatArray>::New();
	st->Allocate(3 * numPoints);
	st->SetNumberOfComponents(3);
	st->SetNumberOfTuples(numPoints);
	st->SetName("st");
	float *st_f = st->WritePointer(0, 3 * numPoints);
	for (vtkIdType id = 0; id < numPoints; ++id) {
		vtkMath::Subtract(centroids_f + 3 * id, points_f+ 3 * id, s_f + 3 * id);
		float tmp = vtkMath::Dot(s_f + 3 * id, normals_f + 3 * id);
		vtkMath::MultiplyScalar(sn_f + 3 * id, tmp);
		vtkMath::Subtract(s_f + 3 * id, sn_f + 3 * id, st_f + 3 * id);
	}
//////////////////////////////////////// u1 ////////////////////////////////////////
	constexpr float f1 = 0.5f;
	vtkSmartPointer<vtkFloatArray> u1 =
		vtkSmartPointer<vtkFloatArray>::New();
	u1->Allocate(3 * numPoints);
	u1->SetNumberOfComponents(3);
	u1->SetNumberOfTuples(numPoints);
	u1->SetName("u1");
	float *u1_f = u1->WritePointer(0, 3 * numPoints);
	memcpy(u1_f, st_f, 3 * numPoints * sizeof(float));
	for (vtkIdType id = 0; id < numPoints; ++id) {
		vtkMath::MultiplyScalar(u1_f + id * 3, f1);
	}
//////////////////////////////////////// mean distance from vertex to neighboring ////////////////////////////////////////	
	float l = 0;
	vtkSmartPointer<vtkFloatArray> mediumDistanceOfNeighbour =
		vtkSmartPointer<vtkFloatArray>::New();
	mediumDistanceOfNeighbour->Allocate(numPoints);
	mediumDistanceOfNeighbour->SetNumberOfComponents(1);
	mediumDistanceOfNeighbour->SetNumberOfTuples(numPoints);
	mediumDistanceOfNeighbour->SetName("mediumDistanceOfNeighbour");
	float *mediumDistanceOfNeighbour_f =
		mediumDistanceOfNeighbour->WritePointer(0, numPoints);
	std::fill_n(
		mediumDistanceOfNeighbour_f,
		numPoints,
		0.0f);
	vtkSmartPointer<vtkIdTypeArray> counter =
		vtkSmartPointer<vtkIdTypeArray>::New();
	//counter->Allocate(numPoints);
	counter->SetNumberOfComponents(1);
	counter->SetNumberOfTuples(numPoints);
	vtkIdType *counter_id = counter->WritePointer(0, numPoints);
	polyData->BuildLinks();
	polyData->BuildCells();
	
	vtkCellArray *polys = polyData->GetPolys();
	polys->InitTraversal();
	for (vtkIdType cid = 0; cid < polys->GetNumberOfCells(); ++cid) {
		vtkIdType npts;
		vtkIdType *pointIds;
		polys->GetNextCell(npts, pointIds);
		for (vtkIdType pid0 = 0; pid0 < npts; ++pid0) {
			counter_id[pointIds[pid0]] += npts - 1;
			for (vtkIdType pid1 = 0; pid1 < npts; ++pid1) {
				// escape the same point
				if (pointIds[pid0] != pointIds[pid1]) {
					float diff[3];
					vtkMath::Subtract(points_f + 3 * pointIds[pid0], points_f + 3 * pointIds[pid1], diff);
					float norm = vtkMath::Norm(diff);
					mediumDistanceOfNeighbour_f[pointIds[pid0]] += norm;
				}
			}
		}
	}
	for (vtkIdType id = 0; id < numPoints; ++id) {
		mediumDistanceOfNeighbour_f[id] /= counter_id[id];
		l += mediumDistanceOfNeighbour_f[id];
	}
	l /= numPoints;
//////////////////////////////////////// u2 ////////////////////////////////////////
	const double rmin = 3.33 * smoothArg;
	const double rmax = 10 * smoothArg;
	const double E = (1 / rmin + 1 / rmax) / 2.;
	const double F = 6. / (1 / rmin - 1 / rmax);
	vtkSmartPointer<vtkFloatArray> u2 =
		vtkSmartPointer<vtkFloatArray>::New();
	u2->Allocate(3 * numPoints);
	u2->SetNumberOfComponents(3);
	u2->SetNumberOfTuples(numPoints);
	u2->SetName("u2");
	for (vtkIdType id = 0; id < numPoints; ++id) {
		// @todo
		// calculation in source code
		//  rinv = (2 * fabs(sn|n))/(l*l);
		// means: 
		// r = l^2 / (2 * sqqrt(norm(sn)))
		// calculation in paper
		// r = l^2 / (2 * norm(sn))
		// following source code for now
		double r_inv = (2 * fabs(vtkMath::Dot(sn_f + id * 3, normals_f + id * 3))) / l * l;

	}
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
	stepOfComputation(image, sphere, bp, 1.0);
	//auto polyDataWriter =
	//	vtkSmartPointer<vtkPolyDataWriter>::New();
	//polyDataWriter->SetInputData(sphere);
	//polyDataWriter->SetFileName("sphere.vtk");
	//polyDataWriter->Write();
}

int main(int argc, char **argv) {

	calculate();
	cin.get();
	return 0;
}
