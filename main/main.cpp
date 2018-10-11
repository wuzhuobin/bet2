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
#include <vtkResliceImageViewer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
// std 
#include <array>
//using namespace bet2;
/**
 * @struct	BET_Parameters
 * @brief	Brain Extraction Tool Parameters
 * @author	WUZHUOBIN
 */
struct BET_Parameters{
	double min;					///< The minimum intensity of the image.
	double max;					///< The maximum intensity of the image.
	double t98;					///< Caculated by looking at the intensity histogram, intensity above 98%.
	double t2;					///< Caculated by looking at the intensity histogram, intensity below 2%.
	double t;					///< Attemptes to distinguish between brain matter and background, 10% between #t2 and #t98.
	double tm;					///< The median intensity of all points with a sphere of the estimated #centerOfMass and #radius.
	double radius;				///< The estimated raduis of the sphere.
	typedef std::array<double, 3> CenterOfMass;
	CenterOfMass centerOfMass;	///< The estimated center of the sphere.
};
/**
 * @fn				std::ostream& operator<<(std::ostream &os, BET_Parameters &bp)	
 * @brief			Serializing output of BET_Parameters
 * @param[in]		os The std::ostream.
 * @param[in]		bp The BET_Parameters need to be print. 
 * @return			The input std::ostream
 */
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

void renderPolyData(vtkPolyData *polyData, vtkResliceImageViewer *viewer) {
	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
	actor->GetMapper()->SetInputDataObject(polyData);
	actor->GetMapper()->SetScalarVisibility(false);
	viewer->GetRenderer()->AddActor(actor);
	viewer->Render();
}

void writePolyData(vtkPolyData *data, const char *name)
{
	vtkSmartPointer<vtkPolyDataWriter> writer =
		vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInputData(data);
	writer->SetFileName(name);
	writer->Write();
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
		subdivision->SetInputData(sphere);
		subdivision->SetNumberOfSubdivisions(1);
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
	// calculate bp.max, bp.min, bp.t, bp.t2, bp.t98;
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
	// calculate bp.centerOfMass
	const double *spacing = image->GetSpacing();
	double counter = 0;
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
				}
			}
		}
	}
	vtkMath::MultiplyScalar(centerOfMass, 1 / counter);
	std::copy(centerOfMass, centerOfMass + 3, bp.centerOfMass.begin());
	// calculate bp.radius
	size_t numIsBrain = 0;
	for (vtkIdType i = image->GetExtent()[0]; i <= image->GetExtent()[1]; ++i) {
		for (vtkIdType j = image->GetExtent()[2]; j <= image->GetExtent()[3]; ++j) {
			for (vtkIdType k = image->GetExtent()[4]; k <= image->GetExtent()[5]; ++k) {
				double voxel = image->GetScalarComponentAsDouble(i, j, k, 0);
				if (voxel > bp.t) {
					++numIsBrain;
				}
			}
		}
	}
	// using the brain's volume as the sphere's volume
	// basing on sphere volume equation, find radius.
	// sphere volume equation: 
	// volume = 4/3 * PI * radius^3
	double brainVolume = numIsBrain * spacing[0] * spacing[1] * spacing[2];
	bp.radius = pow(brainVolume * 0.75 / vtkMath::Pi(), 1.0 / 3.0);
	// transform
	vtkSmartPointer<vtkTransform> transform =
		vtkSmartPointer<vtkTransform>::New();
	transform->Identity();
	transform->Translate(centerOfMass);
	transform->Scale(bp.radius * 0.5, bp.radius * 0.5, bp.radius * 0.5);
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

void stepOfComputation(
	vtkImageData *image,
	vtkPolyData *polyData,
	const int pass,
	const double increaseSmoothing, 
	const int iteraction_number, 
	double &l,
	const BET_Parameters &bp,
	double smoothArg) {
//////////////////////////////////////// normals and centroids////////////////////////////////////////
	//std::cerr  << "normals and centroids\n";
	vtkSmartPointer<vtkPolyDataNormalsCentroids> normalsCentroid =
		vtkSmartPointer<vtkPolyDataNormalsCentroids>::New();
	normalsCentroid->SetInputData(polyData);
	normalsCentroid->SetComputePointNormals(true);
	normalsCentroid->SetComputePointCentroids(true);
	normalsCentroid->SetSplitting(false);
	normalsCentroid->SetConsistency(true);
	normalsCentroid->Update();

	//vtkSmartPointer<vtkCleanPolyData> cleanPolyData =
	//	vtkSmartPointer<vtkCleanPolyData>::New();
	//cleanPolyData->SetInputConnection(normalsCentroid->GetOutputPort());
	//cleanPolyData->SetPointMerging(true);
	//cleanPolyData->Update();

	polyData->ShallowCopy(normalsCentroid->GetOutput());
	const vtkIdType numPoints = polyData->GetNumberOfPoints();
	vtkFloatArray *points = vtkFloatArray::FastDownCast(polyData->GetPoints()->GetData());
	float *points_f = points->WritePointer(0, 3 * numPoints);
	vtkFloatArray *normals = vtkFloatArray::FastDownCast(polyData->GetPointData()->GetNormals());
	float *normals_f = normals->WritePointer(0, 3 * numPoints);
	for (vtkIdType id = 0; id < numPoints; ++id) {
		vtkMath::MultiplyScalar(normals_f + 3 * id, -1.0f);
	}
	vtkFloatArray *centroids = vtkFloatArray::FastDownCast(polyData->GetPointData()->GetArray("Centroids"));
	float *centroids_f = centroids->WritePointer(0, 3 * numPoints);
//////////////////////////////////////// s, st, sn ////////////////////////////////////////
	//std::cerr  << "s, st, sn\n";
	vtkSmartPointer<vtkFloatArray> s =
		vtkSmartPointer<vtkFloatArray>::New();
	s->Allocate(3 * numPoints);
	s->SetNumberOfComponents(3);
	s->SetNumberOfTuples(numPoints);
	s->SetName("s");
	polyData->GetPointData()->AddArray(s);
	float *s_f = s->WritePointer(0, 3 * numPoints);
	vtkSmartPointer<vtkFloatArray> sn =
		vtkSmartPointer<vtkFloatArray>::New();
	sn->Allocate(3 * numPoints);
	sn->SetNumberOfComponents(3);
	sn->SetNumberOfTuples(numPoints);
	sn->SetName("sn");
	polyData->GetPointData()->AddArray(sn);
	float *sn_f = sn->WritePointer(0, 3 * numPoints);
	memcpy(sn_f, normals_f, sizeof(float) * 3 * numPoints);
	vtkSmartPointer<vtkFloatArray> st =
		vtkSmartPointer<vtkFloatArray>::New();
	st->Allocate(3 * numPoints);
	st->SetNumberOfComponents(3);
	st->SetNumberOfTuples(numPoints);
	st->SetName("st");
	polyData->GetPointData()->AddArray(st);
	float *st_f = st->WritePointer(0, 3 * numPoints);
	for (vtkIdType id = 0; id < numPoints; ++id) {
		vtkMath::Subtract(centroids_f + 3 * id, points_f+ 3 * id, s_f + 3 * id);
		float tmp = vtkMath::Dot(s_f + 3 * id, normals_f + 3 * id);
		vtkMath::MultiplyScalar(sn_f + 3 * id, tmp);
		vtkMath::Subtract(s_f + 3 * id, sn_f + 3 * id, st_f + 3 * id);
	}
//////////////////////////////////////// u1 ////////////////////////////////////////
	//std::cerr  << "u1\n";
	constexpr float f1 = 0.5f;
	vtkSmartPointer<vtkFloatArray> u1 =
		vtkSmartPointer<vtkFloatArray>::New();
	u1->Allocate(3 * numPoints);
	u1->SetNumberOfComponents(3);
	u1->SetNumberOfTuples(numPoints);
	u1->SetName("u1");
	polyData->GetPointData()->AddArray(u1);
	float *u1_f = u1->WritePointer(0, 3 * numPoints);
	memcpy(u1_f, st_f, 3 * numPoints * sizeof(float));
	for (vtkIdType id = 0; id < numPoints; ++id) {
		vtkMath::MultiplyScalar(u1_f + id * 3, f1);
	}
//////////////////////////////////////// mean distance from vertex to neighboring ////////////////////////////////////////	
	//std::cerr  << "mean distance from vertext to neighboring\n";
	if (true) {
	//if (iteraction_number == 50 || iteraction_number % 100 == 0) {
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
		counter->Allocate(numPoints);
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
				if (npts != 3) {
					cout << "counter_id: " << counter_id[pointIds[pid0]] << '\n';
					cin.get();
				}
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
		cout << "l: " << l << '\n';
		//cin.get();
	}
	l = 3.34;
	//////////////////////////////////////// u2 ////////////////////////////////////////
		//std::cerr  << "u2\n";
	//const int pass = 0;
	const double increase_smoothing = 0;
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
	polyData->GetPointData()->AddArray(u2);
	float *u2_f = u2->WritePointer(0, 3 * numPoints);
	memcpy(u2_f, sn_f, 3 * numPoints * sizeof(float));
	for (vtkIdType id = 0; id < numPoints; ++id) {
		// @todo
		// calculation in source code
		//  rinv = (2 * fabs(sn|n))/(l*l);
		// means: 
		// r = l^2 / (2 * sqqrt(norm(sn)))
		// calculation in paper
		// r = l^2 / (2 * norm(sn))
		// paper's way is the following: 
		//float r_inv = (2 * vtkMath::Norm(sn_f + 3 * numPoints)) / (l * l);
		// source code is the following: 
		float r_inv = (2 * fabs(vtkMath::Dot(sn_f + id * 3, normals_f + id * 3))) / (l * l);
		float f2 = (1 + tanh(F*(r_inv - E)))*0.5;
		if (pass > 0) {
			if (vtkMath::Dot(s_f + 3 * id, normals_f + 3 * id) > 0) {
				f2 *= increase_smoothing;
				f2 = vtkMath::Min(f2, 1.0f);
			}
		}
		vtkMath::MultiplyScalar(u2_f + 3 * id, f2);
	}
	//////////////////////////////////////// u3 ////////////////////////////////////////
		//std::cerr  << "u3\n";
	float bt = powf(0.5f, 0.275f);
	//if (bt != 0.0)
	//{
	//	bt = Min(1., Max(0., bet_main_parameter + local_th*((*i)->get_coord().Z - zcog) / radius));
	//}
	const int d1 = 7;
	const int d2 = 3;
	const float normal_max_update_fraction = 0.5f;
	const float lambda_fit = 0.1;
	const double *bounds = image->GetBounds();
	const double *origin = image->GetOrigin();
	const double *spacing = image->GetSpacing();
	const int *extent = image->GetExtent();
	vtkSmartPointer<vtkFloatArray> u3 =
		vtkSmartPointer<vtkFloatArray>::New();
	u3->Allocate(3 * numPoints);
	u3->SetNumberOfComponents(3);
	u3->SetNumberOfTuples(numPoints);
	u3->SetName("u3");
	polyData->GetPointData()->AddArray(u3);
	float *u3_f = u3->WritePointer(0, 3 * numPoints);
	memcpy(u3_f, normals_f, 3 * numPoints * sizeof(float));
	double dscale = vtkMath::Min(vtkMath::Min(vtkMath::Min(1.0, spacing[0]), spacing[1]), spacing[2]);
	for (vtkIdType id = 0; id < numPoints; ++id) {
		float Imin = bp.tm;
		float Imax = bp.t;
		float f3 = 0;
		float p[3];
		vtkMath::Subtract(points_f + 3 * id, normals_f + 3 * id, p);
		float iv = (p[0] - origin[0]) / spacing[0] + 0.5;
		float jv = (p[1] - origin[1]) / spacing[1] + 0.5;
		float kv = (p[2] - origin[2]) / spacing[2] + 0.5;
		if (extent[0] <= (int)iv && (int)iv <= extent[1] &&
			extent[2] <= (int)jv && (int)jv <= extent[3] &&
			extent[4] <= (int)kv && (int)kv <= extent[5]) {
			float im = image->GetScalarComponentAsFloat(iv, jv, kv, 0);
			Imin = vtkMath::Min(Imin, im);
			Imax = vtkMath::Max(Imax, im);
			float nxv = normals_f[id * 3 + 0] / spacing[0];
			float nyv = normals_f[id * 3 + 1] / spacing[1];
			float nzv = normals_f[id * 3 + 2] / spacing[2];
			int i2 = iv - (d1 - 1)*nxv;
			int j2 = jv - (d1 - 1)*nyv;
			int k2 = kv - (d1 - 1)*nzv;
			if (extent[0] <= i2 && i2 <= extent[1] &&
				extent[2] <= j2 && j2 <= extent[3] &&
				extent[4] <= k2 && k2 <= extent[5]) {
				im = image->GetScalarComponentAsFloat(i2, j2, k2, 0);
				Imin = vtkMath::Min(Imin, im);
				nxv *= dscale;
				nyv *= dscale;
				nzv *= dscale;
				for (double gi = 2.0; gi < d1; gi += dscale)
				{
					//cout << gi << " " << endl;
					// the following is a quick calc of Pt p = (*i)->get_coord() + (-gi)*n;
					iv -= nxv; jv -= nyv; kv -= nzv;
					im = image->GetScalarComponentAsFloat(iv, jv, kv, 0);
					Imin = vtkMath::Min(Imin, im);

					if (gi < d2) {
						Imax = vtkMath::Max(Imax, im);
					}
				}

				Imin = vtkMath::Max((float)bp.t2, Imin);
				Imax = vtkMath::Min((float)bp.tm, Imax);

				const float tl = (Imax - bp.t2) * bt + bp.t2;
				if (Imax - bp.t2 > 0) {
					f3 = 2 * (Imin - tl) / (Imax - bp.t2);
				}
				else {
					f3 = (Imin - tl) * 2;
				}
			}
		}
		f3 *= (normal_max_update_fraction * lambda_fit * l);
		//cout << "f3: " << f3 << '\n';
		vtkMath::MultiplyScalar(u3_f + 3 * id, f3);
		//if (bounds[0] < p[0] && p[0] < bounds[1] &&
		//	bounds[2] < p[1] && p[1] < bounds[3] &&
		//	bounds[4] < p[2] && p[2] < bounds[5]) {

		//	float pImageCoordinate[3];
		//	memcpy(pImageCoordinate, p, sizeof(pImageCoordinate));
		//	pImageCoordinate[0] -= origin[0];
		//	pImageCoordinate[0] /= spacing[0];
		//	pImageCoordinate[1] -= origin[1];
		//	pImageCoordinate[1] /= spacing[1];
		//	pImageCoordinate[2] -= origin[2];
		//	pImageCoordinate[2] /= spacing[2];
		//	float pixelValue1 = image->GetScalarComponentAsFloat(
		//		pImageCoordinate[0], pImageCoordinate[1], pImageCoordinate[2], 0);
		//	Imin = vtkMath::Min(Imin, pixelValue1);
		//	Imax = vtkMath::Max(Imin, pixelValue1);
		//}
		//double local_t = bet_main_parameter;
		//if (local_th != 0.0)
		//  {
		//    local_t = Min(1., Max(0., bet_main_parameter + local_th*((*i)->get_coord().Z - zcog)/radius));
		//  }

	}
//////////////////////////////////////// u ////////////////////////////////////////
	//std::cerr  << "u\n";
	vtkSmartPointer<vtkFloatArray> u =
		vtkSmartPointer<vtkFloatArray>::New();
	u->Allocate(3 * numPoints);
	u->SetNumberOfComponents(3);
	u->SetNumberOfTuples(numPoints);
	u->SetName("u");
	polyData->GetPointData()->AddArray(u);
	float *u_f = u->WritePointer(0, 3 * numPoints);
	//vtkSmartPointer<vtkPoints> newPoints =
	//	vtkSmartPointer<vtkPoints>::New();
	//newPoints->SetNumberOfPoints(numPoints);
	//float *newPoints_f = vtkFloatArray::FastDownCast(newPoints->GetData())->WritePointer(0, 3 * numPoints);
	//memcpy(newPoints_f, points_f, 3 * numPoints * sizeof(float));
	for (vtkIdType id = 0; id < numPoints; ++id) {
		vtkMath::Add(u1_f + 3 * id, u2_f + 3 * id, u_f + 3 * id);
		vtkMath::Add(u_f + 3 * id, u3_f + 3 * id, u_f + 3 * id);
		if (id == 550) {
			cout << "l: " << l << '\n';
			cout << "u1: " << '[' << u1_f[3 * id + 0] << ',' <<
				u1_f[3 * id + 1] << ',' << u1_f[3 * id + 2] << ',' <<
				']' << '\n';
			cout << "u2: " << '[' << u2_f[3 * id + 0] << ',' <<
				u2_f[3 * id + 1] << ',' << u2_f[3 * id + 2] << ',' <<
				']' << '\n';
			cout << "u3: " << '[' << u3_f[3 * id + 0] << ',' <<
				u3_f[3 * id + 1] << ',' << u3_f[3 * id + 2] << ',' <<
				']' << '\n';
			//if (vtkMath::IsNan(u2_f[3 * id + 0])) {
			//	cin.get();
			//}
		}
		//vtkMath::Add(newPoints_f + 3 * id, u_f + 3 * id, newPoints_f + 3 * id);
		//polyData->GetPoints()->SetPoint(id, point);
		//vtkMath::Add(points_f + 3 * id, u_f + 3 * id, newPoints_f + 3 * id);
		//if (vtkMath::Norm(u_f + 3 * id) > 5) {
		//	cout << "id: " << id << '\n';
		//	cout << "norm: " << vtkMath::Norm(u_f + 3 * id) << '\n';
		//	cin.get();
		//}
		vtkMath::Add(points_f + 3 * id, u_f + 3 * id, points_f + 3 * id);
	}
	polyData->GetPointData()->AddArray(u);
	//polyData->SetPoints(newPoints);
	//polyData->Modified();
}

void calculate()
{
	auto reader =
		vtkSmartPointer<vtkNIFTIImageReader>::New();
	reader->SetFileName("T2.nii");
	reader->Update();
	auto polyDataWriter =
		vtkSmartPointer<vtkPolyDataWriter>::New();
	auto sphere =
		vtkSmartPointer<vtkPolyData>::New();
	auto image =
		vtkSmartPointer<vtkImageData>::New();
	image->ShallowCopy(reader->GetOutput());
	generateSphere(5, sphere);
	auto bp = adjustInitialPolyData(image, sphere);
	std::cout << bp << '\n';
	polyDataWriter->SetInputData(sphere);
	polyDataWriter->SetFileName("sphere0.vtk");
	polyDataWriter->Write();
	for(size_t i = 0; i < 1000; ++i){
		//stepOfComputation(image, sphere, bp, 1.0);
	}
	polyDataWriter->SetInputData(sphere);
	polyDataWriter->SetFileName("sphere1.vtk");
	polyDataWriter->Write();
}



int main(int argc, char **argv) {

	//calculate();
	vtkSmartPointer<vtkImageData> image;
	vtkSmartPointer<vtkNIFTIImageReader> imageReader =
		vtkSmartPointer<vtkNIFTIImageReader>::New();
	imageReader->SetFileName("T2.nii");
	imageReader->Update();
	image = imageReader->GetOutput();

	vtkSmartPointer<vtkPolyData> sphere;
	vtkSmartPointer<vtkPolyDataReader> polyDataReader =
		vtkSmartPointer<vtkPolyDataReader>::New();
	polyDataReader->SetFileName("sphere.vtk");
	polyDataReader->Update();
	sphere = polyDataReader->GetOutput();
	//BET_Parameters bp{ 0, 864, 408, 0, 40.8, 220, 88.3824, {115.153, 111.177, 77.4181} };
	BET_Parameters bp{ 0, 864, 372.384, 0, 37.2384, 192, 88.3824, {115.153, 111.177, 77.4181} };

	cout << bp;
	double l = 0;
	for (size_t i = 0; i < 1000; ++i) {
		stepOfComputation(image, sphere, 0, 0, i, l, bp, 1.0);
		std::cerr << i << '\n';
	}
	writePolyData(sphere, "sphere1.vtk");
	vtkSmartPointer<vtkResliceImageViewer> viewer =
		vtkSmartPointer<vtkResliceImageViewer>::New();
	viewer->SetInputData(image);
	viewer->SetupInteractor(vtkSmartPointer<vtkRenderWindowInteractor>::New());
	viewer->Render();
	viewer->GetRenderer()->ResetCamera();
	viewer->Render();
	renderPolyData(sphere, viewer);


	viewer->GetInteractor()->Start();
	//cin.get();
	return 0;
}
