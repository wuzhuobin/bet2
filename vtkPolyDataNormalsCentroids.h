#ifndef __VTK_POLY_DATA_NORMALS_CENTROID_H__
#define __VTK_POLY_DATA_NORMALS_CENTROID_H__
#pragma once
#include "vtkPolyDataAlgorithm.h"

class vtkFloatArray;
class vtkIdList;
class vtkPolyData;

class vtkPolyDataNormalsCentroids : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkPolyDataNormalsCentroids, vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

	/**
	* Construct with feature angle=30, splitting and consistency turned on,
	* flipNormals turned off, and non-manifold traversal turned on.
	* ComputePointNormals is on and ComputeCellNormals is off.
	*/
	static vtkPolyDataNormalsCentroids *New();

	//@{
	/**
	* Specify the angle that defines a sharp edge. If the difference in
	* angle across neighboring polygons is greater than this value, the
	* shared edge is considered "sharp".
	*/
	vtkSetClampMacro(FeatureAngle, double, 0.0, 180.0);
	vtkGetMacro(FeatureAngle, double);
	//@}

	//@{
	/**
	* Turn on/off the splitting of sharp edges.
	*/
	vtkSetMacro(Splitting, int);
	vtkGetMacro(Splitting, int);
	vtkBooleanMacro(Splitting, int);
	//@}

	//@{
	/**
	* Turn on/off the enforcement of consistent polygon ordering.
	*/
	vtkSetMacro(Consistency, int);
	vtkGetMacro(Consistency, int);
	vtkBooleanMacro(Consistency, int);
	//@}

	//@{
	/**
	* Turn on/off the automatic determination of correct normal
	* orientation. NOTE: This assumes a completely closed surface
	* (i.e. no boundary edges) and no non-manifold edges. If these
	* constraints do not hold, all bets are off. This option adds some
	* computational complexity, and is useful if you don't want to have
	* to inspect the rendered image to determine whether to turn on the
	* FlipNormals flag. However, this flag can work with the FlipNormals
	* flag, and if both are set, all the normals in the output will
	* point "inward".
	*/
	vtkSetMacro(AutoOrientNormals, int);
	vtkGetMacro(AutoOrientNormals, int);
	vtkBooleanMacro(AutoOrientNormals, int);
	//@}

	//@{
	/**
	* Turn on/off the computation of point normals.
	*/
	vtkSetMacro(ComputePointNormals, int);
	vtkGetMacro(ComputePointNormals, int);
	vtkBooleanMacro(ComputePointNormals, int);
	//@}

	//@{
	/**
	* Turn on/off the computation of cell normals.
	*/
	vtkSetMacro(ComputeCellNormals, int);
	vtkGetMacro(ComputeCellNormals, int);
	vtkBooleanMacro(ComputeCellNormals, int);
	//@}

	//@{
	/**
	* Turn on/off the computation of point centroids.
	*/
	vtkSetMacro(ComputePointCentroids, int);
	vtkGetMacro(ComputePointCentroids, int);
	vtkBooleanMacro(ComputePointCentroids, int);
	//@}

	//@{
	/**
	* Turn on/off the computation of cell centroids.
	*/
	vtkSetMacro(ComputeCellCentroids, int);
	vtkGetMacro(ComputeCellCentroids, int);
	vtkBooleanMacro(ComputeCellCentroids, int);
	//@}

	//@{
	/**
	* Turn on/off the global flipping of normal orientation. Flipping
	* reverves the meaning of front and back for Frontface and Backface
	* culling in vtkProperty.  Flipping modifies both the normal
	* direction and the order of a cell's points.
	*/
	vtkSetMacro(FlipNormals, int);
	vtkGetMacro(FlipNormals, int);
	vtkBooleanMacro(FlipNormals, int);
	//@}

	//@{
	/**
	* Turn on/off traversal across non-manifold edges. This will prevent
	* problems where the consistency of polygonal ordering is corrupted due
	* to topological loops.
	*/
	vtkSetMacro(NonManifoldTraversal, int);
	vtkGetMacro(NonManifoldTraversal, int);
	vtkBooleanMacro(NonManifoldTraversal, int);
	//@}

	//@{
	/**
	* Set/get the desired precision for the output types. See the documentation
	* for the vtkAlgorithm::DesiredOutputPrecision enum for an explanation of
	* the available precision settings.
	*/
	vtkSetClampMacro(OutputPointsPrecision, int, SINGLE_PRECISION, DEFAULT_PRECISION);
	vtkGetMacro(OutputPointsPrecision, int);
	//@}

protected:
	vtkPolyDataNormalsCentroids();
	~vtkPolyDataNormalsCentroids() VTK_OVERRIDE {}

	// Usual data generation method
	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) VTK_OVERRIDE;

	double FeatureAngle;
	int Splitting;
	int Consistency;
	int FlipNormals;
	int AutoOrientNormals;
	int NonManifoldTraversal;
	int ComputePointNormals;
	int ComputeCellNormals;
	int ComputePointCentroids;
	int ComputeCellCentroids;
	int NumFlips;
	int OutputPointsPrecision;

private:
	vtkIdList *Wave;
	vtkIdList *Wave2;
	vtkIdList *CellIds;
	vtkIdList *Map;
	vtkPolyData *OldMesh;
	vtkPolyData *NewMesh;
	int *Visited;
	vtkFloatArray *PolyNormals;
	vtkFloatArray *PolyCentroids;
	double CosAngle;

	// Uses the list of cell ids (this->Wave) to propagate a wave of
	// checked and properly ordered polygons.
	void TraverseAndOrder(void);

	// Check the point id give to see whether it lies on a feature
	// edge. If so, split the point (i.e., duplicate it) to topologically
	// separate the mesh.
	void MarkAndSplit(vtkIdType ptId);

private:
	vtkPolyDataNormalsCentroids(const vtkPolyDataNormalsCentroids&) VTK_DELETE_FUNCTION;
	void operator=(const vtkPolyDataNormalsCentroids&) VTK_DELETE_FUNCTION;
};

#endif // !__VTK_POLY_DATA_NORMALS_CENTROID_H__
