// me
#include "vtkPolyDataNormalsCentroids.h"
// vtk
#include <vtkObjectFactory.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkMath.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkTriangleStrip.h>
#include <vtkPriorityQueue.h>
#include <vtkNew.h>
vtkStandardNewMacro(vtkPolyDataNormalsCentroid);
void vtkPolyDataNormalsCentroid::PrintSelf(ostream & os, vtkIndent indent)
{
	Superclass::PrintSelf(os, indent);
}

int vtkPolyDataNormalsCentroid::RequestData(vtkInformation *vtkNotUsed(info), vtkInformationVector ** inputVector, vtkInformationVector * outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType npts = 0;
  vtkIdType *pts = 0;
  vtkIdType numNewPts;
  double flipDirection=1.0;
  vtkIdType numPolys, numStrips;
  vtkIdType cellId;
  vtkIdType numPts;
  vtkPoints *inPts;
  vtkCellArray *inPolys, *inStrips, *polys;
  vtkPoints *newPts = NULL;
  vtkFloatArray *newNormals;
  vtkPointData *pd, *outPD;
  vtkDataSetAttributes* outCD = output->GetCellData();
  double n[3];
  vtkCellArray *newPolys;
  vtkIdType ptId, oldId;

  vtkDebugMacro(<<"Generating surface normals");

  numPolys=input->GetNumberOfPolys();
  numStrips=input->GetNumberOfStrips();
  if ( (numPts=input->GetNumberOfPoints()) < 1 )
  {
    vtkDebugMacro(<<"No data to generate normals for!");
    return 1;
  }


  // If there is nothing to do, pass the data through
  if ( (this->ComputePointNormals == 0 && this->ComputeCellNormals == 0) ||
       (numPolys < 1 && numStrips < 1) )
  { //don't do anything! pass data through
    output->CopyStructure(input);
    output->GetPointData()->PassData(input->GetPointData());
    output->GetCellData()->PassData(input->GetCellData());
    return 1;
  }

  if (numStrips < 1)
  {
    output->GetCellData()->PassData(input->GetCellData());
  }

  // Load data into cell structure.  We need two copies: one is a
  // non-writable mesh used to perform topological queries.  The other
  // is used to write into and modify the connectivity of the mesh.
  //
  inPts = input->GetPoints();
  inPolys = input->GetPolys();
  inStrips = input->GetStrips();

  this->OldMesh = vtkPolyData::New();
  this->OldMesh->SetPoints(inPts);
  if ( numStrips > 0 ) //have to decompose strips into triangles
  {
    vtkDataSetAttributes* inCD = input->GetCellData();
    // When we have triangle strips, make sure to create and copy
    // the cell data appropriately. Since strips are broken into
    // triangles, cell data cannot be passed as it is and needs to
    // be copied tuple by tuple.
    outCD->CopyAllocate(inCD);
    if ( numPolys > 0 )
    {
      polys = vtkCellArray::New();
      polys->DeepCopy(inPolys);
      vtkNew<vtkIdList> ids;
      ids->SetNumberOfIds(numPolys);
      for (vtkIdType i=0; i<numPolys; i++)
      {
        ids->SetId(i, i);
      }
      outCD->CopyData(inCD, ids.GetPointer(), ids.GetPointer());
    }
    else
    {
      polys = vtkCellArray::New();
      polys->Allocate(polys->EstimateSize(numStrips,5));
    }
    vtkIdType inCellIdx = numPolys;
    vtkIdType outCellIdx = numPolys;
    for ( inStrips->InitTraversal(); inStrips->GetNextCell(npts,pts); inCellIdx++)
    {
      vtkTriangleStrip::DecomposeStrip(npts, pts, polys);
      // Copy the cell data for the strip to each triangle.
      for (vtkIdType i=0; i<npts-2; i++)
      {
        outCD->CopyData(inCD, inCellIdx, outCellIdx++);
      }
    }
    this->OldMesh->SetPolys(polys);
    polys->Delete();
    numPolys = polys->GetNumberOfCells();//added some new triangles
  }
  else
  {
    this->OldMesh->SetPolys(inPolys);
    polys = inPolys;
  }
  this->OldMesh->BuildLinks();
  this->UpdateProgress(0.10);

  pd = input->GetPointData();
  outPD = output->GetPointData();

  this->NewMesh = vtkPolyData::New();
  this->NewMesh->SetPoints(inPts);
  // create a copy because we're modifying it
  newPolys = vtkCellArray::New();
  newPolys->DeepCopy(polys);
  this->NewMesh->SetPolys(newPolys);
  this->NewMesh->BuildCells(); //builds connectivity

  // The visited array keeps track of which polygons have been visited.
  //
  if ( this->Consistency || this->Splitting || this->AutoOrientNormals )
  {
    this->Visited = new int[numPolys];
    memset(this->Visited, VTK_CELL_NOT_VISITED, numPolys*sizeof(int));
    this->CellIds = vtkIdList::New();
    this->CellIds->Allocate(VTK_CELL_SIZE);
  }
  else
  {
    this->Visited = NULL;
  }

  //  Traverse all polygons insuring proper direction of ordering.  This
  //  works by propagating a wave from a seed polygon to the polygon's
  //  edge neighbors. Each neighbor may be reordered to maintain consistency
  //  with its (already checked) neighbors.
  //
  this->NumFlips = 0;
  if (this->AutoOrientNormals)
  {
    // No need to check this->Consistency. It's implied.

    // Ok, here's the basic idea: the "left-most" polygon should
    // have its outward pointing normal facing left. If it doesn't,
    // reverse the vertex order. Then use it as the seed for other
    // connected polys. To find left-most polygon, first find left-most
    // point, and examine neighboring polys and see which one
    // has a normal that's "most aligned" with the X-axis. This process
    // will need to be repeated to handle all connected components in
    // the mesh. Report bugs/issues to cvolpe@ara.com.
    int foundLeftmostCell;
    vtkIdType leftmostCellID=-1, currentPointID, currentCellID;
    vtkIdType *leftmostCells;
    unsigned short nleftmostCells;
    vtkIdType *cellPts;
    vtkIdType nCellPts;
    int cIdx;
    double bestNormalAbsXComponent;
    int bestReverseFlag;
    vtkPriorityQueue *leftmostPoints = vtkPriorityQueue::New();
    this->Wave = vtkIdList::New();
    this->Wave->Allocate(numPolys/4+1,numPolys);
    this->Wave2 = vtkIdList::New();
    this->Wave2->Allocate(numPolys/4+1,numPolys);

    // Put all the points in the priority queue, based on x coord
    // So that we can find leftmost point
    leftmostPoints->Allocate(numPts);
    for (ptId=0; ptId < numPts; ptId++)
    {
      leftmostPoints->Insert(inPts->GetPoint(ptId)[0],ptId);
    }

    // Repeat this while loop as long as the queue is not empty,
    // because there may be multiple connected components, each of
    // which needs to be seeded independently with a correctly
    // oriented polygon.
    while (leftmostPoints->GetNumberOfItems())
    {
      foundLeftmostCell = 0;
      // Keep iterating through leftmost points and cells located at
      // those points until I've got a leftmost point with
      // unvisited cells attached and I've found the best cell
      // at that point
      do {
        currentPointID = leftmostPoints->Pop();
        this->OldMesh->GetPointCells(currentPointID, nleftmostCells, leftmostCells);
        bestNormalAbsXComponent = 0.0;
        bestReverseFlag = 0;
        for (cIdx = 0; cIdx < nleftmostCells; cIdx++)
        {
          currentCellID = leftmostCells[cIdx];
          if (this->Visited[currentCellID] == VTK_CELL_VISITED)
          {
            continue;
          }
          this->OldMesh->GetCellPoints(currentCellID, nCellPts, cellPts);
          vtkPolygon::ComputeNormal(inPts, nCellPts, cellPts, n);
          // Ok, see if this leftmost cell candidate is the best
          // so far
          if (fabs(n[0]) > bestNormalAbsXComponent)
          {
            bestNormalAbsXComponent = fabs(n[0]);
            leftmostCellID = currentCellID;
            // If the current leftmost cell's normal is pointing to the
            // right, then the vertex ordering is wrong
            bestReverseFlag = (n[0] > 0);
            foundLeftmostCell = 1;
          } // if this normal is most x-aligned so far
        } // for each cell at current leftmost point
      } while (leftmostPoints->GetNumberOfItems() && !foundLeftmostCell);
      if (foundLeftmostCell)
      {
        // We've got the seed for a connected component! But do
        // we need to flip it first? We do, if it was pointed the wrong
        // way to begin with, or if the user requested flipping all
        // normals, but if both are true, then we leave it as it is.
        if (bestReverseFlag ^ this->FlipNormals)
        {
          this->NewMesh->ReverseCell(leftmostCellID);
          this->NumFlips++;
        }
        this->Wave->InsertNextId(leftmostCellID);
        this->Visited[leftmostCellID] = VTK_CELL_VISITED;
        this->TraverseAndOrder();
        this->Wave->Reset();
        this->Wave2->Reset();
      } // if found leftmost cell
    } // Still some points in the queue
    this->Wave->Delete();
    this->Wave2->Delete();
    leftmostPoints->Delete();
    vtkDebugMacro(<<"Reversed ordering of " << this->NumFlips << " polygons");
  } // automatically orient normals
  else
  {
    if ( this->Consistency )
    {
      this->Wave = vtkIdList::New();
      this->Wave->Allocate(numPolys/4+1,numPolys);
      this->Wave2 = vtkIdList::New();
      this->Wave2->Allocate(numPolys/4+1,numPolys);
      for (cellId=0; cellId < numPolys; cellId++)
      {
        if ( this->Visited[cellId] == VTK_CELL_NOT_VISITED)
        {
          if ( this->FlipNormals )
          {
            this->NumFlips++;
            this->NewMesh->ReverseCell(cellId);
          }
          this->Wave->InsertNextId(cellId);
          this->Visited[cellId] = VTK_CELL_VISITED;
          this->TraverseAndOrder();
        }

        this->Wave->Reset();
        this->Wave2->Reset();
      }

      this->Wave->Delete();
      this->Wave2->Delete();
      vtkDebugMacro(<<"Reversed ordering of " << this->NumFlips << " polygons");
    }//Consistent ordering
  } // don't automatically orient normals

  this->UpdateProgress(0.333);

  //  Initial pass to compute polygon normals without effects of neighbors
  //
  this->PolyNormals = vtkFloatArray::New();
  this->PolyNormals->SetNumberOfComponents(3);
  this->PolyNormals->Allocate(3*numPolys);
  this->PolyNormals->SetName("Normals");
  this->PolyNormals->SetNumberOfTuples(numPolys);

  for (cellId=0, newPolys->InitTraversal(); newPolys->GetNextCell(npts,pts);
       cellId++ )
  {
    if ((cellId % 1000) == 0)
    {
      this->UpdateProgress (0.333 + 0.333 * (double) cellId / (double) numPolys);
      if (this->GetAbortExecute())
      {
        break;
      }
    }
    vtkPolygon::ComputeNormal(inPts, npts, pts, n);
    this->PolyNormals->SetTuple(cellId,n);
  }

  // Split mesh if sharp features
  if ( this->Splitting )
  {
    //  Traverse all nodes; evaluate loops and feature edges.  If feature
    //  edges found, split mesh creating new nodes.  Update polygon
    // connectivity.
    //
      this->CosAngle = cos( vtkMath::RadiansFromDegrees( this->FeatureAngle) );
    //  Splitting will create new points.  We have to create index array
    // to map new points into old points.
    //
    this->Map = vtkIdList::New();
    this->Map->SetNumberOfIds(numPts);
    for (vtkIdType i=0; i < numPts; i++)
    {
      this->Map->SetId(i,i);
    }

    for (ptId=0; ptId < numPts; ptId++)
    {
      this->MarkAndSplit(ptId);
    }//for all input points

    numNewPts = this->Map->GetNumberOfIds();

    vtkDebugMacro(<<"Created " << numNewPts-numPts << " new points");

    //  Now need to map attributes of old points into new points.
    //
    outPD->CopyNormalsOff();
    outPD->CopyAllocate(pd,numNewPts);

    newPts = vtkPoints::New();

    // set precision for the points in the output
    if(this->OutputPointsPrecision == vtkAlgorithm::DEFAULT_PRECISION)
    {
      vtkPointSet *inputPointSet = vtkPointSet::SafeDownCast(input);
      if(inputPointSet)
      {
        newPts->SetDataType(inputPointSet->GetPoints()->GetDataType());
      }
      else
      {
        newPts->SetDataType(VTK_FLOAT);
      }
    }
    else if(this->OutputPointsPrecision == vtkAlgorithm::SINGLE_PRECISION)
    {
      newPts->SetDataType(VTK_FLOAT);
    }
    else if(this->OutputPointsPrecision == vtkAlgorithm::DOUBLE_PRECISION)
    {
      newPts->SetDataType(VTK_DOUBLE);
    }

    newPts->SetNumberOfPoints(numNewPts);
    for (ptId=0; ptId < numNewPts; ptId++)
    {
      oldId = this->Map->GetId(ptId);
      newPts->SetPoint(ptId,inPts->GetPoint(oldId));
      outPD->CopyData(pd,oldId,ptId);
    }
    this->Map->Delete();
  } //splitting

  else //no splitting, so no new points
  {
    numNewPts = numPts;
    outPD->CopyNormalsOff();
    outPD->PassData(pd);
  }

  if ( this->Consistency || this->Splitting )
  {
    delete [] this->Visited;
    this->CellIds->Delete();
  }

  this->UpdateProgress(0.80);

  //  Finally, traverse all elements, computing polygon normals and
  //  accumulating them at the vertices.
  //
  if ( this->FlipNormals && ! this->Consistency )
  {
    flipDirection = -1.0;
  }

  newNormals = vtkFloatArray::New();
  newNormals->SetNumberOfComponents(3);
  newNormals->SetNumberOfTuples(numNewPts);
  newNormals->SetName("Normals");
  float *fNormals = newNormals->WritePointer(0, 3 * numNewPts);
  std::fill_n(fNormals, 3 * numNewPts, 0);

  float *fPolyNormals = this->PolyNormals->WritePointer(0, 3 * numPolys);

  if (this->ComputePointNormals)
  {
    for (cellId=0, newPolys->InitTraversal(); newPolys->GetNextCell(npts, pts);
         ++cellId)
    {
      for (vtkIdType i = 0; i < npts; ++i)
      {
        fNormals[3 * pts[i]] += fPolyNormals[3 * cellId];
        fNormals[3 * pts[i] + 1] += fPolyNormals[3 * cellId + 1];
        fNormals[3 * pts[i] + 2] += fPolyNormals[3 * cellId + 2];
      }
    }

    for (vtkIdType i = 0; i < numNewPts; ++i)
    {
      const double length = sqrt(fNormals[3 * i] * fNormals[3 * i] +
                                 fNormals[3 * i + 1] * fNormals[3 * i + 1] +
                                 fNormals[3 * i + 2] * fNormals[3 * i + 2]
                                 ) * flipDirection;
      if (length != 0.0)
      {
        fNormals[3 * i] /= length;
        fNormals[3 * i + 1] /= length;
        fNormals[3 * i + 2] /= length;
      }
    }
  }

  //  Update ourselves.  If no new nodes have been created (i.e., no
  //  splitting), we can simply pass data through.
  //
  if ( ! this->Splitting )
  {
    output->SetPoints(inPts);
  }

  //  If there is splitting, then have to send down the new data.
  //
  else
  {
    output->SetPoints(newPts);
    newPts->Delete();
  }

  if (this->ComputeCellNormals)
  {
    outCD->SetNormals(this->PolyNormals);
  }
  this->PolyNormals->Delete();

  if (this->ComputePointNormals)
  {
    outPD->SetNormals(newNormals);
  }
  newNormals->Delete();

  output->SetPolys(newPolys);
  newPolys->Delete();

  // copy the original vertices and lines to the output
  output->SetVerts(input->GetVerts());
  output->SetLines(input->GetLines());

  this->OldMesh->Delete();
  this->NewMesh->Delete();

  return 1;
}
