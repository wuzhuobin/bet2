//#include "bet2.h"
#include <vtkNIFTIImageReader.h>
#include <vtkNIFTIImageWriter.h>
#include "vtkBet2Filter.h"
#include <vtkSmartPointer.h>
//using namespace bet2;
int main(int argc, char **argv) {
	argv[1] = "C:/Users/jieji/Desktop/PlaqueQuant/JackyData/nifti_corrected/CUBE T1 corrected.nii";
	argv[2] = "bbb.nii";
	auto reader = 
		vtkSmartPointer<vtkNIFTIImageReader>::New();
	reader->SetFileName(argv[1]);
	reader->Update();
	auto f = vtkSmartPointer<vtkBet2Filter>::New();
	f->SetInputConnection(reader->GetOutputPort());
	f->Update();
	auto writer = vtkSmartPointer<vtkNIFTIImageWriter>::New();
	writer->SetInputConnection(f->GetOutputPort());
	writer->SetFileName(argv[2]);
	writer->Write();
	std::cerr << "Finished. \n";
	cin.get();
	return 0;
}
