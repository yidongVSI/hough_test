#include "vtkAutoInit.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2);

#include <vtkQuadric.h>
#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkOutlineFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkImageData.h>
#include <vtkSphereSource.h>

#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <iostream>
#include <fstream>

void VTPShow(std::string const& inputFileName, std::string const& outFileName);

int main(int argc, char *argv[])
{
  if(argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " input.vtp output.vtp" << std::endl;
    exit(0);
  }

  // parse command line arguments
  std::string inputFileName = argv[1];
  std::string outFileName = argv[2];
  // visualize
  VTPShow(inputFileName, outFileName);
  return 0;
}

void VTPShow(std::string const& inputFileName, std::string const& outFileName)
{
  // Read the input file
  vtkSmartPointer<vtkXMLPolyDataReader> reader =
    vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(inputFileName.c_str());
  reader->Update();

  vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();

  // map the contours to graphical primitives
  vtkPolyDataMapper *contMapper = vtkPolyDataMapper::New();
  contMapper->SetInputData(polydata);

   // create an actor for the contours
  vtkActor *contActor = vtkActor::New();
  contActor->SetMapper(contMapper);

    // a renderer and render window
  vtkRenderer *ren1 = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren1);

    // an interactor
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);

  // add the actors to the scene
  ren1->AddActor(contActor);
  ren1->SetBackground(1,1,1); // Background color white

    // render an image (lights and cameras are created automatically)
  renWin->Render();

    // begin mouse interaction
  iren->Start();

}