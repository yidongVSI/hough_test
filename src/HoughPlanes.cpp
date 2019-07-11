/**
 * This program accepts a vtp file containing a point cloud and attempts to fit many planes
 * to the points. The output is a vtp file containing the input points colored by which plane,
 * if any, they belong to. The colors of the planes are randomly assigned.
 */
#include <fstream>

#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include "vtkHoughPlanes.h"

int main(int argc, char *argv[])
{
  // Verify command line arguments
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " input.vtp output.vtp" << endl;
    exit(0);
  }

  // Parse command line arguments
  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];

  // Read the input file
  vtkSmartPointer<vtkXMLPolyDataReader> reader =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(inputFileName.c_str());
  reader->Update();

  vtkSmartPointer<vtkHoughPlanes> houghPlanes =
      vtkSmartPointer<vtkHoughPlanes>::New();
  houghPlanes->SetInputConnection(reader->GetOutputPort());
  houghPlanes->SetMaxDist(2.00);
  houghPlanes->SetMinDist(0.10);
  houghPlanes->SetAccumulatorMax(50);
  houghPlanes->SetMinSizeAllPoints(5);
  houghPlanes->SetRhoNum(200);
  houghPlanes->SetThetaNum(360);
  houghPlanes->SetPhiNum(176);
  houghPlanes->SetRhoMax(5.00);
  houghPlanes->SetMaxPointPlaneDist(0.025);
  houghPlanes->SetMaxPlanes(10);
  houghPlanes->SetMinPlaneSize(100);
  houghPlanes->SetMinPlanarity(0.03);
  houghPlanes->SetPlaneRatio(0.5);
  houghPlanes->SetPointDist(0.025);
  houghPlanes->SetPeakWindow(false);
  houghPlanes->SetWindowSize(8);
  houghPlanes->SetTrashMax(20);
  houghPlanes->SetAccumulatorType(1);
  houghPlanes->SetHoughAlgorithm(vtkHoughPlanes::Randomized);
  houghPlanes->Update();

  // Write output points colored by plane
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputConnection(houghPlanes->GetOutputPort());
  writer->SetFileName(outputFileName.c_str());
  writer->Write();

  // write out xyzRGB
  std::string outputFileNameXYZ = outputFileName + ".xyz";
  ofstream ofs;
  ofs.open(outputFileNameXYZ.c_str());
  ofs << "x y z r g b\n";
  vector<Point> out_pts = houghPlanes->get_Hough().coloredPoints;
  for (vector<Point>::iterator it = out_pts.begin(); it != out_pts.end(); ++it)
  {
    ofs << std::setprecision(10)
        << it->x << ' ' << it->y << ' ' << it->z << ' '
        << (int)it->rgb[0] << ' ' << (int)it->rgb[1] << ' ' << (int)it->rgb[2]
        << std::endl;
    std::cout << std::setprecision(8)
              << "Point: (" << it->x << ' ' << it->y << ' ' << it->z << ')'
              << " Color: (" << (int)it->rgb[0] << ' ' << (int)it->rgb[1] << ' ' << (int)it->rgb[2] << ')'
              << std::endl;
  }
  ofs.close();
  // Write output as xyzRGB

  return EXIT_SUCCESS;
}
