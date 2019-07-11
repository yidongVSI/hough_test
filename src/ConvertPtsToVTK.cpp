/**
 * This program accepts a xyz, xyzI, xyzRGB, xyzRGBI file containing a point cloud without any header and attempts to
   to fit many planes to the points.  The output is a vtp file containing the input points colored by which plane,
   if any, they belong to.  The colors of the planes are randomly assigned.
   input of this executable:
    1. input point cloud csv files
    2. output filename
    3. configuration file for Hough parameters
 */
#include <map>
#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>

#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkPTSReader.h>
#include <vtkExecutionTimer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPointData.h>

#include "vtkHoughPlanes.h"
#include "hough.h"
#include "slam6d/point.h"


int main(int argc, char *argv[])
{
  // Verify command line arguments
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " input_point_cloud.csv output.vtp config." << std::endl;
    return -1;
  }

  // Parse command line arguments
  std::string inputFilename  = argv[1];
  std::string outputFileName = argv[2];
  std::string configFileName = argv[3];

  // load the configuration
  YAML::Node config = YAML::LoadFile(configFileName.c_str());
  std::string AlgorithmType = config["AlgorithmType"].as<std::string>();
  double MaxDist = config["MaxDist"].as<double>();
  double MinDist = config["MinDist"].as<double>();
  unsigned int AccumulatorMax   = config["AccumulatorMax"].as<unsigned int>();
  unsigned int MinSizeAllPoints = config["MinSizeAllPoints"].as<unsigned int>();
  unsigned int RhoNum   = config["RhoNum"].as<unsigned int>();
  unsigned int ThetaNum = config["ThetaNum"].as<unsigned int>();
  unsigned int PhiNum   = config["PhiNum"].as<unsigned int>();
  unsigned int RhoMax   = config["RhoMax"].as<unsigned int>();
  double MaxPointPlaneDist    = config["MaxPointPlaneDist"].as<double>();
  unsigned int MaxNumPlanes   = config["MaxNumPlanes"].as<unsigned int>();
  unsigned int MinPlaneSize   = config["MinPlaneSize"].as<unsigned int>();
  double MinPlanarity = config["MinPlanarity"].as<double>();
  double PlaneRatio   = config["PlaneRatio"].as<double>();
  double PointDist    = config["PointDist"].as<double>();
  bool   PeakWindow   = config["PeakWindow"].as<bool>();
  unsigned int WindowSize = config["WindowSize"].as<unsigned int>();
  unsigned int TrashMax   = config["TrashMax"].as<unsigned int>();
  unsigned int AccumulatorType = config["AccumulatorType"].as<unsigned int>();
  bool   IsDebug = config["IsDebug"].as<bool>();

  // read the point cloud file format
  vtkSmartPointer<vtkPTSReader> reader = vtkSmartPointer<vtkPTSReader>::New();
  reader->SetFileName(inputFilename.c_str());
  reader->SetLimitToMaxNumberOfPoints(false);
  reader->Update();

  // Debug screen output
  if (IsDebug)
  {
    std::cout << "-------- Hough Parameters ---------" << std::endl;

    std::cout << "AlgorithmType: " << AlgorithmType << '\n';
    std::cout << "MaxDist: " << MaxDist << ", MinDist: " << MinDist << '\n';
    std::cout << "AccumulatorMax: " << AccumulatorMax << '\n';
    std::cout << "MinSizeAllPoints: " << MinSizeAllPoints << '\n';
    std::cout << "Resolution of Hough cell: [" << RhoNum << ','
                                               << ThetaNum << ','
                                               << PhiNum << "]\n";
    std::cout << "RhoMax: " << RhoMax << '\n';
    std::cout << "MaxPointPlaneDist: " << MaxPointPlaneDist << '\n';
    std::cout << "MaxNumPlanes: " << MaxNumPlanes << '\n';
    std::cout << "MinPlaneSize: " << MinPlaneSize << '\n';
    std::cout << "MinPlanarity: " << MinPlaneSize << '\n';
    std::cout << "PlaneRatio: " << PlaneRatio << '\n';
    std::cout << "PointDist: "  << PointDist << '\n';
    std::cout << "PeakWindow: " << PeakWindow << '\n';
    std::cout << "WindowSize: " << WindowSize << '\n';
    std::cout << "TrashMax: "   << TrashMax << '\n';
    std::cout << "AccumulatorType: " << AccumulatorType << '\n';
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "-------- Input PointCloud ---------" << std::endl;
    reader->Print(std::cout);
    std::cout << "Number of input points: " << reader->GetOutput()->GetNumberOfPoints()
              << std::endl;
  }



  // Config the Hough Transformer
  vtkSmartPointer<vtkHoughPlanes> hP = vtkSmartPointer<vtkHoughPlanes>::New();
  hP->SetInputConnection(reader->GetOutputPort());
  hP->SetMaxDist(MaxDist);
  hP->SetMinDist(MinDist);
  hP->SetAccumulatorMax(AccumulatorMax);
  hP->SetMinSizeAllPoints(MinSizeAllPoints);
  hP->SetRhoNum(RhoNum);
  hP->SetPhiNum(PhiNum);
  hP->SetThetaNum(ThetaNum);
  hP->SetRhoMax(RhoMax);
  hP->SetMaxPointPlaneDist(MaxPointPlaneDist);
  hP->SetMaxPlanes(MaxNumPlanes);
  hP->SetMinPlaneSize(MinPlaneSize);
  hP->SetMinPlanarity(MinPlanarity);
  hP->SetPlaneRatio(PlaneRatio);
  hP->SetPointDist(PointDist);
  hP->SetPeakWindow(PeakWindow);
  hP->SetWindowSize(WindowSize);
  hP->SetTrashMax(TrashMax);
  hP->SetAccumulatorType(AccumulatorType);

  if (AlgorithmType == "Randomized")
    hP->SetHoughAlgorithm(vtkHoughPlanes::Randomized);
  else if (AlgorithmType == "Standard")
    hP->SetHoughAlgorithm(vtkHoughPlanes::Standard);
  else if (AlgorithmType == "Probabilistic")
    hP->SetHoughAlgorithm(vtkHoughPlanes::Probabilistic);
  else if (AlgorithmType == "Progressive")
    hP->SetHoughAlgorithm(vtkHoughPlanes::Progressive);
  else if (AlgorithmType == "Adaptive")
    hP->SetHoughAlgorithm(vtkHoughPlanes::Adaptive);
  else {
    std::cerr << "Unknown AlgorithmType: " << AlgorithmType << std::endl;
    return -1;
  }
  // hP->Print(std::cout);
  // Set a CPU timer
  vtkSmartPointer<vtkExecutionTimer> timer = vtkSmartPointer<vtkExecutionTimer>::New();
  timer->SetFilter(hP);

  // run
  hP->Update();

  // Done
  std::cout << "TestExecutionTimer: Filter under inspection ("
            << timer->GetFilter()->GetClassName() << ") execution time: "
            << std::fixed << std::setprecision(12)
            << timer->GetElapsedCPUTime() << " sec (CPU), "
            << timer->GetElapsedWallClockTime() << " sec (wall clock)\n";

  // Write output points colored by plane
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputConnection(hP->GetOutputPort());
  writer->SetFileName(outputFileName.c_str());
  writer->Write();

  // Write output points as xyzRGB
  std::string outputFileNameXYZ = outputFileName + ".xyz";
  ofstream ofs;
  ofs.open(outputFileNameXYZ.c_str());
  ofs << "x y z r g b\n";
  std::vector<Point> out_pts = hP->get_Hough().coloredPoints;
  std::cout << "number of output points: " << out_pts.size() << std::endl;
  for (vector<Point>::iterator it = out_pts.begin(); it != out_pts.end(); ++it)
  {
    ofs << std::setprecision(10)
        << it->x << ' ' << it->y << ' ' << it->z << ' '
        << (int)it->rgb[0] << ' ' << (int)it->rgb[1] << ' ' << (int)it->rgb[2]
        << std::endl;
    // std::cout << std::setprecision(8)
    //           << "Point: (" << it->x << ' ' << it->y << ' ' << it->z << ')'
    //           << " Color: (" << (int)it->rgb[0] << ' ' << (int)it->rgb[1] << ' ' << (int)it->rgb[2] << ')'
    //           << std::endl;
  }
  ofs.close();

// // Get the number of points in the polydata
//   vtkIdType idNumPointsInFile = polydata->GetNumberOfPoints();

//   vtkSmartPointer<vtkDoubleArray> array =
//     vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray(arrayName.c_str()));

//   if(array)
//     {
//       std::cout << "Got array " << arrayName
//                 << " with " << idNumPointsInFile << " values"
//                 << std::endl;
//     for(int i = 0; i < idNumPointsInFile; i++)
//       {
//       double value;
//       value = array->GetValue(i);
//       std::cout << i << ": " << value << std::endl;
//       }
//     }
//   else
//     {
//     std::cout << "The file " << filename
//               << " does not have a PointData array named " << arrayName
//               << std::endl;
//     }

//   return EXIT_SUCCESS;


  return EXIT_SUCCESS;

}
