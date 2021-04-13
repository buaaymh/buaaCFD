// Copyright 2021 Minghao Yang

#ifndef INCLUDE_BUAA_MESH_VTK_READER_HPP_
#define INCLUDE_BUAA_MESH_VTK_READER_HPP_

// C++ system headers:
#include <cassert>
#include <string>
#include <memory>
#include <stdexcept>
#include <utility>
// For `.vtk` files:
#include <vtkDataSet.h>
#include <vtkDataSetReader.h>
// For `.vtu` files:
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
// DataSetAttributes:
#include <vtkFieldData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
// Cells:
#include <vtkCellType.h>  // define types of cells
#include <vtkCellTypes.h>
#include <vtkCell.h>
#include <vtkTriangle.h>
// Helpers:
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>

namespace buaa {
namespace mesh {
namespace vtk {

template <class Mesh>
class Reader {
 private:
  std::unique_ptr<Mesh> mesh_;
  using IdType = typename Mesh::NodeId;

 public:
  bool ReadFromFile(const std::string& file_name) {
    auto vtk_data_set = FileNameToDataSet(file_name.c_str());
    if (vtk_data_set) {
      auto vtk_data_set_owner = vtkSmartPointer<vtkDataSet>();
      vtk_data_set_owner.TakeReference(vtk_data_set);
      mesh_.reset(new Mesh());
      ReadNodes(vtk_data_set);
      ReadCells(vtk_data_set);
      ReadNodeData(vtk_data_set);
      ReadCellData(vtk_data_set);
    } else {
      throw std::runtime_error("Unable to read \"" + file_name + "\".");
    }
    return true;
  }
  std::unique_ptr<Mesh> GetMesh() {
    auto temp = std::make_unique<Mesh>();
    std::swap(temp, mesh_);
    return temp;
  }

 private:
  void ReadNodes(vtkDataSet* vtk_data_set) {
    int n = vtk_data_set->GetNumberOfPoints();
    mesh_->SetNodesNum(n);
    for (int i = 0; i < n; i++) {
      auto xyz = vtk_data_set->GetPoint(i);
      mesh_->EmplaceNode(i, xyz[0], xyz[1]);
    }
  }
  void ReadNodeData(vtkDataSet* vtk_data_set) {
  }
  void ReadCells(vtkDataSet* vtk_data_set) {
    int n = vtk_data_set->GetNumberOfCells();
    mesh_->SetCellsNum(n);
    for (int i = 0; i < n; i++) {
      auto cell = vtk_data_set->GetCell(i);
      auto type = vtk_data_set->GetCellType(i);
      auto id_list = cell->GetPointIds();
      if (type == 5) {
        IdType a = id_list->GetId(0);
        IdType b = id_list->GetId(1);
        IdType c = id_list->GetId(2);
        mesh_->EmplaceCell(i, {a, b, c});
      } else {
        continue;
      }
    }  // for each cell
  }
  void ReadCellData(vtkDataSet* vtk_data_set) {
  }
  vtkDataSet* FileNameToDataSet(const char* file_name) {
    vtkDataSet* vtk_data_set{nullptr};
    auto extension = std::string(file_name + strlen(file_name) - 4);
    // Dispatch based on the file extension
    if (extension == ".vtk") {
      BindReader<vtkDataSetReader>(file_name, &vtk_data_set);
    } else if (extension == ".vtu") {
      BindReader<vtkXMLUnstructuredGridReader>(file_name, &vtk_data_set);
    } else {
      throw std::invalid_argument("Only `.vtk` and `.vtu` are supported!");
    }
    return vtk_data_set;
  }
  template <class Reader>
  void BindReader(const char* file_name, vtkDataSet** vtk_data_set) {
    auto reader = vtkSmartPointer<Reader>::New();
    reader->SetFileName(file_name);
    reader->Update();
    *vtk_data_set = vtkDataSet::SafeDownCast(reader->GetOutput());
    if (*vtk_data_set) {
      (*vtk_data_set)->Register(reader);
    }
  }
};

template <class Mesh>
std::unique_ptr<Mesh> Read(const std::string& file_name) {
  auto reader = Reader<Mesh>();
  reader.ReadFromFile(file_name);
  return reader.GetMesh();
}

}  // namespace vtk
}  // namespace mesh
}  // namespace buaa

#endif  // INCLUDE_BUAA_MESH_VTK_READER_HPP_
