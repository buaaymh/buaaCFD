// Copyright 2021 Minghao Yang

#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "buaa/mesh/data.hpp"
#include "buaa/mesh/dim2.hpp"
#include "buaa/mesh/vtk/reader.hpp"
#include "buaa/mesh/vtk/writer.hpp"
#include "buaa/data/path.hpp"  // defines TEST_DATA_DIR

namespace buaa {
namespace mesh {
namespace vtk {

class ReaderTest : public ::testing::Test {
 protected:
  using MeshType = Mesh<0>;
  using CellType = MeshType::Cell;
  Reader<MeshType> reader;
  const std::string test_data_dir_{TEST_DATA_DIR};
};
TEST_F(ReaderTest, ReadFromFile) {
  EXPECT_TRUE(reader.ReadFromFile(test_data_dir_ + "/tiny.vtk"));
  EXPECT_TRUE(reader.ReadFromFile(test_data_dir_ + "/tiny.vtu"));
}
TEST_F(ReaderTest, GetMesh) {
  for (auto suffix : {".vtk", ".vtu"}) {
    reader.ReadFromFile(test_data_dir_ + "/tiny.vtk");
    auto mesh = reader.GetMesh();
    ASSERT_TRUE(mesh);
    EXPECT_EQ(mesh->CountNodes(), 7);
    EXPECT_EQ(mesh->CountEdges(), 12);
    EXPECT_EQ(mesh->CountCells(), 6);
    // sum of each face's area
    double area = 0.0;
    auto visitor = [&area](const CellType& d) { area += d.Measure(); };
    mesh->ForEachCell(visitor);
    EXPECT_EQ(area, 8.0);
  }
}

class WriterTest : public ::testing::Test {
 public:
  static const char* mesh_name;
 protected:
  const std::string test_data_dir_{TEST_DATA_DIR};
};
const char* WriterTest::mesh_name;
TEST_F(WriterTest, TinyMesh) {
  auto reader = Reader<Mesh<0>>();
  auto writer = Writer<Mesh<0>>();
  for (auto suffix : {".vtk", ".vtu"}) {
    reader.ReadFromFile(test_data_dir_ + "/tiny" + suffix);
    auto mesh_old = reader.GetMesh();
    ASSERT_TRUE(mesh_old);
    // Write the mesh just read:
    writer.SetMesh(mesh_old.get());
    auto filename = test_data_dir_ + std::string("/tiny") + suffix;
    ASSERT_TRUE(writer.WriteToFile(filename));
    // Read the mesh just written:
    reader.ReadFromFile(test_data_dir_ + "/tiny" + suffix);
    auto mesh_new = reader.GetMesh();
    // Check consistency:
    EXPECT_EQ(mesh_old->CountNodes(),
              mesh_new->CountNodes());
    EXPECT_EQ(mesh_old->CountEdges(),
              mesh_new->CountEdges());
    EXPECT_EQ(mesh_old->CountCells(),
              mesh_new->CountCells());
  }
}
TEST_F(WriterTest, MeshWithData) {
  using CellData = mesh::Data<2, 2, 2>;
  using MeshType = Mesh<0, Empty, CellData>;
  // using MeshType = Mesh<0, Empty, Empty>;
  auto reader = Reader<MeshType>();
  auto writer = Writer<MeshType>();
  reader.ReadFromFile(test_data_dir_ + "/medium.vtk");
  auto mesh_old = reader.GetMesh();
  ASSERT_TRUE(mesh_old);
  // Create some data on it:
  MeshType::Cell::scalar_names.at(0) = "X + Y";
  MeshType::Cell::scalar_names.at(1) = "X - Y";
  MeshType::Cell::vector_names.at(0) = "(X, Y)";
  MeshType::Cell::vector_names.at(1) = "(-X, -Y)";
  mesh_old->ForEachCell([](MeshType::Cell& cell) {
    auto center = cell.Center();
    auto y = center.Y();
    auto x = center.X();
    cell.data.scalars.at(0) = x + y;
    cell.data.scalars.at(1) = x - y;
    cell.data.vectors.at(0) = {x, y};
    cell.data.vectors.at(1) = {-x, -y};
  });
  // Write the mesh just read:
  writer.SetMesh(mesh_old.get());
  auto filename = test_data_dir_ + std::string("/medium.vtu");
  ASSERT_TRUE(writer.WriteToFile(filename));
  // Read the mesh just written:
  reader.ReadFromFile(test_data_dir_ + "/medium.vtu");
  auto mesh_new = reader.GetMesh();
  // Check consistency:
  EXPECT_EQ(mesh_old->CountNodes(),
            mesh_new->CountNodes());
  EXPECT_EQ(mesh_old->CountEdges(),
            mesh_new->CountEdges());
  EXPECT_EQ(mesh_old->CountCells(),
            mesh_new->CountCells());
}

}  // namespace vtk
}  // namespace mesh
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  if (argc == 1) {
    buaa::mesh::vtk::WriterTest::mesh_name = "tiny";
  } else {
    buaa::mesh::vtk::WriterTest::mesh_name = argv[1];
  }
  return RUN_ALL_TESTS();
}