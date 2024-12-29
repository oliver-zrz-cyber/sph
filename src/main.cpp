#include "func.h"
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include <igl/readOBJ.h>

int main()
{
	float deltat = 0.002;
	WaterSource source;
	source.pour_deltaT = 300;
	source.initial_velocity = 2;
	source.sample_num = 60;
	source.sourcemesh.load_mesh("E:/research/SPH/MyImplementation/watersource3.obj");
	Eigen::MatrixXd meshV;
	Eigen::MatrixXi meshF;
	igl::readOBJ("E:/research/SPH/MyImplementation/cup2.obj", meshV, meshF);
	ParticleSystem PS;
	std::vector<position>points;
	PS.InitSystem(0.5, 0.5, 100, 300);// -,100,20
	int voxelize_n = 20;
	polyscope::init();
	polyscope::options::groundPlaneHeightMode = polyscope::GroundPlaneHeightMode::Manual;
	polyscope::options::groundPlaneHeight = 0.0f;
	polyscope::options::buildGui = false;
	polyscope::view::setUpDir(polyscope::UpDir::ZUp);
	polyscope::view::farClipRatio = 100.0;
	polyscope::view::lookAt(glm::vec3{ 3,0, 1 }, glm::vec3{ 0., 0., 0.8 });
	auto my_cup = polyscope::registerSurfaceMesh("input mesh", meshV, meshF);
	my_cup->setTransparency(0.35);
	my_cup->setSurfaceColor({ 0,0,0 });
	auto my_cloud = polyscope::registerPointCloud("cloud", PS.get_particles());
	my_cloud->setPointRadius(0.005);
	for (int k = 0; k < 2000; ++k)
	{
		step(k, deltat, source, PS, 0.05,voxelize_n);
		if (k < source.pour_deltaT)
			my_cloud = polyscope::registerPointCloud("cloud", PS.get_particles());
		else
			my_cloud->updatePointPositions(PS.get_particles());
		//my_cloud->setTransparency(0.114);
		my_cloud->setPointColor({ 28/255.f,106/255.f,227/255.f });
		polyscope::frameTick();
	}

}

