#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <random>

struct VoxelIndices {
	int i = 0;
	int j = 0;
	int k = 0;
	int voxelToIndex( int n) {
		return i * n * n + j * n + k;
	}
};



struct position
{
	float x = 0.0, y = 0.0, z = 0.0;
	position(float tx, float ty, float tz) :x(tx),y(ty),z(tz){}
	position(){}
	position operator*(float p)
	{
		return position(x * p, y * p, z * p);
	}

	const position operator*(float p) const
	{
		return position(x * p, y * p, z * p);
	}

	position operator/(float p)
	{
		return position(x / p, y / p, z / p);
	}

	const position operator/(float p)const
	{
		return position(x / p, y / p, z / p);
	}

	float& operator[](int index) {
		switch (index) {
		case 0: return x;
		case 1: return y;
		case 2: return z;
		default: throw std::out_of_range("Index out of range: index must be 0, 1, or 2");
		}
	}
	const float& operator[](int index) const {
		switch (index) {
		case 0: return x;
		case 1: return y;
		case 2: return z;
		default: throw std::out_of_range("Index out of range: index must be 0, 1, or 2");
		}
	}
	// Vector subtraction
	position operator-(const position& other) const {
		return position(x - other.x, y - other.y, z - other.z);
	}
	position operator+(const position& other) const {
		return position(x + other.x, y + other.y, z + other.z);
	}
	// Cross product
	position cross(const position& other) const {
		return position(y * other.z - z * other.y,
			z * other.x - x * other.z,
			x * other.y - y * other.x);
	}

	// Normalize the vector
	void normalize() {
		float length = sqrt(x * x + y * y + z * z);
		if (length == 0) throw std::runtime_error("Zero length vector cannot be normalized.");
		x /= length;
		y /= length;
		z /= length;
	}

	float length() const {
		return std::sqrt(x * x + y * y + z * z);
	}


	float square_length() const {
		return x * x + y * y + z * z;
	}

	int mapCoord(float coord, int n) {

		float normalized = (coord + 0.5f) * n;

		int index = static_cast<int>(std::floor(normalized));

		index = std::clamp(index, 0, n - 1);

		return index;
	}
	VoxelIndices mapToVoxel(int n) {
		VoxelIndices voxel;
		voxel.i = mapCoord(x, n);
		voxel.j = mapCoord(y, n);
		voxel.k = mapCoord(z, n);
		return voxel;
	}

};

using velocity = position;
using acceleration = position;
using realvec3 = position;

struct Particle
{
	Particle(){}
	position pos;
	float density;
	float pressure;
	velocity velocity;
	acceleration acceleration;
};

class ParticleSystem
{
public:
	position get_pos(int id)
	{
		if (id < 0 || id > pos_vec.size())
		{
			std::cout << "get pos ff\n";
			return position(-1, -1, -1);
		}
		return pos_vec[id];
	}

	float get_density(int id)
	{
		if (id < 0 || id > density_vec.size())
		{
			std::cout << "get density ff\n";
			return -1;
		}
		return density_vec[id];
	}

	float get_pressure(int id)
	{
		if (id < 0 || id > pressure_vec.size())
		{
			std::cout << "get pressure ff\n";
			return -1;
		}
		return pressure_vec[id];
	}


	velocity get_velocity(int id)
	{
		if (id < 0 || id > velocity_vec.size())
		{
			std::cout << "get velocity ff\n";
			return velocity(-1, -1, -1);
		}
		return velocity_vec[id];
	}
	acceleration get_acceleration(int id)
	{
		if (id < 0 || id > acceleration_vec.size())
		{
			std::cout << "get acceleration ff\n";
			return acceleration(-1, -1, -1);
		}
		return acceleration_vec[id];
	}

	realvec3 get_fp(int id)
	{
		return fp_vec[id];
	}
	realvec3 get_fv(int id)
	{
		return fv_vec[id];
	}
	realvec3 get_fe(int id)
	{
		return fe_vec[id];
	}
	realvec3 get_fs(int id)
	{
		return fs_vec[id];
	}
	realvec3 get_n(int id)
	{
		return n_vec[id];
	}
	realvec3 get_f(int id)
	{
		return fp_vec[id] + fv_vec[id] + fe_vec[id] + fs_vec[id];
	}
	std::vector<position>& get_particles()
	{
		return pos_vec;
	}
	// -------------------------------------------------------------- 
	void set_pos(int id,const position& pos)
	{
		if (id < 0 || id > pos_vec.size())
		{
			std::cout << "set pos ff\n";
		}
		pos_vec[id] = pos;
	}

	void set_density(int id,const float& density)
	{
		if (id < 0 || id > density_vec.size())
		{
			std::cout << "set density ff\n";
		}
		density_vec[id] = density;
	}

	void set_pressure(int id,const float& pressure)
	{
		if (id < 0 || id > pressure_vec.size())
		{
			std::cout << "set pressure ff\n";
		}
		pressure_vec[id] = pressure;
	}


	void set_velocity(int id,const velocity& v)
	{
		if (id < 0 || id > velocity_vec.size())
		{
			std::cout << "set velocity ff\n";
		}
		velocity_vec[id] = v;
	}

	void set_acceleration(int id,const acceleration& a)
	{
		if (id < 0 || id > acceleration_vec.size())
		{
			std::cout << "set acceleration ff\n";
		}
		acceleration_vec[id] = a;
	}

	void set_fp(int id, const realvec3& fp)
	{
		if (id < 0 || id > fp_vec.size())
		{
			std::cout << "set fp ff\n";
		}
		fp_vec[id] = fp;
	}

	void set_fv(int id, const realvec3& fv)
	{
		if (id < 0 || id > fv_vec.size())
		{
			std::cout << "set fv ff\n";
		}
		fv_vec[id] = fv;
	}

	void set_fe(int id, const realvec3& fe)
	{
		if (id < 0 || id > fe_vec.size())
		{
			std::cout << "set fe ff\n";
		}
		fe_vec[id] = fe;
	}

	void set_fs(int id, const realvec3& fs)
	{
		if (id < 0 || id > fs_vec.size())
		{
			std::cout << "set fs ff\n";
		}
		fs_vec[id] = fs;
	}

	void set_n(int id, const realvec3& n)
	{
		if (id < 0 || id > n_vec.size())
		{
			std::cout << "set n ff\n";
		}
		n_vec[id] = n;
	}

	void InitSystem(float BR,float GC,float vis, float TC)
	{
		BoundRadius = BR;
		GasConstant = GC;
		Viscosity = vis;
		TensionCoefficient = TC;
	}

	void saveAsXYZ(const std::string& filename) {
		std::ofstream file(filename);
		if (!file.is_open()) {
			throw std::runtime_error("Failed to open file for writing.");
		}
		for (int id = 0; id < pos_vec.size(); ++id)
		{
			if (n_vec[id].length() < 0.0)continue;
			const position& pos = pos_vec[id];
			const position& n = n_vec[id];
			file << pos.x << " " << pos.y << " " << pos.z << " " << -n.x << " " << -n.y << " " << -n.z << "\n";
		}
		//for (const position& pos : pos_vec) {
		//	file << pos.x << " " << pos.y << " " << pos.z << "\n";
		//}

		file.close();
	}
	float get_BoundRadius()
	{
		return BoundRadius;
	}
	float get_GasConstant()
	{
		return GasConstant;
	}
	float get_Viscosity()
	{
		return Viscosity;
	}
	float get_TensionCoefficient()
	{
		return TensionCoefficient;
	}
	void add_particle(Particle& particle)
	{
		pos_vec.push_back(particle.pos);
		density_vec.push_back(particle.density);
		pressure_vec.push_back(particle.pressure);
		acceleration_vec.push_back(particle.acceleration);
		velocity_vec.push_back(particle.velocity);
		fp_vec.push_back({ 0, 0, 0 });
		fv_vec.push_back({ 0, 0, 0 });
		fe_vec.push_back({ 0, 0, 0 });
		fs_vec.push_back({ 0, 0, 0 });
		n_vec.push_back({ 0,0,0 });
	}

	int get_size()
	{
		return pos_vec.size();
	}

	void renew_grid(std::vector<std::vector<int>>& grid)
	{
		acc_grid = grid;
	}

	std::vector<int> getAdjacentPoints(int flatIndex, int n) {
		std::vector<int> adjacentPoints;
		adjacentPoints.insert(adjacentPoints.end(), acc_grid[flatIndex].begin(), acc_grid[flatIndex].end());
		return adjacentPoints;
		int i = flatIndex / (n * n);
		int j = (flatIndex % (n * n)) / n;
		int k = flatIndex % n;

		const int di[6] = { 1, -1, 0,  0, 0, 0 };
		const int dj[6] = { 0,  0, 1, -1, 0, 0 };
		const int dk[6] = { 0,  0, 0,  0, 1,-1 };

		for (int d = 0; d < 6; ++d) {
			int ni = i + di[d];
			int nj = j + dj[d];
			int nk = k + dk[d];

			if (ni >= 0 && ni < n && nj >= 0 && nj < n && nk >= 0 && nk < n) {
				int neighborFlatIndex = ni * n * n + nj * n + nk;
				const std::vector<int>& pointsInNeighbor = acc_grid[neighborFlatIndex];
				adjacentPoints.insert(adjacentPoints.end(), pointsInNeighbor.begin(), pointsInNeighbor.end());
			}
		}
		return adjacentPoints;
	}


private:
	std::vector<position> pos_vec;
	std::vector<float> density_vec;
	std::vector<float> pressure_vec;
	std::vector<velocity> velocity_vec;
	std::vector<acceleration> acceleration_vec;
	std::vector<realvec3> fp_vec;
	std::vector<realvec3> fv_vec;
	std::vector<realvec3> fe_vec;
	std::vector<realvec3> fs_vec;
	std::vector<realvec3> n_vec;
	std::vector<std::vector<int>> acc_grid;
	float BoundRadius = -1;
	float GasConstant = -1;
	float Viscosity = 1.0;
	float TensionCoefficient = 30.0;
};

struct tri
{
	int p0, p1, p2;
	tri(int a, int b, int c) :p0(a), p1(b), p2(c) {}
	int& operator[](int index) {
		switch (index) {
		case 0: return p0;
		case 1: return p1;
		case 2: return p2;
		default: throw std::out_of_range("Index out of range: index must be 0, 1, or 2");
		}
	}
	const int& operator[](int index) const {
		switch (index) {
		case 0: return p0;
		case 1: return p1;
		case 2: return p2;
		default: throw std::out_of_range("Index out of range: index must be 0, 1, or 2");
		}
	}
};

class Mesh
{
public:
	std::vector<tri> triangles;
	std::vector<position> points;
	std::vector<position> normals;
	std::vector<std::pair<position,position>> sample(int n) {
		std::vector<std::pair<position, position>> sampled_points;
		std::vector<float> areas(triangles.size());
		float total_area = 0.0;

		// Compute the area of each triangle and total area
		for (size_t i = 0; i < triangles.size(); ++i) {
			const auto& tri = triangles[i];
			const auto& p0 = points[tri[0]];
			const auto& p1 = points[tri[1]];
			const auto& p2 = points[tri[2]];

			// Calculate the area using the cross product (area of parallelogram) / 2
			float area = 0.5 * std::sqrt(
				std::pow((p1[1] - p0[1]) * (p2[2] - p0[2]) - (p2[1] - p0[1]) * (p1[2] - p0[2]), 2) +
				std::pow((p1[2] - p0[2]) * (p2[0] - p0[0]) - (p2[2] - p0[2]) * (p1[0] - p0[0]), 2) +
				std::pow((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]), 2)
			);
			areas[i] = area;
			total_area += area;
		}

		// Random sampling from each triangle based on area
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(0, 1);

		for (size_t i = 0; i < triangles.size(); ++i) {
			const auto& tri = triangles[i];
			const auto& p0 = points[tri[0]];
			const auto& p1 = points[tri[1]];
			const auto& p2 = points[tri[2]];
			int num_samples = std::round(areas[i] / total_area * n);  // Proportional to area

			for (int j = 0; j < num_samples; ++j) {
				float u = dis(gen);
				float v = dis(gen);

				// Make sure u + v <= 1; if not, reflect them inside the triangle
				if (u + v > 1.0) {
					u = 1.0 - u;
					v = 1.0 - v;
				}

				// Barycentric coordinates for the point inside the triangle
				float a = 1 - u - v;
				float b = u;
				float c = v;

				// The sampled point
				position sampled_point = {
					a * p0[0] + b * p1[0] + c * p2[0],
					a * p0[1] + b * p1[1] + c * p2[1],
					a * p0[2] + b * p1[2] + c * p2[2]
				};

				sampled_points.push_back(std::make_pair(sampled_point,normals[i]));
			}
		}

		return sampled_points;
	}
	void load_mesh(const std::string& path) {
		std::ifstream file(path);
		std::string line;

		if (!file.is_open()) {
			std::cerr << "Failed to open file: " << path << std::endl;
			return;
		}

		while (getline(file, line)) {
			std::istringstream iss(line);
			char type;
			iss >> type;
			if (type == 'v') {  // Vertex
				float x, y, z;
				iss >> x >> y >> z;
				points.push_back({ x, y, z });
			}
			else if (type == 'f') {  // Face
				int a, b, c;
				iss >> a >> b >> c;
				// .obj files are 1-based index, we convert them to 0-based for C++ arrays
				triangles.push_back({ a - 1, b - 1, c - 1 });
			}
		}

		for (int id = 0; id < triangles.size(); ++id) {
			normals.push_back(computeNormal(id));
		}


		file.close();
	}
	position computeNormal(int triangleIndex) {
		if (triangleIndex >= triangles.size())
			throw std::out_of_range("Triangle index is out of range.");

		const auto& tri = triangles[triangleIndex];
		position p1 = points[tri[0]];
		position p2 = points[tri[1]];
		position p3 = points[tri[2]];

		position v1 = p2 - p1;
		position v2 = p3 - p1;

		position normal = v1.cross(v2);
		normal.normalize();
		return normal;
	}
};

class WaterSource
{
public:
	Mesh sourcemesh;
	float initial_velocity;
	int sample_num = -1;
	int pour_deltaT = -1;
};