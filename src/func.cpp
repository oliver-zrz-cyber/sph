#include "func.h"

void step(int tstep,const float& dt, WaterSource& source, ParticleSystem& PS, float h, int n)
{
	std::string result_folder = "E:/research/SPH/MyImplementation/result/";
	// check source
	if (source.pour_deltaT > tstep)
		add_source(source, PS);

	renew_accgrid(PS, n);

	// renew density & so pressure at particles
	renew_density(PS, h, n);

	// compute acceleration
	renew_fe(PS);
	renew_fp(PS, h, n);
	renew_fv(PS, h, n);
	renew_fs(PS, h, n);
	renew_acceleration(PS);

	// renew v & pos
	renew_velocity_position(PS, dt);

	// save as pc
	//if (tstep % 10 == 0)
	//	PS.saveAsXYZ(result_folder + std::to_string(tstep) + ".xyz");

}

void add_source(WaterSource& source, ParticleSystem& PS)
{
	std::vector<std::pair<position,position>>sourceparticles = source.sourcemesh.sample(source.sample_num);
	for (int i = 0; i < sourceparticles.size(); ++i)
	{
		auto pos = sourceparticles[i].first;
		auto v = sourceparticles[i].second * source.initial_velocity;
		Particle pc;
		pc.density = 1;
		pc.pos = pos;
		pc.velocity = v;
		pc.pressure = pc.density * PS.get_GasConstant();
		pc.acceleration = acceleration(0,0,0);
		PS.add_particle(pc);
	}

}


float poly6kernel(const float r2, const float h)
{
	float h2 = h * h;
	if (r2 >= 0 && r2 <= h2) {
		float coef = 315.0f / (64.0f * 3.1415926 * pow(h, 9));
		return coef * pow(h2 - r2, 3);
	}
	return 0.0f;
}

realvec3 spikykernelGrad(const position ri_j, const float h)
{
	float r = ri_j.length();
	if (r > 0 && r <= h) {
		float coef = -45.0f / (3.1415926 * pow(h, 6));
		float term = pow(h - r, 2);
		return ri_j * (coef * term / r);
	}
	return realvec3{ 0, 0, 0 };
}

float viscosityLaplacian(float r, float h)
{
	if (r >= 0 && r <= h) {
		float coef = 45.0f / (3.1415926 * pow(h, 6));
		return coef * (h - r);
	}
	return 0.0f;
}

void renew_density(ParticleSystem& PS, const float h, int n)
{
	std::vector<float> new_density; 
	int size = PS.get_size();
	for (int i = 0; i < size; ++i)
	{
		float density = 0;

		if (n == -1) {
			for (int j = 0; j < size; ++j)
			{
				realvec3 ri = PS.get_pos(i);
				realvec3 rj = PS.get_pos(j);
				realvec3 ri_j = ri - rj;
				density += poly6kernel(ri_j.square_length(), h) * 1.0f;// *  mass
			}
		}
		else
		{
			position p = PS.get_pos(i);
			VoxelIndices vid= p.mapToVoxel(n);
			int flatindex = vid.voxelToIndex(n);
			std::vector<int>addpoints = PS.getAdjacentPoints(flatindex,n);
			realvec3 ri = PS.get_pos(i);
			for (int k = 0; k < addpoints.size(); ++k)
			{
				int jid = addpoints[k];
				realvec3 rj = PS.get_pos(jid);
				realvec3 ri_j = ri - rj;
				density += poly6kernel(ri_j.square_length(), h) * 1.0f;// *  mass
			}
		}

		new_density.push_back(density);
	}
	for (int k = 0; k < size; ++k)
	{
		PS.set_density(k, new_density[k]);
		PS.set_pressure(k, new_density[k] * PS.get_GasConstant());
	}
}


void renew_fp(ParticleSystem& PS, float h, int n)
{
	int size = PS.get_size();
	for (int i = 0; i < size; ++i) 
	{
		realvec3 ans(0, 0, 0);
		if (n == -1) {
			for (int j = 0; j < size; ++j)
			{
				float mj = 1.0f;
				float pi = PS.get_pressure(i);
				float pj = PS.get_pressure(j);
				float dj = PS.get_density(j);
				realvec3 ri = PS.get_pos(i);
				realvec3 rj = PS.get_pos(j);
				realvec3 ri_j = ri - rj;
				ans = ans - spikykernelGrad(ri_j, h) * mj * (pi + pj) / 2 / dj;
			}
		}
		else
		{
			position p = PS.get_pos(i);
			VoxelIndices vid = p.mapToVoxel(n);
			int flatindex = vid.voxelToIndex(n);
			std::vector<int>addpoints = PS.getAdjacentPoints(flatindex, n);
			realvec3 ri = PS.get_pos(i);
			float pi = PS.get_pressure(i);
			for (int k = 0; k < addpoints.size(); ++k)
			{
				int jid = addpoints[k];
				float mj = 1.0f;
				float pj = PS.get_pressure(jid);
				float dj = PS.get_density(jid);
				realvec3 rj = PS.get_pos(jid);
				realvec3 ri_j = ri - rj;
				ans = ans - spikykernelGrad(ri_j, h) * mj * (pi + pj) / 2 / dj;
			}
		}
		//std::cout << ans.length() << "fp\n";
		PS.set_fp(i,ans);
	}

}


void renew_fv(ParticleSystem& PS, float h, int n)
{
	int size = PS.get_size();
	float u = PS.get_Viscosity();
	for (int i = 0; i < size; ++i)
	{
		realvec3 ans(0, 0, 0);
		if (n == -1) {
			for (int j = 0; j < size; ++j)
			{
				float mj = 1.0f;
				realvec3 vi = PS.get_velocity(i);
				realvec3 vj = PS.get_velocity(j);
				float dj = PS.get_density(j);
				realvec3 ri = PS.get_pos(i);
				realvec3 rj = PS.get_pos(j);
				realvec3 ri_j = ri - rj;
				ans = ans + (vj - vi) * viscosityLaplacian(ri_j.length(), h) * mj * u / dj;
			}
		}
		else
		{
			position p = PS.get_pos(i);
			VoxelIndices vid = p.mapToVoxel(n);
			int flatindex = vid.voxelToIndex(n);
			std::vector<int>addpoints = PS.getAdjacentPoints(flatindex, n);

			realvec3 ri = PS.get_pos(i);
			realvec3 vi = PS.get_velocity(i);
			float pi = PS.get_pressure(i);

			for (int k = 0; k < addpoints.size(); ++k)
			{
				int jid = addpoints[k];
				float mj = 1.0f;
				realvec3 vj = PS.get_velocity(jid);
				float dj = PS.get_density(jid);
				realvec3 rj = PS.get_pos(jid);
				realvec3 ri_j = ri - rj;
				ans = ans + (vj - vi) * viscosityLaplacian(ri_j.length(), h) * mj * u / dj;
			}
		}
		PS.set_fv(i, ans);
	}
}

void renew_fe(ParticleSystem& PS)
{
	realvec3 g(0, 0, -9.8);
	int size = PS.get_size();
	for (int i = 0; i < size; ++i) {
		PS.set_fe(i, g * PS.get_density(i));
	}
}

void renew_fs(ParticleSystem& PS, float h, int n)
{
	int size = PS.get_size();
	for (int i = 0; i < size; ++i)
	{
		realvec3 fs(0, 0, 0);
		realvec3 normal(0, 0, 0);
		float fsc = 0.0;
		if (n == -1) 
		{
			for (int j = 0; j < size; ++j)
			{
				float mj = 1.0f;
				float dj = PS.get_density(j);
				realvec3 ri = PS.get_pos(i);
				realvec3 rj = PS.get_pos(j);
				realvec3 ri_j = ri - rj;
				normal = normal + poly6KernelGradient(ri_j, h) * mj / dj;
				fsc += poly6KernelLaplacian(ri_j.length(), h) * mj / dj;
			}
		}
		else
		{
			position p = PS.get_pos(i);
			VoxelIndices vid = p.mapToVoxel(n);
			int flatindex = vid.voxelToIndex(n);
			std::vector<int>addpoints = PS.getAdjacentPoints(flatindex, n);
			realvec3 ri = PS.get_pos(i);
			realvec3 vi = PS.get_velocity(i);
			float pi = PS.get_pressure(i);
			for (int k = 0; k < addpoints.size(); ++k)
			{
				int jid = addpoints[k];
				float mj = 1.0f;
				float dj = PS.get_density(jid);
				realvec3 rj = PS.get_pos(jid);
				realvec3 ri_j = ri - rj;
				normal = normal + poly6KernelGradient(ri_j, h) * mj / dj;
				fsc += poly6KernelLaplacian(ri_j.length(), h) * mj / dj;
			}
		}
		PS.set_n(i, normal);
		if (normal.length() < 3.0)PS.set_fs(i, { 0,0,0 });
		else {
			fs = normal * -1.0f * PS.get_TensionCoefficient() * fsc / normal.length();
			PS.set_fs(i, fs);
		}
	}
}

void renew_acceleration(ParticleSystem& PS)
{
	int size = PS.get_size();
	for (int i = 0; i < size; ++i) 
	{
		realvec3 a = PS.get_f(i) / PS.get_density(i);
		PS.set_acceleration(i, a);
	}
}

void renew_velocity_position(ParticleSystem& PS, float dt)
{
	int size = PS.get_size();
	float Br= PS.get_BoundRadius();
	for (int i = 0; i < size; ++i)
	{
		realvec3 a = PS.get_acceleration(i);
		realvec3 prev = PS.get_velocity(i) * 0.998;
		realvec3 newv = a * dt + prev;
		position prepos = PS.get_pos(i);
		position newpos = prepos + newv * dt;

		//check boundary 
		float rold = std::sqrtf(prepos.x * prepos.x + prepos.y * prepos.y);
		float rnew = std::sqrtf(newpos.x * newpos.x + newpos.y * newpos.y);
		bool levelcollision = false;
		if (rnew > Br - 0.01)
		{
			levelcollision = true;
			float prerest = Br - rold;
			float radd = rnew - rold;
			float ratio = prerest / radd;
			newpos = prepos ;
			//rnew = std::sqrtf(newpos.x * newpos.x + newpos.y * newpos.y);
			//if (rnew > Br) {
			//	std::cout << "boundary ff\n";
			//}
		}
		bool verticalcollision = false;
		if (newpos.z < 0)
		{
			verticalcollision = true;
			newpos.z = 0;
		}
		PS.set_pos(i, newpos);
		if (levelcollision) { // compute velocity
			newv = {-newv.x,-newv.y,newv.z};
		}
		if (verticalcollision) {
			newv.z = -newv.z * 0.8;
		}
		PS.set_velocity(i, newv);
	}
}

realvec3 poly6KernelGradient(const realvec3& rVec, float h)
{
	const float rx = rVec[0];
	const float ry = rVec[1];
	const float rz = rVec[2];
	const float r2 = rx * rx + ry * ry + rz * rz;
	const float r = std::sqrt(r2);

	if (r >= h) {
		return { 0.0f, 0.0f, 0.0f };
	}

	static const float invPi = 1.0f / float(3.1415926);
	const float h2 = h * h;
	const float term2 = (h2 - r2);

	const float factor = -945 * invPi / (32.0f * std::pow(h, 9)) * term2 * term2;

	return { factor * rx, factor * ry, factor * rz };
}

float poly6KernelLaplacian(float r, float h)
{
	if (r >= h) {
		return 0.0f;
	}
	static const float invPi = 1.0f / float(3.1415926);

	const float h2 = h * h;
	const float r2 = r * r;

	const float lap = -945.0f * invPi / (32.0f * std::pow(h, 9))
		* (h2 - r2) * (3.0f * h2 - 7.0f * r2);
	return lap;
}

void renew_accgrid(ParticleSystem& PS, int n)
{
	if (n == -1)return;
	std::vector<std::vector<int>> voxelGrid(n * n * n);
	int size = PS.get_size();
//#pragma omp parallel for
	for (int id = 0; id < size; ++id)
	{
		position p = PS.get_pos(id);
		VoxelIndices voxel = p.mapToVoxel(n);
		int flatindex = voxel.voxelToIndex(n);
		voxelGrid[flatindex].push_back(static_cast<int>(id));
	}
	PS.renew_grid(voxelGrid);
}
