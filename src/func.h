#pragma once
#include"dfn.h"

void step(int tstep,const float& dt, WaterSource& source, ParticleSystem& PS,float h,int n);

void add_source(WaterSource& source, ParticleSystem& PS);

float poly6kernel(const float r2, const float h);

realvec3  spikykernelGrad(const position p, const float h);

float viscosityLaplacian(float r, float h);

void renew_density(ParticleSystem& PS,const float h,int n);

void renew_fp(ParticleSystem& PS, float h, int n);

void renew_fv(ParticleSystem& PS, float h, int n);

void renew_fe(ParticleSystem& PS);

void renew_fs(ParticleSystem& PS, float h, int n);

void renew_acceleration(ParticleSystem& PS);

void renew_velocity_position(ParticleSystem& PS, float dt);

realvec3 poly6KernelGradient(const realvec3& rVec, float h);

float poly6KernelLaplacian(float r, float h);

void renew_accgrid(ParticleSystem& PS, int n);

