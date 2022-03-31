#pragma once
#include "Dependencies/implicit_marching_cubes-master/mc.h"
#include <autodiff/forward/dual.hpp>

class Function_T
{

public:

	//std::function<double(Geometry::Vector3D)> sphere = [](Geometry::Vector3D p) {return sin(p[0] * p[2]) - p[1]; };

	Function_T(std::function<double(Geometry::Vector3D)> f, 
		std::function<autodiff::dual(autodiff::dual x, autodiff::dual y, autodiff::dual z)> f_dual,
		std::function<autodiff::dual2nd(autodiff::dual2nd x, autodiff::dual2nd y, autodiff::dual2nd z)> dual2nd
		) 
	{
	
		func = f;
		func_dual = f_dual;
		func_dual_2nd = dual2nd;
	}

	Function_T() {
	
		func = [](Geometry::Vector3D p) {return sin(p[0] * p[2]) - p[1]; };
	
	}

	std::function<double(Geometry::Vector3D)> func;
	std::function<autodiff::dual(autodiff::dual x, autodiff::dual y, autodiff::dual z)> func_dual;
	std::function<autodiff::dual2nd(autodiff::dual2nd x, autodiff::dual2nd y, autodiff::dual2nd z)> func_dual_2nd;


	double evaluate(Geometry::Vector3D at) {
	
		return func(at);
	
	}

	autodiff::dual evalute(Geometry::Vector3D at) {
	
		return func_dual(at[0],at[1],at[2]);
	
	}
	

};

