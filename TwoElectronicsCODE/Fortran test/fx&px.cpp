//第一个核，三个坐标的一阶导
__device__ double3 gx_gy_gz_first_nucleus(const nucleus& first, const nucleus& second)
{
	double Q_squre = pow(Q_constant, 2);
	//坐标平方和
	double loc_squre_sum_first = pow(first.x, 2) + pow(first.y, 2) + pow(first.z, 2);
	
	//一阶导平方和
	double px_py_pz_squre_sum_first = pow(first.px, 2) + pow(first.py, 2) + pow(first.pz, 2);



	//第一个核 一阶导三个计算公式 对应 g1 g3 g5
	double3 gx_gy_gz;
	gx_gy_gz.x = first.px * (1.0 - 1.0 / Q_squre * loc_squre_sum_first * px_py_pz_squre_sum_first
						* exp(A_hardness * (1.0 - pow(loc_squre_sum_first * px_py_pz_squre_sum_first/ Q_squre, 2))));

	gx_gy_gz.y = first.py * (1.0 - 1.0 / Q_squre * loc_squre_sum_first * px_py_pz_squre_sum_first
						* exp(A_hardness * (1.0 - pow(loc_squre_sum_first * px_py_pz_squre_sum_first / Q_squre, 2))));
	
	gx_gy_gz.z = first.pz * (1.0 - 1.0 / Q_squre * loc_squre_sum_first * px_py_pz_squre_sum_first
						* exp(A_hardness * (1.0 - pow(loc_squre_sum_first * px_py_pz_squre_sum_first / Q_squre, 2))));
	
	return gx_gy_gz;
}

//第二个核，三个坐标的一阶导
__device__ double3 gx_gy_gz_second_nucleus(const nucleus& first, const nucleus& second)
{
	const double Q_squre = pow(Q_constant, 2);
	
	//坐标平方和
	const double loc_squre_sum_second = pow(second.x, 2) + pow(second.y, 2) + pow(second.z, 2);

	//一阶导平方和
	const double px_py_pz_squre_sum_second = pow(second.px, 2) + pow(second.py, 2) + pow(second.pz, 2);



	//第二个核 一阶导三个计算公式 对应 g2 g4 g6
	double3 gx_gy_gz;
	gx_gy_gz.x = second.px * (1.0 - 1.0 / Q_squre * loc_squre_sum_second * px_py_pz_squre_sum_second
		* exp(A_hardness * (1.0 - pow(loc_squre_sum_second * px_py_pz_squre_sum_second / Q_squre, 2))));

	gx_gy_gz.y = second.py * (1.0 - 1.0 / Q_squre * loc_squre_sum_second * px_py_pz_squre_sum_second
		* exp(A_hardness * (1.0 - pow(loc_squre_sum_second * px_py_pz_squre_sum_second / Q_squre, 2))));

	gx_gy_gz.z = second.pz * (1.0 - 1.0 / Q_squre * loc_squre_sum_second * px_py_pz_squre_sum_second
		* exp(A_hardness * (1.0 - pow(loc_squre_sum_second * px_py_pz_squre_sum_second / Q_squre, 2))));

	return gx_gy_gz;
}

//第一个核，三个坐标的二阶导
__device__ double3 fx_fy_fz_first_nucleus(const nucleus& first, const nucleus& second)
{
	const double Q_squre = pow(Q_constant, 2);

	//坐标平方和
	const double loc_squre_sum_first = pow(first.x, 2) + pow(first.y, 2) + pow(first.z, 2);
	const double loc_squre_sum_second = pow(second.x, 2) + pow(second.y, 2) + pow(second.z, 2);
	//一阶导平方和
	const double px_py_pz_squre_sum_first = pow(first.px, 2) + pow(first.py, 2) + pow(first.pz, 2);
	const double px_py_pz_squre_sum_second = pow(second.px, 2) + pow(second.py, 2) + pow(second.pz, 2);
	//两核距离平方
	const double distance_squre = nucleus_distance(first, second);
	//两核距离平方的1.5 次方,对应 sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)
	const double distance_1_5_power = pow(distance_squre, 1.5);

	//坐标平方和的1.5 次方 ，对应 sqrt((z1**2.d0+y1**2.d0+x1**2.d0)**3.d0)
	const double loc_1_5_power_first = pow(loc_squre_sum_first, 1.5);
	const double loc_1_5_power_second = pow(loc_squre_sum_second, 1.5);

	//临时变量1 ,第一个粒子（一阶导平方和的平方 / Q方）
	const double temp1 = pow(px_py_pz_squre_sum_first, 2) / Q_squre;
	//临时变量2，第一个粒子（坐标平方和 * 一阶导平方和 / Q方）的平方
	const double temp2 = pow((loc_squre_sum_first * px_py_pz_squre_sum_first / Q_squre), 2);


	double3 fx_fy_fz;
	fx_fy_fz.x = first.x * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_first, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_first)
				+ (first.x - second.x) / distance_1_5_power;

	fx_fy_fz.y = first.y * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_first, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_first)
				+ (first.y - second.y) / distance_1_5_power;

	fx_fy_fz.z = first.z * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_first, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_first)
				+ (first.z - second.z) / distance_1_5_power;

	return  fx_fy_fz;
}


//第二个核，三个坐标的二阶导
__device__ double3 fx_fy_fz_second_nucleus(const nucleus& first, const nucleus& second)
{
	const double Q_squre = pow(Q_constant, 2);

	//坐标平方和
	const double loc_squre_sum_first = pow(first.x, 2) + pow(first.y, 2) + pow(first.z, 2);
	const double loc_squre_sum_second = pow(second.x, 2) + pow(second.y, 2) + pow(second.z, 2);
	//一阶导平方和
	const double px_py_pz_squre_sum_first = pow(first.px, 2) + pow(first.py, 2) + pow(first.pz, 2);
	const double px_py_pz_squre_sum_second = pow(second.px, 2) + pow(second.py, 2) + pow(second.pz, 2);
	//两核距离平方
	const double distance_squre = nucleus_distance(first, second);
	//两核距离平方的1.5 次方,对应 sqrt(((z1-z2)**2.d0+(y1-y2)**2.d0+(x1-x2)**2.d0)**3.d0)
	const double distance_1_5_power = pow(distance_squre, 1.5);

	//坐标平方和的1.5 次方 ，对应 sqrt((z1**2.d0+y1**2.d0+x1**2.d0)**3.d0)
	const double loc_1_5_power_first = pow(loc_squre_sum_first, 1.5);
	const double loc_1_5_power_second = pow(loc_squre_sum_second, 1.5);

	//临时变量1 ,第二个粒子（一阶导平方和的平方 / Q方）
	//对应(pz1**2.d0+px1**2.d0+py1**2.d0)**2.d0/q**2.d0
	const double temp1 = pow(px_py_pz_squre_sum_second, 2) / Q_squre;
	//临时变量2，第二个粒子（坐标平方和 * 一阶导平方和 / Q方）的平方
	const double temp2 = pow((loc_squre_sum_second * px_py_pz_squre_sum_second / Q_squre), 2);


	double3 fx_fy_fz;
	fx_fy_fz.x = second.x * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_second, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_second)
				- (first.x - second.x) / distance_1_5_power;

	fx_fy_fz.y = second.y * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_second, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_second)
				- (first.y - second.y) / distance_1_5_power;

	fx_fy_fz.z = second.z * ((Q_squre / 2.0 / A_hardness / pow(loc_squre_sum_second, 2) + temp1)
							* exp(A_hardness * (1.0 - temp2)) - 2.0 / loc_1_5_power_second)
				- (first.z - second.z) / distance_1_5_power;

	return  fx_fy_fz;

}



