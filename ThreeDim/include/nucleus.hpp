#ifndef NUCLEUS_HPP
#define	NUCLEUS_HPP



struct nucleus{
	double x;
	double y;
	double z;
	
	double px;
	double py;
	double pz;
};

struct  nuclei{
	//初始化
	nucleus init_first;
	nucleus init_second;
	//第一步完事
	nucleus step_one_first;
	nucleus step_one_second;
	//第二步完事
	nucleus step_two_first;
	nucleus step_two_second;
	
};






#endif //NUCLEUS_HPP