#include "bet2.h"
using namespace bet2;
int main(int argc, char **argv) {
	argv[1] = "C:/Users/jieji/Desktop/T2/20130610_144057T2AXTE80SENSEs301a1003.nii";
	argv[2] = "aaa.nii";
	cerr << argv[1] << '\n';
	cerr << argv[2] << '\n';
	volume<float> testvol;
	if (read_volume(testvol, argv[1]) < 0) {
		cerr << "error" << '\n';
		return -1;
	}
	double xarg = 0, yarg = 0, zarg = 0;
	const double bet_main_parameter = pow(fractional_threshold.value(), .275);
	// 2D kludge (worked for bet, but not here in bet2, hohum)
	if (testvol.xsize()*testvol.xdim()<20) testvol.setxdim(200);
	if (testvol.ysize()*testvol.ydim()<20) testvol.setydim(200);
	if (testvol.zsize()*testvol.zdim()<20) testvol.setzdim(200);

	Mesh m;
	make_mesh_from_icosa(5, m);
	bet_parameters bp = adjust_initial_mesh(testvol, m, radiusarg.value(), xarg, yarg, zarg);
	Mesh moriginal = m;

	const double rmin = 3.33 * smootharg.value();
	const double rmax = 10 * smootharg.value();
	const double E = (1 / rmin + 1 / rmax) / 2.;
	const double F = 6. / (1 / rmin - 1 / rmax);
	const int nb_iter = 1000;
	const double self_intersection_threshold = 4000;

	double l = 0;
	for (int i = 0; i<nb_iter; i++)
	{
		step_of_computation(testvol, m, bet_main_parameter, 0, 0, i, l, bp.t2, bp.tm, bp.t, E, F, bp.cog.Z, bp.radius, gradient_threshold.value());
	}

	double tmp = m.self_intersection(moriginal);
	if (verbose.value() && !generate_mesh.value())
		cout << "self-intersection total " << tmp << " (threshold=4000.0) " << endl;

	bool self_intersection;
	if (!generate_mesh.value()) self_intersection = (tmp > self_intersection_threshold);
	else (self_intersection = m.real_self_intersection());
	int pass = 0;

	if (verbose.value() && generate_mesh.value() && self_intersection)
		cout << "the mesh is self-intersecting " << endl;


	//self-intersection treatment
	while (self_intersection)
	{
		if (self_intersection && verbose.value()) { cout << "thus will rerun with higher smoothness constraint" << endl; };
		m = moriginal;
		l = 0;
		pass++;
		for (int i = 0; i<nb_iter; i++)
		{
			double incfactor = pow(10.0, (double)pass + 1);
			if (i > .75 * (double)nb_iter)
				incfactor = 4.*(1. - i / (double)nb_iter) * (incfactor - 1.) + 1.;
			step_of_computation(testvol, m, bet_main_parameter, pass, incfactor, i, l, bp.t2, bp.tm, bp.t, E, F, bp.cog.Z, bp.radius, gradient_threshold.value());
		}
		double tmp = m.self_intersection(moriginal);

		self_intersection = (tmp > self_intersection_threshold);
		if (!generate_mesh.value()) self_intersection = (tmp > self_intersection_threshold);
		else (self_intersection = m.real_self_intersection());

		if (verbose.value() && !generate_mesh.value())
			cout << "self-intersection total " << tmp << " (threshold=4000.0) " << endl;

		if (verbose.value() && generate_mesh.value() && self_intersection)
			cout << "the mesh is self-intersecting" << endl;

		if (pass == 10) // give up
			self_intersection = 0;
	}
	cerr << "adsfasd";

	//display
	volume<short> brainmask = make_mask_from_mesh(testvol, m);
	save_volume(static_cast<short>(1)-brainmask, argv[2]);
	cin.get();
	return 0;
}
