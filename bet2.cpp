
/*  BET - Brain Extraction Tool

    BETv1 Steve Smith
    BETv2 Mickael Pechaud, Mark Jenkinson, Steve Smith
    FMRIB Image Analysis Group

    Copyright (C) 1999-2006 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/9564. */
#include "bet2.h"
using namespace bet2;
int main(int argc, char *argv[]) {

  //parsing options
  OptionParser options(title, examples);
  options.add(outline);
  options.add(generate_mask);
  options.add(skull);
  options.add(no_output);
  options.add(fractional_threshold);
  options.add(gradient_threshold);
  options.add(radiusarg);
  options.add(smootharg);

  options.add(centerarg);
  options.add(apply_thresholding);
  options.add(generate_mesh);
  options.add(verbose);
  options.add(help);
  
  if (argc < 3) {
    if (argc == 1) {options.usage(); exit(EXIT_FAILURE);};
    if (argc>1)
      {
	string s = argv[1];
	if (s.find("-h")!=string::npos | s.find("--help")!=string::npos ) 
	  {options.usage(); exit (0);}
      }
    cerr<<"error: not enough arguments, use bet -h for help"<<endl;
    exit (-1);
  }
  
  vector<string> strarg;
  for (int i=0; i<argc; i++)
    strarg.push_back(argv[i]);
  
  string inputname=strarg[1];
  string outputname=strarg[2];
  
  if (inputname.find("-")==0 || outputname.find("-")==0 )
    {cerr<<"error : two first arguments should be input and output names, see bet -h for help"<<endl; exit(-1);};
 
  /*
  int c=0;
  for (int i=0; i<argc; i++)
    if (i!=1 & i!=2)
      {
	strcpy(argv[c], strarg[i].c_str());
	c++;
      }
  
  argc -=2;
  */

  try {
    options.parse_command_line(argc, argv, 2);  // skip first 2 argv
  }
  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
  } 
  
  if (help.value()) {options.usage(); return 0;};

  string out = outputname;
  if (out.find(".hdr")!=string::npos) out.erase(out.find(".hdr"), 4);
  if (out.find(".img")!=string::npos) out.erase(out.find(".hdr"), 4);
  string in = inputname;
  if (in.find(".hdr")!=string::npos) in.erase(in.find(".hdr"), 4);
  if (in.find(".img")!=string::npos) in.erase(in.find(".hdr"), 4);
  if (out == "default__default") {out=in+"_brain";}
  

  //set a memory hanlder that displays an error message
  set_new_handler(noMoreMemory);


  //the real program

  volume<float> testvol;
  
  if (read_volume(testvol,in.c_str())<0)  return -1;

  double xarg = 0, yarg = 0, zarg = 0;
  if (centerarg.set())
    {
      /*
      if (centerarg.value().size()!=3)
	{
	  cerr<<"three parameters expected for center option !"<<endl;
	  cerr<<"please check there is no space after commas."<<endl;
	  exit (-1);
	}
      else
      */
	{
	  xarg = (double) centerarg.value(0);
	  yarg = (double) centerarg.value(1);
	  zarg = (double) centerarg.value(2);
	  ColumnVector v(4);
	  v << xarg << yarg << zarg << 1.0;
	  v = testvol.niftivox2newimagevox_mat() * v;
	  xarg = v(1);  yarg = v(2);  zarg = v(3);
	}
    }
  
  const double bet_main_parameter = pow(fractional_threshold.value(), .275);
  

  // 2D kludge (worked for bet, but not here in bet2, hohum)
  if (testvol.xsize()*testvol.xdim()<20) testvol.setxdim(200);
  if (testvol.ysize()*testvol.ydim()<20) testvol.setydim(200);
  if (testvol.zsize()*testvol.zdim()<20) testvol.setzdim(200);
  
  Mesh m;
  make_mesh_from_icosa(5, m); 
  
  
  bet_parameters bp = adjust_initial_mesh(testvol, m, radiusarg.value(), xarg, yarg, zarg);
  
  if (verbose.value())
    {
      cout<<"min "<<bp.min<<" thresh2 "<<bp.t2<<" thresh "<<bp.t<<" thresh98 "<<bp.t98<<" max "<<bp.max<<endl;
      cout<<"c-of-g "<<bp.cog.X<<" "<<bp.cog.Y<<" "<<bp.cog.Z<<" mm"<<endl;
      cout<<"radius "<<bp.radius<<" mm"<<endl;
      cout<<"median within-brain intensity "<<bp.tm<<endl;
    }
  

  Mesh moriginal=m;
  
  const double rmin=3.33 * smootharg.value();
  const double rmax=10 * smootharg.value();
  const double E = (1/rmin + 1/rmax)/2.;
  const double F = 6./(1/rmin - 1/rmax);
  const int nb_iter = 1000;
  const double self_intersection_threshold = 4000;
  
  double l = 0;
  for (int i=0; i<nb_iter; i++)
    {
      step_of_computation(testvol, m, bet_main_parameter, 0, 0, i, l, bp.t2, bp.tm, bp.t, E, F, bp.cog.Z, bp.radius, gradient_threshold.value());

      // display progress, revised by liangfu
      if ((i%int(nb_iter*.01f)) == 1){
        fprintf(stdout, "[%3.0f%%] BET: Performing %d of %d iterations.\t\t\r",
          i*100.f / nb_iter, i, nb_iter);
        fflush(stdout);
      }
    }
  
  double tmp = m.self_intersection(moriginal);
  if (verbose.value() && !generate_mesh.value())
    cout<<"self-intersection total "<<tmp<<" (threshold=4000.0) "<<endl;
  
  bool self_intersection;
  if (!generate_mesh.value()) self_intersection = (tmp > self_intersection_threshold);  
  else (self_intersection = m.real_self_intersection());
  int pass = 0;  

  if (verbose.value() && generate_mesh.value() && self_intersection)
    cout<<"the mesh is self-intersecting "<<endl;


  //self-intersection treatment
  while (self_intersection)
    {
      if (self_intersection && verbose.value()) {cout<<"thus will rerun with higher smoothness constraint"<<endl;};
      m = moriginal;
      l = 0;
      pass++;
      for (int i=0; i<nb_iter; i++)
	{
	  double incfactor = pow (10.0,(double) pass + 1);
	  if (i > .75 * (double)nb_iter)
	    incfactor = 4.*(1. - i/(double)nb_iter) * (incfactor - 1.) + 1.;
	  step_of_computation(testvol, m, bet_main_parameter, pass, incfactor, i, l, bp.t2, bp.tm, bp.t, E, F, bp.cog.Z, bp.radius, gradient_threshold.value());
	}
      double tmp = m.self_intersection(moriginal);
  
      self_intersection = (tmp > self_intersection_threshold);
      if (!generate_mesh.value()) self_intersection = (tmp > self_intersection_threshold);  
      else (self_intersection = m.real_self_intersection());
      
      if (verbose.value() && !generate_mesh.value())
	cout<<"self-intersection total "<<tmp<<" (threshold=4000.0) "<<endl;
      
      if (verbose.value() && generate_mesh.value() && self_intersection)
	cout<<"the mesh is self-intersecting"<<endl;
	
      if (pass==10) // give up
	self_intersection=0;
    }
  

  //display
  volume<short> brainmask = make_mask_from_mesh(testvol, m);
  
  if (apply_thresholding.value())
    {
      int xsize = testvol.xsize();
      int ysize = testvol.ysize();
      int zsize = testvol.zsize();
      for (int k=0; k<zsize; k++)
	for (int j=0; j<ysize; j++)
	  for (int i=0; i<xsize; i++)
	    if (testvol.value(i, j, k) < bp.t) brainmask.value(i, j, k) = 1;
    }
  
  if (!(no_output.value()))
    {
      int xsize = testvol.xsize();
      int ysize = testvol.ysize();
      int zsize = testvol.zsize();
      volume<float> output = testvol;
      for (int k=0; k<zsize; k++)
	for (int j=0; j<ysize; j++)
	  for (int i=0; i<xsize; i++)
	    output.value(i, j, k) = (1-brainmask.value(i, j, k)) * output.value(i, j, k);
      if (save_volume(output,out.c_str())<0)  return -1;
    }  
  
  if (generate_mask.value().size()>0){
      if (save_volume((short)1 - brainmask, generate_mask.value().c_str())<0)  return -1;
  }else{
    string maskstr = out + "_mask";
    if (save_volume((short)1 - brainmask, maskstr.c_str())<0)  return -1;
  }
  
  if (outline.value())
    {
      string outlinestr = out+"_overlay";
      volume<float> outline = draw_mesh_bis(testvol, m);
      if (save_volume(outline, outlinestr.c_str())<0)  return -1;
    }

  if (generate_mesh.value())
    {
      string meshstr = out+"_mesh.vtk";
      m.save(meshstr.c_str(),3);
    }

  if (skull.value())
    {
      string skullstr = out+"_skull";
      volume<float> skull = find_skull(testvol, m, bp.t2, bp.t, bp.t98);

      volume<short> bskull;
      copyconvert(skull,bskull);

      if (save_volume(bskull, skullstr.c_str())<0)  return -1;
    }

  return 0;
  
}

