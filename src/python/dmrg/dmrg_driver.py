import calendar
import time
import datetime

class DMRG_Driver:

  def __init__( self ):

    import hydrogen_ext_hubbard;
    self.hydrogen_ext = hydrogen_ext_hubbard.Hydrogen_Extended_Hubbard();
    self.date         = datetime.datetime.now();
    self.seed         = calendar.timegm( self.date.timetuple() );
    self.model        = "H";
    self.M            = 100;
    self.tii          = 2.0;
    self.tij          = -1.0;
    self.u            = 4.0;
    self.h            = 0.0e0;
    self.norb         = 8;
    self.nelec        = self.norb;
    self.mult         = 0;
    self.point_group  = "c1";
    self.sweep_tol    = 1.0e-6;
    self.irrep        = 1;
    self.nroot        = 1;
    self.twodot_to_onedot = 100;
    self.maxiter          = 20;
    self.stochastic_ratio = 0.5e0;
    #self.stochastic_sweep_end = self.maxiter/3;
    self.stochastic_sweep_end = 8;
    self.prefix               = "./"
    self.dot_algorithm        = "twodot";
    self.onepdm               = False;
    self.twopdm               = False;
    self.job_name      = "dmrg.job." + self.date.strftime( "%Y-%m-%d-%H-%M-%S" );
    self.integral_name = self.job_name + ".qcdmrg.FCIDUMP";
    self.integral_path = "";
    self.alpha         = 1.0e0;
    self.pbs_script    = False;
    self.nersc	       = False;

    self.nnode         = 1;
    self.np            = 1;
    self.ppn           = 16;
    self.home          = "/home/";
    self.current_path  = "";
    self.input_name    = self.job_name + ".conf";
    self.input_path    = "";
    self.scratch       = "/scratch/"
    self.nrow          = 3;
    self.ncol          = 8;


  def print_obj( self ):

    print "date:\t\t\t",         self.date;
    print "seed:\t\t\t",         self.seed;
    print "M:\t\t\t",            self.M;
    print "tii:\t\t\t",          self.tii;
    print "tij:\t\t\t",          self.tij;
    print "u:\t\t\t",            self.u;
    print "norb:\t\t\t",         self.norb;
    print "nelec:\t\t\t",        self.nelec;
    print "mult:\t\t\t",         self.mult;
    print "point_group:\t\t",    self.point_group;
    print "sweep_tol:\t\t",      self.sweep_tol;
    print "irrep:\t\t\t",        self.irrep;
    print "nroot:\t\t\t",        self.nroot;
    print "twodot_to_onedot:\t", self.twodot_to_onedot;
    print "maxiter:\t\t",        self.maxiter;
    print "stochastic_ratio:\t", self.stochastic_ratio;
    print "stochastic_last_sweep:\t", self.stochastic_sweep_end;
    print "prefix:\t\t\t",       self.prefix;
    print "dot_algorithm:\t\t",  self.dot_algorithm;
    print "job_name:\t\t",       self.job_name;
    print "integral_name:\t\t",  self.integral_name;
    print "integral_path:\t\t",  self.integral_path;
    print "alpha(for ext. hubbard\t\t", self.alpha;
    print "pbs_script:\t\t\t",   self.pbs_script;
    print "np:\t\t\t",           self.np;
    print "current path:",       self.current_path;


  def make_working_dir( self ):

    import os;
    self.current_path = os.getcwd();
    self.prefix  = self.current_path + "/" + self.job_name;
    import subprocess;
    mkdir_command_line = "mkdir -p " + self.prefix;
    subprocess.check_call( mkdir_command_line, shell=True );


  def run( self ):
    import subprocess;
    self.make_working_dir();
    if self.model == "H":
      self.write_hubbard_int();
    elif self.model == "LH":
      self.write_ladder_hubbard_int();
    elif self.model == "Q2H":
      self.write_quasi_2d_hubbard_int();
    elif self.model == "EXTH":
      self.write_extended_hubbard_int();
    elif self.model == "2DH":
      self.write_2d_hubbard_int();
    elif self.model == "H_PPP":
      self.write_extended_hubbard_hydrogen_chain_int();
    elif self.model == "QC":
      self.write_qc_int();

    if self.pbs_script == False:
      self.output_file = self.prefix + "/" + self.job_name + "." + "out";
      if self.np != 1:
        #command_line = "mpirun -n " + str(self.np) + " dmrg.stochastic " + self.write_input() + " > " + self.output_file;
        command_line = "mpirun -n " + str(self.np) + " " + "/" + self.home + "/weifeng/local/dmrg/src/c++/dmrg.stochastic " + self.write_input() + " > " + self.output_file;
      else: 
        command_line = "/" + self.home + "/weifeng/local/dmrg/src/c++/dmrg.stochastic " + self.write_input() + " > " + self.output_file;
      subprocess.check_call( command_line, shell=True );

    else:

      pbs_path        = "/" + self.prefix + "/" + self.job_name + "." + "sh";
      pbs_input_path  = "/" + self.scratch + "/weifeng/" + self.job_name + "/" + self.input_name;
      pbs_output_path = "/" + self.scratch + "/weifeng/" + self.job_name + "/" + self.job_name + ".out";

      pbs_current_path = self.current_path;

      if self.nersc == True:
# SRUN for cori
        f_pbs = open( pbs_path, "wt" );
        f_pbs.write( "#!/bin/bash\n" );
        f_pbs.write( "#SBATCH -N " + str(self.nnode) + "\n" );
        f_pbs.write( "#SBATCH -p regular\n"  );
        f_pbs.write( "#SBATCH -L SCRATCH\n"  );
        f_pbs.write( "#SBATCH -t 48:00:00\n" );
        f_pbs.write( "#SBATCH -C haswell\n"  );
        f_pbs.write( "#SBATCH -J " + self.job_name + "\n" );

# These are not that useful given that the openmpi is compiled by myself
#   f_pbs.write( "module load impi\n" );
#   f_pbs.write( "export I_MPI_FABRICS=ofi\n" );
#   f_pbs.write( "export I_MPI_OFI_PROVIDER=gni\n" );
#   f_pbs.write( "export I_MPI_OFI_LIBRARY=/usr/common/software/libfabric/1.5.0/gnu/lib/libfabric.so\n" );
#   f_pbs.write( "export I_MPI_PMI_LIBRARY=/usr/lib64/slurmpmi/libpmi.so\n" );
#   f_pbs.write( "\n" );
#   f_pbs.write( "export LD_LIBRARY_PATH=" + "/" + self.home + "/weifeng/local/lib:$LD_LIBRARY_PATH\n" );
#   f_pbs.write( "export LD_LIBRARY_PATH=" + "/" + self.home + "/weifeng/local/lib64:$LD_LIBRARY_PATH\n" );
#   f_pbs.write( "\n" );

# We also don't need these since the NERSC machines have very small home dirs so copying 10000+ files into it then post-processing is not a good option
# just directly put them into the cscratch1
#   f_pbs.write( "mkdir -p " + "/" + self.scratch + "/weifeng/\n" );
#   f_pbs.write( "mkdir -p " + "/" + self.scratch + "/weifeng/" + self.job_name + "\n" );
#   f_pbs.write( "cp -r " + "/" + self.home + "/weifeng/jobs/" + self.job_name + "/* " + "/" + self.scratch + "/weifeng/" + self.job_name + "\n" );
#   f_pbs.write( "\n" );

# We always use the binary in local/bin, need to make sure that the binary is up-to-date
        f_pbs.write( "mpirun -n " + str(self.np) + " " + "/" + self.home + "/weifeng/local/bin/dmrg.stochastic " + pbs_input_path + " > " + pbs_output_path + "\n" );

# We don't delete anything if we run jobs in the cscratch1
#   f_pbs.write( "rm -rf " + self.scratch + "/weifeng/" + self.job_name + "/*.tmp" + "\n" );
#   f_pbs.write( "cp " + self.scratch + "/weifeng/" + self.job_name + "/* " + "/" + self.home + "/weifeng/jobs/" + self.job_name + "/" + "\n" );
#   f_pbs.write( "rm -rf " + self.scratch + "/weifeng/" + self.job_name + "\n" );

        start_path = "/" + self.current_path + "/" + "job_start.sh";
        f_start = open( start_path, "a" );
        f_start.write( "#!/bin/bash\n" );
#  f_start.write( "cp -r "  + self.job_name + "/ " + " " + "/" + self.home + "/weifeng/jobs/\n" );
#  f_start.write( "sbatch " + "/" + self.home + "/weifeng/jobs/" + self.job_name + "/" + self.job_name + ".sh" + "\n" );
        f_start.write( "cd " + self.job_name + "\n" );
        f_start.write( "sh " + self.job_name + ".sh\n");
        f_start.write( "cd ..;\n ");
        f_start.close();
     
#  end_path = self.current_path + "/" + "job_end.sh";
#  f_end = open( end_path, "a" );
#  f_end.write( "#!/bin/bash\n" );
#  f_end.write( "cp -r "  + "/" + self.home + "/weifeng/jobs/" + self.job_name + "/" + " ./\n" );
#  f_end.write( "rm -rf " + "/" + self.home + "/weifeng/jobs/" + self.job_name + "/\n" );
#  f_end.close();

      else:
        f_pbs = open( pbs_path, "wt" );
        f_pbs.write( "#!/bin/bash\n" );
        f_pbs.write( "#PBS -N " + self.job_name + "\n" );
        f_pbs.write( "#PBS -l nodes=" + str(self.np) + ":ppn=" + str(self.ppn) + "\n" );
        f_pbs.write( "\n" );
        f_pbs.write( "export LD_LIBRARY_PATH=/" + self.home + "/weifeng/local/lib64:$LD_LIBRARY_PATH\n" );
        f_pbs.write( "\n" );
        f_pbs.write( "mkdir -p " + "/" + self.scratch + "/weifeng/\n" );
        f_pbs.write( "mkdir -p " + "/" + self.scratch + "/weifeng/"      + self.job_name + "\n" );
        f_pbs.write( "cp -r "    + "/" + self.home    + "/weifeng/jobs/" + self.job_name + "/* " + "/" + self.scratch + "/weifeng/" + self.job_name + "\n" );
        f_pbs.write( "\n" );
 
        if self.np != 1:
          f_pbs.write( "mpirun -n " + str(self.np) + " /" + str( self.home ) + "/weifeng/local/bin/dmrg.stochastic " + pbs_input_path + " > " + pbs_output_path + "\n" );
        else:
          f_pbs.write( "/" + self.home + "/weifeng/local/bin/dmrg.stochastic " + pbs_input_path + " > " + pbs_output_path + "\n" );
 
        f_pbs.write( "cp "     + "/" + self.scratch + "/weifeng/" + self.job_name + "/* " + "/" + self.home + "/weifeng/jobs/" + self.job_name + "/" + "\n" );
        f_pbs.write( "rm -rf " + "/" + self.scratch + "/weifeng/" + self.job_name + "\n" );
        f_pbs.close();
        pbs_job_line = "echo qsub " + pbs_path + " >> pbs_script.sh ";
        subprocess.check_call( pbs_job_line, shell = True );
  
        start_path = "/" + self.current_path + "/" + "job_start.sh";
        f_start = open( start_path, "a" );
        f_start.write( "#!/bin/bash\n" );
        f_start.write( "cp -r " + self.job_name + "/ " + " " + "/" + self.home + "/weifeng/jobs/\n" );
        f_start.write( "qsub " + "/" + self.home + "/weifeng/jobs/" + self.job_name + "/" + self.job_name + ".sh" + "\n" );
        f_start.close();
  
        end_path = "/" + self.current_path + "/" + "job_end.sh";
        f_end = open( end_path, "a" );
        f_end.write( "#!/bin/bash\n" );
        f_end.write( "cp -r "  + "/" + self.home + "/weifeng/jobs/" + self.job_name + "/" + " ./\n" );
        f_end.write( "rm -rf " + "/" + self.home + "/weifeng/jobs/" + self.job_name + "/\n" );
        f_end.close();
 
      self.write_input();


  def write_hubbard_int( self ):
    self.integral_path = "/" + self.prefix + "/" + self.integral_name;
    f_integral = open( self.integral_path, "wt" );
    f_integral.write( " &FCI  NORB= " + str( self.norb ) + ",NELEC= " + str( self.nelec ) + ",MS2= " + str( self.mult) + ",\n" );
    f_integral.write( "  ORBSYM=" );
    for i in range( 0, self.norb ):
      f_integral.write( "1," );
    f_integral.write( "\n" );
    f_integral.write( "  ISYM=1\n");
    f_integral.write( " &END\n" );
    for i in range( 0, self.norb ):
      f_integral.write( str( self.u ) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(i+1) + "\n" );
    for i in range( 0, self.norb ):
      f_integral.write( str( self.tii ) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(0) + "    " + str(0) + "\n" );
    for i in range( 0, self.norb - 1 ):
      f_integral.write( str( self.tij ) + "    " + str(i+1) + "    " + str(i+2) + "    " + str(0) + "    " + str(0) + "\n" );
      f_integral.write( str( self.tij ) + "    " + str(i+2) + "    " + str(i+1) + "    " + str(0) + "    " + str(0) + "\n" );
    f_integral.write( str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + "\n" );
    f_integral.close();


  def write_ladder_hubbard_int( self ):
    self.integral_path = "/" + self.prefix + "/" + self.integral_name;
    f_integral = open( self.integral_path, "wt" );
    f_integral.write( " &FCI  NORB= " + str( self.norb ) + ",NELEC= " + str( self.nelec ) + ",MS2= " + str( self.mult) + ",\n" );
    f_integral.write( "  ORBSYM=" );
    for i in range( 0, self.norb ):
      f_integral.write( "1," );
    f_integral.write( "\n" );
    f_integral.write( "  ISYM=1\n");
    f_integral.write( " &END\n" );
    for i in range( 0, self.norb ):
      f_integral.write( str( self.u ) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(i+1) + "\n" );
    for i in range( 0, self.norb ):
      f_integral.write( str( self.tii ) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(0) + "    " + str(0) + "\n" );

    for i in range( 0, self.norb - 2, 2 ):
      f_integral.write( str( self.tij ) + "    " + str(i+1) + "    " + str(i+3) + "    " + str(0) + "    " + str(0) + "\n" );
      f_integral.write( str( self.tij ) + "    " + str(i+3) + "    " + str(i+1) + "    " + str(0) + "    " + str(0) + "\n" );
    for i in range( 1, self.norb - 2, 2 ):
      f_integral.write( str( self.tij ) + "    " + str(i+1) + "    " + str(i+3) + "    " + str(0) + "    " + str(0) + "\n" );
      f_integral.write( str( self.tij ) + "    " + str(i+3) + "    " + str(i+1) + "    " + str(0) + "    " + str(0) + "\n" );
    for i in range( 0, self.norb - 1, 2 ):
      f_integral.write( str( self.tij ) + "    " + str(i+1) + "    " + str(i+2) + "    " + str(0) + "    " + str(0) + "\n" );
      f_integral.write( str( self.tij ) + "    " + str(i+2) + "    " + str(i+1) + "    " + str(0) + "    " + str(0) + "\n" );

    f_integral.write( str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + "\n" );
    f_integral.close();


  def write_quasi_2d_hubbard_int( self ):
    norb = self.norb;
    import math;

    self.integral_path = "/" + self.prefix + "/" + self.integral_name;
    f_integral = open( self.integral_path, "wt" );
    f_integral.write( " &FCI  NORB= " + str( self.norb ) + ",NELEC= " + str( self.nelec ) + ",MS2= " + str( self.mult) + ",\n" );
    f_integral.write( "  ORBSYM=" );
    for i in range( 0, self.norb ):
      f_integral.write( "1," );
    f_integral.write( "\n" );
    f_integral.write( "  ISYM=1\n");
    f_integral.write( " &END\n" );
    for i in range( 0, self.norb ):
      f_integral.write( str( self.u ) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(i+1) + "\n" );

    for i in range( 0, self.norb ):
      f_integral.write( str( self.tii ) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(0) + "    " + str(0) + "\n" );

    # now print the hopping integrals
    ncomplete = 0;
    m = self.nrow;
    n = self.ncol;
    if m * n != norb:
      print "nrow x ncol must equal norb";
      quit();
    if m == 2:
      print "this quasi 2d does not work for 2 x l ladders";
      quit();

    arr = [];
    for i in range( 0, m ):
      subarr = [];
      for j in range( 0, n ):
        subarr.append(0);
      arr.append( subarr );
        
    i = 0;
    j = 0;
    ncomplete = 1;
    while ncomplete <= norb and i < m and j < n:
#      print i, j, ncomplete
      if i == 0 or j == n - 1:
        if i == 0 and j < n - 1:
          arr[i][j] = ncomplete;
          j = j + 1;
          ncomplete += 1;
          if ncomplete > norb :
            break; 
        elif j == n - 1 and i < m - 1:
          arr[i][j] = ncomplete;
          i = i + 1;
          ncomplete += 1;
          if ncomplete > norb :
            break;
        else:
          arr[i][j] = ncomplete;
          ncomplete += 1;
        while i < m - 1 and j > 0 and ncomplete <= norb:
          arr[i][j] = ncomplete;
          i += 1; j -= 1;
          ncomplete += 1;

      elif j == 0 or i == m - 1:
        if j == 0 and i < m - 1:
          arr[i][j] = ncomplete;
          i += 1;
          ncomplete += 1;
          if ncomplete > norb:
            break;
        elif i == m - 1 and j < n - 1:
          arr[i][j] = ncomplete;
          j += 1;
          ncomplete += 1;
          if ncomplete > norb:
            break;
        else:
          arr[i][j] = ncomplete;
          ncomplete += 1;
        while i > 0 and j < n - 1 and ncomplete <= norb:
          arr[i][j] = ncomplete;
          i -= 1; j += 1;
          ncomplete += 1;

    for i in range(0, m):
      for j in range(0, n):
        if i + 1 < m:
          f_integral.write( str( self.tij ) + "    " + str( arr[i][j] ) + "    " + str( arr[i+1][j] ) + "    " + str(0) + "    " + str(0) + "\n" );
          f_integral.write( str( self.tij ) + "    " + str( arr[i+1][j] ) + "    " + str( arr[i][j] ) + "    " + str(0) + "    " + str(0) + "\n" );
        if j + 1 < n:
          f_integral.write( str( self.tij ) + "    " + str( arr[i][j] ) + "    " + str( arr[i][j+1] ) + "    " + str(0) + "    " + str(0) + "\n" );
          f_integral.write( str( self.tij ) + "    " + str( arr[i][j+1] ) + "    " + str( arr[i][j] ) + "    " + str(0) + "    " + str(0) + "\n" );
    f_integral.close();


  def write_2d_hubbard_int( self ):
    norb = self.norb;
    import math;
    sqrt_root_of_norb = int( math.sqrt(norb) );
    self.norb = sqrt_root_of_norb * sqrt_root_of_norb;

    self.integral_path = "/" + self.prefix + "/" + self.integral_name;
    f_integral = open( self.integral_path, "wt" );
    f_integral.write( " &FCI  NORB= " + str( self.norb ) + ",NELEC= " + str( self.nelec ) + ",MS2= " + str( self.mult) + ",\n" );
    f_integral.write( "  ORBSYM=" );
    for i in range( 0, self.norb ):
      f_integral.write( "1," );
    f_integral.write( "\n" );
    f_integral.write( "  ISYM=1\n");
    f_integral.write( " &END\n" );
    for i in range( 0, self.norb ):
      f_integral.write( str( self.u ) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(i+1) + "\n" );

    for i in range( 0, self.norb ):
      f_integral.write( str( self.tii ) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(0) + "    " + str(0) + "\n" );

    # a small algorithm problem DP
    ind_group = [];
    n_layer = 2 * sqrt_root_of_norb - 1;
    forward = False;
    increase = True;
    current_length = 0;
    start = 0;
    for i_layer in range( 0, n_layer ):
      start += current_length ;
      if i_layer  < sqrt_root_of_norb :
        current_length += 1;
      else:
        current_length -= 1;

      #print i_layer, "-" , start, "--", current_length
      ind_subgroup = [];
      if forward == True:
        for i in range( 0, current_length ):
          new_ind = start + i + 1;
          ind_subgroup.append( new_ind );

          if i_layer < sqrt_root_of_norb:
            if i < ( current_length - 1 ):
              prev_i = i;
              f_integral.write( str( self.tij ) + "    " + str( new_ind ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str(0) + "    " + str(0) + "\n" );
              f_integral.write( str( self.tij ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str( new_ind ) + "    " + str(0) + "    " + str(0) + "\n" );

            if i > 0:
              prev_i = i - 1;
              f_integral.write( str( self.tij ) + "    " + str( new_ind ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str(0) + "    " + str(0) + "\n" );
              f_integral.write( str( self.tij ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str( new_ind ) + "    " + str(0) + "    " + str(0) + "\n" );

          else:
            prev_i = i + 1;
            f_integral.write( str( self.tij ) + "    " + str( new_ind ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str(0) + "    " + str(0) + "\n" );
            f_integral.write( str( self.tij ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str( new_ind ) + "    " + str(0) + "    " + str(0) + "\n" );
            prev_i = i;
            f_integral.write( str( self.tij ) + "    " + str( new_ind ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str(0) + "    " + str(0) + "\n" );
            f_integral.write( str( self.tij ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str( new_ind ) + "    " + str(0) + "    " + str(0) + "\n" );

#          print str( new_ind ) + ", ",

        forward = False;

      else:
        for i in range( 0, current_length ):
          new_ind = start + current_length - i;
          ind_subgroup.append( new_ind );

          if i_layer < sqrt_root_of_norb:
            if i < ( current_length - 1 ):
              prev_i = i;
              f_integral.write( str( self.tij ) + "    " + str( new_ind ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str(0) + "    " + str(0) + "\n" );
              f_integral.write( str( self.tij ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str( new_ind ) + "    " + str(0) + "    " + str(0) + "\n" );

            if i > 0:
              prev_i = i - 1;
              f_integral.write( str( self.tij ) + "    " + str( new_ind ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str(0) + "    " + str(0) + "\n" );
              f_integral.write( str( self.tij ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str( new_ind ) + "    " + str(0) + "    " + str(0) + "\n" );

          else:
            prev_i = i + 1;
            f_integral.write( str( self.tij ) + "    " + str( new_ind ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str(0) + "    " + str(0) + "\n" );
            f_integral.write( str( self.tij ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str( new_ind ) + "    " + str(0) + "    " + str(0) + "\n" );
            prev_i = i;
            f_integral.write( str( self.tij ) + "    " + str( new_ind ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str(0) + "    " + str(0) + "\n" );
            f_integral.write( str( self.tij ) + "    " + str( ind_group[i_layer-1][prev_i] ) + "    " + str( new_ind ) + "    " + str(0) + "    " + str(0) + "\n" );


#          print str( new_ind ) + ", ",
        forward = True;

      ind_group.append( ind_subgroup );
#      print 

    f_integral.write( str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + "\n" );
    f_integral.close();


  def write_extended_hubbard_int( self ):
    self.integral_path = "/" + self.prefix + "/" + self.integral_name;
    f_integral = open( self.integral_path, "wt" );
    f_integral.write( " &FCI  NORB= " + str( self.norb ) + ",NELEC= " + str( self.nelec ) + ",MS2= " + str( self.mult) + ",\n" );
    f_integral.write( "  ORBSYM=" );
    for i in range( 0, self.norb ):
      f_integral.write( "1," );
    f_integral.write( "\n" );
    f_integral.write( "  ISYM=1\n");
    f_integral.write( " &END\n" );
    for i in range( 0, self.norb ):
      f_integral.write( str( self.u ) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(i+1) + "\n" );
      for j in range( i + 1, self.norb ):
        dist = ( j - i ) * 1.0e0;
        value = self.u/( dist * self.alpha );
        f_integral.write( str( value ) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(j+1) + "    " + str(j+1) + "\n" );
        f_integral.write( str( value ) + "    " + str(j+1) + "    " + str(j+1) + "    " + str(i+1) + "    " + str(i+1) + "\n" );

    for i in range( 0, self.norb ):
      f_integral.write( str( self.tii ) + "    " + str(i+1) + "    " + str(i+1) + "    " + str(0) + "    " + str(0) + "\n" );
    for i in range( 0, self.norb - 1, 1 ):
      f_integral.write( str( self.tij ) + "    " + str(i+1) + "    " + str(i+2) + "    " + str(0) + "    " + str(0) + "\n" );
      f_integral.write( str( self.tij ) + "    " + str(i+2) + "    " + str(i+1) + "    " + str(0) + "    " + str(0) + "\n" );

    f_integral.write( str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + "\n" );
    f_integral.close();


  def write_extended_hubbard_hydrogen_chain_int( self ):

    self.integral_path = "/" + self.prefix + "/" + self.integral_name;
    f_integral = open( self.integral_path, "wt" );
    f_integral.write( " &FCI  NORB= " + str( self.norb ) + ",NELEC= " + str( self.nelec ) + ",MS2= " + str( self.mult) + ",\n" );
    f_integral.write( "  ORBSYM=" );
    for i in range( 0, self.norb ):
      f_integral.write( "1," );
    f_integral.write( "\n" );
    f_integral.write( "  ISYM=1\n");
    f_integral.write( " &END\n" );

    # U is 0
    # V splits into ri and rb
    rb = self.hydrogen_ext.rb;
    ri = self.hydrogen_ext.ri;
    v  = self.hydrogen_ext.v;
    p  = self.hydrogen_ext.p;
    for i in range( 0, self.norb ):
      coord_i = self.hydrogen_ext.get_1d_coord(i);
      for j in range( 0, self.norb ):
        coord_j = self.hydrogen_ext.get_1d_coord(j);
        R = abs( coord_i - coord_j );
        if R != 0:
          value = float(v) * pow( float(R), -p );
          f_integral.write( str( value ) + "    " + str( i + 1 ) + "    " + str( i + 1 ) + "    " +  str( j + 1 ) + "    " + str( j + 1 ) + "\n" );

    # tii is 0
    # but we have chemical potential term mu * sum( ci^+ ci - 1/2 )
    # nu is defined by the number of particles
    for i in range( 0, self.norb ):
      mu = float( self.nelec );
    # also now we add the PH symmetry ( ni - 1/2 ) ( nj - 1/2 ), so finally these one-body term will be addted to tii
#      coord_i = self.hydrogen_ext.get_1d_coord(i);
      tii = mu;
#      for j in range( 0, self.norb ):
#        coord_j = self.hydrogen_ext.get_1d_coord(j);
#        R = abs( coord_i - coord_j );
#        if R != 0:
#          value = -0.5e0 * float(v) * pow( float(R), -p );
#          tii += value;
      f_integral.write( str( tii ) + "     " + str( i + 1 ) + "     " + str( i + 1 ) + "     " + str(0) + "     " + str(0) + "\n" );

    # tij splits into t1 and t2
    t1 = self.hydrogen_ext.t1;
    t2 = self.hydrogen_ext.t2;
    for i in range( 0, self.norb - 1, 2 ):
      f_integral.write( str( t1 ) + "    " + str( i + 1 ) + "    " + str( i + 2 ) + "    " +  str(0) + "    " + str(0) + "\n" );
      f_integral.write( str( t1 ) + "    " + str( i + 2 ) + "    " + str( i + 1 ) + "    " +  str(0) + "    " + str(0) + "\n" );
    for i in range( 1, self.norb - 1, 2 ):
      f_integral.write( str( t2 ) + "    " + str( i + 1 ) + "    " + str( i + 2 ) + "    " +  str(0) + "    " + str(0) + "\n" );
      f_integral.write( str( t2 ) + "    " + str( i + 2 ) + "    " + str( i + 1 ) + "    " +  str(0) + "    " + str(0) + "\n" );

    core_energy = 0.0e0;
#    for i in range( 0, self.norb ):
#      value = -0.5e0 * float( self.nelec );
#      coord_i = self.hydrogen_ext.get_1d_coord(i);
#      for j in range( 0, self.norb ):
#        coord_j = self.hydrogen_ext.get_1d_coord(j);
#        R = abs( coord_i - coord_j );
#        if R != 0:
#          value += 0.125e0 * float(v) * pow( float(R), -p );
#      core_energy += value;

    f_integral.write( str(core_energy) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + "\n" );
    f_integral.close();


  def write_qc_int( self ):
    command_cp_integral = "cp " + self.integral_path + " " + self.job_name + "/" + self.job_name + ".qcdmrg.FCIDUMP";
    #print command_cp_integral;
    self.integral_path = self.job_name + "/" + self.job_name + ".qcdmrg.FCIDUMP";
    from subprocess import check_call;
    check_call( command_cp_integral, shell = True );
    return 0; 


#  def write_extended_hubbard_hydrogen_chain_int( self ):
#
#    self.integral_path = self.prefix + "/" + self.integral_name;
#    f_integral = open( self.integral_path, "wt" );
#    f_integral.write( " &FCI  NORB= " + str( self.norb ) + ",NELEC= " + str( self.nelec ) + ",MS2= " + str( self.mult) + ",\n" );
#    f_integral.write( "  ORBSYM=" );
#    for i in range( 0, self.norb ):
#      f_integral.write( "1," );
#    f_integral.write( "\n" );
#    f_integral.write( "  ISYM=1\n");
#    f_integral.write( " &END\n" );
#
#    # U is 0
#    # V splits into ri and rb
#    rb = self.hydrogen_ext.rb;
#    ri = self.hydrogen_ext.ri;
#    v  = self.hydrogen_ext.v;
#    p  = self.hydrogen_ext.p;
#    for i in range( 0, self.norb ):
#      coord_i = self.hydrogen_ext.get_1d_coord(i);
##      f_integral.write( "6.0e0" + "    " + str( i + 1 ) + "    " + str( i + 1 ) + "    " +  str( i + 1 ) + "    " + str( i + 1 ) + "\n" );
#      for j in range( 0, self.norb ):
#        coord_j = self.hydrogen_ext.get_1d_coord(j);
#        R = abs( coord_i - coord_j );
#        if R != 0:
#          value = float(v) * pow( float(R), -p );
#          f_integral.write( str( value ) + "    " + str( i + 1 ) + "    " + str( i + 1 ) + "    " +  str( j + 1 ) + "    " + str( j + 1 ) + "\n" );
#          #f_integral.write( str( value ) + "    " + str( j + 1 ) + "    " + str( j + 1 ) + "    " +  str( i + 1 ) + "    " + str( i + 1 ) + "\n" );
#
##          f_integral.write( str( 0.1e0 ) + "    " + str( i + 1 ) + "    " + str( j + 1 ) + "    " +  str( i + 1 ) + "    " + str( j + 1 ) + "\n" );
##          f_integral.write( str( 1.0e-6 ) + "    " + str( j + 1 ) + "    " + str( i + 1 ) + "    " +  str( j + 1 ) + "    " + str( i + 1 ) + "\n" );
#
#    # tii is 0
#    # but now we add the PH symmetry ( ni - 1/2 ) ( nj - 1/2 ), so finally these one-body term will be addted to tii
#    for i in range( 0, self.norb ):
#      coord_i = self.hydrogen_ext.get_1d_coord(i);
#      tii = 0.0e0;
#      for j in range( 0, self.norb ):
#        coord_j = self.hydrogen_ext.get_1d_coord(j);
#        R = abs( coord_i - coord_j );
#        if R != 0:
#          value = -0.5e0 * float(v) * pow( float(R), -p );
#          tii += value;
#      f_integral.write( str( tii ) + "     " + str( i + 1 ) + "     " + str( i + 1 ) + "     " + str(0) + "     " + str(0) + "\n" );
#
#    # tij splits into t1 and t2
#    t1 = self.hydrogen_ext.t1;
#    t2 = self.hydrogen_ext.t2;
##    for i in range( 0, self.norb ):
##      f_integral.write( str( 1.0e-6 ) + "    " + str( i + 1 ) + "    " + str( i + 1 ) + "    " +  str(0) + "    " + str(0) + "\n" );
#    for i in range( 0, self.norb - 1, 2 ):
#      f_integral.write( str( t1 ) + "    " + str( i + 1 ) + "    " + str( i + 2 ) + "    " +  str(0) + "    " + str(0) + "\n" );
#      f_integral.write( str( t1 ) + "    " + str( i + 2 ) + "    " + str( i + 1 ) + "    " +  str(0) + "    " + str(0) + "\n" );
#    for i in range( 1, self.norb - 1, 2 ):
#      f_integral.write( str( t2 ) + "    " + str( i + 1 ) + "    " + str( i + 2 ) + "    " +  str(0) + "    " + str(0) + "\n" );
#      f_integral.write( str( t2 ) + "    " + str( i + 2 ) + "    " + str( i + 1 ) + "    " +  str(0) + "    " + str(0) + "\n" );
#
#    core_energy = 0.0e0;
#    for i in range( 0, self.norb ):
#      coord_i = self.hydrogen_ext.get_1d_coord(i);
#      for j in range( 0, self.norb ):
#        coord_j = self.hydrogen_ext.get_1d_coord(j);
#        R = abs( coord_i - coord_j );
#        if R != 0:
#          value = 0.125e0 * float(v) * pow( float(R), -p );
#          core_energy += value;
#
#    f_integral.write( str(core_energy) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + "\n" );
##    f_integral.write( str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + "\n" );
#
#    f_integral.close();


  # write input file for DMRG program
  def write_input( self ):
    self.input_path = "/" + self.prefix + "/" + self.input_name;

#    print "write dmrg input to ", input_name;
    f_dmrg = open( self.input_path, "wt" );
    f_dmrg.write( "orbitals " + "/" + self.scratch + "/" + self.job_name + "/" + self.integral_name + "\n" );
    f_dmrg.write( "prefix "   + "/" + self.scratch + "/" + self.job_name + "/" + "\n" );

    f_dmrg.write( "nelec "       + str( self.nelec )   + "\n" );
    f_dmrg.write( "spin "        + str( self.mult )    + "\n" );
    f_dmrg.write( "point_group " + self.point_group    + "\n" );

    # HF occupancy is not implemented in the interfaced DMRG code
    #    f_dmrg.write( "hf_occ ");
    #    for iorb in range( 0, self.norb ):
    #      f_dmrg.write( "1 0 " );
    #    for iorb in range(nelec, norb):
    #     f_dmrg.write( "0 " );
    #    f_dmrg.write("\n");
    #    f_dmrg.write( "noreorder\n");

    f_dmrg.write( "sweep_tol " + str( self.sweep_tol ) + "\n");
    f_dmrg.write( "irrep " + str( self.irrep ) + "\n" );
    if self.stochastic_ratio > 0.01e0 :
      f_dmrg.write( "stoch_ratio " + str( self.stochastic_ratio ) + "\n" );
      f_dmrg.write( "stoch_seed " + str( self.seed ) + "\n" );
      f_dmrg.write( "stoch_sweep_end " + str( self.stochastic_sweep_end ) + "\n" );
      f_dmrg.write( "schedule\n" );
      f_dmrg.write( "0 "  + str( self.M ) + " 1.0e-7  0.0e0\n" );
      f_dmrg.write( "end\n");
    else :
      f_dmrg.write( "schedule\n" );
      f_dmrg.write( "0 "  + str( self.M / 4) + " 1.0e-5  1.0e-4\n" );
      f_dmrg.write( "4 "  + str( self.M / 2) + " 1.0e-6  1.0e-5\n" );
      f_dmrg.write( "8 "  + str( self.M / 1) + " 1.0e-7  1.0e-7\n" );
      f_dmrg.write( "12 " + str( self.M / 1) + " 5.0e-8  0.0e0\n" );
      f_dmrg.write( "end\n");

#    f_dmrg.write( "twodot_to_onedot " + str( twodot_to_onedot ) + "\n" );
    f_dmrg.write( self.dot_algorithm + "\n" );
    f_dmrg.write( "maxiter " + str(self.maxiter) + "\n" );
    f_dmrg.write( "nroots "  + str(self.nroot) + "\n" );

# We don't do excited states here so this functionally is neglected
#  f_dmrg.write( "weights " );
#  for iweight in range( 0, len( weights ) ):
#    f_dmrg.write( str( weights[iweight] ) + " " );
#  f_dmrg.write( "\n" );

    f_dmrg.write( "outputlevel 5\n" );
    if self.onepdm == True:
      f_dmrg.write( "onepdm\n" );
    if self.twopdm == True:
      f_dmrg.write( "twopdm\n" );
# f_dmrg.write( "docd\n" );
    f_dmrg.close();
    return self.input_path;


  def get_energy( self ):
    import decimal;
    retval = 0.0e0;
    dmrg_energy_file = self.prefix + "/dmrg.e";
    f_energy = open( dmrg_energy_file, "rt" );
    line = f_energy.readline();
    retval = decimal.Decimal( line );
    f_energy.close();
    return retval;


  def get_onepdm( self ):
    from dmrg import dmrg_data;
    import decimal;
    retval = dmrg_data.ONEPDM();
    retval.resize( self.norb );
    dmrg_onepdm_file = self.prefix + "/" + "spatial_onepdm.0.0";
    f_onepdm = open( dmrg_onepdm_file, "rt" );
    line = f_onepdm.readline();
    for i in range( self.norb ):
      for j in range( self.norb ):
        line = f_onepdm.readline();
        fields = line.split( ' ' );
        ind_i = int( fields[0] );
        ind_j = int( fields[1] );
        value = decimal.Decimal( fields[2] );
        retval.store[ind_i][ind_j] = value;
    f_onepdm.close();
    return retval;

  def compute_sigma( self, hole_density ):
    xs = [];
    for i in range( 0, self.norb ):
      coord_i = self.hydrogen_ext.get_1d_coord(i);
      xs.append( coord_i );

#    mid = ( -xs[0] + xs[ self.norb - 1 ] ) / 2.0e0;
#    for i in range( 0, self.norb ):
#      xs[i] -= mid;
#      print xs[i];
#    print xs;
    x_square_exp = 0.0e0;
    for i in range( 0, self.norb ):
      value = float(xs[i]) * float(xs[i]) * float(hole_density[i]);
      x_square_exp += value;

    x_exp_square = 0.0e0
    for i in range( 0, self.norb ):
      value = float(xs[i]) * float(hole_density[i]);
      x_exp_square += value;
    x_exp_square *= x_exp_square;

    print x_square_exp, x_exp_square;
    import math;
    sigma = math.sqrt( x_square_exp - x_exp_square );

    return sigma;


  def compute_sigma_qc( self, hole_density ):
    xs = [];
    for i in range( 0, self.norb ):
      coord_i = self.coordinate[i][3];
      xs.append( coord_i );

#    mid = ( -xs[0] + xs[ self.norb - 1 ] ) / 2.0e0;
#    for i in range( 0, self.norb ):
#      xs[i] -= mid;
#      print xs[i];
#    print xs;
    x_square_exp = 0.0e0;
    for i in range( 0, self.norb ):
      value = float(xs[i]) * float(xs[i]) * float(hole_density[i]);
      x_square_exp += value;

    x_exp_square = 0.0e0
    for i in range( 0, self.norb ):
      value = float(xs[i]) * float(hole_density[i]);
      x_exp_square += value;
    x_exp_square *= x_exp_square;

    #print x_square_exp, x_exp_square;
    import math;
    sigma = math.sqrt( x_square_exp - x_exp_square );

    return sigma;

  def compute_cdf( self, hole_density ):
    x = [];
    for i in range( 0, self.norb ):
      coord_i = self.coordinate[i][3];
      x.append( coord_i );
  
    cdf = [];
    integral_rho = 0.0e0;
    for ind in range( 0, len( hole_density ) ):
      if( ind > 0 ) :
        new_value = 0.5e0 * ( hole_density[ ind ] + hole_density[ ind - 1] ) * ( x[ind] - x[ind-1] );
        integral_rho += new_value;
      cdf.append( integral_rho );
  
    return cdf
