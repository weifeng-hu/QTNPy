
class DriverObj:

  def __init__( self ):
    self.molecule_name = "";
    self.coordinate = [];
    self.natom = 0;
    self.basis_set = "sto-3g";
    self.dmrg_pbs = False;
    self.dmrg_M = 100;
    self.dmrg_norb = 0;
    self.dmrg_nelec = 0;
    self.dmrg_mult = 0;
    self.dmrg_point_group = "c1";
    self.dmrg_sweep_tol = 1.0e-6;
    self.dmrg_irrep = 1;
    self.dmrg_nroot = 1;
    self.dmrg_twodot_to_onedot = 100;
    self.dmrg_maxiter = 20;
    self.dmrg_prefix = "./";
    self.dmrg_dot_algorithm = "twodot_to_onedot 14";
    self.dmrg_onepdm = False;
    self.dmrg_twopdm = False;
    self.dmrg_np = 1;
    self.dmrg_nnode = 1;
    self.dmrg_ppn = 1;
    self.dmrg_home = "/home/"
    self.dmrg_scratch = "/scratch/"

  def run_trsf( self ):
    print "generating integral transformed ... ";
    from qc import qc_interface;
    new_interface_obj = qc_interface.QC_Interface();
    new_interface_obj.pgm_name = "/home/weifeng/local/6000.integrated_quantum_dynamics/src/integral/trsf_norma/norma";
    new_interface_obj.natom = self.natom;
    new_interface_obj.norb = self.dmrg_norb;
    new_interface_obj.nelec = self.dmrg_nelec;
    new_interface_obj.coordinate = self.coordinate;
    new_interface_obj.basis_set = self.basis_set;
    new_interface_obj.dir_name = "ov";
    new_interface_obj.run_trsf();


  def run_int( self ):
    print "generating integral ... ";
    from qc import qc_interface;
    new_interface_obj = qc_interface.QC_Interface();
    new_interface_obj.pgm_name = "/home/weifeng/local/7013.quantum_chemistry.orca_test_pkg/x86_exe/orca";
    new_interface_obj.natom = self.natom;
    new_interface_obj.norb = self.dmrg_norb;
    new_interface_obj.nelec = self.dmrg_nelec;
    new_interface_obj.coordinate = self.coordinate;
    new_interface_obj.basis_set = self.basis_set;
    new_interface_obj.dir_name = "int";
    new_interface_obj.run();


  def run_ov( self ):
    print "generating C and S ... ";
    from qc import qc_interface;
    new_interface_obj = qc_interface.QC_Interface();
    new_interface_obj.pgm_name = "/home/weifeng/local/7014.quantum_chemistry.orca_lno_pkg/x86_exe/orca";
    new_interface_obj.natom = self.natom;
    new_interface_obj.norb = self.dmrg_norb;
    new_interface_obj.nelec = self.dmrg_nelec;
    new_interface_obj.coordinate = self.coordinate;
    new_interface_obj.basis_set = self.basis_set;
    new_interface_obj.dir_name = "ov";
    new_interface_obj.run();

  def run_s( self ):

    print "generating S ... ";
    from qc import qc_interface;
    new_interface_obj = qc_interface.QC_Interface();
    new_interface_obj.pgm_name = "/home/weifeng/local/7014.quantum_chemistry.orca_lno_pkg/x86_exe/orca";
    new_interface_obj.natom = self.natom;
    new_interface_obj.norb = self.dmrg_norb;
    new_interface_obj.nelec = self.dmrg_nelec;
    new_interface_obj.coordinate = self.coordinate;
    new_interface_obj.basis_set = self.basis_set;
    new_interface_obj.dir_name = "ov";
    new_interface_obj.run();

    print "computing S^(-1/2) ... ";
    from utility import ph_analysis;
    from utility.ph_analysis import python;
    from utility.ph_analysis.python import interface;
    new_ph_analysis_obj           = interface.Interface();
    new_ph_analysis_obj.norb      = self.dmrg_norb;
    new_ph_analysis_obj.job_name  = "h" + str( self.natom );
    new_ph_analysis_obj.dir_name  = "ov";
    new_ph_analysis_obj.pgm_name  = "/home/weifeng/local/6000.integrated_quantum_dynamics/src/utility/ph_analysis/ph_trans";
    new_ph_analysis_obj.run_s_only();


  def run_dmrg( self ):
    from dmrg import dmrg_driver;
    new_dmrg_driver_obj_neutral = dmrg_driver.DMRG_Driver();
    new_dmrg_driver_obj_neutral.norb            = self.dmrg_norb;
    new_dmrg_driver_obj_neutral.nelec           = self.dmrg_nelec;
    new_dmrg_driver_obj_neutral.mult            = 0;
    new_dmrg_driver_obj_neutral.maxiter         = self.dmrg_maxiter;
    new_dmrg_driver_obj_neutral.sweeptol        = self.dmrg_sweep_tol;
    new_dmrg_driver_obj_neutral.M               = self.dmrg_M;
    new_dmrg_driver_obj_neutral.model           = "QC";
    new_dmrg_driver_obj_neutral.stochastic_ratio= 0.0;
    new_dmrg_driver_obj_neutral.pbs_script      = self.dmrg_pbs;
    new_dmrg_driver_obj_neutral.dot_algorithm   = self.dmrg_dot_algorithm;
    new_dmrg_driver_obj_neutral.onepdm          = self.dmrg_onepdm;
    new_dmrg_driver_obj_neutral.twopdm          = self.dmrg_twopdm;
    new_dmrg_driver_obj_neutral.home            = self.dmrg_home;
    new_dmrg_driver_obj_neutral.scratch            = self.dmrg_scratch;
    new_dmrg_driver_obj_neutral.np              = self.dmrg_np;
    new_dmrg_driver_obj_neutral.ppn             = self.dmrg_ppn;
    new_dmrg_driver_obj_neutral.nnode           = self.dmrg_np/self.dmrg_ppn;
    new_dmrg_driver_obj_neutral.integral_path   = "int/h" + str( self.natom ) + ".qcdmrg.FCIDUMP";
    print "runing dmrg for neutral ... ", new_dmrg_driver_obj_neutral.job_name;
    new_dmrg_driver_obj_neutral.run();

    from time import sleep;
    sleep(2);

    new_dmrg_driver_obj_polarized = dmrg_driver.DMRG_Driver();
    new_dmrg_driver_obj_polarized.norb          = self.dmrg_norb;
    new_dmrg_driver_obj_polarized.nelec         = self.dmrg_nelec - 1;
    new_dmrg_driver_obj_polarized.mult          = 1;
    new_dmrg_driver_obj_polarized.maxiter       = self.dmrg_maxiter;
    new_dmrg_driver_obj_polarized.sweeptol      = self.dmrg_sweep_tol;
    new_dmrg_driver_obj_polarized.M             = self.dmrg_M;
    new_dmrg_driver_obj_polarized.model         = "QC";
    new_dmrg_driver_obj_polarized.stochastic_ratio= 0.0;
    new_dmrg_driver_obj_polarized.pbs_script    = self.dmrg_pbs;
    new_dmrg_driver_obj_polarized.dot_algorithm = self.dmrg_dot_algorithm;
    new_dmrg_driver_obj_polarized.onepdm        = self.dmrg_onepdm;
    new_dmrg_driver_obj_polarized.twopdm        = self.dmrg_twopdm;
    new_dmrg_driver_obj_polarized.home          = self.dmrg_home;
    new_dmrg_driver_obj_polarized.scratch          = self.dmrg_scratch;
    new_dmrg_driver_obj_polarized.np            = self.dmrg_np;
    new_dmrg_driver_obj_polarized.ppn           = self.dmrg_ppn;
    new_dmrg_driver_obj_polarized.nnode         = self.dmrg_np/self.dmrg_ppn;
    new_dmrg_driver_obj_polarized.integral_path   = "int/h" + str( self.natom ) + ".qcdmrg.FCIDUMP";
    print "runing dmrg for polarized ... ", new_dmrg_driver_obj_polarized.job_name;
    new_dmrg_driver_obj_polarized.run();

    return [ new_dmrg_driver_obj_neutral.job_name, new_dmrg_driver_obj_polarized.job_name ];


  def run_transform( self, job_pair ):
    from utility import ph_analysis;
    from utility.ph_analysis import python;
    from utility.ph_analysis.python import interface;
    
#    new_umat_obj = interface.Interface();
#    new_umat_obj.norb = self.dmrg_norb;
#    new_umat_obj.job_name = "h" + str(self.natom);
#    new_umat_obj.dir_name = "umat";
#    new_umat_obj.pgm_name = "/local_scratch/weifeng/local/6000.integrated_quantum_dynamics/src/utility/ph_analysis/ph_trans";
#    new_umat_obj.run_ctsh();

    new_pdm_transform_neutral_obj           = interface.Interface();
    new_pdm_transform_neutral_obj.norb      = self.dmrg_norb;
    new_pdm_transform_neutral_obj.job_name  = "h" + str( self.natom );
    new_pdm_transform_neutral_obj.mo_pdm_dir_name  = job_pair[0];
    new_pdm_transform_neutral_obj.dir_name  = job_pair[0] + ".pdm.ao";
    new_pdm_transform_neutral_obj.pgm_name  = "/home/weifeng/local/6000.integrated_quantum_dynamics/src/utility/ph_analysis/ph_trans";
    new_pdm_transform_neutral_obj.run_transform();

    new_pdm_transform_polarized_obj           = interface.Interface();
    new_pdm_transform_polarized_obj.norb      = self.dmrg_norb;
    new_pdm_transform_polarized_obj.job_name  = "h" + str( self.natom );
    new_pdm_transform_polarized_obj.mo_pdm_dir_name  = job_pair[1];
    new_pdm_transform_polarized_obj.dir_name  = job_pair[1] + ".pdm.ao";
    new_pdm_transform_polarized_obj.pgm_name  = "/home/weifeng/local/6000.integrated_quantum_dynamics/src/utility/ph_analysis/ph_trans";
    new_pdm_transform_polarized_obj.run_transform();

    return [ new_pdm_transform_neutral_obj.dir_name, new_pdm_transform_polarized_obj.dir_name ];


  def run_analysis( self, job_pair ):

    from dmrg import dmrg_driver;
    new_dmrg_analysis_obj_neutral          = dmrg_driver.DMRG_Driver();
    new_dmrg_analysis_obj_neutral.coordinate = self.coordinate;
    new_dmrg_analysis_obj_neutral.norb     = self.dmrg_norb;
    new_dmrg_analysis_obj_neutral.job_name = job_pair[0];
    new_dmrg_analysis_obj_neutral.prefix   = job_pair[0];
    onepdm_neutral = new_dmrg_analysis_obj_neutral.get_onepdm();
    energy_neutral = new_dmrg_analysis_obj_neutral.get_energy();

    new_dmrg_analysis_obj_polarized          = dmrg_driver.DMRG_Driver();
    new_dmrg_analysis_obj_polarized.coordinate = self.coordinate;
    new_dmrg_analysis_obj_polarized.norb     = self.dmrg_norb;
    new_dmrg_analysis_obj_polarized.job_name = job_pair[1];
    new_dmrg_analysis_obj_polarized.prefix   = job_pair[1];
    onepdm_polarized = new_dmrg_analysis_obj_polarized.get_onepdm();
    energy_polarized = new_dmrg_analysis_obj_polarized.get_energy();

    from dmrg import dmrg_data;
    hole_density = dmrg_data.compute_hole_density( onepdm_neutral, onepdm_polarized );
    for i in range( 0, len( hole_density ) ):
      print i, hole_density[i];

    print "CDF:";

    cdf = new_dmrg_analysis_obj_neutral.compute_cdf( hole_density );
    for i in range( 0, len( hole_density ) ):
      print i, cdf[i];

    print "sigma: ", new_dmrg_analysis_obj_polarized.compute_sigma_qc( hole_density ), "  energy gap: ", str( energy_polarized - energy_neutral );

  def run( self ):
    self.read_coordinates();
    self.clean();
    self.run_int();
    self.run_ov();
    job_pair = self.run_dmrg();
    if( self.dmrg_pbs == False ):
      pdm_pair = self.run_transform( job_pair );
      self.run_analysis( pdm_pair );

  def run_ao( self ):
    self.read_coordinates();
    self.clean();
    self.run_s();
    self.run_trsf();
    job_pair = self.run_dmrg();
    if( self.dmrg_pbs == False ):
      self.run_analysis( job_pair );

  def read_coordinates( self ):
    xyz_file = self.molecule_name + ".xyz";
    fxyz = open( xyz_file, "rt" );
    line = fxyz.readline();
    self.natom = int( line );
    for iline in range( 0, self.natom ):
      line = fxyz.readline();
      fields = line.split();
      atom = fields[0];
      x = float( fields[1] );
      y = float( fields[2] );
      z = float( fields[3] );
      new_coord = [ atom, x, y, z ];
      self.coordinate.append( new_coord );
    fxyz.close();

  def clean( self ):
    print "pre-cleaning ... "
    from subprocess import check_call;
    command_clean = "rm -rf int/ ov/ umat/ dmrg* *sh";
    check_call( command_clean, shell = True ); 
