import subprocess;

class QC_Interface():

  def __init__( self ):
    self.natom = 0;
    self.norb = 0;
    self.nelec = 0;
    self.job_name = "";
    self.coordinate = [];
    self.basis_set = "sto-3g";
    self.dir_name = "";
    self.pgm_name = "";
    self.input_name = "";

  def write_input_norma( self ):
    self.job_name = "h" + str( self.natom );
    self.input_name = self.dir_name + "/" + self.job_name + ".norma.conf";
    print "  writing input file ", self.input_name;
    fcidump_infile = self.input_name;
    f_dump = open( fcidump_infile, "wt" );
    f_dump.write( "oe_input\n" );
    f_dump.write( self.job_name + ".orca.ha\n" );
    f_dump.write( "te_input\n" );
    f_dump.write( "orca.tei\n" );
    f_dump.write( "nuclear_energy_input\n" );
    f_dump.write( "ne.out\n" );
    f_dump.write( "coeff_input\n" );
    f_dump.write( self.job_name + ".s_inversehalf.txt\n" );
    f_dump.write( "fcidump\n" );
    f_dump.write( self.job_name + ".qcdmrg.FCIDUMP\n" );
    f_dump.write( "nbas\n" );
    f_dump.write( str( self.natom ) + "\n" );
    f_dump.write( "norb\n" );
    f_dump.write( str( self.natom ) + "\n" );
    f_dump.close();

  def write_input( self ):
    self.job_name = "h" + str( self.natom );
    self.input_name = self.dir_name + "/" + self.job_name + ".inp";
    print "  writing inputfile ", self.input_name;
    # write input for orca
    fcidump_infile = self.input_name;
    f_dump = open( fcidump_infile, "wt" );
    f_dump.write('%coords\n');
    f_dump.write(' CTyp xyz\n');
    f_dump.write(' Units Angs\n');
    f_dump.write(' coords\n');
    for icoord in range( 0, self.natom ):
     f_dump.write( str(self.coordinate[icoord][0]) + " " +
                   str(self.coordinate[icoord][1]) + " " +
                   str(self.coordinate[icoord][2]) + " " +
                   str(self.coordinate[icoord][3]) ); f_dump.write('\n');
    f_dump.write(' end\n');
    f_dump.write('end\n');
    f_dump.write('\n');
    f_dump.write("!rhf " + self.basis_set + " extremescf\n");
    f_dump.write('\n');
#    f_dump.write('%output\n');
#    f_dump.write('  print[P_ReducedOrbPopMO_L] 1\n');
#    f_dump.write('  pi_threshold  40.0e0\n');
#    f_dump.write('end\n');
#    f_dump.write('\n');
    f_dump.write('$new_job\n');
    f_dump.write('%coords\n');
    f_dump.write('  CTyp xyz\n');
    f_dump.write('  Units Angs\n');
    f_dump.write('  coords\n');
    for icoord in range( 0, self.natom ):
     f_dump.write( str(self.coordinate[icoord][0]) + " " +
                   str(self.coordinate[icoord][1]) + " " +
                   str(self.coordinate[icoord][2]) + " " +
                   str(self.coordinate[icoord][3]) ); f_dump.write('\n');
    f_dump.write('  end\n');
    f_dump.write('end\n');
    f_dump.write('\n');

    f_dump.write("!casscf " + self.basis_set + " extremescf moread\n");
    f_dump.write('%moinp '); f_dump.write('"'); f_dump.write( self.job_name ); f_dump.write('.copy.gbw"\n');
    f_dump.write('%casscf nel '); f_dump.write( str( self.nelec ) ); f_dump.write('\n');
    f_dump.write('        norb '); f_dump.write( str( self.norb ) ); f_dump.write('\n');
    f_dump.write('        mult 1\n');
    f_dump.write('        nroots 1\n');
    f_dump.write('        weights[0] = 1\n');
    f_dump.write('        cistep dmrg\n');
    f_dump.write('        maxiter 1\n');
    f_dump.write('        etol 1.0e-14\n');
    f_dump.write('        doct 1\n');
    f_dump.write('        dmrg_para\n');
    f_dump.write('           SweepTol 1.0e-12\n');
    f_dump.write('           new_version true\n');
    f_dump.write('           Nroots 2\n');
    f_dump.write('           weights = 0.0e0, 1.0e0\n');
    f_dump.write('           irrep 0\n');
    f_dump.write('           use_mpi true\n');
    f_dump.write('           use_host false\n');
    f_dump.write('           np 8\n');
    f_dump.write('           pathname "/home/wh288/Block_stracking/block.spin_adapted"\n ');
    f_dump.write('           nschedule 4\n');
    f_dump.write('           m = 200, 600, 1000, 1000\n');
    f_dump.write('           iteration = 0, 4, 8, 12\n');
    f_dump.write('           david_tol = 1.0e-7, 1.0e-7, 1.0e-7, 1.0e-8\n');
    f_dump.write('           noise = 1.0e-4, 1.0e-5, 1.0e-8, 0.0e0\n');
    f_dump.write('           twodot true\n');
    f_dump.write('           Twodot_to_Onedot 14\n');
    f_dump.write('           Maxiter 20\n')
    f_dump.write('           Restart true\n');
    f_dump.write('           Reset_iter false\n');
    f_dump.write('           Savetransf false\n');
    f_dump.write('           MOM false\n');
    f_dump.write('           prefix "./"' ); f_dump.write('\n');
    f_dump.write('           hostfile "hosts"\n');
    f_dump.write(' end\n');
    f_dump.write('end\n');
    f_dump.close();

  def run_trsf( self ):
    import subprocess;
    from subprocess import check_call;
    self.write_input_norma();
    command_run = "cd " + self.dir_name + "; ulimit -s unlimited;" + self.pgm_name + " " + self.job_name + ".norma.conf 1 ";
    print "  running: ", command_run;
    check_call( command_run, shell = True );

    print " copying FCIDUMP to int/";
    command_mkdir_and_cp = "mkdir int/; cp " + self.dir_name + "/*.qcdmrg.FCIDUMP int/;";
    print "  running: ", command_mkdir_and_cp;
    check_call( command_mkdir_and_cp, shell = True );

    print "  done";


  def run( self ):
    import subprocess;
    from subprocess import check_call;
#    print "  pre-cleaning ...";
#    command_clean = "rm -rf " + "./" + self.dir_name;
#    check_call( command_clean, shell = True );
    print "  create subdir ", self.dir_name;
    command_mkdir = "mkdir -p " + "./" + self.dir_name;
    check_call( command_mkdir, shell = True );

    self.write_input();

    command_run = "cd " + self.dir_name + "; ulimit -s unlimited;" + self.pgm_name + " " + self.job_name + ".inp" + " > " + self.job_name + ".out" ;
    print "  running: ", command_run;
    check_call( command_run, shell = True );
    print "  done";

