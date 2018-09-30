import dmrg_driver
import sys

nsites = [ 8 ];
M      = [ 100 ];
nroot  = 1;
e      = -2.0e0;
t      = -0.5e0;
u      = 6.0e0;
repeat = 10;

for ir in range( 0, repeat ):
  for i_size in range( 0, len(nsites) ):
    for iM in range( 0, len(M) ):
      dmrg_driver_obj = dmrg_driver.DMRG_Driver();
      dmrg_driver_obj.M = M[iM];
      dmrg_driver_obj.tii = e;
      dmrg_driver_obj.tij = t;
      dmrg_driver_obj.u   = u;
      dmrg_driver_obj.norb = nsites[ i_size ];
      dmrg_driver_obj.job_name  = "dmrg." + ".n." + str( nsites[isite] ) + ".M." + str( M[iM] ) + ".e." + str( e ) + ".t." + str( t ) + ".u." + str( u ) + ".repeat." + str( ir );
      dmrg_driver_obj.run();
      energy = dmrg_driver_obj.get_energy();
      print "n:\t", str( nsites[isite] ), "\t\tM:\t", str( M[iM] ), "\t\te:\t", str( e ), "\t\tt:\t", str( t ), "\t\tu:\t", str( u ), "\t\tEn:\t", energy
