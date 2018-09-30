import dmrg
from dmrg import dmrg_driver

#new_driver_obj = dmrg_driver.DMRG_Driver()
#new_driver_obj.norb = 12
#new_driver_obj.nelec = 12
#new_driver_obj.stochastic_ratio = 0.5
#new_driver_obj.maxiter = 6
#new_driver_obj.sweeptol = 1.0e-4
repeat = 1;
#new_driver_obj.print_obj()

#new_driver_obj.run()
for ir in range( 0, repeat ):
  new_driver_obj = dmrg_driver.DMRG_Driver()
  new_driver_obj.norb = 40
  new_driver_obj.nelec = 40
  new_driver_obj.mult = 0
  new_driver_obj.stochastic_ratio = 0.0
  new_driver_obj.maxiter = 20
  new_driver_obj.sweeptol = 1.0e-4
  new_driver_obj.M = 50;
  new_driver_obj.model = "H_PPP";
  new_driver_obj.pbs_script = False;
  new_driver_obj.dot_algorithm = "twodot_to_onedot 10";
  new_driver_obj.onepdm = True;
  new_driver_obj.twopdm = False;
  new_driver_obj.run();
#  new_driver_obj.make_working_dir();
#  new_driver_obj.write_2d_hubbard_int();
#  new_driver_obj.write_input();
  print new_driver_obj.get_energy()
