
class CONVERGENCE:

   def __init__( self ):

     self.Ms  = [];
     self.deterministic_energy = [];
     self.stochastic_energy = [];
     self.one_over_hilbertspace_size = [];

def plot( data, plot_name ):

  import matplotlib.pyplot as plt;
  import numpy as np;

  if plot_name == "":
    plot_name = "sample";
  png_name = plot_name + ".png";

  fig = plt.figure();
  plt.xlabel( "1/(Hilbert space size)", fontsize = 12 );
  plt.ylabel( "Energy", fontsize = 12 );
  plt.xlim( 0, data.one_over_hilbertspace_size[0] );
  plt.ylim( data.deterministic_energy[-1] - 0.01, data.deterministic_energy[0] + 0.01 );
  plt.scatter( data.one_over_hilbertspace_size, data.deterministic_energy, s = 30, marker = '^', color = 'c', label = "deterministic", alpha = 0.7 );
  plt.scatter( data.one_over_hilbertspace_size, data.stochastic_energy, s = 40, marker = '8', color = 'r', label = "stochastic", alpha = 0.5 );
  plt.legend( loc = 'upper right', fontsize = 12 );
  plt.savefig( png_name );

