
class Hydrogen_Extended_Hubbard:

  def __init__( self ):

    self.t1 = -1.0e0;
    self.t2 = 0.1e0;
    self.rb = 1;
    self.ri = 3;
    self.v  = 2;
    self.p  = 1;


  def get_1d_coord( self, ind ):
    n_site = ind + 1;
    nb_short = 0;
    nb_long  = 0;
    if  n_site % 2  == 0 :
      nb_short = n_site/2;
      nb_long = nb_short - 1;
    else:
      nb_short = n_site/2;
      nb_long = nb_short;
    value = self.rb * nb_short + self.ri * nb_long;
    return value;

