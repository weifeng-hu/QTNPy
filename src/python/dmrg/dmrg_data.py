
class ONEPDM:

  def __init__( self ):

    self.norb = 0;
    self.store = [];


  def resize( self, norb ):

    self.norb = norb;
    for i in range( 0, norb ):
      col = [];
      for j in range( 0, norb ):
        col.append( 0.0e0 );
      self.store.append( col );


def compute_hole_density( onepdm_neutral, onepdm_polarized ):

  if onepdm_neutral.norb != onepdm_polarized.norb:
    print "Two onepdms do not share the same norb";
    quit();

  result = [];
  for i in range( 0, onepdm_neutral.norb ):
    value_i = onepdm_neutral.store[i][i] - onepdm_polarized.store[i][i];
    result.append(value_i);

  # symmetrize
  for i in range( 0, onepdm_neutral.norb ):
    sum = float(result[i])+ float(result[ onepdm_neutral.norb - i - 1 ]);
    average = sum/2.0;
    result[i] = average;
    result[ onepdm_neutral.norb - i - 1 ] = average;

  return result;


def compute_hole_density_by_bond( onepdm_neutral, onepdm_polarized ):

  if onepdm_neutral.norb != onepdm_polarized.norb:
    print "Two onepdms do not share the same norb";
    quit();

  result = [];
  for i in range( 0, onepdm_neutral.norb, 2 ):
    value_i = onepdm_neutral.store[i][i] + onepdm_neutral.store[i+1][i+1] - ( onepdm_polarized.store[i][i] + onepdm_polarized.store[i+1][i+1] );
    result.append(value_i);

  return result;
