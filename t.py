import chromopredict as cp

# original woodward-fieser rules
abs_max, description, image = cp.predict(
  smiles='O=C(OC)c(c(N)ccc1)c1',
  solvent=None,
  verbose=True, # return increments of all structural features
  chromlib='woodward')

#refined rules by us
abs_max, description, image = cp.predict(
  smiles='O=C(OC)c(c(N)ccc1)c1',
  solvent=None,
  verbose=True,
  chromlib='woodward_refine')
