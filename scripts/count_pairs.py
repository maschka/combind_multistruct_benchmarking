from glob import glob
from shared_paths import shared_paths
from density_estimate import DensityEstimate


total, mcss = 0, 0
for fname in glob('{}/*/stats/{}/*mcss-reference.de'.format(shared_paths['data'], shared_paths['stats']['version'])):
	de = DensityEstimate.read(fname)
	total += 1
	mcss += de.n_samples > 0
print(total, mcss)