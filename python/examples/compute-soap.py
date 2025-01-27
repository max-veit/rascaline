import sys
import chemfiles

from rascaline import SoapPowerSpectrum

# read structures using chemfiles
with chemfiles.Trajectory(sys.argv[1]) as trajectory:
    frames = [f for f in trajectory]
# or using ASE, at your own convenience
# frames = ase.io.read(sys.argv[1], ":")

# define hyper parameters for the calculation
HYPER_PARAMETERS = {
    "cutoff": 5.0,
    "max_radial": 6,
    "max_angular": 4,
    "atomic_gaussian_width": 0.3,
    "gradients": False,
    "radial_basis": {
        "Gto": {},
    },
    "cutoff_function": {
        "ShiftedCosine": {"width": 0.5},
    },
}

calculator = SoapPowerSpectrum(**HYPER_PARAMETERS)

# run the actual calculation
descriptor = calculator.compute(frames)

# Transform the descriptor to dense representation,
# with one sample for each atom-centered environment
descriptor.densify(["species_neighbor_1", "species_neighbor_2"])

# you can now use descriptor.values as the
# input of a machine learning algorithm
print(descriptor.values.shape)
