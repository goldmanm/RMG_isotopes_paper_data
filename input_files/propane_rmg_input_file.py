database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='propane_ooo',
    reactive=True,
    structure=SMILES("CCC"),
)
species(
    label='He',
    reactive=False,
    structure=SMILES("[He]"),
)


# Reaction systems
simpleReactor(
    temperature=(1273,'K'),
    pressure=(2.0,'bar'),
    initialMoleFractions={
        "propane_ooo": 1,
	"He": 99,
    },
    terminationTime=(1e-1,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000,
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=False,
    saveSimulationProfiles=True,
)

generatedSpeciesConstraints(
    maximumCarbonAtoms=6,
    maximumRadicalElectrons=1,
)

