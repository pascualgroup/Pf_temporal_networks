RANDOM_SEED = 1

N_GENES_INITIAL = 12000 # pool
N_GENES_PER_STRAIN = 60
N_LOCI = 2
N_ALLELES_INITIAL = N_LOCI * [1200]

SELECTION_MODE = 'SPECIFIC_IMMUNITY'
N_INFECTIONS_FOR_GENERAL_IMMUNITY = 0
GENERAL_IMMUNITY_PARAMS = [-1, -1, -1, -1]
CLEARANCE_RATE_IMMUNE = 0.041628833086078301 # For generalized immunity

TRANSITION_RATE_NOT_IMMUNE = 1.0 / 6.0 # switching rate. every 6 days
TRANSITION_RATE_IMMUNE = 1000

PRINT_FUNCTION_TRACE = False
PRINT_DEBUG_LEVEL = 0

T_YEAR = 360.0
T_BURNIN = 0 # in days
T_END = 50.0 * T_YEAR

SAMPLE_DB_FILENAME = '"test_02.sqlite"'

PRINT_INFO_PERIOD = 120.0

VERIFICATION_ON = True
VERIFICATION_PERIOD = 50.0 * T_YEAR

SAVE_TO_CHECKPOINT = False
CHECKPOINT_SAVE_FILENAME = '""'
CHECKPOINT_SAVE_PERIOD = T_YEAR * 20

LOAD_FROM_CHECKPOINT = False
CHECKPOINT_LOAD_FILENAME = '""'

OUTPUT_HOSTS = True
OUTPUT_GENES = False # Set to false to save space
OUTPUT_STRAINS = False # Set to false to save space

HOST_SAMPLING_ON = True
HOST_SAMPLING_PERIOD = 30.0
HOST_SAMPLE_SIZE = 100

GENE_TRANSMISSIBILITY = 0.5
COINFECTION_REDUCES_TRANSMISSION = True

ECTOPIC_RECOMBINATION_RATE = 1.8e-07
P_ECTOPIC_RECOMBINATION_IS_CONVERSION = 0

IMMUNITY_LOSS_RATE = 0.00001

MUTATION_RATE = 1.42e-08

T_LIVER_STAGE = 14.0


MEAN_HOST_LIFETIME = 30.0 * T_YEAR
MAX_HOST_LIFETIME = 80.0 * T_YEAR

N_POPULATIONS = 4
DISTANCE_MAT=[[1,100,100,100],[100,1,100,100],[100,100,1,100],[100,100,100,1]]
DIST_POWER=1

N_HOSTS = [10000]
N_INITIAL_INFECTIONS = [20]

BITING_RATE_MEAN = [1]
BITING_RATE_RELATIVE_AMPLITUDE = [0.0]
BITING_RATE_PEAK_PHASE = [0.0]
DAILY_BITING_RATE_DISTRIBUTION = [1] # the "raw" values of mosquito numbers from Mathematica. needs to be of length 360. Otherwise will run without it

IRS_ON = False
IRS_START_TIMES = [] # A vector with start times for the IRS interventions. Each element is the starting time of an IRS
BITING_RATE_FACTORS = [] # An array in which each element is a vector of values of mosquito numbers as the ones in DAILY_BITING_RATE_DISTRIBUTION. These are obtained from the model externally run in Mathematica
IRS_IMMIGRATION_RATE_FACTORS = [] # A vector with length as number of IRS events. Each element is a factor to multiply the IMMIGRATION_RATE. For example, 0.3 will reduce the usual immigration rate to 30% of its original size.

MDA_ON = False
MDA_START_TIMES = []
HOST_FAIL_RATE = [] # % of hosts that did not take the drug
DRUG_EFF_DURATION = [] # How long the drug remains effective in the body
MDA_IMMIGRATION_RATE_FACTORS = []


IMMIGRATION_ON = True
IMMIGRATION_RATE = [1.0]
P_IMMIGRATION_INCLUDES_NEW_GENES = 0.5
N_IMMIGRATION_NEW_GENES = 0
POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION=False # If this is False, we assume a regional intervention, whereby the gene pool size is reduced by the same proportion as the local population gene pool size after the intervention is lifted. For example, if IRS or MDA reduced the local gene diversity to 30% of its original size before the intervention (where original size it taken at the beginning point of the intervention), then when the intervention stops, the general gene pool is also reduced to 30% of its size (e.g., for N_GENES_INITIAL = 12000 the pool size will be 3600). A value of True means no change and the regional pool size remains as if no intervention were done (i.e., 12000 in our example).
