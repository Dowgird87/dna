# Path to DNA files.
FILES_PATH = r"D:\Genetyka\raw_dna"

# Path to .xlsx file.
WORKING_DIRECTORY = r"D:\Genetyka\vp_tab"

# Path to min_map.txt file.
MAP_PATH = r"D:\Genetyka"

# SIBLINGS to be compared. Make sure that no two files share the same name.
SIBLINGS = ['TERESA','ALICJA','BOGUSLAWAN']

# COUSINS to be compared with SIBLINGS in a pre-existing file. No two files 
COUSINS = []

# Phased files to be compared to each other.
PHASED_FILES = []
# Chromosome selected. Leave empty to select all the chromosomes.
CHROMOSOMES = []

# Excel file name. Leave ".xlsx" out.
EXCEL_FILE_NAME = "vp_tab"

# Suppress no-matches. Set to True if display of no-matches is desired.
SHOW_NO_MATCHES = True

# Chromosome true size. Set to False for normalized size.
CHROM_TRUE_SIZE = True

# Linearize the chromosome.
LINEAR_CHROMOSOME = True

# Resolution. If CHROM_TRUE_SIZE is set to False, keep under 10. Set to 100 if 
# full resolution is desired. If CHROM_TRUE_SIZE is set to True, RESOLUTION has 
# the opposite behavior to False. Chromosome lengths DECREASE with increasing 
# RESOLUTION value. Keep between 3 and 30 (30 (True) is approximately the 
# equivalent of 1 (False). If LINEAR_CHROMOSOME is set to True, RESOLUTION will
# be ignored, unless it is set to 10 (10x resolution).
RESOLUTION = 1

# Set AUTO_REC_PNTS to True if calculation of RPs is desired. ARP is not 
# activated in LINEAR_CHROMOSOME mode or when COUSINS is not empty.
AUTO_REC_PNTS = False

# When AUTO_REC_PNTS is activated, Columns with pixel numbers less than this 
# value will be deleted. Set to minimum desired column width (pixels). 
# Default = 4.
ARP_TOLERANCE = 4

# HIR Minimum segment length (cM). The default is 7.
HIR_CUTOFF = 3

# FIR cutoff. FIRs less than 1cM in length are probably not significant.
FIR_CUTOFF = 1

# Display Fir tables.
FIR_TABLES = True

# Turn scale on and off. Set to False if not required.
SCALE_ON = True

# Column to freeze. Set to "A" if freezing not required.
FREEZE_COLUMN = 'A'

# Linux font string. An alternative is:
# "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf" 
LINUX_FONT_STRING = "*/fonts/truetype/family/DejaVuSerif-Bold.ttf"

# Elapsed times are shown for each step.
SHOW_TIMES = True

# Notifies the completion of each step. Set to False if you don't want to see 
# this.
SHOW_MATCH_PAIR_PROGRESS = True

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
You shouldn't have to change the parameters below.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# The column width per pixel factor.
SCALE_FACTOR = .1355

# Minimum number of HIR SNPs default 200.
HIR_SNP_MIN = 100

# Minimum number of FIR SNPs.
FIR_SNP_MIN = 75

# Number of Kbs between mismatches to end segment. Default 1000
MM_DIST = 500

# Character assigned to no calls.
NO_CALL = 'X'