"""
Configuration settings for Visual Phaser application.
Contains default values and user-configurable parameters.
"""

# Default paths - will be overridden by user uploads in the Streamlit app
DATA_DIR = "data"
MAP_DIR = "maps"
DEFAULT_MAP_FILE = "maps/min_map.txt"

# Default settings (can be modified in the Streamlit UI)
DEFAULT_CONFIG = {
    # Core comparison settings
    "SIBLINGS": [],
    "COUSINS": [],
    "PHASED_FILES": [],
    "CHROMOSOMES": [],
    
    # Display options
    "SHOW_NO_MATCHES": True,
    "CHROM_TRUE_SIZE": True,
    "LINEAR_CHROMOSOME": True,
    "RESOLUTION": 1,
    "AUTO_REC_PNTS": False,
    "ARP_TOLERANCE": 4,
    "HIR_CUTOFF": 3,
    "FIR_CUTOFF": 1,
    "FIR_TABLES": True,
    "SCALE_ON": True,
    
    # Processing settings
    "SHOW_TIMES": True,
    "SHOW_MATCH_PAIR_PROGRESS": True,
    
    # Technical settings
    "MM_DIST": 500,
    "HIR_SNP_MIN": 100,
    "FIR_SNP_MIN": 75,
    "NO_CALL": 'X',
}

# Font settings - Linux paths
LINUX_FONT_STRING = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"

# Add this to DEFAULT_CONFIG
DEFAULT_CONFIG["LINUX_FONT_STRING"] = LINUX_FONT_STRING
