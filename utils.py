"""
Utility functions for the Visual Phaser application.
Includes data loading, processing, and image generation.
"""

import os
import pandas as pd
import numpy as np
import time
from PIL import Image, ImageDraw, ImageFont
import platform
from io import BytesIO
import streamlit as st
import base64

def get_dna_files(data_dir):
    """List all DNA files in the data directory with full paths."""
    if os.path.exists(data_dir):
        # Get files with .txt or .csv extensions, filter out any hidden files (starting with .)
        files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) 
                if f.endswith(('.txt', '.csv')) and not f.startswith('.')]
        return files
    return []

def extract_individual_names(files):
    """Extract individual names from DNA filenames."""
    names = []
    for filepath in files:
        # Get just the filename without the directory path
        filename = os.path.basename(filepath)
        
        # Check for different common DNA file naming patterns
        if "_raw" in filename:
            name = filename.split("_raw")[0]
        elif "-raw" in filename:
            name = filename.split("-raw")[0]
        elif "_data" in filename:
            name = filename.split("_data")[0]
        elif "-data" in filename:
            name = filename.split("-data")[0]
        # Extract name from .txt or .csv extension
        elif filename.endswith(".txt") or filename.endswith(".csv"):
            name = os.path.splitext(filename)[0]
        else:
            # If no pattern matches, use the whole filename
            name = filename
            
        if name and name not in names:
            names.append(name)
            
    return names

def cm_calc(start_pos, end_pos, map_data):
    """
    Calculate genetic distance in centiMorgans between two positions
    using the genetic map data.
    """
    if map_data is None or map_data.empty:
        # Fallback without map: rough approximation (1Mb ≈ 1cM)
        return (end_pos - start_pos) / 1_000_000
    
    # Check column names - map file might use 'position' or 'pos' for positions
    # and 'cm', 'cM', or 'centimorgans' for genetic distance
    pos_col = None
    cm_col = None
    
    # Try to find the position column
    for possible_pos_col in ['position', 'pos', 'Position', 'POSITION']:
        if possible_pos_col in map_data.columns:
            pos_col = possible_pos_col
            break
    
    # Try to find the centimorgan column
    for possible_cm_col in ['cm', 'cM', 'centimorgans', 'centiMorgans', 'CM']:
        if possible_cm_col in map_data.columns:
            cm_col = possible_cm_col
            break
    
    # If we couldn't find appropriate columns, use the fallback
    if pos_col is None or cm_col is None:
        return (end_pos - start_pos) / 1_000_000
    
    # Find closest positions in the map
    start_row = map_data[map_data[pos_col] >= start_pos].iloc[0] if not map_data[map_data[pos_col] >= start_pos].empty else None
    end_row = map_data[map_data[pos_col] >= end_pos].iloc[0] if not map_data[map_data[pos_col] >= end_pos].empty else None
    
    if start_row is None or end_row is None:
        # Fallback if positions are outside map range
        return (end_pos - start_pos) / 1_000_000
    
    # Return cM difference
    return end_row[cm_col] - start_row[cm_col]

def load_map_data(map_path):
    """Load genetic map data from file."""
    try:
        if os.path.exists(map_path):
            # Try to load with different separators and options
            try:
                map_data = pd.read_csv(map_path, sep="\t", header=0)
            except:
                try:
                    map_data = pd.read_csv(map_path, sep=",", header=0)
                except:
                    # Last attempt, try with different separator detection
                    map_data = pd.read_csv(map_path, sep=None, engine='python', header=0)
            
            # Print loaded columns for debugging
            st.write(f"Map data columns: {map_data.columns.tolist()}")
            
            # If we have a chromosome column, filter to make sure it's numeric
            if 'chromosome' in map_data.columns:
                # Replace 'X' with '23' if present
                if map_data['chromosome'].dtype == 'object':
                    map_data.replace({'chromosome': {'X': '23', 'XY': '23'}}, inplace=True)
                    # Convert to numeric if possible
                    map_data['chromosome'] = pd.to_numeric(map_data['chromosome'], errors='coerce')
                    # Drop any rows with non-numeric chromosomes
                    map_data = map_data.dropna(subset=['chromosome'])
            
            return map_data
        return None
    except Exception as e:
        st.error(f"Error loading map data: {e}")
        return None

def save_uploaded_file(file, directory, filename=None):
    """Save an uploaded file to the specified directory."""
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    if filename is None:
        filename = file.name
    
    file_path = os.path.join(directory, filename)
    with open(file_path, "wb") as f:
        f.write(file.getbuffer())
    
    return file_path

def image_to_base64(img):
    """Convert a PIL image to base64 for embedding in HTML."""
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    return base64.b64encode(buffered.getvalue()).decode()

def get_chromosome_image(dplot, pair_name, chrom, config):
    """Generate a PIL Image representing the chromosome matches."""
    # Handle empty data case
    if dplot is None or len(dplot) == 0:
        # Return a placeholder image with error message when no data
        img = Image.new("RGB", (400, 50), color="white")
        img1 = ImageDraw.Draw(img)
        img1.text((10, 20), "No data to display for this chromosome", fill="black")
        return img
        
    # Ensure minimum image width
    width = max(len(dplot), 10)  # Minimum width of 10 pixels
    
    # Create image with 50px height as requested
    img = Image.new("RGB", (width, 50), color="white")
    img1 = ImageDraw.Draw(img)

    pos = 0
    for i in range(len(dplot)):
        try:
            # Get the color values from the dataframe
            match_color = dplot.loc[i, 'match']
            bar_color = dplot.loc[i, 'bar']
            
            # Draw lines for the top and bottom halves
            # Adjust line heights to fill the 50px height exactly
            img1.line([(pos, 0), (pos, 24)], fill=match_color, width=0)
            img1.line([(pos, 25), (pos, 49)], fill=bar_color, width=0)
            pos += 1
        except Exception as e:
            # Skip errors in individual pixels
            pass

    return img

def get_scale_image(dplot, chrom, config):
    """Generate a PIL Image with a scale bar for the chromosome."""
    # Handle empty data case
    if dplot is None or len(dplot) == 0:
        # Return a placeholder image when no data
        img = Image.new("RGB", (400, 50), color="white")
        img1 = ImageDraw.Draw(img)
        
        # Load font
        try:
            if platform.system() == 'Windows':
                fnt = ImageFont.truetype("arial.ttf", 10)
            elif platform.system() == 'Darwin':
                fnt = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 11)
            else:  # Linux
                fnt = ImageFont.truetype(config["LINUX_FONT_STRING"], 11)
                
            img1.text((10, 10), "No scale data available", fill="black", font=fnt)
        except:
            # If font loading fails, still draw text without specified font
            img1.text((10, 10), "No scale data available", fill="black")
            
        return img
    
    # Ensure minimum image width
    width = max(len(dplot), 100) + 30  # Minimum width plus margin
    
    img = Image.new("RGB", (width, 50), color="white")
    img1 = ImageDraw.Draw(img)

    # Font handling - wrap in try/except to handle font issues
    try:
        if platform.system() == 'Windows':
            fnt = ImageFont.truetype("arial.ttf", 10)
            fnt1 = ImageFont.truetype("arial.ttf", 8)
        elif platform.system() == 'Darwin':
            fnt = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 11)
            fnt1 = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 8)
        else:  # Linux
            fnt = ImageFont.truetype(config["LINUX_FONT_STRING"], 11)
            fnt1 = ImageFont.truetype(config["LINUX_FONT_STRING"], 8)
    except Exception as e:
        # Default to no font if there's an issue with font loading
        fnt = None
        fnt1 = None

    try:
        # Find position column if it exists
        position_col = 'position'
        if position_col not in dplot.columns and 'pos' in dplot.columns:
            position_col = 'pos'
            
        # Check if the position column exists
        if position_col not in dplot.columns:
            img1.text((10, 10), "No position data in dataframe", fill="black", font=fnt if fnt else None)
            return img
            
        max_position = int(dplot[position_col].max()) if position_col in dplot.columns else 0
        if max_position == 0:
            img1.text((10, 10), "Zero maximum position detected", fill="black", font=fnt if fnt else None)
            return img
            
        pixels_per_base = len(dplot) / max_position if max_position != 0 else 0

        # ticks every 0.5 Mb
        for position in range(0, max_position + 1, int(0.5 * 10**6)):
            pixel_position = int(position * pixels_per_base)
            if position % (2 * 10**6) == 0:
                # major tick
                if fnt:
                    img1.text((pixel_position - 7, 0), f"{position/1e6:.0f}", font=fnt, fill="black")
                img1.line([(pixel_position, 20), (pixel_position, 48)], fill="black", width=1)
            elif position % (10**6) == 0:
                # medium tick
                img1.line([(pixel_position, 30), (pixel_position, 48)], fill="black", width=1)
            else:
                # short tick
                img1.line([(pixel_position, 35), (pixel_position, 48)], fill="black", width=1)
    except Exception as e:
        # If any error occurs, return a basic image with error message
        img = Image.new("RGB", (400, 50), color="white")
        img1 = ImageDraw.Draw(img)
        img1.text((10, 10), f"Error generating scale: {str(e)}", fill="black", font=fnt if fnt else None)

    return img

def conditions(al1x, al2x, al1y, al2y, no_call):
    """Return color code depending on matches / mismatches / no-calls."""
    if al1x == no_call or al1y == no_call:
        return 'grey'  # Grey for no-call SNPs
    elif al1x == al2x and al1y == al2y and al1x != al1y:
        return 'crimson'  # Crimson for Half Identical Regions (HIR) - homozygous but different alleles
    elif (al1x == al1y and al2x == al2y) or (al1x == al2y and al2x == al1y):
        return 'limegreen'  # Limegreen for Fully Identical Regions (FIR) - shared alleles
    else:
        return 'yellow'  # Yellow for "half match" scenarios
        
def save_analysis_results(results, chromosome, pair, config, timestamp=None):
    """
    Save analysis results to the results directory.
    
    Args:
        results (dict): Analysis results containing HIR/FIR data and images
        chromosome (int): Chromosome number
        pair (tuple): Pair of individuals (name1, name2)
        config (dict): Configuration parameters used for the analysis
        timestamp (str, optional): Timestamp to use. If None, generate a new one.
        
    Returns:
        str: Analysis ID (timestamp) used to save the results
    """
    import json
    import pickle
    from datetime import datetime
    
    # Create results directory if it doesn't exist
    results_dir = "results"
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
        
    # Create a timestamp if not provided
    if timestamp is None:
        # Use a format that's file-safe for directories but will be displayed with full date and time
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
    # Create a unique analysis ID
    analysis_id = timestamp
    pair_str = f"{pair[0]}-{pair[1]}"
    
    # Create directory for this analysis
    analysis_dir = os.path.join(results_dir, f"{pair_str}_chr{chromosome}_{analysis_id}")
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
        
    # Save configuration as JSON
    with open(os.path.join(analysis_dir, "config.json"), 'w') as f:
        # Convert any NumPy or non-serializable types to standard Python types
        config_copy = {}
        for key, value in config.items():
            # NumPy 2.0 compatibility - check for numeric types more generically
            if hasattr(value, 'item') and callable(getattr(value, 'item')):
                # This handles both np.integer and np.floating types
                try:
                    config_copy[key] = value.item()
                except:
                    config_copy[key] = float(value) if '.' in str(value) else int(value)
            elif isinstance(value, (list, dict, str, int, float, bool, type(None))):
                config_copy[key] = value
            else:
                # Convert non-serializable types to string
                config_copy[key] = str(value)
        
        json.dump(config_copy, f, indent=2)
        
    # Save metadata (pair, chromosome, timestamp)
    metadata = {
        "pair": pair,
        "chromosome": chromosome,
        "timestamp": timestamp,
        "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }
    with open(os.path.join(analysis_dir, "metadata.json"), 'w') as f:
        json.dump(metadata, f, indent=2)
        
    # Save table data as CSV if available
    if results.get("hir_data") is not None:
        results["hir_data"].to_csv(os.path.join(analysis_dir, "hir_data.csv"), index=False)
    if results.get("fir_data") is not None:
        results["fir_data"].to_csv(os.path.join(analysis_dir, "fir_data.csv"), index=False)
        
    # Save images (using PNG format)
    if "chr_image" in results and results["chr_image"] is not None:
        results["chr_image"].save(os.path.join(analysis_dir, "chromosome.png"))
    if "scale_image" in results and results["scale_image"] is not None:
        results["scale_image"].save(os.path.join(analysis_dir, "scale.png"))
        
    # Save complete results using pickle (for internal use)
    with open(os.path.join(analysis_dir, "full_results.pkl"), 'wb') as f:
        # Remove Pillow images from pickled data (we saved them separately)
        pickle_results = {k: v for k, v in results.items() 
                         if k not in ["chr_image", "scale_image"]}
        pickle.dump(pickle_results, f)
        
    return analysis_id

def load_analysis_results(analysis_path):
    """
    Load saved analysis results from the specified path.
    
    Args:
        analysis_path (str): Path to the analysis directory
        
    Returns:
        tuple: (metadata, config, results) - All data for the analysis
    """
    import json
    import pickle
    from PIL import Image
    
    metadata = None
    config = None
    results = None
    
    try:
        # Load metadata
        with open(os.path.join(analysis_path, "metadata.json"), 'r') as f:
            metadata = json.load(f)
            
        # Load configuration
        with open(os.path.join(analysis_path, "config.json"), 'r') as f:
            config = json.load(f)
            
        # Load full results
        with open(os.path.join(analysis_path, "full_results.pkl"), 'rb') as f:
            results = pickle.load(f)
            
        # Load images
        if os.path.exists(os.path.join(analysis_path, "chromosome.png")):
            results["chr_image"] = Image.open(os.path.join(analysis_path, "chromosome.png"))
        
        if os.path.exists(os.path.join(analysis_path, "scale.png")):
            results["scale_image"] = Image.open(os.path.join(analysis_path, "scale.png"))
            
        # Load CSV data if needed (normally loaded from pickle, this is a fallback)
        if "hir_data" not in results and os.path.exists(os.path.join(analysis_path, "hir_data.csv")):
            results["hir_data"] = pd.read_csv(os.path.join(analysis_path, "hir_data.csv"))
        if "fir_data" not in results and os.path.exists(os.path.join(analysis_path, "fir_data.csv")):
            results["fir_data"] = pd.read_csv(os.path.join(analysis_path, "fir_data.csv"))
            
    except Exception as e:
        st.error(f"Error loading analysis results: {str(e)}")
        return None, None, None
        
    return metadata, config, results

def get_saved_analyses():
    """
    Get list of all saved analyses from the results directory.
    
    Returns:
        dict: Dictionary with keys as chromosome+pair combinations and values as lists of analysis metadata
    """
    import json
    import glob
    
    results_dir = "results"
    if not os.path.exists(results_dir):
        return {}
        
    analyses = {}
    
    # Find all analysis directories
    analysis_dirs = glob.glob(os.path.join(results_dir, "*"))
    
    for analysis_dir in analysis_dirs:
        try:
            # Load metadata
            metadata_path = os.path.join(analysis_dir, "metadata.json")
            if os.path.exists(metadata_path):
                with open(metadata_path, 'r') as f:
                    metadata = json.load(f)
                    
                # Create key for chromosome-pair combination
                pair = metadata.get("pair", ["Unknown", "Unknown"])
                chromosome = metadata.get("chromosome", "Unknown")
                key = f"chr{chromosome}_{pair[0]}-{pair[1]}"
                
                # Add directory path to metadata
                metadata["path"] = analysis_dir
                
                # Add to analyses dictionary
                if key not in analyses:
                    analyses[key] = []
                analyses[key].append(metadata)
        except Exception as e:
            # Skip directories with errors
            continue
            
    # Sort analyses by timestamp (newest first)
    for key in analyses:
        analyses[key].sort(key=lambda x: x.get("timestamp", ""), reverse=True)
        
    return analyses

def delete_analysis(analysis_path):
    """
    Delete a saved analysis.
    
    Args:
        analysis_path (str): Path to the analysis directory
        
    Returns:
        bool: True if deletion was successful, False otherwise
    """
    import shutil
    
    try:
        if os.path.exists(analysis_path) and os.path.isdir(analysis_path):
            shutil.rmtree(analysis_path)
            return True
        else:
            return False
    except Exception as e:
        st.error(f"Error deleting analysis: {str(e)}")
        return False
