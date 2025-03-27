"""
Visual Phaser core functionality.
Provides DNA comparison and visualization between siblings and cousins.

Based on original code by Mick Jolley (mickj1948@gmail.com)
Adapted for web use.
"""

import numpy as np
import pandas as pd
import time
import streamlit as st
from io import BytesIO
from PIL import Image, ImageDraw, ImageFont
import os
from itertools import combinations
from utils import cm_calc, conditions, get_chromosome_image, get_scale_image, image_to_base64

def load_dna_files(data_dir, pair, chrom, no_call='X'):
    """
    For a given pair (e.g., ("TERESA","ALICJA")) and chromosome, load 
    and merge their raw .txt or .csv files into a single DataFrame.
    """
    file_names = os.listdir(data_dir)
    dm = pd.DataFrame()

    for ind in list(pair):
        for filname in file_names:
            if ind + "_raw" in filname:
                this_file = os.path.join(data_dir, filname)

                # Handle different file-type patterns:
                if "Ancestry" in filname:
                    df = pd.read_csv(
                        this_file, sep="\t", skip_blank_lines=True,
                        comment="#", header=0,
                        names=["rsid", "chromosome", "position",
                              f"{ind}-allele1", f"{ind}-allele2"],
                    )
                else:
                    df = pd.read_csv(
                        this_file, sep="\t", skip_blank_lines=True,
                        comment="#", header=0, low_memory=False,
                        names=["rsid","chromosome","position",
                              f"{ind}_allele_pair"],
                    )
                    df[f"{ind}-allele1"] = df[f"{ind}_allele_pair"].str[0]
                    df[f"{ind}-allele2"] = df[f"{ind}_allele_pair"].str[1]
                    df.drop([f"{ind}_allele_pair"], axis=1, inplace=True)

                    # Replace X, XY with 23, remove Y, MT
                    df.replace("X","23",inplace=True)
                    df.replace("XY","23",inplace=True)
                    df = df.drop(df[df["chromosome"] == "Y"].index)
                    df = df.drop(df[df["chromosome"] == "MT"].index)
                    df = df.astype({"chromosome":"int64"})

                df = df[df["chromosome"] == chrom]
                df = df[
                    df[f"{ind}-allele1"].isin(['A','C','G','T',no_call])
                    & df[f"{ind}-allele2"].isin(['A','C','G','T',no_call])
                ]

                df = df.reset_index(drop=True)

                # Merge this person's DNA df into common df
                if len(dm.index) > 1:
                    dm = pd.merge(dm, df, on=("rsid","chromosome","position"))
                else:
                    dm = df
    return dm


def scan_genomes(dm, pair, chrom, map_data, config):
    """
    Evaluate the 'match' column (crimson/yellow/limegreen/grey) for HIR and FIR segments.
    Returns two DataFrames: dx (HIR info) and ds (FIR info).
    """
    dm["match"] = np.vectorize(conditions)(
        dm[f"{pair[0]}-allele1"], dm[f"{pair[0]}-allele2"],
        dm[f"{pair[1]}-allele1"], dm[f"{pair[1]}-allele2"],
        config["NO_CALL"]
    )

    # Add visualization bar - initially black (will be updated later)
    dm["bar"] = "black"

    segflag = False
    stpos = 0
    pos = 0
    nmms = 0
    fflag = False
    fsnps = 0
    
    dx = pd.DataFrame()  # Will collect HIR segments
    ds = pd.DataFrame()  # Will collect FIR segments

    # If chromosome 23, use cutoff=15 cM for HIR if not specified.
    if chrom == 23:
        cutoff = 15
    else:
        cutoff = config["HIR_CUTOFF"]

    length = len(dm)    
    for i in range(length):
        if i == 0 and (dm.loc[i, "match"] in ['yellow','limegreen','grey']):
            nsnps = 1
            segflag = True
            stpos = dm.loc[i, "position"]
            if dm.loc[i,'match'] == 'limegreen' or dm.loc[i,'match'] == 'grey':
                fsnps = 1
                fstpos = dm.loc[i, "position"]
                fflag = True
        elif not segflag and (dm.loc[i, "match"] in ['yellow','limegreen','grey']):
            nsnps = 1
            segflag = True
            stpos = dm.loc[i, "position"]
            if not fflag and (dm.loc[i,'match'] == 'limegreen' or dm.loc[i,'match'] == 'grey'):
                fsnps = 1
                fstpos = dm.loc[i, "position"]
                fflag = True
        elif not fflag and (dm.loc[i,'match'] == 'limegreen' or dm.loc[i,'match'] == 'grey'):
            fsnps = 1
            fstpos = dm.loc[i, "position"]
            fflag = True
        elif segflag and (dm.loc[i, "match"] in ['yellow','limegreen','grey']):
            nsnps += 1
            pos = dm.loc[i, "position"]
            if fflag:
                if dm.loc[i,'match'] == 'limegreen' or dm.loc[i,'match'] == 'grey':
                    fsnps += 1
                    fpos = dm.loc[i, "position"]
                else: 
                    # we just ended a FIR segment
                    fflag = False
                    if fsnps > config["FIR_SNP_MIN"]:
                        dcm = cm_calc(fstpos, fpos, map_data)
                        if dcm > config["FIR_CUTOFF"]:
                            addn = pd.DataFrame({
                                "Chr": chrom,
                                "Start Mb": [fstpos],
                                "Finish Mb": [fpos],
                                "No. SNPs": [fsnps],
                                "Length (cM)": [round(dcm,1)]
                            })
                            ds = pd.concat([ds, addn], ignore_index=True)
                    fsnps = 0
        elif segflag and dm.loc[i, "match"] == 'crimson':
            # End of a FIR segment if we were in one
            fflag = False
            if fsnps > config["FIR_SNP_MIN"]:
                dcm = cm_calc(fstpos, fpos, map_data)
                if dcm > config["FIR_CUTOFF"]:
                    addn = pd.DataFrame({
                        "Chr": chrom,
                        "Start Mb": [fstpos],
                        "Finish Mb": [fpos],
                        "No. SNPs": [fsnps],
                        "Length (cM)": [round(dcm,1)]
                    })
                    ds = pd.concat([ds, addn], ignore_index=True)
            fsnps = 0

            nmms += 1
            if nmms == 1:
                mmpos = dm.loc[i, "position"]
            else:
                if dm.loc[i, "position"] - mmpos < config["MM_DIST"]*1000:
                    segflag = False
                    nmms = 0
                    if nsnps > config["HIR_SNP_MIN"]:
                        dcm = cm_calc(stpos, pos, map_data)
                        if dcm > cutoff:
                            addn = pd.DataFrame({
                                "Chr": chrom,
                                "Start Mb": [stpos],
                                "Finish Mb": [pos],
                                "No. SNPs": [nsnps],
                                "Length (cM)": [round(dcm,1)]
                            })
                            dx = pd.concat([dx, addn], ignore_index=True)
                    nsnps = 0
                else:
                    nmms = 1
                    mmpos = dm.loc[i, "position"]

        if i == length - 1:
            # Reached end of chromosome
            if nsnps > config["HIR_SNP_MIN"]:
                dcm = cm_calc(stpos, pos, map_data)
                if dcm > cutoff:
                    addn = pd.DataFrame({
                        "Chr": chrom,
                        "Start Mb": [stpos],
                        "Finish Mb": [pos],
                        "No. SNPs": [nsnps],
                        "Length (cM)": [round(dcm,1)]
                    })
                    dx = pd.concat([dx, addn], ignore_index=True)

            if fsnps > config["FIR_SNP_MIN"]:
                dcm = cm_calc(fstpos, fpos, map_data)
                if dcm > 3:
                    addn = pd.DataFrame({
                        "Chr": chrom,
                        "Start Mb": [fstpos],
                        "Finish Mb": [fpos],
                        "No. SNPs": [fsnps],
                        "Length (cM)": [round(dcm,1)]
                    })
                    ds = pd.concat([ds, addn], ignore_index=True)
    
    # Update the bar colors to show HIR and FIR segments
    # Initially all positions are black (default), now we color HIR blue and FIR orange
    for index, row in dx.iterrows():  # For each HIR segment
        start_pos = row['Start Mb']
        end_pos = row['Finish Mb']
        # Mark all positions within this HIR segment as blue in the bar column
        dm.loc[(dm['position'] >= start_pos) & (dm['position'] <= end_pos), 'bar'] = 'blue'
        
    for index, row in ds.iterrows():  # For each FIR segment
        start_pos = row['Start Mb']
        end_pos = row['Finish Mb']
        # Mark all positions within this FIR segment as orange in the bar column
        dm.loc[(dm['position'] >= start_pos) & (dm['position'] <= end_pos), 'bar'] = 'orange'
    
    # For positions with no match, use grey in both top and bottom
    dm.loc[dm['match'] == 'grey', 'bar'] = 'grey'
    
    # Handle LINEAR_CHROMOSOME mode
    if config["LINEAR_CHROMOSOME"]:
        # Create an initial linear grid with positions along the entire chromosome
        try:
            # These are approximate chromosome lengths - adjust as needed
            chr_lens = [249e6, 243e6, 199e6, 191e6, 182e6, 171e6, 159e6, 146e6, 
                        141e6, 136e6, 135e6, 134e6, 115e6, 107e6, 102e6, 90e6, 
                        84e6, 80e6, 59e6, 64e6, 47e6, 51e6, 157e6] # 23 is X
            
            length = chr_lens[chrom - 1]
            
            # Create linear position grid
            if config["RESOLUTION"] == 10:
                # Higher resolution grid (10001 points)
                linear_grid = pd.DataFrame(index=range(10001))
                linear_grid['position'] = np.linspace(0, length, 10001)
                linear_grid['match'] = 'grey'  # Default color for unfilled positions
                linear_grid['bar'] = 'grey'    # Default color for unfilled positions
            else:
                # Standard resolution grid (1001 points)
                linear_grid = pd.DataFrame(index=range(1001))
                linear_grid['position'] = np.linspace(0, length, 1001)
                linear_grid['match'] = 'grey'  # Default color for unfilled positions
                linear_grid['bar'] = 'grey'    # Default color for unfilled positions
                
            # Map actual data positions to the grid
            # For each data point in dm, find the closest grid position and update colors
            for idx, row in dm.iterrows():
                # Find the closest position in the grid
                pos = row['position']
                grid_idx = int(round(pos / length * (len(linear_grid) - 1)))
                
                # Ensure index is within bounds
                if 0 <= grid_idx < len(linear_grid):
                    # Copy the match and bar colors to the grid
                    linear_grid.loc[grid_idx, 'match'] = row['match']
                    linear_grid.loc[grid_idx, 'bar'] = row['bar']
            
            # Replace original data with the linear grid
            dm = linear_grid
            
        except Exception as e:
            # If there's an error in LINEAR_CHROMOSOME processing, log it but continue
            print(f"Error in LINEAR_CHROMOSOME processing: {str(e)}")
            # We'll still use the original data
    
    return dx, ds, dm


def process_chromosome(data_dir, pair, chrom, map_data, config):
    """Process a single chromosome comparison for a pair of individuals."""
    start_time = time.time()
    
    # Load DNA files for this pair and chromosome
    dm = load_dna_files(data_dir, pair, chrom, config["NO_CALL"])
    if dm.empty:
        # Create placeholder images for empty data with 50px height
        placeholder_chr = Image.new("RGB", (400, 50), color="white")
        placeholder_draw = ImageDraw.Draw(placeholder_chr)
        placeholder_draw.text((10, 20), f"No data found for chromosome {chrom}", fill="black")
        
        placeholder_scale = Image.new("RGB", (400, 50), color="white")
        placeholder_scale_draw = ImageDraw.Draw(placeholder_scale)
        placeholder_scale_draw.text((10, 20), "No scale data available", fill="black")
        
        return {
            "hir_data": None,
            "fir_data": None,
            "match_data": None,
            "chr_image": placeholder_chr,
            "scale_image": placeholder_scale,
            "loading_time": 0,
            "processing_time": 0,
            "status": "No data found for this chromosome"
        }
    
    loading_time = time.time() - start_time
    
    # Scan genomes to find matching segments
    hir_data, fir_data, match_data = scan_genomes(dm, pair, chrom, map_data, config)
    
    # Generate visualization images
    try:
        chr_image = get_chromosome_image(match_data, f"{pair[0]}-{pair[1]}", chrom, config)
    except Exception as e:
        # Create error placeholder image if image generation fails
        chr_image = Image.new("RGB", (400, 50), color="white")
        chr_draw = ImageDraw.Draw(chr_image)
        chr_draw.text((10, 20), f"Error generating chromosome image: {str(e)}", fill="black")
    
    scale_image = None
    if config["SCALE_ON"]:
        try:
            scale_image = get_scale_image(match_data, chrom, config)
        except Exception as e:
            # Create error placeholder image if scale generation fails
            scale_image = Image.new("RGB", (400, 50), color="white")
            scale_draw = ImageDraw.Draw(scale_image)
            scale_draw.text((10, 20), f"Error generating scale: {str(e)}", fill="black")
    
    processing_time = time.time() - start_time - loading_time
    
    # Ensure the images have valid dimensions
    if chr_image.width == 0 or chr_image.height == 0:
        chr_image = Image.new("RGB", (400, 50), color="white")
        chr_draw = ImageDraw.Draw(chr_image)
        chr_draw.text((10, 20), "Invalid chromosome image dimensions", fill="black")
    
    if scale_image and (scale_image.width == 0 or scale_image.height == 0):
        scale_image = Image.new("RGB", (400, 50), color="white")
        scale_draw = ImageDraw.Draw(scale_image)
        scale_draw.text((10, 20), "Invalid scale image dimensions", fill="black")
    
    return {
        "hir_data": hir_data,
        "fir_data": fir_data,
        "match_data": match_data,
        "chr_image": chr_image,
        "scale_image": scale_image,
        "loading_time": loading_time,
        "processing_time": processing_time,
        "status": "success"
    }


def process_all_pairs(data_dir, individuals, chromosomes, map_data, config):
    """Process all pairs of individuals across specified chromosomes."""
    results = {}
    
    # Generate all possible pairs of individuals
    pairs = list(combinations(individuals, 2))
    
    for pair in pairs:
        pair_key = f"{pair[0]}-{pair[1]}"
        results[pair_key] = {}
        
        for chrom in chromosomes:
            if config["SHOW_MATCH_PAIR_PROGRESS"]:
                st.write(f"Processing {pair_key} for chromosome {chrom}...")
                
            chr_results = process_chromosome(data_dir, pair, chrom, map_data, config)
            results[pair_key][chrom] = chr_results
    
    return results
