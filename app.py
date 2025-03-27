"""
Visual Phaser Web Application

This Streamlit app provides a web interface for the Visual Phaser DNA analysis tool.
It allows users to upload DNA files, configure analysis parameters, and view results.

Based on original code by Mick Jolley (mickj1948@gmail.com)
Adapted for web use.
"""

import streamlit as st
import pandas as pd
import numpy as np
import os
import time
import platform
from itertools import combinations
from PIL import Image, ImageDraw, ImageFont
from io import BytesIO
import base64

from config import DEFAULT_CONFIG, DATA_DIR, MAP_DIR, DEFAULT_MAP_FILE
from utils import (
    get_dna_files, extract_individual_names, load_map_data, 
    save_uploaded_file, image_to_base64, save_analysis_results,
    load_analysis_results, get_saved_analyses, delete_analysis
)
from visual_phaser import process_all_pairs, process_chromosome

# Page configuration
st.set_page_config(
    page_title="Visual Phaser DNA Analysis",
    page_icon="🧬",
    layout="wide",
)

# Create necessary directories if they don't exist
for directory in [DATA_DIR, MAP_DIR]:
    if not os.path.exists(directory):
        os.makedirs(directory)

# Initialize session state for configuration
if 'config' not in st.session_state:
    st.session_state.config = DEFAULT_CONFIG.copy()

if 'map_data' not in st.session_state:
    # Try to load default map file if it exists
    if os.path.exists(DEFAULT_MAP_FILE):
        st.session_state.map_data = load_map_data(DEFAULT_MAP_FILE)
    else:
        st.session_state.map_data = None

# Check for existing DNA files
if 'dna_files_uploaded' not in st.session_state:
    # Check if there are already DNA files in the data directory
    existing_files = get_dna_files(DATA_DIR)
    st.session_state.dna_files_uploaded = len(existing_files) > 0
    
    # Store the list of existing DNA files
    if st.session_state.dna_files_uploaded:
        st.session_state.existing_dna_files = existing_files

if 'analysis_complete' not in st.session_state:
    st.session_state.analysis_complete = False

if 'results' not in st.session_state:
    st.session_state.results = {}
    
# Initialize saved analyses state
if 'saved_analyses' not in st.session_state:
    st.session_state.saved_analyses = get_saved_analyses()
    
# For toggling between current and saved analyses
if 'view_saved' not in st.session_state:
    st.session_state.view_saved = False
    
# Store current tab selection for saved analyses
if 'saved_tab_index' not in st.session_state:
    st.session_state.saved_tab_index = {}

# Title and introduction
st.title("Visual Phaser DNA Analysis")
st.markdown("""
This application performs DNA comparison between siblings and cousins,
displaying HIR/FIR tables and chromosome visualizations.

**Top-Bar Colors:**
- **Limegreen**: Fully Identical Regions (FIR) - matching segments on both chromosomes or at least one sibling has a "no call"
- **Crimson**: Half Identical Regions (HIR) - both siblings are homozygous but for different alleles
- **Yellow**: Half-match scenario - siblings share one allele but differ in the other
- **Grey**: Placeholder color for unfilled positions in linear chromosome mode

**Bottom-Bar Colors:**
- **Blue**: Position found within an HIR (half-identical region) segment
- **Orange**: Position found within an FIR (fully-identical region) segment
- **Black**: Default background if not within an HIR or FIR segment
- **Grey**: Placeholder color for unfilled positions in linear chromosome mode

© 2024 Mick Jolley (mickj1948@gmail.com) - Original code
Web adaptation - Streamlit version
""")

# Sidebar for configuration
with st.sidebar:
    st.header("Configuration")
    
    # File Upload Section
    st.subheader("DNA Files")
    
    # Show current DNA files if any exist
    if st.session_state.dna_files_uploaded:
        available_files = get_dna_files(DATA_DIR)
        
        if available_files:
            st.info(f"{len(available_files)} DNA files currently loaded:")
            file_list = "<ul>"
            for file in available_files:
                # Extract just the filename from the path
                filename = os.path.basename(file)
                file_list += f"<li>{filename}</li>"
            file_list += "</ul>"
            st.markdown(file_list, unsafe_allow_html=True)
    
    # File uploader for new files
    st.subheader("Upload Additional DNA Files")
    dna_files = st.file_uploader(
        "Upload raw DNA files (.txt or .csv)",
        type=["txt", "csv"],
        accept_multiple_files=True
    )
    
    if dna_files:
        for file in dna_files:
            save_uploaded_file(file, DATA_DIR)
            st.session_state.dna_files_uploaded = True
        st.success(f"{len(dna_files)} DNA files uploaded successfully")
        # Force refresh of available files
        st.rerun()
    
    # Map file upload
    st.subheader("Upload Genetic Map")
    
    # Display status of current map
    if st.session_state.map_data is not None:
        if os.path.exists(DEFAULT_MAP_FILE):
            st.info(f"Using default genetic map: {DEFAULT_MAP_FILE}")
        else:
            st.info("Using uploaded genetic map")
    else:
        st.warning("No genetic map loaded. Please upload one or a default 1cM/1Mb approximation will be used.")
        
    map_file = st.file_uploader(
        "Upload min_map.txt file",
        type=["txt"]
    )
    
    if map_file:
        map_path = save_uploaded_file(map_file, MAP_DIR, "min_map.txt")
        st.session_state.map_data = load_map_data(map_path)
        if st.session_state.map_data is not None:
            st.success("Genetic map uploaded successfully")
    
    # Individual Selection
    if st.session_state.dna_files_uploaded:
        available_files = get_dna_files(DATA_DIR)
        individual_names = extract_individual_names(available_files)
        
        st.subheader("Individual Selection")
        selected_siblings = st.multiselect(
            "Select siblings to compare",
            options=individual_names,
            default=individual_names[:min(3, len(individual_names))]
        )
        
        selected_cousins = st.multiselect(
            "Select cousins to compare (optional)",
            options=individual_names,
            default=[]
        )
        
        st.session_state.config["SIBLINGS"] = selected_siblings
        st.session_state.config["COUSINS"] = selected_cousins
    
    # Chromosome Selection
    st.subheader("Chromosome Selection")
    all_chroms = list(range(1, 24))  # Chromosomes 1-23 (23 is X)
    selected_chroms = st.multiselect(
        "Select chromosomes to analyze (empty = all)",
        options=all_chroms,
        default=[]
    )
    
    if not selected_chroms:  # If empty, select all
        selected_chroms = all_chroms
    
    st.session_state.config["CHROMOSOMES"] = selected_chroms
    
    # Analysis Parameters
    st.subheader("Analysis Parameters")
    
    st.session_state.config["SHOW_NO_MATCHES"] = st.checkbox(
        "Show No Matches", 
        value=st.session_state.config["SHOW_NO_MATCHES"]
    )
    
    st.session_state.config["CHROM_TRUE_SIZE"] = st.checkbox(
        "Chromosome True Size", 
        value=st.session_state.config["CHROM_TRUE_SIZE"]
    )
    
    st.session_state.config["LINEAR_CHROMOSOME"] = st.checkbox(
        "Linear Chromosome", 
        value=st.session_state.config["LINEAR_CHROMOSOME"]
    )
    
    st.session_state.config["RESOLUTION"] = st.slider(
        "Resolution", 
        min_value=1, 
        max_value=100, 
        value=st.session_state.config["RESOLUTION"]
    )
    
    st.session_state.config["HIR_CUTOFF"] = st.slider(
        "HIR Cutoff (cM)", 
        min_value=1, 
        max_value=15, 
        value=st.session_state.config["HIR_CUTOFF"]
    )
    
    st.session_state.config["FIR_CUTOFF"] = float(st.slider(
        "FIR Cutoff (cM)", 
        min_value=0.1, 
        max_value=5.0, 
        value=float(st.session_state.config["FIR_CUTOFF"]),
        step=0.1
    ))
    
    st.session_state.config["FIR_TABLES"] = st.checkbox(
        "Show FIR Tables", 
        value=st.session_state.config["FIR_TABLES"]
    )
    
    st.session_state.config["SCALE_ON"] = st.checkbox(
        "Show Scale", 
        value=st.session_state.config["SCALE_ON"]
    )
    
    # Run Analysis Button
    if st.session_state.dna_files_uploaded and st.session_state.config["SIBLINGS"]:
        if st.button("Run Analysis"):
            with st.spinner("Running DNA analysis..."):
                # Process all pairs across all selected chromosomes
                results = process_all_pairs(
                    DATA_DIR, 
                    st.session_state.config["SIBLINGS"],
                    st.session_state.config["CHROMOSOMES"],
                    st.session_state.map_data,
                    st.session_state.config
                )
                
                st.session_state.results = results
                st.session_state.analysis_complete = True
                st.success("Analysis complete!")
                st.rerun()

# Main content area
if not st.session_state.dna_files_uploaded:
    st.info("Please upload DNA files using the sidebar to begin.")
elif not st.session_state.config["SIBLINGS"]:
    st.warning("Please select at least two siblings to compare in the sidebar.")
elif st.session_state.analysis_complete or len(st.session_state.saved_analyses) > 0:
    # Toggle between current and saved analyses
    view_options = ["Current Analysis"]
    if len(st.session_state.saved_analyses) > 0:
        view_options.append("Saved Analyses")
        
    col1, col2 = st.columns([3, 1])
    with col1:
        st.header("Analysis Results")
    
    with col2:
        if len(view_options) > 1:
            selected_view = st.radio("View", view_options, horizontal=True)
            st.session_state.view_saved = selected_view == "Saved Analyses"
        else:
            st.session_state.view_saved = False
    
    # Function to display current analysis results
    def display_current_analysis():
        # Get all unique chromosomes across all pairs
        all_chromosomes = set()
        for pair in st.session_state.results.keys():
            all_chromosomes.update(st.session_state.results[pair].keys())
        
        # Sort chromosomes numerically
        all_chromosomes = sorted(all_chromosomes, key=lambda x: int(x))
        
        # Create tabs for each chromosome
        chrom_tabs = st.tabs([f"Chromosome {chrom}" for chrom in all_chromosomes])
        
        for i, chrom in enumerate(all_chromosomes):
            with chrom_tabs[i]:
                st.subheader(f"Chromosome {chrom} Analysis")
                
                # Add Save button for each chromosome
                st.button(f"💾 Save Chromosome {chrom} Analysis", key=f"save_chr_{chrom}", 
                          on_click=lambda c=chrom: save_current_analysis(c))
                
                # Create combined tables for all pairs for this chromosome
                all_hir_data = []
                all_fir_data = []
                
                # Combine data from all pairs for this chromosome
                for pair in st.session_state.results.keys():
                    if chrom in st.session_state.results[pair]:
                        results = st.session_state.results[pair][chrom]
                        
                        if results["status"] == "success":
                            # Add HIR data to combined table
                            if results["hir_data"] is not None and not results["hir_data"].empty:
                                # Add pair info to the HIR data
                                hir_with_pair = results["hir_data"].copy()
                                hir_with_pair["Match-pair"] = pair
                                hir_with_pair["Type"] = "HIR"
                                all_hir_data.append(hir_with_pair)
                            
                            # Add FIR data to combined table if enabled
                            if st.session_state.config["FIR_TABLES"]:
                                if results["fir_data"] is not None and not results["fir_data"].empty:
                                    # Add pair info to the FIR data
                                    fir_with_pair = results["fir_data"].copy()
                                    fir_with_pair["Match-pair"] = pair
                                    fir_with_pair["Type"] = "FIR"
                                    all_fir_data.append(fir_with_pair)
                
                # Create combined dataframe for all segments
                all_segments = []
                if all_hir_data:
                    all_segments.extend(all_hir_data)
                if all_fir_data:
                    all_segments.extend(all_fir_data)
                
                # Display combined table
                if all_segments:
                    combined_df = pd.concat(all_segments, ignore_index=True)
                    # Reorder columns to put Match-pair first
                    cols = ["Match-pair", "Chr", "Start Mb", "Finish Mb", "No. SNPs", "Length (cM)", "Type"]
                    combined_df = combined_df[cols]
                    st.subheader("Segments")
                    st.dataframe(combined_df, use_container_width=True)
                else:
                    st.info(f"No segments found for chromosome {chrom}")
                
                # Display scale if enabled (before chromosome visualizations)
                if st.session_state.config["SCALE_ON"]:
                    st.subheader("Scale")
                    # Just show one scale (from the first pair with valid data)
                    for pair in st.session_state.results.keys():
                        if chrom in st.session_state.results[pair] and st.session_state.results[pair][chrom]["status"] == "success":
                            results = st.session_state.results[pair][chrom]
                            if results["scale_image"] is not None:
                                try:
                                    scale_byte_arr = BytesIO()
                                    results["scale_image"].save(scale_byte_arr, format='PNG')
                                    scale_byte_arr.seek(0)
                                    st.image(scale_byte_arr, caption="Scale (Mb)")
                                    break
                                except Exception as e:
                                    st.error(f"Error displaying scale: {str(e)}")
                                    continue
                
                # Display chromosome visualization for each pair
                st.subheader("Chromosome Visualization")
                
                for pair in st.session_state.results.keys():
                    if chrom in st.session_state.results[pair]:
                        results = st.session_state.results[pair][chrom]
                        
                        if results["status"] == "success":
                            try:
                                # Display pair name
                                st.markdown(f"**{pair}**")
                                
                                # Convert PIL Image to bytes for Streamlit
                                if results["chr_image"] is not None:
                                    # Create a BytesIO object
                                    img_byte_arr = BytesIO()
                                    # Save the image to BytesIO object as PNG
                                    results["chr_image"].save(img_byte_arr, format='PNG')
                                    # Set the position to the beginning of BytesIO object
                                    img_byte_arr.seek(0)
                                    # Display the image from bytes
                                    st.image(img_byte_arr, caption=f"{pair} - Chromosome {chrom}")
                            except Exception as e:
                                st.error(f"Error displaying images for {pair}: {str(e)}")
                        else:
                            st.warning(f"{pair}: {results['status']}")
    
    # Function to save current analysis for a chromosome
    def save_current_analysis(chrom):
        saved_pairs = []
        
        # Save each pair's results for this chromosome
        for pair in st.session_state.results.keys():
            if chrom in st.session_state.results[pair]:
                results = st.session_state.results[pair][chrom]
                
                if results["status"] == "success":
                    # Save this pair's analysis
                    analysis_id = save_analysis_results(
                        results, 
                        int(chrom), 
                        [pair.split('-')[0], pair.split('-')[1]], 
                        st.session_state.config
                    )
                    saved_pairs.append(pair)
        
        # Refresh saved analyses list
        st.session_state.saved_analyses = get_saved_analyses()
        
        if saved_pairs:
            st.success(f"Analysis for Chromosome {chrom} saved successfully for pairs: {', '.join(saved_pairs)}")
        else:
            st.warning(f"No successful analyses to save for Chromosome {chrom}")
            
    # Function to display saved analyses
    def display_saved_analyses():
        # Refresh the saved analyses list
        st.session_state.saved_analyses = get_saved_analyses()
        
        if not st.session_state.saved_analyses:
            st.info("No saved analyses found.")
            return
            
        # Group saved analyses by chromosome
        chr_groups = {}
        for key, analyses in st.session_state.saved_analyses.items():
            # Key format is chr{num}_{pair1}-{pair2}
            parts = key.split('_')
            chr_num = parts[0].replace('chr', '')
            
            if chr_num not in chr_groups:
                chr_groups[chr_num] = {}
                
            chr_groups[chr_num][key] = analyses
                
        # Sort chromosomes numerically
        sorted_chroms = sorted(chr_groups.keys(), key=lambda x: int(x))
        
        # Create tabs for each chromosome
        chrom_tabs = st.tabs([f"Chromosome {chrom}" for chrom in sorted_chroms])
        
        for i, chrom in enumerate(sorted_chroms):
            with chrom_tabs[i]:
                # For each chromosome, group analyses by date
                analyses_by_date = {}
                
                # Collect all analyses for this chromosome across all pairs
                for pair_key, pair_analyses in chr_groups[chrom].items():
                    for analysis in pair_analyses:
                        # Use full timestamp instead of just date
                        timestamp = analysis.get("date", "Unknown Date")
                        
                        if timestamp not in analyses_by_date:
                            analyses_by_date[timestamp] = []
                            
                        analyses_by_date[timestamp].append({
                            "metadata": analysis,
                            "pair_key": pair_key
                        })
                
                # Sort dates with newest first
                sorted_dates = sorted(analyses_by_date.keys(), reverse=True)
                
                # Skip if no analyses
                if not sorted_dates:
                    st.info(f"No saved analyses for Chromosome {chrom}")
                    continue
                    
                # Create tabs for each analysis timestamp
                # Format tabs to show date and time (e.g., "Analysis 2025-03-27 14:30")
                date_tabs = st.tabs([f"Analysis {timestamp}" for timestamp in sorted_dates])
                
                for j, date in enumerate(sorted_dates):
                    with date_tabs[j]:
                        analyses_data = analyses_by_date[date]
                        
                        # Add a delete button for this analysis date
                        if st.button("🗑️ Delete This Analysis", key=f"delete_chr{chrom}_date{date}"):
                            # Delete all analyses for this date
                            deleted_count = 0
                            for analysis_info in analyses_data:
                                if delete_analysis(analysis_info["metadata"]["path"]):
                                    deleted_count += 1
                                    
                            if deleted_count > 0:
                                st.success(f"Deleted {deleted_count} analyses successfully.")
                                # Refresh saved analyses
                                st.session_state.saved_analyses = get_saved_analyses()
                                st.rerun()
                            else:
                                st.error("Error deleting analyses.")
                        
                        # Load the first analysis to get common parameters
                        first_analysis = analyses_data[0]
                        metadata, config, _ = load_analysis_results(first_analysis["metadata"]["path"])
                        
                        if metadata and config:
                            # Display analysis date and parameters once at the top
                            st.markdown(f"**Analysis Date:** {metadata.get('date', 'Unknown')}")
                            st.markdown("**Analysis Parameters:**")
                            
                            param_cols = st.columns(2)
                            with param_cols[0]:
                                st.markdown(f"- HIR Cutoff: {config.get('HIR_CUTOFF', 'N/A')} cM")
                                st.markdown(f"- FIR Cutoff: {config.get('FIR_CUTOFF', 'N/A')} cM")
                                st.markdown(f"- Resolution: {config.get('RESOLUTION', 'N/A')}")
                                
                            with param_cols[1]:
                                st.markdown(f"- Linear Chromosome: {config.get('LINEAR_CHROMOSOME', False)}")
                                st.markdown(f"- Chromosome True Size: {config.get('CHROM_TRUE_SIZE', False)}")
                                st.markdown(f"- Show No Matches: {config.get('SHOW_NO_MATCHES', False)}")
                        
                        # Collect all segments for this chromosome
                        all_hir_data = []
                        all_fir_data = []
                        
                        # Load data from each analysis
                        for analysis_info in analyses_data:
                            _, _, results = load_analysis_results(analysis_info["metadata"]["path"])
                            
                            if results:
                                # Get the full match-pair name from the pair_key
                                # Format: chr1_Name1-Name2 -> we want to preserve the full Name1-Name2 format
                                pair_name = analysis_info["pair_key"].split('_', 1)[1] if '_' in analysis_info["pair_key"] else analysis_info["pair_key"]
                                
                                # Add HIR data
                                if results.get("hir_data") is not None and not results["hir_data"].empty:
                                    hir_with_pair = results["hir_data"].copy()
                                    hir_with_pair["Match-pair"] = pair_name
                                    hir_with_pair["Type"] = "HIR"
                                    all_hir_data.append(hir_with_pair)
                                    
                                # Add FIR data
                                if config and config.get("FIR_TABLES", False):
                                    if results and results.get("fir_data") is not None and not results["fir_data"].empty:
                                        fir_with_pair = results["fir_data"].copy()
                                        fir_with_pair["Match-pair"] = pair_name
                                        fir_with_pair["Type"] = "FIR"
                                        all_fir_data.append(fir_with_pair)
                        
                        # Create combined dataframe for all segments
                        all_segments = []
                        if all_hir_data:
                            all_segments.extend(all_hir_data)
                        if all_fir_data:
                            all_segments.extend(all_fir_data)
                        
                        # Display combined table
                        if all_segments:
                            combined_df = pd.concat(all_segments, ignore_index=True)
                            # Reorder columns to put Match-pair first
                            cols = ["Match-pair", "Chr", "Start Mb", "Finish Mb", "No. SNPs", "Length (cM)", "Type"]
                            cols = [col for col in cols if col in combined_df.columns]
                            combined_df = combined_df[cols]
                            st.subheader("Segments")
                            st.dataframe(combined_df, use_container_width=True)
                        else:
                            st.info(f"No segments found for chromosome {chrom}")
                        
                        # Display scale if enabled (only once at the top)
                        if config and config.get("SCALE_ON", False):
                            st.subheader("Scale")
                            for analysis_info in analyses_data:
                                _, _, results = load_analysis_results(analysis_info["metadata"]["path"])
                                if results and results.get("scale_image") is not None:
                                    try:
                                        scale_byte_arr = BytesIO()
                                        results["scale_image"].save(scale_byte_arr, format='PNG')
                                        scale_byte_arr.seek(0)
                                        st.image(scale_byte_arr, caption="Scale (Mb)")
                                        break  # Show just one scale
                                    except Exception as e:
                                        continue  # Try the next one if there's an error
                        
                        # Display chromosome visualization for each pair
                        st.subheader("Chromosome Visualization")
                        
                        for analysis_info in analyses_data:
                            # Get the full match-pair name from the pair_key
                            pair_name = analysis_info["pair_key"].split('_', 1)[1] if '_' in analysis_info["pair_key"] else analysis_info["pair_key"]
                            _, _, results = load_analysis_results(analysis_info["metadata"]["path"])
                            
                            if results and results.get("chr_image") is not None:
                                try:
                                    # Display pair name
                                    st.markdown(f"**{pair_name}**")
                                    
                                    # Create a BytesIO object
                                    img_byte_arr = BytesIO()
                                    # Save the image to BytesIO object as PNG
                                    results["chr_image"].save(img_byte_arr, format='PNG')
                                    # Set the position to the beginning of BytesIO object
                                    img_byte_arr.seek(0)
                                    # Display the image from bytes
                                    st.image(img_byte_arr, caption=f"{pair_name} - Chromosome {chrom}")
                                except Exception as e:
                                    st.error(f"Error displaying image for {pair_name}: {str(e)}")
                            else:
                                st.warning(f"No visualization available for {pair_name}")
    
    # Display either current or saved analyses based on the toggle
    if st.session_state.view_saved:
        display_saved_analyses()
    else:
        if st.session_state.analysis_complete:
            display_current_analysis()
        else:
            st.info("No current analysis results. Please run an analysis or view saved results.")
            # Automatically switch to saved analyses if available
            if len(st.session_state.saved_analyses) > 0:
                st.session_state.view_saved = True
                st.rerun()
else:
    st.info("Configure your analysis parameters in the sidebar, then click 'Run Analysis'")

# Footer
st.markdown("---")
st.markdown("Visual Phaser DNA Analysis Tool - Web Edition")
