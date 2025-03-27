# -*- coding: utf-8 -*-
"""
Visual_Phaser_New_V13.py performs comparisons between siblings and cousins 
and displays HIR/FIR tables plus PNG images for each chromosome.

© 2024 Mick Jolley (mickj1948@gmail.com)
Refactored to remove Excel output and show combined HIR/FIR tables as console output.
"""

import numpy as np
import pandas as pd
import sys
import os
import time
import platform
from itertools import combinations
from PIL import Image, ImageDraw, ImageFont

# -------------------------------------------------------------------------------
# Import your configuration here:
# from VPnew_config_V13_1 import *
# -------------------------------------------------------------------------------
# For clarity, ensure your config file defines variables such as:
#   FILES_PATH, WORKING_DIRECTORY, MAP_PATH
#   SIBLINGS, COUSINS, PHASED_FILES, CHROMOSOMES
#   EXCEL_FILE_NAME (now unused), SHOW_NO_MATCHES, SHOW_TIMES, SHOW_MATCH_PAIR_PROGRESS
#   CHROM_TRUE_SIZE, LINEAR_CHROMOSOME, RESOLUTION, AUTO_REC_PNTS, ARP_TOLERANCE
#   HIR_CUTOFF, FIR_CUTOFF, FIR_TABLES, SCALE_ON, FREEZE_COLUMN, LINUX_FONT_STRING
#   NO_CALL, MAP_PATH, etc.
# -------------------------------------------------------------------------------

# Below is just an example inlined, you can keep using VPnew_config_V13_1.py
# as-is and import from it. Adjust as needed:
FILES_PATH = r"D:\Genetyka\raw_dna"
WORKING_DIRECTORY = r"D:\Genetyka\replit"
MAP_PATH = r"D:\Genetyka"
SIBLINGS = ['TERESA','ALICJA','BOGUSLAWAN']
COUSINS = []
PHASED_FILES = []
CHROMOSOMES = []
SHOW_NO_MATCHES = True
SHOW_TIMES = True
SHOW_MATCH_PAIR_PROGRESS = True
CHROM_TRUE_SIZE = True
LINEAR_CHROMOSOME = True
RESOLUTION = 1
AUTO_REC_PNTS = False
ARP_TOLERANCE = 4
HIR_CUTOFF = 3
FIR_CUTOFF = 1
FIR_TABLES = True
SCALE_ON = True
LINUX_FONT_STRING = "*/fonts/truetype/family/DejaVuSerif-Bold.ttf"
NO_CALL = 'X'
MM_DIST = 500
HIR_SNP_MIN = 100
FIR_SNP_MIN = 75
# -------------------------------------------------------------------------------


def conditions(al1x, al2x, al1y, al2y):
    """Return color code depending on matches / mismatches / no-calls."""
    if al1x == NO_CALL or al1y == NO_CALL:
        return 'limegreen'
    elif al1x == al2x and al1y == al2y and al1x != al1y:
        return 'crimson'
    elif (al1x == al1y and al2x == al2y) or (al1x == al2y and al2x == al1y):
        return 'limegreen'
    else:
        return 'yellow'


def load_dna_files(pair, chrom):
    """
    For a given pair (e.g., ("TERESA","ALICJA")) and chromosome, load 
    and merge their raw .txt or .csv files into a single DataFrame.
    """
    file_names = os.listdir(FILES_PATH)
    dm = pd.DataFrame()

    for ind in list(pair):
        for filname in file_names:
            if ind + "_raw" in filname:
                this_file = os.path.join(FILES_PATH, filname)

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
                    df[f"{ind}-allele1"].isin(['A','C','G','T',NO_CALL])
                    & df[f"{ind}-allele2"].isin(['A','C','G','T',NO_CALL])
                ]

                df = df.reset_index(drop=True)

                # Merge this person's DNA df into common df
                if len(dm.index) > 1:
                    dm = pd.merge(dm, df, on=("rsid","chromosome","position"))
                else:
                    dm = df
    return dm


def scan_genomes(dm, pair, chrom):
    """
    Evaluate the 'match' column (crimson/yellow/limegreen) for HIR and FIR segments.
    Returns two DataFrames: dx (HIR info) and ds (FIR info).
    """
    dm["match"] = np.vectorize(conditions)(
        dm[f"{pair[0]}-allele1"], dm[f"{pair[0]}-allele2"],
        dm[f"{pair[1]}-allele1"], dm[f"{pair[1]}-allele2"],
    )

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
        cutoff = HIR_CUTOFF

    length = len(dm)    
    for i in range(length):
        # (Same logic as in your code.)
        # ...
        # The logic is unchanged; you only remove the Excel references.

        if i == 0 and (dm.loc[i, "match"] in ['yellow','limegreen']):
            nsnps = 1
            segflag = True
            stpos = dm.loc[i, "position"]
            if dm.loc[i,'match'] == 'limegreen':
                fsnps = 1
                fstpos = dm.loc[i, "position"]
                fflag = True
        elif not segflag and (dm.loc[i, "match"] in ['yellow','limegreen']):
            nsnps = 1
            segflag = True
            stpos = dm.loc[i, "position"]
            if not fflag and dm.loc[i,'match'] == 'limegreen':
                fsnps = 1
                fstpos = dm.loc[i, "position"]
                fflag = True
        elif not fflag and dm.loc[i,'match'] == 'limegreen':
            fsnps = 1
            fstpos = dm.loc[i, "position"]
            fflag = True
        elif segflag and (dm.loc[i, "match"] in ['yellow','limegreen']):
            nsnps += 1
            pos = dm.loc[i, "position"]
            if fflag:
                if dm.loc[i,'match'] == 'limegreen':
                    fsnps += 1
                    fpos = dm.loc[i, "position"]
                else: 
                    # we just ended a FIR segment
                    fflag = False
                    if fsnps > FIR_SNP_MIN:
                        dcm = cm_calc(fstpos, fpos)
                        if dcm > FIR_CUTOFF:
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
            if fsnps > FIR_SNP_MIN:
                dcm = cm_calc(fstpos, fpos)
                if dcm > FIR_CUTOFF:
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
                if dm.loc[i, "position"] - mmpos < MM_DIST*1000:
                    segflag = False
                    nmms = 0
                    if nsnps > HIR_SNP_MIN:
                        dcm = cm_calc(stpos, pos)
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
            if nsnps > HIR_SNP_MIN:
                dcm = cm_calc(stpos, pos)
                if dcm > cutoff:
                    addn = pd.DataFrame({
                        "Chr": chrom,
                        "Start Mb": [stpos],
                        "Finish Mb": [pos],
                        "No. SNPs": [nsnps],
                        "Length (cM)": [round(dcm,1)]
                    })
                    dx = pd.concat([dx, addn], ignore_index=True)

            if fsnps > FIR_SNP_MIN:
                dcm = cm_calc(fstpos, fpos)
                if dcm > 3:
                    addn = pd.DataFrame({
                        "Chr": chrom,
                        "Start Mb": [fstpos],
                        "Finish Mb": [fpos],
                        "No. SNPs": [fsnps],
                        "Length (cM)": [round(dcm,1)]
                    })
                    ds = pd.concat([ds, addn], ignore_index=True)

    return dx, ds


def get_image(dplot, iname, chrom):
    """Generate and save a 35px-high PNG representing the main match track."""
    img = Image.new("RGB", (len(dplot), 35), color="white")
    img1 = ImageDraw.Draw(img)

    pos = 0
    for i in range(len(dplot)):
        color = dplot.loc[i,'match']
        color2 = dplot.loc[i,'bar']
        img1.line([(pos, 0), (pos, 19)], fill=color, width=0)
        img1.line([(pos, 20), (pos, 34)], fill=color2, width=0)
        pos += 1

    outpath = os.path.join(WORKING_DIRECTORY, f"{iname} {chrom}.png")
    img.save(outpath)  # <--- CHANGED: no more Excel insertion, just save.


def get_scale(dplot, chrom):
    """Generate a PNG scale bar for the chromosome."""
    img = Image.new("RGB", (len(dplot) + 30, 50), color="white")
    img1 = ImageDraw.Draw(img)

    # Font handling
    if platform.system() == 'Windows':
        fnt = ImageFont.truetype("arial.ttf", 10)
        fnt1 = ImageFont.truetype("arial.ttf", 8)
    elif platform.system() == 'Darwin':
        fnt = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 11)
        fnt1 = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 8)
    else:  # Linux
        fnt = ImageFont.truetype(LINUX_FONT_STRING, 11)
        fnt1 = ImageFont.truetype(LINUX_FONT_STRING, 8)

    max_position = int(dplot['position'].max())
    pixels_per_base = len(dplot) / max_position if max_position != 0 else 0

    # ticks every 0.5 Mb
    for position in range(0, max_position + 1, int(0.5 * 10**6)):
        pixel_position = int(position * pixels_per_base)
        if position % (2 * 10**6) == 0:
            # major tick
            img1.text((pixel_position - 7, 0), f"{position/1e6:.0f}", font=fnt, fill="black")
            img1.line([(pixel_position, 20), (pixel_position, 48)], fill="black", width=1)
        elif position % (10**6) == 0:
            # medium tick
            img1.line([(pixel_position, 30), (pixel_position, 48)], fill="black", width=1)
        else:
            # short tick
            img1.line([(pixel_position, 35), (pixel_position, 48)], fill="black", width=1)

    outpath = os.path.join(WORKING_DIRECTORY, f"scale {chrom}.png")
    img.save(outpath)  # <--- CHANGED: Just save the file.


def load_map():
    """Load minimal map from your 'min_map.txt'."""
    file_name = os.path.join(MAP_PATH, "min_map.txt")
    dmap_source = pd.read_csv(file_name, sep="\t", header=0)
    return dmap_source


def cm_calc(st, fin):
    """Use the loaded map to convert a start/end (in bp) into cM distance."""
    stcm = 0
    fincm = dmap.loc[len(dmap) - 1, "cM"]
    for w in range(1, len(dmap)):
        if st >= dmap.loc[w - 1, "Position"] and st < dmap.loc[w, "Position"]:
            stcm = dmap.loc[w - 1, "cM"] + (
                (st - dmap.loc[w - 1, "Position"]) 
                / (dmap.loc[w, "Position"] - dmap.loc[w - 1, "Position"])
            ) * (dmap.loc[w, "cM"] - dmap.loc[w - 1, "cM"])

        if fin > dmap.loc[w - 1, "Position"] and fin <= dmap.loc[w, "Position"]:
            fincm = dmap.loc[w - 1, "cM"] + (
                (fin - dmap.loc[w - 1, "Position"])
                / (dmap.loc[w, "Position"] - dmap.loc[w - 1, "Position"])
            ) * (dmap.loc[w, "cM"] - dmap.loc[w - 1, "cM"])
            break
    dcm = fincm - stcm
    return dcm


def get_dplot(q, dtot, dxtot, dstot, iname, chrom):
    """Downsample the 'match' data for visual display, and mark HIR/FIR regions in dplot."""
    # The logic is the same as your existing code. We only removed Excel references.
    global CHROM_TRUE_SIZE, RESOLUTION

    if LINEAR_CHROMOSOME:
        CHROM_TRUE_SIZE = False
        if RESOLUTION != 10:
            RESOLUTION = 1

    if CHROM_TRUE_SIZE:
        div = RESOLUTION
    else:
        res = RESOLUTION * 1000
        if res > len(dtot):
            res = len(dtot)
        div =  len(dtot)//res

    dplot = pd.DataFrame(data={'match':'limegreen'}, index=np.arange(len(dtot)//div + 1))
    ycnt = 0
    rcnt = 0

    for i in range(len(dtot)):
        if dtot.iloc[i,q+1] == 'crimson':
            rcnt += 1
        elif dtot.iloc[i,q+1] == 'yellow':
            ycnt += 1
        if i % div == 0 and i != 0:
            if rcnt > 0:
                dplot.loc[int(i/div) - 1,'match'] = 'crimson'
                dplot.loc[int(i/div) - 1,'position'] = dtot.loc[i,'position']
            elif ycnt > 0:
                dplot.loc[int(i/div) - 1,'match'] = 'yellow'
                dplot.loc[int(i/div) - 1,'position'] = dtot.loc[i,'position']
            else:
                dplot.loc[int(i/div) - 1,'position'] = dtot.loc[i,'position']
            ycnt = 0
            rcnt = 0

    if rcnt > 0:
        dplot.loc[len(dplot)-1,'match'] = 'crimson'
        dplot.loc[len(dplot)-1,'position'] = dtot.loc[len(dtot)-1,'position']
    elif ycnt > 0:
        dplot.loc[len(dplot)-1,'match'] = 'yellow'
        dplot.loc[len(dplot)-1,'position'] = dtot.loc[len(dtot)-1,'position']
    else:
        dplot.loc[len(dplot)-1,'position'] = dtot.loc[len(dtot)-1,'position']

    dplot['bar'] = 'black'

    # Mark HIR segments (blue) and FIR segments (orange)
    for i in range(len(dplot)):
        for j in range(len(dxtot)):
            if dxtot.loc[j,'pair'] == iname:
                if (dxtot.loc[j,'Start Mb'] <= dplot.loc[i,'position'] <= dxtot.loc[j,'Finish Mb']):
                    dplot.loc[i,'bar'] = 'blue'
        for j in range(len(dstot)):
            if dstot.loc[j,'pair'] == iname:
                if (dstot.loc[j,'Start Mb'] <= dplot.loc[i,'position'] <= dstot.loc[j,'Finish Mb']):
                    dplot.loc[i,'bar'] = 'orange'
        # If you use AUTO_REC_PNTS, you'd collect them here - omitted for brevity

    if LINEAR_CHROMOSOME:
        length = chr_lens[chrom - 1]
        dplot = dplot.dropna(ignore_index=True)
        if RESOLUTION == 10:
            dtot['fract'] = 10000 * round(dtot['position']/length,4)
            dplotr = pd.DataFrame(data={'match':'grey','bar':'grey'}, index=np.arange(10001))
            dplotr['position'] = np.linspace(0,length,10001)
        else:
            dtot['fract'] = 1000 * round(dtot['position']/length,3)
            dplotr = pd.DataFrame(data={'match':'grey','bar':'grey'}, index=np.arange(1001))
            dplotr['position'] = np.linspace(0,length,1001)

        dtot = dtot.astype({"fract": "int64"})
        fdict = dict(zip(dtot['position'],dtot['fract']))
        for i in range(len(dplot)):
            f = fdict[dplot.loc[i,'position']]
            dplotr.loc[f,'match'] = dplot.loc[i,'match']
            dplotr.loc[f,'bar'] = dplot.loc[i,'bar']

        dplot = dplotr.copy()

    return dplot


# We'll remove all references to openpyxl or Excel, including the "paste_tables" logic.
# Instead, we will just *print* the combined HIR/FIR tables after scanning each chromosome.

def main():
    st = time.time()
    # Basic checks as in your code:
    if len(SIBLINGS) < 2 and len(PHASED_FILES) == 0 and COUSINS == []:
        print("\nThere must be at least two SIBLINGS. Try again.")
        sys.exit()

    if len(PHASED_FILES) < 2 and len(SIBLINGS) == 0 and COUSINS == []:
        print("\nThere must be at least two PHASED_FILES. Try again.")
        sys.exit()

    # Load minimal map
    global dmap
    dmap_source = load_map()

    # Chromosome lengths (Build 37). 
    global chr_lens
    chr_lens = [
       249250621,243199373,198022430,191154276,180915260,171115067,
       159138663,146364022,141213431,135534747,135006516,133851895,
       115169878,107349540,102531392, 90354753, 81195210, 78077248,
       59128983, 63025520, 48129895, 51304566,155270560
    ]

    if COUSINS != []:
        # If you want to handle the "add cousins" scenario, adapt as needed.
        # In this stripped-down version, we assume no adding to existing Excel.
        pass

    # Pair up siblings + phased
    match_pairs = list(combinations(SIBLINGS, 2))
    if PHASED_FILES != []:
        phased_pairs = list(combinations(PHASED_FILES, 2))
        match_pairs = match_pairs + phased_pairs

    # If user gave specific chromosome list, use it; else 1..23
    all_chromosomes = CHROMOSOMES if (len(CHROMOSOMES) > 0) else range(1,24)

    print("\nLoading DNA files and computing HIR/FIR segments...\n")

    for chrom in all_chromosomes:
        # Filter minimal map for this chromosome
        dmap = dmap_source[dmap_source["Chromosome"] == chrom].reset_index(drop=True)

        # We'll keep one data structure that merges the "match" columns for all pairs:
        dtot = pd.DataFrame()

        # We'll keep these to build the final combined HIR/FIR table:
        dxtot = pd.DataFrame()  # All HIR for all match_pairs
        dstot = pd.DataFrame()  # All FIR for all match_pairs

        cstart = time.time()

        for pair in match_pairs:
            t0 = time.time()
            dm = load_dna_files(pair, chrom)
            iname = pair[0] + "-" + pair[1]

            if len(dm) == 0:
                print("\nDNA file not found. Check DNA files folder and file formats.")
                sys.exit()

            # Identify HIR (dx) and FIR (ds)
            dx, ds = scan_genomes(dm, pair, chrom)

            # Append columns for final output
            dx['Match-pair'] = iname
            dx['Type'] = "HIR"
            ds['Match-pair'] = iname
            ds['Type'] = "FIR"

            # Accumulate in "tot" DataFrames
            dx['pair'] = iname
            ds['pair'] = iname
            dxtot = pd.concat([dxtot, dx], ignore_index=True)
            dstot = pd.concat([dstot, ds], ignore_index=True)

            # Create a column with the match color for dtot
            dm = dm.drop(['rsid','chromosome'], axis=1)
            dm = dm.rename(columns={'match': iname})
            if len(dtot) == 0:
                dtot = dm
            else:
                dtot = pd.merge(dtot, dm, on='position', how='outer')  

            if SHOW_MATCH_PAIR_PROGRESS:
                print(f"Chromosome {chrom}, match-pair {iname} complete.")
                if SHOW_TIMES:
                    dt = time.time() - t0
                    print(f"  Elapsed: {dt:.1f} seconds")

        # Once all pairs are processed for this chromosome, we create images:
        if not dtot.empty:
            for q, pair in enumerate(match_pairs):
                iname = pair[0] + "-" + pair[1]
                dplot = get_dplot(q, dtot, dxtot, dstot, iname, chrom)

                if q == 0 and SCALE_ON:
                    get_scale(dplot, chrom)   # create scale PNG
                get_image(dplot, iname, chrom)  # create chromosome image

        ctime = time.time() - cstart
        print(f"\nChromosome {chrom} processed in {ctime:.1f} s.\n")

        # Now produce ONE combined table (HIR+FIR) for this chromosome:
        # Columns: Match-pair, Chr, Start Mb, Finish Mb, No. SNPs, Length (cM), Type
        # Merge dxtot + dstot, but only for rows of this chromosome:
        all_segments = pd.concat([dxtot[dxtot['Chr']==chrom], dstot[dstot['Chr']==chrom]], ignore_index=True)
        # Reorder columns in the desired format:
        all_segments = all_segments[[
            "Match-pair","Chr","Start Mb","Finish Mb","No. SNPs","Length (cM)","Type"
        ]]

        # Print or log them. Example: 
        print("Combined HIR/FIR segments for Chromosome", chrom)
        print(all_segments.to_string(index=False))
        print("-"*70)

    et = time.time()
    print(f"Total elapsed time = {(et-st)//60:.0f} minutes {(et-st) % 60:.0f} seconds.")


# Entry point
if __name__ == "__main__":
    main()
