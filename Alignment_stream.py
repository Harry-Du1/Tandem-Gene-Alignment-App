from collections import defaultdict
import subprocess
import tempfile
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import random
import os
import argparse
import shutil
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib as mpl

global chunk_size
chunk_size = 500
exon_height = 0.05
spacing = 0.02

#######################################################HELPER FUNCTIONS####################################################
def gradient_colors(start_color, end_color, n):
    def hex_to_rgb(hex_color):
        hex_color = hex_color.lstrip('#')
        return [int(hex_color[i:i+2], 16) for i in (0, 2, 4)]

    def rgb_to_hex(rgb_color):
        return '#{:02X}{:02X}{:02X}'.format(*rgb_color)

    start_rgb = hex_to_rgb(start_color)
    end_rgb = hex_to_rgb(end_color)

    gradient = []
    for i in range(n):
        interpolated = [
            int(start_rgb[j] + (end_rgb[j] - start_rgb[j]) * i / (n - 1))
            for j in range(3)
        ]
        gradient.append(rgb_to_hex(interpolated))
    return gradient

def group_mappings(mappings, threshold=100):  # threshold: bp distance to be considered "same"
    groups = []
    used = [False] * len(mappings)

    for i, m in enumerate(mappings):
        if used[i]:
            continue
        group = [i]
        used[i] = True

        for j in range(i + 1, len(mappings)):
            if used[j]:
                continue

            # Compare both query and ref start/end
            q_similar = abs(m[0] - mappings[j][0]) < threshold and abs(m[1] - mappings[j][1]) < threshold
            r_similar = abs(m[2] - mappings[j][2]) < threshold and abs(m[3] - mappings[j][3]) < threshold

            if q_similar or r_similar:
                group.append(j)
                used[j] = True

        groups.append(group)
    return groups

def readFasta(fileName):
    with open(fileName) as lines:
        sequences = []
        currSeq = ""
        for line in lines:
            if line[0] == ">":
                if currSeq:
                    sequences.append(currSeq)
                    currSeq = ""
            else:
                currSeq += line.strip()
        if currSeq:
            sequences.append(currSeq)
        return sequences

def revComp(string):
    result = string[::-1]
    a = 'ATCG'
    b = 'TAGC'
    table = str.maketrans(a,b)
    return result.translate(table)

def writeFasta(seqs, fileName):
    with open(fileName, "w") as f:
        for i, seq in enumerate(seqs):
            f.write(f">query_{i}\n")
            f.write(f"{seq}\n")
def chop_seq(seq):
    c = []
    for i in range(0, len(seq), chunk_size):
        c.append(seq[i:i + chunk_size])
    return c

#######################################################DEFINE SEQUENCE CLASS####################################################
class AdhSeq:
    exon_47 = [
        (126, 8839, "47 dup1"),
        (1754, 3478, "Adh1 gene"),
        (8845, 17599, "47 dup2"),
        (10570, 12431, "Adh1 gene")
    ]

    exon_160i = [
        (657, 7682, "160i dup1"),
        (7683, 14689, "160i dup2"),
        (1689, 3413, "Adh1 gene"),
        (8759, 10620, "Adh1 gene")
    ]

    exon_nov14 = [
        (2206, 4058, "Adh1 mRNA")
    ]
    exon_SC10C = [
        (466, 7519, "SC10C dup1"),
        (1534, 3256, "Adh1 mRNA"),
        (7520, 14573, "SC10C dup2"),
        (8596, 10457, "Adh1 mRNA")
    ]
    exon_09 = [
        (345, 5479, "09 dup1"),
        (1472, 3196, "Adh1 mRNA"),
        (5480, 10573, "09 dup2"),
        (6579, 8441, "Adh1 mRNA")
    ]
    exon_49 = [
        (464, 7519, "49 dup1"),
        (1532, 3256, "Adh1 mRNA"),
        (7520, 14576, "49 dup2"),
        (8596, 10457, "Adh1 mRNA")
    ]
    exon_longNov = [
        (7001, 8853, "Adh1 mRNA")
    ]
    def __init__(self, file, color, tp):
        self.file = file
        self.name = self.file[:-6]
        self.seq = readFasta(self.file)[0]
        self.type = tp
        
        # Generate a list of colors, one per chunk
        if self.type == 'ref':
            self.color = gradient_colors(color[0], color[1], int(len(self.seq) / chunk_size * 20) + 1)
        else:
            self.color = ['#808080'] * (int(len(self.seq) / chunk_size * 20) + 1)
        
        if '47' in self.name:
            self.ex = AdhSeq.exon_47
            if self.type == 'ref':
                for i in range(len(self.color)):
                    if chunk_size * i /20 > 8845 and chunk_size * i /20 < 17599:
                        self.color[i] = self.color[i - int(8845/chunk_size * 20)]
        elif '160i' in self.name:
            self.ex = AdhSeq.exon_160i
            if self.type == 'ref':
                for i in range(len(self.color)):
                    if chunk_size * i / 20 > 7683 and chunk_size * i / 20 < 14689:
                        self.color[i] = self.color[i - int(7683/chunk_size * 20)]
        elif 'nov14' in self.name:
            self.ex = AdhSeq.exon_nov14
        elif 'SC10C' in self.name:
            self.ex = AdhSeq.exon_SC10C
            if self.type == 'ref':
                for i in range(len(self.color)):
                    if chunk_size * i /20 > 7520 and chunk_size * i /20 < 14576 and chunk_size * i /20 < 7520 + 7519 - 466:
                        self.color[i] = self.color[i - int(7520/chunk_size * 20)]
        elif '09' in self.name:
            self.ex = AdhSeq.exon_09
            if self.type == 'ref':
                for i in range(len(self.color)):
                    if chunk_size * i / 20 > 5539 and chunk_size * i /20 < 12555 and chunk_size * i /20 < 5539 + 5538 - 404:
                        self.color[i] = self.color[i - int(5539/chunk_size * 20)]
        elif '49' in self.name:
            self.ex = AdhSeq.exon_49
            if self.type == 'ref':
                for i in range(len(self.color)):
                    if chunk_size * i /20 > 7520 and chunk_size * i /20 < 14576 and chunk_size * i /7520 + 7519 - 464:
                        self.color[i] = self.color[i - int(7520/chunk_size * 20)]
        elif 'long' in self.name:
            self.ex = AdhSeq.exon_longNov
        else:
            self.ex = None

        
    def getSeq(self):
        return self.seq
    def getName(self):
        return self.name
    def getExons(self):
        return self.ex
    def getFile(self):
        return self.file
    def getColor(self):
        return self.color
    def changeColor(self, i, newColor):
        self.color[i] = newColor
        return

    

#To be continued...

#######################################################Run Program####################################################

def _find_lastz():
    # 1) env var
    p = os.environ.get("LASTZ_PATH")
    if p and Path(p).exists():
        return p
    # 2) system PATH
    p = shutil.which("lastz")
    if p:
        return p
    # 3) bundled binary in repo
    here = Path(__file__).resolve().parent
    bundled = here / "bin" / "lastz"
    if bundled.exists():
        return str(bundled)
    return None

def run_lastz(query_file, ref_file, output_file="out.txt"):
    lastz_path = _find_lastz()
    if not lastz_path:
        raise RuntimeError(
            "LASTZ not found. Set LASTZ_PATH env var, install it on PATH, "
            "or include a Linux binary at Alignment_stream/bin/lastz."
        )
    with open(output_file, "w") as out:
        subprocess.run([
            lastz_path,
            f"{ref_file}[multiple]",
            f"{query_file}",
            "--format=general:name1,start1,end1,name2,start2,end2,identity",
            "--strand=both",
            "--maxwordcount=2"
        ], stdout=out, check=True)



def parse_lastz(lastz_file, chunk_size):
    mappings = []
    with open(lastz_file) as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 7 or fields[1] == 'start1':
                continue
            ref_name = fields[0]
            rstart = int(fields[1])
            rend = int(fields[2])
            query_name = fields[3]
            qstart = int(fields[4])
            qend = int(fields[5])
            
            # Fix: parse identity as fraction
            num, denom = map(int, fields[6].split('/'))
            identity = float(num) / float(denom) if denom != 0 else 0

            chunk_idx = int(query_name.split('_')[1])
            qstart_global = qstart + chunk_idx * chunk_size
            qend_global = qend + chunk_idx * chunk_size

            mappings.append((qstart_global, qend_global, rstart, rend, identity))
    return mappings


def stack_exons(exons):
    """
    exons entries can be:
      (start, end, label) or (start, end, label, color_hex)
    Yields: (start, end, label, level, color_or_None)
    """
    sorted_exons = sorted(exons, key=lambda x: x[0])
    levels = []
    for exon in sorted_exons:
        start, end = exon[0], exon[1]
        label = exon[2] if len(exon) > 2 else ""
        color = exon[3] if len(exon) > 3 else None
        placed = False
        for level, track in enumerate(levels):
            if all(start >= t[1] or end <= t[0] for t in track):
                track.append((start, end))
                yield (start, end, label, level, color)
                placed = True
                break
        if not placed:
            levels.append([(start, end)])
            yield (start, end, label, len(levels) - 1, color)



def plot_all(query, ref, query_coords, ref_coords, query_len, ref_len,
             exons_query=None, exons_ref=None, colors=None, identity=None,
             window_size=100):
    
    prev_pos = 0
    rpPos = -1000
    fig, ax = plt.subplots(figsize=(14, 6), dpi=200)

    ax.plot([0, query_len], [1, 1], color='black', linewidth=2, label='Query (Gene A)')
    ax.plot([0, ref_len], [0, 0], color='black', linewidth=2, label='Reference (Gene B)')

    combined = list(zip(query_coords, ref_coords, identity))

    combined.sort(key=lambda x: x[1][0])
    
    for i, (q, r, idt) in enumerate(combined):
        q_start, q_end = q
        r_start, r_end = r

        # Use a random color or pre-defined color
        refColor = ref.getColor()
        queryColor = query.getColor()

        # Adjust transparency based on alignment type
        alpha = 0.8

        rcPos = (r_start + r_end)/2
        qcPos = (q_start + q_end)/2

        currRefColor = refColor[int(rcPos/chunk_size * 20)]

        for i in range(int(q_start/chunk_size * 20), int(q_end/chunk_size * 20)):
            query.changeColor(i, refColor[int(r_start/chunk_size*20) + i - int(q_start/chunk_size*20)])
        
        rpPos = rcPos
        polygon = Polygon([
            (q_start, 1), (q_end, 1),
            (r_end, 0), (r_start, 0)
        ], closed=True, facecolor=currRefColor,
                          edgecolor='black',
                          alpha=alpha)
        ax.add_patch(polygon)
    for i in range(max(len(query.getColor()), len(ref.getColor()))):
        if i < len(query.getColor()):
            if i == len(query.getColor()) -1:
                rect_1 = plt.Rectangle((chunk_size*i / 20, 1), query_len %(chunk_size)/20 + 20, exon_height,
                                     facecolor= query.getColor()[i], #edgecolor='black',
                                 )
            else:
                rect_1 = plt.Rectangle((chunk_size*i / 20, 1), chunk_size / 20, exon_height,
                                     facecolor= query.getColor()[i], #edgecolor='black',
                                 )
        if i < len(ref.getColor()):
            if i == len(ref.getColor()) -1:
                rect_2 = plt.Rectangle((chunk_size * i / 20, -exon_height), ref_len % (chunk_size / 20) + 20, exon_height,
                                         facecolor=ref.getColor()[i], #edgecolor='black',
                                     )
            else:
                rect_2 = plt.Rectangle((chunk_size * i/ 20, -exon_height), chunk_size / 20, exon_height,
                                         facecolor=ref.getColor()[i], #edgecolor='black',
                                     )
            
            
            ax.add_patch(rect_2)
        
        ax.add_patch(rect_1)


    # Query exons
    if exons_query:
        for start, end, label, level, color in stack_exons(exons_query):
            y_base = 1.05 + level * (exon_height + spacing)
            rect = plt.Rectangle((start, y_base + 0.1), end - start, exon_height,
                                 facecolor=(color if color else 'orange'),
                                 edgecolor='black')
            ax.add_patch(rect)
            ax.text((start + end) / 2, y_base + exon_height + 0.105, label,
                    ha='center', va='bottom', fontsize=8)


    # Reference exons
    if exons_ref:
        if exons_ref:
            for start, end, label, level, color in stack_exons(exons_ref):
                y_base = -0.15 - level * (exon_height + spacing)
                rect = plt.Rectangle((start, y_base), end - start, exon_height,
                                     facecolor=(color if color else 'skyblue'),
                                     edgecolor='black')
                ax.add_patch(rect)
                ax.text((start + end) / 2, y_base - 0.01, label,
                        ha='center', va='top', fontsize=8)



    #Divergence Area Graph
    ref_bins = int(ref_len / window_size) + 1
    identity_sums = np.zeros(ref_bins)
    counts = np.zeros(ref_bins)

    # Bin alignment identity
    for (r_start, r_end), idt in zip(ref_coords, identity):
        start_bin = int(r_start / window_size)
        end_bin = int(r_end / window_size)
        for b in range(start_bin, end_bin + 1):
            bin_start = b * window_size
            bin_end = bin_start + window_size
            overlap_start = max(r_start, bin_start)
            overlap_end = min(r_end, bin_end)
            overlap_len = max(0, overlap_end - overlap_start)
            if overlap_len > 0:
                identity_sums[b] += idt * overlap_len
                counts[b] += overlap_len
    

    avg_identity = np.divide(identity_sums, counts, out=np.full_like(identity_sums, np.nan), where=counts != 0)
    divergence = 1 - avg_identity

    divergence_filled = np.nan_to_num(divergence, nan=0.0)

    divergence_smooth = gaussian_filter1d(divergence_filled, sigma=1, mode='nearest', truncate=2.0)

    valid_bins = counts > 0  # numpy boolean array

    divergence_smooth[~valid_bins] = np.nan

    bin_centers = np.arange(ref_bins) * window_size + window_size // 2

    # Define color map
    cmap = plt.cm.plasma
    norm = mpl.colors.Normalize(vmin=0.0, vmax=1-min(identity))
    scale = 1
    if 1-min(identity) < 0.1:
        scale = 3
    # Plot color-mapped segments
    for i in range(1, len(divergence_smooth)):
        if not np.isnan(divergence_smooth[i]) and not np.isnan(divergence_smooth[i - 1]):
            x_vals = [bin_centers[i - 1], bin_centers[i]]
            y_vals = [-0.5 + divergence_smooth[i - 1] * scale,
                      -0.5 + divergence_smooth[i] * scale]
            ax.fill_between(x_vals, [-0.5, -0.5], y_vals,
                            color=cmap(norm((divergence_smooth[i] + divergence_smooth[i - 1]) / 2)),
                            alpha=0.8,
                            step='pre')

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.015, pad=0.01)
    cbar.set_label('Sequence Divergence (1 - Identity)', fontsize=8)

        ### --- Query Divergence Track --- ###
    query_bins = int(query_len / window_size) + 1
    q_identity_sums = np.zeros(query_bins)
    q_counts = np.zeros(query_bins)

    for (q_start, q_end), idt in zip(query_coords, identity):
        q_start_bin = int(q_start / window_size)
        q_end_bin = int(q_end / window_size)
        for b in range(q_start_bin, q_end_bin + 1):
            q_identity_sums[b] += idt
            q_counts[b] += 1

    q_avg_identity = np.divide(q_identity_sums, q_counts, out=np.full_like(q_identity_sums, np.nan), where=q_counts != 0)
    q_divergence = 1 - q_avg_identity
    q_divergence_filled = np.nan_to_num(q_divergence, nan=0.0)
    q_divergence_smooth = gaussian_filter1d(q_divergence_filled, sigma=1, mode='nearest', truncate=2.0)
    q_valid_bins = q_counts > 0 
    q_divergence_smooth[~q_valid_bins] = np.nan
    q_bin_centers = np.arange(query_bins) * window_size + window_size // 2
    q_baseline = 1.5  # slightly above top of query track
    q_scale = 1
    if 1-min(identity) < 0.1:
        q_scale = 3
    for i in range(1, len(q_divergence_smooth)):
        if not np.isnan(q_divergence_smooth[i]) and not np.isnan(q_divergence_smooth[i - 1]):
            x_vals = [q_bin_centers[i - 1], q_bin_centers[i]]
            y_vals = [q_baseline - q_divergence_smooth[i - 1] * q_scale,
                      q_baseline - q_divergence_smooth[i] * q_scale]
            ax.fill_between(x_vals, [q_baseline, q_baseline], y_vals,
                            color=cmap(norm((q_divergence_smooth[i] + q_divergence_smooth[i - 1]) / 2)),
                            alpha=0.9,
                            step='pre')
            
    # Axes and labels
    ax.set_ylim(-1.2, 2.0)
    ax.set_xlim(0, 18000)
    ax.set_yticks([0, 1, -0.5, 1.5])
    ax.set_yticklabels([ref.getName(), query.getName(), "Referece Divergence", "Query Divergence"])
    ax.set_xlabel("Position (bp)")
    ax.set_title("Sequence Alignment")
    ax.grid(True, axis='x', linestyle='--', linewidth=0.5)
    plt.tight_layout()

    return fig

def getGraph(query, ref, group_threshold=100, window_size=100):
    query_fasta = query.getFile()
    ref_fasta = ref.getFile()

    temp_query_file = "chopped_query.fasta"
    paf_file = "alignments_lastz.txt"

    # Prepare data
    seq = query.getSeq()
    chopped_seqs = chop_seq(seq)

    writeFasta(chopped_seqs, temp_query_file)
    run_lastz(temp_query_file, ref_fasta, paf_file)
    mappings = parse_lastz(paf_file, chunk_size)

    print(f"Found {len(mappings)} mappings.")

    exons_query = query.getExons()
    exons_ref = ref.getExons()

    # Extract coordinates
    # Extract coordinates
    query_coords = [(m[0], m[1]) for m in mappings]
    ref_coords = [(m[2], m[3]) for m in mappings]
    identity = [m[4] for m in mappings]

    # Group mappings by similarity
    groups = group_mappings(mappings, threshold=group_threshold)
    group_colors = gradient_colors("#FF6B6B", "#4D96FF", len(groups))

    # Assign a color to each mapping based on group
    color_map = [None] * len(mappings)
    for i, group in enumerate(groups):
        for idx in group:
            color_map[idx] = group_colors[i]



    return plot_all(query, ref, query_coords, ref_coords, len(seq), len(readFasta(ref_fasta)[0]),
                exons_query, exons_ref, color_map, identity, window_size=window_size)
#######################################################MAIN EXECUTION####################################################

