import io
from pathlib import Path
import streamlit as st
import Alignment_stream as AS
from Alignment_stream import AdhSeq, getGraph
import pandas as pd

st.set_page_config(page_title="Tandem Duplicate Alignment Visualizer", layout="wide")

# ---------------- Session State Init ----------------
st.session_state.setdefault("ann_q", [])  # list of dicts: {name,start,end,color}
st.session_state.setdefault("ann_r", [])

# --- Color Picker Column compatibility (older Streamlit versions may not have it) ---
HAS_COLOR_PICKER_COLUMN = hasattr(st.column_config, "ColorPickerColumn")
if HAS_COLOR_PICKER_COLUMN:
    COLOR_COL_Q = st.column_config.ColorPickerColumn("Color", help="Pick or paste hex", default="#FFA500")
    COLOR_COL_R = st.column_config.ColorPickerColumn("Color", help="Pick or paste hex", default="#87CEEB")
else:
    COLOR_COL_Q = st.column_config.TextColumn("Color (hex)", help="Enter hex like #FFA500")
    COLOR_COL_R = st.column_config.TextColumn("Color (hex)", help="Enter hex like #87CEEB")
    st.info("This Streamlit version doesn't support ColorPickerColumn; using hex text fields instead.")

# ---------------- Title & Intro ----------------
st.title("Tandem Duplicate Alignment Visualizer")
st.markdown(
    "Upload or paste FASTA for **Query** and **Reference**, then explore the alignment with adjustable colors and parameters. "
    "Use **Add Annotation** to draw custom blocks (e.g., exons) with labels and colors."
)

# ---------------- Sidebar: Appearance ----------------
st.sidebar.header("Appearance")
palette = st.sidebar.selectbox(
    "Color gradient",
    ["Sunset (red→blue)", "Teal→Purple", "Green→Orange", "Custom"],
    key="palette_select"
)
PRESETS = {
    "Sunset (red→blue)": ("#FF6B6B", "#4D96FF"),
    "Teal→Purple":       ("#00BFA6", "#7E57C2"),
    "Green→Orange":      ("#00C853", "#FF6D00"),
}
if palette == "Custom":
    start_hex = st.sidebar.color_picker("Start color", "#FF6B6B", key="color_start")
    end_hex   = st.sidebar.color_picker("End color",   "#4D96FF", key="color_end")
    GRADIENT = (start_hex, end_hex)
else:
    GRADIENT = PRESETS[palette]

# DPI & render method
fig_dpi = st.sidebar.slider("Figure DPI", 100, 400, 200, 50, key="fig_dpi")
use_png_buffer = st.sidebar.checkbox("Render via high-DPI PNG (sharper)", value=False, key="use_png_buffer")

# ---------------- Sidebar: Parameters ----------------
st.sidebar.header("Parameters")
# set module global so chopping & color indexing stay consistent
AS.chunk_size = st.sidebar.slider("Chunk size (bp)", 100, 2000, AS.chunk_size, 50, key="chunk_size")
group_threshold = st.sidebar.slider("Grouping threshold (bp)", 10, 1000, 100, 10, key="group_thresh")
window_size = st.sidebar.slider("Divergence window (bp)", 50, 1000, 100, 50, key="window_size")

# ---------------- Input Method ----------------
st.header("Input Sequences")
method = st.selectbox("Input method", ["Upload", "Paste FASTA text"], key="input_method")

query_file = None
ref_file = None
query_text = ""
ref_text = ""

if method == "Upload":
    query_file = st.file_uploader("Upload Query FASTA", type=["fasta", "fa", "fna"], key="query_upload")
    ref_file   = st.file_uploader("Upload Reference FASTA", type=["fasta", "fa", "fna"], key="ref_upload")
else:
    query_text = st.text_area("Paste Query FASTA", key="query_text")
    ref_text   = st.text_area("Paste Reference FASTA", key="ref_text")

# ---------------- Add Annotation ----------------
with st.expander("➕ Add Annotation", expanded=False):
    # defaults for names (used in dropdowns before ready)
    default_q_name = Path(query_file.name).stem if (method == "Upload" and query_file is not None) else "Query"
    default_r_name = Path(ref_file.name).stem if (method == "Upload" and ref_file is not None) else "Reference"

    seq_choice = st.selectbox("Annotate which sequence?", [default_q_name, default_r_name], key="ann_seq_choice")
    ann_name = st.text_input("Annotation label", value="exon", key="ann_label")
    c1, c2 = st.columns(2)
    with c1:
        ann_start = st.number_input("Start (bp)", min_value=0, value=0, step=1, key="ann_start")
    with c2:
        ann_end = st.number_input("End (bp)", min_value=1, value=100, step=1, key="ann_end")
    ann_color = st.color_picker("Color", value="#FFA500", key="ann_color")

    if st.button("Add annotation", key="add_ann_btn") and ann_end > ann_start:
        ann = {"name": ann_name.strip() or "exon", "start": int(ann_start), "end": int(ann_end), "color": ann_color}
        if seq_choice == default_q_name:
            st.session_state.ann_q.append(ann)
        else:
            st.session_state.ann_r.append(ann)
        st.success("Annotation added.")

# ---------------- Edit / Delete Annotations ----------------
with st.expander("✏️ Edit Annotations", expanded=False):
    st.caption("Edit labels, positions, and colors; add or delete rows. Click **Save** to apply.")
    left, right = st.columns(2)

    with left:
        ann_q_df = st.data_editor(
            st.session_state.ann_q if st.session_state.ann_q else [{"name":"exon","start":0,"end":100,"color":"#FFA500"}],
            num_rows="dynamic",
            key="ann_q_editor",
            use_container_width=True,
            column_config={
                "name":  st.column_config.TextColumn("Label"),
                "start": st.column_config.NumberColumn("Start (bp)", min_value=0, step=1),
                "end":   st.column_config.NumberColumn("End (bp)",   min_value=1, step=1),
                "color": COLOR_COL_Q,
            },
        )
        if st.button("Save Query Annotations", key="save_q_anns"):
            st.session_state.ann_q = [
                {
                    "name":  str(row.get("name", "exon")).strip(),
                    "start": int(row.get("start", 0)),
                    "end":   int(row.get("end", 0)),
                    "color": str(row.get("color", "#FFA500")).strip() or "#FFA500",
                }
                for row in ann_q_df
                if str(row.get("name", "")).strip() != ""
            ]
            st.success("Saved query annotations.")

    with right:
        ann_r_df = st.data_editor(
            st.session_state.ann_r if st.session_state.ann_r else [{"name":"exon","start":0,"end":100,"color":"#87CEEB"}],
            num_rows="dynamic",
            key="ann_r_editor",
            use_container_width=True,
            column_config={
                "name":  st.column_config.TextColumn("Label"),
                "start": st.column_config.NumberColumn("Start (bp)", min_value=0, step=1),
                "end":   st.column_config.NumberColumn("End (bp)",   min_value=1, step=1),
                "color": COLOR_COL_R,
            },
        )
        if st.button("Save Reference Annotations", key="save_r_anns"):
            st.session_state.ann_r = [
                {
                    "name":  str(row.get("name", "exon")).strip(),
                    "start": int(row.get("start", 0)),
                    "end":   int(row.get("end", 0)),
                    "color": str(row.get("color", "#87CEEB")).strip() or "#87CEEB",
                }
                for row in ann_r_df
                if str(row.get("name", "")).strip() != ""
            ]
            st.success("Saved reference annotations.")

# ---------------- CSV Import / Export & Figure Download ----------------
with st.expander("📥 Upload / 📤 Download", expanded=False):
    st.markdown("**Annotations CSV** — [See CSV format guide](#annotation-csv-format)")

    colA, colB = st.columns(2)
    with colA:
        # Download annotations as CSV (sequence,name,start,end,color)
        all_rows = (
            [{"sequence":"query", **a} for a in st.session_state.ann_q] +
            [{"sequence":"ref",   **a} for a in st.session_state.ann_r]
        )
        csv_df = pd.DataFrame(all_rows, columns=["sequence","name","start","end","color"]) if all_rows else pd.DataFrame(columns=["sequence","name","start","end","color"]) 
        csv_bytes = csv_df.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download Annotations CSV",
            data=csv_bytes,
            file_name="annotations.csv",
            mime="text/csv",
            key="dl_ann_csv",
            use_container_width=True,
        )
    with colB:
        # Upload annotations CSV
        up = st.file_uploader("Upload Annotations CSV", type=["csv"], key="ann_csv_upload")
        if up is not None:
            try:
                df = pd.read_csv(up)
                required = {"sequence","name","start","end","color"}
                if not required.issubset(set(df.columns.str.lower())):
                    st.error("CSV missing required columns: sequence, name, start, end, color")
                else:
                    # normalize column names (lowercase)
                    df.columns = [c.lower() for c in df.columns]
                    q_rows = df[df["sequence"].str.lower().isin(["q","query"])].to_dict("records")
                    r_rows = df[df["sequence"].str.lower().isin(["r","ref","reference"])].to_dict("records")
                    st.session_state.ann_q = [
                        {"name": str(r.get("name","exon")).strip(), "start": int(r.get("start",0)), "end": int(r.get("end",0)), "color": str(r.get("color","#FFA500")).strip() or "#FFA500"}
                        for r in q_rows
                        if str(r.get("name",""))
                    ]
                    st.session_state.ann_r = [
                        {"name": str(r.get("name","exon")).strip(), "start": int(r.get("start",0)), "end": int(r.get("end",0)), "color": str(r.get("color","#87CEEB")).strip() or "#87CEEB"}
                        for r in r_rows
                        if str(r.get("name",""))
                    ]
                    st.success("Annotations loaded from CSV.")
            except Exception as e:
                st.error(f"Failed to parse CSV: {e}")

# Small in-page format guide anchor
st.markdown("""
#### Annotation CSV format
- Columns (**required**): `sequence`, `name`, `start`, `end`, `color`  
- `sequence`: one of `query`/`q` or `ref`/`r`/`reference`  
- `start`/`end`: integer basepair coordinates (inclusive start, exclusive end recommended)  
- `color`: hex like `#FFA500`  
- Example:
```
sequence,name,start,end,color
query,exon1,1200,1600,#FFA500
ref,dupA,7500,8200,#87CEEB
```
""", help=None)

# ---------------- Ready & Naming ----------------
ready = (
    (method == "Upload" and query_file is not None and ref_file is not None)
    or
    (method == "Paste FASTA text" and query_text.strip() and ref_text.strip())
)

ready = (
    (method == "Upload" and query_file is not None and ref_file is not None)
    or
    (method == "Paste FASTA text" and query_text.strip() and ref_text.strip())
)

if ready:
    # Default names based on inputs
    default_q_name = Path(query_file.name).stem if (method == "Upload" and query_file is not None) else "Query"
    default_r_name = Path(ref_file.name).stem if (method == "Upload" and ref_file is not None) else "Reference"
    q_name = st.text_input("Query name", value=default_q_name, key="q_name")
    r_name = st.text_input("Reference name", value=default_r_name, key="r_name")

    run_now = st.button("Run Alignment / Update Plot", key="run_btn")
    if run_now:
        # Grab sequences
        if method == "Upload":
            query_seq = query_file.getvalue().decode("utf-8")
            ref_seq   = ref_file.getvalue().decode("utf-8")
        else:
            query_seq = query_text
            ref_seq   = ref_text
        # Write to temp paths used by pipeline
        with open("/tmp/query.fasta", "w") as f: f.write(query_seq)
        with open("/tmp/ref.fasta", "w") as f: f.write(ref_seq)

        # Build objects with chosen gradient
        query_obj = AdhSeq("/tmp/query.fasta", GRADIENT, "query")
        ref_obj   = AdhSeq("/tmp/ref.fasta",   GRADIENT, "ref")
        # Override display names
        query_obj.name = q_name
        ref_obj.name   = r_name

        # ---------------- Use ONLY user annotations ----------------
        q_len = len(query_obj.getSeq())
        r_len = len(ref_obj.getSeq())
        add_q_ex = []
        for a in st.session_state.ann_q:
            s = max(0, int(a['start']))
            e = max(s + 1, int(a['end']))
            s = min(s, max(0, q_len - 1))
            e = min(e, q_len)
            if e > s:
                add_q_ex.append((s, e, a['name'], a['color']))
        add_r_ex = []
        for a in st.session_state.ann_r:
            s = max(0, int(a['start']))
            e = max(s + 1, int(a['end']))
            s = min(s, max(0, r_len - 1))
            e = min(e, r_len)
            if e > s:
                add_r_ex.append((s, e, a['name'], a['color']))
        query_obj.ex = add_q_ex
        ref_obj.ex   = add_r_ex

        # Generate figure
        with st.spinner("Running LASTZ and generating plot..."):
            fig = getGraph(query_obj, ref_obj, group_threshold=group_threshold, window_size=window_size)
            fig.set_dpi(fig_dpi)

        st.subheader("Alignment Result")
        if use_png_buffer:
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=fig_dpi, bbox_inches="tight")
            st.image(buf, use_container_width=True)
        else:
            # Download figure as PNG
            dl_buf = io.BytesIO()
            fig.savefig(dl_buf, format="png", dpi=fig_dpi, bbox_inches="tight")
            st.download_button("Download Figure (PNG)", data=dl_buf.getvalue(), file_name="alignment.png", mime="image/png", key="dl_fig_png")

            st.pyplot(fig, dpi=fig_dpi, use_container_width=True)
else:
    st.info("Choose an input method and provide both sequences.")
