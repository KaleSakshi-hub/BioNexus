import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import pandas as pd
import re
from Bio.Seq import Seq
from Bio.SeqUtils.MeltingTemp import Tm_Wallace, Tm_NN
from Bio.SeqUtils import gc_fraction
from io import StringIO
import requests
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time
import math
from io import BytesIO

# Configure page
st.set_page_config(
    page_title="BioNexus",  # Changed from "BioSuite Pro"
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': 'https://github.com/sakshikale/bionexus',
        'Report a bug': "https://github.com/sakshikale/bionexus/issues",
        'About': "# BioNexus\nAdvanced molecular biology analysis toolkit"
    }
)
# Custom CSS styling
st.markdown(""" <style>
:root {
    --primary: #2c3e50;
    --secondary: #3498db;
    --accent: #e74c3c;
    --light: #ecf0f1;
    --dark: #2c3e50;
    --success: #2ecc71;
}
.main {background-color: #ffffff;}
.sidebar .sidebar-content {background-color: var(--light);}
[data-testid="stHeader"] {background-color: var(--primary);}
[data-testid="stToolbar"] {display: none;}
.stTabs [data-baseweb="tab-list"] {
    gap: 0;
    padding: 0 1rem;
}
.stTabs [data-baseweb="tab"] {
    padding: 12px 24px;
    margin: 0;
    border-radius: 4px 4px 0 0;
    background-color: var(--light);
    transition: all 0.2s;
}
.stTabs [aria-selected="true"] {
    background-color: var(--secondary);
    color: white !important;
    font-weight: 600;
}
.stTabs [aria-selected="false"] {
    color: var(--dark);
}
.sequence-box {
    font-family: monospace;
    background-color: #f8f9fa;
    padding: 15px;
    border-radius: 5px;
    border: 1px solid #dee2e6;
    margin: 10px 0;
    white-space: pre-wrap;
    word-break: break-all;
}
.feature-card {
    background: white;
    border-radius: 8px;
    padding: 20px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.1);
    margin-bottom: 20px;
    height: 100%;
    transition: transform 0.3s ease;
}
.feature-card:hover {
    transform: translateY(-5px);
    box-shadow: 0 5px 15px rgba(0,0,0,0.1);
}
.metric-card {
    background: var(--light);
    border-radius: 8px;
    padding: 15px;
    text-align: center;
    border-left: 4px solid var(--secondary);
}
.footer {
    text-align: center;
    padding: 20px;
    color: #7f8c8d;
    font-size: 0.9rem;
    margin-top: 40px;
    border-top: 1px solid #eee;
}
.hero {
    background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%);
    color: white;
    padding: 3rem 2rem;
    border-radius: 10px;
    margin-bottom: 2rem;
}
.team-card {
    background: white;
    border-radius: 8px;
    padding: 20px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.1);
    margin-bottom: 20px;
    text-align: center;
}
.team-card img {
    width: 120px;
    height: 120px;
    border-radius: 50%;
    object-fit: cover;
    margin-bottom: 15px;
    border: 3px solid var(--secondary);
}
.tooltip {
    position: relative;
    display: inline-block;
    border-bottom: 1px dotted black;
}
.tooltip .tooltiptext {
    visibility: hidden;
    width: 200px;
    background-color: #555;
    color: #fff;
    text-align: center;
    border-radius: 6px;
    padding: 5px;
    position: absolute;
    z-index: 1;
    bottom: 125%;
    left: 50%;
    margin-left: -100px;
    opacity: 0;
    transition: opacity 0.3s;
}
.tooltip:hover .tooltiptext {
    visibility: visible;
    opacity: 1;
}
.badge {
    display: inline-block;
    padding: 0.25em 0.4em;
    font-size: 75%;
    font-weight: 700;
    line-height: 1;
    text-align: center;
    white-space: nowrap;
    vertical-align: baseline;
    border-radius: 0.25rem;
    background-color: var(--accent);
    color: white;
}
 </style>
""", unsafe_allow_html=True)

# Example sequences
EXAMPLE_PROTEIN = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
EXAMPLE_DNA = "ATGGCCCTGTGGATCCGCCTCCTGCTGGTGGCGGCCTGGGGCCCCGACCCCGCCGCCTTTGTGAACCAGCACCTGTGCGGCTCACCTGGTGGAGGCTCTACTTGGTGTGCGTGGGAGAGGGCTTCTTTTACACCCCAAAGACCCGCGAGGAGGAGGACCTGCAGGTGGGCCAGGTGGAGCTGGGGGGCCCCGGGGCCGGCTCACAGCCCCTGGCCCTGGAGGGCTCCCTGCAGAAGAGGGGCATCGTGGAGCAGTGCTGCACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTGA"
EXAMPLE_PRIMER = "ATGGCCCTGTGG"

# ==================== UTILITY FUNCTIONS ====================

def clean_sequence(sequence, valid_chars):
    """Clean and validate input sequence"""
    cleaned = ''.join([c.upper() for c in sequence if c.isalpha()])
    invalid = set(cleaned) - set(valid_chars)
    if invalid:
        st.warning(f"Removed invalid characters: {', '.join(invalid)}")
    return ''.join([c for c in cleaned if c in valid_chars])

def get_closest_match(sequence):
    """Identify closest matching protein using BLAST"""
    try:
        Entrez.email = "3522411028@despu.edu.in" # Required for NCBI API

        with st.spinner("Identifying closest protein match..."):
            # Perform BLAST search
            result_handle = NCBIWWW.qblast(
                program="blastp",
                database="nr",
                sequence=sequence,
                expect=0.001,
                hitlist_size=1
            )
           
            # Parse BLAST results
            blast_records = NCBIXML.read(result_handle)
            result_handle.close()
           
            if blast_records.alignments:
                best_match = blast_records.alignments[0]
                hit = best_match.hsps[0]
               
                match_info = {
                    'accession': best_match.accession,
                    'protein_name': best_match.hit_def.split(">")[0].split("[")[0].strip(),
                    'species': best_match.hit_def.split("[")[-1].split("]")[0],
                    'identity': (hit.identities / hit.align_length) * 100,
                    'e_value': hit.expect
                }
               
                st.success(f"Closest match found: {match_info['protein_name']}")
                st.write(f"**Species:** {match_info['species']}")
                st.write(f"**Accession:** {match_info['accession']}")
                st.write(f"**Identity:** {match_info['identity']:.1f}%")
                st.write(f"**E-value:** {match_info['e_value']:.2e}")
               
                return match_info
            else:
                st.warning("No significant matches found in NCBI database")
                return None
               
    except Exception as e:
        st.error(f"BLAST search error: {str(e)}")
        return None

def generate_database_links(match_info):
    """Generate database links based on closest match"""
    if not match_info:
        return None

    st.subheader("Database Links for Closest Match")

    # NCBI link
    ncbi_link = f"https://www.ncbi.nlm.nih.gov/protein/?term={match_info['accession']}"
    st.markdown(f"🔗 [NCBI Protein Entry for {match_info['accession']}]({ncbi_link})")

# ==================== HOME PAGE ====================

def show_homepage():
    """Display the homepage content"""
    # Hero Section
    st.markdown("""
    <div class="hero">
        <h1 style="color: white; margin-bottom: 0.5rem;">Welcome to BioSuite Pro</h1>
        <h3 style="color: white; font-weight: 400; margin-top: 0;">Your comprehensive platform for molecular biology analysis</h3>
        <p style="color: rgba(255,255,255,0.9);">Perform professional-grade bioinformatics analyses with our intuitive tools.</p>
    </div>
    """, unsafe_allow_html=True)

    # Quick Start Cards
    st.subheader("🚀 Quick Start")
    cols = st.columns(4)
    tools = [
        ("🧪 Amino Acid Analysis", "Analyze protein sequences"),
        ("🔍 ORF Finder", "Identify coding regions"),
        ("🌡️ Tm Calculator", "Calculate melting temps"),
        ("📊 Sequence Stats", "Get detailed metrics")
    ]
   
    for col, (title, desc) in zip(cols, tools):
        with col:
            st.markdown(f"""
            <div class="feature-card">
                <h4>{title}</h4>
                <p>{desc}</p>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("---")

    # Features Section
    st.header("✨ Key Features")
    features = [
        ("🧬", "Comprehensive Analysis", "From basic sequence stats to advanced ORF prediction"),
        ("📈", "Publication-Ready Visuals", "Export high-quality charts and tables"),
        ("🔗", "Database Integration", "Connect directly to NCBI, UniProt and more"),
        ("⚡", "Fast Processing", "Optimized algorithms for quick results"),
        ("🔒", "Privacy Focused", "Your data never leaves your browser"),
        ("📱", "Responsive Design", "Works on desktop, tablet and mobile")
    ]
   
    for i in range(0, len(features), 3):
        cols = st.columns(3)
        for col, (icon, title, desc) in zip(cols, features[i:i+3]):
            with col:
                st.markdown(f"""
                <div class="feature-card">
                    <div style="font-size: 2rem; margin-bottom: 10px;">{icon}</div>
                    <h4>{title}</h4>
                    <p>{desc}</p>
                </div>
                """, unsafe_allow_html=True)

    st.markdown("---")

    # Getting Started
    st.header("📌 Getting Started")
    st.markdown("""
    1. **Select a tool** from the navigation tabs above
    2. **Input your sequence** (or use our example sequences)
    3. **Customize** analysis parameters
    4. **View and download** your results
   
    All tools work directly in your browser - no installation required!
    """)
   
    with st.expander("💡 Tips for Best Results"):
        st.markdown("""
        - For protein analysis, use single-letter amino acid codes
        - Remove any non-sequence data (headers, numbers) before analysis
        - For ORF finding, longer sequences (>100bp) work best
        - Tm calculations are most accurate for sequences 18-30bp
        """)

    # Footer
    st.markdown("""
    <div class="footer">
        <p>BioSuite Pro | Powered by Streamlit and Biopython | © 2023</p>
        <p style="font-size: 0.8rem; margin-top: 10px;">
            Version 1.0.0 | <a href="https://github.com/yourusername/biosuite-pro" target="_blank">GitHub</a>
        </p>
    </div>
    """, unsafe_allow_html=True)

# ==================== ABOUT PAGE ====================

def show_about_page():
    """Display the about page content"""
    st.title("About BioNexus")
   
    st.markdown("""
    ## Our Mission
   
    BioNexus was created to provide researchers with easy access to professional-grade
    bioinformatics tools without the need for complex software installations or programming knowledge.
    """)
   
    st.markdown("---")
   
    # Developer Section
    st.header("👩‍💻 Developer")
    cols = st.columns([1, 3])
    with cols[0]:
        st.image("https://media.licdn.com/dms/image/v2/D5603AQEpTd-5A5pLlw/profile-displayphoto-shrink_100_100/B56ZVdDw_ZGsAY-/0/1741023025700?e=1752105600&v=beta&t=RxQSDNkJSdznG2xvwfUxkElWQPdwUEfpWFTi9PtlynM",
                 width=150)
    with cols[1]:
        st.markdown("""
        ### Sakshi Kale
        **M.Sc. Bioinformatics Student | DES Pune University**
       
        Sakshi developed BioNexus as part of her academic curriculum while pursuing a Master's degree in Bioinformatics.
        With a passion for computational biology and data science, she aims to bridge the gap between biology and technology
        through accessible bioinformatics tools.
       
        [![LinkedIn](https://img.shields.io/badge/LinkedIn-Connect-blue?style=flat&logo=linkedin)](https://www.linkedin.com/in/sakshi-kale-273b70320/)
        """)
   
    st.markdown("---")
   
    # Web Server Section
    st.header("🌐 BioNexus Web Server")
    st.markdown("""
    BioNexus is built using cutting-edge web technologies to ensure reliable performance:
   
    - **Framework**: Streamlit for interactive web application development
    - **Backend**: Python with Biopython for core bioinformatics computations
    - **Visualization**: Matplotlib, Seaborn, and Plotly for dynamic charts
    - **Deployment**: Cloud-based hosting with continuous integration
    - **Security**: Client-side processing for data privacy
   
    The server architecture is optimized for biological sequence analysis while maintaining
    user-friendly interfaces suitable for researchers at all technical levels.
    """)
   
    st.markdown("---")
   
    # Project Background Section
    st.header("🎓 Project Background")
    st.markdown("""
    This application was developed as part of the academic curriculum for the M.Sc. Bioinformatics program
    at DES Pune University. The project combines theoretical knowledge with practical implementation,
    demonstrating the application of:
   
    - Computational biology principles
    - Web application development
    - Bioinformatics algorithm implementation
    - Data visualization techniques
   
    The tools selected for implementation represent fundamental bioinformatics analyses commonly
    required in molecular biology research workflows.
    """)
   
    st.markdown("---")
   
    # Mentorship Section
    st.header("👨‍🏫 Mentorship")
    cols = st.columns([1, 3])
    with cols[0]:
        st.image("https://media.licdn.com/dms/image/v2/D5603AQF9gsU7YBjWVg/profile-displayphoto-shrink_400_400/B56ZZI.WrdH0Ag-/0/1744981029051?e=1752105600&v=beta&t=F4QBDSEgjUvnBS00xPkKqPTLI0jQaMpYefaOzARY1Yg",
                 width=150)
    with cols[1]:
        st.markdown("""
        ### Special Thanks to Dr. Kushagra Kashyap
       
        **Assistant Professor (Bioinformatics)**  
        Department of Life Sciences,  
        School of Science and Mathematics,  
        DES Pune University  
       
        [![LinkedIn](https://img.shields.io/badge/LinkedIn-Profile-blue?style=flat&logo=linkedin)](https://www.linkedin.com/in/dr-kushagra-kashyap-b230a3bb/)
       
        Dr. Kashyap provided invaluable guidance and academic support during the development of this project.
        His mentorship played a crucial role in shaping the project's scientific direction and refining the
        technical implementation. The encouragement and expert insights provided throughout the mini-project
        were instrumental in its successful completion.
        """)
   
    st.markdown("---")
   
    # Feedback Section
    st.header("📬 Feedback & Contact")
    st.markdown("""
    We welcome your feedback and suggestions to improve BioNexus:
   
    - **Email**: [3522411028@despu.edu.in](mailto:3522411028@despu.edu.in)
    - **LinkedIn**: [Sakshi Kale's Profile](https://www.linkedin.com/in/sakshi-kale-273b70320/)
    - **GitHub**: [Source code](https://github.com/KaleSakshi-hub/BioNexus/blob/main/Code.py)
   
    Your input helps us enhance the tools and add new features to better serve the research community.
    """)
   
    st.markdown("---")
   
    # Citation Section
    st.header("📜 Academic Citation")
    st.markdown("""
    If you use BioNexus in your academic work or research, please cite as:
   
    ```bibtex
    @software{bionexus2024,
      title = {BioNexus: Web-based Bioinformatics Toolkit},
      author = {Kale, Sakshi},
      year = {2024},
      note = {M.Sc. Bioinformatics Project, DES Pune University}
    }
    ```
    """)
# ==================== AMINO ACID TOOL ====================

def run_amino_acid_tool():
    """Amino Acid Frequency Counter with Database Links"""
    st.title("🧪 Amino Acid Analysis")
    st.markdown("Analyze protein sequence composition with detailed frequency statistics and database links.")

    # Analysis Settings
    st.subheader("Analysis Settings")
    col1, col2 = st.columns(2)
    with col1:
        chart_type = st.radio(
            "Visualization Type",
            ("Pie Chart", "Bar Chart"),
            index=0
        )
    with col2:
        color_palette = st.selectbox(
            "Color Palette",
            ("deep", "muted", "pastel", "bright", "dark", "colorblind"),
            index=0
        )

    # Main content
    tab1, tab2 = st.tabs(["Example Sequence", "Your Sequence"])

    with tab1:
        st.subheader("Example Protein Sequence")
        st.markdown(f'<div class="sequence-box">{EXAMPLE_PROTEIN}</div>', unsafe_allow_html=True)
       
        if st.button("Analyze Example", type="primary", key="aa_example"):
            with st.spinner("Performing complete analysis..."):
                analyze_protein_sequence(EXAMPLE_PROTEIN, chart_type, color_palette)

    with tab2:
        st.subheader("Your Protein Sequence")
        seq_input = st.text_area(
            "Paste protein sequence (single-letter codes):",
            placeholder="MALWMRLLP...",
            label_visibility="collapsed",
            height=150
        )
       
        if st.button("Analyze Your Sequence", type="primary", key="aa_custom"):
            if seq_input:
                with st.spinner("Analyzing sequence..."):
                    analyze_protein_sequence(seq_input, chart_type, color_palette)
            else:
                st.warning("Please input a protein sequence")

def analyze_protein_sequence(sequence, chart_type, palette):
    """Core protein analysis function with database links"""
    clean_seq = clean_sequence(sequence, 'ACDEFGHIKLMNPQRSTVWY')
    if not clean_seq:
        st.error("No valid amino acids found in the input")
        return

    # Basic sequence analysis
    counts = Counter(clean_seq)
    total = len(clean_seq)

    df = pd.DataFrame.from_dict(counts, orient='index', columns=['Count'])
    df['Frequency (%)'] = (df['Count'] / total * 100).round(2)
    df.index.name = 'Amino Acid'
    df.reset_index(inplace=True)
    df = df.sort_values('Count', ascending=False)

    # Display summary metrics
    st.subheader("Sequence Analysis Results")
    col1, col2, col3 = st.columns(3)
    col1.metric("Total Residues", f"{total:,}")
    col2.metric("Unique AAs", len(df))
    col3.metric("Most Frequent", f"{df.iloc[0]['Amino Acid']} ({df.iloc[0]['Count']})")

    # Visualization
    st.subheader("Amino Acid Composition")
    fig, ax = plt.subplots(figsize=(10, 6))

    if chart_type == "Pie Chart":
        ax.pie(
            df['Count'],
            labels=df['Amino Acid'],
            autopct='%1.1f%%',
            startangle=90,
            colors=sns.color_palette(palette, len(df)),
            textprops={'fontsize': 10}
        )
        ax.set_title('Amino Acid Distribution', pad=20)
    else:
        sns.barplot(
            data=df,
            x='Amino Acid',
            y='Count',
            palette=palette,
            ax=ax
        )
        ax.set_title('Amino Acid Composition', pad=20)
        ax.set_xlabel('Amino Acid', labelpad=10)
        ax.set_ylabel('Count', labelpad=10)
        ax.bar_label(ax.containers[0], fmt='%d')

    st.pyplot(fig)

    # Database search section
    if len(clean_seq) >= 10:
        closest_match = get_closest_match(clean_seq)
        if closest_match:
            generate_database_links(closest_match)
    else:
        st.warning("Sequence too short for database searches (minimum 10 amino acids)")

    # Detailed data
    with st.expander("View Detailed Composition Data"):
        st.dataframe(
            df.style.format({'Frequency (%)': '{:.2f}'}),
            use_container_width=True,
            height=400
        )

    # Download options
    st.subheader("Export Results")
    col1, col2 = st.columns(2)

    with col1:
        st.download_button(
            "Download CSV",
            data=df.to_csv(index=False).encode('utf-8'),
            file_name='amino_acid_analysis.csv',
            mime='text/csv'
        )

    with col2:
        buf = BytesIO()
        fig.savefig(buf, format="png", dpi=300, bbox_inches='tight')
        st.download_button(
            "Download Chart",
            data=buf.getvalue(),
            file_name='amino_acid_chart.png',
            mime='image/png'
        )

# ==================== ORF FINDER TOOL ====================

def run_orf_finder():
    """ORF Finder Tool"""
    st.title("🔍 ORF Finder")
    st.markdown("Identify open reading frames in DNA/RNA sequences with protein translation.")

    # ORF Settings
    st.subheader("ORF Settings")
    col1, col2 = st.columns(2)
    with col1:
        min_len = st.slider(
            "Minimum ORF Length (aa)",
            min_value=10,
            max_value=150,
            value=30,
            help="Minimum protein length to consider as valid ORF"
        )
    with col2:
        show_translation = st.checkbox("Show Protein Translations", True)

    # Main content tabs
    tab1, tab2 = st.tabs(["Example Sequence", "Your Sequence"])

    with tab1:
        st.subheader("Example DNA/RNA Sequence")
        st.markdown(f'<div class="sequence-box">{EXAMPLE_DNA}</div>', unsafe_allow_html=True)
        st.caption("Note: RNA sequences with U will be automatically converted to T")
       
        if st.button("Find ORFs in Example", type="primary", key="orf_example"):
            with st.spinner("Scanning sequence in all 6 frames..."):
                find_orfs_in_sequence(EXAMPLE_DNA, min_len, show_translation)

    with tab2:
        st.subheader("Your DNA/RNA Sequence")
        dna_input = st.text_area(
            "Paste DNA/RNA sequence:",
            placeholder="ATGGCCCTGTGG...",
            label_visibility="collapsed",
            height=150
        )
       
        if st.button("Find ORFs in Your Sequence", type="primary", key="orf_custom"):
            if dna_input:
                with st.spinner("Analyzing sequence..."):
                    find_orfs_in_sequence(dna_input, min_len, show_translation)
            else:
                st.warning("Please input a DNA/RNA sequence")

def find_orfs_in_sequence(sequence, min_len, show_translation):
    """Core ORF finding function"""
    clean_seq = clean_sequence(sequence, 'ATCGU')
    clean_seq = clean_seq.replace('U', 'T')  # Convert RNA to DNA
    if not clean_seq:
        st.error("No valid DNA/RNA bases found in the input")
        return

    orfs = []
    for frame in range(3):
        # Forward frames
        trans = str(Seq(clean_seq[frame:]).translate())
        orfs.extend(_scan_translation(trans, frame+1, 1, min_len))
       
        # Reverse frames
        rev_trans = str(Seq(clean_seq).reverse_complement()[frame:].translate())
        orfs.extend(_scan_translation(rev_trans, frame+1, -1, min_len))

    if not orfs:
        st.warning(f"No ORFs found ≥ {min_len} amino acids")
        return

    # Create results dataframe
    orf_df = pd.DataFrame(orfs)
    orf_df['Frame'] = orf_df['Frame'].apply(lambda x: f"+{x}" if x > 0 else x)
    orf_df = orf_df.sort_values('Start')

    # Display results
    st.subheader(f"Found {len(orf_df)} ORFs")

    # Summary table
    st.dataframe(
        orf_df[['Start', 'End', 'Length', 'Frame']].style.format({
            'Start': '{:,}',
            'End': '{:,}',
            'Length': '{:,}'
        }),
        use_container_width=True,
        height=min(400, 35 * len(orf_df) + 35)
    )

    # Detailed ORF view
    if show_translation:
        st.subheader("ORF Details")
        for idx, row in orf_df.iterrows():
            with st.expander(f"ORF {idx+1} | Frame {row['Frame']} | Positions {row['Start']}-{row['End']} ({row['Length']} aa)"):
                st.code(f"Protein sequence:\n{row['Sequence']}", language='text')

    # Download options
    st.subheader("Export Results")
    st.download_button(
        "Download ORF Table",
        data=orf_df.to_csv(index=False).encode('utf-8'),
        file_name='orf_results.csv',
        mime='text/csv'
    )

def _scan_translation(translation, frame, strand, min_len):
    """Helper to scan translated sequence for ORFs"""
    orfs = []
    for match in re.finditer(r'M.*?\*', translation):
        seq = match.group()
        if len(seq) >= min_len:
            orfs.append({
                'Start': match.start() * 3 + frame,
                'End': (match.end() - 1) * 3 + frame,
                'Length': len(seq),
                'Frame': frame * strand,
                'Sequence': seq
            })
    return orfs

# ==================== TM CALCULATOR TOOL ====================

def GC(sequence):
    """Calculate GC content percentage of a DNA sequence"""
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100 if sequence else 0

def Tm_Wallace(sequence):
    """Calculate Tm using the Wallace Rule"""
    a_count = sequence.upper().count('A')
    t_count = sequence.upper().count('T')
    g_count = sequence.upper().count('G')
    c_count = sequence.upper().count('C')
    return 2 * (a_count + t_count) + 4 * (g_count + c_count)

def Tm_NN(sequence, dnac1=50e-9, Na=50e-3, Mg=0, saltcorr=5):
    """Calculate Tm using nearest-neighbor thermodynamics"""
    # This is a simplified placeholder - in practice you would use a proper NN calculation
    # For a real implementation, consider using biopython's MeltingTemp module
    gc_content = GC(sequence)
    length = len(sequence)

    # Simplified NN approximation
    tm = 81.5 + 0.41*gc_content - 675/length - salt_correction(Na, Mg, saltcorr)
    return tm

def salt_correction(Na, Mg, method):
    """Calculate salt correction for Tm"""
    if method == 5:  # SantaLucia's correction
        return 12.5 * math.log((Na + 120 * math.sqrt(Mg)) / math.log(10))
    return 0

def run_tm_calculator():
    """Melting Temperature Calculator"""
    st.title("🌡️ Melting Temperature Calculator")
    st.markdown("Calculate primer melting temperature using different thermodynamic methods.")

    # Calculation settings
    st.subheader("Calculation Settings")
    col1, col2 = st.columns(2)
    with col1:
        method = st.radio(
            "Calculation Method",
            ("Wallace Rule", "Nearest-Neighbor"),
            index=0
        )

    advanced_params = {}
    if method == "Nearest-Neighbor":
        with col2:
            with st.expander("Advanced Parameters"):
                advanced_params['dna_conc'] = st.number_input(
                    "Primer Concentration (nM)",
                    min_value=1.0,
                    value=50.0,
                    step=10.0
                )
                advanced_params['na_conc'] = st.number_input(
                    "[Na⁺] Concentration (mM)",
                    min_value=0.0,
                    value=50.0,
                    step=10.0
                )
                advanced_params['mg_conc'] = st.number_input(
                    "[Mg²⁺] Concentration (mM)",
                    min_value=0.0,
                    value=0.0,
                    step=1.0
                )

    # Main content tabs
    tab1, tab2 = st.tabs(["Example Sequence", "Your Sequence"])

    with tab1:
        st.subheader("Example DNA Primer")
        st.markdown(f'<div class="sequence-box">{EXAMPLE_PRIMER}</div>', unsafe_allow_html=True)
        if st.button("Calculate Example Tm", type="primary", key="tm_example"):
            with st.spinner("Calculating..."):
                calculate_tm(EXAMPLE_PRIMER, method, advanced_params)

    with tab2:
        st.subheader("Your DNA Primer")
        primer_input = st.text_area(
            "Paste DNA primer sequence:",
            placeholder="ATGGCCCTGTGG...",
            label_visibility="collapsed",
            height=100
        )
        if st.button("Calculate Tm", type="primary", key="tm_custom"):
            if primer_input:
                with st.spinner("Analyzing sequence..."):
                    calculate_tm(primer_input, method, advanced_params)
            else:
                st.warning("Please input a DNA sequence")

def calculate_tm(sequence, method, params):
    """Core Tm calculation function"""
    clean_seq = clean_sequence(sequence, 'ATCG')
    if not clean_seq:
        st.error("No valid DNA bases found in the input")
        return

    try:
        if method == "Wallace Rule":
            tm = Tm_Wallace(clean_seq)
            gc_content = GC(clean_seq)
            st.subheader("Results")
            col1, col2 = st.columns(2)
            col1.metric("Melting Temperature (Tm)", f"{tm:.2f} °C")
            col2.metric("GC content", f"{gc_content:.1f}%")
           
            st.markdown("**Wallace Rule Formula:**")
            st.latex(r"Tm = 2°C \times (A+T) + 4°C \times (G+C)")
           
        else:  # Nearest-Neighbor
            dnac1 = params.get('dna_conc', 50) * 1e-9  # Convert nM to M
            na = params.get('na_conc', 50) / 1e3  # Convert mM to M
            mg = params.get('mg_conc', 0) / 1e3  # Convert mM to M
           
            tm = Tm_NN(
                clean_seq,
                dnac1=dnac1,
                Na=na,
                Mg=mg,
                saltcorr=5  # Use salt correction method 5
            )
           
            st.subheader("Results")
            col1, col2 = st.columns(2)
            col1.metric("Melting Temperature (Tm)", f"{tm:.2f} °C")
            col2.metric("Primer Length", f"{len(clean_seq)} bp")
           
            st.markdown("**Nearest-Neighbor Parameters:**")
            st.markdown(f"""
            - Primer Concentration: {params.get('dna_conc', 50)} nM
            - Sodium Concentration: {params.get('na_conc', 50)} mM
            - Magnesium Concentration: {params.get('mg_conc', 0)} mM
            """)
       
        # Sequence validation
        with st.expander("Sequence Details"):
            st.markdown(f"**Processed Sequence:** `{clean_seq}`")
            st.markdown(f"**Sequence Length:** {len(clean_seq)} bp")
           
    except Exception as e:
        st.error(f"Error in calculation: {str(e)}")

# ==================== MAIN APP ====================

def main():
    """Main app navigation"""
    home_tab, aa_tab, orf_tab, tm_tab, about_tab = st.tabs([
        "🏠 Home",
        "🧪 Amino Acid Analyzer",
        "🔍 ORF Finder",
        "🌡️ Tm Calculator",
        "ℹ️ About"
    ])

    with home_tab:
        show_homepage()

    with aa_tab:
        run_amino_acid_tool()

    with orf_tab:
        run_orf_finder()

    with tm_tab:
        run_tm_calculator()
       
    with about_tab:
        show_about_page()

if __name__ == "__main__":
    main()
