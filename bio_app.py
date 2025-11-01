import streamlit as st
import os
import io
import ssl
import subprocess
from Bio import Entrez, SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight

# --- SSL Certificate Workaround ---
# This is a workaround for SSL certificate verification errors that can occur when
# fetching data from NCBI on some networks or systems. It tells Python not to
# verify the server's SSL certificate.
try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    # Python < 2.7.9 / 3.4 does not have this attribute
    pass
else:
    ssl._create_default_https_context = _create_unverified_https_context


# --- Core Functions ---

def fetch_genbank_record(accession_id):
    """Fetches a GenBank record using the provided accession ID."""
    if not accession_id:
        st.warning("Please enter a valid Accession ID.")
        return None
    try:
        st.info(f"Fetching record for {accession_id} from NCBI...")
        Entrez.email = "your.email@example.com"  # IMPORTANT: Always tell Entrez who you are
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        st.success(f"Successfully fetched record for {record.id} ({len(record.seq)} bp)")
        return record
    except Exception as e:
        st.error(f"An error occurred while fetching the record: {e}")
        return None


def parse_fasta_input(fasta_input):
    """Parses FASTA data from either an uploaded file or a text area."""
    if not fasta_input:
        return []
    try:
        # Use StringIO to treat the string data as a file
        stringio = io.StringIO(fasta_input)
        records = list(SeqIO.parse(stringio, "fasta"))
        if not records:
            st.warning("Could not parse any sequences. Please ensure the input is in valid FASTA format.")
            return []
        st.success(f"Successfully parsed {len(records)} sequence(s).")
        return records
    except Exception as e:
        st.error(f"An error occurred during FASTA parsing: {e}")
        return []


def display_features(record):
    """Displays the features of a GenBank record in a structured way."""
    st.subheader("Genomic Features")
    if not record.features:
        st.info("No features found in this record.")
        return

    features_data = []
    for feature in record.features:
        features_data.append({
            "Type": feature.type,
            "Location": str(feature.location),
            "Qualifier": ", ".join(f"{k}={v}" for k, v in feature.qualifiers.items())
        })
    st.dataframe(features_data)


def translate_dna_to_protein(sequence):
    """Translates a DNA sequence to a protein sequence."""
    st.subheader("DNA to Protein Translation")
    try:
        # Ensure the sequence length is a multiple of 3
        if len(sequence) % 3 != 0:
            st.warning("Sequence length is not a multiple of 3. Translation might be incomplete.")
            trimmed_len = len(sequence) - (len(sequence) % 3)
            sequence = sequence[:trimmed_len]

        protein_seq = sequence.translate()
        st.write("**DNA Sequence:**")
        st.code(str(sequence), language="text")
        st.write("**Translated Protein Sequence:**")
        st.code(str(protein_seq), language="text")
    except Exception as e:
        st.error(f"An error occurred during translation: {e}")


def find_orfs(sequence, min_len=50):
    """Finds Open Reading Frames (ORFs) in a DNA sequence."""
    orfs = []
    strands = {
        "Forward": sequence,
        "Reverse": sequence.reverse_complement()
    }

    for strand_name, seq in strands.items():
        for frame in range(3):
            # Trim sequence to be a multiple of 3
            trimmed_seq = seq[frame: len(seq) - (len(seq) - frame) % 3]
            protein_seq = trimmed_seq.translate()

            # Split by stop codon '*'
            for potential_orf in str(protein_seq).split('*'):
                # Find the first start codon 'M'
                start_index = potential_orf.find('M')
                if start_index != -1:
                    # The actual ORF starts from 'M'
                    orf = potential_orf[start_index:]
                    if len(orf) >= min_len:
                        orfs.append((orf, len(orf), strand_name, frame + 1))
    return orfs


def detailed_sequence_analysis(record):
    """Performs a detailed analysis of a sequence record."""
    st.subheader(f"Detailed Analysis of {record.id}")

    # --- Basic Stats ---
    st.markdown("---")
    st.write(f"**ID:** {record.id}")
    st.write(f"**Length:** {len(record.seq)} bp")
    st.write(f"**GC Content:** {gc_fraction(record.seq) * 100:.2f}%")
    st.write(f"**Molecular Weight (DNA):** {molecular_weight(record.seq, seq_type='DNA'):.2f}")
    st.markdown("---")

    # --- ORF Finding ---
    st.subheader(f"Annotating {record.id} by finding ORFs")
    orfs = find_orfs(record.seq)

    if orfs:
        st.write(f"Found {len(orfs)} potential ORFs (length >= 50 aa).")
        # Sort ORFs by length, descending
        orfs.sort(key=lambda x: x[1], reverse=True)
        for i, (orf_seq, orf_len, strand, frame) in enumerate(orfs):
            st.write(f"**ORF {i + 1} (Length: {orf_len} aa, Strand: {strand}, Frame: {frame})**")
            snippet = orf_seq[:60] + "..." if len(orf_seq) > 60 else orf_seq
            st.code(snippet, language="text")
    else:
        st.info("No potential ORFs found with a minimum length of 50 amino acids.")


def perform_pairwise_alignment(seq1_str, seq2_str, align_type):
    """Performs pairwise global or local alignment on two sequences."""
    st.subheader(f"Pairwise {align_type} Alignment")
    if not seq1_str or not seq2_str:
        st.warning("Please provide two sequences for alignment.")
        return

    seq1 = Seq(seq1_str)
    seq2 = Seq(seq2_str)
    aligner = Align.PairwiseAligner()

    if align_type == "Global":
        alignments = aligner.align(seq1, seq2)
        st.write(f"**Global Alignment Score:** {aligner.score(seq1, seq2)}")
    else:  # Local
        alignments = aligner.align(seq1, seq2)
        st.write(f"**Local Alignment Score:** {aligner.score(seq1, seq2)}")

    st.write("**Alignment:**")
    # Show the first and best alignment
    if alignments:
        st.code(str(alignments[0]), language="text")
    else:
        st.info("No alignment could be produced.")


def perform_msa(records):
    """Performs Multiple Sequence Alignment using ClustalW."""
    st.subheader("Multiple Sequence Alignment (MSA)")
    st.info("This feature requires ClustalW to be installed and accessible in your system's PATH.")

    # Create a temporary file to store the input sequences
    temp_fasta_file = "msa_input.fasta"
    with open(temp_fasta_file, "w") as f:
        SeqIO.write(records, f, "fasta")

    # Define output file
    output_aln_file = "msa_output.aln"

    try:
        # Check if clustalw executable is available
        result = subprocess.run(["clustalw", "-version"], check=True, capture_output=True, text=True)
        st.write(f"Using ClustalW: {result.stdout.strip()}")

        # Run ClustalW directly using subprocess
        st.write("Running ClustalW...")
        cmd = ["clustalw", f"-INFILE={temp_fasta_file}", f"-OUTFILE={output_aln_file}"]
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)

        # Display stdout if available
        if result.stdout:
            st.text(result.stdout)

        # Read and display the alignment
        if os.path.exists(output_aln_file):
            with open(output_aln_file, "r") as f:
                alignment_content = f.read()
            st.write("**MSA Result (Clustal Format):**")
            st.code(alignment_content, language="text")
        else:
            st.error("Alignment output file was not created.")

    except FileNotFoundError:
        st.error("ClustalW not found. Please make sure it is installed and in your PATH.")
        st.info(
            "On Windows, download ClustalW from http://www.clustal.org/download/current/ and add it to your system PATH.")
    except subprocess.CalledProcessError as e:
        st.error(f"ClustalW encountered an error: {e}")
        if e.stderr:
            st.text(f"Error output: {e.stderr}")
    except Exception as e:
        st.error(f"An error occurred during MSA: {e}")
    finally:
        # Clean up temporary files
        if os.path.exists(temp_fasta_file):
            os.remove(temp_fasta_file)
        if os.path.exists(output_aln_file):
            os.remove(output_aln_file)
        if os.path.exists("msa_input.dnd"):  # ClustalW guide tree file
            os.remove("msa_input.dnd")


# --- Streamlit App Layout ---

st.set_page_config(page_title="Bioinformatics Power Tool", layout="wide")
st.title("🧬 Bioinformatics Power Tool 🔬")
st.write("A web app to perform common bioinformatics tasks using Biopython.")

# --- Sidebar for Input Selection ---
st.sidebar.title("Input Options")
input_method = st.sidebar.radio("Choose your input method:", ("Enter Accession ID", "Upload or Paste FASTA"))

records = []

if input_method == "Enter Accession ID":
    accession_id = st.sidebar.text_input("GenBank Accession ID", placeholder="e.g., NC_000852")
    if st.sidebar.button("Fetch and Process"):
        record = fetch_genbank_record(accession_id)
        if record:
            records = [record]  # Store as a list for consistency
            st.session_state['records'] = records  # Use session state to store data

elif input_method == "Upload or Paste FASTA":
    uploaded_file = st.sidebar.file_uploader("Upload a FASTA file", type=["fasta", "fa", "fna"])
    pasted_sequence = st.sidebar.text_area("Or paste FASTA sequence(s) here", height=200)

    if st.sidebar.button("Process FASTA"):
        if uploaded_file:
            # To read file as string:
            string_data = uploaded_file.getvalue().decode("utf-8")
            records = parse_fasta_input(string_data)
        elif pasted_sequence:
            records = parse_fasta_input(pasted_sequence)
        else:
            st.sidebar.warning("Please upload a file or paste a sequence.")

        if records:
            st.session_state['records'] = records

# --- Main Panel for Operations and Output ---
if 'records' in st.session_state and st.session_state['records']:
    st.header("Choose Operation")

    # Use the first record for single-sequence operations
    main_record = st.session_state['records'][0]

    # Let user select which record to use if multiple are available
    if len(st.session_state['records']) > 1:
        record_ids = [rec.id for rec in st.session_state['records']]
        selected_id = st.selectbox("Select a sequence for single-record analysis:", record_ids)
        main_record = next((rec for rec in st.session_state['records'] if rec.id == selected_id), main_record)

    operation = st.selectbox(
        "Select an analysis to perform:",
        [
            "--- Select ---",
            "Analyze Features (GenBank input only)",
            "Detailed Sequence Analysis",
            "Translate DNA to Protein",
            "Sequence Alignment",
        ]
    )

    if operation == "Analyze Features (GenBank input only)":
        if main_record.features:
            display_features(main_record)
        else:
            st.warning(
                "Feature analysis requires a GenBank file with annotations. The provided FASTA file does not contain features.")

    elif operation == "Detailed Sequence Analysis":
        detailed_sequence_analysis(main_record)

    elif operation == "Translate DNA to Protein":
        translate_dna_to_protein(main_record.seq)

    elif operation == "Sequence Alignment":
        st.header("Alignment Options")
        alignment_mode = st.radio("Select alignment type:", ("Pairwise Alignment", "Multiple Sequence Alignment (MSA)"))

        if alignment_mode == "Pairwise Alignment":
            st.subheader("Input Sequences for Pairwise Alignment")
            seq1_default = str(main_record.seq)
            seq2_default = ""

            col1, col2 = st.columns(2)
            with col1:
                seq1_input = st.text_area("Sequence 1", value=seq1_default, height=200)
            with col2:
                seq2_input = st.text_area("Sequence 2", value=seq2_default, height=200,
                                          placeholder="Enter second sequence here...")

            pairwise_type = st.radio("Select Pairwise Alignment method:", ("Global", "Local"))

            if st.button(f"Run {pairwise_type} Alignment"):
                perform_pairwise_alignment(seq1_input, seq2_input, pairwise_type)

        elif alignment_mode == "Multiple Sequence Alignment (MSA)":
            if len(st.session_state['records']) > 1:
                st.write(f"**Found {len(st.session_state['records'])} sequences for MSA.**")
                if st.button("Run MSA with ClustalW"):
                    perform_msa(st.session_state['records'])
            else:
                st.warning(
                    "MSA requires multiple sequences. Please upload a FASTA file containing at least two sequences.")
