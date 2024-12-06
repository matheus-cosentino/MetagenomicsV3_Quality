#MsC. Matheus Cosentino
#analise de dados Illumina, com base no output do arquivo Slurm gerado
import re
import argparse
import csv

#Parse Arguments parser
parser = argparse.ArgumentParser(description="Process a Slurm output file of Metagenomics.")
parser.add_argument(
    "--slurm",
    type=str,
    required=True,
    help="Path to the Slurm output file (.out) generated in MetagenomicV3.\nExample: MetaGen3_173518.out"
    )
parser.add_argument(
    "--output",
    type=str,
    required=True,
    help="Path to the output file (.csv) with processed data by this script .\nExample: MetaGen3_173518.csv"
    )

args = parser.parse_args()

#Function to get Library Name
def get_names(section):
    # Extract and clean the sample name
    sample_name = section.split('\n')[0].strip()
    sample_name = sample_name.strip()  # Final cleanup
    sample_data[sample_name] = {}
    return sample_name

# Function to get the reads passed filter (Q30 value)
def get_Q30_R1(section):
    lines = section.split('\n')
    R1_Q30 = None
    R1_Q30_percent = None
    for i, line in enumerate(lines):
        if 'Read1 after filtering' in line:
            value=lines[i+1].split(':')[1].strip()
            R1_Q30 = value
            R1_Q30_percent = re.search(r'\((\d+\.\d+)%\)', lines[i+4]).group(1)
    if R1_Q30 is not None and R1_Q30_percent is not None:
        result = f"{R1_Q30} ({R1_Q30_percent}%)"
        return (result)
    else:
        return "No acha"
def get_Q30_R2(section):
    lines = section.split('\n')
    R2_Q30 = None
    R2_Q30_percent = None
    for i, line in enumerate(lines):
        if 'Read1 after filtering' in line:
            value=lines[i+1].split(':')[1].strip()
            R2_Q30 = value
            R2_Q30_percent = re.search(r'\((\d+\.\d+)%\)', lines[i+4]).group(1)
    if R2_Q30 is not None and R2_Q30_percent is not None:
        result = f"{R2_Q30} ({R2_Q30_percent}%)"
        return (result)
    else:
        return "No acha"

# Function to get the raw reads quantity
def get_raw_r1(section):
    lines = section.split('\n')
    R1_RAW = None
    R1_RAW_percent = None
    for i, line in enumerate(lines):
        if 'Read1 before filtering:' in line:
            value=lines[i+1].split(':')[1].strip()
            R1_RAW = value
            R1_RAW_percent = re.search(r'\((\d+\.\d+)%\)', lines[i+4]).group(1)
    if R1_RAW is not None and R1_RAW_percent is not None:
        result = f"{R1_RAW} ({R1_RAW_percent}%)"
        return (result)
    else:
        return "No acha"
def get_raw_r2(section):
    lines = section.split('\n')
    R2_RAW = None
    R2_RAW_percent = None
    for i, line in enumerate(lines):
        if 'Read1 before filtering:' in line:
            value=lines[i+1].split(':')[1].strip()
            R2_RAW = value
            R2_RAW_percent = re.search(r'\((\d+\.\d+)%\)', lines[i+4]).group(1)
    if R2_RAW is not None and R2_RAW_percent is not None:
        result = f"{R2_RAW} ({R2_RAW_percent}%)"
        return (result)
    else:
        return "No acha"

# Function to get the paired reads quality
def paired_quality(section):
    lines = section.split('\n')
    reads_passed_filter = None
    reads_failed_lq = None
    reads_failed_short = None
    reads_adapter = None
    duplication_rate = None
    insert_peak = None
    for i, line in enumerate(lines):
        if 'Filtering result:' in line:
            reads_passed_filter = lines[i+1].split(':')[1].strip()
            reads_failed_lq = lines[i+2].split(':')[1].strip()
            reads_failed_short = lines[i+4].split(':')[1].strip()
            reads_adapter = lines[i+5].split(':')[1].strip()
        if 'Duplication rate:' in line:
            duplication_rate = line.split(':')[1].strip()
            duplication_rate = f"{float(line.split(':')[1].strip()[:-1]):,.2f}%"
        if 'Insert size peak' in line:
            insert_peak = line.split(':')[-1].strip()
        if reads_passed_filter is not None and reads_failed_lq is not None and reads_failed_short is not None and reads_adapter is not None and duplication_rate is not None and insert_peak is not None:
            result = f"{reads_passed_filter}, {reads_failed_lq}, {reads_failed_short}, {reads_adapter}, {duplication_rate}, {insert_peak}" 
            return result

#Function to get reads mapped to host genome    
def get_mapped_paired(section):
    lines = section.split('\n')
    mapped_reads = None
    for i, line in enumerate(lines):
        if 'Creating MAPPED and UNMAPPED from PAIRED SORTED BAM files' in line:
            mapped_reads = lines[i+2].split()[-2].strip()
            unmapped_reads = lines[i+4].split()[-2].strip()
    if mapped_reads is not None and unmapped_reads is not None:
        result = f"{mapped_reads} {unmapped_reads}"
        return result
def get_mapped_unpaired(section):
    lines = section.split('\n')
    mapped_reads = None
    for i, line in enumerate(lines):
        if 'Creating MAPPED and UNMAPPED from UNPAIRED SORTED BAM files' in line:
            mapped_reads = lines[i+2].split()[-2].strip()
            unmapped_reads = lines[i+4].split()[-2].strip()
    if mapped_reads is not None and unmapped_reads is not None:
        result = f"{mapped_reads} {unmapped_reads}"
        return result
def get_mapped_contigs(section):
    lines = section.split('\n')
    mapped_contigs = None
    for i, line in enumerate(lines):
        if 'Creating MAPPED and UNMAPPED from CONTIGS SORTED BAM files' in line:
            mapped_contigs = lines[i+2].split()[-2].strip()
            unmapped_contigs = lines[i+4].split()[-2].strip()
    if mapped_contigs is not None and unmapped_contigs is not None:
        result = f"{mapped_contigs} {unmapped_contigs}"
        return result


#Function to get reads information taxonomically iddentified by Kraken
def get_kraken_reads(section):
    lines = section.split('\n')
    kraken_reads = None
    for i, line in enumerate(lines):
        if 'Performing Kraken2 taxonomic assignment (reads)' in line:
            try:
                kraken_reads = lines[i+4].split()[0].strip()
                kraken_classified = lines[i+5].split()[0].strip()
                kraken_classified_percent = lines[i+5].split()[-1].strip()
                kraken_unclassified = lines[i+6].split()[0].strip()
                kraken_unclassified_percent = lines[i+6].split()[-1].strip()
            except IndexError:
                print("Error: unexpected format in Kraken2 output")
                return None
    if kraken_reads is not None and kraken_classified is not None and kraken_classified_percent is not None and kraken_unclassified is not None and kraken_unclassified_percent is not None:
        result = f"{kraken_reads} {kraken_classified} {kraken_classified_percent} {kraken_unclassified} {kraken_unclassified_percent}"
        return result
    else:
        print("Error: unable to extract Kraken2 reads")
        return None

def get_kraken_contigs(section):
    lines = section.split('\n')
    kraken_contigs = None
    for i, line in enumerate(lines):
        if 'Performing Kraken2 taxonomic assignment (contigs)' in line:
            try:
                kraken_contigs = lines[i+4].split()[0].strip()
                kraken_classified = lines[i+5].split()[0].strip()
                kraken_classified_percent = lines[i+5].split()[-1].strip()
                kraken_unclassified = lines[i+6].split()[0].strip()
                kraken_unclassified_percent = lines[i+6].split()[-1].strip()
            except IndexError:
                print("Error: unexpected format in Kraken2 output")
                return None
    if kraken_contigs is not None and kraken_classified is not None and kraken_classified_percent is not None and kraken_unclassified is not None and kraken_unclassified_percent is not None:
        result = f"{kraken_contigs} {kraken_classified} {kraken_classified_percent} {kraken_unclassified} {kraken_unclassified_percent}"
        return result
    else:
        print("Error: unable to extract Kraken2 reads")
        return None

# Function to get reads information taxonomically iddentified by Diamond
def get_Diamond(section):
    lines = section.split('\n')
    diamond_reads = None
    diamond_contigs = None
    for i, line in enumerate(lines):
        if 'Creating DIAMOND diversity plot (reads) for library' in line:
            try:
                diamond_reads = lines[i-2].split()[0].strip()
            except IndexError:
                print("Error: unable to extract DIAMOND reads")
                return None
        if 'Creating DIAMOND diversity plot (contigs) for library' in line:
            try:
                diamond_contigs = lines[i-2].split()[0].strip()
            except IndexError:
                print("Error: unable to extract DIAMOND contigs")
                return None
    if diamond_reads is not None and diamond_contigs is not None:
        result = f"{diamond_reads} {diamond_contigs}"
        return result
    else:
        print("Error: unable to extract DIAMOND reads or contigs")
        return None

# Check if the file path was provided and process it
try:
    with open(args.slurm, 'r') as file:
        # Skip the first 60 lines
        for _ in range(60):
            next(file)
        
        # Read the rest of the file
        content = file.read()
        
        # Split the content into sections
        sections = content.split('Sample directory:')
        
        # Skip the first section
        sections = sections[1:]
except FileNotFoundError:
    print(f"Error: The file '{args.slurm}' was not found.")
    kill
except Exception as e:
    print(f"An error occurred: {e}")
    kill
    
# Check if the file path was provided and process it
try:
    with open(args.slurm, 'r') as file:
        # Skip the first 60 lines
        for _ in range(60):
            next(file)
        
        # Read the rest of the file
        content = file.read()
        
        # Split the content into sections
        sections = content.split('Sample directory:')
        
        # Skip the first section
        sections = sections[1:]
except FileNotFoundError:
    print(f"Error: The file '{args.slurm}' was not found.")
    kill
except Exception as e:
    print(f"An error occurred: {e}")
    kill
    
#Separet file per library
sample_data = {}

#process
for section in sections:
    # Get the sample name
    sample_name = get_names(section)

    r1_raw = get_raw_r1(section)
    r2_raw = get_raw_r2(section)
    r1_passed = get_Q30_R1(section)
    r2_passed = get_Q30_R2(section)
    paired_quality_reads = paired_quality(section)
    mapped_paired = get_mapped_paired(section)
    mapped_unpaired = get_mapped_unpaired(section)
    mapped_contigs = get_mapped_contigs(section)
    kraken_reads = get_kraken_reads(section)
    kraken_contigs = get_kraken_contigs(section)
    diamond = get_Diamond(section)

    if len(paired_quality_reads.split(',')) != 6:
        print(f"Error: unexpected format for paired quality reads in sample {sample_name}")
        continue

    if len(mapped_paired.split()) != 2:
        print(f"Error: unexpected format for mapped reads in sample {sample_name}")
        continue

    sample_data[sample_name]['r1_raw'] = r1_raw
    sample_data[sample_name]['r2_raw'] = r2_raw 
    sample_data[sample_name]['r1_passed'] = r1_passed
    sample_data[sample_name]['r2_passed'] = r2_passed
    sample_data[sample_name]['reads_passed_filter'] = paired_quality_reads.split(',')[0]
    sample_data[sample_name]['reads_failed_lq'] = paired_quality_reads.split(',')[1]
    sample_data[sample_name]['reads_failed_short'] = paired_quality_reads.split(',')[2]
    sample_data[sample_name]['reads_adapter'] = paired_quality_reads.split(',')[3]
    sample_data[sample_name]['duplication_rate'] = paired_quality_reads.split(',')[4]
    sample_data[sample_name]['insert_peak'] = paired_quality_reads.split(',')[5]
    sample_data[sample_name]['paired_mapped_reads'] = mapped_paired.split()[0]
    sample_data[sample_name]['paired_unmapped_reads'] = mapped_paired.split()[1]
    sample_data[sample_name]['unpaired_mapped_reads'] = mapped_unpaired.split()[0]
    sample_data[sample_name]['unpaired_unmapped_reads'] = mapped_unpaired.split()[1]
    sample_data[sample_name]['contigs_mapped'] = mapped_contigs.split()[0]
    sample_data[sample_name]['contigs_unmapped'] = mapped_contigs.split()[1]
    sample_data[sample_name]['kraken_reads_processed'] = kraken_reads.split()[0]
    sample_data[sample_name]['kraken_reads_classified'] = f"{kraken_reads.split()[1]} {kraken_reads.split()[2]}"
    sample_data[sample_name]['kraken_reads_unclassified'] = f"{kraken_reads.split()[3]} {kraken_reads.split()[4]}"
    sample_data[sample_name]['kraken_contigs_processed'] = kraken_contigs.split()[0]
    sample_data[sample_name]['kraken_contigs_classified'] = f"{kraken_contigs.split()[1]} {kraken_contigs.split()[2]}"
    sample_data[sample_name]['kraken_contigs_unclassified'] = f"{kraken_contigs.split()[3]} {kraken_contigs.split()[4]}"
    sample_data[sample_name]['diamond_reads'] = diamond.split()[0]
    sample_data[sample_name]['diamond_contigs'] = diamond.split()[1]


# Open the output file in write mode
with open(args.output, 'w', newline='') as f:
    # Create a CSV writer object
    writer = csv.writer(f, delimiter=',')

    # Write the header row
    header_row = ['Sample Name'] + list(sample_data[sample_name].keys())
    writer.writerow(header_row)

    # Write the data rows
    for sample_name, sample_data in sample_data.items():
        row = [sample_name] + list(sample_data.values())
        writer.writerow(row)