import os
import argparse
import subprocess
import pandas as pd
import datetime
import glob
import multiprocessing

################################################################################
#################           FUNCTIONS           ################################
################################################################################

################################################################################
# 1. READS TRIMMING
def ListReadsSamples(folder):
  try:
    # Check if the path is a directory
    if not os.path.isdir(folder):
      print(f"The path '{folder}' is not a directory.")
      return

    # List all entries in the directory
    samples = glob.glob(os.path.join(folder, '*/*R1_001*.f*q.gz'))
    samples = [sample for sample in samples if 'control' not in sample]
    samples = [sample for sample in samples if 'Nc' not in sample]
    return(samples)

  except Exception as e:
    print(f"An error occurred: {e}")
    return([])

def RunTrimGalore(out_folder, reads_forward, reads_reverse, cpus);
  args = f'trim_galore --fastqc --paired -j {str(cpus)} -o {os.path.join(out_folder, 'trimmed_reads')} {reads_forward} {reads_reverse}'
  subprocess.call(args, shell = True)
  
  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(' '.join(args) + '\n\n')
    
def ProcessTrimGalore(out_folder, short_reads_samples, cpus):
  for reads_forward in short_reads_samples:
    reads_reverse = reads_forward.replace('R1_001.fastq.gz', 'R2_001.fastq.gz')
    RunTrimGalore(out_folder, reads_forward, reads_reverse, cpus)

def RunPorechop(out_folder, reads_in, reads_out, cpus):
  args = f'porechop --threads {str(cpus)} -i {reads_in} -o {reads_out}'
  subprocess.call(args, shell = True)

  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(args + '\n\n')

def RunChopper(out_folder, reads_in, reads_out, cpus):
  args = f'gunzip -c {reads_in} | chopper -q {str(12)} --threads {str(cpus)} | gzip > {reads_out}'
  subprocess.call(args, shell = True)

  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(args + '\n\n')
    
def ProcessPorechopChopper(out_folder, long_reads_samples, cpus):
  for reads_in in long_reads_samples:
    reads_porechop = os.path.join(out_folder, 'trimmed_reads', 'long_reads', reads_in.split('/')[-1].replace('.fastq.gz', '_porechopped.fastq.gz'))
    reads_chopper = os.path.join(out_folder, 'trimmed_reads', 'long_reads', reads_in.split('/')[-1].replace('.fastq.gz', '_chopped.fastq.gz'))
    RunPorechop(out_folder, reads_in, reads_porechop, cpus)
    RunChopper(out_folder, reads_porechop, reads_chopper, cpus)
    
################################################################################
# 2. ASSEMBLIES (one co-assembly, several smaller co-assemblies, individual assemblies)
def LoadMetadata():

def EvaluateSamples(type, metadata):
  if type == 'co':
    
  elif type == 'sub:
    for samples, i in enumerate(list_of_list_samples):
      name = f'sub_assembly_{str(i)}'
      RunMegahit(os.path.join(out_folder, name), samples, cpus)
  elif type == 'single:
    for samples, i in enumerate(list_of_list_samples):
      name = f'single_assembly_{samples[0]}'
      RunMegahit(os.path.join(out_folder, name), samples, cpus)
  trimmed_reads_forward = glob.glob()

def RunMegahit(folder_path_name, samples, cpus):
  forward_list = [f'{sample}_val_1.fq.gz' for sample in samples]
  reverse_list = [f'{sample}_val_2.fq.gz' for sample in samples]
  args = f'megahit --presets meta-large --k-min 27 --k-max 127 --k-step 10 --min-contig-len 1000 -t {cpus} -1 {forward_list} -2 {reverse_list} -o {folder_path_name}'
  subprocess.call(args, shell = True)

def ProcessAssembly(out_folder, list_of_list_samples, type, cpus):
  if type == 'co':
    name = 'coassembly'
    RunMegahit(os.path.join(out_folder, name), list_of_list_samples[0], cpus)
  elif type == 'sub:
    for samples, i in enumerate(list_of_list_samples):
      name = f'sub_assembly_{str(i)}'
      RunMegahit(os.path.join(out_folder, name), samples, cpus)
  elif type == 'single:
    for samples, i in enumerate(list_of_list_samples):
      name = f'single_assembly_{samples[0]}'
      RunMegahit(os.path.join(out_folder, name), samples, cpus)
    

################################################################################
#################           MAIN                ################################
################################################################################    
def Main():
  # Create an argument parser
  parser = argparse.ArgumentParser(description="List files in a folder")

  # Add the folder path argument
  parser.add_argument("-illu", "--illumina_folder", type=str,
                        help="Path to the folder as a string", required=True)
  parser.add_argument("-nano", "--nanopore_folder", type=str,
                        help="Path to the folder as a string")
  parser.add_argument("-n", "--name", type=str,
                        help="Name of the results folder (_results will be added at the end)", required=True)
  parser.add_argument("-t", "--type", type=str,
                        help="co=coassembly, sub=subcoassemblies, single=singleassemblies", required=True)
  parser.add_argument("-m", "--metadata_file", type=str,
                        help="Path to the metadata tsv file", required=True)
  parser.add_argument("-c", "--cpus", type=int,
                        help="Number of cpus to use for multiprocessing-compatible tasks", required=True)
  parser.add_argument("--skippreprocessing", action='store_true',
                        help="To add if you want to skip preprocessing")

  # Parse argument, create results folder:
  args = parser.parse_args()
  out_folder = f'{args.name}_results'
  os.makedirs(out_folder)

  # Initiate the log file
  with open(f'{results_folder_name}/log.txt', 'w') as log:
    log.write(f"Log file for the run {results_folder_name}, time and date: {datetime.datetime.now().strftime('%I:%M%p on %B %d, %Y')}" + '\n\n')

  # Check if it is an hybrid assembly or short reads only
  is_hybrid = False
  if args.nanopore_folder is not None
    is_hybrid = True

  ################################################################################
  # 1. READS TRIMMING
  os.makedirs(f'{out_folder}/trimmed_reads')
  os.makedirs(f'{out_folder}/trimmed_reads/short_reads')
  short_reads_samples = ListReadsSamples(args.illumina_folder)
  ProcessTrimGalore(out_folder, short_reads_samples)

  long_reads_samples = []
  if is_hybrid:
    os.makedirs(f'{out_folder}/trimmed_reads/long_reads')
    long_reads_samples = ListReadsSamples(args.nanopore_folder)
    ProcessPorechopChopper(out_folder, long_reads_samples)

  ################################################################################
  # 2. ASSEMBLIES (one co-assembly, several smaller co-assemblies, individual assemblies)
  assemblies_folder = f'{out_folder}/assemblies'
  os.makedirs(assemblies_folder)
  metadata = LoadMetadata(args.metadata_file)

  if 'co' in args.type:
    samples_to_process = EvaluateSamples('co', short_reads_samples, long_reads_samples, metadata)
    ProcessAssembly(assemblies_folder, samples_to_process['coassembly'], 'co', args.cpus)
  if 'sub' in args.type:
    samples_to_process = EvaluateSamples('sub', short_reads_samples, long_reads_samples, metadata)
    ProcessAssembly(assemblies_folder, samples_to_process['subcoassemblies'], 'sub', args.cpus)
  if 'single' in args.type:
    samples_to_process = EvaluateSamples('single', short_reads_samples, long_reads_samples, metadata)
    ProcessAssembly(assemblies_folder, samples_to_process['singleassemblies'], 'single', args.cpus)

  ################################################################################
  # 3. Abundance and Functional annotation
  os.makedirs(f'{out_folder}/metagenome')

  ################################################################################
  # 4. MAGs recovery
  os.makedirs(f'{out_folder}/mags')










