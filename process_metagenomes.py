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
  parser.add_argument("-c", "--cpus", type=str,
                        help="Number of cpus to use for multiprocessing-compatible tasks", required=True)
  parser.add_argument("--skippreprocessing", action='store_true',
                        help="To add if you want to skip preprocessing")

  # Parse argument, create results folder:
  args = parser.parse_args()
  out_folder = f'{args.name}_results'

  # Initiate the log file
  with open(f'{results_folder_name}/log.txt', 'w') as log:
    log.write(f"Log file for the run {results_folder_name}, time and date: {datetime.datetime.now().strftime('%I:%M%p on %B %d, %Y')}" + '\n\n')

  # Check if it is an hybrid assembly or short reads only
  is_hybrid = False
  if args.nanopore_folder is not None
    is_hybrid = True

  ################################################################################
  # 1. READS TRIMMING
  short_reads_samples = ListReadsSamples(args.illumina_folder)
  ProcessTrimGalore(out_folder, short_reads_samples)

  long_reads_samples = []
  if is_hybrid:
    long_reads_samples = ListReadsSamples(args)
    ProcessPorcheopChopper(out_folder, long_reads_samples)

  ################################################################################
  # 2. ASSEMBLIES (one co-assembly, several smaller co-assemblies, individual assemblies)
  metadata = LoadMetadata(args.metadata_file)
  samples_to_process = EvaluateSamples(short_reads_samples, long_reads_samples, metadata)

  if 'co' in args.type:
    ProcessAssembly(out_folder, samples_to_process['coassembly'])
  if 'sub' in args.type:
    ProcessAssembly(out_folder, samples_to_process['subcoassemblies'])
  if 'single' in args.type:
    ProcessAssembly(out_folder, samples_to_process['singleassemblies'])

  ################################################################################
  # 3. Abundance and Functional annotation

  ################################################################################
  # 4. MAGs recovery









