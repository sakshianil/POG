nextflow.enable.dsl=2

// Define input channels


params {
  input_dir = "data/proteomes"
  output_dir = "results"
  signalp_script = "scripts/run_signalp_pipeline.py"
  deeptm_script = "scripts/run_deeptm_pipeline.py"
  wideeffhunter_script = "scripts/run_wideeffhunter.sh"
  pooe_script = "scripts/run_pooe_kaggle.sh"
  ortholog_script = "scripts/fetch_orthologs.sh"
  upstream_script = "scripts/fetch_upstream_protein.sh"
  motif_script = "scripts/motif_analysis.sh"
}

workflow {

  downloadProteomes()
  runSignalP()
  runDeepTM()
  runWideEffHunter()
  runPOOE()
  fetchOrthologs()
  extractSequences()
  runMotifAnalysis()
}

process downloadProteomes {
  output:
    path "${params.input_dir}/"
  """
  echo "Proteomes should be placed in ${params.input_dir}"
  """
}

process runSignalP {
  input:
    path params.input_dir
  output:
    path "signalp_output"
  script:
    """
    python3 ${params.signalp_script}
    """
}

process runDeepTM {
  input:
    path "signalp_output"
  output:
    path "deeptm_output"
  script:
    """
    python3 ${params.deeptm_script}
    """
}

process runWideEffHunter {
  input:
    path params.input_dir
  output:
    path "${params.input_dir}/WideEff_Hunter_outputs"
  script:
    """
    bash ${params.wideeffhunter_script}
    """
}

process runPOOE {
  output:
    path "POOE_output"
  script:
    """
    bash ${params.pooe_script}
    """
}

process fetchOrthologs {
  output:
    path "orthologs.csv"
  script:
    """
    bash ${params.ortholog_script}
    """
}

process extractSequences {
  output:
    path "upstream_sequences"
    path "protein_sequences"
  script:
    """
    bash ${params.upstream_script}
    """
}

process runMotifAnalysis {
  input:
    path "upstream_sequences"
    path "protein_sequences"
  output:
    path "motif_results"
  script:
    """
    bash ${params.motif_script}
    """
}

