cwlVersion: v1.2

$graph:
  # CommandLineTool: protein-quest search uniprot -> uniprot_accs.txt
  - class: CommandLineTool
    id: search_uniprot
    label: protein-quest search uniprot
    baseCommand: [protein-quest, search, uniprot]
    requirements:
      InlineJavascriptRequirement: {}
      NetworkAccess:
        networkAccess: true
    inputs:
      taxon_id:
        type: ["null", string]
        inputBinding:
          prefix: --taxon-id
          position: 1
      reviewed:
        type: ["null", boolean]
        inputBinding:
          prefix: --reviewed
          position: 2
      no_reviewed:
        type: ["null", boolean]
        inputBinding:
          prefix: --no-reviewed
          position: 3
      subcellular_location_uniprot:
        type: ["null", string]
        inputBinding:
          prefix: --subcellular-location-uniprot
          position: 4
      subcellular_location_go:
        type:
          - "null"
          - type: array
            items: string
        inputBinding:
          prefix: --subcellular-location-go
          position: 5
      molecular_function_go:
        type:
          - "null"
          - type: array
            items: string
        inputBinding:
          prefix: --molecular-function-go
          position: 6
      limit:
        type: int
        default: 10000
        inputBinding:
          prefix: --limit
          position: 7
      timeout:
        type: int
        default: 1800
        inputBinding:
          prefix: --timeout
          position: 8
      output_filename:
        type: string
        default: uniprot_accs.txt
        inputBinding:
          position: 9
    outputs:
      uniprot_accessions:
        type: File
        outputBinding:
          glob: $(inputs.output_filename)

  # CommandLineTool: protein-quest search pdbe uniprot_accs.txt -> pdbe.csv
  - class: CommandLineTool
    id: search_pdbe
    label: protein-quest search pdbe
    baseCommand: [protein-quest, search, pdbe]
    requirements:
      InlineJavascriptRequirement: {}
      NetworkAccess:
        networkAccess: true
    inputs:
      accs_file:
        type: File
        inputBinding:
          position: 1
      limit:
        type: int
        default: 10000
        inputBinding:
          prefix: --limit
          position: 2
      timeout:
        type: int
        default: 1800
        inputBinding:
          prefix: --timeout
          position: 3
      output_filename:
        type: string
        default: pdbe.csv
        inputBinding:
          position: 4
    outputs:
      pdb_ids_csv:
        type: File
        outputBinding:
          glob: $(inputs.output_filename)

  # CommandLineTool: protein-quest retrieve pdbe pdbe.csv -> downloads-pdbe/
  - class: CommandLineTool
    id: retrieve_pdbe
    label: protein-quest retrieve pdbe
    baseCommand: [protein-quest, retrieve, pdbe]
    requirements:
      InlineJavascriptRequirement: {}
      NetworkAccess:
        networkAccess: true
    inputs:
      pdbe_csv:
        type: File
        inputBinding:
          position: 1
      output_dirname:
        type: string
        default: downloads-pdbe
        inputBinding:
          position: 2
      max_parallel_downloads:
        type: int
        default: 5
        inputBinding:
          prefix: --max-parallel-downloads
          position: 3
    outputs:
      out_dir:
        type: Directory
        outputBinding:
          glob: $(inputs.output_dirname)

  # CommandLineTool: protein-quest filter chain pdbe.csv downloads-pdbe -> filtered-chains/
  - class: CommandLineTool
    id: filter_chain
    label: protein-quest filter chain
    baseCommand: [protein-quest, filter, chain]
    requirements:
      InlineJavascriptRequirement: {}
    inputs:
      chains:
        type: File
        inputBinding:
          position: 1
      input_dir:
        type: Directory
        inputBinding:
          position: 2
      output_dirname:
        type: string
        default: filtered-chains
        inputBinding:
          position: 3
      scheduler_address:
        type: ["null", string]
        inputBinding:
          prefix: --scheduler-address
          position: 4
    outputs:
      out_dir:
        type: Directory
        outputBinding:
          glob: $(inputs.output_dirname)

  # CommandLineTool: protein-quest filter residue --min-residues 100 --max-residues 1000 filtered-chains -> filtered/
  - class: CommandLineTool
    id: filter_residue
    label: protein-quest filter residue
    baseCommand: [protein-quest, filter, residue]
    requirements:
      InlineJavascriptRequirement: {}
    inputs:
      min_residues:
        type: int
        default: 0
        inputBinding:
          prefix: --min-residues
          position: 1
      max_residues:
        type: int
        default: 10000000
        inputBinding:
          prefix: --max-residues
          position: 2
      input_dir:
        type: Directory
        inputBinding:
          position: 3
      output_dirname:
        type: string
        default: filtered
        inputBinding:
          position: 4
      write_stats:
        type: ["null", string]
        inputBinding:
          prefix: --write-stats
          position: 5
    outputs:
      out_dir:
        type: Directory
        outputBinding:
          glob: $(inputs.output_dirname)

  # Workflow: chain all steps ensuring each output feeds the next one
  - class: Workflow
    id: main
    label: protein-quest search + retrieve + filter workflow
    inputs: {}
    steps:
      search_uniprot:
        run: '#search_uniprot'
        in:
          taxon_id:
            default: "9606"
          reviewed:
            default: true
          subcellular_location_uniprot:
            default: nucleus
          subcellular_location_go:
            default: ["GO:0005634"]
          molecular_function_go:
            default: ["GO:0003677"]
          limit:
            default: 100
          output_filename:
            default: uniprot_accs.txt
        out: [uniprot_accessions]

      search_pdbe:
        run: '#search_pdbe'
        in:
          accs_file: search_uniprot/uniprot_accessions
          limit:
            default: 100
          timeout:
            default: 1800
          output_filename:
            default: pdbe.csv
        out: [pdb_ids_csv]

      retrieve_pdbe:
        run: '#retrieve_pdbe'
        in:
          pdbe_csv: search_pdbe/pdb_ids_csv
          output_dirname:
            default: downloads-pdbe
          max_parallel_downloads:
            default: 5
        out: [out_dir]

      filter_chain:
        run: '#filter_chain'
        in:
          chains: search_pdbe/pdb_ids_csv
          input_dir: retrieve_pdbe/out_dir
          output_dirname:
            default: filtered-chains
        out: [out_dir]

      filter_residue:
        run: '#filter_residue'
        in:
          input_dir: filter_chain/out_dir
          min_residues:
            default: 100
          max_residues:
            default: 1000
          output_dirname:
            default: filtered
        out: [out_dir]

    outputs:
      uniprot_accs:
        type: File
        outputSource: search_uniprot/uniprot_accessions
      pdbe_csv:
        type: File
        outputSource: search_pdbe/pdb_ids_csv
    #   downloads_pdbe:
    #     type: Directory
    #     outputSource: retrieve_pdbe/out_dir
    #   filtered_chains:
    #     type: Directory
    #     outputSource: filter_chain/out_dir
      filtered:
        type: Directory
        outputSource: filter_residue/out_dir

# Original commands for reference:
# protein-quest search uniprot \
#     --taxon-id 9606 \
#     --reviewed \
#     --subcellular-location-uniprot nucleus \
#     --subcellular-location-go GO:0005634 \
#     --molecular-function-go GO:0003677 \
#     --limit 100 \
#     uniprot_accs.txt
# protein-quest search pdbe uniprot_accs.txt pdbe.csv
# protein-quest retrieve pdbe pdbe.csv downloads-pdbe/
# protein-quest filter chain \
#     pdbe.csv \
#     ./downloads-pdbe ./filtered-chains
# protein-quest filter residue  \
#     --min-residues 100 \
#     --max-residues 1000 \
#     ./filtered-chains ./filtered

# Run with

# cwltool --prov cwlprov workflow.cwl
# runcrate convert cwlprov -o crate
# Browse crate/ro-crate-metadata.json in https://ro-crate.ldaca.edu.au/explorer
