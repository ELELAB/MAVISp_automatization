import importlib.util
import os
import re
from pathlib import Path


ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX'
SCRIPT_DIR = Path(__file__).resolve().parent


def get_chains(pdb):
  from Bio.PDB import PDBParser

  parser = PDBParser(QUIET=True)
  structure = parser.get_structure('', pdb)
  chains = [c.id for c in structure.get_chains()]
  return chains


def read_pdb_numbering(pdb_path):
    numbering = {}
    with Path(pdb_path).open() as handle:
        for line in handle:
            if line[:6] != 'ATOM  ' and not (line[:6] == 'HETATM' and line[17:20] == 'MSE'):
                continue
            # Same fixed-width chain/residue fields and MSE handling as the model parser.
            numbering.setdefault(line[21:22], {})[line[22:27].strip()] = None
    return numbering


def infer_offset(pdb_path, chain, offset=None, metadata_only=False):
    if len(chain) != 1:
        raise ValueError('Postprocessing requires one selected PDB chain.')
    if metadata_only:
        residues = list(read_pdb_numbering(pdb_path).get(chain, {}))
    else:
        from protein_mpnn_utils import alt_parse_PDB

        parsed = alt_parse_PDB(str(pdb_path), chain)[0]
        residues = parsed['resn_list'] if parsed['seq'] else []
    if not residues:
        raise ValueError('No protein residues found for chain %s.' % chain)
    try:
        numbers = [int(number) for number in residues]
    except ValueError:
        raise ValueError('Residue insertion codes require an explicit mapping; one offset is insufficient.') from None
    if offset is not None:
        return offset

    match = re.search(r'_(\d+)-(\d+)$', Path(pdb_path).stem)
    if match is None:
        raise ValueError('Use a PDB filename ending in _START-END.pdb or provide --offset.')
    start, end = map(int, match.groups())
    if start < 1 or end < start:
        raise ValueError('Invalid residue range in the PDB filename.')
    first, last = min(numbers), max(numbers)
    # The parser keeps gaps; either original or 1-based domain numbering can use this offset.
    if (first == start and last <= end) or (first == 1 and last <= end - start + 1):
        return start
    raise ValueError(
        'Filename range %d-%d does not match chain %s numbering %d-%d; check the input or provide --offset.'
        % (start, end, chain, first, last)
    )


def main(cfg, args):
    import pandas as pd
    from omegaconf import OmegaConf
    from datasets import Mutation
    from protein_mpnn_utils import alt_parse_PDB
    from thermompnn_benchmarking import get_trained_model
    from SSM import get_ssm_mutations
    from postprocess_thermompnn import postprocess

    # define config for model loading
    config = {
        'training': {
            'num_workers': 8,
            'learn_rate': 0.001,
            'epochs': 100,
            'lr_schedule': True,
        },
        'model': {
            'hidden_dims': [64, 32],
            'subtract_mut': True,
            'num_final_layers': 2,
            'freeze_weights': True,
            'load_pretrained': True,
            'lightattn': True,
            'lr_schedule': True,
        }
    }

    cfg = OmegaConf.merge(config, cfg)

    # load the chosen model and dataset
    models = {
        "ThermoMPNN": get_trained_model(model_name=args.model_path,
                                        config=cfg, override_custom=True)
    }

    input_pdb = args.pdb
    pdb_id = os.path.basename(input_pdb).rstrip('.pdb')

    datasets = {
        pdb_id: args.pdb
    }

    raw_pred_df = pd.DataFrame(columns=['Model', 'Dataset', 'ddG_pred', 'position', 'wildtype', 'mutation',])
    row = 0
    for name, model in models.items():
        model = model.eval()
        model = model.cuda()
        for dataset_name, dataset in datasets.items():
            if len(args.chain) < 1:  # if unspecified, take first chain
                chain = get_chains(input_pdb)[0]
            else:
                chain = args.chain
            mut_pdb = alt_parse_PDB(input_pdb, chain)
            mutation_list = get_ssm_mutations(mut_pdb[0])
            final_mutation_list = []

            # build into list of Mutation objects
            for n, m in enumerate(mutation_list):
                if m is None:
                    final_mutation_list.append(None)
                    continue
                m = m.strip()  # clear whitespace
                wtAA, position, mutAA = str(m[0]), int(str(m[1:-1])), str(m[-1])

                assert wtAA in ALPHABET, f"Wild type residue {wtAA} invalid, please try again with one of the following options: {ALPHABET}"
                assert mutAA in ALPHABET, f"Wild type residue {mutAA} invalid, please try again with one of the following options: {ALPHABET}"
                mutation_obj = Mutation(position=position, wildtype=wtAA, mutation=mutAA,
                                        ddG=None, pdb=mut_pdb[0]['name'])
                final_mutation_list.append(mutation_obj)

            pred, _ = model(mut_pdb, final_mutation_list)

            for mut, out in zip(final_mutation_list, pred):
                if mut is not None:
                    col_list = ['ddG_pred', 'position', 'wildtype', 'mutation', 'pdb', 'chain']
                    val_list = [out["ddG"].cpu().item(), mut.position, mut.wildtype,
                                mut.mutation, mut.pdb.strip('.pdb'), chain]
                    for col, val in zip(col_list, val_list):
                        raw_pred_df.loc[row, col] = val

                    raw_pred_df.loc[row, 'Model'] = name
                    raw_pred_df.loc[row, 'Dataset'] = dataset_name
                    row += 1

    print(raw_pred_df)
    output_path = args.output if args.output is not None else Path("ThermoMPNN_inference_%s.csv" % pdb_id)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    raw_pred_df.to_csv(output_path)
    if args.postprocess_output is not None:
        args.postprocess_output.parent.mkdir(parents=True, exist_ok=True)
        postprocess(output_path, args.postprocess_output, args.offset)


def parse_args(argv=None):
    import argparse

    model_paths = {'default': SCRIPT_DIR / "thermoMPNN_default.pt"}
    parser = argparse.ArgumentParser(description='Run ThermoMPNN and write raw and MAVISp CSVs.')
    parser.add_argument('--pdb', required=True, help='Trimmed input PDB to use for inference')
    parser.add_argument('--chain', type=str, default='A', help='Chain in input PDB to use.')
    parser.add_argument(
        '--gpu', type=int, metavar='INDEX',
        help='CUDA index among visible GPUs (e.g. 0 or 1); omitted uses the current CUDA default'
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help='Show the plan and check files/package availability without loading a model or writing outputs'
    )
    parser.add_argument(
        '--output', type=Path,
        help='Raw prediction CSV path (default: ThermoMPNN_inference_<pdb_id>.csv in the working directory)'
    )
    processing_options = parser.add_mutually_exclusive_group()
    processing_options.add_argument(
        '--postprocess-output', type=Path,
        help='MAVISp CSV path (default: raw filename with _postprocessed.csv)'
    )
    processing_options.add_argument(
        '--raw-only', action='store_true', help='Write only the raw CSV and skip postprocessing'
    )
    parser.add_argument(
        '--offset', type=int, help='Override the residue offset read from the input PDB filename'
    )
    model_options = parser.add_mutually_exclusive_group()
    model_options.add_argument(
        '--model', choices=model_paths, default='default',
        help='Named ThermoMPNN model (default: default)'
    )
    model_options.add_argument(
        '--model-path', '--model_path', type=str,
        help='Path to a ThermoMPNN checkpoint instead of a named model'
    )

    args = parser.parse_args(argv)
    if args.gpu is not None and args.gpu < 0:
        parser.error('--gpu must be a non-negative CUDA device index')
    if args.model_path is None:
        args.model_path = str(model_paths[args.model])
    args.model_path = str(Path(args.model_path).expanduser().resolve())
    args.pdb = str(Path(args.pdb).expanduser().resolve())
    if not Path(args.pdb).is_file():
        parser.error('Input PDB not found: %s' % args.pdb)
    pdb_id = os.path.basename(args.pdb).rstrip('.pdb')
    args.output = (args.output or Path('ThermoMPNN_inference_%s.csv' % pdb_id)).expanduser().resolve()
    if args.raw_only:
        if args.offset is not None:
            parser.error('--offset cannot be used with --raw-only')
        if args.dry_run:
            chains = list(read_pdb_numbering(args.pdb))
            if not args.chain and chains:
                args.chain = chains[0]
            if args.chain not in chains:
                parser.error('No protein residues found for chain %s.' % args.chain)
    else:
        args.postprocess_output = (
            args.postprocess_output or args.output.with_name(args.output.stem + '_postprocessed.csv')
        ).expanduser().resolve()
        if not args.chain:
            chains = list(read_pdb_numbering(args.pdb)) if args.dry_run else get_chains(args.pdb)
            if not chains:
                parser.error('No chains found in the input PDB.')
            args.chain = chains[0]
        try:
            args.offset = infer_offset(args.pdb, args.chain, args.offset, metadata_only=args.dry_run)
        except ValueError as error:
            parser.error(str(error))
    if args.output == args.postprocess_output:
        parser.error('Raw and postprocessed outputs must use different paths.')
    protected_paths = {
        Path(args.pdb), Path(args.model_path).expanduser().resolve(),
        SCRIPT_DIR / 'local.yaml', SCRIPT_DIR / 'vanilla_model_weights' / 'v_48_020.pt',
    }
    for output in (args.output, args.postprocess_output):
        if output in protected_paths:
            parser.error('Output paths must not overwrite the input PDB, configuration or model weights.')
        if output is not None and output.is_dir():
            parser.error('Output path is a directory: %s' % output)
    return args


def dry_run(args):
    print('DRY RUN: no model loaded, no GPU computation, no files or folders written.')
    print('Shared code/config folder: %s' % SCRIPT_DIR)
    print('Input PDB: %s' % args.pdb)
    print('Chain: %s' % (args.chain or '(first chain)'))
    if args.gpu is not None:
        print('Requested GPU: cuda:%d (index among visible GPUs; not checked)' % args.gpu)
    else:
        print('Requested GPU: current CUDA default (visibility inherited; not checked)')
    print('ThermoMPNN checkpoint: %s' % args.model_path)
    print('Raw CSV would be written to: %s' % args.output)
    if args.postprocess_output is not None:
        print('Postprocessed CSV would be written to: %s' % args.postprocess_output)
        print('Residue offset: %d' % args.offset)
    else:
        print('Postprocessing: skipped (--raw-only)')

    files = [
        Path(args.pdb), Path(args.model_path), SCRIPT_DIR / 'local.yaml',
        SCRIPT_DIR / 'vanilla_model_weights' / 'v_48_020.pt',
    ]
    helpers = [
        'datasets.py', 'protein_mpnn_utils.py', 'thermompnn_benchmarking.py',
        'train_thermompnn.py', 'transfer_model.py', 'SSM.py', 'model_utils.py', 'cache.py',
        'postprocess_thermompnn.py',
    ]
    files.extend(SCRIPT_DIR / name for name in helpers)
    errors = []
    for path in files:
        if not path.is_file() or not os.access(path, os.R_OK):
            errors.append('Required file missing or unreadable: %s' % path)
        elif path.stat().st_size == 0:
            errors.append('Required file is empty: %s' % path)
    for output in (args.output, args.postprocess_output):
        if output is None:
            continue
        parent = output.parent
        while not parent.exists() and parent != parent.parent:
            parent = parent.parent
        if not parent.is_dir() or not os.access(parent, os.W_OK | os.X_OK):
            errors.append('Cannot create output under: %s' % parent)
        elif output.exists() and not os.access(output, os.W_OK):
            errors.append('Existing output is not writable: %s' % output)
        elif output.exists():
            print('Existing output would be overwritten: %s' % output)

    packages = ['numpy', 'pandas', 'torch', 'omegaconf', 'Bio', 'tqdm', 'pytorch_lightning', 'torchmetrics', 'wandb']
    missing = []
    for package in packages:
        try:
            available = importlib.util.find_spec(package) is not None
        except (ImportError, ValueError):
            available = False
        if not available:
            missing.append(package)
    if missing:
        errors.append('Python packages missing in this environment: %s' % ', '.join(missing))
    for error in errors:
        print('ERROR: %s' % error)
    print('GPU availability, checkpoint/config compatibility and prediction success are not checked.')
    print('Preflight checks: %s' % ('FAILED' if errors else 'OK'))
    return 1 if errors else 0


if __name__ == "__main__":
    args = parse_args()
    if args.dry_run:
        raise SystemExit(dry_run(args))
    import torch
    from omegaconf import OmegaConf

    if args.gpu is not None:
        torch.cuda.set_device(args.gpu)
        print('Selected GPU: cuda:%d (index among visible GPUs)' % args.gpu)
    cfg = OmegaConf.load(str(SCRIPT_DIR / "local.yaml"))
    cfg.platform.thermompnn_dir = str(SCRIPT_DIR)
    if args.postprocess_output is not None:
        print('Postprocessing offset: %d' % args.offset)
    with torch.no_grad():
        main(cfg, args)
