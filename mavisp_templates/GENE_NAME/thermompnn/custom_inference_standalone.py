import os
import sys
from pathlib import Path


ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX'


def get_chains(pdb):
  from Bio.PDB import PDBParser

  parser = PDBParser(QUIET=True)
  structure = parser.get_structure('', pdb)
  chains = [c.id for c in structure.get_chains()]
  return chains


def main(cfg, args):
    # Add the external installation before importing its helper modules.
    sys.path.append(str(args.thermompnn_dir))
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
            if not mut_pdb[0]['seq']:
                raise ValueError('No protein residues found for chain %s.' % chain)
            if args.postprocess_output is not None:
                try:
                    for number in mut_pdb[0]['resn_list']:
                        int(number)
                except ValueError:
                    raise ValueError('Residue insertion codes cannot be mapped with one offset.') from None
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

    parser = argparse.ArgumentParser(description='Run ThermoMPNN and write raw and MAVISp CSVs.')
    parser.add_argument('--pdb', required=True, help='Trimmed input PDB to use for inference')
    parser.add_argument(
        '--thermompnn-dir', required=True, type=Path,
        help='External installation containing helper modules, postprocessor and model weights'
    )
    parser.add_argument('--chain', type=str, default='A', help='Chain in input PDB to use.')
    parser.add_argument(
        '--gpu', type=int, metavar='INDEX',
        help='CUDA index among visible GPUs; omitted uses the current CUDA default'
    )
    parser.add_argument(
        '--output', type=Path,
        help='Raw CSV path (default: ThermoMPNN_inference_<pdb_id>.csv in the working directory)'
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
        '--offset', type=int,
        help='Starting residue added to zero-based positions; required unless --raw-only'
    )
    model_options = parser.add_mutually_exclusive_group()
    model_options.add_argument(
        '--model', choices=['default'], default='default',
        help='Named ThermoMPNN model (default: default)'
    )
    model_options.add_argument(
        '--model-path', '--model_path', type=str,
        help='Path to a ThermoMPNN checkpoint instead of a named model'
    )

    args = parser.parse_args(argv)
    args.thermompnn_dir = args.thermompnn_dir.expanduser().resolve()
    if not args.thermompnn_dir.is_dir():
        parser.error('ThermoMPNN installation directory not found: %s' % args.thermompnn_dir)
    if args.gpu is not None and args.gpu < 0:
        parser.error('--gpu must be a non-negative CUDA device index')
    if args.model_path is None:
        args.model_path = str(args.thermompnn_dir / 'thermoMPNN_default.pt')
    args.model_path = str(Path(args.model_path).expanduser().resolve())
    if not Path(args.model_path).is_file():
        parser.error('ThermoMPNN checkpoint not found: %s' % args.model_path)
    args.pdb = str(Path(args.pdb).expanduser().resolve())
    if not Path(args.pdb).is_file():
        parser.error('Input PDB not found: %s' % args.pdb)
    pdb_id = os.path.basename(args.pdb).rstrip('.pdb')
    args.output = (args.output or Path('ThermoMPNN_inference_%s.csv' % pdb_id)).expanduser().resolve()
    if args.raw_only:
        if args.offset is not None:
            parser.error('--offset cannot be used with --raw-only')
    else:
        if args.offset is None:
            parser.error('--offset is required for postprocessing; it is not inferred from the PDB filename')
        if args.chain and len(args.chain) != 1:
            parser.error('Postprocessing requires one selected PDB chain')
        args.postprocess_output = (
            args.postprocess_output or args.output.with_name(args.output.stem + '_postprocessed.csv')
        ).expanduser().resolve()
    if args.output == args.postprocess_output:
        parser.error('Raw and postprocessed outputs must use different paths.')
    protected_paths = {
        Path(args.pdb), Path(args.model_path),
        args.thermompnn_dir / 'vanilla_model_weights' / 'v_48_020.pt',
        Path(__file__).resolve(),
    }
    for output in (args.output, args.postprocess_output):
        if output in protected_paths:
            parser.error('Output paths must not overwrite the input PDB, inference script or model weights.')
        if output is not None and output.is_dir():
            parser.error('Output path is a directory: %s' % output)
    return args


if __name__ == "__main__":
    args = parse_args()
    import torch
    from omegaconf import OmegaConf

    if args.gpu is not None:
        torch.cuda.set_device(args.gpu)
        print('Selected GPU: cuda:%d (index among visible GPUs)' % args.gpu)
    cfg = OmegaConf.create({
        'platform': {'accel': 'gpu', 'thermompnn_dir': str(args.thermompnn_dir)}
    })
    if args.postprocess_output is not None:
        print('Postprocessing offset: %d' % args.offset)
    with torch.no_grad():
        main(cfg, args)
