import argparse
from logging import info
import os
import tempfile


def wrap(string):
    import re
    return re.sub("\\s{2,}", " ", string).strip()


def create_snapshot_dataset(args):
    create_dataset_cmd = f"bq --apilog=false --project_id={args.project} mk {args.dataset_name}"
    pipe = os.popen(create_dataset_cmd)
    info(pipe.read())
    ret = pipe.close()

    if ret:
        raise RuntimeError(f"Unexpected exit code from dataset creation: {ret}")

    fd, schema_file = tempfile.mkstemp(suffix=".json")

    try:
        with open(fd, 'w') as f:
            schema = """[
              {
                "name": "snapshot_name",
                "type": "STRING",
                "mode": "REQUIRED"
              },
              {
                "name": "integration_workflow_name",
                "type": "STRING",
                "mode": "REQUIRED"
              },
              {
                "name": "gvs_workflow_name",
                "type": "STRING",
                "mode": "REQUIRED"
              }
            ]
            """
            f.write(schema)

        create_table_cmd = wrap(f"""

        bq
          --apilog=false
          --project_id={args.project}
          mk
          --table {args.project}:{args.dataset_name}.integration_snapshot_metadata
          {schema_file}

        """)

        pipe = os.popen(create_table_cmd)
        info(pipe.read())
        ret = pipe.close()

        if ret:
            raise RuntimeError(f"Unexpected exit code from snapshot metadata table creation: {ret}")
    finally:
        os.remove(schema_file)


def build_argument_parser():
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Create integration test snapshot metadata dataset and table')
    parser.add_argument('--project', type=str, default='gvs-internal',
                        help='GCP project in which to create the snapshot metadata dataset')
    parser.add_argument('--dataset-name', type=str, help='Name of the snapshot metadata dataset to create',
                        required=True)
    return parser


def configure_logging():
    import logging
    import sys
    # https://stackoverflow.com/a/14058475
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root.addHandler(handler)


if __name__ == '__main__':
    configure_logging()
    arg_parser = build_argument_parser()
    parsed_args = arg_parser.parse_args()
    create_snapshot_dataset(parsed_args)
