"""Filter local structured diagnostics without loading complete logs into memory."""
import argparse
from collections import deque
import json
import logging
from pathlib import Path
import re
from .diagnostics import log_directory


def read_events(directory, *, limit=30, level='WARNING', **filters):
    threshold = logging.getLevelName(level)
    if not isinstance(threshold, int) or limit < 1:
        raise ValueError('Invalid log level or limit')
    selected = deque(maxlen=limit)
    root = Path(directory)
    # Rotated files are oldest first, followed by the active file.
    paths = [root/f'events.jsonl.{i}' for i in (3, 2, 1)] + [root/'events.jsonl']
    for path in paths:
        try:
            with path.open(encoding='utf-8') as stream:
                for line in stream:
                    try: record = json.loads(line)
                    except ValueError: continue  # tolerate a final interrupted write
                    if not isinstance(record, dict): continue
                    priority = logging.getLevelName(record.get('level', 'INFO'))
                    if not isinstance(priority, int) or priority < threshold: continue
                    if any(record.get(key) != value for key, value in filters.items() if value is not None): continue
                    selected.append(record)
        except FileNotFoundError:
            continue
    return list(selected)


def main(argv=None):
    parser = argparse.ArgumentParser(description='Inspect BioMolExplorer diagnostics; technical messages retain their original language.')
    parser.add_argument('--directory', type=Path, help='Log root (default: BIOMOL_LOG_DIR or logs)')
    parser.add_argument('--job', help='Inspect one job directory')
    parser.add_argument('--run', help='Filter run_id')
    parser.add_argument('--stage', help='Filter stage_id')
    parser.add_argument('--operation')
    parser.add_argument('--level', choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], default='WARNING')
    parser.add_argument('--limit', type=int, default=30)
    parser.add_argument('--json', action='store_true', help='Return a JSON array')
    args = parser.parse_args(argv)
    if args.limit < 1 or args.limit > 10000: parser.error('--limit must be between 1 and 10000')
    directory = args.directory or log_directory()
    if args.job:
        if not re.fullmatch(r'[A-Za-z0-9_-]{1,100}', args.job): parser.error('Invalid job identifier')
        directory = directory/'jobs'/args.job
    events = read_events(directory, limit=args.limit, level=args.level,
                         job_id=args.job, run_id=args.run, stage_id=args.stage, operation=args.operation)
    if args.json:
        print(json.dumps(events, ensure_ascii=False, indent=2))
    else:
        print(f'Diagnostics: {directory}')
        if not events: print('No matching events. Check the directory, filters and log level.')
        for record in events:
            identity = ' '.join(f'{k}={record[k]}' for k in ('run_id','stage_id','job_id','pair','command_id') if record.get(k))
            print(f"\n{record['timestamp']} {record['level']} {record.get('error_code',record.get('event','message'))} {identity}")
            print(record.get('message', ''))
            if record.get('action'): print('Action: '+record['action'])
            if record.get('configuration'): print('Configuration: '+record['configuration'])
            if record.get('diagnostic_path'): print('Details: '+record['diagnostic_path'])
            source=record.get('source', {})
            if source: print(f"Source: {source.get('file')}:{source.get('line')}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
