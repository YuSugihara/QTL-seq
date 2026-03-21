
import os
import sys
import glob
import shutil
from datetime import datetime


def time_stamp():
    return '[QTL-seq:{}]'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

def clean_cmd(cmd):
    return ' '.join(cmd.split())

def check_required_commands(commands, context):
    missing = [cmd for cmd in commands if shutil.which(cmd) is None]
    if len(missing) == 0:
        return

    sys.stderr.write(('{} Error: Missing required command(s) for {}: {}.\n').format(time_stamp(),
                                                                                      context,
                                                                                      ', '.join(missing)))
    sys.stderr.write(('{} Please check your PATH and installation.\n').format(time_stamp()))
    sys.stderr.write(('{} Example checks: {}\n\n').format(time_stamp(),
                                                          ' ; '.join(['which {}'.format(cmd) for cmd in missing])))
    sys.exit(1)


def call_log(out_dir, name, cmd, returncode=None, hint=None):
    log_path = '{}/log/{}.log'.format(out_dir, name)

    print(time_stamp(), 
          '!!ERROR!! {}\n'.format(cmd), 
          flush=True)

    if returncode is not None:
        print(time_stamp(),
              'Exit status: {}'.format(returncode),
              flush=True)
        if returncode == 127:
            print(time_stamp(),
                  'It looks like a command in this step is not available in PATH.',
                  flush=True)

    print(time_stamp(),
          'Expected log file: {}'.format(log_path),
          flush=True)

    print('please check below:\n')

    if os.path.exists(log_path):
        with open(log_path) as log:
            for line in log:
                print(line, end='')
    else:
        print(time_stamp(),
              'Log file was not found.',
              flush=True)
        related_logs = sorted(glob.glob('{}/log/{}*.log'.format(out_dir, name)))
        if len(related_logs) > 0:
            print(time_stamp(),
                  'Related log files:',
                  flush=True)
            for related_log in related_logs:
                print(' - {}'.format(related_log), flush=True)

    if hint is not None:
        print('\n{} Hint: {}'.format(time_stamp(), hint), flush=True)
