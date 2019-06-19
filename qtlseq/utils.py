
import os
import sys
from datetime import datetime


def time_stamp():
    return '[QTL-seq:{}]'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

def clean_cmd(cmd):
    return ' '.join(cmd.split())

def call_log(out_dir, name, cmd):
    print(time_stamp(), 
          '!!ERROR!! {}\n'.format(cmd), 
          flush=True)

    print('please check below:\n')

    with open('{}/log/{}.log'.format(out_dir, name)) as log:
        for line in log:
            print(line, end='')
