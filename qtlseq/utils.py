
import os
import sys
from datetime import datetime


def time_stamp():
    return '[QTL-seq:{}]'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

def clean_cmd(cmd):
    return ' '.join(cmd.split())
