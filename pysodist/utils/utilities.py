from datetime import datetime as dt
import sys


def log(msg, outfile=None):
    msg = '{} --> {}'.format(dt.now().strftime('%Y-%m-%d %H:%M:%S'), msg)
    print(msg)
    sys.stdout.flush()
    if outfile is not None:
        try:
            with open(outfile, 'a') as f:
                f.write(msg + '\n')
        except Exception as e:
            log(e)


def vlog(msg, outfile=None):
    msg = '{} -v-> {}'.format(dt.now().strftime('%Y-%m-%d %H:%M:%S'), msg)
    print(msg)
    sys.stdout.flush()
    if outfile is not None:
        try:
            with open(outfile, 'a') as f:
                f.write(msg + '\n')
        except Exception as e:
            vlog(e)
