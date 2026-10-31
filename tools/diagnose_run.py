"""Run a check with GIL-safe Python stack diagnostics, without signal handlers.

Some Sage/cysignals builds crash when faulthandler's watchdog inspects frames
asynchronously. This optional CI helper reads frames in an ordinary Python
thread instead. It neither terminates the check nor changes its exit status.
"""
import argparse
from pathlib import Path
import runpy
import sys
import threading
import traceback


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=float, default=120)
    parser.add_argument('kind', choices=('module', 'script'))
    parser.add_argument('target')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.interval <= 0:
        parser.error('--interval must be positive')
    stopped = threading.Event()
    main_thread = threading.get_ident()

    def report():
        while not stopped.wait(args.interval):
            frame = sys._current_frames().get(main_thread)
            print('\nDiagnostic snapshot (check continues):', file=sys.stderr, flush=True)
            if frame is not None:
                traceback.print_stack(frame, limit=18, file=sys.stderr)
            sys.stderr.flush()
            del frame

    reporter = threading.Thread(target=report, name='check-diagnostics', daemon=True)
    reporter.start()
    sys.argv = [args.target, *args.arguments]
    try:
        if args.kind == 'module':
            sys.path.insert(0, str(Path.cwd()))
            runpy.run_module(args.target, run_name='__main__', alter_sys=True)
        else:
            sys.path.insert(0, str(Path(args.target).resolve().parent))
            runpy.run_path(args.target, run_name='__main__')
    finally:
        stopped.set()
        reporter.join(timeout=1)


if __name__ == '__main__':
    main()
