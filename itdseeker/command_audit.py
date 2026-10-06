"""Durable audit of external processes; argv is authoritative, display is POSIX quoted."""
import json
import logging
import os
from pathlib import Path
import shlex
import shutil
from datetime import datetime, timezone

_log_path = None


def configure_command_log(path):
    """Append across reruns; timestamps distinguish invocations. None disables file output."""
    global _log_path
    _log_path = Path(path) if path is not None else None
    if _log_path is not None:
        _log_path.parent.mkdir(parents=True, exist_ok=True)


def record_command(argv, *, stdin, stdout, stderr):
    """Record the actual argument vector and stream routing before process launch."""
    argv = [os.fspath(arg) for arg in argv]
    record = dict(
        timestamp=datetime.now(timezone.utc).isoformat(),
        argv=argv, executable=shutil.which(argv[0]), cwd=os.getcwd(),
        shell=False, command=shlex.join(argv),
        stdin=stdin, stdout=stdout, stderr=stderr,
    )
    logging.getLogger(__name__).info(
        "External command: %s [stdin=%s; stdout=%s; stderr=%s]",
        record['command'], stdin, stdout, stderr,
    )
    if _log_path is not None:
        with _log_path.open('a', encoding='utf-8') as handle:
            handle.write(json.dumps(record, ensure_ascii=False) + '\n')
    return record
