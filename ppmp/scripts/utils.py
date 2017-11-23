import string
import os
import errno


def alpha_label(iterable):
    unique = set(iterable)
    LABELS = [i+n for i in string.ascii_uppercase for n in string.digits]
    mapping = {v: LABELS[i] for i, v in enumerate(unique)}
    label = [mapping[i] for i in iterable]
    return label, mapping


def handle_missing_dirs(path):
    try:
        os.makedirs(os.path.dirname(path))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    return os.path.abspath(path)
