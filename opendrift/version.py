__version__ = "1.11.11"


def git_describe():
    """
    Return git version if available.
    """
    import os.path
    from subprocess import check_output, DEVNULL

    path = os.path.dirname(__file__)
    args = [
        "git", "-C", path, "describe", "--tags", "--abbrev=7", "--dirty",
        "--broken"
    ]

    try:
        version = check_output(args, cwd=path, stderr=DEVNULL).decode().strip()
        return version
    except:
        return None

def version_or_git():
    v = git_describe()
    if v is None:
        return __version__
    else:
        return "%s / %s" % (__version__, v)

