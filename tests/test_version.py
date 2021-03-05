import opendrift

def test_git_describe():
    v = opendrift.version.git_describe()
    print("git:", v)

def test_version_or_git():
    v = opendrift.version.version_or_git()
    print("full version:", v)
    assert isinstance(v, str)

