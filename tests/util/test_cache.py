from opendrift.util.cache import file_cache

def test_file_cache(tmp_path):
    global side
    side = 1

    @file_cache(path = tmp_path / "test", timeout = None)
    def fun():
        global side

        side += 1

        return 2

    assert fun() == 2
    assert side == 2

    assert fun() == 2
    assert side == 2

    @file_cache(path = tmp_path / "test", timeout = None)
    def fun2():
        return 3

    assert fun2() == 2

